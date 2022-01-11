# 07 Jan 22
# Classifier finds the optimal minor and major parent for recombinant.

import argparse
import collections
import os
from pathlib import Path
import pandas as pd
import numpy as np
import math
from collections import defaultdict
import ast
from Bio import SeqIO, AlignIO
import distance
from intervaltree import Interval, IntervalTree
import re

class classifier:

    # Lists of paths for the three file types we're parsing
    rec_events_path = Path
    seq_events_path_path = Path
    alignment_path = Path

    # Recombination events and sequence events files
    alignment = dict
    rec_events = pd.DataFrame
    seq_events = pd.DataFrame

    # The longest genome in the alignment files.
    maxGenomeLength = 0
    numberOfSeqs = 0

    # Generation matrix is a numpy array that is [number of alignments x max genome length] ([row x columns])
    generationMatrix = np.array

    # Dictionaries
    seqmap_dict = dict
    events_dict = dict
    events_map = dict
    inv_seqmap_dict = defaultdict(set)

    # Gap dict
    gaps = defaultdict(set)

    def __init__(self, alig, rec, seq):

        # Get relavent files that will be used in the parsing
        self.alignment_path = Path(alig)
        self.rec_events_path = Path(rec)
        self.seq_events_path = Path(seq)        
        self.major_parents = {}
        self.minor_parents = {}

        # Read in files function
        self.readFiles()
        # Create dictionaries used in generation matrix
        self.create_dictionaries()
        # Find posistion of Gap characters in the sequences
        self.getGaps()
        # Create generation count matrix
        self.createGenerationMatrix()
        # Calc Parents
        self.calcParents()
        # Output to csv.
        self.output()


    def readFiles(self):
        # JOSH: Trying the AlignIO feature from BioPython as they have get max length and number of Seq Fnc.
        # Useful if faster than my function for this.
        self.alignment = AlignIO.read(self.alignment_path, 'fasta')
        self.maxGenomeLength = self.alignment.get_alignment_length()
        self.numberOfSeqs = self.alignment.__len__()
        self.alignment = SeqIO.to_dict(self.alignment)
        #self.alignment = SeqIO.to_dict(SeqIO.parse(self.alignment_path, 'fasta'))

        #changing to just store sequence, dont need other entries that biopython stores in dictionary        
        for k, v in self.alignment.items():
            self.alignment[k] = v.seq

        # Read in Recombination events file
        self.rec_events = pd.read_csv(
            self.rec_events_path,
            sep=r"*",
            usecols=["EventNum", "Breakpoints", "Generation"],
        )

        # Remove the brackets surrounding breakpoints
        self.rec_events.Breakpoints = self.rec_events.Breakpoints.str.strip("[]")

        # Split breakpoints into start and end, drop breakpoints
        self.rec_events[["Start", "End"]] = self.rec_events.Breakpoints.str.split(",", expand=True,)

        # Read in sequence events map
        self.seq_events = pd.read_csv(self.seq_events_path, delimiter="*", index_col="Sequence")

    def create_dictionaries(self):
        # generating dictionaries from dataframes
        self.seqmap_dict = {
            i: ast.literal_eval(v)
            for i, v in enumerate(self.seq_events["Events"].to_numpy(), 1)
        }
        self.events_dict = {
            event: ast.literal_eval(bp)
            for event, bp in zip(
                self.rec_events["EventNum"], self.rec_events["Breakpoints"]
            )
        }
        # Creating an inverted seqmap dictionary (event:sequences instead of sequence:events)
        # with key,value pairs: event, [sequences containing event]
        for key, value in self.seqmap_dict.items():
            for eventnum in value:
                self.inv_seqmap_dict[eventnum].add(key)

    def getGaps(self):
        for key, seq in self.alignment.items():
            currentGaps = ([pos for pos, char in enumerate(str(seq)) if char == '-'])
            self.gaps[int(key)] = currentGaps

    def createGenerationMatrix(self):
        # Generation matrix is a numpy array that is [number of alignments x max genome length] ([row x columns])

        #self.generationMatrix = np.zeros(shape=(self.numberOfSeqs, self.maxGenomeLength))
        self.generationMatrix = np.full(shape=(self.numberOfSeqs, self.maxGenomeLength), dtype='O', fill_value=0)
        for key, gaps in self.gaps.items():
            self.generationMatrix[key-1, [pos for pos in gaps]] = '-'

        for event, seqs in self.inv_seqmap_dict.items():

            start = int(self.rec_events.Start[self.rec_events.EventNum == int(event)])

            end = int(self.rec_events.End[self.rec_events.EventNum == int(event)])

            # Fixes indexing error involving maximum genome length
            fix = 1
            if end == self.maxGenomeLength:
                fix = 0

            box = []
            for x in iter(seqs):
                box.append(x - 1)

            self.generationMatrix[box, start : (end + 1)] = np.full(
                (len(seqs), end + 1 * fix - start), int(event), dtype=np.int32
                )

    def calcHammingDistance(self, seq1, seq2):
        
        # Set for instances where seq1 and seq2 are 0 after dropping -
        normalisedDistance = None

        # Find the index of the gap characters in both sequences
        seq1Gaps = ({pos for pos, char in enumerate(str(seq1)) if char == '-'})
        seq2Gaps = ({pos for pos, char in enumerate(str(seq2)) if char == '-'})

        # Find the union of these Gap characters
        union = seq1Gaps | seq2Gaps

        # Remove the characters at the indexes of the union. 
        # 'A-TT-G' (Gaps at 1 and 4)
        # 'GG-TAA' (Gap at 3)
        # Union is {1,3,4}
        # Output sequences will be 'ATG' and 'GTA'. Hamming distance will be 2. 
        # Normalised over the total sites where neither sequence has a gap (3)
        # 2/3 = 0.66

        for count, indx in enumerate(union):
            seq1 = seq1[:indx-count] + seq1[indx+1-count:]
            seq2 = seq2[:indx-count] + seq2[indx+1-count:]
        
        #returns the hamming distance between two sequences
        d = distance.hamming(seq1, seq2)

        #Length of seqs
        totalNucleotides = len(seq1)

        if totalNucleotides > 0:
            # Normalise over number of chars without gap characters
            normalisedDistance = d/totalNucleotides

        return normalisedDistance

    def findEventPositions(self):
        #we need to know where the "recombination event blocks" are, i.e. which sections of the alignment we need to compare sequences within to find parents
        #these sections are different for every event, and arent neccesarily continuous: these blocks can be overwritten by later events, creating fragmented blocks

        #will make a dictionary to store this information
        #format is as follows {key = event number, value = dictionary of sequences}
        #where the dictionary of sequences stores the nucleotide positions for each sequence
        #for example if positions 100-200 for sequence 5 and 10 are written as event 2543 in the generation matrix, the dictionary entry will be as follows:
        #{2543: {5: [100, 200], 10: [100, 200]}}

        #The point is to have a way to efficiently know where any given recombination block in the generation matrix is, 
        #without having to search it for all entries of a given event number
        block_dict = {x:{} for x in self.events_dict.keys()}       

        #extracting raw array
        gen_matrix = self.generationMatrix

        #iterate through generation matrix, where index = a tuple (sequence name, nucleotide position) 
        for index, entry in np.ndenumerate(gen_matrix):
            #entry of 0 means no recombination event touched this nucleotide, so we don't have to do anything
            if not entry == 0:
                seq = index[0]
                nucleotide_pos = index[1]
                
                if seq in block_dict[entry]:
                    #if entry exists, change as appropiate
                    min_max = block_dict[entry][seq]
                    max_pos = min_max[-1][-1]                    

                    if nucleotide_pos == max_pos:
                        #continuous, uninterrupted block of same event, thus just add +1 to maximum range
                        min_max[-1][-1] += 1
                        block_dict[entry][seq] = min_max
                    else:
                        #discontinuity, so need to make a new range of nucleotides
                        min_max.append([nucleotide_pos, nucleotide_pos+1])
                        block_dict[entry][seq] = min_max
                else:   
                    #if dictionary entry doesnt exist yet, create a new one          
                    block_dict[entry][seq] = [[nucleotide_pos, nucleotide_pos+1]]

        return block_dict

    def calculateParents(self, block_dict, range_flipper):
        #range_flipper = flips ranges, so True => major parents are calculated instead of minor
               
        #returns a dictionary with {key = event number, value = set(major parents)}
        #TODO HERE: how to find best major parent out of multiple options?
        #It is not guaranteed that all sequences will have the same best major parent

        #this variable (deleted_nucleotides) will keep track of nucleotides with higher event numbers than all current events under consideration,
        #thus if a nucleotide falls into this range, it shouldnt be considered for parent calculations
        #Format is {sequence_name: interval_tree}
        #interval tree holds all intervals of deleted nucleotides
        deleted_nucleotides = {}             
        parents = {} 
        #we traverse block_dict in reverse order, calculating parents and then adding the ranges traversed to deleted nucleotides
        #i.e. we start at the highest event number, calculate parents with the recombinant region,
        #then that recombinant region is added to deleted nucleotides. Since all future events will have a smaller event number (earlier generation),
        #these nucleotides shouldnt be considered for any parent calculations.
        for event_number, sequence_ranges_dict in reversed(block_dict.items()):

            #divide sequences into recombinant and potential parents
            sequences_in_block = set(sequence_ranges_dict.keys())
            sequences_not_in_block = set(range(self.numberOfSeqs)) - sequences_in_block      
            best_parents = []

            #calculate best parents for all sequences of the current recombination event
            for sequence, ranges in sequence_ranges_dict.items():  
                recombinant_seq = str(self.alignment[str(sequence+1)])  
                hamming_distances = {}               

                #calculate hamming distances for all potential parents
                for parent in sequences_not_in_block:
                    parent_seq = str(self.alignment[str(parent+1)])
                    seq1 = '' 
                    seq2 = ''

                    #if we do minor parents, take recombinant region and then remove intervals that have been deleted
                    if not range_flipper:
                        ranges_tree = IntervalTree.from_tuples(ranges)
                        if parent in deleted_nucleotides.keys():                                               
                            for j in deleted_nucleotides[parent]:
                                ranges_tree.chop(j.begin, j.end) 
                    #for major parents, its the same except we take the complement of the recombinant region (all regions not in the recombinant region)
                    else:
                        recombinant_region = IntervalTree.from_tuples(ranges)
                        ranges_tree = IntervalTree.from_tuples([[0, self.maxGenomeLength]])
                        for r in recombinant_region:
                            ranges_tree.chop(r.begin, r.end)
                        if parent in deleted_nucleotides.keys():                                               
                            for j in deleted_nucleotides[parent]:
                                ranges_tree.chop(j.begin, j.end) 

                    
                    #can extract sequence now, since we have the final region:
                    #recombinant region with all deleted nucleotides removed
                    for k in ranges_tree:
                        seq1 = seq1 + recombinant_seq[k.begin:k.end]
                        seq2 = seq2 + parent_seq[k.begin:k.end]

                    #calculating hamming distance
                    hamming_distance = self.calcHammingDistance(seq1, seq2)
                    
                    if hamming_distance != None:
                        hamming_distances[parent] = hamming_distance

                    #JOSH EDITED: Normalising in the function now due to dropping characters
                    # total_nucleotides = (len(seq1))
                        
                    # #now add the sequence that has been compared to, together with normalised total hamming distance
                    # if total_nucleotides > 0:
                    #     hamming_distances[parent] = hamming_distance/total_nucleotides

                #now all the distances have been calculated for this particular sequence, need to find minimum
                minimum_seq = min(hamming_distances, key=hamming_distances.get)
                #add this minimum hamming distance sequence, together with its hamming distance to best parents list 
                #so this output is: recombinant sequence, best parent, hamming distance               
                best_parents.append((sequence+1, minimum_seq+1, hamming_distances[minimum_seq]))

                #now add the nucleotides we have traversed to deleted nucleotides, these won't be considered in future events
                if sequence in deleted_nucleotides.keys():                                                            
                    deleted_nucleotides[sequence] = deleted_nucleotides[sequence].union(IntervalTree.from_tuples(ranges))                     
                    deleted_nucleotides[sequence].merge_overlaps()
                else:
                    deleted_nucleotides[sequence] = IntervalTree.from_tuples(ranges)

            parents[event_number] = best_parents  

        return parents

    def calcParents(self):
        #This function uses the generation matrix to calculate the best minor and major parents for each recombination event

        #we need to know where the "recombination event blocks" are, i.e. which sections of the alignment to compare to find parents
        #will make a dictionary to store this information, see function for more details on dictionary
        block_dict = self.findEventPositions()         
        #now we can use this dictionary to find the major parents
        print("Calculating minor parents...")
        self.minor_parents = self.calculateParents(block_dict, False)      
        print("Calculating major parents...")
        self.major_parents = self.calculateParents(block_dict, True)
        print("Done")       

    def output(self):
        print("Creating output file...")
        
        # Create unique key for the file name
        key = re.search(r'(?<=alignment_).*', self.alignment_path.name).group()[:-3]
        fileName = ("output/RPD_Output_" + key + '.csv')
        filePath = Path(fileName)
        
        try:
            os.makedirs('output')
        except FileExistsError:
            pass

        with open(fileName, "w") as g:
            header = ['SantaEventNumber', 'StartBP', 'EndBP', 'Recombinant', 'MinorParent', 'MajorParent'] 
            g.write('\t'.join(str(s) for s in header) + '\n')

        for events in self.minor_parents.keys():
            startBP, EndBP = self.events_dict[events]
            for minorTup, MajorTup in zip(self.minor_parents[events], self.major_parents[events]):
                recom = minorTup[0]
                minor = minorTup[1]
                major = MajorTup[1]
                
                with open(fileName, "a") as f:
                    content = [events, startBP, EndBP, recom, minor, major]
                    f.write('\t'.join(str(s) for s in content) + '\n')
        print('Done')

def getFilePaths():
    # This function is used to get the folder path from command line.
    # For future pipeline use.

    # Define command line argument parser
    parser = argparse.ArgumentParser(
        description="Parse Recombination Information from SantaSim"
    )

    # Add arguments for command line
    parser.add_argument(
        "-a",
        dest="alignment_path",
        type=str,
        help="recombination events file",
        required=True,
    )
    parser.add_argument(
        "-f",
        dest="recombination_path",
        type=str,
        help="recombination events file",
        required=True,
    )
    parser.add_argument(
        "-f",
        dest="sequence_path",
        type=str,
        help="sequence events map file",
        required=True,
    )

    # Parse Events
    args = parser.parse_args()

    # Returns the root folder to parse
    return (args.alignment_path, args.recombination_path, args.sequence_path)


if __name__ == "__main__":

    # The line below will be used to get fileNames from pipeline later
    # alignment_path, recombination_path, sequence_path = getFilePaths()

    # Currently used for testing purposes.
    alignment_path = "data/alignment_XML1-2500-0.01-12E-5-100-13.fa"
    recombination_path = "data/recombination_events_XML1-2500-0.01-12E-5-100-13.txt"
    sequence_path = "data/sequence_events_map_XML1-2500-0.01-12E-5-100-13.txt"

    # Create classifier class by initialising file paths
    parser = classifier(alignment_path, recombination_path, sequence_path)

    # XML1-2500-0.01-12E-5-100-13
    # XML5-4000-0.02-12E-5-50-4-3