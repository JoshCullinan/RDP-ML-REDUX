# 03 Jan 22
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


class classifier:

    # Lists of paths for the three file types we're parsing
    rec_events_path = Path()
    seq_events_path_path = Path()
    alignment_path = Path()

    # Recombination events and sequence events files
    alignment = pd.DataFrame
    rec_events = pd.DataFrame
    seq_events = pd.DataFrame

    # The longest genome in the alignment files.
    maxGenomeLength = 0
    numberOfAlignments = 0

    # Generation matrix is a pandas dataframe (for now) that is [number of alignments x max genome length] ([row x columns])
    generationMatrix = pd.DataFrame

    # Dictionaries
    seqmap_dict = {}
    events_dict = {}
    events_map = {}
    inv_seqmap_dict = defaultdict(set)

    def __init__(self, alig, rec, seq):

        # Get relavent files that will be used in the parsing
        self.alignment_path = Path(alig)
        self.rec_events_path = Path(rec)
        self.seq_events_path = Path(seq)

        # Read in files function
        self.readFiles()
        # Find the maximum genome length from the alignment file.
        self.get_max_genome_length()
        # Create dictionaries used in generation matrix
        self.create_dictionaries()
        # Create generation count matrix
        self.createGenerationMatrix()

    def get_max_genome_length(self):
        # Find the longest genome amongst the alignments
        self.maxGenomeLength = max(self.alignment.iloc[:, 0].apply(len))

        # Whilst I have the alignment file open find the number of alignments in it (between 100 and 200 I think)
        # Divide by two because pandas interprets the .fa files weirdly
        self.numberOfAlignments = math.ceil(len(self.alignment.index) / 2)

    def readFiles(self):
        # Load in the files as pandas dataframes
        self.alignment = pd.read_csv(self.alignment_path)

        # Read in Recombination events file
        self.rec_events = pd.read_csv(
            self.rec_events_path,
            sep=r"*",
            usecols=["EventNum", "Breakpoints", "Generation"],
        )

        # Remove the brackets surrounding breakpoints
        self.rec_events.Breakpoints = self.rec_events.Breakpoints.str.strip("[]")
        # Split breakpoints into start and end, drop breakpoints
        self.rec_events[["Start", "End"]] = self.rec_events.Breakpoints.str.split(
            ",",
            expand=True,
        )

        # Read in sequence events map
        self.seq_events = pd.read_csv(
            self.seq_events_path, delimiter="*", index_col="Sequence"
        )

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

    def createGenerationMatrix(self):
        # Generation matrix is a pandas dataframe (for now) that is [number of alignments x max genome length] ([row x columns])
        self.generationMatrix = pd.DataFrame(
            np.zeros(
                shape=(self.numberOfAlignments, self.maxGenomeLength), dtype=np.int8
            )
        )

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

            self.generationMatrix.iloc[box, start : (end + 1)] = np.full(
                (len(seqs), end + 1 * fix - start), event, dtype=int
            )

        # Just for testing purposes.
        print(self.generationMatrix)


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
