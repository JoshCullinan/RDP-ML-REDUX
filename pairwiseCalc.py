import argparse
import textdistance
import os
import psutil
import logging
import numpy as np
from multiprocessing import Pool
from Bio import AlignIO
from scipy import stats
from pathlib import Path


def hamSim(i, x, j, x_prime):
    '''
    Calculates the normalised hamming similarity between supplied sequences.
    '''
    sim = textdistance.hamming.normalized_similarity(x, x_prime)
    return [i, j, sim]

def ReadData(alignmentName):
    '''
    Reads in the sequences. 
    '''
    alignment = AlignIO.read(alignmentName, 'fasta')
    seqs = np.array([str(alig.seq) for alig in alignment])

    return seqs

def multicoreSimilarity(seqs, threads):
    '''
    Calculates the pariwise similarity between all the sequences in an alignment. Quick 'n dirty multicore implemention because vectorising is hard and this is quick.
    '''
    length =  len(seqs)

    pairwiseSim = np.empty((length,length))
    boolarr = np.full((length,length), True)
    np.fill_diagonal(boolarr, False)

    pairwiseArr = []
    for i, x in enumerate(seqs):
        for j, x_prime in enumerate(seqs):
            pairwiseArr.append((i, x, j, x_prime))

    with Pool(threads) as pool:
        results = pool.starmap(hamSim, pairwiseArr)

    for res in results:
        pairwiseSim[res[0], res[1]] = res[2]

    avSim = pairwiseSim.mean(where=boolarr)
    maximum = np.amax(pairwiseSim, initial = 0, where=boolarr)
    minimum = np.amin(pairwiseSim, initial = 1, where= boolarr)

    return avSim, minimum, maximum, pairwiseSim


def cullData(pairwiseSim):
    """
    Aim to remove highest similarity and lowest similarity sequences from the alignment until pairwise similarities reached.
    """
    length = len(pairwiseSim)
    boolarr = np.full((length,length), True)
    np.fill_diagonal(boolarr, False)

    avSim = pairwiseSim.mean(where=boolarr)
    maximum = np.amax(pairwiseSim, initial = 0, where=boolarr)
    minimum = np.amin(pairwiseSim, initial = 1, where= boolarr)

    print(f'Initial values\nAverage: {avSim:.4f}, Max: {maximum:.4f}, Min: {minimum:.4f}.')

    output = pairwiseSim

    while output.mean(where=boolarr) < 0.80:

        if np.amin(output, initial = 1, where= boolarr) < 0.75:
            np.fill_diagonal(output, 1)
            locMin = stats.mode(np.argmin(output, axis = 0, keepdims=True), axis=None)[0]
            output = np.delete(output, int(locMin), axis = 0)
            output = np.delete(output, int(locMin), axis = 1)

            length = len(output)
            boolarr = np.full((length,length), True)
            np.fill_diagonal(boolarr, False)


        if np.amax(output, initial = 0, where=boolarr) > 0.97:
            np.fill_diagonal(output, 0)
            locMax = stats.mode(np.argmax(output, axis = 0, keepdims=True), axis=None)[0]
            output = np.delete(output, int(locMax), axis = 0)
            output = np.delete(output, int(locMax), axis = 1)
        
            length = len(output)
            boolarr = np.full((length,length), True)
            np.fill_diagonal(boolarr, False)

        avSim = output.mean(where=boolarr)
        maximum = np.amax(output, initial = 0, where=boolarr)
        minimum = np.amin(output, initial = 1, where= boolarr)

        print(f'Average Similarity: {avSim:.4f}, Maximum Similarity: {maximum:.4f}, Minimum similarity: {minimum:.4f}.')
    
    return output

def getFileNames(folderToParse = ''):
    """
     Walk through the folder and find all alignment files ending with .fa and store them in lists.
    """
    # List to fill
    rdpFiles = []

    # Folder to walk through and find the recom
    targetFolder = Path(folderToParse)

    # Walk through the folder and find all of the files in it.
    # Then determine if the file matches an aligment, recombination event or sequence event file
    # If matching add to the relavent list as type Path
    for paths in os.walk(targetFolder):
        for files in paths[2]:
            if files.endswith('.fa'):
                rdpFiles.append(
                    Path(paths[0] + '/' + files))

    return rdpFiles

def findPairwise(rdpFiles, folder, threads):
    """
    Uses the multicoreSimilarity function to find the average, min, max similarities for an entire directory of alignments.
    """

    print(f'Started now. {len(rdpFiles)} to complete.')

    for i, files in enumerate(rdpFiles):
        seqs = ReadData(files)
        
        print(f'Current CPU utilistaion: {psutil.cpu_percent(interval=None, percpu=False)}')
        print(psutil.cpu_percent(interval=None, percpu=True))
        avSim, minimum, maximum, _ = multicoreSimilarity(seqs, threads)
        
        fileName =  Path(f'Pairwise_Similarity_{folder.name}.csv')

        try:
            with open(fileName, 'x') as f:
                f.write(f'Name,AverageSimilarity,MinimumSimilarity,MaximumSimilarity\n')
        except FileExistsError:
            logging.info('File already exists')
    
        with open(fileName, 'a') as f:
            f.write(f'{files.name},{avSim},{minimum},{maximum}\n')

        print(f'Completed {i+1} out of {len(rdpFiles)}. Completed {files.name}.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate Average Normalised Hamming similarity of an alignment file using -f and of a directory using -d')
    parser.add_argument('-f', type=str, help='Alignment file.')
    parser.add_argument('-d', type=str, help='Directory of Alignments.')
    parser.add_argument('-t', type=str, help='Threads to use.')
    args = parser.parse_args()

    if args.t:
        threads = int(args.t)
    else:
        threads = int(os.cpu_count())

    if args.f:
        logging.info('Parsing an alignment file')
        f = Path(args.f)
        if f.exists() and not f.is_dir():
            seqs = ReadData(f)
            avsim, minimum, maximum, pairwiseSim = multicoreSimilarity(seqs, threads)
            print(f"Average Similarity: {avsim:.04f}\nMinimum similarity: {minimum:.04f}\nMax Similartiy: {maximum:.04f}")
        else:
            logging.error('The path given is a directory or the file does not exist.')

    if args.d:
        logging.info('Parsing an entire directory of alignments')
        folder = Path(args.d)
        if folder.exists() and folder.is_dir():
            rdpFiles = getFileNames(folderToParse=folder)
            findPairwise(rdpFiles, folder, threads)
        else:
            logging.error('The path given is a file or does not exist.')