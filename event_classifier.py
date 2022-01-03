# 03 Jan 22
# Classifier finds the optimal minor and major parent for recombinant.

import argparse
import os
from pathlib import PurePath


class classifier():

    folderToParse = ''

    rec_events_files = []
    seq_events_files = []
    alignment_files = []

    maxGenomeLength = 0

    def __init__(self, rootPath):

        # Get folder name to parse
        self.folderToParse = rootPath

        # Walk through the folder and find all files recombination event files, sequence event files and
        # alignment files and store them in lists.
        self.rec_events_files, self.seq_events_files, self.alignment_files = self.getFileNames()

        # Find the maximum genome length from the alignment files.
        # self.maxGenomeLength = self.get_max_genome_length()

    def getFileNames(self):

        # Folder to walk through and find the recom
        targetFolder = PurePath(self.folderToParse)

        filesToParse = []
        for paths in os.walk(targetFolder):
            for files in paths[2]:
                if files.endswith('.fa'):
                    filesToParse.append(PurePath(
                        r"C:\Users\joshc\OneDrive\University\Masters\RDP-ML", paths[0] + '/' + files)
                    )

    # def get_max_genome_length(self):
    #     # Find the longest genome in the folder being parsed.
    #     try:
    #         alignment = open(self.alignment_fileName, "r")
    #         alignment.seek(0)
    #         _ = alignment.readline()
    #         data = alignment.readline()
    #     finally:
    #         self.maximum_genome_length = len(data)


def getPath():
    # This function is used to get the folder path from command line.
    # For future pipeline use.

    # Define command line argument parser
    parser = argparse.ArgumentParser(
        description='Parse Recombination Information from SantaSim'
    )

    # Add arguments for command line
    parser.add_argument('-f', dest='folderToParse', type=str,
                        help='The root folder to parse', required=True)

    # Parse Events
    args = parser.parse_args()

    # Returns the root folder to parse
    return(args.folderToParse)


if __name__ == '__main__':

    # The line below will be used to get fileNames from pipeline later
    # folderToParse = getPath()

    # Currently used for testing purposes.
    folderToParse = 'data/'

    # Create classifier class by initialising file paths
    parser = classifier(folderToParse)

    # Gets maximum genome length.
    parser.get_max_genome_length()
