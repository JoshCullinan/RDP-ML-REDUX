import os
from pathlib import PurePath

folderToParse = ''

alignment_files = []
rec_events_files = []
seq_events_files = []


def getFileNames():
    # Walk through the folder and find all recombination event files, sequence event files and
    # alignment files and store them in lists.

    # Folder to walk through and find the recom
    targetFolder = PurePath(folderToParse)

    # Walk through the folder and find all of the files in it.
    # Then determine if the file matches an aligment, recombination event or sequence event file
    # If matching add to the relavent list as type Path
    for paths in os.walk(targetFolder):
        for files in paths[2]:
            if files.endswith('.fa'):
                alignment_files.append(
                    PurePath(paths[0] + '/' + files))

            if files.startswith('recombination_events'):
                rec_events_files.append(
                    PurePath(paths[0] + '/' + files))

            if files.startswith('sequence_events_map'):
                seq_events_files.append(
                    PurePath(paths[0] + '/' + files))
