import os
from pathlib import Path, Path
import event_classifier
import re
import gc
from multiprocessing import Pool

alignment_files = []

def getFileNames(folderToParse = ''):
    # Walk through the folder and find all recombination event files, sequence event files and
    # alignment files and store them in lists.

    # Folder to walk through and find the recom
    targetFolder = Path(folderToParse)

    # Walk through the folder and find all of the files in it.
    # Then determine if the file matches an aligment, recombination event or sequence event file
    # If matching add to the relavent list as type Path
    for paths in os.walk(targetFolder):
        for files in paths[2]:
            if files.endswith('.fa'):
                alignment_files.append(
                    Path(paths[0] + '/' + files))
    
def parsing_loop(alig):
    # total = len(alignment_files)

    # for count, alig in enumerate(alignment_files):
    key = re.search(r'(?<=alignment_).*', alig.name).group()[:-3]
    rec = Path(alig.parents[0] / ('recombination_events_' + key + '.txt'))
    seq = Path(alig.parents[0] / ('sequence_events_map_' + key + '.txt'))

    if alig.exists() and rec.exists() and seq.exists():
        # print(f'Parsing {count+1} out of {total}.')
        parse = event_classifier.classifier(alig, rec, seq)
        
        #Remove parse after use and create new.
        del parse
        gc.collect()
        
    else:
        print("The requested files don't exist")
        print('Alig: ' + str(alig.exists()))
        print('Rec: ' + str(rec.exists()))
        print('Seq: ' + str(seq.exists()))
        pass


if __name__ == '__main__':
    # Change the path below to your target path.
    getFileNames(folderToParse=Path(r'dataRaw/Test'))
    with Pool() as p:
        p.map(parsing_loop, alignment_files)
    