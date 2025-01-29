#09.01.2022
#Joshua Cullinan
#SantaSim Simulation ease of access.
#Parses santaSim XMLS and runs all the combinations of variables that have been given.
from decimal import localcontext
import os
import logging
import argparse
import time
import shutil
import re
import ast
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from multiprocessing import Pool
from itertools import product

# java -jar /home/clljos001/santa.jar -generationCount=$GC -recombinationProbability=$RP -mutationRate=$MR -sampleSize=$SS /home/clljos001/xml/$file.xml

def folderfnc(key):

    alig =  local.parent / Path(f'alignment_{key}.fa')
    recom = local.parent / Path(f'recombination_events_alignment_{key}.txt')
    seq = local.parent / Path(f'sequence_events_map_alignment_{key}.txt')

    try:
        alig_dst = Path(outputs / Path(f'alignment_{key}.fa'))
        shutil.move(alig, alig_dst)

        recom_dst = Path(outputs / Path(f'recombination_events_alignment_{key}.txt'))
        shutil.move(recom, recom_dst)

        seq_dst = Path(outputs / Path(f'sequence_events_map_alignment_{key}.txt'))
        shutil.move(seq, seq_dst)
    except:
        logging.error(outputs / Path(f'alignment_{key}.fa'))
        logging.error(outputs / Path(f'recombination_events_alignment_{key}.txt'))
        logging.error(outputs / Path(f'sequence_events_map_alignment_{key}.txt'))
        logging.error('Failed to move alignment files')

def execute(input):

    cmds, variableNames = input
    maxHeapSize = 4

    javaCMD = [
        'java',
        '-jar',
        '-Xms1G',
        f'-Xmx{maxHeapSize}G',
        'santa.jar',
    ]

    localCount = cmds[-1]

    key = f'Number:{cmds[-1]}_{now}'
    for value, varName in zip(cmds[0:-1], variableNames[0:-1]):
        key = key + str(value) + '_'
        javaCMD.append(f'-{varName}={value}')
    
    localCount = cmds[-1]

    #Create the logfile key and alignment details.
    key = key.rstrip('_')
    an = f'alignment_{key}.fa'
    logFile = outputs/Path(f'Alignment_Log_{key}.txt')
    
    #Add alignment name to the javaCMD.
    javaCMD.append(f'-alignmentName={an}')

    #Last edits to work.
    javaCMD.append(XML)

    completed = False
    while completed != True:
        with open(logFile, 'a+') as lf:
            lf.write(f'Working on {localCount}.\n')
            print(f'Working on {localCount}.')
            try:
                subprocess.run(javaCMD, shell=False, stdout=lf, stderr=subprocess.STDOUT, check=True)
                
                try:
                    folderfnc(key)
                except:
                    logging.error('Folder Function Error')
                    lf.write('\nError in folder function.')

                completed = True
                lf.write('\nFin.')

            except subprocess.CalledProcessError as e:
                logging.error(f'Error occured in JAVA for {localCount} -- {an}. Reason logged to file.')
                
                #If error exists write to file
                if e.cmd:
                    lf.write(str(e.cmd))
                if e.output:
                    lf.wite(str(e.output))
                if e.stderr:
                    lf.write(str(e.stderr))
                if e.stdout:
                    lf.write(str(e.stdout))

                    
                #Delete alignment and start again.
                an = Path(an)
                if an.exists():
                    os.remove(an)

                #Increase max heap size up to a specific limit.
                if maxHeapSize < 15:
                    maxHeapSize += 1
                
                javaCMD[3] = f'-Xmx{maxHeapSize}G'
                print(javaCMD)
                print(f'Max Heap Size increased to {maxHeapSize}GBs')

def parseXML():
    with open(XML, 'r+') as f:
        start = f.tell()
        data = list()
        #Read the file for the comments containing the variables.
        while f.readline().startswith('<!--'):
            f.seek(start)
            data.append(f.readline().strip('<!->\n= '))
            start = f.tell()

    variables = {}
    #Create a dictionary of all the variables and their values that we are going to use.
    for values in data:
        value = re.search(r'\[.*\]', values).group()
        varName = re.search(r'.*(?= \=)', values).group()
        variables[varName] = ast.literal_eval(value)

    return variables


if __name__ == '__main__':
    """
    This file takes in an XML and parses it to get the variables from the XML file and sets up a
    multiprocessing approach to simulate all of them.
    """
    #Parse arguments.
    argParser = argparse.ArgumentParser(description='Simulate SantaSim (Cluster or Local)')
    argParser.add_argument('-o', dest='output', help='Where is the output to be saved. Optional flag usually to denote cluster usage.')
    argParser.add_argument('-t', dest='threads', help='How many simulations to run.')
    argParser.add_argument('-xml', dest='xml', help='XML file to use.', required=True)
    args = argParser.parse_args()
    
    #Set some global variables
    #Used for unique file name. 
    global now
    now = datetime.now().strftime('%Y-%m-%d-%H_%M_%S_')
    
    #Find local file path.
    global local
    local = Path(__file__).resolve()

    #Set output directory. Local or specified.
    global outputs
    if args.output:
        outputs = Path(args.output)
    else:
        outputs = Path(local.parent / Path('outputs'))

    # Make outputs folder
    if not outputs.exists():
        print(f'Making outputs folder at {outputs}.')
        outputs.mkdir(parents=True)

    #Setting XML as global based on argParse
    global XML 
    XML = args.xml
    variables = parseXML()

    #Set number of threads to use. Default is 1.
    global threads
    if args.threads:
        threads = int(args.threads)
    else:
        threads = 1

    Store = variables.copy()

    #Set up lists of permutations to use.
    VariableNames = list()
    Vals = list()
    while True:
        try:
            currVarName, currVal = variables.popitem()
            #I'm doing it this way because dictionary pops from the bottom (weird).
            VariableNames.append(currVarName)
            Vals.append(currVal)
        except KeyError:
            break
    
    #Add the count variable name for later use
    VariableNames.append('count')

    #Get all the combinations of the variables given. 
    #Will adjust to whatever size is given. 
    expectedRuns = list(product(*Vals))
    count = 0
    for i in range(len(expectedRuns)):
        count+=1
        expectedRuns[i] = expectedRuns[i] + (i+1,)

    with Pool(threads) as p:
        print(f'Starting now. {len(expectedRuns)} runs to complete.')
        p.map(execute, zip(expectedRuns, [VariableNames] * len(expectedRuns)))
    
    # When finished copy and rename the XML file into the output directory to store that runs data.
    shutil.copyfile(XML, Path(outputs/Path(r'#' + XML.name)))