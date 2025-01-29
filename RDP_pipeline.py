# 28/12/2021
# To run through the output data using RDP and parsing the results.
# Need to fix naming in output parser.

import os
import subprocess
from subprocess import run
from pathlib import Path
import sys
import re
import argparse


class parsing_script():

    def __init__(self):
        pass

    def script(self, folder):
        
        print(f"Does this folder exist? {folder.exists()}")

        filesToParse = []
        for paths in os.walk(folder):
            for files in paths[2]:
                if files.endswith('.fa'):
                    filesToParse.append(Path("C:/Users/joshc/OneDrive - University of Cape Town/University/Masters/RDP-ML-REDUX", paths[0] + '/' + files))
                # if files.endswith("SimVSRealCompare.csv"):
                #     alreadyDone.append(Path("C:/Users/joshc/OneDrive - University of Cape Town/University/Masters/RDP-ML-REDUX", paths[0] + '/' + files))


        for prog, f in enumerate(filesToParse):
            print(f"Currently parsing number {prog+1} out of {len(filesToParse)} - {f.name}")

            nameToWrite = f.parents[0] / (f.name + '_log.txt')
            with open(nameToWrite, 'a+') as b:
                b.write('Starting Now\n')

            # Look for files from previous simulation runs.
            key = re.search(r'(?<=alignment_).*', f.name).group()[:-3]

            # It took me 45 minutes to debug a spelling error. Instead of RDP it was RPD. JFC.
            rdp5ml = Path(
                f.parents[0] / (r"RPD_Output_" + key + r".rdp5ML"
                ))

            Done = Path(
                f.parents[0] / (r"alignment_" + key + r".faSimVsRealCompare.csv"
                ))

            if (rdp5ml.exists() and f.exists() and not Done.exists()):
                with open(nameToWrite, 'a') as b:
                    b.write(f'The rdp5ml file is:{str(rdp5ml)}.\nThe alignment file is: {str(f)}\n')  

                # newWD = Path("C:/Users/joshc/OneDrive - University of Cape Town/University/Masters/RDP-ML-REDUX/" / folder)
                
                cmds = "RDP5CL.exe" 
                try:

                    execute(["cd", folder, "&&", cmds, "-f", f.name, "-rdp5ml", rdp5ml.name, "-ds"])
                    
                    # result = run([cmds, "-f", f, "-rdp5ml", rdp5ml, "-ds"],
                    #                 capture_output=True,
                    #                 text=True,
                    #                 universal_newlines=True,
                    #                 bufsize=-1,
                    #                 check=True,
                    #                 timeout=600,
                    #                 cwd=newWD
                    #                 )

                    # with open(nameToWrite, 'a') as b:
                    #     b.write(result.stdout)
                    #     b.write(result.stderr)
                    
                except:
                    with open(nameToWrite, 'a') as b:
                        b.write("\nRDP threw an error")

                    print("Error in RDP execution for file: " + f.name)
            else:
                print("RDP5ML file present: " + str(rdp5ml.exists()))
                print("Alignment file present: " + str(f.exists()))
                print("Already parsed: " + str(Done.exists()))

            _ = os.system(
                'del 3seqTable, BinProbs, LastSave.rdp5, PairsScores, SCF, RDP5FSSRDP, tempfile, RDP5Redolist* > nul 2>&1& cd ..')

def execute(command):
    subprocess.check_call(command, shell=True, stdout=sys.stdout, stderr=subprocess.STDOUT)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process RDP5 files')
    parser.add_argument('-f', dest='folder', help='file path to parse')
    args = parser.parse_args()

    if not args.folder:
        print("Grrr give me a file...")
        raise FileNotFoundError

    folder = Path(args.folder)

    script = parsing_script()
    script.script(folder)