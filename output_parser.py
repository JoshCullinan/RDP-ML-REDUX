import os
import re
import argparse
import pandas as pd
from pathlib import Path, Path


rdpStatsFiles = []
rdpSimVReal = []

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
            if files.endswith('.faRecombIdentifyStats.csv'):
                rdpStatsFiles.append(
                    Path(paths[0] + '/' + files))
            if files.endswith('.faSimVSRealCompare.csv'):
                rdpSimVReal.append(
                    Path(paths[0] + '/' + files))

def validate_and_clean_triplets(processed_df, rdpFiles):
    """
    Validate triplets and remove those with more than one recombinant.
    
    Parameters:
    processed_df (pandas.DataFrame): The processed DataFrame
    
    Returns:
    pandas.DataFrame: Cleaned DataFrame with invalid triplets removed
    dict: Statistics about the cleaning process
    """
    rows_to_keep = []
    stats = {
        'original_triplets': len(processed_df) // 3,
        'removed_triplets': 0,
        'remaining_triplets': 0
    }
    
    for i in range(0, len(processed_df), 3):
        triplet = processed_df.iloc[i:i+3]
        recomb_count = triplet['is_recombinant'].sum()
        
        if recomb_count == 1:
            rows_to_keep.extend([i, i+1, i+2])
        else:
            stats['removed_triplets'] += 1
            print(f"Removing triplet at rows {i}-{i+2}: {recomb_count} recombinants found, in files {rdpFiles[0]} & {rdpFiles[1]}")
    
    cleaned_df = processed_df.iloc[rows_to_keep].copy()
    stats['remaining_triplets'] = len(cleaned_df) // 3
    
    return cleaned_df, stats


# def parser(filesToParse):
#     #Load in file
#     RecDF = pd.read_csv(filesToParse[0], index_col=False)
#     SimVRealDF = pd.read_csv(filesToParse[1], index_col=False)

#     # Strip leading whitespace
#     RecDF.columns = RecDF.columns.str.strip()
#     SimVRealDF.columns = SimVRealDF.columns.str.strip()

#     # # Drop un-needed data
#     # RecDF.drop(["Event","StartBP", "EndBP"], axis = 1, inplace = True)

#     file_name = f"output_test/ml_input_{folder.name}.txt"
#     file_path = Path(file_name)

#     # if file doesnt exist yet, create it and write header
#     if (not (file_path.exists())):
#         os.makedirs(os.path.dirname(file_name), exist_ok=True)
#         with open(file_name, "w") as g:
#             header = ['Recombinant'] + RecDF.columns.to_list()
#             g.write('\t'.join(str(s) for s in header) + '\n')

#     with open(file_name, "a") as f:
#         for col in range(0, len(RecDF), 3):
#             out = [RecDF.iloc[col+s].tolist() for s in range(3) ]
            
            
#             zipped = zip(out[0], out[1], out[2])
#             recombinant_vector = [1.0, 0.0, 0.0]

#             final_output = [tuple(recombinant_vector)] + [x for x in zipped]

#             f.write('\t'.join(str(y) for y in final_output) + '\n')

def save_processed_data(processed_df, output_path):
    """
    Save the processed DataFrame to a CSV file.
    
    Parameters:
    processed_df (pandas.DataFrame): The processed DataFrame to save
    output_path (str): Path where the CSV should be saved
    """
    if not output_path.exists():
        processed_df.to_csv(output_path, index=False)
    else:
        processed_df.to_csv(output_path, mode='a', header=False, index=False)

def process_recombination_data(recomb_stats_path, sim_compare_path):
    """
    Process recombination statistics and simulation comparison data to create a merged dataset
    with binary labels indicating recombinant (1) or parent (0) sequences.
    Each triplet in RecombIdentifyStats corresponds to one row in SimVSRealCompare.
    
    Parameters:
    recomb_stats_path (str): Path to the RecombIdentifyStats CSV file
    sim_compare_path (str): Path to the SimVSRealCompare CSV file
    
    Returns:
    pandas.DataFrame: Merged dataset with binary labels
    """
    # Read the CSV files
    recomb_stats = pd.read_csv(recomb_stats_path, index_col=False)
    sim_compare = pd.read_csv(sim_compare_path, index_col=False)

    #Strip white space from column headings
    sim_compare.columns = sim_compare.columns.str.strip()
    recomb_stats.columns = recomb_stats.columns.str.strip()
    

    recomb_stats.drop(["Event","StartBP", "EndBP"], axis = 1, inplace = True)

    # Create a new column to store the binary labels
    recomb_stats['is_recombinant'] = 0
    
    # Process data in triplets
    for i in range(0, len(sim_compare)):
        # Get the actual recombinant ID for this triplet
        actual_recomb = sim_compare.iloc[i]['ActualRecomb']
        
        # Get the corresponding three rows from recomb_stats
        start_idx = i * 3
        for j in range(3):  # Process each row in the triplet
            current_idx = start_idx + j
            
            if current_idx >= len(recomb_stats):
                print(f"Warning: Reached end of recomb_stats at index {current_idx}")
                break
                
            # Get the sequence IDs from ISeqs(A)
            if pd.isna(recomb_stats.iloc[current_idx]['ISeqs(A)']):
                continue
                
            seq_ids = str(recomb_stats.iloc[current_idx]['ISeqs(A)']).split('$')
            seq_ids = [int(id_) for id_ in seq_ids if id_.strip()]
            
            # Check if the actual recombinant is in this sequence
            if actual_recomb in seq_ids:
                recomb_stats.at[current_idx, 'is_recombinant'] = 1
            else:
                recomb_stats.at[current_idx, 'is_recombinant'] = 0

    return recomb_stats

def parsing_loop():
    #Does len of recombIdent == SimVCompare?
    total = len(rdpStatsFiles)
    if len(rdpStatsFiles) != len(rdpSimVReal):
        print("The number of files in the recombination event folder does not match the number of files in the sequence event folder")
        print(f"Recombination event files: {total}")
        print(f"Sequence event files: {len(rdpSimVReal)}")
        return

    for count, rdpFiles in enumerate(zip(rdpStatsFiles,rdpSimVReal)):
        # key = re.search(r'(?<=alignment_).*', alig.name).group()[:-3]

        if rdpFiles[0].exists():
            print(f'Parsing {count+1} out of {total}.')
            # parser(rdpFiles)
            processed_data = process_recombination_data(rdpFiles[0], rdpFiles[1])
            
            # Clean the data and get statistics
            cleaned_data, cleaning_stats = validate_and_clean_triplets(processed_data, rdpFiles)
            
            # Print cleaning statistics
            print(f"Original number of triplets: {cleaning_stats['original_triplets']}")
            print(f"Removed triplets: {cleaning_stats['removed_triplets']}")
            print(f"Remaining triplets: {cleaning_stats['remaining_triplets']}")

            # Save the processed data
            file_name = f"output_test/ml_input_{folder.name}.txt"
            file_path = Path(file_name)
            cleaned_data.drop(["ISeqs(A)"], axis = 1, inplace = True)  
            save_processed_data(cleaned_data, file_path)
        
        else:
            print("The requested file don't exist")
            print('Alig: ' + str(rdpFiles[0].exists()))
            pass

if __name__ == '__main__':
    argParser = argparse.ArgumentParser(description='Process RDP5 files')
    argParser.add_argument('-f', dest='folder', help='file path to parse')
    args = argParser.parse_args()

    global folder
    folder = Path(args.folder)

    # Change the path below to your target path.
    getFileNames(folderToParse=Path(folder))
    parsing_loop()

# folder = Path('output_test')
# rdpStatsFiles = [Path('dataRaw/UnseenTestSet/alignment_TestSet_3500-0.075-0.000118-200.faRecombIdentifyStats.csv')]
# rdpSimVReal = [Path('dataRaw/UnseenTestSet/alignment_TestSet_3500-0.075-0.000118-200.faSimVSRealCompare.csv')]
