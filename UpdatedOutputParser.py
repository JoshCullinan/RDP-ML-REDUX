import pandas as pd
import numpy as np

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
            if pd.isna(recomb_stats.iloc[current_idx][' ISeqs(A)']):
                continue
                
            seq_ids = str(recomb_stats.iloc[current_idx][' ISeqs(A)']).split('$')
            print(seq_ids)
            seq_ids = [int(id_) for id_ in seq_ids if id_.strip()]
            
            # Check if the actual recombinant is in this sequence
            if actual_recomb in seq_ids:
                recomb_stats.at[current_idx, 'is_recombinant'] = 1
            else:
                recomb_stats.at[current_idx, 'is_recombinant'] = 0
    
    return recomb_stats



def save_processed_data(processed_df, output_path):
    """
    Save the processed DataFrame to a CSV file.
    
    Parameters:
    processed_df (pandas.DataFrame): The processed DataFrame to save
    output_path (str): Path where the CSV should be saved
    """
    processed_df.to_csv(output_path, index=False)

def validate_and_clean_triplets(processed_df):
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
            print(f"Removing triplet at rows {i}-{i+2}: {recomb_count} recombinants found")
    
    cleaned_df = processed_df.iloc[rows_to_keep].copy()
    stats['remaining_triplets'] = len(cleaned_df) // 3
    
    return cleaned_df, stats

if __name__ == '__main__':
    recomb_stats_path = 'dataRaw/UnseenTestSet/alignment_TestSet_2022-08-20-15_55_16_3300-0.07-0.000114-100.faRecombIdentifyStats.csv'
    sim_compare_path = 'dataRaw/UnseenTestSet/alignment_TestSet_2022-08-20-15_55_16_3300-0.07-0.000114-100.faSimVSRealCompare.csv'

    processed_data = process_recombination_data(recomb_stats_path, sim_compare_path)

    # Clean the data and get statistics
    cleaned_data, cleaning_stats = validate_and_clean_triplets(processed_data)

    # Print cleaning statistics
    print(f"Original number of triplets: {cleaning_stats['original_triplets']}")
    print(f"Removed triplets: {cleaning_stats['removed_triplets']}")
    print(f"Remaining triplets: {cleaning_stats['remaining_triplets']}")

    # Save the processed data
    output_path = 'processed_recombination_data.csv'
    save_processed_data(cleaned_data, output_path)