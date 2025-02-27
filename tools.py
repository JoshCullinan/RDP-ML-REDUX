# RDP Tensorflow Impl.
import pandas as pd
import numpy as np
from sklearn.metrics import classification_report
from collections import defaultdict

def class_report(y_true, y_preds):
    print(
        classification_report(y_true,
                              y_preds,
                              digits= 4,
                              target_names=["Parent", "Recombinant"]))
    
def flip(series):
    split = series.str.split(",", expand=True).values 

    return split.reshape(len(series) * 3)

def ingestor(recom_path):
    
    recom = pd.read_csv(recom_path, sep="\t")

    # def f(x): return x.strip("()")
    # for (colname, _) in recom.iteritems():
    #     recom.loc[:, colname] = recom.loc[:, colname].apply(f)
    
    allData = recom.apply(flip, axis=0) 

    return allData

def combine_three_rows(input_file, output_file):
    """
    Reads a CSV file and combines every three rows into one, renaming columns appropriately.
    Particularily used in the posistion selection NN to combine the original data into data that has all features from all viruses.
    
    Parameters:
    input_file (str): Path to the input CSV file or a Pandas DataFrame
    output_file (str): Path to save the output CSV file
    
    Returns:
    pd.DataFrame: The transformed DataFrame
    """
    # Read the CSV file
    if (type(input_file) == pd.DataFrame):
        df = input_file
    else:
        df = pd.read_csv(input_file)
    
    # Initialize lists to store the new rows
    new_rows = []
    
    # Process rows in groups of three
    for i in range(0, len(df), 3):
        if i + 2 < len(df):  # Make sure we have 3 rows to combine
            new_row = {}
            new_row['id'] = i // 3  # Add an ID column
            
            # Add recombinant values
            new_row['Recombinant1'] = df.iloc[i]['is_recombinant']
            new_row['Recombinant2'] = df.iloc[i+1]['is_recombinant']
            new_row['Recombinant3'] = df.iloc[i+2]['is_recombinant']
            
            # Process each column (except 'is_recombinant')
            for col in df.columns:
                if col != 'is_recombinant':
                    # Add values for each of the three rows with numbered suffix
                    new_row[f'{col}1'] = df.iloc[i][col]
                    new_row[f'{col}2'] = df.iloc[i+1][col]
                    new_row[f'{col}3'] = df.iloc[i+2][col]
            
            new_rows.append(new_row)
    
    # Create new DataFrame
    result_df = pd.DataFrame(new_rows)
    
    # Save to CSV
    result_df.to_csv(output_file, index=False)
    
    return result_df




def verify_triplet_positives(data):
    """
    Verify that each triplet has exactly one positive case.
    
    Args:
        data: numpy array of shape (n_triplets, 3, n_features) or DataFrame
        
    Returns:
        bool: True if valid, raises ValueError if invalid
    """
    if isinstance(data, pd.DataFrame):
        # Convert DataFrame to numpy array and reshape
        n_rows = len(data)
        if n_rows % 3 != 0:
            raise ValueError(f"Number of rows ({n_rows}) is not divisible by 3")
        values = data.values
        data = values.reshape(n_rows // 3, 3, -1)
    
    for triplet_idx, triplet in enumerate(data):
        positive_count = sum(row[-1] == 1 for row in triplet)
        if positive_count != 1:
            raise ValueError(
                f"Triplet {triplet_idx} has {positive_count} positive cases, expected exactly 1.\n"
                f"Triplet values:\n{triplet}"
            )
    return True

def balance_triplet_positions(input_file, output_file = None, random_seed=42):
    """
    Shuffle triplets in a dataset to ensure balanced positioning of positive examples,
    with thorough redistribution across positions.
    
    Args:
        input_file (str): Path to input CSV file
        output_file (str): Path to output CSV file
        random_seed (int, optional): Seed for random number generator
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    
    # Read the CSV file
    df = pd.read_csv(input_file)
    
    # Verify input data
    print("Verifying input data...")
    verify_triplet_positives(df)
    
    # Ensure the number of rows is divisible by 3
    if len(df) % 3 != 0:
        raise ValueError("Number of rows must be divisible by 3")
    
    # Convert dataframe to numpy array for easier manipulation
    data = df.values
    n_triplets = len(data) // 3
    
    # Reshape into triplets
    triplets = data.reshape(n_triplets, 3, -1)
    
    # First, identify the original position of positive examples in each triplet
    original_positions = []
    for triplet in triplets:
        pos_idx = np.where([row[-1] == 1 for row in triplet])[0][0]
        original_positions.append(pos_idx)
    
    # Count original distribution
    position_counts = defaultdict(int)
    for pos in original_positions:
        position_counts[pos] += 1
    
    # For each original position, redistribute a portion of its triplets
    shuffled_triplets = [None] * n_triplets  # Pre-allocate list with correct size
    processed_indices = set()
    
    for orig_pos in range(3):
        # Get triplets that originally had positive in this position
        pos_triplets = [(idx, trip) for idx, trip in enumerate(triplets) 
                       if original_positions[idx] == orig_pos and idx not in processed_indices]
        
        n_triplets_in_pos = len(pos_triplets)
        
        # Keep ~1/3 in original position, redistribute rest equally
        n_keep = n_triplets_in_pos // 3
        n_redistribute = n_triplets_in_pos - n_keep
        n_per_other_pos = n_redistribute // 2
        
        # Shuffle the triplets for random selection
        np.random.shuffle(pos_triplets)
        
        # Process triplets in three groups
        for i, (idx, triplet) in enumerate(pos_triplets):
            new_triplet = np.copy(triplet)
            
            if i < n_keep:
                # Keep in original position
                new_pos = orig_pos
            elif i < n_keep + n_per_other_pos:
                # Move to next position (cyclically)
                new_pos = (orig_pos + 1) % 3
            else:
                # Move to remaining position
                new_pos = (orig_pos + 2) % 3
            
            # Swap positive example to new position if needed
            if new_pos != orig_pos:
                new_triplet[[orig_pos, new_pos]] = new_triplet[[new_pos, orig_pos]]
            
            shuffled_triplets[idx] = (idx, new_triplet)
            processed_indices.add(idx)
    
    # Ensure all triplets were processed
    if None in shuffled_triplets:
        raise ValueError("Some triplets were not processed")
        
    # Extract just the triplets from the list
    final_triplets = [t[1] for t in shuffled_triplets]
    
    # Verify we have the correct number of triplets
    if len(final_triplets) != n_triplets:
        raise ValueError(f"Expected {n_triplets} triplets but got {len(final_triplets)}")
    
    # Combine shuffled triplets back into a single array
    shuffled_data = np.vstack(final_triplets)
    
    # Convert back to dataframe and save
    shuffled_df = pd.DataFrame(shuffled_data, columns=df.columns)
    
    # Print statistics
    print("\nOriginal distribution of positive examples:")
    for pos in range(3):
        count = position_counts[pos]
        print(f"Position {pos}: {count} ({count/n_triplets*100:.1f}%)")
    
    print("\nFinal distribution of positive examples:")
    final_positions = defaultdict(int)
    for triplet in final_triplets:
        pos_idx = np.where([row[-1] == 1 for row in triplet])[0][0]
        final_positions[pos_idx] += 1
    
    for pos in range(3):
        count = final_positions[pos]
        print(f"Position {pos}: {count} ({count/n_triplets*100:.1f}%)")
        
    # Calculate and print transition matrix
    print("\nTransition matrix (Original -> Final positions):")
    transition_matrix = np.zeros((3, 3))
    for orig_pos, triplet in zip(original_positions, final_triplets):
        final_pos = np.where([row[-1] == 1 for row in triplet])[0][0]
        transition_matrix[orig_pos][final_pos] += 1
    
    print("    Final Pos:")
    print("         0      1      2")
    for i in range(3):
        row_sum = sum(transition_matrix[i])
        percentages = [f"{(x/row_sum*100):6.1f}%" for x in transition_matrix[i]]
        print(f"Orig {i}: {' '.join(percentages)}")
    
    # Verify output data before saving
    print("\nVerifying output data...")
    verify_triplet_positives(shuffled_df)
    
    if output_file != None:
        # Save to new CSV file
        shuffled_df.to_csv(output_file, index=False)
        print(f"\nShuffled data saved to {output_file}")

    return shuffled_df