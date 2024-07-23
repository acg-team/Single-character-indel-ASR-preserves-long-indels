## Python 3
# Select subset of samples that have matching topologies
# Created by:  Gholamhossein Jowkar <xjok@zhaw.ch>
# ACGT ZHAW
# Created date: July 2024

import os, sys
import pandas as pd



if __name__ == '__main__':
    os.chdir(os.path.dirname(sys.argv[0]))  # Change the working directory to the script directory
    # Step 1: Load the IDs from the first text file into a DataFrame
    ids_df = pd.read_csv('../../mammals_data/00_agreeing_topology_IDs.txt', header=None, names=['id'])
    ids = ids_df['id'].tolist()  # Convert to a list if necessary

    # Step 2: Load the second DataFrame from a CSV file
    data_df = pd.read_csv('../../mammals_data/01_new_summary_stat_global_df.csv', index_col=0)
    # Display the first few rows of the DataFrame
    print(data_df.head())

    # Step 3: Filter the second DataFrame using the IDs from the first DataFrame
    filtered_data_df = data_df[data_df['OMA_group'].isin(ids)]

    # Display the filtered DataFrame
    print(filtered_data_df)
    # Step 4: Save the filtered DataFrame to a new CSV file
    # filtered_data_df.to_csv('../../mammals_data/02_filtered_summary_stat_global_df.csv')

    # Step 5: Number of samples with is gap flagged
    print("Number of samples with is_gap flag: {}".format(filtered_data_df['is_gappy'].sum()))
    # Step 6: Number of samples without is gap flagged
    print("Number of samples without is_gap flag: {}".format(filtered_data_df.shape[0] - filtered_data_df['is_gappy'].sum())    )


    # Step 7: select join table
    join_df = pd.read_csv('../../mammals_data/02_join_b_seq_len_per_OMA_per_species_gappy.csv')

    # Step 8: Filter the second DataFrame using the IDs from the first DataFrame
    filtered_data_df = join_df[join_df['OMA_group'].isin(ids)]

    # Display the filtered DataFrame
    print(filtered_data_df)
    # Step 9: Save the filtered DataFrame to a new CSV file
    filtered_data_df.to_csv('../../mammals_data/02_sub_join_b_seq_len_per_OMA_per_species_gappy.csv')
