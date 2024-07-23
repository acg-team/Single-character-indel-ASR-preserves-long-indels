## Python 3
# Compute average gap length of taxa and root
# Created by:  Gholamhossein Jowkar <jowk@zhaw.ch>
# ACGT ZHAW
# Created date: Sep 2023
# Modified by:
# Modified date:

import os, re, sys
from Bio import AlignIO
from statistics import mean

import pandas as pd



def extract_gap_len(align, spc_list_name, dic_species={}):
    # making a dictionary of species
    for elmn in spc_list_name:
        dic_species[elmn] = []

    # indel/gap size dist per specie
    for record in align:
        tmp_rc_id = ''
        matches = list(re.finditer('-+', str(record.seq)))
        if len(matches) != 0:
            if tmp_rc_id != record.id:
                tmp_rc_id = record.id
            for region_number, match in enumerate(matches, 1):
                # store the indel length
                indel_length = match.end() - match.start()
                if record.id == 'Mus':
                    dic_species['Mus'].append(indel_length)
                elif record.id == 'Rattus':
                    dic_species['Rattus'].append(indel_length)
                elif record.id == 'Pan':
                    dic_species['Pan'].append(indel_length)
                elif record.id == 'Homo':
                    dic_species['Homo'].append(indel_length)
                elif record.id == 'Gorilla':
                    dic_species['Gorilla'].append(indel_length)
                elif record.id == 'Macaca':
                    dic_species['Macaca'].append(indel_length)
                elif record.id == 'root':
                    dic_species['root'].append(indel_length)
    return dic_species


def read_fasta_file(file):
    return AlignIO.read(file, "fasta")

if __name__=="__main__":
    os.chdir(os.path.dirname(sys.argv[0]))  # Change the working directory to the script directory
    path = "../../mammals_data/new_output/"
    os.chdir(path)
    print(f"Current working directory: {path}")

    list_node_names = ["Mus", "Rattus", "Pan", "Homo", "Gorilla", "Macaca", "root"]
    selected_taxa = ["Mus", "Rattus", "Pan", "Homo", "Gorilla", "Macaca"]

    global_dict_mean_gap_len = {}

    file_count = 0
    non_gap_taxa_count = 0
    non_gap_root_count = 0

    for file in os.listdir():
        if file.endswith('.arpipasr.fasta'):
            try:
                oma_name = file.split('_')[3].split('.')[0]

                # Read the alignment
                aln = read_fasta_file(file)

                file_count += 1
                dict_local_gap_len_pr_sp = {}
                dict_local_gap_len_pr_sp = extract_gap_len(aln, list_node_names)

                # For taxa
                lst_taxa_gap_len = []
                for i_key in selected_taxa:
                    if i_key in dict_local_gap_len_pr_sp:
                        lst_taxa_gap_len.extend(list(dict_local_gap_len_pr_sp[i_key]))

                if not lst_taxa_gap_len:
                    non_gap_taxa_count += 1
                    lst_taxa_gap_len = [0]

                lst_root_gap_len = list(dict_local_gap_len_pr_sp['root'])
                if not lst_root_gap_len:
                    non_gap_root_count += 1
                    lst_root_gap_len = [0]

                global_dict_mean_gap_len[oma_name] = {'taxa': mean(lst_taxa_gap_len), 'root': mean(lst_root_gap_len)}



            except Exception as e:
                print(f"An error occurred: {str(e)}")

    print("The total number of %s file were processed" %file_count)
    print("Total number of %s taxa without gap" %non_gap_taxa_count)
    print("Total number of %s root without gap" %non_gap_root_count)

    # load the stat data frame to select gappy samples

    data = pd.read_csv("../01_sub_summary_stat_global_df.csv")
    data.head()

    df_mean_gap_len = pd.DataFrame.from_dict(global_dict_mean_gap_len, orient='index')
    df_mean_gap_len.reset_index(inplace=True)
    df_mean_gap_len.rename(columns={'index': 'OMA_group'}, inplace=True)
    df_mean_gap_len['OMA_group'] = df_mean_gap_len['OMA_group'].astype(int)
    # Print the resulting DataFrame
    print(df_mean_gap_len.head())

    condition = data[data['is_gappy'] == 1]
    # condition['OMA_group'] = condition['OMA_group'].astype(int)
    # condition.head()
    print(len(condition))

    df_means = pd.merge(df_mean_gap_len, condition, on='OMA_group', how='inner')  # OMA_group
    df_means = df_means[['OMA_group', 'taxa', 'root']]
    df_means.head()
    print(len(df_means))
    df_means.to_csv('../02_sub_df_gap_length_taxa_root.csv')