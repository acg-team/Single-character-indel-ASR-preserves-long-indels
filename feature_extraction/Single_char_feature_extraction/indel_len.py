## Python 3
# Indel length for mammalian data
# Created by:  Gholamhossein Jowkar <jowk@zhaw.ch>
# ACGT ZHAW
# Created date: Feb 2024
# Modified by:
# Modified date:


import numpy as np
import pandas as pd
import os, sys, time, re, operator, statistics, math, json
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter
from Bio import SeqIO, AlignIO, Seq, SeqRecord



############## functions ##################
def get_indel_events(event_file_address, list_names, lst_gappy_col_numbs=0):
    """
    param: seq_len: length of the sequence
    param: events: event file in text fomat
    param: list_names: name of all nodes
    return: puting the event in 3 files: mat_in/del matrix of events seperately
    mat_insrtion/deletion of size(m*n) and events of size(2*m)
    """

    # read the events file
    try:
        f = open(event_file_address, 'r')
        events = f.readlines()
        seq_len = len(events)
    #         print("seq_len", seq_len)
        if seq_len == 0:
            print("The input file is empty")

        dic_event = {}
        mat_insertion = np.zeros([seq_len, len(list_names)])
        mat_insertion
        mat_deletion = np.zeros([seq_len, len(list_names)])
        mat_deletion
        dic_name_index = {}
        counter = 0

        # making a dictionary of name and index of tree's nodes
        for element in list_names:
            dic_name_index[element] = counter
            counter += 1
        # check  to see if there is a column with gap in the alignment file
        if lst_gappy_col_numbs:
            for i, line in enumerate(events):
                if i in lst_gappy_col_numbs:
                    try:
                        if line.strip():
                            dic_event[i + 1] = line.strip()
                            s = line.split(';')
                            for sub in s:
                                if ':X' in sub:
                                    strd = sub.split(':X')
                                    for sub_str_id in strd[:-1]:
                                        mat_deletion[i, dic_name_index[sub_str_id]] = 1
                                elif ':I' in sub:
                                    stri = sub.split(':I')
                                    for sub_str_id in stri[:-1]:
                                        mat_insertion[i, dic_name_index[sub_str_id]] = 1
                    except KeyError:
                        if not (sub_str_id in list_names):
                            print("Error: Node with name [%s]  could not be found in the tree." % sub_str_id)
                        else:
                            print("ERROR")

        #     print("event dictionary", dic_event)
        df_event = pd.Series(dic_event)
        f.close()
        return df_event, mat_insertion, mat_deletion, dic_name_index
    except KeyError:
        print("Error: We cannot read the input file")


def extract_gappy_cols(msa):
    """
    :param msa: alignment file in AlionIO format (see BioPython for more)
    return: list of column numbers containing gap charcter
    """
    gappy_col = []
    # row counter
    r_counter = 0
    for record in msa:
        # col counter
        c_counter = 0
        for i in record:
            if i == '-':
                gappy_col.append(c_counter)
            c_counter += 1
        r_counter += 1
    return list(set(gappy_col))

def read_fasta_file(file):
    # read the file
    return AlignIO.read(file, "fasta")

# Function to count consecutive 1s in a list
def count_consecutive_ones(column):
    count = 0
    consecutive_counts = []

    for value in column:
        if value == 1:
            count += 1
        else:
            if count > 0:
                consecutive_counts.append(count)
                count = 0

    if count > 0:
        consecutive_counts.append(count)

    return consecutive_counts

if __name__=="__main__":
    ############## main code ##################
    path = "../../mammals_data/new_output/"
    os.chdir(path)

    start_time = time.time()
    count = 0

    dict_gap_colmns_lst = {}

    list_insertion = []
    list_deletion = []

    list_node_names = ["Mus", "Rattus", "Pan", "Homo", "Gorilla", "Macaca", "V2", "V5", "V7", "V9", "root"]
    lst_human_specific_indel_names = []

    for file in os.listdir():
        if file.endswith('.arpipasr.fasta'):
            oma_name = file.split('_')[3].split('.')[0]
            aln = read_fasta_file(file)
            aln_len = aln.get_alignment_length()
            #         print(oma_name)
            count += 1
            dict_gap_colmns_lst[oma_name] = extract_gappy_cols(aln)
            #         print("Length of gap elements are:", len(dict_gap_colmns_lst[oma_name] ))
            if dict_gap_colmns_lst[oma_name]:
                if os.path.isfile(path + "/" + file.split('.')[0] + ".arpipindel.txt"):
                    sub_file = path + "/" + file.split('.')[0] + ".arpipindel.txt"
                    local_data_lst_indel = []
                    lst_gappy_col_numbs = []
                    lst_gappy_col_numbs = dict_gap_colmns_lst[oma_name]
                    if lst_gappy_col_numbs:
                        is_gappy = 1
                        df_event, mat_ins, mat_del, dic_name_index = get_indel_events(sub_file, list_node_names,
                                                                                      lst_gappy_col_numbs)
                        if mat_del[:, :-1].any():
                            # Make dataframe from the deletion points:
                            df_del = pd.DataFrame(mat_del[:, :-1], columns=list_node_names[:-1])
                            #                     print(df_del)
                            consecutive_counts_per_column_dels = df_del.apply(count_consecutive_ones)
                            # Display the resulting dataframe
                            for col, col_counts in consecutive_counts_per_column_dels.items():
                                #                         print(f'{col}: {col_counts}')
                                if col_counts:  # Check if counts are not empty
                                    list_insertion.append(list(col_counts))
                        if mat_ins[:, :-1].any():
                            df_ins = pd.DataFrame(mat_ins[:, :-1], columns=list_node_names[:-1])
                            #                     print(df_ins)
                            # Apply the function to each column
                            consecutive_counts_per_column_ins = df_ins.apply(count_consecutive_ones)

                            # Display the resulting dataframe
                            for col, col_counts in consecutive_counts_per_column_ins.items():
                                #                         print(f'{col}: {col_counts}')
                                if col_counts:  # Check if counts are not empty
                                    list_deletion.append(list(col_counts))

    #             print('*******************', count)
    #         if count==12: break

    print("The total number of processed is %s" %count)
    print("--- %s seconds --- " %(time.time() - start_time))

    # Specify the output file
    output_file_1 = '../03_sub_deletion_length_list.json'

    # Write the list to JSON file
    with open(output_file_1, 'w') as json_file:
        json.dump(list_deletion, json_file)

    output_file_2 = '../03_sub_insertion_length_list.json'

    # Write the list to JSON file
    with open(output_file_2, 'w') as json_file:
        json.dump(list_insertion, json_file)