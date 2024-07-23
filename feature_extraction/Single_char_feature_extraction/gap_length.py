## Python 3
# Extracting gap length from the alignment file
# Created by:  Gholamhossein Jowkar <xjok@zhaw.ch>
# ACGT ZHAW
# Created date: November 2022

from collections import defaultdict
import numpy as np
import pandas as pd
import os, sys, time, re, operator, statistics, json

from ete3 import Tree

from Bio import SeqIO, AlignIO, Seq, SeqRecord


def aln_split(aln, asr_pref="V"):
    """
    Splite the alignment file into msa and maa using the prefix provided by user
    The default value is "V".
    :param aln: alignment file from fasta file converted to AlignIO (see BioPython for more).
    :param asr_pref: The prefix which the default value is "V".
    :return msa file and maa file seperately.
    """
    msa = []
    asr = []
    for rec in aln:
        try:
            if asr_pref in rec.id or "root" in rec.id:
                asr.append(rec)
            else:
                msa.append(rec)
        except KeyError as ker:
            print("ERROR: In the alignment file %s", str(ker))
    MSA = AlignIO.MultipleSeqAlignment(msa)
    ASR = AlignIO.MultipleSeqAlignment(asr)
    return MSA, ASR


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


def extract_gap_len(align, spc_list_name):
    dic_species = defaultdict(list)

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
                if record.id:
                    dic_species[record.id].append(indel_length)
    return dic_species


def gap_count(aln):
    """
    Count number of gap charchter denoted by '-' in the alignment file
    param: aln: AlignIO file (see BioPython for more)
    return: gap_count: The exact number of observed gap chars.
    """
    gap_count = 0
    for record in aln:
        for i in record:
            if i == '-':
                gap_count += 1
    return gap_count


def gap_count_pr_spc(aln):
    """
    return the dictinary of nodes with the number of gaps
    """
    dic_gap_count = {}
    for rec in aln:
        gap_count = 0
        for i in rec:
            if i == '-':
                gap_count += 1
        dic_gap_count[rec.id] = gap_count
    return dic_gap_count


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
        mat_deletion = np.zeros([seq_len, len(list_names)])
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



def merge_dic_gap_len_pr_spc(dic_global_gep_len, dic_local_gap_len):
    """
    Merge two dictionaries
    """
    overlapping_keys = dic_global_gep_len.keys()
    for key in overlapping_keys:
        dic_global_gep_len[key] = dic_global_gep_len[key] + dic_local_gap_len[key]
    return dic_global_gep_len


def merge_dic_gap_count_pr_spc(dic_global_gap_count_pr_spc, dic_local_gap_count_pr_spc):
    """
    sum two values with the same key
    :param dic_global_gap_count_pr_spc:
    :param dic_local_gap_count_pr_spc:
    :return:
    """
    overlapping_keys = dic_local_gap_count_pr_spc.keys()
    for key in overlapping_keys:
        dic_global_gap_count_pr_spc[key] = dic_global_gap_count_pr_spc[key]+dic_local_gap_count_pr_spc[key]
    return dic_global_gap_count_pr_spc



def read_fasta_file(file):
    # read the file
    return AlignIO.read(file, "fasta")


def read_file(path, file):
    """
    :param: path to the file
    :param: file: the name of the file
    """
    file_path = f"{path}/{file}"
    f = open(file_path, 'r')
    return len(f.readlines()) - 1, f.readlines()  # -1: since the last line is a blank line


def summary_gap_count(dic_spc):
    max_val_spc = max(dic_spc.values())
    argmax_spc = max(dic_spc, key=dic_spc.get)

    min_val_spc = min(dic_spc.values())
    argmin_spc = min(dic_spc, key=dic_spc.get)

    #     avg_spc = statistics.mean(dic_spc.values())
    return (argmin_spc, min_val_spc), (argmax_spc, max_val_spc)





def tree_len(file):
    """
    Compute total tree length
    :param file: tree file in newick format
    """
    tr_len = 0
    tree = Tree(file)
    for n in tree.traverse(strategy="postorder"):
        tr_len += n.dist
    return tree, tr_len


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # path = os.getcwd()
    path = 'D:\\paper\\2022.10.26.04_arpip_prank_phyml_codeml_final\\04_arpip_prank_phyml_codeml_final\\output\\'
    os.chdir(path)
    print(f"Current working directory: {path}")

    count = 0

    start_time = time.time()

    # header for the output files
    header_asr = ["OMA_group", "aln_len", "total_gap_count", "msa_gap_count", "asr_gap_count",  "gap_ratio", "avg_gap_len"]

    header_indel = ["OMA_group", "num_ins", "num_del", "indel_ratio", "is_gappy"]

    header_tree = ["OMA_group", "tree_len"]

    dic_local_gap_len_pr_sp = {} #{"Mus": [], "Rattus": [], "Pan": [], "Homo": [], "Gorilla": [], "Macaca": [], "V2": [],
                               # "V5": [], "V7": [], "V9": [], "root": []}

    dic_global_gap_len_pr_sp ={}#= dic_local_gap_len_pr_sp

    dic_global_gap_count_pr_spc = {"Mus": 0, "Rattus": 0, "Pan": 0, "Homo": 0, "Gorilla": 0, "Macaca": 0, "V2": 0,
                                   "V5": 0, "V7": 0, "V9": 0, "root": 0}
    dic_avg_gap_len = {}

    dic_count_ins_pr_sp = {"Mus": 0, "Rattus": 0, "Pan": 0, "Homo": 0, "Gorilla": 0, "Macaca": 0, "V2": 0, "V5": 0,
                           "V7": 0, "V9": 0, "root": 0}
    dic_count_del_pr_sp = {"Mus": 0, "Rattus": 0, "Pan": 0, "Homo": 0, "Gorilla": 0, "Macaca": 0, "V2": 0, "V5": 0,
                           "V7": 0, "V9": 0, "root": 0}

    list_node_names = ["Mus", "Rattus", "Pan", "Homo", "Gorilla", "Macaca", "V2", "V5", "V7", "V9", "root"]
    dic_gap_colmns_lst = {}

    lst_oma_names = []
    # defining lists and array of lists for storing the data
    global_data_lst_asr = []
    local_data_lst_asr = []
    global_data_lst_param = []
    local_data_lst_param = []
    global_data_lst_tree = []
    local_data_lst_tree = []
    global_data_lst_indel = []
    local_data_lst_indel = []

    global_num_ins = np.zeros([1, len(list_node_names)])
    global_num_del = np.zeros([1, len(list_node_names)])

    for file in os.listdir():
        if file.endswith(".arpipasr.fasta"):
            try:
                count += 1
                local_data_lst_asr = []
                oma_name = file.split('_')[3].split('.')[0]

                lst_oma_names.append(oma_name)
                aln = read_fasta_file(file)
                aln_len = aln.get_alignment_length()

                msa, asr = aln_split(aln, 'V')
                dic_local_gap_count_pr_spc = gap_count_pr_spc(aln)

                dic_global_gap_count_pr_spc = merge_dic_gap_count_pr_spc(dic_global_gap_count_pr_spc,
                                                                         dic_local_gap_count_pr_spc)

                dic_gap_colmns_lst[oma_name] = extract_gappy_cols(msa)
                num_gappy_sites = len(dic_gap_colmns_lst[oma_name])
                # print("These are column containg gap:",dic_gap_colmns_lst)

                dic_local_gap_len_pr_sp = extract_gap_len(aln, list_node_names)
                # store gap length per species with oma_name in nested dictionary
                dic_global_gap_len_pr_sp[oma_name] = dic_local_gap_len_pr_sp
                # dic_global_gap_len_pr_sp = merge_dic_gap_len_pr_spc(dic_global_gap_len_pr_sp, dic_local_gap_len_pr_sp)


                # local_data_lst_asr = [oma_name, aln_len,  num_gappy_sites ]
                # global_data_lst_asr.append(local_data_lst_asr)
            except KeyError as er:
                print("Error: There is a problem with file [%s]: %s " % (oma_name, str(er)))

    print("The total number of files that has been processed is %s" % count)
    print("--- %s seconds --- " % (time.time() - start_time))

        # store dic_global_gap_len_pr_sp in json file
    with open('gap_len.json', 'w') as f:
        json.dump(dic_global_gap_len_pr_sp, f, indent=4 )
