## Python 3
# Extracting feature vector for the current dataset
# Created by:  Gholamhossein Jowkar <jowk@zhaw.ch>
# ACGT ZHAW
# Created date: November 2022
# Modified by:
# Modified date:
from collections import defaultdict

import numpy as np
import pandas as pd
import os, sys, time, re, operator, statistics, json

from ete3 import Tree

from Bio import SeqIO, AlignIO, Seq, SeqRecord


# path to the data
os.chdir(os.path.dirname(sys.argv[0]))

# function definition
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


def extract_pip_param(file):
    with open(file) as fp:
        for line in fp:
            s_list = [s.strip() for s in line.split(":")]
            if "Mu" in s_list:
                Mu = s_list[1]
            elif "Lambda" in s_list:
                Lam = s_list[1]
            elif "Loglikelihood" in s_list:
                Log = s_list[1]
    return Mu, Lam, Log


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


def gap_ratio(aln):
    """
    gap ratio = number of gaps/total number of AAs in the sequence
    """
    gap_cnt = gap_count(aln)
    return gap_cnt, gap_cnt / (aln.get_alignment_length() * len(aln))


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



def merge_dic_gap_len_pr_spc(dic_global_gep_len, dic_local_gap_len):
    """
    Merge two dictionaries
    """
    #     new _dic = {}
    overlapping_keys = dic_global_gep_len.keys()
    for key in overlapping_keys:
        dic_global_gep_len[key] = dic_global_gep_len[key] + dic_local_gap_len[key]
    # print("The local version", dic_local_gap_len)
    # print("The global version", dic_global_gep_len)
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
    # print("The local v", dic_local_gap_count_pr_spc)
    # print("The global", dic_global_gap_count_pr_spc)
    return dic_global_gap_count_pr_spc

def mode_count(lst_nums):
    lst_count = [[x, lst_nums.count(x)] for x in set(lst_nums)]
    # remove count <= 1
    lst_count = [x for x in set(lst_nums) if lst_nums.count(x) > 1]
    if len(lst_count) > 1:
        print([lst_count[0], lst_count[1]])
        x = lst_count[0]
    else:
        x = lst_count[0]
        print(x)
    return x

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


def summary_gap_length(dic_species):
    """
    :param dic_species: dictionary of name and list of gap lengths
    """
    dic_avg_spc = {}
    dic_min_spc = {}
    dic_max_spc = {}
    lst_len = []
    for spc_name, lst_spc_len in dic_species.items():
        lst_len.extend(lst_spc_len)
        if (len(lst_spc_len) != 0):
            dic_avg_spc[spc_name] = sum(lst_spc_len) / len(lst_spc_len)
            dic_max_spc[spc_name] = max(lst_spc_len)
            dic_min_spc[spc_name] = min(lst_spc_len)
        else:
            dic_avg_spc[spc_name] = 0
            dic_max_spc[spc_name] = 0
            dic_min_spc[spc_name] = 0
    if lst_len:
        # print("the list of lens", lst_len)
        avg_val = statistics.mean(lst_len)
        # mode_val = mode_count(lst_len)
        # print("custom", mode_val)
        mode_val = statistics.mode(lst_len)
        # print("stat package", mode_val)
        median_val = statistics.median(lst_len)
    else:
        print("Warning: the gap length list is empty")
        avg_val, mode_val, median_val = 0, 0, 0
    return avg_val, mode_val, median_val, dic_avg_spc, dic_min_spc, dic_max_spc


def tree_len(file):
    """
    Compute total tree length
    :param file: tree file in newick format
    """
    tr_len = 0
    tree = Tree(file, format=1)
    for n in tree.traverse(strategy="postorder"):
        tr_len += n.dist
    return tree, tr_len

if __name__=='__main__':
    count = 0
    ungap_count = 0

    start_time = time.time()
    path = "../../mammals_data/new_output/"
    os.chdir(path)

    header_asr = ["OMA_group", "aln_len", "total_gap_count", "msa_gap_count", "asr_gap_count", "min_gap_count"
        , "spc_wmin_gap_count", "max_gap_count", "spc_wmax_gap_count", "gap_ratio", "num_gappy_sites", "avg_gap_len"
        , "median_gap_len", "mode_gap_len", "min_gap_len", "spc_wmin_gap_len", "max_gap_len", "spc_wmax_gap_len"]
    header_indel = ["OMA_group", "num_ins", "num_del", "indel_ratio", "is_gappy"]
    header_tree = ["OMA_group", "tree_len"]
    header_param = ["OMA_group", "mu", "lambda", "loglikelihood", "intensity", "expected_len"]

    dic_local_gap_len_pr_sp = {"Mus": [], "Rattus": [], "Pan": [], "Homo": [], "Gorilla": [], "Macaca": [], "V2": [],
                               "V5": [], "V7": [], "V9": [], "root": []}
    dic_global_gap_len_pr_sp = dic_local_gap_len_pr_sp
    dic_global_gap_count_pr_spc = {"Mus": 0, "Rattus": 0, "Pan": 0, "Homo": 0, "Gorilla": 0, "Macaca": 0, "V2": 0,
                                   "V5": 0, "V7": 0, "V9": 0, "root": 0}
    dic_avg_gap_len = {}
    dic_count_ins_pr_sp = {"Mus": 0, "Rattus": 0, "Pan": 0, "Homo": 0, "Gorilla": 0, "Macaca": 0, "V2": 0, "V5": 0,
                           "V7": 0, "V9": 0, "root": 0}
    dic_count_del_pr_sp = {"Mus": 0, "Rattus": 0, "Pan": 0, "Homo": 0, "Gorilla": 0, "Macaca": 0, "V2": 0, "V5": 0,
                           "V7": 0, "V9": 0, "root": 0}
    # three data 1st:msa(multiple sequence alignment) 2nd:mas(multiple ancestral sequences)
    # 3rd:masa(multiple ancestral and sequence alignment)
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
                # print("file name", oma_name)
                lst_oma_names.append(oma_name)
                aln = read_fasta_file(file)
                aln_len = aln.get_alignment_length()

                msa, asr = aln_split(aln, 'V')
                msa_gap_count, msa_gap_ratio = gap_ratio(msa)
                asr_gap_count, asr_gap_ratio = gap_ratio(asr)
                aln_gap_count, aln_gap_ratio = gap_ratio(aln)
                # print("Alignment gap count", aln_gap_count)

                dic_local_gap_count_pr_spc = gap_count_pr_spc(aln)
                tpl_min_gap_count_pr_spc, tpl_max_gap_count_pr_spc = summary_gap_count(dic_local_gap_count_pr_spc)
                # print("Min gap count", tpl_min_gap_count_pr_spc)
                # print("Max gap count", tpl_max_gap_count_pr_spc)
                dic_global_gap_count_pr_spc = merge_dic_gap_count_pr_spc(dic_global_gap_count_pr_spc,
                                                                         dic_local_gap_count_pr_spc)

                dic_gap_colmns_lst[oma_name] = extract_gappy_cols(msa)
                num_gappy_sites = len(dic_gap_colmns_lst[oma_name])
                # print("These are column containg gap:",dic_gap_colmns_lst)

                dic_local_gap_len_pr_sp = extract_gap_len(aln, list_node_names)
                dic_global_gap_len_pr_sp = merge_dic_gap_len_pr_spc(dic_global_gap_len_pr_sp, dic_local_gap_len_pr_sp)
                # print("for the file %s is %s" %(oma_name, dic_local_gap_len_pr_sp))
                # print("the commulative is", dic_global_gap_len_pr_sp)

                avg_len_val, mode_len_val, median_len_val, dic_avg_gap_len, dic_min_gap_len, dic_max_gap_len \
                    = summary_gap_length(dic_local_gap_len_pr_sp)
                # print("dictionary of avg lengths", dic_avg_gap_len)
                # print("dictionary of min lengths", dic_min_gap_len)
                # print("dictionary of max lengths", dic_max_gap_len)
                # avg_spc_gap_len = statistics.mean(list(map(int, dic_local_gap_len_pr_sp.values())))
                # print("mean length:", avg_spc_gap_len)
                arg_max_specie_gap_len = max(dic_max_gap_len.items(), key=operator.itemgetter(1))[0]
                max_spc_gap_len = max(dic_max_gap_len.values())
                # print("max length:", max_spc_gap_len)
                # print("spc with max gap len", arg_max_specie_gap_len)
                arg_min_specie_gap_len = min(dic_min_gap_len.items(), key=operator.itemgetter(1))[0]
                min_spc_gap_len = min(dic_min_gap_len.values())
                # print("min length:", min_spc_gap_len)
                # print("scp with min gap len", arg_min_specie_gap_len, avg_len_val)

                local_data_lst_asr = [oma_name, aln_len, aln_gap_count, msa_gap_count, asr_gap_count
                    , tpl_min_gap_count_pr_spc[1], tpl_min_gap_count_pr_spc[0], tpl_max_gap_count_pr_spc[1]
                    , tpl_max_gap_count_pr_spc[0], aln_gap_ratio, num_gappy_sites, avg_len_val, median_len_val
                    , mode_len_val, min_spc_gap_len, arg_min_specie_gap_len, max_spc_gap_len, arg_max_specie_gap_len]
                global_data_lst_asr.append(local_data_lst_asr)
            except KeyError as er:
                print("Error: There is a problem with file [%s]: %s " % (oma_name, str(er)))

            #     elif file.endswith(".pipparam.txt"):
            if os.path.isfile(path + "/" + file.split('.')[0] + ".pipparam.txt"):
                try:
                    sub_file = path + "/" + file.split('.')[0] + ".pipparam.txt"
                    local_data_lst_param = []
                    # oma_name = file.split('/')[:-1].split('_')[3].split('.')[0]
                    mu, lam, log = extract_pip_param(sub_file)
                    # print("The value of%s ,%s and %s" %(mu, lam, log))
                    local_data_lst_param = [oma_name, mu, lam, log, round(float(lam) / float(mu), 7),
                                            round(float(lam) * float(mu), 7)]  # I=lam*mu, E=lam/mu
                    global_data_lst_param.append(local_data_lst_param)
                except KeyError as kr:
                    print("Error: Something is wrong with the  file [%s]" % (oma_name, str(kr)))
            #     elif file.endswith(".arpiptree.nwk"):
            if os.path.isfile(path + "/" + file.split('.')[0] + ".arpiptree.nwk"):
                try:
                    sub_file = path + "/" + file.split('.')[0] + ".arpiptree.nwk"
                    local_data_lst_tree = []
                    # oma_name = file.split('_')[3].split('.')[0]
                    tree, tr_len = tree_len(sub_file)
                    # print("Total tree length is", tr_len)
                    local_data_lst_tree = [oma_name, round(tr_len, 5)]
                    global_data_lst_tree.append(local_data_lst_tree)
                except KeyError as kr:
                    print("Error: Something is wrong with the file [%s]" % (oma_name, str(kr)))
            #     elif file.endswith(".arpipindel.txt"):
            if os.path.isfile(path + "/" + file.split('.')[0] + ".arpipindel.txt"):
                try:
                    sub_file = path + "/" + file.split('.')[0] + ".arpipindel.txt"
                    local_data_lst_indel = []
                    # oma_name = file.split('_')[3].split('.')[0]
                    #         print(oma_name)
                    #         events_count, events = read_file(path, file)
                    #         print(events)
                    #         print("This is the gappy column list", dic_gap_colmns_lst[oma_name])
                    # event_file_address = f"{path}/{file}"
                    lst_gappy_col_numbs = []
                    lst_gappy_col_numbs = dic_gap_colmns_lst[oma_name]
                    if lst_gappy_col_numbs:
                        is_gappy = 1
                        df_event, mat_ins, mat_del, dic_name_index = get_indel_events(sub_file, list_node_names,
                                                                                      lst_gappy_col_numbs)

                        num_ins = int(mat_ins.sum())
                        num_del = int(mat_del.sum())
                        # print("Total number of insertion is %s and deletion is %s" %(num_ins, num_del))

                        num_ins_pr_spc = np.sum(mat_ins, axis=0)
                        num_del_pr_spc = np.sum(mat_del, axis=0)
                        # print("Total number of insertion per spc", num_ins_pr_spc)
                        # print("Total number of deletion per spc", num_del_pr_spc)
                        global_num_ins = np.add(num_ins_pr_spc, global_num_ins)
                        global_num_del = np.add(num_del_pr_spc, global_num_del)

                        indel_ratio = num_del / num_ins
                        local_data_lst_indel = [oma_name, num_ins, num_del, indel_ratio, is_gappy]
                        global_data_lst_indel.append(local_data_lst_indel)
                    else:
                        is_gappy = 0
                        local_data_lst_indel = [oma_name, num_ins, num_del, indel_ratio, is_gappy]
                        global_data_lst_indel.append(local_data_lst_indel)
                        ungap_count += 1
                        # print(
                        #     "Warning: There was no gap in the file [%s] since the gappy column list is empty." % oma_name)
                except KeyError as kr:
                    print("Error: Something is wrong with the file [%s]" % (oma_name, str(kr)))
    #         display(df_event)
    #     if count == 500:
    #         break
    print("The total number of files that has been processed is %s" % count)
    print("The total number of files without gap is %s" % ungap_count)
    print("--- %s seconds --- " % (time.time() - start_time))

    # Building a dataFrames from each ARPIP file:
    asr_df = pd.DataFrame(global_data_lst_asr, columns=header_asr)
    print("Length of ars dataframe", len(asr_df))
    param_df = pd.DataFrame(global_data_lst_param, columns=header_param)
    indel_df = pd.DataFrame(global_data_lst_indel, columns=header_indel)
    print("Length of indel dataframe", len(indel_df))
    tree_df = pd.DataFrame(global_data_lst_tree, columns=header_tree)
    global_df = pd.merge(pd.merge(pd.merge(asr_df, indel_df, on='OMA_group'), param_df, on='OMA_group'), tree_df,
                         on='OMA_group')
    global_df.set_index('OMA_group', inplace=True)
    print(global_df)

    # store the data to file
    global_df.to_csv("../01_sub_summary_stat_global_df.csv")

    # Building a dataFrame for species view:
    lst_avg_gap_len_pr_spc = list(dic_avg_gap_len.values())
    lst_min_gap_len_pr_spc = list(dic_min_gap_len.values())
    lst_max_gap_len_pr_spc = list(dic_max_gap_len.values())

    lst_global_num_ins = global_num_ins.tolist()[0]
    lst_global_num_del = global_num_del.tolist()[0]
    lst_global_gap_count = list(dic_global_gap_count_pr_spc.values())

    header_pr_spc = ['number_insertion', 'number_deletion', 'number_gap', 'avg_gap_len'
        , 'min_gap_len', 'max_gap_len']

    d = {'number_insertion': lst_global_num_ins, 'number_deletion': lst_global_num_del
        , 'number_gap': lst_global_gap_count, 'avg_gap_len': lst_avg_gap_len_pr_spc
        , 'min_gap_len': lst_min_gap_len_pr_spc, 'max_gap_len': lst_max_gap_len_pr_spc}

    # global_spc_df
    spc_df = pd.DataFrame(d, index=list_node_names, columns=header_pr_spc)
    spc_df.reset_index()

    print(spc_df)
    # spc_df.to_csv('../01_sub_summary_stat_spc_df.csv')

    # Write down the gap length dictionary to a file
    with open("../01_sub_dic_global_gap_len_pr_sp.txt", 'w') as length_file:
        length_file.write(json.dumps(dic_global_gap_len_pr_sp))
    # print the dictionary
    selected_taxa = ["Mus", "Rattus", "Pan", "Homo", "Gorilla", "Macaca", "V2", "V5", "V7", "V9", "root"]
    dict_taxa_gap_len = {}
    dict_taxa_tmp = {}
    for i_key in selected_taxa:
        if i_key in dic_global_gap_len_pr_sp:
            dict_taxa_tmp[i_key] = dic_global_gap_len_pr_sp[i_key]
            print(i_key, len(dic_global_gap_len_pr_sp[i_key]))