## Python 3
# long indel situation
# Created by:  Gholamhossein Jowkar <jowk@zhaw.ch>
# ACGT ZHAW
# Created date: October 2023
# Modified by: jowk
# Modified date: Nov 2023

import os, sys, time, statistics, json, re, fnmatch
from statistics import mean

from Bio import AlignIO, SeqIO

from ete3 import Tree

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score, auc



def extract_relationship_nodes(lines):
    dict_node_relationship = {}
    for line in lines[1:-1]:
        s = line.split('\t')
        dict_node_relationship[s[0]] = s[1].strip()
    print("Relation file (kid:father)", dict_node_relationship)
    return dict_node_relationship

def extract_branch_len(tree, dict_node_relationship):
    list_node_names = []
    dict_br_len = {}
    for node in tree.traverse('postorder'):
        if node not in tree.iter_leaves():
            a = node.children[0]
            node.name = dict_node_relationship.get(a.name)
        dict_br_len[node.name] = node.dist
    return dict_br_len

def extract_branch_age(tree, dict_node_relationship):
    list_node_names = []
    dict_br_len = {}
    for node in tree.traverse('postorder'):
        if node not in tree.iter_leaves():
            a = node.children[0]
            node.name = dict_node_relationship.get(a.name)
        dict_br_len[node.name] = tr.get_distance(node.name)
    return dict_br_len



def extract_ungap_sequence_len(alignment):
    dict_seq_len={}
    for record in alignment:
        #print("%s without gaps is %s" % (record.id, record.seq.ungap("-")))
        #print("%s with length of %s" % (record.id, len(record.seq.ungap("-"))))
        dict_seq_len[record.id] = len(record.seq.ungap("-"))
    return dict_seq_len

def extract_gap_len_pr_spc(alignment, spc_list_name, dic_species={}):
    # making a dictionary of species
    for elmn in spc_list_name:
        dic_species[elmn] = []

    # indel/gap size dist per specie
    for record in alignment:
        tmp_rc_id = ''
        matches = list(re.finditer('-+', str(record.seq)))
        if len(matches) != 0:
            if tmp_rc_id != record.id:
                tmp_rc_id = record.id
            for region_number, match in enumerate(matches, 1):
                # store the indel length
                indel_length = match.end() - match.start()
                dic_species[record.id].append(indel_length)
    return dic_species


if __name__=="__main__":

    path = "../../simulation_data/indelibl_arpip_prank_output/"
    os.chdir(path)
    
    dict_val_index = {"Mus": [], "Rattus": [], "Pan": [], "Homo": [], "Gorilla": [], "Macaca": [], "V2": [],
                      "V5": [], "V7": [], "V9": [], "root": []}
    list_node_taxa = ["Mus", "Rattus", "Pan", "Homo", "Gorilla", "Macaca"]
    list_node_internal = ["V2", "V5", "V7", "V9", "root"]
    selected_taxa = ["Mus", "Rattus", "Pan", "Homo", "Gorilla", "Macaca"]

    dict_data_seq = {}
    dict_data_brln = {}
    dict_data_age = {}

    dict_data_gap_len = {}

    global_dict_taxa_gap_len = {}
    lst_taxa_gap_len = []

    global_dict_root_gap_len = {}
    lst_root_gap_len = []

    # extract relationship between nodes for tree traverse:
    rel_name = './312335.arpipnode_rel.txt'
    f = open(rel_name, 'r')
    lines = f.readlines()
    dict_node_relationship = extract_relationship_nodes(lines)

    oma_count = 0
    start_time = time.time()

    for file in os.listdir():
        if file.endswith('.arpipasr.fasta'):
            oma_name = file.split('.')[0]
            # print(oma_name)
            oma_count += 1

            # alignment information
            arpip_msa_aln = AlignIO.read(f"./{oma_name}.arpipmsa.fasta", "fasta")
            arpip_asr_aln = AlignIO.read(f"./{oma_name}.arpipasr.fasta", "fasta")

            # tree information
            tr = Tree(f"./{oma_name}.arpiptree.nwk", format=1)
            dict_data_brln[oma_name] = extract_branch_len(tr, dict_node_relationship)
            # print(dict_data_brln)

            # compute the ungappy sequence length:
            dict_seq = {}
            dict_msa = extract_ungap_sequence_len(arpip_msa_aln)
            dict_seq.update(dict_msa)
            dict_asr = extract_ungap_sequence_len(arpip_asr_aln)
            dict_seq.update(dict_asr)
            dict_data_seq[oma_name] = dict_seq
            # print(dict_data_seq)

            # compute gap length
            dict_gap = {}
            dict_gap_msa = extract_gap_len_pr_spc(arpip_msa_aln, list_node_taxa)
            dict_gap.update(dict_gap_msa)
            for i_key in selected_taxa:
                if i_key in dict_gap_msa:
                    lst_taxa_gap_len.extend(list(dict_gap_msa[i_key]))

            dict_gap_asr = extract_gap_len_pr_spc(arpip_asr_aln, list_node_internal)
            lst_root_gap_len.extend(list(dict_gap_asr['root']))
            dict_gap.update(dict_gap_asr)
            dict_data_gap_len[oma_name] = dict_gap
    #         print("The length of root elements is %s and the list is:"%(len(lst_root_gap_len)), lst_root_gap_len)
    #         print("The list of taxa gap lengths:", lst_taxa_gap_len)
    #         print("The dictionary of taxa gap lengths:", dict_data_gap_len)
    #         if oma_count > 4:
    #             break

    print("The total number of %s files that has been processed." % oma_count)
    print("--- %s seconds --- " % (time.time() - start_time))

    # making dataframes for all three dictionaries
    df_seq = pd.DataFrame(dict_data_seq).transpose()
    df_seq.to_csv('../04_seq_nogap_len_per_OMA_per_species_prank.csv')
    print("The sequence file has been written successflly.")

    df_brl = pd.DataFrame(dict_data_brln).transpose()
    df_brl.to_csv('../04_br_len_per_OMA_per_species.csv')
    print("The branch length file has been written successflly.")

    # Writing down the dictionary of gap lengths
    global_dict_taxa_gap_len['taxa'] = lst_taxa_gap_len
    global_dict_root_gap_len['root'] = lst_root_gap_len

    print("The length of taxa elements:", len(lst_taxa_gap_len))
    print("The length of root elements:", len(lst_root_gap_len))

    # Create DataFrames from the dictionaries
    df_taxa_gap_len = pd.DataFrame.from_dict({'Gap length': global_dict_taxa_gap_len['taxa']})
    df_root_gap_len = pd.DataFrame.from_dict({'Gap length': global_dict_root_gap_len['root']})
    df_taxa_gap_len.to_json('../04_dynamic_taxa_gap_length_prank.json', indent=2)
    df_root_gap_len.to_json('../04_dynamic_root_gap_length_prank.json', indent=2)

    # Write the dictionary to a JSON file
    with open('../04_gap_len_pr_oma_pr_spc_prank.json', 'w') as json_file:
        json.dump(dict_data_gap_len, json_file, indent=2)

    print("Json file has written sucessfully.")
