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

    path = "../../simulation_data/indelible_arpip_output_filtered/"
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
    rel_name = 'indelible_312335.arpipnode_rel.txt'
    f = open(rel_name, 'r')
    lines = f.readlines()
    dict_node_relationship = extract_relationship_nodes(lines)

    oma_count = 0
    start_time = time.time()

    for file in os.listdir():
        if file.endswith('.arpipasr.fasta'):
            oma_name = file.split('_')[1].split('.')[0]
            # print(oma_name)
            oma_count += 1
            tr = Tree(f"indelible_{oma_name}.arpiptree.nwk", format=1)
            dict_data_brln[oma_name] = extract_branch_len(tr, dict_node_relationship)
            dict_data_age[oma_name] = extract_branch_age(tr, dict_node_relationship)

    print("The total number of %s files that has been processed." %oma_count)
    print("--- %s seconds --- " %(time.time() - start_time))
    df_brl = pd.DataFrame(dict_data_brln).transpose()
    df_brl.to_csv('../01_br_len_per_OMA_per_species.csv')
    print("The branch length file has been written successflly.")



