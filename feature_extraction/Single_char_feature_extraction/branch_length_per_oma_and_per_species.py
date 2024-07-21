## Python 3
# Branch length extraction per lineage:
# Created by:  Gholamhossein Jowkar <jowk@zhaw.ch>
# ACGT ZHAW
# Created date: March 2023
# Modified by:
# Modified date:

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import os, sys, time, re, operator, statistics

from ete3 import Tree

from Bio import SeqIO, AlignIO, Seq, SeqRecord

import json
from pathlib import Path


def tree_read(file):
    """
    Compute total tree length
    :param file: tree file in newick format
    """
    tr_len = 0
    tree = Tree(file)
    for n in tree.traverse(strategy="postorder"):
        tr_len += n.dist
    return tree, tr_len


def read_file(path, file):
    """
    reading the input file
    """
    file_path = f"{path}{file}"
    f = open(file_path, 'r')
    return f.readlines()


def extract_rel_nodes(lines):
    rel_dic = {}
    for line in lines[1:-1]:
        s = line.split('\t')
        rel_dic[s[0]] = s[1].strip()
    #     print("Relation file (kid:father)",rel_dic)
    return rel_dic

if __name__=="__main__":
    ###### main code #######
    path = "../../mammals_data/new_output/"
    os.chdir(path)

    start_time = time.time()
    rel_name = 'prank_codeml_OMAGroup_311903.arpipnode_rel.txt'
    lines = read_file(path, rel_name)
    dic_relation = extract_rel_nodes(lines)

    dic_node_branches_tmp = {"Mus": [], "Rattus": [], "Pan": [], "Homo": [], "Gorilla": [], "Macaca": [], "V2": [],
                             "V5": [], "V7": [], "V9": [], "root": []}


    mat_lineage_branches = np.zeros([3906, len(dic_node_branches_tmp)])

    # making a dictionary of name and index of tree's nodes
    dic_name_index = {}
    counter = 0
    for element in dic_node_branches_tmp.keys():
        dic_name_index[element] = counter
        counter += 1
    global_lst_br_ln = []

    count = 0
    os.chdir(path)
    for file in os.listdir():
        if file.endswith(".arpiptree.nwk"):
            asr_file = Path(f"{path}{file[:-14]}.arpipasr.fasta")
            if asr_file.is_file():
                oma_name = file.split('_')[3].split('.')[0]
                # print(oma_name)
                try:
                    count += 1
                    tree = Tree(file, format=1)
                    list_br = []
                    dic_node_branches = dic_node_branches_tmp
                    for node in tree.traverse("postorder"):
                        if node not in tree.iter_leaves():
                            a = node.children[0]
                            node.name = dic_relation.get(a.name)
                        dic_node_branches[node.name] = node.dist
                    # print(dic_node_branches.values())
                    list_br = list(dic_node_branches.values())
                    list_br.append(oma_name)
                    # print(list_br)
                    global_lst_br_ln.append(list_br)
                except KeyError:
                    print(file, count)

    print("The total number of files that has been processed is %s" % count)
    print("--- %s seconds --- " % (time.time() - start_time))

    # Building the dataFrame for each OMA file
    header_br_len = ["Mus", "Rattus", "Pan", "Homo", "Gorilla", "Macaca", "V2", "V5", "V7", "V9", "root", "OMA_group"]
    print(global_lst_br_ln[0])
    branch_length_df = pd.DataFrame(global_lst_br_ln, columns=header_br_len)
    branch_length_df = branch_length_df[
        ['OMA_group', 'Mus', 'Rattus', 'Pan', 'Homo', 'Gorilla', 'Macaca', 'V2', 'V5', 'V7', 'V9', 'root']]
    branch_length_df
    # store the data to file
    branch_length_df.to_csv("../01_sub_br_len_pr_oma_pr_spc_df.csv")