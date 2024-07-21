## Python 3
#  Branch length as divergence proxy
#  Distribution of branches per scpecie:
# Created by:  Gholamhossein Jowkar <jowk@zhaw.ch>
# ACGT ZHAW
# Created date: March 2023
# Modified by:
# Modified date:

import json
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import os, sys, time, re, operator, statistics

from ete3 import Tree

from Bio import SeqIO, AlignIO, Seq, SeqRecord


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
    print("Relation file (kid:father)", rel_dic)
    return rel_dic

if __name__=="__main__":
    ###### main code #######
    path = "../../mammals_data/new_output/"
    os.chdir(path)

    start_time = time.time()
    rel_name = 'prank_codeml_OMAGroup_311903.arpipnode_rel.txt'
    lines = read_file(path, rel_name)
    dic_relation = extract_rel_nodes(lines)

    dic_node_branches = {"Mus": [], "Rattus": [], "Pan": [], "Homo": [], "Gorilla": [], "Macaca": [], "V2": [],
                         "V5": [], "V7": [], "V9": [], "root": []}
    count = 0

    os.chdir(path)
    for file in os.listdir():
        if file.endswith(".arpiptree.nwk"):
            asr_file = Path(f"{path}{file[:-14]}.arpipasr.fasta")
            if asr_file.is_file():
                try:
                    count += 1
                    tree = Tree(file, format=1)
                    #             print(count)
                    list_names = []
                    for node in tree.traverse("postorder"):
                        if node not in tree.iter_leaves():
                            a = node.children[0]
                            node.name = dic_relation.get(a.name)
                        list_names.append(node.name)
                        dic_node_branches[node.name].append(node.dist)
                #         print(dic_node_branches)
                except KeyError:
                    print(file, count)

    print("The total number of files that has been processed is %s" % count)
    print("--- %s seconds --- " % (time.time() - start_time))

    # Write down the gap length dictionary to a file
    with open("../01_sub_dic_tree_branches_len_pr_specie.txt", 'w') as br_file:
        br_file.write(json.dumps(dic_node_branches))
    print("done")
    # print avergate per specie
    for key, value in dic_node_branches.items():
        print(key, statistics.mean(value))