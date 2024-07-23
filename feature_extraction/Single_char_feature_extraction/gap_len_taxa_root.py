## Python 3
# The empirical gap length distribution and compare average gap length of taxa and root
# Created by:  Gholamhossein Jowkar <jowk@zhaw.ch>
# ACGT ZHAW
# Created date: Sep 2023
# Modified by:
# Modified date:

import numpy as np
import pandas as pd
import os, sys, time, re, operator, statistics, json



if __name__=="__main__":
    os.chdir(os.path.dirname(sys.argv[0]))  # Change the working directory to the script directory

    # Reading the indel length file
    with open('../../mammals_data/01_sub_dic_global_gap_len_pr_sp.txt') as f:
        gap_len = f.read()

    dict_gap_len_js = json.loads(gap_len)

    selected_taxa = ["Mus", "Rattus", "Pan", "Homo", "Gorilla", "Macaca"]
    dict_taxa_gap_len = {}
    dict_taxa_tmp = {}

    # For taxa
    for i_key in selected_taxa:
        if i_key in dict_gap_len_js:
            dict_taxa_tmp[i_key] = dict_gap_len_js[i_key]
            print(i_key, len(dict_gap_len_js[i_key]))

    # lst_taxa_gap_values = list(dict_taxa_tmp.values())
    lst_taxa_gap_values = [value for values in dict_taxa_tmp.values() for value in values]

    lst_taxa_gap_values_filtered = [value for value in lst_taxa_gap_values if value <= 100]
    # lst_taxa_gap_values_filtered = lst_taxa_gap_values

    dict_taxa_gap_len["taxa"] = lst_taxa_gap_values_filtered
    print("The length of taxa dictionary", len(dict_taxa_gap_len["taxa"]))

    # For root
    dict_root_gap_len = {}
    lst_root_gap_values = list(dict_gap_len_js["root"])
    lst_root_gap_values_filtered = [value for value in lst_root_gap_values if value <= 100]
    # lst_root_gap_values_filtered = lst_root_gap_values

    dict_root_gap_len["root"] = lst_root_gap_values_filtered

    # print(dict_root_gap_len)
    # Create DataFrames from the dictionaries
    df_taxa = pd.DataFrame.from_dict({'Gap length': dict_taxa_gap_len['taxa']})
    df_root = pd.DataFrame.from_dict({'Gap length': dict_root_gap_len['root']})

    df_taxa.to_json('../../mammals_data/02_sub_json_taxa_gap_length.json', indent=2)
    df_root.to_json('../../mammals_data/02_sub_json_root_gap_length.json', indent=2)
    print("Json file wrote sucessfully.")