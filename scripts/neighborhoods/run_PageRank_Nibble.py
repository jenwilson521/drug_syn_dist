def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def get_all_targs(targ_combos_list):
    """
        Converts all target combinations to single targets to be used as seeds in PRN

        :param targ_combos_list: list
            [T1__T2, T3__T4, ...]

        :returns all_targs: list
            [T1, T2, T3, T4, ...]
    """
    all_targs = set()
    for targ_combo in targ_combos_list:
        for targ in targ_combo.split('__'):
            for t in targ.split('_'): #for cases of multi-target drugs
                all_targs.add({t}) #MyGraph takes in a set of genes

    return list(all_targs)             

def run_prn(restart_proba, interactome_dir, all_targs):
    """
        Calls MyGraph, initializing mygraph class in order to run PageRank-Nibble 
        for all input targets using a given restart probability

        :param restart_proba: float
            Restart probability used for lazy random walk calculation of approx. PageRank vector
        
        :param interactome_dir: string
            Directory where interactome is stored
        
        :param all_targs: list
            List of all input targets from given experimental source
            [T1, T2, T3, ...]
        
        :retrun nbhd_dict: dictionary
            [Target] = [neighbor 1, neighbor 2, ...]
        
        :reutrn conduct_dict: dictionary
            [Target] = neighborhood conductance
    
    """
    my_graph = MyGraph(interactome_dir)

    nbhd_dict = defaultdict(list)
    conduct_dict = dict()
    for targ in all_targs:
        nbhd, cond = my_graph.local_community_detection_pagerank_nibble(targ, alpha = restart_proba, 
                                                                        epsilon = 1e-7, window = 3)

        nbhd_dict[list(targ)[0]] = list(nbhd)
        conduct_dict[list(targ)[0]] = cond
    
    return nbhd_dict, conduct_dict

def main():

    #Load interactome
    interactome_source = 'PathFX'
    threshold = 0.50
    exper_source = 'GI'
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    interactome_path = '../../data/interactomes'
    interactome_file = str(threshold) + 'thr_filtered_' + interactome_source + '_scored_interactome.txt'
    interactome_dir = os.path.join(interactome_path, interactome_file)
    interactome_df = pd.read_csv(interactome_dir, sep = '\t', na_filter = False)
    interactome = nx.from_pandas_edgelist(interactome_df, 'protein1', 'protein2', edge_attr = 'cost')


    #Load experimental data, extracting target combos
    if exper_source == 'GI':
        data_path = '../../data/GI_screen_inputs'
        file_name = 'GT_scores_dict.pkl'
        targ_to_syn_dict = pickle.load(open(os.path.join(data_path, file_name), 'rb')) #[T1__T2] = GT score
        targ_combos_list = list(targ_to_syn_dict.keys())
    elif exper_source == 'DREAM':
        data_path = '../../data/DREAM_inputs'
        file_name = 'targset_to_syn_dict.pkl'
        targ_to_syn_dict_per_cl = pickle.load(open(os.path.join(data_path, file_name), 'rb')) #[T1_T2_T3__TT1_TT2] = syn
        targ_combos_set = set()
        for cl in targ_to_syn_dict_per_cl.keys():
            for targset_combo in targ_to_syn_dict_per_cl[cl].keys():
                targ_combos_set.add(targset_combo)
            
        targ_combos_list = list(targ_combos_set)

    all_targs = get_all_targs(targ_combos_list)


    #Run PRN for all targets using different restart probabilities
    restart_proba_all = np.concatenate((np.arange(0.02, 0.15, 0.02), np.arange(0.05, 0.80, 0.05))) #22 values in the range 0.02-0.75
    dd_list = functools.partial(defaultdict, list)
    dd_float = functools.partial(defaultdict, float)
    nbhd_dict_all_alpha = defaultdict(dd_list)
    conduct_dict_all_alpha = defaultdict(dd_float)
    for restart_proba in restart_proba_all:
        rp = np.round(restart_proba)
        nbhd_dict, conduct_dict = run_prn(rp, interactome_dir, all_targs)
        nbhd_dict_all_alpha[rp] = nbhd_dict
        conduct_dict_all_alpha[rp] = conduct_dict

        print(rp, 'done')
    
    #Save results
    prn_results_path = '../../results/PRN_neighborhoods'
    prn_results_dir = os.path.join(prn_results_path, interactome_aname, exper_source)
    check_dir(prn_results_dir)

    file_name = 'PageRank_Nibble_nbhds_0.02-0.75_restart_proba.pkl'
    pickle.dump(nbhd_dict_all_alpha, open(os.path.join(prn_results_dir, file_name), 'wb'))

    file_name = 'PageRank_Nibble_conductances_0.02-0.75_restart_proba.pkl'
    pickle.dump(conduct_dict_all_alpha, open(os.path.join(prn_results_dir, file_name), 'wb'))


if __name__ == '__main__':
    import sys, pickle, os, functools
    import pandas as pd
    import networkx as nx
    import numpy as np
    from collections import defaultdict
    from os import path, curdir
    #sys.path.append(path.abspath('..') + "/src/")
    from myGraph import MyGraph
    main()