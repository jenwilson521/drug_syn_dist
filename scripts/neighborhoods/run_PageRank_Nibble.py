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
    targs_set = set()
    for targ_combo in targ_combos_list:
        for targ in targ_combo.split('__'):
            for t in targ.split('_'): #for cases of multi-target drugs
                targs_set.add(t)
    
    all_targs = list()
    for targ in targs_set:
        all_targs.append({targ})
         

    return all_targs    

def get_all_unique_NCI_targs(dbid_to_targs_dict, GI_nbhds, DREAM_nbhds):
    """
        Extracts all targets from the NCI ALMANAC dataset that don't overlap with targets from other screens

        :param dbid_to_targs_dict: dictionary
            [DrugBank ID] = [T1, T2, T3, ...]
        
        :param GI_nbhds: dictionary
            [Target] = [nbhd_protein_1, ...]

        :param DREAM_nbhds: dictionary
            [Target] = [nbhd_protein_1, ...]     
            
        :returns all_targs: list
            [T1, T2, T3, T4, ...]
    """ 

    targs_set = set()
    for dbid in dbid_to_targs_dict.keys():
        for targ in dbid_to_targs_dict[dbid]:
            if targ not in list(GI_nbhds.keys()) and targ not in list(DREAM_nbhds.keys()):
                targs_set.add(targ)
    
    all_targs = list()
    for targ in targs_set:
        all_targs.append({targ})

    return all_targs
            

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
        
        :return nbhd_dict: dictionary
            [Target] = [neighbor 1, neighbor 2, ...]
        
        :return conduct_dict: dictionary
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
    exper_source = 'NCI_ALMANAC'
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    interactome_path = '../../data'
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
        all_targs = get_all_targs(targ_combos_list)
    
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

    elif exper_source == 'NCI_ALMANAC':
        data_path = '../../data/NCI_ALMANAC_inputs'
        file_name = 'dbid_to_targs_dict_db_022825.pkl'
        dbid_to_targs_dict = pickle.load(open(os.path.join(data_path, file_name), 'rb')) #DBID = target list
         
        #Load GI and DREAM optimal neighborhoods, only applying PRN to UNIQUE targets from the NCI dataset
        nbhd_path = '../../results/neighborhoods'
        nbhd_file = 'PRN_nbhd_opt_alpha.pkl'
        nbhd_exper_source = 'GI'
        GI_nbhds = pickle.load(open(os.path.join(nbhd_path, interactome_aname, nbhd_exper_source, nbhd_file), 'rb'))
        nbhd_exper_source = 'DREAM'
        DREAM_nbhds = pickle.load(open(os.path.join(nbhd_path, interactome_aname, nbhd_exper_source, nbhd_file), 'rb'))

        all_targs = get_all_unique_NCI_targs(dbid_to_targs_dict, GI_nbhds, DREAM_nbhds)


    #Run PRN for all targets using different restart probabilities
    restart_proba_all = np.concatenate((np.arange(0.02, 0.15, 0.02), np.arange(0.05, 0.80, 0.05))) #22 values in the range 0.02-0.75
    dd_list = functools.partial(defaultdict, list)
    dd_float = functools.partial(defaultdict, float)
    nbhd_dict_all_alpha = defaultdict(dd_list)
    conduct_dict_all_alpha = defaultdict(dd_float)
    for restart_proba in restart_proba_all:
        rp = np.round(restart_proba, 3)
        nbhd_dict, conduct_dict = run_prn(rp, interactome_dir, all_targs)
        nbhd_dict_all_alpha[rp] = nbhd_dict
        conduct_dict_all_alpha[rp] = conduct_dict

        print(rp, 'restart probability done')
    
    #Save results
    prn_results_path = '../../results/neighborhoods'
    prn_results_dir = os.path.join(prn_results_path, interactome_aname, exper_source)
    check_dir(prn_results_dir)

    if exper_source == 'NCI_ALMANAC':
        file_name = 'Unique_PageRank_Nibble_nbhds_0.02-0.75_restart_proba.pkl'
    else:
        file_name = 'PageRank_Nibble_nbhds_0.02-0.75_restart_proba.pkl'
    pickle.dump(nbhd_dict_all_alpha, open(os.path.join(prn_results_dir, file_name), 'wb'))

    if exper_source == 'NCI_ALMANAC':
        file_name = 'Unique_PageRank_Nibble_conductances_0.02-0.75_restart_proba.pkl'
    else:
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