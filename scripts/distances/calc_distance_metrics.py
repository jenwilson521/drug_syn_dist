def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def calc_asp_dist(combo):

    """
        Calculates average shortest path (ASP) distance between two target neighborhoods
        Saves result to sp_dict

        :param combo: string
            target or target set pair
    """

    T1, T2 = combo.split('__')
    T1_nbhd = nbhd_dict[T1]
    T2_nbhd = nbhd_dict[T2]
    
    #Discarding hub protein
    if 'UBC' in T1_nbhd:
        T1_nbhd.remove('UBC')
    if 'UBC' in T2_nbhd:
        T2_nbhd.remove('UBC')
    
    spl = list()
    nopath = 0
    for n1 in T1_nbhd:
        for n2 in T2_nbhd:
            try:
                sp = nx.shortest_path(interactome, source = n1, target = n2)
                spl.append(len(sp) - 1)
            except:
                nopath += 1
    
    if len(spl) != 0:
        spl_max = np.max(spl)

        #Adds penalty (i.e., increases overall avg. distance) for all neighborhood pairs for which there is no shortest path
        nopath_list = [spl_max + 2] * nopath 
        sp = np.mean(spl + nopath_list)
        sp_dict[combo] =  sp
    else:
        sp_dict[combo] = float('nan')

def get_all_combos(targ_list):

    """
        Generates all pairwise combinations of targets

        :parm targ_list: list
            all input targets to be combined
        
        :pairwise_combos: list
            all pairwise combos
    
    """
    pairwise_combos = set()
    for t1 in targ_list:
        for t2 in targ_list:
            combo = '__'.join(sorted((t1, t2)))
            pairwise_combos.add(combo)
    
    return list(pairwise_combos)

def proba_t2t(source, target):
    """
        Calculates T1 -> T2 and T2 ->T1 probabilities, no neighborhood inclusion

        :param source: string
            Target whose diffusion profile is used, i.e., perturbation effects are from the perspective of source
        
        :parm target: string
            Target who is affected by source's perturbation
        
        :return cumm_proba / denom:
            Average probability of reaching target
        
    """
    source_dp = targ_dp_dict[source]
    targ_split = target.split('_') #Included to account for multi-target target sets

    cumm_proba = 0
    denom = 0
    for targ in targ_split:
        cumm_proba += source_dp[int_nodes_sorted.index(targ)]
        denom += 1

    return cumm_proba / denom

def avg_proba_t2n(source, target):
    """
        Calculates T1 -> N2 and T2 -> N1 probabilities, partial neighborhood inclusion

        :param source: string
            Target whose diffusion profile is used, i.e., perturbation effects are from the perspective of source
        
        :parm target: string
            Target whose neighborhood is affected by source's perturbation
        
        :return cumm_proba / denom:
            Average probability of reaching target neighborhood
    """

    source_dp= targ_dp_dict[source]
    target_nbhd = nbhd_dict[target]
    
    cumm_proba = 0
    denom = 0
    for n in target_nbhd:
        if n != 'UBC': #Discarding hub protein
            cumm_proba += source_dp[int_nodes_sorted.index(n)]
            denom += 1
    
    return cumm_proba / denom

def avg_proba_n2n(source, target):
    """
        Calculates N1 -> N2 and N2 -> N1 probabilities, complete neighborhood inclusion

        :param source: string
            Target whose neighborhood diffusion profile is used, i.e., perturbation 
            effects are from the perspective of source neighborhood
        
        :parm target: string
            Target whose neighborhood is affected by source neighborhood's perturbation
        
        :return cumm_proba / denom:
            Average probability of reaching target neighborhood
    """

    source_dp = nbhd_dict[source]
    target_nbhd = nbhd_dict[target]

    cumm_proba = 0
    denom = 0
    for n in target_nbhd:
        if n != 'UBC':
            cumm_proba += source_dp[int_nodes_sorted.index(n)]
            denom += 1
    
    return cumm_proba / denom

def calc_diffusion_dist(combo):

    """
    Calculates -log(avg probablity) of reaching one target (or target neighborhood) from the other and vice versa
    Saves distances from both directions in a sorted tuple so first value is from target with the larger neighborhood

    :param combo: string
        target or target set pair

    """

    #Determines which target has the larger neighborhood
    T1, T2 = combo.split('__')
    nbhd_size_dict = {T: len(nbhd_dict[T]) for T in [T1, T2]}
 
    T_ordered = sorted(nbhd_size_dict, key = nbhd_size_dict.get, reverse = True)
    if len(T_ordered) > 1:
        T1_ordered = T_ordered[0]
        T2_ordered = T_ordered[1]
    else: #Sham combos, e.g., T1__T1
        T1_ordered = T_ordered[0]
        T2_ordered = T_ordered[0]
    
    #Calculates probability of reaching T2 from T1 and vice versa
    proba_t12t2 = proba_t2t(T1_ordered, T2_ordered)
    proba_t22t1 = proba_t2t(T2_ordered, T1_ordered)

    #Calculates average probability of reaching T2's neighborhood from T1 and vice versa
    proba_t12n2 = avg_proba_t2n(T1_ordered, T2_ordered)
    proba_t22n1 = avg_proba_t2n(T2_ordered, T1_ordered)

    #Calculates average probability of reaching T2's neighborhood from T1's neighborhood and vice versa
    proba_n12n2 = avg_proba_n2n(T1_ordered, T2_ordered)
    proba_n22n1 = avg_proba_n2n(T2_ordered, T1_ordered)

    #Converting probabilities to distances and saving results
    t2t_dist_dict[combo] = (-np.log(proba_t12t2), -np.log(proba_t22t1))
    t2n_dist_dict[combo] = (-np.log(proba_t12n2), -np.log(proba_t22n1))
    n2n_dist_dict[combo] = (-np.log(proba_n12n2), -np.log(proba_n22n1))

def main():

    #Global variables
    global nbhd_dict
    global sp_dict
    global t2t_dist_dict
    global t2n_dist_dict
    global n2n_dist_dict

    global interactome
    global int_nodes_sorted

    global targ_dp_dict
    global nbhd_dp_dict


    #Load interactome
    interactome_source = 'PathFX'
    threshold = 0.50
    exper_source = 'DREAM'
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    interactome_path = '../../data/interactomes'
    interactome_file = str(threshold) + '_filtered_' + interactome_source + '_scored_interactome.txt'
    interactome_df = pd.read_csv(os.path.join(interactome_path, interactome_file), sep = '\t', na_filter = False)
    interactome = nx.from_pandas_edgelist(interactome_df, 'protein1', 'protein2', edge_attr = 'cost')
    int_nodes = list(interactome.nodes())
    int_nodes_sorted = sorted(int_nodes)
    
    global nbhd_dict
    prn_nbhd_path = '../../results/PRN_neighborhoods'
    prn_nbhd_dir = os.path.join(prn_nbhd_path, interactome_aname, exper_source)

    target_type = 'single' #multi
    if target_type == 'single': #all_targs = [T1, T2, T3, ...]
        #Load PageRank-Nibble neighborhoods
        nbhd_filename = 'PRN_nbhd_opt_alpha.pkl'

    elif target_type == 'multi': #all_targs = [T1_T2, TT1_TT2_TT3, ...]
        #Load merged PageRank-Nibble neighborhoods
        nbhd_filename = 'targset_to_merged_nbhd_dict.pkl'
    
    nbhd_dict = pickle.load(open(os.path.join(prn_nbhd_dir, nbhd_filename), 'rb'))
    all_targs = [targ for targ, nbhd in nbhd_dict.items() if len(nbhd) > 0] 
    pairwise_combos = get_all_combos(all_targs)

    #Diffusion profile path/directory
    dp_results_path = '../../results/RWR_diffusion_profiles'
    dp_results_dir = os.path.join(dp_results_path, exper_source, interactome_aname)


    manager = Manager()
    num_cores = mp.cpu_count()
    dd_float = functools.partial(defaultdict, float)

    alpha_list = [0.06, 0.15, 0.30, 0.45, 0.60, 0.75] #[0.06] #[0.04, 0.20, 0.50, 0.70]
    starting_list = ['targ', 'nbhd']
    t2t_dist_all_alpha = defaultdict(dd_float)
    t2n_dist_all_alpha = defaultdict(dd_float)
    n2n_dist_all_alpha = defaultdict(dd_float)

    for alphaD in alpha_list:
         
        #Load target- and target neighborhood-specific diffusion profiles
        for starting_nodes_choice in starting_list:
            dp_filename = target_type + '_target_RWR_diffusion_profiles_' + starting_nodes_choice + \
                        '_starting_nodes_alphaD_' + str(alphaD) + '.npz'
            
            dp_load_dict = np.load(os.path.join(dp_results_dir, dp_filename))
            if starting_nodes_choice == 'targ':
                #Something weird happens when opening dictionary, need to create a new one
                targ_dp_dict = dict()
                for targ in all_targs:
                    targ_dp_dict[targ] = dp_load_dict[targ]
            elif starting_nodes_choice == 'nbhd':
                #Something weird happens when opening dictionary, need to create a new one
                nbhd_dp_dict = dict()
                for targ in all_targs:
                    nbhd_dp_dict[targ] = dp_load_dict[targ]
             
        #Run multi-processing to compute all pairwise distances
        sp_dict = manager.dict()
        t2t_dist_dict = manager.dict() #Target to target distances
        t2n_dist_dict = manager.dict() #Target to neighborhood distances
        n2n_dist_dict = manager.dict() #Neighborhood to neighborhood distances

        with Pool(processes = num_cores) as pool:
            pool.map(calc_asp_dist, pairwise_combos)
            pool.map(calc_diffusion_dist, pairwise_combos)
        
        #Pool process does something weird with dictionaries, need to convert before using again
        t2t_dist_dict2 = dict(t2t_dist_dict)
        t2n_dist_dict2 = dict(t2n_dist_dict)
        n2n_dist_dict2 = dict(n2n_dist_dict)
        
        #Combines distances across all alphaDs
        for combo in pairwise_combos:
            t2t_dist_all_alpha[alphaD][combo] = t2t_dist_dict2[combo]
            t2n_dist_all_alpha[alphaD][combo] = t2n_dist_dict2[combo]
            n2n_dist_all_alpha[alphaD][combo] = n2n_dist_dict2[combo]
    
    #Save results
    results_path = '../../results/distance_dictionaries'
    results_dir = os.path.join(results_path, exper_source, interactome_aname)
    check_dir(results_dir)

    sp_filename = target_type + '_target_SP_dict.pkl'
    t2t_dist_filename = target_type + '_target_t2t_diffusion_dist_0.06-0.75_alphaD_dict.pkl'
    t2n_dist_filename = target_type + '_target_t2n_diffusion_dist_0.06-0.75_alphaD_dict.pkl'
    n2n_dist_filename = target_type + '_target_n2n_diffusion_dist_0.06-0.75_alphaD_dict.pkl'

    pickle.dump(dict(sp_dict), open(os.path.join(results_dir, sp_filename), 'wb'))
    pickle.dump(t2t_dist_all_alpha, open(os.path.join(results_dir, t2t_dist_filename), 'wb'))
    pickle.dump(t2n_dist_all_alpha, open(os.path.join(results_dir, t2n_dist_filename), 'wb'))
    pickle.dump(n2n_dist_all_alpha, open(os.path.join(results_dir, n2n_dist_filename), 'wb'))


if __name__ == "__main__":
    import pickle, os, functools
    import networkx as nx
    import pandas as pd
    import numpy as np
    from collections import defaultdict
    from multiprocessing import Pool, Manager
    import multiprocessing as mp
    main()