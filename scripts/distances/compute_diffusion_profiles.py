def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def build_transition_matrix():
    """
        Builds interactome's transition matrix (T = D^(-1)A)

        :param: None

        :return T: matrix
            Interactome transition matrix
    
    """
    #Adjacanecy matrix
    A = nx.adj_matrix(interactome, weight = 'cost') 
    A = A.todense()
    A = np.array(A, dtype = np.float64)

    #Diagonal matrix
    D = np.diag(np.sum(A, axis=0))
    d = np.sum(A, axis = 0) 
    d = 1/np.sqrt(np.sum(A, axis = 0))
    D_sqrt_inv = np.diag(d)

    #Transition matrix
    T = np.dot(np.linalg.inv(D), A) #Row normalized

    return T

def process_diffusion_profile(targ_str):

    """
        Runs RWR algorithm according to the starting point choice (target or neighborhood)
        Generalized to consider both single and multi-target RWRs

        :param targ_str: string
            if target_type = 'single', targ_str = T1
            if target_type = 'multi', targ_str = T1_T2_T3_...

    """
    targ_str
    targs = [t for t in targ_str.split('_')] #list of targets

    #Merges all target neighborhoods
    merged_nbhd = set()
    for t in targs:
        for n in prn_nbhd_dict[t]:
            if n != 'UBC':
                merged_nbhd.add(n)

    if len(merged_nbhd) != 0:
        print(targ_str, 'nbhd size:', len(merged_nbhd))

        if starting_nodes_choice == 'nbhd':
            p = rwr(list(merged_nbhd))
        else:
            p = rwr(targs)

        p_node_combined = zip(int_nodes, p)
        p_sorted = np.array([proba for _, proba in sorted(p_node_combined)])

    else:
        print(targ_str, 'No neighborhood detected')
        p_sorted = np.array(len(int_nodes))
    
    dp_dict[targ_str] = p_sorted

def rwr(starting_nodes):
    """
    Computes diffusion profiles using the random walk with restart (RWR) algorithm, starting from:
        All neighborhood proteins with equally probability if starting_nodes = nbhd
        The target if starting_nodes = targ
        All targets within a target set with equal probability if starting_nodes = targ and target_type = 'multi'

    :param starting_nodes: list
        All nodes where the RWR is initialized
    
    :return p: vector
        Steady-state probability distribution
 
    """
    
    s = np.zeros(len(int_nodes)) #starting position vector (constant for all iterations)
    p = np.zeros(len(int_nodes)) #probability distribution vector (updated at each iteration
    
    starting_prob = 1 / len(starting_nodes)

    #Distributing walker to all starting nodes with equal probability
    for n in starting_nodes:
        s[int_nodes.index(n)] = starting_prob
    
    p_t_1 = s #initializing p at t = t - 1 (p(t = 0) = s)
    p = (alphaD * s) + (1 - alphaD) * np.dot(p_t_1, T) #RWR standard equation

    n_steps = 1
    #Continue random walk until steady state is reached (defined by L1 norm <= tol)
    while np.linalg.norm((p - p_t_1), ord = 1) >= tol:
        p_t_1 = p
        p = (alphaD * s) + (1 - alphaD) * np.dot(p_t_1, T)
        n_steps += 1

    print('STEPS', n_steps)
    return p

def main():

    #Global variables
    global interactome
    global int_nodes
    global alphaD
    global tol
    global T
    global starting_nodes_choice
    global dp_dict

    #Load interactome
    interactome_source = 'PathFX'
    threshold = 0.50
    exper_source = 'GI'
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    interactome_path = '../../data/interactomes'
    interactome_file = str(threshold) + '_filtered_' + interactome_source + '_scored_interactome.txt'
    interactome_df = pd.read_csv(os.path.join(interactome_path, interactome_file), sep = '\t', na_filter = False)
    interactome = nx.from_pandas_edgelist(interactome_df, 'protein1', 'protein2', edge_attr = 'cost')
    int_nodes = list(interactome.nodes())

    #Load PageRank-Nibble neighborhoods
    prn_nbhd_path = '../../results/PRN_neighborhoods'
    prn_nbhd_dir = os.path.join(prn_nbhd_path, interactome_aname, exper_source)
    prn_filename = 'PRN_nbhd_opt_alpha.pkl'
    global prn_nbhd_dict
    prn_nbhd_dict = pickle.load(open(os.path.join(prn_nbhd_dir, prn_filename), 'rb'))

    target_type = 'single' #multi

    tol = 1e-7 #Tolerance used to determine RWR convergence
    T = build_transition_matrix()
    manager = Manager()
    num_cores = mp.cpu_count()

    starting_list = ['targ', 'nbhd']
    alpha_list = [0.06, 0.15, 0.30, 0.45, 0.60, 0.75] 
    for starting_nodes_choice in starting_list:
        for alphaD in alpha_list:

            #Load target starting points for RWR
            if target_type == 'single':
                targs = list(prn_nbhd_dict.keys())
            elif target_type == 'multi': #DREAM challenge multi-target sets
                DREAM_inputs_dir = '../../data/DREAM_inputs'
                targset2syn_file = 'targset_to_syn_dict.pkl'
                targset_to_syn_dict = pickle.load(open(os.path.join(DREAM_inputs_dir, targset2syn_file), 'rb')) #[cell line][Targetset1__Targetset2] = synergy
                multi_targ_set = set()
                for cl in targset_to_syn_dict.keys():
                    for combo in targset_to_syn_dict[cl].keys():
                        for targset in combo.split('__'):
                            multi_targ_set.add(targset) #{T1_T2_T3, ...}
                
                targs = list(multi_targ_set)

            #Multi-processing for full vector computation
            dp_dict = manager.dict()
            with Pool(processes = num_cores) as pool:
                pool.map(process_diffusion_profile, list(targs)) #Sorted by interactome node for consistency when computing distances

                pool.close() 
                pool.join()
            
            #Saving results
            #Something weird happens with dp_dict during pool process, need to convert back to dictionary
            dp_dict_save = dict(dp_dict)
            dp_results_path = '../../results/RWR_diffusion_profiles'
            dp_results_dir = os.path.join(dp_results_path, exper_source, interactome_aname)
            check_dir(dp_results_dir)

            filename = target_type + '_target_RWR_diffusion_profiles_' + starting_nodes_choice + \
                        '_starting_nodes_alphaD_' + str(alphaD) + '.npz'
            np.savez(os.path.join(dp_results_dir, filename), **dp_dict_save)
            print(starting_nodes_choice, alphaD, 'DONE')


if __name__ == "__main__":
    import pickle, os, functools
    import networkx as nx
    import pandas as pd
    import numpy as np
    from collections import defaultdict
    from multiprocessing import Pool, Manager
    import multiprocessing as mp

    main()