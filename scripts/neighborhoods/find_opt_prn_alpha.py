def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def rank_alpha(targ_to_tuple_dict):
    """
        Sorts all target neighborhoods (n = 22 alphas) according to size and conductance,
        retaining the alpha that yields the largest and most compact neighborhood per target

        :param targ_to_tuple_dict: dictionary
            [Target] = [(alpha_1, nbhd_size, conductance), (alpha_2, nbhd_size, conductance), ...]

        :return targ_to_opt_alpha_tuple_dict: dictionary
            [Target] = (optimal alpha, nbhd_size, conductance)
    """
    targ_to_opt_alpha_tuple_dict = dict()
    for targ in targ_to_tuple_dict.keys():
        #Sort alphas and neighborhood sizes according to size (bigger size = lower rank)
        alpha_by_size = list(zip(*sorted(targ_to_tuple_dict[targ], key = lambda x:x[1], reverse = True)))[0]
        size_by_size = list(zip(*sorted(targ_to_tuple_dict[targ], key = lambda x:x[1], reverse = True)))[1]
        inv_size_by_size = [1/size for size in size_by_size] #For rankdata function -> larger neighborhood = smaller rank

        #Sort alphas and neighborhood conductances according to conductance (lower conductance = lower rank)
        alpha_by_conduct = list(zip(*sorted(targ_to_tuple_dict[targ], key = lambda x:x[2], reverse = False)))[0]
        conduct_by_conduct = list(zip(*sorted(targ_to_tuple_dict[targ], key = lambda x:x[2], reverse = False)))[2]

        #Give sorted lists ranks (e.g., 1, 2, ...)
        size_rank = {alpha: rank for (alpha, rank) in zip(alpha_by_size, rankdata(inv_size_by_size, method = 'min'))}
        conduct_rank = {alpha: rank for (alpha, rank) in zip(alpha_by_conduct, rankdata(conduct_by_conduct, method = 'min'))}

        avg_rank = {alpha: np.mean([(0.5 * size_rank[alpha]), conduct_rank[alpha]]) for alpha in size_rank.keys()} #Prioritizing size ranking x2

        min_rank_alpha = min(avg_rank, key = avg_rank.get)

        targ_to_opt_alpha_tuple_dict[targ] = (min_rank_alpha, len(nbhd_dict_all_alpha[min_rank_alpha][targ]), conduct_dict_all_alpha[min_rank_alpha][targ])

        return targ_to_opt_alpha_tuple_dict
    
def make_hist(ax, all_list, list_label):
    """
        Plots histogram for various neighborhood metrics

        :param ax: object
            Axis object for all subplots
        
        :param all_list: list
            List of given metric (list_label) for all targets
        
        :param list_label: string
            Metric being plotted

    """
    
    n_bins = np.round(len(all_list) / 6, 0)
    plt.hist(all_list, bins = int(n_bins), color = 'r', alpha = 0.75, label = 'n = ' + str(len(all_list)))
    plt.axvline(x = np.median(all_list), color = 'r', linestyle = 'dashed', label = 'median = ' + str(np.round(np.median(all_list), 3))) 
    ax.set_xlabel(list_label)
    plt.legend(loc = 'upper right', fontsize = 14)

def main():

    #Global variables
    global nbhd_dict_all_alpha
    global conduct_dict_all_alpha

    #Load interactome
    interactome_source = 'PathFX'
    threshold = 0.50
    exper_source = 'GI'
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    interactome_path = '../../data/interactomes'
    interactome_file = str(threshold) + 'thr_filtered_' + interactome_source + '_scored_interactome.txt'
    interactome_df = pd.read_csv(os.path.join(interactome_path, interactome_file), sep = '\t', na_filter = False)
    interactome = nx.from_pandas_edgelist(interactome_df, 'protein1', 'protein2', edge_attr = 'cost')
    int_nodes = list(interactome.nodes())

    #Load neighborhood results
    prn_results_path = '../../results/PRN_neighborhoods'
    prn_results_dir = os.path.join(prn_results_path, interactome_aname, exper_source)
    check_dir(prn_results_dir)

    file_name = 'PageRank_Nibble_nbhds_0.02-0.75_restart_proba.pkl'
    nbhd_dict_all_alpha = pickle.load(open(os.path.join(prn_results_dir, file_name), 'rb'))

    file_name = 'PageRank_Nibble_conductances_0.02-0.75_restart_proba.pkl'
    conduct_dict_all_alpha = pickle.load(open(os.path.join(prn_results_dir, file_name), 'rb'))

    #Make dictionary that stores (alpha, neighborhood size, conductance) as tuple
    targ_to_tuple_dict = defaultdict(list)
    for alpha in nbhd_dict_all_alpha.keys():
        for targ in nbhd_dict_all_alpha[alpha].keys():
            if len(nbhd_dict_all_alpha[alpha][targ]) > 0:
                conduct = conduct_dict_all_alpha[alpha][targ]
                nbhd_size = nbhd_dict_all_alpha[alpha][targ]

                targ_to_tuple_dict[targ].append((alpha, nbhd_size, conduct))

    #Rank alphas according to neighborhood size and conductance, save optimal
    targ_to_opt_alpha_tuple_dict = rank_alpha(targ_to_tuple_dict)     

    #Extract alpha, size, and conductance lists across targets for histogram plots
    all_alpha = list(zip(*list(targ_to_opt_alpha_tuple_dict.values())))[0]
    all_size = list(zip(*list(targ_to_opt_alpha_tuple_dict.values())))[1]
    all_conduct = list(zip(*list(targ_to_opt_alpha_tuple_dict.values())))[2]

    #Define first neighbors neighborhoods, size and conductance
    firstn_nbhd_dict = defaultdict(list)
    firstn_conduct_dict = dict()
    firstn_size = list()
    firstn_conduct = list()
    for targ in targ_to_opt_alpha_tuple_dict.keys():
        firstn_nbhd = [fn for fn in interactome.neighbors(targ)]
        firstn_nbhd = firstn_nbhd + [targ]
        firstn_nbhd_dict[targ] = firstn_nbhd
        firstn_size.append(len(firstn_nbhd))

        conduct = nx.conductance(interactome, firstn_nbhd, T = None, weight = 'cost')
        firstn_conduct_dict[targ] = conduct
        firstn_conduct.append(conduct)
    
    #Plot stats for PRN w/ optimal alpha vs. first neighbors neighborhoods
    fig = plt.figure(figsize = (20, 7), dpi = 300)

    ax = plt.subplot(2, 3, 1)
    make_hist(ax, all_alpha, 'Optimal PRN alpha')

    ax = plt.subplot(2, 3, 2)
    make_hist(ax, all_size, 'Optimal PRN size')

    ax = plt.subplot(2, 3, 3)
    make_hist(ax, firstn_size, 'Firstn size')

    ax = plt.subplot(2, 3, 4)
    make_hist(ax, all_conduct, 'Optimal conductance')

    ax = plt.subplot(2, 3, 5)
    make_hist(ax, firstn_conduct, 'Firstn conductance')

    #Save figure
    figure_path = '../../figures/PRN_neighborhoods'
    check_dir(figure_path)
    figure_name = 'PRN_optimal_alpha_vs_firstn.png'
    plt.savefig(os.path.join(figure_path, figure_name))

    #Save optimal alpha and firstn neighborhood dictionaries
    opt_alpha_dict = dict()
    opt_nbhd_dict = dict()
    opt_conduct_dict = dict()
    for targ in targ_to_opt_alpha_tuple_dict:
        alpha = targ_to_opt_alpha_tuple_dict[targ][0]
        opt_alpha_dict[targ] = alpha
        opt_nbhd_dict[targ] = nbhd_dict_all_alpha[alpha][targ]
        opt_conduct_dict[targ] = conduct_dict_all_alpha[alpha][targ]

    file_name = 'PRN_alpha_opt_alpha.pkl'
    outpath = os.path.join(prn_results_dir, file_name)
    pickle.dump(opt_alpha_dict, open(outpath, 'wb'))

    file_name = 'PRN_nbhd_opt_alpha.pkl'
    outpath = os.path.join(prn_results_dir, file_name)
    pickle.dump(opt_nbhd_dict, open(outpath, 'wb'))

    file_name = 'PRN_conduct_opt_alpha.pkl'
    outpath = os.path.join(prn_results_dir, file_name)
    pickle.dump(opt_conduct_dict, open(outpath, 'wb'))

    file_name = 'firstn_nbhd.pkl'
    outpath = os.path.join(prn_results_dir, file_name)
    pickle.dump(firstn_nbhd_dict, open(outpath, 'wb'))

    file_name = 'firstn_conduct.pkl'
    outpath = os.path.join(prn_results_dir, file_name)
    pickle.dump(firstn_conduct_dict, open(outpath, 'wb'))

if __name__ == '__main__':
    import os, pickle, functools
    import pandas as pd
    import networkx as nx
    import numpy as np
    from collections import defaultdict
    import matplotlib.pyplot as plt
    from scipy.stats import rankdata
    main()