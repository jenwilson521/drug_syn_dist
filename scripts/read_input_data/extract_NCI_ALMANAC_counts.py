""""
02/27/2025 - initial search for NCI ALMANAC synergy distributions, target stats, etc

"""
def plot_syn_dist(cl_to_dbid_syn_dict):

    fig = plt.figure(dpi = 300)
    ax = plt.subplot()
    fig.tight_layout() 

    i_subplot = 1
    for cl in cl_to_dbid_syn_dict.keys():
        ax = plt.subplot(2, 3, i_subplot)
        syn_sum = [np.sum(syn_list) for syn_list in cl_to_dbid_syn_dict[cl].values()]
        syn_mean = [np.mean(syn_list) for syn_list in cl_to_dbid_syn_dict[cl].values()]
        syn_std = [np.std(syn_list) for syn_list in cl_to_dbid_syn_dict[cl].values()]

        l = 'mean = ' + str(np.round(np.mean(syn_std), 3)) + '\n' + 'std = ' + str(np.round(np.std(syn_std), 3))
        n_bins = np.round(len(syn_sum) / 30, 0)
        plt.hist(syn_std, bins = int(n_bins), alpha = 0.5, label = ' n = ' + str(len(syn_std)))
        plt.axvline(x = np.mean(syn_std), color = 'black', linestyle = 'dashed', label = l) 

        plt.title(cl, fontsize = 8)

        if i_subplot in [4, 5, 6]:
            plt.xlabel('ComboScore', fontsize = 6)
        
        plt.ylabel('# of combos', fontsize = 6)

        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
        for axis in ['top','right']:
            ax.spines[axis].set_linewidth(0)

        #plt.xticks(np.arange(-600, 800, 200), fontsize = 5)
        #plt.yticks(np.arange(0, 400, 100), fontsize = 5)
        plt.legend(fontsize = 5)  
        
        i_subplot += 1
    
    figure_name = 'NCI_ALMANAC_synergy_std_per_combo_distribtuion_6_breast.png'
    plt.savefig(os.path.join(figure_path, figure_name), bbox_inches = 'tight', pad_inches = 0.1)

def main():

    global data_path
    global results_path
    global figure_path
    data_path = '../../data'
    results_path = '../../results'
    figure_path = '../../results/figures'

    #Load NCI ALMANAC target/synergy dictionaries
    outpath = os.path.join(data_path, 'NCI_ALMANAC_inputs')
    dbid_to_targs_dict = pickle.load(open(os.path.join(outpath, 'dbid_to_targs_dict_db_022825.pkl'), 'rb')) #[DBID] = [list of targets]
    cl_to_dbid_syn_dict = pickle.load(open(os.path.join(outpath, 'cl_to_dbid_syn_dict_db_022825.pkl'), 'rb')) #[Cell line][DBID1__DBID2] = [list of synergy scores for dif doses]
    cl_to_dbid_targs_dict = pickle.load(open(os.path.join(outpath, 'cl_to_dbid_targs_dict_db_022825.pkl'), 'rb')) #[Cell line][DBID1__DBID2] = [[DBID1 target list], [DBID2 target list]]

    #Target frequency across DBIDs
    targ_to_dbid_dict = defaultdict(list)
    for dbid in dbid_to_targs_dict.keys():
        for targ in dbid_to_targs_dict[dbid]:
            targ_to_dbid_dict[targ].append(dbid)
    
    print('Number of unique targets:', len(targ_to_dbid_dict.keys()))
    #print(sorted(targ_to_dbid_dict.items(), key = lambda x: len(x[1]), reverse = True))
    three_count = 0
    two_count = 0
    one_count = 0
    for targ in targ_to_dbid_dict.keys():
        #print(targ)
        if len(targ_to_dbid_dict[targ]) == 3:
            three_count += 1
        if len(targ_to_dbid_dict[targ]) == 2:
            two_count += 1
        if len(targ_to_dbid_dict[targ]) == 1:
            one_count += 1
    print(list(targ_to_dbid_dict.keys()))
    print('# of targets with DBID frequency = 3:', three_count)
    print('# of targets with DBID frequency = 2:', two_count)
    print('# of targets with DBID frequency = 1:', one_count)

    #Load DREAM and GI target neighborhoods
    nbhd_file = 'PRN_nbhd_opt_alpha.pkl'
    exper_source = 'GI'
    GI_nbhds = pickle.load(open(os.path.join(results_path, 'neighborhoods', 'PathFX_0.5thr_int', exper_source, nbhd_file), 'rb'))
    exper_source = 'DREAM'
    DREAM_nbhds = pickle.load(open(os.path.join(results_path, 'neighborhoods', 'PathFX_0.5thr_int', exper_source, nbhd_file), 'rb'))
    
    targs_with_nbhds = set()
    for targ in targ_to_dbid_dict.keys():
        if targ in list(GI_nbhds.keys()):
            targs_with_nbhds.add(targ)
        if targ in list(DREAM_nbhds.keys()):
            targs_with_nbhds.add(targ)
    print('Targets with GI and/or DREAM neighborhoods:', targs_with_nbhds, len(targs_with_nbhds))

    targs_without_nbhds = set()
    for targ in targ_to_dbid_dict.keys():
        if targ not in list(GI_nbhds.keys()) and targ not in list(DREAM_nbhds.keys()):
            targs_without_nbhds.add(targ)
    print('Targets without GI and/or DREAM neighborhoods:', targs_without_nbhds, len(targs_without_nbhds))


    #Plot synergy distributions
    plot_syn_dist(cl_to_dbid_syn_dict)




if __name__ == "__main__":
    import os, pickle, functools
    import pandas as pd
    from collections import defaultdict
    import numpy as np
    import networkx as nx
    import matplotlib.pyplot as plt
    main()