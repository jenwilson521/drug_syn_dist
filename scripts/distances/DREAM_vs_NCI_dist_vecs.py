"""
03/05/2025 - comparing distance vector length distributions for DREAM vs. NCI ALAMANAC combos

"""
def get_pw_combos_per_cl(cl_to_combo_to_targset_dict, cl_to_combo_to_syn_dict, screen_name, syn_analysis):
    
    cl_to_dist_vec_dict = defaultdict(list)
    for cl in cl_to_combo_to_targset_dict.keys():
        for combo in cl_to_combo_to_targset_dict[cl].keys():
            pw_combo_set = set()
            if screen_name == 'NCI_ALMANAC':
                syn = np.sum(cl_to_combo_to_syn_dict[cl][combo])
            if screen_name == 'DREAM':
                syn = cl_to_combo_to_syn_dict[cl][combo]
                
            if syn_analysis == 'posyn' and syn > 0:
                ts1_list = cl_to_combo_to_targset_dict[cl][combo][0]
                ts2_list = cl_to_combo_to_targset_dict[cl][combo][1]

                for t1 in ts1_list:
                    for t2 in ts2_list:
                        pw_combo_set.add('__'.join(sorted((t1, t2))))
            
                cl_to_dist_vec_dict[cl].append(list(pw_combo_set))
                #cl_to_dist_vec_dict[cl].append((combo, list(pw_combo_set)))
            elif syn_analysis == 'negsyn' and syn <= 0:
                ts1_list = cl_to_combo_to_targset_dict[cl][combo][0]
                ts2_list = cl_to_combo_to_targset_dict[cl][combo][1]

                for t1 in ts1_list:
                    for t2 in ts2_list:
                        pw_combo_set.add('__'.join(sorted((t1, t2))))
            
                cl_to_dist_vec_dict[cl].append(list(pw_combo_set))
                #cl_to_dist_vec_dict[cl].append((combo, list(pw_combo_set)))
            elif syn_analysis == 'allsyn':
                ts1_list = cl_to_combo_to_targset_dict[cl][combo][0]
                ts2_list = cl_to_combo_to_targset_dict[cl][combo][1]

                for t1 in ts1_list:
                    for t2 in ts2_list:
                        pw_combo_set.add('__'.join(sorted((t1, t2))))
            
                cl_to_dist_vec_dict[cl].append(list(pw_combo_set))
                #cl_to_dist_vec_dict[cl].append((combo, list(pw_combo_set)))
    
    return cl_to_dist_vec_dict

def plot_dist_vec_len_distribution(cl_to_dist_vec_dict, screen_name, syn_analysis):

    if screen_name == 'DREAM':
        screen_color = 'dodgerblue'
    else:
        screen_color = 'orangered'

    fig = plt.figure(dpi = 300)
    ax = plt.subplot()
    fig.tight_layout() 

    i_subplot = 1
    for cl in cl_to_dist_vec_dict.keys():
        ax = plt.subplot(2, 3, i_subplot)

        dist_vec_len_list = [len(dist_vec) for dist_vec in cl_to_dist_vec_dict[cl]]
        print(max(dist_vec_len_list))

        dist_vec_len_mean = np.mean(dist_vec_len_list)
        dist_vec_len_median = np.median(dist_vec_len_list)

        l_mean = 'mean = ' + str(np.round(dist_vec_len_mean, 3)) + '\n' + 'std = ' + str(np.round(np.std(dist_vec_len_list), 3))
        l_median = 'median = ' + str(np.round(dist_vec_len_median, 3))

        #n_bins = np.round(len(dist_vec_len_list) / 30, 0)
        plt.hist(dist_vec_len_list, bins = 10, color = screen_color, alpha = 0.5, label = ' n = ' + str(len(dist_vec_len_list)))
        plt.axvline(x = dist_vec_len_mean, color = 'black', linestyle = 'dashed', label = l_mean) 
        plt.axvline(x = dist_vec_len_median, color = 'black', label = l_median) 

        plt.title(cl, fontsize = 8)

        plt.ylabel('# of combos', fontsize = 6)
        plt.yticks(fontsize = 6)
        plt.xticks(fontsize = 6)
        
        if i_subplot in [4, 5, 6]:
            plt.xlabel('Distance vector length', fontsize = 6)

        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
        for axis in ['top','right']:
            ax.spines[axis].set_linewidth(0)
        
        plt.legend(fontsize = 5)  
    
        i_subplot += 1
    
    figure_name = screen_name + '_' + syn_analysis + '_dist_vec_length_distribtuion_' + str(i_subplot - 1) + '_breast.png'
    plt.savefig(os.path.join(figure_path, figure_name), bbox_inches = 'tight', pad_inches = 0.1)


def main():

    global figure_path

    data_path = '../../data'
    results_path = '../../results'
    figure_path = os.path.join(results_path, 'figures')

    #Load DREAM target set to synergy score dictionary
    DREAM_inputs_dir = '../../data/DREAM_inputs'
    targset2syn_file = 'targset_to_syn_dict.pkl'
    DREAM_targset_to_syn_dict = pickle.load(open(os.path.join(DREAM_inputs_dir, targset2syn_file), 'rb')) #[cl][target set] = synergy

    #Reformatting so it is the same as the NCI data
    dd_list = functools.partial(defaultdict, list)
    DREAM_targset_to_targset_dict = defaultdict(dd_list)
    for cl in ['BT-20', 'BT-474', 'CAMA-1', 'MCF7', 'MDA-MB-157']: #DREAM_targset_to_syn_dict.keys():
        for targset_combo in DREAM_targset_to_syn_dict[cl]:
            [targset_1, targset_2] = targset_combo.split('__')
            ts1_list = list()
            for targ in targset_1.split('_'):
                ts1_list.append(targ)
            ts2_list = list()
            for targ in targset_2.split('_'):
                ts2_list.append(targ)
            DREAM_targset_to_targset_dict[cl][targset_combo] = [ts1_list, ts2_list]
    
            

    #Load NCI ALMANAC DB ID combos to target set
    outpath = os.path.join(data_path, 'NCI_ALMANAC_inputs')
    NCI_dbid_to_targset_dict = pickle.load(open(os.path.join(outpath, 'cl_to_dbid_targs_dict_db_022825.pkl'), 'rb')) #[cl][DB ID combo] = [[DBID1 target set], [DBID2 target set]]
    NCI_cl_to_dbid_syn_dict = pickle.load(open(os.path.join(outpath, 'cl_to_dbid_syn_dict_db_022825.pkl'), 'rb'))

    syn_analysis = 'negsyn'
    DREAM_cl_to_dist_vec_dict = get_pw_combos_per_cl(DREAM_targset_to_targset_dict, DREAM_targset_to_syn_dict, 'DREAM', syn_analysis)
    NCI_cl_to_dist_vec_dict = get_pw_combos_per_cl(NCI_dbid_to_targset_dict, NCI_cl_to_dbid_syn_dict, 'NCI_ALMANAC', syn_analysis)

    #dist_vec_dict_sorted = list(sorted(NCI_cl_to_combo_dist_vec_dict['BT-549'], key=lambda item: len(item[1]), reverse = True))
    #print(dist_vec_dict_sorted[:10])


    plot_dist_vec_len_distribution(DREAM_cl_to_dist_vec_dict, 'DREAM', syn_analysis)
    plot_dist_vec_len_distribution(NCI_cl_to_dist_vec_dict, 'NCI_ALMANAC', syn_analysis)

   


if __name__ == "__main__":
    import os, pickle, functools
    import pandas as pd
    from collections import defaultdict
    import numpy as np
    import networkx as nx
    import matplotlib.pyplot as plt
    main()