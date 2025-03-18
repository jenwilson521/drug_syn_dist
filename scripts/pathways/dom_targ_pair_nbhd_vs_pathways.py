def plot_np_overlap_vs_syn_corr(np_overlap_list, syn_list, cui, cl, m, a):
    fig = plt.figure(dpi = 300)
    ax = plt.subplot()
    fig.tight_layout()

    n_combos = 'n = ' + str(len(np_overlap_list))

    pearson_corr, pearson_pval = pearsonr(np_overlap_list, syn_list)
    spearman_corr, spearman_pval = spearmanr(np_overlap_list, syn_list)
    if pearson_pval <= 0.05:
        p_label = 'p = ' + str(np.round(pearson_corr, 4)) + '*'
    else:
        p_label = 'p = ' + str(np.round(pearson_corr, 4)) + 'ns'
    if spearman_pval <= 0.05:
        s_label = 's = ' + str(np.round(spearman_corr, 4)) + '*'
    else:
        s_label = 's = ' + str(np.round(spearman_corr, 4)) + 'ns'

    plt.scatter(np_overlap_list, syn_list, alpha = 0.5, c = 'dodgerblue', label = n_combos + '\n' + p_label + '\n' + s_label)

    plt.xlabel('nbhd-pathway overlap / nbhd size') 
    plt.ylabel('Synergy')
    plt.legend()

    figure_name = cl + '_' + m + '_' + str(a) + '_dom_pair_nbhd_' + cui + '_overlap_vs_synergy.png'
    plt.savefig(os.path.join(figure_path, 'pathways', figure_name), bbox_inches = 'tight', pad_inches = 0.1)




def main():
    global data_path
    global results_path
    global figure_path
    data_path = '../../data'
    results_path = '../../results'
    figure_path = '../../results/figures'
    resources_path = '../../resources'
    gda_path = os.path.join(data_path, 'GDA_files')

    #Load CUIs to genes dictionary
    c2g_dict = pickle.load(open(os.path.join(resources_path, 'Pfx050120_merged_unique_cuis2genes.pkl'), 'rb'))

    #Load CUIs of interest
    breast_cuis = pickle.load(open(os.path.join(resources_path, 'curated_breast_cuis.pkl'), 'rb'))
    
    #Open GDA files and save CUI -> genes dictionary
    gda_dict = defaultdict(list)
    all_gda_files = glob.glob(os.path.join(gda_path, "*.tsv"))
    for f in all_gda_files:
        gda_df = pd.read_csv(f, delimiter='\t')
        gda_dict[gda_df['Disease_id'].iloc[0]] = gda_df['Gene'].tolist()


    #Load DREAM dominant target pair dictionary
    opt_targ_pair_dict = pickle.load(open(os.path.join(results_path, 'curve_fitting', 'opt_targ_pair_dict_37_metrics_031424.pkl'), 'rb'))
    print(opt_targ_pair_dict['BT-20'].keys()) # [Cell line] = [(target set combo, dom target pair, distance)]

    #Load target set to synergy score dictionary
    DREAM_inputs_dir = '../../data/DREAM_inputs'
    targset2syn_file = 'targset_to_syn_dict.pkl'
    targset_to_syn_dict = pickle.load(open(os.path.join(DREAM_inputs_dir, targset2syn_file), 'rb'))

    #Load neighborhoods
    nbhd_file = 'PRN_nbhd_opt_alpha.pkl'
    exper_source = 'DREAM'
    DREAM_nbhds = pickle.load(open(os.path.join(results_path, 'neighborhoods', exper_source, nbhd_file), 'rb'))

    #cui = 'C0006142'
    #print(len(c2g_dict[cui]))

    #cl = 'MCF7'

    for cui in ['C0860580']: #gda_dict.keys():
        print(cui, len(gda_dict[cui]))
        for cl in ['MDA-MB-157']: #opt_targ_pair_dict.keys():
            for m in ['N1_to_N2']: #opt_targ_pair_dict[cl].keys():
                for a in [0.15]: #opt_targ_pair_dict[cl][m].keys():
                    syn_list = list()
                    merged_nbhd = set()
                    merged_nbhd_all_combos = list()
                    merged_nbhd_pathway_overlap_all_combos = list()
                    for (targset, dom_targ_pair, dist) in opt_targ_pair_dict[cl][m][a]:
                        syn_list.append(targset_to_syn_dict[cl][targset])

                        [T1, T2] = dom_targ_pair.split('__')
                        for n in DREAM_nbhds[T1]:
                            merged_nbhd.add(n)
                        for n in DREAM_nbhds[T2]:
                            merged_nbhd.add(n)
                        
                        #print(dom_targ_pair, merged_nbhd)
        
                        nbhd_pathway_overlap = [n for n in merged_nbhd if n in gda_dict[cui]]
                        merged_nbhd_pathway_overlap_all_combos.append(len(nbhd_pathway_overlap) / len(list(merged_nbhd)))

        
                    #merged_nbhd_all_combos.append(list(merged_nbhd))
                    pearson_corr, pearson_pval = pearsonr(merged_nbhd_pathway_overlap_all_combos, syn_list)
                    spearman_corr, spearman_pval = spearmanr(merged_nbhd_pathway_overlap_all_combos, syn_list)

                    if pearson_pval <= 0.05 or spearman_pval <= 0.05:
                        print(cui, cl, m, a, pearson_corr, pearson_pval, spearman_corr, spearman_pval)

    plot_np_overlap_vs_syn_corr(merged_nbhd_pathway_overlap_all_combos, syn_list, cui, cl, m, a)





if __name__ == "__main__":
    import os, pickle, functools, glob
    import pandas as pd
    from collections import defaultdict
    import numpy as np
    import networkx as nx
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr, spearmanr
    main()