"""

07/23/2025 - takes in pairwise GO term Wang semantic similarity and calculates similarity per dominant target pair.
Correlates semantic similarity with dominant target pair distance and synergy

"""
def calc_dom_targ_pair_sem_sim(wang_sem_sim_dict, dom_targ_pair_set, depth_cutoff):
    dom_targ_pair_sem_sim_dict = dict()
    for targ_pair in dom_targ_pair_set:

        T1, T2 = targ_pair.split('__')
        if T1 in sym2id_dict.keys() and T2 in sym2id_dict.keys():
            T1_id = sym2id_dict[T1]
            T2_id = sym2id_dict[T2]
            if T1_id in ns2assoc[NS].keys() and T2_id in ns2assoc[NS].keys():
                T1_terms = [term for term in ns2assoc[NS][T1_id] if godag[term].depth >= depth_cutoff]
                T2_terms = [term for term in ns2assoc[NS][T2_id] if godag[term].depth >= depth_cutoff]
                num_terms_combined = len(T1_terms) + len(T2_terms)

                if len(T1_terms) > 0 and len(T2_terms) > 0:

                    T1_to_T2_term_sim = list()
                    for T1_term_i in T1_terms:
                        T1_to_T2_term_sim.append([(T1_term_i, godag[T1_term_i].name, T2_term_j, godag[T2_term_j].name, wang_sem_sim_dict[T1_term_i][T2_term_j]) for T2_term_j in T2_terms])
                    
                    T2_to_T1_term_sim = list()
                    for T2_term_j in T2_terms:
                        T2_to_T1_term_sim.append([(T2_term_j, godag[T2_term_j].name, T1_term_i, godag[T1_term_i].name, wang_sem_sim_dict[T2_term_j][T1_term_i]) for T1_term_i in T1_terms])
                    
                    T1_best_match_sum = np.sum(sorted([np.max(list(zip(*term_pair_sim_list))[4]) for term_pair_sim_list in T1_to_T2_term_sim], reverse = True))
                    T2_best_match_sum = np.sum(sorted([np.max(list(zip(*term_pair_sim_list))[4]) for term_pair_sim_list in T2_to_T1_term_sim], reverse = True))

                    T1_best_match_avg = np.mean(sorted([np.max(list(zip(*term_pair_sim_list))[4]) for term_pair_sim_list in T1_to_T2_term_sim], reverse = True))
                    T2_best_match_avg = np.mean(sorted([np.max(list(zip(*term_pair_sim_list))[4]) for term_pair_sim_list in T1_to_T2_term_sim], reverse = True))
                    #num_terms_combined = 4
                    T1_T2_sim = (T1_best_match_sum + T2_best_match_sum) / num_terms_combined
                    T1_T2_sim_2 = (T1_best_match_avg + T2_best_match_avg) / 2

                    if targ_pair in ['ATR__CHEK1', 'BRAF__MAP2K2', 'BRAF__MAP2K7', 'AKT1__BCL2L1', 'AKT2__BCL2L1']: # and depth_cutoff == 7: #['ATR__CHEK1', 'BRAF__MAP2K2', 'BRAF__MAP2K7', 'AKT1__BCL2L1', 'AKT2__BCL2L1']:
                    #    print('Number of terms', len(T1_terms), len(T2_terms))
                        #print(targ_pair, T1_to_T2_term_sim, T2_to_T1_term_sim, T1_T2_sim)
                        print(targ_pair, sorted([np.max(list(zip(*term_pair_sim_list))[4]) for term_pair_sim_list in T1_to_T2_term_sim], reverse = True), sorted([np.max(list(zip(*term_pair_sim_list))[4]) for term_pair_sim_list in T2_to_T1_term_sim], reverse = True))

                    dom_targ_pair_sem_sim_dict[targ_pair] = T1_T2_sim
                else:
                    continue
                    #print(targ_pair, 'has no terms at this depth cutoff')
    return dom_targ_pair_sem_sim_dict

def sem_sim_across_depth_cutoffs(wang_sem_sim_dict, depth_cutoff_all, NS):
    cl_to_depth_to_sem_dist_corr_dict = defaultdict(dd_float)
    cl_to_depth_to_sem_syn_corr_dict = defaultdict(dd_float)
    cl_to_depth_to_num_combos_dict = defaultdict(dd_float)
    for depth_cutoff in depth_cutoff_all:
        print('Depth cutoff:', depth_cutoff)
        dom_targ_pair_sem_sim_dict = calc_dom_targ_pair_sem_sim(wang_sem_sim_dict, dom_targ_pair_set, depth_cutoff)

        
        for cl in dom_targ_pair_dict.keys():
            dist_list = list()
            sem_sim_list = list()
            syn_list = list()
            print(cl)
            for (combo_targset, dom_targ_pair, dom_targ_pair_dist) in dom_targ_pair_dict[cl][m][a]:
                
                if dom_targ_pair in dom_targ_pair_sem_sim_dict.keys():
                    sem_sim = dom_targ_pair_sem_sim_dict[dom_targ_pair]
                    syn = targset_to_syn_dict[cl][combo_targset]
                    dist_list.append(dom_targ_pair_dist)
                    sem_sim_list.append(sem_sim)
                    syn_list.append(syn)
            
            num_combos = len(dist_list)
            cl_to_depth_to_num_combos_dict[cl][depth_cutoff] = num_combos
            if num_combos <= 1:
                cl_to_depth_to_sem_dist_corr_dict[cl][depth_cutoff] = 0
                cl_to_depth_to_sem_syn_corr_dict[cl][depth_cutoff] = 0
            else:
                ds_pearson_corr, pval = pearsonr(dist_list, sem_sim_list)
                ss_pearson_corr, pval = pearsonr(syn_list, sem_sim_list)

                ds_spearman_corr, pval = spearmanr(dist_list, sem_sim_list)
                ss_spearman_corr, pval = spearmanr(syn_list, sem_sim_list)

                cl_to_depth_to_sem_dist_corr_dict[cl][depth_cutoff] = ds_spearman_corr
                cl_to_depth_to_sem_syn_corr_dict[cl][depth_cutoff] = ss_spearman_corr

                print(len(dist_list), pearsonr(dist_list, sem_sim_list), spearmanr(dist_list, sem_sim_list))
                print(len(dist_list), pearsonr(syn_list, sem_sim_list), spearmanr(syn_list, sem_sim_list))
    
    return cl_to_depth_to_num_combos_dict, cl_to_depth_to_sem_dist_corr_dict, cl_to_depth_to_sem_syn_corr_dict


def plot_corr_across_depth_cutoff(cl_to_depth_dict, y_title, title):
    fig = plt.figure(dpi = 300)
    fig.tight_layout()
    ax = plt.subplot()

    for cl in cl_to_depth_dict.keys():
        depth_cutoff_all = list(cl_to_depth_dict[cl].keys())
        val_all = list(cl_to_depth_dict[cl].values())

        plt.plot(depth_cutoff_all, val_all, marker = 'o', color = cl_to_color_dict[cl], label = cl)
    
    #Figure setting specifications
    for axis in ['left', 'bottom']:
        ax.spines[axis].set_linewidth(1)
    for axis in ['top','right']:
        ax.spines[axis].set_linewidth(0)                            

    plt.grid()
    plt.legend()
    ax.tick_params(axis = 'y', width = 1)
    ax.tick_params(axis = 'x', width = 1)
    #plt.setp(ax, xlim = [0, 3])
    #plt.xticks([0.5, 1.5, 2.5], labels = [])
    #plt.yticks(np.arange(-0.55, -0.35, 0.05), labels = [])
    plt.ylabel(y_title)
    plt.xlabel('Depth cutoff')

    figure_name = title + '_all_cl.png'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches = 'tight', pad_inches = 0.1) 


def first_degree_poly(dist, a, b):
    return a * dist + b

def fit_curve(x_list, y_list):

    """
        Uses selected curve type (e.g., first-degree polynomial) to generate curve parameters

        :param dist_syn_list: list of tuples
            [(distance, synergy),...]

        :return popt: array
            Optimal parameters from curve fitting
        
        :return pcov: 2-D array
            Covariance matrix of the optimal parameters from curve fitting
    """


    popt, pcov = curve_fit(first_degree_poly, x_list, y_list)
    #popt, pcov = curve_fit(second_degree_poly, dist, syn)
    #popt, pcov = curve_fit(third_degree_poly, dist, syn)
    #popt, pcov = curve_fit(exp_func, dist, syn)

    return popt, pcov

def make_corr_plots(x_list, y_list, x_title, y_title, cl):
    fig = plt.figure(dpi = 300)
    fig.tight_layout() 
    ax = plt.subplot()

    popt, _ = fit_curve(x_list, y_list)
    x = np.arange(min(x_list), max(x_list), 1)
    y = [(x_i * popt[0]) + popt[1] for x_i in x]

    plt.scatter(x_list, y_list, 110, color = 'dodgerblue', alpha = 0.5, zorder = 1)
    plt.plot(x, y, color = 'dodgerblue', linewidth = 1.5, zorder = 0)

    for axis in ['bottom','left']:
        ax.spines[axis].set_linewidth(3)
    for axis in ['top', 'right']:
        ax.spines[axis].set_linewidth(0)
    
    ax.tick_params(width = 3, length = 6)
    if 'dist' in x_title:
        plt.xticks(np.arange(2, 18, 2), labels = [])
    elif 'syn' in x_title:
        plt.xticks(np.arange(0, 120, 10), labels = [])
    
    plt.yticks(np.arange(0, 1.2, 0.2), labels = [])

    figure_name = '_'.join(['__'.join([m, str(a)]), cl, x_title, 'vs', y_title, 'corr.png'])
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches = 'tight', pad_inches = 0.1) 



def main():

    global figure_dir
    global m
    global a

    global dd_float
    dd_float = functools.partial(defaultdict, float)
    figure_dir = '../../results/figures/functionality/network_pharma'

    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    exper_source = 'DREAM'

    global cl_to_color_dict
    cl_to_color_dict = {'BT-474': 'navy', 
                            'CAMA-1': 'blue',
                            'MCF7': 'dodgerblue', 
                            'MDA-MB-157': 'lightblue',
                            'BT-20': 'indigo'}

    #LOAD SYNERGY SCORES
    global targset_to_syn_dict
    DREAM_inputs_dir = '../../data/DREAM_inputs'
    targset2syn_file = 'targset_to_syn_dict.pkl'
    targset_to_syn_dict = pickle.load(open(os.path.join(DREAM_inputs_dir, targset2syn_file), 'rb'))

    #LOAD DISTANCES
    distances_results_dir = '../../results/distance_dictionaries'
    merged_nbhd_dist_filename = 'multi_target_37_dist_metric_dicts.pkl'
    merged_nbhd_results_path = os.path.join(distances_results_dir, interactome_aname, exper_source, merged_nbhd_dist_filename)
    merged_dist_dict = pickle.load(open(merged_nbhd_results_path, 'rb')) # All combos

    global dom_targ_pair_dict
    model_name = 'posyn'
    syn_type = 'per_targset'
    curve_fitting_results_path = os.path.join('../../results/Multi_Targ_Curve_Fitting', interactome_aname, exper_source)
    dom_targ_pair_results_dir = os.path.join(curve_fitting_results_path, 'opt_targ_pair_dict_37_metrics_' + model_name + '_' + syn_type + '.pkl')
    dom_targ_pair_dict = pickle.load(open(dom_targ_pair_results_dir, 'rb')) #[cl][dist][a] = [(combo_targset, best_targ_pair, best_targ_pair_dist)]

    global dom_targ_pair_set
    m = 'N2_to_N1'
    a = 0.75
    dom_targ_pair_set = set()
    for cl in dom_targ_pair_dict.keys():
        for (_, dom_targ_pair, _) in dom_targ_pair_dict[cl][m][a]:
            dom_targ_pair_set.add(dom_targ_pair)
    
    #LOAD SEMANTIC SIMILARITIES
    global NS
    NS = 'BP'
    sem_sim_results_dir = os.path.join('../../results/functionality', interactome_aname, exper_source)
    filename = m + '_' + str(a) + '_dom_targ_pair_' + NS + '_wang_sem_sim_dict.pkl'
    wang_sem_sim_dict = pickle.load(open(os.path.join(sem_sim_results_dir, filename), 'rb'))
    

    #INITIALIZE GOEA
    global obodag
    global godag
    global ns2assoc
    fin_gene2go = download_ncbi_associations() #download associations
    obodag = GODag("../../resources/go-basic.obo") #load ontologies

    #Not sure if this is the same as the line above...
    godag = get_godag("../../resources/go-basic.obo", optional_attrs={'relationship'})

    # Read NCBI's gene2go. Store annotations in a list of namedtuples
    objanno = Gene2GoReader(fin_gene2go, taxids=[9606]) #homo sapien taxid = 9606
    ns2assoc = objanno.get_ns2assc() #Dictionary with MF, CC, BP keys and gene -> GO term annotations ex: 109617024: {'GO:0006396'}
    
    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated homo sapien genes". format(NS = nspc, N = len(id2gos)))

    #Load background genes (PathFX 0.5 threshold interactome)
    sym2id_df = pd.read_csv('../../data/0.5thr_PathFX_int_GeneSymbol2GeneID.csv', sep = ',')

    #Create dictionary that maps gene ID to symbol in order to convert target neighborhoods to gene IDs
    global id2sym_dict
    global sym2id_dict
    id2sym_dict = dict()
    sym2id_dict = dict()
    for i, row in sym2id_df.iterrows():
        if row['converted_alias'] != 'None':
            id = float(row['converted_alias'])
            sym = row['initial_alias']
            id2sym_dict[id] = sym
            sym2id_dict[sym] = id

    """
    global goeaobj
    goeaobj = GOEnrichmentStudyNS(
        id2sym_dict.keys(), # List of PathFX protein-coding genes with Gene IDs
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    """

    #Calculate semantic similarity between proteins in all dominant target pairs
    depth_cutoff_all = np.arange(0, 11, 1)
    #cl_to_depth_to_num_combos_dict, cl_to_depth_to_sem_dist_corr_dict, cl_to_depth_to_sem_syn_corr_dict  = sem_sim_across_depth_cutoffs(wang_sem_sim_dict, depth_cutoff_all, NS)
     
    
    #print(cl_to_depth_to_sem_dist_corr_dict)
    #print(cl_to_depth_to_sem_syn_corr_dict)

    #plot_corr_across_depth_cutoff(cl_to_depth_to_sem_dist_corr_dict, 'Spearman corr', NS + '_sem_sim_vs_' + m + '_' + str(a) + '_dist_spearman_corr_pes_avg_method')
    #plot_corr_across_depth_cutoff(cl_to_depth_to_sem_syn_corr_dict, 'Spearman corr', NS + '_sem_sim_vs_' + m + '_' + str(a) + '_syn_spearman_corr_pes_avg_method')
    #plot_corr_across_depth_cutoff(cl_to_depth_to_num_combos_dict, 'Num combos', NS + '_num_combos_pes_avg_method')


    #Depth cutoff chosen from parameter serach for each NS
    ns_to_depth_cutoff_dict = {'MF': 5,
                            'BP': 4}
    
    dom_targ_pair_sem_sim_dict = calc_dom_targ_pair_sem_sim(wang_sem_sim_dict, dom_targ_pair_set, ns_to_depth_cutoff_dict[NS])

    for cl in dom_targ_pair_dict.keys():
        dist_list = list()
        sem_sim_list = list()
        syn_list = list()
        print(cl)
        for (combo_targset, dom_targ_pair, dom_targ_pair_dist) in dom_targ_pair_dict[cl][m][a]:
            
            if dom_targ_pair in dom_targ_pair_sem_sim_dict.keys():
                sem_sim = dom_targ_pair_sem_sim_dict[dom_targ_pair]
                syn = targset_to_syn_dict[cl][combo_targset]
                if syn < 120:
                    dist_list.append(dom_targ_pair_dist)
                    sem_sim_list.append(sem_sim)
                    syn_list.append(syn)
        
        print(len(dist_list), 'Dist-semsim corr:', pearsonr(dist_list, sem_sim_list), spearmanr(dist_list, sem_sim_list))
        print(len(syn_list), 'Syn-semsim corr:', pearsonr(syn_list, sem_sim_list), spearmanr(syn_list, sem_sim_list))
        
        make_corr_plots(dist_list, sem_sim_list, '_'.join([m, str(a), 'dist']), NS + '_wang_sem_sim_wang_avg_method', cl)
        make_corr_plots(syn_list, sem_sim_list, '_'.join([m, str(a), 'syn']), NS + '_wang_sem_sim_wang_avg_method', cl)
            




if __name__ == "__main__":
    import pickle, os, functools, math
    from collections import defaultdict
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import networkx as nx
    from scipy.stats import ttest_ind, pearsonr, spearmanr, rankdata
    from scipy.optimize import curve_fit

    #Download ontologies
    from goatools.base import download_go_basic_obo
    #Download associations
    from goatools.base import download_ncbi_associations
    #Load ontologies
    from goatools.obo_parser import GODag
    #Load associations
    from goatools.anno.genetogo_reader import Gene2GoReader
    #Initialize GOEA object
    from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

    #For Wang similarity
    from goatools.base import get_godag
    from goatools.semsim.termwise.wang import SsWang
    
    main()