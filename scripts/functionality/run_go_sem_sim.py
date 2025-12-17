"""
07/22/2025 - calculates Wang's semantic similarity between DREAM's dominant target pairs

"""
def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)


def main():

    global interactome_path
    global figure_dir
    global exper_color
    global m
    global a

    dd_float = functools.partial(defaultdict, float)
    figure_dir = '../../results/figures/functionality'

    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    exper_source = 'DREAM'

    #LOAD SYNERGY SCORES
    global targset_to_syn_dict
    DREAM_inputs_dir = '../../data/DREAM_inputs'
    targset2syn_file = 'targset_to_syn_dict.pkl'
    targset_to_syn_dict = pickle.load(open(os.path.join(DREAM_inputs_dir, targset2syn_file), 'rb'))

    distances_results_dir = '../../results/distance_dictionaries'
    merged_nbhd_dist_filename = 'multi_target_37_dist_metric_dicts.pkl'
    merged_nbhd_results_path = os.path.join(distances_results_dir, interactome_aname, exper_source, merged_nbhd_dist_filename)
    merged_dist_dict = pickle.load(open(merged_nbhd_results_path, 'rb')) # All combos

    model_name = 'posyn'
    syn_type = 'per_targset'
    curve_fitting_results_path = os.path.join('../../results/Multi_Targ_Curve_Fitting', interactome_aname, exper_source)
    dom_targ_pair_results_dir = os.path.join(curve_fitting_results_path, 'opt_targ_pair_dict_37_metrics_' + model_name + '_' + syn_type + '.pkl')
    dom_targ_pair_dict = pickle.load(open(dom_targ_pair_results_dir, 'rb')) #[cl][dist][a] = [(combo_targset, best_targ_pair, best_targ_pair_dist)]

    m = 'N2_to_N1'
    a = 0.75
    dom_targ_pair_set = set()
    for cl in dom_targ_pair_dict.keys():
        for (_, dom_targ_pair, _) in dom_targ_pair_dict[cl][m][a]:
            dom_targ_pair_set.add(dom_targ_pair)
    
    print(len(dom_targ_pair_set))



    #INITIALIZE GOEA
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
    id2sym_dict = dict()
    sym2id_dict = dict()
    for i, row in sym2id_df.iterrows():
        if row['converted_alias'] != 'None':
            id = int(row['converted_alias'])
            sym = row['initial_alias']
            id2sym_dict[id] = sym
            sym2id_dict[sym] = id

    global goeaobj
    goeaobj = GOEnrichmentStudyNS(
        id2sym_dict.keys(), # List of PathFX protein-coding genes with Gene IDs
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    
    NS = 'BP'
    #Generate alll pairwise GO terms from dominant target pair set
    go_pair_set = set()
    for targ_pair in dom_targ_pair_set:
        T1, T2 = targ_pair.split('__')
        if T1 in sym2id_dict.keys() and T2 in sym2id_dict.keys():
            T1_id = sym2id_dict[T1]
            T2_id = sym2id_dict[T2]
            if T1_id in ns2assoc[NS].keys() and T2_id in ns2assoc[NS].keys():
                T1_terms = list(ns2assoc[NS][T1_id])
                T2_terms = list(ns2assoc[NS][T2_id])
                for T1_term_i in T1_terms:
                    for T2_term_j in T2_terms:
                        go_pair_set.add('__'.join(sorted([T1_term_i, T2_term_j])))
    print(len(go_pair_set))

    #Calculate semantic similarity between all GO term pairs
    relationships = {'part_of', 'is_a', 'regulates', 'negatively_regulates', 'positively_regulates'}
    wang_sem_sim_dict = defaultdict(dd_float)
    for go_pair in go_pair_set:
        term_1, term_2 = go_pair.split('__')
        terms_12 = [term_1, term_2]
        wang = SsWang(terms_12, godag, relationships)
        sim = wang.get_sim(term_1, term_2)
        wang_sem_sim_dict[term_1][term_2] = sim
        wang_sem_sim_dict[term_2][term_1] = sim
    
    results_dir = os.path.join('../../results/functionality', interactome_aname, exper_source)
    check_dir(results_dir)
    filename = m + '_' + str(a) + '_dom_targ_pair_' + NS + '_wang_sem_sim_dict.pkl'
    outpath = os.path.join(results_dir, filename)
    pickle.dump(wang_sem_sim_dict, open(outpath, 'wb'))

    """

    dom_targ_pair_sem_sim_dict = dict()
    dom_targ_pair_sem_sim_2_dict = dict()


    c = 0
    for targ_pair in dom_targ_pair_set:  

        c += 1  
    
        T1, T2 = targ_pair.split('__')
        if T1 in sym2id_dict.keys() and T2 in sym2id_dict.keys():
            T1_id = sym2id_dict[T1]
            T2_id = sym2id_dict[T2]
            if T1_id in ns2assoc[NS].keys() and T2_id in ns2assoc[NS].keys():
                T1_terms = list(ns2assoc[NS][T1_id])
                T2_terms = list(ns2assoc[NS][T2_id])
                T1_to_T2_term_sim = defaultdict(list)
                T2_to_T1_term_sim = defaultdict(list)
                for T1_term_i in T1_terms:
                    for T2_term_j in T2_terms:
                        
                        terms_ij = [T1_term_i, T2_term_j]
                        wang = SsWang(terms_ij, godag, relationships)
                        sim = wang.get_sim(T1_term_i, T2_term_j)

                        T1_term_i_name = godag[T1_term_i].name
                        T2_term_j_name = godag[T2_term_j].name
                        T1_term_i_depth = godag[T1_term_i].depth
                        T2_term_j_depth = godag[T2_term_j].depth
                        #print(T1_term_i, T1_term_i_name, T1_term_i_depth, T2_term_j, T2_term_j_name, T2_term_j_depth, sim)

                        T1_to_T2_term_sim[T1_term_i].append(sim)
                        T2_to_T1_term_sim[T2_term_j].append(sim)

                        wang_sem_sim_dict[T1_term_i].add((T2_term_j, sim))
                        wang_sem_sim_dict[T2_term_j].add((T1_term_i, sim))
                
                T1_best_match_sum = np.sum([max(sim_list) for sim_list in T1_to_T2_term_sim.values()])
                T2_best_match_sum = np.sum([max(sim_list) for sim_list in T2_to_T1_term_sim.values()])
                T1_best_match_avg = np.mean([max(sim_list) for sim_list in T1_to_T2_term_sim.values()])
                T2_best_match_avg = np.mean([max(sim_list) for sim_list in T2_to_T1_term_sim.values()])
                num_terms_combined = len(T1_terms) + len(T2_terms)
                T1_T2_sim = (T1_best_match_sum + T2_best_match_sum) / num_terms_combined
                T1_T2_sim_2 = (T1_best_match_avg + T2_best_match_avg) / 2
                print(targ_pair, c, T1_T2_sim, T1_T2_sim_2)

                dom_targ_pair_sem_sim_dict[targ_pair] = T1_T2_sim
                dom_targ_pair_sem_sim_2_dict[targ_pair] = T1_T2_sim_2

            else:
                print('No association in GO DAG:', T1, T2)
        else:
            print('No symbol to ID mapping:', T1, T2)
    
    results_dir = os.path.join('../../results/functionality', interactome_aname, exper_source)
    check_dir(results_dir)
    filename = m + '_' + str(a) + '_dom_targ_pair_' + NS + '_wang_sem_sim_dict.pkl'
    outpath = os.path.join(results_dir, filename)
    pickle.dump(wang_sem_sim_dict, open(outpath, 'wb'))
    """

    
    """
    for cl in dom_targ_pair_dict.keys():
        dist_list = list()
        sem_sim_list = list()
        sem_sim_2_list = list()
        syn_list = list()
        print(cl)
        for (combo_targset, dom_targ_pair, dom_targ_pair_dist) in dom_targ_pair_dict[cl][m][a]:
            
            if dom_targ_pair in dom_targ_pair_sem_sim_dict.keys():
                sem_sim = dom_targ_pair_sem_sim_dict[dom_targ_pair]
                sem_sim_2 = dom_targ_pair_sem_sim_2_dict[dom_targ_pair]
                if sem_sim < 0.8:
                    syn = targset_to_syn_dict[cl][combo_targset]
                    dist_list.append(dom_targ_pair_dist)
                    sem_sim_list.append(sem_sim)
                    sem_sim_2_list.append(sem_sim_2)
                    syn_list.append(syn)
        
        print('Similarity score 1:')
        print(len(dist_list), pearsonr(dist_list, sem_sim_list), spearmanr(dist_list, sem_sim_list))
        print(len(dist_list), pearsonr(syn_list, sem_sim_list), spearmanr(syn_list, sem_sim_list))

        print('Similarity score 2:')
        print(len(dist_list), pearsonr(dist_list, sem_sim_2_list), spearmanr(dist_list, sem_sim_2_list))
        print(len(dist_list), pearsonr(syn_list, sem_sim_2_list), spearmanr(syn_list, sem_sim_2_list))

        make_corr_plots(dist_list, sem_sim_list, '_'.join([m, str(a), 'dist']), NS + '_wang_sem_sim_wang_avg_method', cl)
        make_corr_plots(dist_list, sem_sim_2_list, '_'.join([m, str(a), 'dist']), NS + '_wang_sem_sim_pesquita_avg_method', cl)
    """

if __name__ == "__main__":
    import pickle, os, functools, math
    from collections import defaultdict
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import networkx as nx
    from scipy.stats import ttest_ind, pearsonr, spearmanr, rankdata

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