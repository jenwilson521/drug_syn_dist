"""
05/29/2024 - cell line specific correlation plots 
            and script to extract neighborhood overlap info

"""

def make_corr_plots(merged_dist_dict, opt_targ_pair_dict):

    for cl in opt_targ_pair_dict.keys():

        fig = plt.figure(dpi = 300)
        fig.tight_layout() 
        ax = plt.subplot(2, 1, 1)

        pre_dist = list(zip(*merged_dist_dict[cl][m][a]))[1]
        pre_syn = list(zip(*merged_dist_dict[cl][m][a]))[2]

        post_dist = list(zip(*opt_targ_pair_dict[cl][m][a]))[2]
        post_syn = [targset_to_syn_dict[cl][ts] for (ts, combo, dist) in opt_targ_pair_dict[cl][m][a]]

        print(cl)
        print('PRE-FIT:', pearsonr(pre_dist, pre_syn), spearmanr(pre_dist, pre_syn))
        print('POST-FIT:', pearsonr(post_dist, post_syn), spearmanr(post_dist, post_syn))

        pre_fit_m = pre_params_dict[cl][m][a][2][0]
        pre_fit_b = pre_params_dict[cl][m][a][2][1]

        post_fit_m = post_params_dict[cl][m][a][2][0]
        post_fit_b = post_params_dict[cl][m][a][2][1]

        x_all = np.arange(3, 18) #np.arange(min(dist), max(dist) + 1)
        y_pre = [(pre_fit_m * x) + pre_fit_b for x in x_all]
        y_post = [(post_fit_m * x) + post_fit_b for x in x_all]

        plt.plot(x_all, y_pre, color = 'black', linewidth = 1.5, zorder = 0)

        for (ts, dist, syn) in merged_dist_dict[cl][m][a]:
            if syn < 120:
                if ts == 'AKT1_AKT2__BCL2L1':
                    plt.scatter(dist, syn, 70, color = 'purple', edgecolors = 'black', zorder = 2)
                elif ts == 'ATR__CHEK1':
                    plt.scatter(dist, syn, 70, color = 'orange', edgecolors = 'black', marker = 'o', zorder = 2)
                elif ts == 'BRAF__MAP2K1_MAP2K2_MAP2K3_MAP2K4_MAP2K5_MAP2K6_MAP2K7':
                    plt.scatter(dist, syn, 70, color = 'forestgreen', edgecolors = 'black', marker = 'o', zorder = 2)
                else:
                    plt.scatter(dist, syn, color = 'black', marker = 'o', alpha = 0.5, zorder = 1)

        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
        for axis in ['top', 'right']:
            ax.spines[axis].set_linewidth(0)
        
        #for label in ax.yaxis.get_ticklabels()[::2]:
        #    label.set_visible(False)   

        ax.tick_params(width = 1.5)
        plt.xticks(np.arange(0, 20, 4), labels = [])
        plt.yticks(np.arange(0, 120, 10), labels = [])  

        ax = plt.subplot(2, 1, 2)          

        plt.plot(x_all, y_post, color = exper_color, linewidth = 1.5, zorder = 0)

        for (ts, combo, dist) in opt_targ_pair_dict[cl][m][a]:
            syn = targset_to_syn_dict[cl][ts]
            if syn < 120:
                if ts == 'AKT1_AKT2__BCL2L1':
                    plt.scatter(dist, syn, 70, color = 'purple', edgecolors = exper_color, marker = 'o', zorder = 2)
                elif ts == 'ATR__CHEK1':
                    plt.scatter(dist, syn, 70, color = 'orange', edgecolors = exper_color, marker = 'o', zorder = 2)
                elif ts == 'BRAF__MAP2K1_MAP2K2_MAP2K3_MAP2K4_MAP2K5_MAP2K6_MAP2K7':
                    plt.scatter(dist, syn, 70, color = 'forestgreen', edgecolors = exper_color, marker = 'o', zorder = 2)
                else:
                    plt.scatter(dist, syn, color = exper_color, marker = 'o', alpha = 0.5, zorder = 1)            

        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
        for axis in ['top', 'right']:
            ax.spines[axis].set_linewidth(0)
        
        #for label in ax.yaxis.get_ticklabels()[::2]:
        #    label.set_visible(False)

        ax.tick_params(width = 1.5)
        plt.xticks(np.arange(0, 20, 4), labels = [])
        plt.yticks(np.arange(0, 120, 10), labels = [])

        figure_name = cl + '_pre_post_corr.png'
        plt.savefig(os.path.join(figure_path, figure_name), bbox_inches = 'tight', pad_inches = 0.1) 

def get_dp_nbhd_overlap(targ_pair):

    local_thresh = math.exp(-(avg_dist - (2 * std_dist)))

    T1, T2 = targ_pair.split('__')

    T1_nbhd = prn_nbhds[T1]

    for n in T1_nbhd:
        print(n)
    print('BREAK')
    T2_nbhd = prn_nbhds[T2]

    print(targ_pair)
    print('T1 neighborhood size', len(T1_nbhd))
    print('T2 neighborhood size', len(T2_nbhd))

    T2_dp = rwr_vecs[T2]
    T1_dp = rwr_vecs[T1]

    T1_proba_int_nodes = zip(int_nodes_sorted, T1_dp)
    T2_proba_int_nodes = zip(int_nodes_sorted, T2_dp)


    T1_proba_int_nodes_dict = {n: p for (n, p) in T1_proba_int_nodes} #if n in T1_nbhd or n in T2_nbhd}
    T2_proba_int_nodes_dict = {n: p for (n, p) in T2_proba_int_nodes} #if n in T1_nbhd or n in T2_nbhd}

    T1_proba_list = [-1/np.log(p) for (n, p) in T1_proba_int_nodes_dict.items()] # if n in T1_nbhd or n in T2_nbhd]
    T2_proba_list = [-1/np.log(p) for (n, p) in T2_proba_int_nodes_dict.items()] # if n in T1_nbhd or n in T2_nbhd]

    #L2_norm = np.linalg.norm(np.array(T1_proba_list)- np.array(T2_proba_list))
    #print('L2 norm:', L2_norm)

    local_pert = [n for (n, p) in T2_proba_int_nodes_dict.items() if p >= local_thresh]
    print('T2 local pert size:', len(local_pert))
    #for n in local_pert:
    #    print(n)
    
    print('BREAK')
    local_pert_nbhd_overlap = [n for n in local_pert if n in T1_nbhd]

    print(len(local_pert_nbhd_overlap), local_pert_nbhd_overlap)
    for n in local_pert_nbhd_overlap:
        print(n)

    return local_pert_nbhd_overlap, T1_nbhd, local_pert


def get_dp_dp_overlap(targ_pair):
    local_thresh = math.exp(-(avg_dist - (2 * std_dist)))

    T1, T2 = targ_pair.split('__')

    T1_nbhd = prn_nbhds[T1]
    T2_nbhd = prn_nbhds[T2]

    print(targ_pair)
    print('T1 neighborhood size', len(T1_nbhd))
    print('T2 neighborhood size', len(T2_nbhd))

    T2_dp = rwr_vecs[T2]
    T1_dp = rwr_vecs[T1]

    T1_proba_int_nodes = zip(int_nodes_sorted, T1_dp)
    T2_proba_int_nodes = zip(int_nodes_sorted, T2_dp)


    T1_proba_int_nodes_dict = {n: p for (n, p) in T1_proba_int_nodes if n in T1_nbhd or n in T2_nbhd}
    T2_proba_int_nodes_dict = {n: p for (n, p) in T2_proba_int_nodes if n in T1_nbhd or n in T2_nbhd}

    T1_proba_list = [-1/np.log(p) for (n, p) in T1_proba_int_nodes_dict.items()] # if n in T1_nbhd or n in T2_nbhd]
    T2_proba_list = [-1/np.log(p) for (n, p) in T2_proba_int_nodes_dict.items()] # if n in T1_nbhd or n in T2_nbhd]

    L2_norm = np.linalg.norm(np.array(T1_proba_list)- np.array(T2_proba_list))
    #print('L2 norm:', L2_norm)

    T2_local_pert = [n for (n, p) in T2_proba_int_nodes_dict.items() if p >= local_thresh]
    T1_local_pert = [n for (n, p) in T1_proba_int_nodes_dict.items() if p >= local_thresh]
    for n in T2_local_pert:
        print(n)
    print('break')
    for n in T1_local_pert:
        print(n)

    local_pert_overlap = [n for n in T2_local_pert if n in T1_lcaol_pert]

    print(len(local_pert_overlap), local_pert_overlap)
    for n in local_pert_overlap:
        print(n)

    return local_pert_overlap


def run_goea(aname, protein_list, depth_cutoff):
        
    """
        Runs GOEA for a target's PageRank-Nibble neighborhood
        Saves ALL (even those with Benjamini-Hochberg corrected p-val > 0.05) in text file (one text file/target)
        Saves significant results (GO term and Benjamini-Hochberg corrected p-value) in dictionary

    
    """

    prn_ids = [id for id, sym in id2sym_dict.items() if sym in protein_list] #Neighborhood gene symbol -> gene ID
    print('Number of genes with IDs:', len(prn_ids))
    goea_results_all = goeaobj.run_study(prn_ids, prt = None) #Run GOEA
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05] #Retain sig. results
    go_pval_tuple = [(r.GO, r.name, r.p_fdr_bh) for r in goea_results_all
                      if r.p_fdr_bh < 0.05 and r.NS == 'BP' and r.depth >= depth_cutoff] #Save [(GO, pval)] for sig. BP results

    #textfile_name = targ + '_PRN_opt_alpha_GOEA.txt'
    #goeaobj.wr_txt(os.path.join(goea_results_path, textfile_name), goea_results_sig)

    goea_dict[aname] = go_pval_tuple


def dist_syn_plot_by_num_go_terms(cl_to_opt_targ_pair_dict, targ_pair_to_marker_size_dict):
    fig = plt.figure(dpi = 300, figsize = (3, 5))
    fig.tight_layout() 
    ax = plt.subplot(2, 1, 1)

    for cl in cl_to_opt_targ_pair_dict.keys():
        for (targ_pair, dist, syn) in cl_to_opt_targ_pair_dict[cl]:

            plt.scatter(dist, syn, targ_pair_to_marker_size_dict[targ_pair], 
            edgecolors = targ_pair_to_color_dict[targ_pair], color = cl_to_edge_color_dict[cl], linewidth = 1, alpha = 0.7)
    
    for axis in ['bottom','left']:
        ax.spines[axis].set_linewidth(1)
    for axis in ['top', 'right']:
        ax.spines[axis].set_linewidth(0)
    
    ax.tick_params(width = 1)
    plt.xticks(np.arange(2, 14, 2), labels = [])
    plt.yticks(np.arange(0, 100, 20), labels = [])
    
    figure_name = '__'.join([m, str(a)]) + '_dist_syn_plot_by_num_go_terms.png'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches = 'tight', pad_inches = 0.1) 

def goea_hm(term_to_pval_dict, combo):
    targs_and_combos = list(term_to_pval_dict.keys())
    terms = list(term_to_pval_dict[targs_and_combos[0]].keys())

    #Clustering terms according to their Wang semantic similarity
    relationships = {'part_of', 'is_a', 'regulates', 'negatively_regulates', 'positively_regulates'}
    sim_matrix = np.empty((len(terms), len(terms)))
    for i, term1 in enumerate(terms):
        for j, term2 in enumerate(terms):
            term1_split = term1.split(' ')
            term2_split = term2.split(' ')
            go1 = term1_split[0]
            go2 = term2_split[0]
            wang = SsWang([go1, go2], godag, relationships)
            sim = wang.get_sim(go1, go2)
            #sim = lin_sim(go1, go2, obodag, termcounts)
            #sim = resnik_sim(go1, go2, obodag, termcounts)
            #mbl = min_branch_length(go1, go2, obodag, None)
            sim_matrix[i, j] = 1 - sim
            #sim_matrix[i, j] = mbl
    
    sim_matrix_condensed = pdist(sim_matrix)
    dendogram = hierarchy.ward(sim_matrix_condensed)
    clusters = hierarchy.fcluster(dendogram, 0.4, criterion = 'distance')
    term_cluster_tuples = [(t, c) for (t, c) in zip(terms, clusters)]
    terms_sorted = list(zip(*sorted(term_cluster_tuples, key = lambda x:x[1])))[0]
    clusters_sorted = list(zip(*sorted(term_cluster_tuples, key = lambda x:x[1])))[1]
    print(terms_sorted)
    print(clusters_sorted)
    
    hm_matrix = np.empty((len(terms), len(targs_and_combos)))
    for i, tc in enumerate(targs_and_combos):
        for j, t in enumerate(terms_sorted):
            pval = term_to_pval_dict[tc][t]
            if pval != 1:
                hm_matrix[j, i] = -math.log10(pval)
            else:
                hm_matrix[j, i] = np.nan
        #hm_matrix[i, j] = [-math.log10(pval) if pval != 1 else np.nan for pval in term_to_pval_dict[tc].values()]

        #pvals_all.append([-math.log10(pval) for pval in term_to_pval_dict[t].values()])
    
    #enrich_per_term = np.sum(hm_matrix, axis=1)
    #sorted_term_indices = np.argsort(enrich_per_term)[::-1]
    #hm_matrix_sorted = hm_matrix[sorted_term_indices]
    #indices_terms_sorted = sorted(zip(sorted_term_indices, terms))
    #terms_sorted = [t for _, t in indices_terms_sorted]

    hm_df = pd.DataFrame(hm_matrix)
    hm_df['terms'] = terms_sorted
    hm_df['non_nan_count'] = hm_df.notna().sum(axis=1)
    sorted_hm_df = hm_df.sort_values(by='non_nan_count', ascending=False)
    print(sorted_hm_df)
    terms_sorted_by_count = sorted_hm_df['terms'].to_list()   
    sorted_hm_df = sorted_hm_df.drop(columns=['non_nan_count', 'terms'])

    fig = plt.figure(dpi = 300, figsize = (3, 5)) #, figsize = (6, 10)) #, figsize = (3, 5))
    fig.tight_layout()
    ax = plt.subplot()
    cmap = combo_to_hm_color_dict[combo]
    cmap.set_bad('gray')
    hm = sns.heatmap(sorted_hm_df, cmap = cmap, linecolor = 'black', linewidths = 1.5) # xticklabels = False, yticklabels = False)
    
    for _, spine in hm.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1.5)
    
    ax.set_yticks(range(0, len(terms)), labels = []) #terms_sorted_by_count)
    ax.set_xticks(range(0, len(targs_and_combos)), labels = []) #targs_and_combos)
    plt.yticks(rotation = 0)
    plt.xticks(rotation = 90)
    ax.tick_params(width = 1.5)

    xticks = ax.get_xticks()
    yticks = ax.get_yticks()

    # Shift tick positions to center
    ax.set_xticks(xticks + 0.5)
    ax.set_yticks(yticks + 0.5)

    figure_name = '__'.join([m, str(a)]) + '_dist_' + combo + '_combo_top_terms_goea_hm.png'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches = 'tight', pad_inches = 0.1)

    
def make_goea_tables(targ, protein_list):
    """
        Runs GOEA for single target neighborhood

        :param targ: string
            Single target or combo target
        
        :param protein_list: list
            List of proteins associated with targ for GOEA
        
        :return rows: list
            List of enrichment results for all relevant GO terms, to be converted to a dataframe
    
    """

    protein_ids = [id for id, sym in id2sym_dict.items() if sym in protein_list] #Neighborhood gene symbol -> gene ID
    goea_results_all = goeaobj.run_study(protein_ids, prt = None) #Run GOEA
    rows = list()
    for r in goea_results_all:
        if r.p_fdr_bh <= 0.05:
            ratio_in_study = r.ratio_in_study
            ratio_in_study_str = str(ratio_in_study[0]) + '/' + str(ratio_in_study[1])
            ratio_in_pop = r.ratio_in_pop
            ratio_in_pop_str = str(ratio_in_pop[0]) + '/' + str(ratio_in_pop[1])

            if len(r.study_items) > 0:

                gene_sym = [sym for id, sym in id2sym_dict.items() if id in r.study_items]
                gene_sym_str = ','.join(gene_sym)
            else:
                gene_sym_str = ''

            row_data = [r.GO, r.NS, r.p_fdr_bh, ratio_in_study_str, ratio_in_pop_str, 
                        r.depth, r.name, gene_sym_str]

            rows.append(row_data)
    
    return rows

def run_make_goea_tables(dp_nbhd_overlap_dict):
    """
        Runs GOEA for list of targets and saves results in one spreadsheet

        :param dp_nbhd_overlap_dict: dictionary
        
    """
    col_names = ['GO term', 'Aspect', 'BH-corrected p-val', 'Ratio of term in protein list',
                     'Ratio of term in interactome', 'Depth', 'Name', 'Term-protein list overlapping proteins']
    targ_df_dict = {targ: pd.DataFrame() for targ in dp_nbhd_overlap_dict.keys()}
    for targ in dp_nbhd_overlap_dict.keys():
        protein_list = dp_nbhd_overlap_dict[targ]
        rows = make_goea_tables(targ, protein_list)
        targ_df = pd.DataFrame(rows, columns = col_names)
        targ_df_dict[targ] = targ_df
    
    filename = 'Example_dominant_target_pair_GOEA_tables.xlsx'
    with pd.ExcelWriter(os.path.join(functionality_results_dir, filename), engine='xlsxwriter') as writer:
            for targ in targ_df_dict.keys():
                #Manually updating spreadsheet names
                if '__' in targ:
                    aname = " ".join([targ, 'overlap'])
                elif targ in ['CHEK1', 'BRAF', 'BCL2L1']:
                    aname = " ".join([targ, 'neighborhood'])
                elif targ in ['ATR', 'MAP2K2', 'MAP2K7', 'AKT1', 'AKT2']:
                    aname = " ".join([targ, 'local perturbation'])

                targ_df_dict[targ].to_excel(writer, sheet_name = aname)

def main():

    global interactome_path
    global figure_dir
    global exper_color
    global m
    global a

    dd_float = functools.partial(defaultdict, float)
    figure_dir = '../../results/figures/functionality/dom_targ_pair'

    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    exper_source = 'DREAM'

    global functionality_results_dir
    functionality_results_dir = os.path.join('../../results', 'functionality', interactome_aname, exper_source)

    #LOAD SYNERGY SCORES
    global targset_to_syn_dict
    DREAM_inputs_dir = '../../data/DREAM_inputs'
    targset2syn_file = 'targset_to_syn_dict.pkl'
    targset_to_syn_dict = pickle.load(open(os.path.join(DREAM_inputs_dir, targset2syn_file), 'rb'))

    #for cl in targset_to_syn_dict.keys():
    #    print(cl)
    #    for (targset, syn) in targset_to_syn_dict[cl].items():
    #        if 'ALK_MET' in targset:
    #            print(targset, syn)

    distances_results_dir = '../../results/distance_dictionaries'
    merged_nbhd_dist_filename = 'multi_target_37_dist_metric_dicts.pkl'
    merged_nbhd_results_path = os.path.join(distances_results_dir, interactome_aname, exper_source, merged_nbhd_dist_filename)
    merged_dist_dict = pickle.load(open(merged_nbhd_results_path, 'rb')) # All combos

    model_name = 'posyn'
    syn_type = 'per_targset'
    curve_fitting_results_path = os.path.join('../../results/Multi_Targ_Curve_Fitting', interactome_aname, exper_source)
    opt_targ_pair_results_dir = os.path.join(curve_fitting_results_path, 'opt_targ_pair_dict_37_metrics_' + model_name + '_' + syn_type + '.pkl')
    opt_targ_pair_dict = pickle.load(open(opt_targ_pair_results_dir, 'rb')) #[cl][dist][a] = [(combo_targset, best_targ_pair, best_targ_pair_dist)]

    global pre_params_dict
    global post_params_dict
    pre_params_filename = 'pre_fit_params_dict_minfunc_' + model_name + '_' + syn_type + '.pkl'
    post_params_filename = 'post_fit_params_dict_minfunc_' + model_name + '_' + syn_type + '.pkl'
    pre_params_dict = pickle.load(open(os.path.join(curve_fitting_results_path, pre_params_filename), 'rb'))
    post_params_dict = pickle.load(open(os.path.join(curve_fitting_results_path, post_params_filename), 'rb'))

    exper_color = 'dodgerblue'
    m = 'T2_to_N1'
    if 'T' in m.split('_')[0]:
        starting_choice = 'targ'
    elif 'N' in m.split('_')[0]:
        starting_choice = 'nbhd'
    a = 0.75
    #make_corr_plots(merged_dist_dict, opt_targ_pair_dict)

    #LOAD NEIGHBORHOODS
    global prn_nbhds
    nbhd_dir = '../../results/neighborhoods'
    nbhd_filename = 'PRN_nbhd_opt_alpha.pkl'
    nbhd_path = os.path.join(nbhd_dir, interactome_aname, exper_source, nbhd_filename)
    prn_nbhds = pickle.load(open(nbhd_path, 'rb'))

    #LOAD RWR VECS
    global rwr_vecs
    rwr_vec_dir = '../../results/RWR_diffusion_profiles'
    rwr_vec_filename = 'single_and_multi_target_RWR_diffusion_profiles_' + starting_choice + '_starting_nodes_alphaD_' + str(a) + '.npz'
    rwr_vec_path = os.path.join(rwr_vec_dir, interactome_aname, exper_source, rwr_vec_filename)
    rwr_vecs = np.load(rwr_vec_path)

    #LOAD FULL DISTANCE DISTRIBUTIONS
    single_targ_dist_filename = 'single_target_37_dist_metric_dicts.pkl'
    single_targ_results_path = os.path.join(distances_results_dir, interactome_aname, exper_source, single_targ_dist_filename)
    combined_pw_distance_lists_load = pickle.load(open(single_targ_results_path, 'rb'))

    
    #07/18/2025 - Initial analysis was done exlcuding targets combined with themselves (n = 7260), need to remove those from the distance lists to ensure distribution is the same
    combined_pw_distance_lists = {m: {a: {targ_pair: combined_pw_distance_lists_load[m][a][targ_pair] for targ_pair in combined_pw_distance_lists_load[m][a].keys() 
                                if targ_pair.split('__')[0] != targ_pair.split('__')[1]} for a in combined_pw_distance_lists_load[m].keys()} 
                                for m in combined_pw_distance_lists_load.keys()}


    global avg_dist
    global std_dist
    print('NUMBER OF DISTANCES:', len(combined_pw_distance_lists['T2_to_T1'][a].values()))
    avg_dist = np.mean(list(combined_pw_distance_lists['T2_to_T1'][a].values()))
    std_dist = np.std(list(combined_pw_distance_lists['T2_to_T1'][a].values()))
    for (targ_pair, dist) in combined_pw_distance_lists['T2_to_T1'][a].items():
        if np.isinf(dist):
            print(targ_pair, dist) 
        if np.isnan(dist):
            print(targ_pair, dist)
        
    print('BRAF__MAP2K7 N2_to_N1__0.75 distance:', combined_pw_distance_lists[m][a]['BRAF__MAP2K7'])
    print('BRAF__MAP2K2 N2_to_N1__0.75 distance:', combined_pw_distance_lists[m][a]['BRAF__MAP2K2'])

    global int_nodes_sorted
    interactome_path = '../../data'
    interactome_file = str(threshold) + 'thr_filtered_' + interactome_source + '_scored_interactome.txt'
    interactome_df = pd.read_csv(os.path.join(interactome_path, interactome_file), sep = '\t', na_filter = False)
    interactome = nx.from_pandas_edgelist(interactome_df, 'protein1', 'protein2', edge_attr = 'cost')
    int_nodes = list(interactome.nodes())
    int_nodes_sorted = sorted(int_nodes)
    int_node_path = os.path.join(interactome_path, str(threshold) + 'thr_filtered_' + interactome_source + '_nodes.txt')
    #with open(int_node_path, 'w') as file:
    #    for node in int_nodes_sorted:
    #       file.write(str(node) + '\n')

    cl_to_opt_targ_pair_dict = defaultdict(list)
    for cl in opt_targ_pair_dict.keys():
        for (combo_targset, dom_targ_pair, dom_targ_pair_dist) in opt_targ_pair_dict[cl][m][a]:
            syn = targset_to_syn_dict[cl][combo_targset]
            if combo_targset == 'BRAF__MAP2K1_MAP2K2_MAP2K3_MAP2K4_MAP2K5_MAP2K6_MAP2K7':
                print(cl, dom_targ_pair, dom_targ_pair_dist)
                cl_to_opt_targ_pair_dict[cl].append((dom_targ_pair, dom_targ_pair_dist, syn))
            if combo_targset == 'AKT1_AKT2__BCL2L1':
                print(cl, dom_targ_pair, dom_targ_pair_dist)
                cl_to_opt_targ_pair_dict[cl].append((dom_targ_pair, dom_targ_pair_dist, syn))
            if combo_targset == 'ATR__CHEK1':
                print(cl, dom_targ_pair, dom_targ_pair_dist)
                cl_to_opt_targ_pair_dict[cl].append((dom_targ_pair, dom_targ_pair_dist, syn))

    dp_nbhd_overlap_dict = dict()
    #combo_titles = ['BRAF__MAP2K2', 'BRAF__MAP2K7']
    #combo_titles = ['BRAF__MAP2K1', 'BRAF__MAP2K3']
    #combo_titles = ['BCL2L1__AKT1', 'BCL2L1__AKT2']
    #combo_titles = ['CHEK1__ATR']
    combo_titles = ['CHEK1__ATR', 'BRAF__MAP2K2', 'BRAF__MAP2K7', 'BCL2L1__AKT1', 'BCL2L1__AKT2']
    #combo_titles = ['CHEK1__ATR', 'BRAF__MAP2K1', 'BRAF__MAP2K3', 'BCL2L1__AKT1', 'BCL2L1__AKT2']
    for targ_pair in combo_titles: #['BRAF__MAP2K2', 'BRAF__MAP2K7']: #['BCL2L1__AKT1', 'BCL2L1__AKT2']: #['BRAF__MAP2K2', 'BRAF__MAP2K7']: #['CHEK1__ATR', 'BRAF__MAP2K2', 'BRAF__MAP2K7', 'BCL2L1__AKT1', 'BCL2L1__AKT2']:
        T1, T2 = targ_pair.split('__')
        dp_nbhd_overlap, T1_nbhd, T2_local_pert = get_dp_nbhd_overlap(targ_pair)
        targ_pair_title = '_'.join([targ_pair, 'overlap'])
        T1_nbhd_title = '_'.join([T1, 'nbhd'])
        T2_local_pert_title = '_'.join([T2, 'local_pert'])
        dp_nbhd_overlap_dict[targ_pair_title] = dp_nbhd_overlap
        dp_nbhd_overlap_dict[T1_nbhd_title] = T1_nbhd
        dp_nbhd_overlap_dict[T2_local_pert_title] = T2_local_pert

    #Saves protein lists for external GOEA
    #print(dp_nbhd_overlap_dict)
    #pickle.dump(dp_nbhd_overlap_dict, open(os.path.join(functionality_results_dir, 'T2_to_N1_0.75_example_combos_local_pert_nbhd_overlap_dict.pkl'), 'wb'))


    #INITIALIZE GOEA
    global obodag
    fin_gene2go = download_ncbi_associations() #download associations
    obodag = GODag("../../resources/go-basic.obo") #load ontologies

    global godag
    #Not sure if this is the same as the line above...
    godag = get_godag("../../resources/go-basic.obo", optional_attrs={'relationship'})

    # Read NCBI's gene2go. Store annotations in a list of namedtuples
    objanno = Gene2GoReader(fin_gene2go, taxids=[9606]) #homo sapien taxid = 9606
    ns2assoc = objanno.get_ns2assc() #Dictionary with MF, CC, BP keys and gene -> GO term annotations ex: 109617024: {'GO:0006396'}
    
    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated homo sapien genes". format(NS = nspc, N = len(id2gos)))

    global termcounts
    termcounts = TermCounts(obodag, ns2assoc['BP'])

    #Load background genes (PathFX 0.5 threshold interactome)
    sym2id_df = pd.read_csv('../../data/0.5thr_PathFX_int_GeneSymbol2GeneID.csv', sep = ',')

    #Create dictionary that maps gene ID to symbol in order to convert target neighborhoods to gene IDs
    global id2sym_dict
    id2sym_dict = dict()
    for i, row in sym2id_df.iterrows():
        if row['converted_alias'] != 'None':
            id = float(row['converted_alias'])
            sym = row['initial_alias']
            id2sym_dict[id] = sym
    
    global goeaobj
    goeaobj = GOEnrichmentStudyNS(
        id2sym_dict.keys(), # List of PathFX protein-coding genes with Gene IDs
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    
    #run_make_goea_tables(dp_nbhd_overlap_dict)
    
    
    global goea_dict
    goea_dict = defaultdict(list)
    depth_cutoff = 0
    for targ_pair in dp_nbhd_overlap_dict.keys():
        print(targ_pair)
        protein_list = dp_nbhd_overlap_dict[targ_pair]
        title = targ_pair
        run_goea(title, protein_list, depth_cutoff)
        print(title, goea_dict[title])

   

    #Are there shared GO terms between protein lists from the same drug combo?

    #combo_titles = ['BRAF', 'MAP2K2', 'MAP2K7', 'BRAF__MAP2K2', 'BRAF__MAP2K7']
    #combo_titles = ['BRAF', 'MAP2K1', 'MAP2K3', 'BRAF__MAP2K1', 'BRAF__MAP2K3']
    #combo_titles = ['BCL2L1', 'AKT1', 'AKT2', 'BCL2L1__AKT1', 'BCL2L1__AKT2']
    combo_titles = ['CHEK1', 'ATR', 'CHEK1__ATR']
    top_terms_dict = dict()
    top_terms_per_combo = set()
    for title in combo_titles: #['BCL2L1__AKT1', 'BCL2L1__AKT2', 'BCL2L1', 'AKT1', 'AKT2']: #['BRAF__MAP2K2', 'BRAF__MAP2K7', 'BRAF', 'MAP2K2', 'MAP2K7']:
        if len(goea_dict[title]) > 3:
            if '__' in title:
                term_limit = 5
            else:
                term_limit = 3
            top_terms = [(go, name, pval) for (go, name, pval) in goea_dict[title] if obodag[go].depth >= 7 and pval > 1e-30][:term_limit]
        
        elif len(goea_dict[title]) <= 3 and len(goea_dict[title]) > 0:
            top_terms = [(go, name, pval) for (go, name, pval) in goea_dict[title]]
        else:
            top_terms = []    

        top_terms_dict[title] = top_terms
        for (go, name, pval) in top_terms:
            top_terms_per_combo.add(' '.join([go, name]))
    
    top_terms_to_pval_dict = defaultdict(dd_float)
    for top_term in top_terms_per_combo:
        for title in combo_titles:
            go_term_to_pval = {' '.join([go, term]): pval for (go, term, pval) in goea_dict[title]}
            if top_term in go_term_to_pval.keys():
                top_terms_to_pval_dict[title][top_term] = go_term_to_pval[top_term]
            else:
                top_terms_to_pval_dict[title][top_term] = 1

    print(len(top_terms_per_combo))
    print(print(top_terms_dict))

    global combo_to_hm_color_dict
    combo_to_hm_color_dict = {'BRAF__MAP2K': plt.cm.Greens,
                          'AKT__BCL2L1': plt.cm.Purples,
                          'ATR__CHEK1': plt.cm.YlOrRd}
    
    #Plots GO BP enrichment for select terms for a given combination
    goea_hm(top_terms_to_pval_dict, 'ATR__CHEK1')
      

    
    #For dot plot (Figure 6A)
    num_terms_ranges = [(0, 0),
                    (1, 10),
                    (11, 20),
                    (21, 30),
                    (31, 40),
                    (41, 50),
                    (51, 60), 
                    (61, 70)]
    targ_pair_to_marker_size_dict = dict()
    for title in goea_dict.keys():
        if '__' in title:
            T1, T2 = title.split('__')
            targ_pair = title
            num_terms = len(goea_dict[title])
            for (i, r) in enumerate(num_terms_ranges):
                if num_terms >= r[0] and num_terms <= r[1]:
                    print(targ_pair, num_terms, r)
                    marker_size = 40 * 1/2 * (i + 1)
                    targ_pair_to_marker_size_dict['__'.join(sorted([T1, T2]))] = marker_size
    
    global targ_pair_to_color_dict
    targ_pair_to_color_dict = {'ATR__CHEK1': 'orange',
                            'AKT1__BCL2L1': 'purple',
                            'AKT2__BCL2L1': 'purple',
                            'BRAF__MAP2K1': 'forestgreen',
                            'BRAF__MAP2K2': 'forestgreen',
                            'BRAF__MAP2K3': 'forestgreen',
                            'BRAF__MAP2K7': 'forestgreen'}

    global cl_to_edge_color_dict
    cl_to_edge_color_dict = {'BT-474': 'darkslategrey', 
                            'CAMA-1': 'teal',
                            'MCF7': 'dodgerblue', 
                            'MDA-MB-157': 'darkturquoise',
                            'BT-20': 'lightslategrey'}
    dist_syn_plot_by_num_go_terms(cl_to_opt_targ_pair_dict, targ_pair_to_marker_size_dict)



if __name__ == "__main__":
    from scipy.cluster import hierarchy
    import pickle, os, functools, math
    from collections import defaultdict, Counter
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import networkx as nx
    from scipy.stats import ttest_ind, pearsonr, spearmanr, rankdata
    from scipy.spatial.distance import pdist

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

    #For other semantic similarities
    from goatools.semantic import *

    main()