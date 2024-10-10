def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def calc_performance_metric(pearson_dict, spearman_dict):
    """
        Calculates metric for evaluating the performance of each distance metric
		in the curve fitting framework
		
		metric = avg([-log(sd_PCC) * mean_PCC * merged_nbhd_PCC, 
		        -log(sd_spearman) * mean_spearman * merged_nbhd_spearman])
		
		:param pearson_dict: nested dictionary
            Pearson's correlation coefficients from curve fitting for 
			merged neighborhood (MN) and 100 random starting points
            [cell line][metric][alphaD] = [(MN_PCC, PCC_sig), (rand_PCC, PCC_sig), ...]

		:param spearman_dict: nested dictionary
            Spearman's correlation coefficients from curve fitting for 
			merged neighborhood (MN) and 100 random starting points
            [cell line][metric][alphaD] = [(MN_spearman, spearman_sig), (rand_spearman, spearman_sig), ...]
		
		:return perform_dict: nested dictionary
            Performance metric results for all cell lines and distance metrics
			[cell line][metric__alphaD] = performance metric
    """

    perform_dict = defaultdict(dd_float)
    for cl in pearson_dict.keys(): #['BT-20', 'BT-474', 'CAMA-1', 'MCF7', 'MDA-MB-157']
        for m in pearson_dict[cl].keys(): #[T1_to_T2, T2_to_T1, T1_to_N2, T2_to_N1, N1_to_N2, N2_to_N1]
            for a in pearson_dict[cl][m].keys(): #[0.06, 0.15, 0.30, 0.45, 0.60, 0.75]

                merged_nbhd_pcc = pearson_dict[cl][m][a][0][0]
                avg_pcc = np.mean(list(zip(*pearson_dict[cl][m][a][1:]))[0])
                std_pcc = np.std(list(zip(*pearson_dict[cl][m][a][1:]))[0])
                pcc_perform = merged_nbhd_pcc * avg_pcc * -np.log(std_pcc)

                merged_nbhd_spearman = spearman_dict[cl][m][a][0][0]
                avg_spearman = np.mean(list(zip(*spearman_dict[cl][m][a][1:]))[0])
                std_spearman = np.std(list(zip(*spearman_dict[cl][m][a][1:]))[0])
                spearman_perform = merged_nbhd_spearman * avg_spearman * -np.log(std_spearman)

                perform_dict[cl]['__'.join([m, str(a)])] = np.mean([pcc_perform, spearman_perform])

    return perform_dict

def count_top_metric_features(top_metrics):
    """
        Decomposes distance metric names (saved as metric__alphaD) into their three network feature assumptions

        :param top_metrics: list
            Distance metrics that surpass user-specified performance threshold
            [m__a, m__a, ...] 
                where m = [T1_to_T2, T2_to_T1, T1_to_N2, T2_to_N1, N1_to_N2, N2_to_N1]
                      a = [0.06, 0.15, 0.30, 0.45, 0.60, 0.75]
        
        :return feature_count_dict: nested dictionary
            Occurence of each feature among metrics in top_metrics
    
    """

    #Initialize counts for all assumptions for network features (diffusion depth, nbhd inclusion, directionality)
    low_alphaD = 0
    mid_alphaD = 0
    high_alphaD = 0

    tt = 0
    tn = 0
    nn = 0

    one_to_two = 0
    two_to_one = 0

    feature_count_dict = defaultdict(dd_float)
    #Total number of times each assumptions occurs in all 37 metrics across 5 cell lines (n = 185)
    feature_count_dict['Low alphaD']['Total'] = 12 * 5 #alphaD = 0.06, 0.15
    feature_count_dict['Mid alphaD']['Total'] = 12 * 5 #alphaD = 0.30, 0.45
    feature_count_dict['High alphaD']['Total'] = 12 * 5 #alphaD = 0.60, 0.75

    feature_count_dict['TT']['Total'] = 12 * 5 #No nbhd info
    feature_count_dict['TN']['Total'] = 12 * 5 #Partial nbhd info
    feature_count_dict['NN']['Total'] = 12 * 5 #All nbhd info

    feature_count_dict['1->2']['Total'] = 18 * 5 #Large to small nbhd
    feature_count_dict['2->1']['Total'] = 18 * 5 #Small to large nbhd

    for m__a in top_metrics:
        m__a_split = m__a.split('_')

        if m__a_split[4] in ['0.06', '0.15']:
            low_alphaD += 1
        elif m__a_split[4] in ['0.3', '0.45']:      
            mid_alphaD += 1
        elif m__a_split[4] in ['0.6', '0.75']:
            high_alphaD += 1
            
        if 'T' in m__a_split[0] and 'T' in m__a_split[2]:
            tt += 1
        elif 'T' in m__a_split[0] and 'N' in m__a_split[2]:
            tn += 1
        elif 'N' in m__a_split[0] and 'N' in m__a_split[2]:
            nn += 1

        if '1' in m__a_split[0]:
            one_to_two += 1
        elif '2' in m__a_split[0]:
            two_to_one += 1
    
    feature_count_dict['Low alphaD']['Top'] = low_alphaD
    feature_count_dict['Mid alphaD']['Top'] = mid_alphaD
    feature_count_dict['High alphaD']['Top'] = high_alphaD

    feature_count_dict['TT']['Top'] = tt
    feature_count_dict['TN']['Top'] = tn
    feature_count_dict['NN']['Top'] = nn

    feature_count_dict['1->2']['Top'] = one_to_two
    feature_count_dict['2->1']['Top'] = two_to_one

    return feature_count_dict

def feature_assumption_enrichment(feature_count_dict, M, n):
    """
        Performs enrichment analysis for all feature assumptions among the top performing metrics

        :param feature_count_dict: nested dictionary
            [Feature assumption][Total] = # of occurences of feat. assump. among all 185 curves 
                                          N in hypergeometric test
            [Feature assumption][Top] = # of occurences of feat. assump. among top performing curves 
                                          k in hypergeometric test
        
        :param M: float
            Total number of curves (M in hypergeometric test)

        :param n: float
            Total number of top performing curves (n in hypergeometric test)                       

        :return enrich_dict: dictionary
            [Feat. assump] = enrichment p-value
    """
    enrich_dict = dict()
    for feat_assump in feature_count_dict.keys():
        N = feature_count_dict[feat_assump]['Total']
        k = feature_count_dict[feat_assump]['Top']

        pval = 1 - hypergeom.cdf(k, M, n, N)
        enrich_dict[feat_assump] = pval
    
    return enrich_dict

def plot_performance_metric(perform_dict_norm):
    """
        Plots sorted bar charts, seperated by cell line

        :param perform_dict_norm: nested dictionary
            Performance metric results for all cell lines and distance metrics
			[cell line][metric__alphaD] = normalized performance metric

    """

    for cl in perform_dict_norm.keys():

        fig = plt.figure(dpi = 300)
        ax = plt.subplot()

        metric_perform_sorted = sorted(perform_dict_norm[cl].items(), key = lambda x:x[1], reverse = True)

        perform_sorted = list(zip(*metric_perform_sorted))[1]
        metric_sorted = list(zip(*metric_perform_sorted))[0]

        x = len(perform_sorted)
        plt.bar(x, perform_sorted, width = 0.75, color = 'navy', alpha = 0.75)

        #Figure setting specifications
        for axis in ['left']:
            ax.spines[axis].set_linewidth(2)
        for axis in ['top','right', 'bottom']:
            ax.spines[axis].set_linewidth(0)
        ax.set_ylim([-0.1, 1])
        plt.xticks(np.arange(len(metric_sorted)), labels = metric_sorted, fontsize = 5, rotation = 270)
        plt.yticks(np.arange(0, 1.25, 0.25), labels = [], fontsize = 8)

        figure_name = cl + '_PCC_Spearman_performance_metric_bars.png'
        plt.savefig(os.path.join(figure_path, figure_name))


def main():
    global dd_float
    dd_float = functools.partial(defaultdict, float)

    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'
    exper_source = 'DREAM'
    model_name = 'posyn'

    #Load curve fitting results
    curve_fitting_results_dir = os.path.join('../../results/Multi_Targ_Curve_Fitting', exper_source, interactome_aname)
    check_dir(curve_fitting_results_dir)

    post_pcc_dict = pickle.load(open(os.path.join(curve_fitting_results_dir, 'post_fit_pcc_and_sig_dict_' + model_name + '.pkl'), 'rb'))
    post_spearman_dict = pickle.load(open(os.path.join(curve_fitting_results_dir, 'post_fit_spearman_and_sig_dict_' + model_name + '.pkl'), 'rb'))

    #Calculate performance score for 37 metrics across 5 cell lines (n = 185)
    perform_dict = calc_performance_metric(post_pcc_dict, post_spearman_dict)

    #Find max performance to normalize across cell lines
    max_perform = 0
    for cl in perform_dict.keys():
        for perform in perform_dict[cl].values():
            if perform > max_perform:
                max_perform = perform
    
    #Normalize performance scores
    perform_dict_norm = defaultdict(dd_float)
    for cl in perform_dict.keys():
        perform_dict_norm[cl] = {m__a: perform / max_perform for m__a, perform in perform_dict[cl].items()}

    #Find top performing distance metrics across all cell lines
    perform_thresh = 0.70
    top_metric_list = list()
    for cl in perform_dict_norm.keys():
        for m__a, perform_norm in perform_dict_norm.items():
            if perform_norm >= perform_thresh:
                top_metric_list.append(m__a)

    #Count number of occurences of each network feature assumption for top performing curves
    feature_count_dict = count_top_metric_features(top_metric_list)         

    #Calculate enrichment p-value for each network feature assumption
    enrich_dict = feature_assumption_enrichment(feature_count_dict, 185, len(top_metric_list))

    file_name = 'top_metric_feature_enrichment_dict.pkl'
    pickle.dump(enrich_dict, open(os.path.join(curve_fitting_results_dir, file_name), 'wb'))

    global figure_path
    figure_path = '../../figures/multi_targ_curve_fitting'
    check_dir(figure_path)

    #Plot metric performance for all cell lines
    plot_performance_metric(perform_dict_norm)


if __name__ == "__main__":
    import pickle, os, functools, math
    import numpy as np
    import pandas as pd
    from collections import defaultdict
    from scipy.stats import iqr, hypergeom
    import matplotlib.pyplot as plt
    main()