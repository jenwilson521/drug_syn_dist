def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def find_robust_metrics(pre_params_dict):

    """
        Extracts metrics per cell line that are robust across curve fitting randomizations

        :param pre_params_dict: nested dictionary
            [cell line][NI + dir (e.g., T1_to_T2)][alphaD] = [(pre-fit slope, pre-fit intercept), ...]

        :return robust_metrics_per_cl: Dictionary of robust metrics
            [cell line] = [T1_to_T2__0.06, ....]
    """
    
    robust_metrics_per_cl = defaultdict(list)

    for cl in pre_params_dict.keys():
        m_all = list(pre_params_dict[cl].keys())
        m_all.remove('SP')
        for m in m_all:
            for a in pre_params_dict[cl][m].keys():
                pre_slope_list = list(zip(*pre_params_dict[cl][m][a]))[0]

                num_neg_slope = len([s for s in pre_slope_list[:100] if s <= 0])
                #print(cl, m, a, num_neg_slope)

                if num_neg_slope >= 75:
                    robust_metrics_per_cl[cl].append('__'.join([m, str(a)]))

    return robust_metrics_per_cl

def dependent_corr_ztest(xy, xz, yz, n, twotailed = True):

    """
        Computes Steiger's asymptotic z-test to compare dependent correlation coefficients
        
        :param xy: float
            Distance 1 (y) - synergy (x) correlation coefficient
        
        :param xz: float
            Distance 1 (z) - synergy (x) correlation coefficient
        
        :param yz: float
            Distance 1 (y) - distance 2 (z) correlation coefficient
        
        :param n: float
            Number of samples used to calculate correlation coefficients

        :param twotailed: boolean
            Indicator for two- vs. one-tailed z-test
            Default: True
        
        :return z: float
            Z-statistic
        
        :return p: float
            P-value from asymptotic z-test

    """
    #Fisher transformation
    xy_trans = np.arctanh(xy)
    xz_trans = np.arctanh(xz)

    #Asymptotic covariance calculations (s)
    rbar = (xy + xz) / 2 #avg. distance-synergy correlations
    rbar_sq = rbar * rbar

    s_num = yz * (1 - rbar_sq - rbar_sq) - 0.5 * rbar_sq * (1 - rbar_sq - rbar_sq - (yz * yz))
    s_den = (1 - rbar_sq) * (1 - rbar_sq)

    #Z-statistic calculations
    z_num = np.sqrt((n - 3)) * (xy_trans - xz_trans)
    z_den = np.sqrt((2 - 2 * (s_num / s_den)))

    z = z_num / z_den

    p = 1 - norm.cdf(abs(z))

    if twotailed:
        p *= 2
    
    return z, p

def diffusion_depth_tests(alphas_to_compare):
    """
        Creates correlation line plots, sorted by increasing alphaD with NI and directionality held constant in each plot
        Runs and saves significance of all alphaD sensitivity tests, calculated with a dependent correlation test

        :param alphas_to_compare: List of tuples
            [(alpha 1, alpha 2, ...)]
        
        :return diffusion_depth_dict: Nested dictionary
            [cell line][alpha 1__alpha 2__metric] = p-value
    """
    x_all = [i + 0.5 for i, _ in enumerate(alphas_to_compare)]

    m_all = list(post_pcc_dict['BT-20'].keys()) #[T1_to_T2, T2_to_T1, T1_to_N2, T2_to_N1, N1_to_N2, N2_to_N1]
    m_all.remove('SP') #ASP not considered in sensitivity analysis

    diffusion_depth_dict = defaultdict(dd_float)
    for m in m_all:

        fig = plt.figure(dpi = 300)
        fig.tight_layout()
        ax = plt.subplot()

        for cl in post_pcc_dict.keys():

            #Collecting robust correlations to be plotted
            x_per_cl = list()
            pcc_per_cl = list()
            for (x, a) in zip(x_all, alphas_to_compare):
                if '__'.join([m, str(a)]) in robust_metrics_per_cl[cl]:
                    x_per_cl.append(x)
                    pcc_per_cl.append(post_pcc_dict[cl][m][a][2])
                else:
                    continue
            
            plt.plot(x_per_cl, pcc_per_cl, marker = 'o', color = cl_to_color_dict[cl], 
                     linewidth = 2.5, markersize = 15, label = cl)
            
            #Correlation comparison tests
            compared = []
            for a1 in alphas_to_compare:
                for a2 in alphas_to_compare:
                    if a1 != a2:
                        if '__'.join([m, str(a1)]) in robust_metrics_per_cl[cl] and '__'.join([m, str(a2)]) in robust_metrics_per_cl[cl]:

                                a1__a2 = '__'.join(sorted([str(a1), str(a2)]))

                                if a1__a2 not in compared:
                                    compared.append(a1__a2)
                                    a1_dist = list(zip(*opt_targ_pair_dict[cl][m][a1]))[2]
                                    a2_dist = list(zip(*opt_targ_pair_dict[cl][m][a2]))[2]

                                    a1_a2_corr, _ = pearsonr(a1_dist, a2_dist) #Distance 1 - distance 2 corr

                                    a1_syn_corr = abs(post_pcc_dict[cl][m][a1][0][0])
                                    a2_syn_corr = abs(post_pcc_dict[cl][m][a2][0][0])

                                    z, p = dependent_corr_ztest(a1_syn_corr, a2_syn_corr, a1_a2_corr, len(a1_dist))

                                    diffusion_depth_dict[cl]['__'.join([str(a1), str(a2), m])] = p

        #Figure setting specifications
        for axis in ['left']:
                ax.spines[axis].set_linewidth(1.5)
        for axis in ['top','right', 'bottom']:
                ax.spines[axis].set_linewidth(0)                            

        plt.grid()
        ax.tick_params(axis = 'y', width = 2)
        ax.tick_params(axis = 'x', width = 0)
        plt.setp(ax, xlim = [0, 3])
        plt.xticks([0.5, 1.5, 2.5], labels = [])
        plt.yticks(np.arange(-0.7, -0.2, 0.1), labels = [])

        figure_name = m + '_diffusion_depth.png'
        plt.savefig(os.path.join(figure_path, 'Diffusion_Depth', figure_name), bbox_inches = 'tight', pad_inches = 0.1)            

        return diffusion_depth_dict                

def nbhd_tests(nbhds_to_compare, alphas_to_compare):
    """
        Creates correlation line plots, sorted by increasing amounts of neighborhood info with alphaD and directionality held constant in each plot
        Runs and saves significance of all NI sensitivity tests, calculated with a dependent correlation test

        :param nbhds_to_compare: List of tuples
            [(NI 1, N1 2, NI 3), ...]
        
        param alphas_to_compare: List of tuples
            [(alpha 1, alpha 2, ...)]
            Included in order to restrict comparisons to only include alphaDs from diffusion depth tests
        
        :return diffusion_depth_dict: Nested dictionary
            [cell line][metric 1__metric 2__alphaD] = p-value
    """
    x_all = [i + 0.5 for i, _ in enumerate(nbhds_to_compare)]


    nbhd_inclusion_dict = defaultdict(dd_float)
    for i, dir in enumerate(['12', '21']):
        for a in alphas_to_compare[0]:

            fig = plt.figure(dpi = 300)
            fig.tight_layout()
            ax = plt.subplot()

            for cl in post_pcc_dict.keys():

                #Collecting robust correlations to be plotted
                x_per_cl = list()
                pcc_per_cl = list()
                for (x, m) in zip(x_all, nbhds_to_compare[i]):
                    if '__'.join([m, str(a)]) in robust_metrics_per_cl[cl]:
                        x_per_cl.append(x)
                        pcc_per_cl.append(post_pcc_dict[cl][m][a][2])
                    else:
                        continue
                
                plt.plot(x_per_cl, pcc_per_cl, marker = 'o', color = cl_to_color_dict[cl], 
                         linewidth = 2.5, markersize = 15, label = cl)
                
                #Correlation comparison tests
                compared = []
                for m1 in nbhds_to_compare[i]:
                    for m2 in nbhds_to_compare[i]:
                        if m1 != m2:

                            if '__'.join([m1, str(a)]) in robust_metrics_per_cl[cl] and '__'.join([m2, str(a)]) in robust_metrics_per_cl[cl]:

                                m1__m2 = '__'.join(sorted([m1, m2]))

                                if m1__m2 not in compared:
                                    compared.append(m1__m2)
                                    m1_dist = list(zip(*opt_targ_pair_dict[cl][m1][a]))[2]
                                    m2_dist = list(zip(*opt_targ_pair_dict[cl][m2][a]))[2]

                                    m1_m2_corr, _ = pearsonr(m1_dist, m2_dist)

                                    m1_syn_corr = abs(post_pcc_dict[cl][m1][a][0][0])
                                    m2_syn_corr = abs(post_pcc_dict[cl][m2][a][0][0])

                                    t, p = dependent_corr_ztest(m1_syn_corr, m2_syn_corr, m1_m2_corr, len(m1_dist))

                                    nbhd_inclusion_dict[cl]['__'.join([m1, m2, str(a)])] = p

            #Figure setting specifications
            for axis in ['left']:
                ax.spines[axis].set_linewidth(1.5)
            for axis in ['top','right', 'bottom']:
                ax.spines[axis].set_linewidth(0)
        
            #plt.legend()
            plt.grid()
            ax.tick_params(axis = 'y', width = 2)
            ax.tick_params(axis = 'x', width = 0)
            plt.setp(ax, xlim = [0, 3])
            plt.xticks([0.5, 1.5, 2.5], labels = [])
            plt.yticks(np.arange(-0.7, -0.2, 0.1), labels = [])


            figure_name = dir + '_' + str(a) + 'aD_nbhd_inclusion.png'
            plt.savefig(os.path.join(figure_path, 'Neighborhood_Inclusion', figure_name), bbox_inches = 'tight', pad_inches = 0.1)
    
    return nbhd_inclusion_dict

def directionality_tests(dir_to_compare, alphas_to_compare):

    """
        Creates correlation line plots, sorted by 1->2 then 2->1 directionality with alphaD and NI held constant in each plot
        Runs and saves significance of all directionality sensitivity tests, calculated with a dependent correlation test

        :param dir_to_compare: List of tuples
            [(1_to_2, 2_to_1), ...]
        
        param alphas_to_compare: List of tuples
            [(alpha 1, alpha 2, ...)]
            Included in order to restrict comparisons to only include alphaDs from diffusion depth tests
        
        :return directionality_dict: Nested dictionary
            [cell line][metric 1__metric 2__alphaD] = p-value
    """
    x_all = [i + 0.5 for i, _ in enumerate(dir_to_compare)]

    directionality_dict = defaultdict(dd_float)
    for i, nbhd in enumerate(['TT', 'TN', 'NN']):
        for a in alphas_to_compare[0]:

            fig = plt.figure(dpi = 300)
            fig.tight_layout()
            ax = plt.subplot()

            for cl in post_pcc_dict.keys():

                #Collecting robust correlations to be plotted
                x_per_cl = list()
                pcc_per_cl = list()
                for (x, m) in zip(x_all, dir_to_compare[i]):
                    if '__'.join([m, str(a)]) in robust_metrics_per_cl[cl]:
                        x_per_cl.append(x)
                        pcc_per_cl.append(post_pcc_dict[cl][m][a][2])
                    else:
                        continue
                
                plt.plot(x_per_cl, pcc_per_cl, marker = 'o', color = cl_to_color_dict[cl], 
                         linewidth = 2.5, markersize = 15, label = cl)

                #Correlation comparison tests
                compared = []
                for m1 in dir_to_compare[i]:
                    for m2 in dir_to_compare[i]:
                        if m1 != m2:

                            m1__m2 = '__'.join(sorted([m1, m2]))

                            if '__'.join([m1, str(a)]) in robust_metrics_per_cl[cl] and '__'.join([m2, str(a)]) in robust_metrics_per_cl[cl]:


                                if m1__m2 not in compared:
                                    compared.append(m1__m2)

                                    m1_dist = list(zip(*opt_targ_pair_dict[cl][m1][a]))[2]
                                    m2_dist = list(zip(*opt_targ_pair_dict[cl][m2][a]))[2]

                                    m1_m2_corr, _ = pearsonr(m1_dist, m2_dist)

                                    m1_syn_corr = abs(post_pcc_dict[cl][m1][a][0][0])
                                    m2_syn_corr = abs(post_pcc_dict[cl][m2][a][0][0])

                                    t, p  = dependent_corr_ztest(m1_syn_corr, m2_syn_corr, m1_m2_corr, len(m1_dist))

                                    directionality_dict[cl]['__'.join([m1, m2, str(a)])] = p

            #Figure setting specifications
            for axis in ['left']:
                ax.spines[axis].set_linewidth(1.5)
            for axis in ['top','right', 'bottom']:
                ax.spines[axis].set_linewidth(0)
        
            plt.grid()
            ax.tick_params(axis = 'y', width = 2)
            ax.tick_params(axis = 'x', width = 0)
            plt.setp(ax, xlim = [0, 2])
            plt.xticks([0.5, 1.5], labels = [])
            plt.yticks(np.arange(-0.7, -0.2, 0.1), labels = [])

            figure_name = nbhd + '_' + str(a) + 'aD_directionality.png'
            plt.savefig(os.path.join(figure_path, 'Directionality', figure_name), bbox_inches = 'tight', pad_inches = 0.1)

    return directionality_dict
    

def main():

    #Global variables
    global figure_path
    global post_pcc_dict
    global robust_metrics_per_cl
    global cl_to_color_dict
    global opt_targ_pair_dict

    cl_to_color_dict = {'BT-474': 'navy', 
                    'CAMA-1': 'blue',
                    'MCF7': 'dodgerblue', 
                    'MDA-MB-157': 'lightblue',
                    'BT-20': 'indigo'}

    #Path where figures are saved
    figure_path = '../../figures/sensitivity_analysis'
    check_dir(figure_path)

    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'
    exper_source = 'DREAM'

    #Load curve fitting results
    curve_fitting_results_path = os.path.join('../../results/Curve_Fitting', exper_source, interactome_aname)
    analysis = 'posyn'
    pre_params_filename = 'pre_fit_params_dict_' + analysis + '.pkl'
    post_pcc_filename = 'post_fit_pcc_and_sig_dict_' + analysis + '.pkl'

    post_pcc_dict = pickle.load(open(os.path.join(curve_fitting_results_path, post_pcc_filename), 'rb')) #[cl][m][a] = [(corr,p-val - merged nbhd), (corr,p-val - min), (corr,p-val - max), ... n=98 rand ]
    pre_params_dict = pickle.load(open(os.path.join(curve_fitting_results_path, pre_params_filename), 'rb')) #[cl][m][a] = [(b,m - merged nbhd), (b,m - min), (b,m - max), ... n=98 rand ] 
    opt_targ_pair_dict = pickle.load(open(os.path.join(curve_fitting_results_path, 'opt_targ_pair_dict_5_cl_37_metrics.pkl'), 'rb')) #[cl][m][a] = [(combo target set, opt target pair, opt distance), ...]

    #Find metrics robust to starting point randomizations in curve fitting framework
    robust_metrics_per_cl = find_robust_metrics(pre_params_dict)

    global dd_float
    dd_float = functools.partial(defaultdict, float)
    all_feature_dict = dict()

    #Diffusion depth (alphaD) sensitivity tests and plots
    alphas_to_compare = [0.06, 0.30, 0.75]
    diffusion_depth_dict = diffusion_depth_tests(alphas_to_compare)
    all_feature_dict['Diffusion depth'] = diffusion_depth_dict

    nbhds_to_compare = [('T1_to_T2', 'T1_to_N2', 'N1_to_N2'), ('T2_to_T1', 'T2_to_N1', 'N2_to_N1')]
    nbhd_inclusion_dict = nbhd_tests(nbhds_to_compare, alphas_to_compare)
    all_feature_dict['Neighborhood inclusion'] = nbhd_inclusion_dict

    dir_to_compare = [('T1_to_T2', 'T2_to_T1'), ('T1_to_N2', 'T2_to_N1'), ('N1_to_N2', 'N2_to_N1')]
    directionality_dict = nbhd_tests(dir_to_compare, alphas_to_compare)
    all_feature_dict['Directionality'] = directionality_dict

    #Saving results
    sens_analysis_results_path = os.path.join('../../results/Sensitivity_Analysis', exper_source, interactome_aname)
    check_dir(sens_analysis_results_path)
    outpath = os.path.join(sens_analysis_results_path, 'sens_analysis_sig_all_features_dict.pkl')
    pickle.dump(all_feature_dict, open(outpath, 'wb'))


if __name__ == "__main__":

    import os, pickle, functools
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr, spearmanr
    from scipy.stats import norm
    from collections import defaultdict

    main()