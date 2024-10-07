def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def get_all_combos(targset_list):
    """
        Generates all target pairs between a combo's 2 target sets

        :param targset_list: list
            List of targets to be paired within a given combo's target sets
        
        :retun set_combo: set
            Set of all target pairs for a given combo's target sets
    
    """
    set_combo = set()
    if len(targset_list) == 1:
        for t1 in targset_list[0]:
            for t2 in targset_list[0]:
                if t1 != t2:
                    combo = '__'.join(sorted((t1, t2)))
                    set_combo.add(combo)
    else:
        D1_targset = targset_list[0]
        D2_targset = targset_list[1]
        for t1 in D1_targset:
            for t2 in D2_targset:
                combo = '__'.join(sorted((t1, t2)))
                set_combo.add(combo)         

    return set_combo

def get_cell_line_combos(dist_dict):

    """
        Creates a dictionary with (pw_combo, distance) tuples, unique to the given cell line

        :param dist_dict: dictionary
            [target 1__target 2] = distance
        
        :return pw_combo_and_dist_dict: dictionary of lists
            Dictionary of all target pairs and their distances for all combos
            [target set 1__target set 2] = [(target 1__target 2, distance), ...]
    """

    pw_combo_and_dist_dict = defaultdict(list)

    for combo_targset in list(combo_targset_list):

        combo_targset_split = combo_targset.split('__')
        D1_targset = combo_targset_split[0]
        D2_targset = combo_targset_split[1]

        D1_targset_split = D1_targset.split('_')
        D2_targset_split = D2_targset.split('_')
        D1D2_targset_split = [D1_targset_split, D2_targset_split]

        D1_pairwise_combos = get_all_combos([D1_targset_split])
        D2_pairwise_combos = get_all_combos([D2_targset_split])
        D1D2_pairwise_combos = get_all_combos(D1D2_targset_split)

        syn = targset_to_syn_dict[cl][combo_targset]
        
        pw_combo_and_dist_dict[combo_targset] = [(pw_combo, dist_dict[pw_combo]) for pw_combo in D1D2_pairwise_combos]
    
    return pw_combo_and_dist_dict

def first_degree_poly(dist, a, b):
    return a * dist + b

def second_degree_poly(dist, a, b, c): 
    return (a * (dist**2)) + (b * dist) + c

def third_degree_poly(dist, a, b, c, d):
    return (a * (dist**3)) + (b * (dist**2)) + (c * dist) + d

def exp_func(dist, a, b, c):
    return a * np.exp(-b * dist) + c

def infer_syn(dist, a, b, c = 0, d = 0):
    """
        Uses fitted parameters (e.g., a and b) to calculate predicted synergy

        :param dist: list
            List of distances used to calculated predicted synergy, sorted according to combo_targset

        :param a: float
            Slope parameter for first-degree poly
            x^2 coefficient for second-degree poly
            x^3 coefficient for third-degree poly
            e^-x coefficient for exponential function
        
        :param b: float
            Intercept parameter for first-degree poly
            x coefficient for second-degree poly
            x^2 coefficient for third-degree poly
            e^-x exponent scaling factor for exponential function
        
        :param c: float
            Intercept parameter for second-degree poly
            x coefficient for third-degree poly
            Intercept parameter for exponential function
        
        :param d: float
            Intercept parameter for third-degree poly
        
        :pred_syn: list
            List of predicted synergy scores


    """
    pred_syn = [(d * a) + b for d in dist] #First-degree polynomial
    #pred_syn = [(a * (d**2)) + (b * d) + c for d in dist] #Second-degree polynomial
    #pred_syn = [(a * (di**3)) + (b * (di**2)) + (c * di) + d for di in dist] #Third-degree polynomial
    #pred_syn = [a * np.exp(-b * d) + c for d in dist] #Exponential function
    return pred_syn

def model_lossfunc(dist, exper_syn, a, b, c = 0, d = 0):

    """
        Predicts synergy according to curve parameters 
        and then calculates predicted vs. true synergy model loss

        :param dist: list
            List of distances used to calculate predicted synergy, sorted according to combo_targset
        
        :param exper_syn: list
            List of experimental synergy scores, sorted according to combo_targset

        :param a: float
            Slope parameter for first-degree poly
            x^2 coefficient for second-degree poly
            x^3 coefficient for third-degree poly
            e^-x coefficient for exponential function
        
        :param b: float
            Intercept parameter for first-degree poly
            x coefficient for second-degree poly
            x^2 coefficient for third-degree poly
            e^-x exponent scaling factor for exponential function
        
        :param c: float
            Intercept parameter for second-degree poly
            x coefficient for third-degree poly
            Intercept parameter for exponential function
        
        :param d: float
            Intercept parameter for third-degree poly

        :param model_loss: float
            Sum of squared errors between predicted and real synergies 
    """

    pred_syn = infer_syn(dist, a, b)

    model_loss = np.nansum([(p - e)**2 for p, e in zip(pred_syn, exper_syn)])

    return model_loss

def optimize_lossfunc(exper_syn, initial_pred, a, b, c = 0, d = 0, maxiter = 1000):

    """
        Updates curve parameters, distances, and predicted synergies in order to minimize model loss

        :param exper_syn: list
            List of experimental synergy scores, sorted according to combo_targset
        
        :param initial_pred: list
            Distances used to initialize optimization, sorted according to combo_targset

        :param a: float
            Slope parameter for first-degree poly
            x^2 coefficient for second-degree poly
            x^3 coefficient for third-degree poly
            e^-x coefficient for exponential function
        
        :param b: float
            Intercept parameter for first-degree poly
            x coefficient for second-degree poly
            x^2 coefficient for third-degree poly
            e^-x exponent scaling factor for exponential function
        
        :param c: float
            Intercept parameter for second-degree poly
            x coefficient for third-degree poly
            Intercept parameter for exponential function
        
        :param d: float
            Intercept parameter for third-degree poly
        
        :param max_iter: float
            Maximum number of iterations to be used in minimize function
        
        :return opt: object
            Optimization results, opt.x = distance solution array

    """
    
    #Minimization is restricted by each combo's min and max pairwise distance
    bnds = [(np.min(pw_dist), np.max(pw_dist)) for combo_targset, pw_dist in pw_dist_dict.items()]
                     
    #opt = least_squares(model_lossfunc, initial_pred, args=(exper_syn, a, b, c))
    #                #method = 'Powell') 

    opt = minimize(model_lossfunc, initial_pred, args=(exper_syn, a, b, c, d), method = 'Powell', bounds = bnds)    

    return opt

def fit_curve(dist_syn_list):

    """
        Uses selected curve type (e.g., first-degree polynomial) to generate curve parameters

        :param dist_syn_list: list of tuples
            [(distance, synergy),...]

        :return popt: array
            Optimal parameters from curve fitting
        
        :return pcov: 2-D array
            Covariance matrix of the optimal parameters from curve fitting
    """

    dist = list(zip(*dist_syn_list))[1]
    syn = list(zip(*dist_syn_list))[2]

    popt, pcov = curve_fit(first_degree_poly, dist, syn)
    #popt, pcov = curve_fit(second_degree_poly, dist, syn)
    #popt, pcov = curve_fit(third_degree_poly, dist, syn)
    #popt, pcov = curve_fit(exp_func, dist, syn)

    return popt, pcov

def find_opt_pw_dist(combo_targset, opt_dist):
    """
        Finds pairwise distance within combo's distance vector
        that is most similar to the optimal distance following fitting

        :param combo_targset: list
            [target set 1__target set 2, ...]
        
        :param opt_dist: list
            Length of and sorted according to combo_targset
        
        :return opt_pw_combo_dict: dictionary
            [target set 1__target set 2] = (target 1__target 2, distance)
    """

    dist_dif_dict = defaultdict(list)
    for i, cts in enumerate(combo_targset):

        #Calculating the difference between each pairwise distance and the optimal from fitting
        dist_dif_dict[cts] = [(pw_combo, pw_d, abs(pw_d - opt_dist[i])) 
                                for pw_combo, pw_d in pw_combo_and_dist_dict[cts]]

    #Finding and saving the target pair with the minimum distance
    opt_pw_combo_dict = dict()
    for cts in dist_dif_dict.keys():
        differences = list(zip(*dist_dif_dict[cts]))[2]
        pw_dist = list(zip(*dist_dif_dict[cts]))[1]
        pw_combos = list(zip(*dist_dif_dict[cts]))[0]
        
        opt_pw_ind = differences.index(min(differences))

        opt_pw_combo_dict[cts] = (pw_combos[opt_pw_ind], pw_dist[opt_pw_ind])

    return opt_pw_combo_dict

def pipeline(dist_syn_list):

    """

    :param dist_syn_list: list of tuples 
        [(target set 1__target set 2, dist, synergy), ...] for curve fitting initialization which is used to:
    
    
    (1) Find the initial line of best fit (fit_curve)
    (2) Optimize the line of best fit by finding the parameters and distances 
        that result in minimal model loss, restricted by each combo's distance vector (optimize_lossfunc)
    (3) Map the optimal distances from curve fitting to those within each combo's distance vector
        (find_opt_pw_dist)
    (4) Find final line of best fit using distances from step (3)
    (5) Calculate final synergy predictions using final line of best fit parameters (infer_syn)
    (6) Calculate performance metrics (Pearson, Spearman, RMSE) and saves results
    
    :return opt_pw_combo_dict: dictionary 
        [target set 1__target set 2] = (target 1__target 2, distance)
    """

    combo_targset = list(zip(*dist_syn_list))[0]
    input_dist = list(zip(*dist_syn_list))[1]
    exper_syn = list(zip(*dist_syn_list))[2]

    pre_popt, pcov = fit_curve(dist_syn_list) #pre_opt = curve parameters, pcov = covariance matrix

    #NOTE - update pre_opt parameters according to curve type (e.g., for second-degree poly, use pre_popt[0], pre_popt[2], pre_popt[3])
    opt = optimize_lossfunc(exper_syn, input_dist, pre_popt[0], pre_popt[1], maxiter = 1000) #opt = object output from scipy.optimize.minimize

    opt_pw_combo_dict = find_opt_pw_dist(combo_targset, opt.x) #

    opt_pw_dist = list(zip(*list(opt_pw_combo_dict.values())))[1]

    opt_popt, pcov = fit_curve(list(zip(combo_targset, opt_pw_dist, exper_syn)))

    #NOTE - update pre_opt parameters according to curve type (e.g., for second-degree poly, use opt_popt[0], opt_popt[2], opt_popt[3])
    final_syn_pred = infer_syn(opt_pw_dist, opt_popt[0], opt_popt[1])

    pcc_corr, pcc_pval = pearsonr(opt_pw_dist, exper_syn)
    spearman_corr, spearman_pval = spearmanr(opt_pw_dist, exper_syn)

    df = len(opt_pw_dist) - len(opt_popt) #degrees of freedom
    rmse = np.sqrt(np.sum((np.array(exper_syn) - np.array(final_syn_pred))**2) / df) #root mean squared error

    #Saving results according to analysis type
    if analysis == 'shuffle synergy':
        post_fit_pcc_shuffle[cl][metric][a].append((pcc_corr, pcc_pval))
        post_fit_spearman_shuffle[cl][metric][a].append((spearman_corr, spearman_pval))
        post_fit_rmse_shuffle[cl][metric][a].append(rmse)
    elif analysis == 'random starting point':
        rand_post_fit_pcc.append((pcc_corr, pcc_pval))
        rand_post_fit_spearman.append((spearman_corr, spearman_pval))
        rand_post_fit_rmse.append(rmse)
        rand_pre_fit_params.append(pre_popt)
        rand_cond_num.append(np.linalg.cond(pcov))
    else:
        post_fit_pcc[cl][metric][a].append((pcc_corr, pcc_pval))
        post_fit_spearman[cl][metric][a].append((spearman_corr, spearman_pval))
        post_fit_rmse[cl][metric][a].append(rmse)
        pre_fit_params[cl][metric][a].append(pre_popt)
        cond_num[cl][metric][a].append(np.linalg.cond(pcov))

    return opt_pw_combo_dict

def run_random_starting_point(num_cores):

    """
    :param num_cores = # of processing units

    :return: None
    
    """

    #Generating 98 lists of random target pair distance initializations
    iter = 98
    rand_dist_syn_list = list()
    for i in range(iter):
        rand_dist_syn_list.append([(combo_targset, choice(list(zip(*pw_combo_and_dist_dict[combo_targset]))[1]), targset_to_syn_dict[cl][combo_targset])
                                   for combo_targset in combo_targset_list])

    #Initialzing lists to be compatible with multi-processing 
    #since nested dictionaries of lists (where all results are stored) cannot be updated during multi-processing
    manager = Manager()
    global rand_post_fit_pcc
    global rand_post_fit_spearman
    global rand_post_fit_rmse
    global rand_cond_num
    global rand_pre_fit_params
    rand_post_fit_pcc = manager.list()
    rand_post_fit_spearman = manager.list()
    rand_post_fit_rmse = manager.list()
    rand_cond_num = manager.list()
    rand_pre_fit_params = manager.list()

    global analysis
    analysis = 'random starting point'
    #Multi-processing to run curve fitting pipeline x98
    with Pool(processes = num_cores) as pool:
                    pool.map(pipeline, rand_dist_syn_list)

                    pool.close() 
                    pool.join()
    
    #Joining results with full dictionaires
    for (pcc_corr, spearman_corr, rmse, cn, pre_params, post_params) in zip(rand_post_fit_pcc, rand_post_fit_spearman, rand_post_fit_rmse, rand_pre_fit_params, rand_cond_num):

                    post_fit_pcc[cl][metric][a].append(pcc_corr)
                    post_fit_spearman[cl][metric][a].append(spearman_corr)
                    post_fit_rmse[cl][metric][a].append(rmse)
                    cond_num[cl][metric][a].append(cn)
                    pre_fit_params[cl][metric][a].append(pre_params)

def run_shuffle_syn(merged_nbhd_dist_syn_list):

    """
        Shuffles synergy scores and runs curve fitting pipeline x 100

        :param merged_nbhd_dist_syn_list: list of tuples
            [(Target set 1__Target set 2, merged nbhd dist, synergy)]
        
        :return: None
    """
    global analysis
    analysis = 'shuffle synergy'

    iter = 100
    for i in range(iter):

        syn_list = list(targset_to_syn_dict[cl].values())

        #Randomly shuffle synergy scores
        random.shuffle(syn_list)

        targset_to_shuffle_syn_dict = {combo_targset: syn for (combo_targset, syn) in zip(combo_targset_list, syn_list)}

        merged_nbhd_dist_shuffle_syn_list = [(combo_targset, dist, targset_to_shuffle_syn_dict[combo_targset])
                            for (combo_targset, dist, real_syn)  in merged_nbhd_dist_syn_list]
        
        pipeline(merged_nbhd_dist_shuffle_syn_list)

        #Re-ordering synergy scores according to combo_targset_list
        syn_list = [syn for combo, syn in targset_to_syn_dict[cl].items() if combo in combo_targset_list]

def main():

    global targset_to_syn_dict
    global post_fit_pcc
    global post_fit_spearman
    global post_fit_rmse
    global cond_num
    global pre_fit_params

    global post_fit_pcc_shuffle
    global post_fit_spearman_shuffle
    global post_fit_rmse_shuffle
	
    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'
    exper_source = 'DREAM'
	
    distance_results_dir = os.path.join('../../../results/Distance_Dictionaries', exper_source, interactome_aname)
    check_dir(distance_results_dir)

    #All target pair distances across all 27 distance metrics
    combined_dist_filename = 'combined_distance_dicts.pkl'
    combined_dist_dict = pickle.load(open(os.path.join(distance_results_dir, combined_dist_filename), 'rb')) #[metric (e.g., T1_to_T2)][alphaD][target pair] = disance

    DREAM_inputs_dir = '../../../data/DREAM_inputs'
    targset2syn_file = 'targset_to_syn_dict.pkl'
    targset_to_syn_dict = pickle.load(open(os.path.join(DREAM_inputs_dir, targset2syn_file), 'rb')) #[cell line][Targetset1__Targetset2] = synergy

    #Used for initializing fitting with each combo's merged neighborhood distance
    merged_nbhd_dist_by_targset_filename = 'combined_distance_dicts_by_targset.pkl'
    merged_nbhd_dist_by_targset_dict = pickle.load(open(os.path.join(distance_results_dir, merged_nbhd_dist_by_targset_filename), 'rb')) #[metric (e.g., T1_to_T2)][alphaD][target set combo] = distance

    cl_all = ['BT-20', 'BT-474', 'CAMA-1', 'MCF7', 'MDA-MB-157']

    #Initializing nested dictionaries to save final results
    pre_fit_params = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}

    post_fit_pcc = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}
    post_fit_spearman = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}
    post_fit_rmse = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}

    post_fit_pcc_shuffle = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}
    post_fit_spearman_shuffle = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}
    post_fit_rmse_shuffle = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}
    cond_num = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}
    
    opt_pw_dist_all_cl = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}

    global cl
    global metric
    global a
    global analysis
    global model_name

    model_name = 'posyn'

    #Looping through all cell lines (n = 5), metrics (n = 6), and alphaDs (n = 6)
    for cl in cl_all:

        global combo_targset_list
        if model_name == 'posyn':
            combo_targset_list = [combo for combo in targset_to_syn_dict[cl].keys() if targset_to_syn_dict[cl][combo] > 0]
        elif model_name == 'negsyn':
            combo_targset_list = [combo for combo in targset_to_syn_dict[cl].keys() if targset_to_syn_dict[cl][combo] < 0]

        #syn_list = [syn for combo, syn in targset_to_syn_dict[cl].items() if combo in combo_targset_list]

        for metric in combined_dist_dict.keys():
            for a in combined_dist_dict[metric].keys():
                 
                global pw_combo_and_dist_dict #[combo_targset] = [(targ_pair, pw_dist)]
                pw_combo_and_dist_dict = get_cell_line_combos(combined_dist_dict[metric][a])

                global pw_dist_dict #[combo_targset] = [pw_dist1, pw_dist2, ...]
                pw_dist_dict = {combo_targset: list(zip(*combo_dist))[1] 
                                for combo_targset, combo_dist in pw_combo_and_dist_dict.items()}

                merged_nbhd_dist_syn_list = [(combo_targset, merged_nbhd_dist_by_targset_dict[metric][a][combo_targset], targset_to_syn_dict[cl][combo_targset])
                        for combo_targset in combo_targset_list] #[(Target set 1__Target set 2, merged nbhd dist, synergy)]
                
                #Run pipeline with merged neighborhood distance initialization
                analysis = 'merged neighborhood'
                opt_pw_dist_dict = pipeline(merged_nbhd_dist_syn_list)

                #Combining results across all cell lines, metrics, and alphaDs
                for combo_targset, pw_dist_tuple in opt_pw_dist_dict.items():
                    opt_pw_dist_all_cl[cl][metric][a].append((combo_targset, pw_dist_tuple[0], pw_dist_tuple[1]))

                #Run pipeline by initializing it with (1) min and (2) max pairwise distances
                min_dist_syn_list = [(combo_targset, np.min(list(zip(*pw_combo_and_dist_dict[combo_targset]))[1]), targset_to_syn_dict[cl][combo_targset])
                      for combo_targset in combo_targset_list]
                analysis = 'min'
                pipeline(min_dist_syn_list)

                analysis = 'max'
                max_dist_syn_list = [(combo_targset, np.max(list(zip(*pw_combo_and_dist_dict[combo_targset]))[1]), targset_to_syn_dict[cl][combo_targset])
                            for combo_targset in combo_targset_list]
                pipeline(max_dist_syn_list)

                #Run pipeline using random initializations (i.e., starting points) x98
                num_cores = mp.cpu_count()
                run_random_starting_point(num_cores)

                #Run pipeline using shuffled synergy scores x100
                run_shuffle_syn(merged_nbhd_dist_syn_list)
        
    
    #Saving all results
    outpath = os.path.join('../../../results/Multi_Targ_Curve_Fitting', exper_source, interactome_aname)
    check_dir(outpath)

    pickle.dump(opt_pw_dist_all_cl, open(os.path.join(outpath, 'opt_targ_pair_dict_5_cl_37_metrics.pkl'), 'wb'))

    pickle.dump(post_fit_pcc, open(os.path.join(outpath, 'post_fit_pcc_and_sig_dict_' + model_name + '.pkl'), 'wb'))
    pickle.dump(post_fit_spearman, open(os.path.join(outpath, 'post_fit_spearman_and_sig_dict_' + model_name + '.pkl'), 'wb'))
    pickle.dump(post_fit_rmse, open(os.path.join(outpath, 'post_fit_rse_dict_' + model_name + '.pkl'), 'wb'))
    pickle.dump(cond_num, open(os.path.join(outpath, 'post_fit_cond_num_' + model_name + '.pkl'), 'wb'))
    pickle.dump(pre_fit_params, open(os.path.join(outpath, 'pre_fit_curve_params_' + model_name + '.pkl'), 'wb'))

    pickle.dump(post_fit_pcc_shuffle, open(os.path.join(outpath, 'post_fit_pcc_and_sig_dict_shuffle_' + model_name + '.pkl'), 'wb'))
    pickle.dump(post_fit_spearman_shuffle, open(os.path.join(outpath, 'post_fit_spearman_and_sig_dict_shuffle_' + model_name + '.pkl'), 'wb'))
    pickle.dump(post_fit_rmse_shuffle, open(os.path.join(outpath, 'post_fit_rse_dict_shuffle_' + model_name + '.pkl'), 'wb'))


if __name__ == "__main__":
    import pickle, os, functools
    import pandas as pd
    import numpy as np
    from collections import defaultdict
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit, minimize, least_squares
    from numpy.random import choice
    from scipy.stats import pearsonr, spearmanr
    from multiprocessing import Pool, Manager
    import multiprocessing as mp
    import random

    main()