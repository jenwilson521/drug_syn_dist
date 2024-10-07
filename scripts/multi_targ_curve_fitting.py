def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)


def main():
	
    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'
    exper_source = 'DREAM'
	
    distance_results_dir = os.path.join('../../../results/Distance_Dictionaries', exper_source, interactome_aname)
    check_dir(distance_results_dir)

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
    post_fit_pcc = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}
    post_fit_spearman = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}
    post_fit_rmse = {cl: {metric: {a: [] for a in combined_dist_dict[metric].keys()} for metric in combined_dist_dict.keys()} for cl in cl_all}


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