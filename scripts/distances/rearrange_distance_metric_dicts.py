""" 
    Rearranging distance dictionaries for better access in the rest of the pipeline

"""

def main():
    dd_float = functools.partial(defaultdict, float)

    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'
    exper_source = 'NCI_ALMANAC'

    target_type = 'multi' #'single' #'multi'

    #Load distance dictionaries
    distance_results_dir = os.path.join('../../results/distance_dictionaries', interactome_aname, exper_source)

    sp_filename = target_type + '_target_SP_dict.pkl'
    t2t_dist_filename = target_type + '_target_t2t_diffusion_dist_0.06-0.75_alphaD_dict.pkl'
    t2n_dist_filename = target_type + '_target_t2n_diffusion_dist_0.06-0.75_alphaD_dict.pkl'
    n2n_dist_filename = target_type + '_target_n2n_diffusion_dist_0.06-0.75_alphaD_dict.pkl'

    sp_dict = pickle.load(open(os.path.join(distance_results_dir, sp_filename), 'rb'))
    t2t_dist_dict = pickle.load(open(os.path.join(distance_results_dir, t2t_dist_filename), 'rb'))
    t2n_dist_dict = pickle.load(open(os.path.join(distance_results_dir, t2n_dist_filename), 'rb'))
    n2n_dist_dict = pickle.load(open(os.path.join(distance_results_dir, n2n_dist_filename), 'rb'))

    #Initializes combined dictionaries
    distances = ['T1_to_T2', 'T2_to_T1', 'T1_to_N2', 'T2_to_N1', 'N1_to_N2', 'N2_to_N1']
    alphaDs = list(t2t_dist_dict.keys())
    #pairwise_combos = t2t_dist_dict[alphaDs[0]].keys()
    pairwise_combos = [combo for combo in sp_dict.keys() if not math.isnan(sp_dict[combo])]
    combined_distance_dicts = {d: {a: {c: [] for c in pairwise_combos} for a in alphaDs} for d in distances}
    combined_distance_dicts.update({'ASP': {0: {c: [] for c in pairwise_combos}}})

    for a in alphaDs:
        for combo in pairwise_combos:

            asp_dist = sp_dict[combo]

            if math.isnan(asp_dist):
                continue
            else:

                combined_distance_dicts['ASP'][0][combo] = sp_dict[combo] #ASP doesn't use alphaD, using 0 as a placeholder

                combined_distance_dicts['T1_to_T2'][a][combo] = t2t_dist_dict[a][combo][0]
                combined_distance_dicts['T2_to_T1'][a][combo] = t2t_dist_dict[a][combo][1]

                combined_distance_dicts['T1_to_N2'][a][combo] = t2n_dist_dict[a][combo][0]
                combined_distance_dicts['T2_to_N1'][a][combo] = t2n_dist_dict[a][combo][1]

                combined_distance_dicts['N1_to_N2'][a][combo] = n2n_dist_dict[a][combo][0]
                combined_distance_dicts['N2_to_N1'][a][combo] = n2n_dist_dict[a][combo][1]

    #Save results
    combined_dist_dicts_filename = target_type + '_target_37_dist_metric_dicts.pkl'
    pickle.dump(combined_distance_dicts, open(os.path.join(distance_results_dir, combined_dist_dicts_filename), 'wb'))


if __name__ == "__main__":
    import os, pickle, functools, math
    from collections import defaultdict
    main()