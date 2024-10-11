def main():

    """
        Merges optimal PageRank-Nibble neighborhoods for multi-target drugs
    
    """
    interactome_source = 'PathFX'
    threshold = 0.50
    exper_source = 'DREAM'
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    #Load optimal alpha PRN neighborhoods
    prn_nbhd_path = '../../results/PRN_neighborhoods'
    prn_nbhd_dir = os.path.join(prn_nbhd_path, interactome_aname, exper_source)
    file_name = 'PRN_nbhd_opt_alpha.pkl'
    prn_nbhd_dict = pickle.load(open(os.path.join(prn_nbhd_dir, file_name), 'rb'))

    #Load DREAM data
    DREAM_inputs_dir = '../../data/DREAM_inputs'
    targset2syn_file = 'targset_to_syn_dict.pkl'
    targset_to_syn_dict = pickle.load(open(os.path.join(DREAM_inputs_dir, targset2syn_file), 'rb')) #[cell line][Targetset1__Targetset2] = synergy

    #Create set of all target sets
    all_targsets = set()
    for cl in targset_to_syn_dict.keys():
        for targset_combo in targset_to_syn_dict[cl].keys():
            for targset in targset_combo.split('__'):
                all_targsets.add(targset)

    #For all targets within a target set, merge PRN neighborhoods 
    merged_nbhd_dict = defaultdict(set)
    for targset in all_targsets:
        for targ in targset.split('_'):
            if targ in list(prn_nbhd_dict.keys()): #Accounting for cases where experimental targets didn't have a detectable PRN neighborhood
                for n in prn_nbhd_dict[targ]:
                    merged_nbhd_dict[targset].add(n)
    
    #Save results
    merged_nbhd_filename = 'targset_to_merged_nbhd_dict.pkl'
    pickle.dump(merged_nbhd_dict, open(os.path.join(prn_nbhd_dir, merged_nbhd_filename), 'wb'))
    

if __name__ == "__main__":
    import os, pickle
    from collections import defaultdict