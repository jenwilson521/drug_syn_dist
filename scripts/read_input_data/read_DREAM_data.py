"""
    Reads and filters DREAM challenge data
    Reconfigures data structure so it is compatible with the rest of the pipeline

"""
def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def main():

	#Load interactome
    interactome_source = 'PathFX'
    threshold = 0.50
    exper_source = 'GI'
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    interactome_path = '../../data/interactomes'
    interactome_file = str(threshold) + 'thr_filtered_' + interactome_source + '_scored_interactome.txt'
    interactome_df = pd.read_csv(os.path.join(interactome_path, interactome_file), sep = '\t', na_filter = False)
    interactome = nx.from_pandas_edgelist(interactome_df, 'protein1', 'protein2', edge_attr = 'cost')
    int_nodes = list(interactome.nodes())

    #Load DREAM data
    rfile = '../../data/OI_combinations_synergy_scores_final.txt'
    DREAM_df = pd.read_csv(rfile,sep='\t')
    DREAM_df = DREAM_df[DREAM_df.QA == 1] #only retaining combos with high quality data
    DREAM_df = DREAM_df.reset_index(drop=True)
    DREAM_df = DREAM_df.replace('Gemcitibine', 'Gemcitabine') 

    #Load DrugBank info
    dint = pickle.load(open('../../resources/Pfx050120_dint.pkl','rb'))
    drugbank_file = '../../resources/Drugbank050120.xlsx'
    drugbank_df = pd.read_excel(drugbank_file, sheet_name = 'DrugsToTargets')
	
    dd_list = functools.partial(defaultdict, list)
    dd_float = functools.partial(defaultdict, float)

    #5 breast cancer cell lines used for model validation
    cell_lines = ['BT-474', 'BT-20', 'CAMA-1', 'MCF7', 'MDA-MB-157']

    #Erroneous paralogs identified in post-processing
    targs_to_remove = ['AKT1S1', 'AKTIP', 'PRKCSH', 'PRKCDBP', 'ERBB2IP', 'FGFR1OP', 
                       'FGFR1OP2', 'JAKMIP1', 'Abl', 'Abacavir', 'Gemcitabine', 'IKBKAP', 'M2']
    
    combo_to_syn_dict = defaultdict(dd_list)
    targset_to_combo_dict = defaultdict(dd_list) #[cell line][T1_T2_T3__TT1_TT2__TT3__TT4] = ['D1__D2', 'D1__D2',...]
    for cl in cell_lines:

        comp_a = DREAM_df[DREAM_df['CELL_LINE'] == cl].COMPOUND_A.tolist()
        comp_b = DREAM_df[DREAM_df['CELL_LINE'] == cl].COMPOUND_B.tolist()
        exper_syn = DREAM_df[DREAM_df['CELL_LINE'] == cl].SYNERGY_SCORE.tolist()

        for a, b, syn in zip(comp_a, comp_b, exper_syn):
            combo_to_syn_dict[cl]['__'.join(sorted([a, b]))].append(syn)

        for combo in combo_to_syn_dict[cl].keys():
            combo_target_set = list()
            combo_split = combo.split('__') #Example, combo_split[0] = 'MAP2K_2'
            for multi_targ in combo_split:
                target_set = list()
                if multi_targ == 'CarboTaxol': #CarboTaxol = Paclitxel + Carboplatin
                    multi_targ_split = ['Paclitaxel', 'Carboplatin']
                else:
                    multi_targ_split = multi_targ.split('_') #Example, multi_targ_split = ['MAP2K', '2']

                for targ in multi_targ_split:
                    

                    if targ in dint: #For drugs with identified names in DrugBank
                        for n in dint[targ]:
                            if 'target' in drugbank_df[(drugbank_df['DrugName'] == targ) & 
                                                       (drugbank_df['geneSym'] == n)].category.tolist():
                                #Only saving targets (not transporters, enzymes, etc)
                                if n in int_nodes:
                                    target_set.append(n)
                    elif targ in int_nodes: #For drugs whose targets are specific paralogs, e.g., MAP2K1
                        target_set.append(targ)
                    elif targ.isnumeric(): #For drugs named with additional numeric identifier, e.g., MAP2K_1
                        continue
                    else: #For drugs whose targets consist of paralagous proteins
                        if targ == 'IAP':
                            targ = 'BIRC' 
                        targ_to_int_node = [n.startswith(targ) for n in int_nodes] #Example, targ_to_int_node =['MAP2K1', 'MAP2K2', ...]
                        int_index = [ind for ind in range(len(targ_to_int_node)) if targ_to_int_node[ind]]
                        targ_int_nodes = [int_nodes[n] for n in int_index]
                        if targ == 'CHK1':
                            targ_int_nodes = ['CHEK1'] 
                        for n in targ_int_nodes:
                            if n not in targs_to_remove:
                                target_set.append(n)


                combo_target_set.append('_'.join(sorted(target_set)))

            combo_target_set_str = '__'.join(combo_target_set)
            targset_to_combo_dict[cl][combo_target_set_str].append(combo)

    #Combining all synergy scores (and averaging) for a given combination of target sets since some target sets corresponded to more than one drug
    targset_to_syn_dict = defaultdict(dd_float)
    for cl in cell_lines:
        for targset_combo in targset_to_combo_dict[cl].keys():
            syn_per_targset_combo = list()
            for combo in targset_to_combo_dict[cl][targset_combo]:
                syn_list = combo_to_syn_dict[cl][combo]

                for syn in syn_list:
                    syn_per_targset_combo.append(syn)
            
            targset_to_syn_dict[cl][targset_combo] = np.mean(syn_per_targset_combo)
    
    #Save results
    outpath = '../../data/DREAM_inputs'
    check_dir(outpath)

    file_name = 'targset_to_syn_dict.pkl'
    pickle.dump(targset_to_syn_dict, open(os.path.join(outpath, file_name), 'wb'))

    file_name = 'targset_to_combo_dict.pkl'
    pickle.dump(targset_to_combo_dict, open(os.path.join(outpath, file_name), 'wb'))

    file_name = 'combo_to_syn_dict.pkl'
    pickle.dump(combo_to_syn_dict, open(os.path.join(outpath, file_name), 'wb'))

    

if __name__ == "__main__":
    import pandas as pd
    import pickle, os, re, functools
    import networkx as nx
    import numpy as np
    from collections import defaultdict

    dd_list = functools.partial(defaultdict, list)
    main()