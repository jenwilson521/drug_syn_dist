def main():

    data_path = '../../data'
    resources_path = '../../resources'
    results_path = '../../results'

    filename = 'Nair_et_al_drug_names.xlsx'
    anchor_drugs_df = pd.read_excel(os.path.join(data_path, filename), sheet_name = 'c_AnchorDrugsConcentration')
    anchor_drugs = anchor_drugs_df['Anchor_Drug_Name'].to_list()

    library_drugs_df = pd.read_excel(os.path.join(data_path, filename), sheet_name = 'd_LibraryDrugsConcentration')
    library_drugs = library_drugs_df['Library_Drug_Name'].to_list()
    all_drugs = set(anchor_drugs).union(set(library_drugs))
    print(all_drugs, len(all_drugs))

    #Load interactome
    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_file = str(threshold) + 'thr_filtered_' + interactome_source + '_scored_interactome.txt'
    interactome_df = pd.read_csv(os.path.join(data_path, interactome_file), sep = '\t', na_filter = False)
    interactome = nx.from_pandas_edgelist(interactome_df, 'protein1', 'protein2', edge_attr = 'cost')
    int_nodes = list(interactome.nodes())

    #New new DrugBank   
    drugbank_df = pd.read_excel(os.path.join(resources_path, 'DrugBank022825.xlsx'))
    targ_list = drugbank_df['gene_symbol'].to_list()
    num_targ_list = [len(targ_str.split(',')) for targ_str in targ_list]


    """
    Option 3 - only extracting drugs with DB IDs and getting their targets from DrugBank
    
    """
    # Consider all common name synonyms (from Jen)
    vocab_file = 'drugbank_vocabulary.csv'
    syns_df = pd.read_csv(os.path.join(resources_path, vocab_file))
    syns_df = syns_df.fillna('none')
    dbid_to_all_names = defaultdict(list)
    syns_to_dbid = {}
    for (dbid,cn,syn_str) in zip(syns_df['DrugBank ID'].to_list(),syns_df['Common name'].to_list(),syns_df['Synonyms'].to_list()):
        dbid_to_all_names[dbid].append(cn.lower())
        syns_to_dbid[cn.lower()] = dbid
        if syn_str == 'none':
                continue
        for dsyn in syn_str.split(' | '):
                dbid_to_all_names[dbid].append(dsyn.lower())
                syns_to_dbid[dsyn.lower()] = dbid
    
    name_to_dbid = {}
    c = 0
    for drug in all_drugs:
        if drug.lower() in syns_to_dbid:
                dbid = syns_to_dbid[drug.lower()]
                name_to_dbid[drug] = dbid
                c += 1
    print('Number of drugs with DB ID:', c)

    #for drug in all_drugs:
    #    if drug not in name_to_dbid.keys():
    #        print(drug)

    dbid_to_targs_dict = defaultdict(list)
    for dbid in name_to_dbid.values():
        all_targs = set()
        targ_list = drugbank_df[drugbank_df['drugbank_id'] == dbid].gene_symbol.tolist()
        for targ_str in targ_list:
            if targ_str == 'No mapping':
                continue

            if 'target' in drugbank_df[(drugbank_df['drugbank_id'] == dbid) & 
                            (drugbank_df['gene_symbol'] == targ_str)].category.tolist():
                
                for targ in targ_str.split(','):
                    if targ in int_nodes:
                        all_targs.add(targ)
        #print(dbid, all_targs)
        
        if len(all_targs) == 0:
            continue
        else:
            dbid_to_targs_dict[dbid] = list(all_targs)
    print(dbid_to_targs_dict, len(dbid_to_targs_dict.keys()))

    drugbank_targs = set()
    for dbid, targs in dbid_to_targs_dict.items():
        for targ in targs:
            drugbank_targs.add(targ)
    #print(drugbank_targs, len(drugbank_targs))

    #Load targets from NCI ALMANAC
    outpath = os.path.join(data_path, 'NCI_ALMANAC_inputs')
    nci_dbid_to_targs_dict = pickle.load(open(os.path.join(outpath, 'dbid_to_targs_dict_db_022825.pkl'), 'rb')) #[DBID] = [list of targets]
    nci_targs = set()
    for dbid, targs in nci_dbid_to_targs_dict.items():
        for targ in targs:
            nci_targs.add(targ)
    
    #Load DREAM and GI target neighborhoods
    nbhd_file = 'PRN_nbhd_opt_alpha.pkl'
    exper_source = 'GI'
    GI_nbhds = pickle.load(open(os.path.join(results_path, 'neighborhoods', exper_source, nbhd_file), 'rb'))
    exper_source = 'DREAM'
    DREAM_nbhds = pickle.load(open(os.path.join(results_path, 'neighborhoods', exper_source, nbhd_file), 'rb'))
    
    targs_with_nbhds = set()
    for targ in drugbank_targs:
        if targ in list(GI_nbhds.keys()):
            targs_with_nbhds.add(targ)
        if targ in list(DREAM_nbhds.keys()):
            targs_with_nbhds.add(targ)
        if targ in nci_targs:
            targs_with_nbhds.add(targ)
    print('DrugBank targets with GI, DREAM, NCI neighborhoods:', targs_with_nbhds, len(targs_with_nbhds))


    """
    Option 1 - taking Nair et al target lists as they are
    
    """
    anchor_drugs_to_targs = dict([x for x in zip(anchor_drugs_df['Anchor_Drug_Name'].to_list(), anchor_drugs_df['all_targets'].to_list())])
    library_drugs_to_targs = dict([x for x in zip(library_drugs_df['Library_Drug_Name'].to_list(), library_drugs_df['all_targets'].to_list())])

    drug_to_targs_dict = defaultdict(list)
    nair_targs = set()
    for drug, targ_str in anchor_drugs_to_targs.items():
        all_targs = set()
        for targ in targ_str.split(','):
            if targ in int_nodes:
                nair_targs.add(targ)
                all_targs.add(targ)
            
        if len(all_targs) == 0:
            continue
        else:
            drug_to_targs_dict[drug] = list(all_targs)
    
    for drug, targ_str in library_drugs_to_targs.items():
        all_targs = set()
        for targ in targ_str.split(','):
            if targ in int_nodes:
                nair_targs.add(targ)
                all_targs.add(targ)
            
        if len(all_targs) == 0:
            continue
        else:
            drug_to_targs_dict[drug] = list(all_targs)
    
    print('Number of Nair drugs with interactome targets:', len(drug_to_targs_dict.keys()))
    print(nair_targs, len(nair_targs))

    targs_with_nbhds = set()
    for targ in nair_targs:
        if targ in list(GI_nbhds.keys()):
            targs_with_nbhds.add(targ)
        if targ in list(DREAM_nbhds.keys()):
            targs_with_nbhds.add(targ)
        if targ in nci_targs:
            targs_with_nbhds.add(targ)
    print('Nair targets with GI, DREAM, NCI neighborhoods:', targs_with_nbhds, len(targs_with_nbhds))
        
    

if __name__ == "__main__":
    import os, pickle, functools
    import pandas as pd
    from collections import defaultdict
    import numpy as np
    import networkx as nx
    import matplotlib.pyplot as plt
    main()