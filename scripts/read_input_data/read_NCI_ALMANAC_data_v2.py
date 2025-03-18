def plot_syn_dist(cl_to_dbid_syn_dict):

    fig = plt.figure(dpi = 300)
    ax = plt.subplot()
    fig.tight_layout() 

    i_subplot = 1
    for cl in cl_to_dbid_syn_dict.keys():
        ax = plt.subplot(2, 3, i_subplot)
        syn_sum = [np.sum(syn_list) for syn_list in cl_to_dbid_syn_dict[cl].values()]
        syn_mean = [np.mean(syn_list) for syn_list in cl_to_dbid_syn_dict[cl].values()]
        syn_std = [np.std(syn_list) for syn_list in cl_to_dbid_syn_dict[cl].values()]

        l = 'mean = ' + str(np.round(np.mean(syn_sum), 3)) + '\n' + 'std = ' + str(np.round(np.std(syn_sum), 3))
        n_bins = np.round(len(syn_sum) / 30, 0)
        plt.hist(syn_sum, bins = int(n_bins), alpha = 0.5, color = 'orangered', label = ' n = ' + str(len(syn_sum)))
        plt.axvline(x = np.mean(syn_sum), color = 'black', linestyle = 'dashed', label = l) 

        plt.title(cl, fontsize = 8)

        #if i_subplot in [4, 5, 6]:
        #    plt.xlabel('ComboScore', fontsize = 6)
        
        #plt.ylabel('# of combos', fontsize = 6)

        for axis in ['bottom','left']:
            ax.spines[axis].set_linewidth(1.5)
        for axis in ['top','right']:
            ax.spines[axis].set_linewidth(0)

        plt.xticks(np.arange(-600, 800, 200), fontsize = 5)
        plt.yticks(np.arange(0, 400, 100), fontsize = 5)
        plt.legend(fontsize = 5)  
        
        i_subplot += 1
    
    figure_name = 'NCI_ALMANAC_synergy_distribtuion_6_breast.png'
    plt.savefig(os.path.join(figure_path, figure_name), bbox_inches = 'tight', pad_inches = 0.1)

def main():

    global figure_path
    data_path = '../../data'
    resources_path = '../../resources'
    figure_path = '../../results/figures'

    # Make cell line -> combo -> combo score dict
    NCI_file = 'ComboDrugGrowth_Nov2017.csv'
    NCI_df = pd.read_csv(os.path.join(data_path, NCI_file), sep=',')
    breast_df = NCI_df[(NCI_df.PANEL == 'Breast Cancer') & (NCI_df.VALID == 'Y')]
    breast_df = breast_df.drop(breast_df[breast_df['SCORE'].isnull()].index) #Removing single drug experiments
    breast_df['COMBONAME'] = breast_df.apply(lambda row: '__'.join(sorted([str(int(row.NSC1)), str(int(row.NSC2))])), axis = 1)
    #print(breast_df)

    cl_to_combo_dict = {k: f.groupby('COMBONAME')['SCORE'].apply(list).to_dict()
     for k, f in breast_df.groupby('CELLNAME')}
    
    single_ncid_set = set()
    for cl in cl_to_combo_dict.keys():
        for combo in cl_to_combo_dict[cl].keys():
            [D1, D2] = combo.split('__')
            single_ncid_set.add(D1)
            single_ncid_set.add(D2)

    # Map NCI identifiers to common names (from Jen)
    common_names_file = 'ComboCompoundNames_all.txt'
    ncid_to_common_name_dict = defaultdict(list)
    data_lines = [l.strip() for l in open(os.path.join(data_path, common_names_file),'r').readlines()]
    for dl in data_lines:
        if '\t' in dl:
            [n, cn] = dl.split('\t')
        elif ' ' in dl:
            n = dl.split(' ')[0]
            cn = dl.replace(n+' ','')
        elif '  ' in dl:
            [n, cn] = dl.split('  ')
        elif '   ' in dl:
            [n, cn] = dl.split('   ')
        else:
            print(dl)
        ncid_to_common_name_dict[n].append(cn)

    # Consider all common name synonyms (from Jen)
    vocab_file = 'drugbank_vocabulary.csv'
    db_df = pd.read_csv(os.path.join(resources_path, vocab_file))
    db_df = db_df.fillna('none')
    dbid_to_all_names = defaultdict(list)
    syns_to_dbid = {}
    for (dbid,cn,syn_str) in zip(db_df['DrugBank ID'].to_list(),db_df['Common name'].to_list(),db_df['Synonyms'].to_list()):
        dbid_to_all_names[dbid].append(cn.lower())
        syns_to_dbid[cn.lower()] = dbid
        if syn_str == 'none':
                continue
        for dsyn in syn_str.split(' | '):
                dbid_to_all_names[dbid].append(dsyn.lower())
                syns_to_dbid[dsyn.lower()] = dbid
    
    nc_to_dbid = {}
    for (ncid, cn_list) in ncid_to_common_name_dict.items():
        for cn in cn_list:
            if cn.lower() in syns_to_dbid:
                nc_dbid = syns_to_dbid[cn.lower()]
                nc_to_dbid[ncid] = nc_dbid
    unmapped = [k for k in ncid_to_common_name_dict.keys() if k not in nc_to_dbid]

    # Manually mapping those in unmapped (from Jen)
    nc_to_dbid['63878'] = 'DB00987'
    nc_to_dbid['138783'] = 'DB06769'
    nc_to_dbid['256439'] = 'DB01177'
    nc_to_dbid['702294'] = 'DB01196'
    nc_to_dbid['750690'] = 'DB01268'
    nc_to_dbid['707389'] = 'DB08871'
    nc_to_dbid['737754'] = 'DB06589'
    nc_to_dbid['753082'] = 'DB08881'
    nc_to_dbid['763371'] = 'DB08877'
    
    # Map DBID to targets (filtering for druggable and within the interactome)
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
    print('# of drugs that target protein groups:', sum(i > 1 for i in num_targ_list))
    print('# of DB IDs:', len(set(drugbank_df['drugbank_id'].to_list())))
   
    c = 0
    dbid_to_targs_dict = defaultdict(list)
    for ncid in single_ncid_set:
        all_targs = set()
        dbid = nc_to_dbid[ncid]
        targ_list = drugbank_df[drugbank_df['drugbank_id'] == dbid].gene_symbol.tolist()
        for targ_str in targ_list:
            if targ_str == 'No mapping':
                continue

            if 'target' in drugbank_df[(drugbank_df['drugbank_id'] == dbid) & 
                            (drugbank_df['gene_symbol'] == targ_str)].category.tolist():
                
                for targ in targ_str.split(','):
                    all_targs.add(targ)
        #print(dbid, all_targs)
        
        if len(all_targs) == 0:
            continue
        else:
            dbid_to_targs_dict[dbid] = list(all_targs)
    
    a = 0
    d = 0
    unique_targs = set()
    for dbid in dbid_to_targs_dict.keys():
        #print(dbid, dbid_to_targs_dict[dbid])
        a += len(dbid_to_targs_dict[dbid])
        d += 1
        for targ in dbid_to_targs_dict[dbid]:
            unique_targs.add(targ)
    print('Average number of targets per drug:', a/d, d)
    print('Unique targets:', unique_targs, len(unique_targs))

    #Filtering for interactome protein targets
    dbid_to_int_targs_dict = defaultdict(list)
    for dbid in dbid_to_targs_dict.keys():
        for targ in dbid_to_targs_dict[dbid]:
            if targ in int_nodes:
                dbid_to_int_targs_dict[dbid].append(targ)
    
    a = 0
    d = 0
    unique_targs = set()
    for dbid in dbid_to_int_targs_dict.keys():
        print(dbid, dbid_to_int_targs_dict[dbid], len(dbid_to_int_targs_dict[dbid]))
        a += len(dbid_to_targs_dict[dbid])
        d += 1
        for targ in dbid_to_int_targs_dict[dbid]:
            unique_targs.add(targ)
    print('Average number of interactome targets per drug:', a/d, d)
    print('Unique interactome targets:', unique_targs, len(unique_targs))

    #print([(dbid, targs) for (dbid, targs) in dbid_to_targs_dict.items() if dbid not in dbid_to_int_targs_dict.keys()])

    #Load old DrugBank results for comparison
    outpath = os.path.join(data_path, 'NCI_ALMANAC_inputs')
    old_db_dbid_to_int_targs_dict = pickle.load(open(os.path.join(outpath, 'dbid_to_targs_dict_db_050120.pkl'), 'rb')) #[DBID] = [list of targets]

    #Which DB IDs are excluded in the old version of DrugBank?
    print([(dbid, targs) for (dbid, targs) in dbid_to_int_targs_dict.items() if dbid not in old_db_dbid_to_int_targs_dict.keys()])

    #Which DB IDs are excluded in the new version of DrugBank?
    print([(dbid, targs) for (dbid, targs) in old_db_dbid_to_int_targs_dict.items() if dbid not in dbid_to_int_targs_dict.keys()])

    #For the overlapping DB IDs, how similar are their targets?
    unique_new_targs_dict = defaultdict(list)
    unique_old_targs_dict = defaultdict(list)
    for dbid in dbid_to_int_targs_dict.keys():
        if dbid in old_db_dbid_to_int_targs_dict.keys():
            unique_new_targs_dict[dbid] = [targ for targ in dbid_to_int_targs_dict[dbid] if targ not in old_db_dbid_to_int_targs_dict[dbid]]
            unique_old_targs_dict[dbid] = [targ for targ in old_db_dbid_to_int_targs_dict[dbid] if targ not in dbid_to_int_targs_dict[dbid]]

    c = 0
    for dbid in unique_new_targs_dict.keys():
        if len(unique_new_targs_dict[dbid]) != 0:
            #print(dbid, unique_new_targs_dict[dbid])
            c += 1
    print('Unique targets in new DrugBank:', c)

    c = 0
    for dbid in unique_old_targs_dict.keys():
        if len(unique_old_targs_dict[dbid]) != 0:
            #print(dbid, unique_old_targs_dict[dbid])
            c += 1
    print('Unique targets in old DrugBank:', c)

    # Final filtering of combos with scores based on whether both drugs have targets in the interactome
    dd_list = functools.partial(defaultdict, list)
    cl_to_dbid_syn_dict = defaultdict(dd_list)
    for cl in cl_to_combo_dict.keys():
        for combo in cl_to_combo_dict[cl].keys():
            [ncid_1, ncid_2] = combo.split('__')
            dbid_1 = nc_to_dbid[ncid_1]
            dbid_2 = nc_to_dbid[ncid_2]

            if dbid_1 in dbid_to_int_targs_dict.keys() and dbid_2 in dbid_to_int_targs_dict.keys():
                combo_dbid = '__'.join(sorted([dbid_1, dbid_2]))
                cl_to_dbid_syn_dict[cl][combo_dbid] = cl_to_combo_dict[cl][combo]
    
    cl_to_dbid_targs_dict = defaultdict(dd_list)
    for cl in cl_to_dbid_syn_dict.keys():
        print(cl, len(cl_to_dbid_syn_dict[cl].keys()))
        for combo in cl_to_dbid_syn_dict[cl].keys():
            [dbid_1, dbid_2] = combo.split('__')
            cl_to_dbid_targs_dict[cl][combo] = [list(dbid_to_int_targs_dict[dbid_1]), list(dbid_to_int_targs_dict[dbid_2])]

    #plot_syn_dist(cl_to_dbid_syn_dict)

    print([n for n in int_nodes if 'PLA2G' in n])

    plot_syn_dist(cl_to_dbid_syn_dict)

    #Saving results
    outpath = os.path.join(data_path, 'NCI_ALMANAC_inputs')
    pickle.dump(dbid_to_int_targs_dict, open(os.path.join(outpath, 'dbid_to_targs_dict_db_022825.pkl'), 'wb'))
    pickle.dump(cl_to_dbid_syn_dict, open(os.path.join(outpath, 'cl_to_dbid_syn_dict_db_022825.pkl'), 'wb'))
    pickle.dump(cl_to_dbid_targs_dict, open(os.path.join(outpath, 'cl_to_dbid_targs_dict_db_022825.pkl'), 'wb'))


    
    #for targ in unique_targs:
    #    print(targ)


if __name__ == "__main__":
    import os, pickle, functools
    import pandas as pd
    from collections import defaultdict
    import numpy as np
    import networkx as nx
    import matplotlib.pyplot as plt
    main()
