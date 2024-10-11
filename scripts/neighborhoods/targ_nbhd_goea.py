def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def make_goea_tables(targ):
    """
        Runs GOEA for single target neighborhood

        :param targ: string
            Single target
        
        :return rows: list
            List of enrichment results for all relevant GO terms, to be converted to a dataframe
    
    """

    prn_ids = [id for id, sym in id2sym_dict.items() if sym in opt_prn_nbhd_dict[targ]] #Neighborhood gene symbol -> gene ID
    goea_results_all = goeaobj.run_study(prn_ids, prt = None) #Run GOEA
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

def run_make_goea_tables(targs):
    """
        Runs GOEA for list of targets and saves results in one spreadsheet

        :param targs: list
        
    """
    col_names = ['GO term', 'Aspect', 'BH-corrected p-val', 'Ratio of term in nbhd',
                     'Ratio of term in interactome', 'Depth', 'Name', 'Term-nbhd overlapping proteins']
    targ_df_dict = {targ: pd.DataFrame() for targ in prn_nbhd_targs}
    for targ in prn_nbhd_targs:
        rows = make_goea_tables(targ)
        targ_df = pd.DataFrame(rows, columns = col_names)
        targ_df_dict[targ] = targ_df
    
    filename = 'Example_GOEA_tables.xlsx'
    with pd.ExcelWriter(os.path.join(goea_results_path, filename), engine='xlsxwriter') as writer:
            for targ in targ_df_dict.keys():
                targ_df_dict[targ].to_excel(writer, sheet_name = targ)

def firstn_goea(targ):
    """
        Runs GOEA for a target's first neighbors neighborhood
        Saves ALL (even those with Benjamini-Hochberg corrected p-val > 0.05) in text file (one text file/target)
        Saves significant results (GO term and Benjamini-Hochberg corrected p-value) in dictionary

        :param targ: string
    
    """
    #Gene symbol -> Gene ID
    firstn_ids = [id for id, sym in id2sym_dict.items() if sym in firstn_nbhd_dict[targ]]
    goea_results_all = goeaobj.run_study(firstn_ids, prt = None)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05] #Retain sig. results
    go_pval_tuple = [(r.GO, r.p_fdr_bh) for r in goea_results_all if r.p_fdr_bh < 0.05 and r.NS == 'BP']

    textfile_name = targ + '_firstn_GOEA.txt'
    goeaobj.wr_txt(os.path.join(goea_results_path, textfile_name), goea_results_sig)

    firstn_go_dict[targ] = go_pval_tuple

def prn_goea(targ):
        
    """
        Runs GOEA for a target's PageRank-Nibble neighborhood
        Saves ALL (even those with Benjamini-Hochberg corrected p-val > 0.05) in text file (one text file/target)
        Saves significant results (GO term and Benjamini-Hochberg corrected p-value) in dictionary

        :param targ: string
    
    """

    prn_ids = [id for id, sym in id2sym_dict.items() if sym in opt_prn_nbhd_dict[targ]] #Neighborhood gene symbol -> gene ID
    goea_results_all = goeaobj.run_study(prn_ids, prt = None) #Run GOEA
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05] #Retain sig. results
    go_pval_tuple = [(r.GO, r.p_fdr_bh) for r in goea_results_all
                      if r.p_fdr_bh < 0.05 and r.NS == 'BP'] #Save [(GO, pval)] for sig. BP results

    textfile_name = targ + '_PRN_opt_alpha_GOEA.txt'
    goeaobj.wr_txt(os.path.join(goea_results_path, textfile_name), goea_results_sig)

    prn_go_dict[targ] = go_pval_tuple


def main():
    #Global variables
    global goea_results_path
    global prn_nbhd_targs
    global opt_prn_nbhd_dict
    global firstn_nbhd_dict
    global firstn_go_dict
    global prn_go_dict

    #Need to make all GO resources global
    global goeaobj
    global id2sym_dict

    num_cores = mp.cpu_count()
    dd_tuple = functools.partial(defaultdict, tuple)

    #Load interactome
    interactome_source = 'PathFX'
    threshold = 0.50
    exper_source = 'DREAM'
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    interactome_path = '../../data/interactomes'
    interactome_file = str(threshold) + '_filtered_' + interactome_source + '_scored_interactome.txt'
    interactome_df = pd.read_csv(os.path.join(interactome_path, interactome_file), sep = '\t', na_filter = False)
    interactome = nx.from_pandas_edgelist(interactome_df, 'protein1', 'protein2', edge_attr = 'cost')
    int_nodes = list(interactome.nodes())
    int_nodes_sorted = sorted(int_nodes)

    fin_gene2go = download_ncbi_associations() #download associations
    obodag = GODag("../../../rscs/go-basic.obo") #load ontologies

    # Read NCBI's gene2go. Store annotations in a list of namedtuples
    objanno = Gene2GoReader(fin_gene2go, taxids=[9606]) #homo sapien taxid = 9606
    ns2assoc = objanno.get_ns2assc() #Dictionary with MF, CC, BP keys and gene -> GO term annotations ex: 109617024: {'GO:0006396'}
    
    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated homo sapien genes". format(NS = nspc, N = len(id2gos)))

    #Load background genes (PathFX 0.5 threshold interactome)
    sym2id_df = pd.read_csv('../../data/0.5thr_PathFX_int_GeneSymbol2GeneID.csv', sep = ',')

    #Create dictionary that maps gene ID to symbol in order to convert target neighborhoods to gene IDs
    id2sym_dict = dict()
    for i, row in sym2id_df.iterrows():
        if row['converted_alias'] != 'None':
            id = int(row['converted_alias'])
            sym = row['initial_alias']
            id2sym_dict[id] = sym
    
    goeaobj = GOEnrichmentStudyNS(
        id2sym_dict.keys(), # List of PathFX protein-coding genes with Gene IDs
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    
    #Load neighborhoods (PRN using optimal alpha and first neighbors)
    prn_results_dir = os.path.join('../../results/PRN_neighborhoods', interactome_aname, exper_source)

    file_name = 'PRN_nbhd_opt_alpha.pkl'
    opt_prn_nbhd_dict = pickle.load(open(os.path.join(prn_results_dir, file_name), 'rb'))
    prn_nbhd_targs = [targ for targ, nbhd in opt_prn_nbhd_dict.items() if len(nbhd) > 0]

    file_name = 'firstn_nbhd.pkl'
    firstn_nbhd_dict = pickle.load(open(os.path.join(prn_results_dir, file_name), 'rb'))

    #Path where results are saved
    goea_results_path = os.path.join(prn_results_dir, 'GOEA')
    check_dir(goea_results_path)

    #For instances when GOEA is only necessary for select targets and want to save them as a spreadsheet
    select_targs = ['ATR', 'PIK3CG', 'MAP2K2', 'BCL2L1', 'EGFR']
    #run_make_goea_tables(select_targs)

    #Run GOEA for all PRN and first neighbors neighborhoods using pool process
    manager = Manager()
    #Save GOEA results across all targets in one place 
    prn_go_dict = manager.dict()
    firstn_go_dict = manager.dict()
    with Pool(processes = num_cores) as pool:
        pool.map(firstn_goea, prn_nbhd_targs)
        pool.map(prn_goea, prn_nbhd_targs)

    #Something weird happens in pool process and need to convert dictionaries back to dictionaries
    firstn_go_dict2 = dict(firstn_go_dict)
    prn_go_dict2 = dict(prn_go_dict)

    #Combine first neighbors and PRN GO results in single dictionary
    nbhd2go_dict = defaultdict(dd_tuple)
    for targ in prn_nbhd_targs:
        nbhd2go_dict['PRN opt. alpha'][targ] = prn_go_dict2[targ]
        nbhd2go_dict['Firstn'][targ] = firstn_go_dict2[targ]

    #Save dictionary
    file_name = 'nbhd2go_dict_opt_alpha_firstn.pkl'
    pickle.dump(nbhd2go_dict, open(os.path.join(goea_results_path, file_name), 'wb'))


if __name__ == "__main__":
    from __future__ import print_function
    import os, pickle, functools
    import pandas as pd
    import networkx as nx
    import numpy as np
    from collections import defaultdict
    from multiprocessing import Pool, Manager
    import multiprocessing as mp

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

    main()