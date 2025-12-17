"""
07/15/2025 - enrichment analysis for target sets across all drug combo screens


"""
def run_goea(aname, protein_list, depth_cutoff):
        
    """
        Runs GOEA for a target's PageRank-Nibble neighborhood
        Saves ALL (even those with Benjamini-Hochberg corrected p-val > 0.05) in text file (one text file/target)
        Saves significant results (GO term and Benjamini-Hochberg corrected p-value) in dictionary

    
    """

    prn_ids = [id for id, sym in id2sym_dict.items() if sym in protein_list] #Neighborhood gene symbol -> gene ID
    goea_results_all = goeaobj.run_study(prn_ids, prt = None) #Run GOEA
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05] #Retain sig. results
    go_pval_tuple = [(r.GO, r.name, r.p_fdr_bh) for r in goea_results_all
                      if r.p_fdr_bh < 0.05 and r.NS == 'MF' and r.depth >= depth_cutoff] #Save [(GO, pval)] for sig. BP results

    #textfile_name = targ + '_PRN_opt_alpha_GOEA.txt'
    #goeaobj.wr_txt(os.path.join(goea_results_path, textfile_name), goea_results_sig)

    goea_dict[aname] = go_pval_tuple

def goea_hm(terms_pvals_tuples, title):

    terms = list(zip(*terms_pvals_tuples))[0]
    pvals_1 = [-math.log10(pval) for pval in list(zip(*terms_pvals_tuples))[1]]
    pvals_2 = [-math.log10(pval) for pval in list(zip(*terms_pvals_tuples))[2]]
    pvals_all = [pvals_1, pvals_2]

    fig = plt.figure(dpi = 300, figsize = (10, 2))
    fig.tight_layout()
    ax = plt.subplot()

    hm = sns.heatmap(pvals_all, vmin = 1.1, cmap = 'Reds', linecolor = 'black', linewidths = 1, xticklabels = False, yticklabels = False)
    ax.set_xticks(range(0, len(terms)))
    ax.set_yticks(range(0, 2))

    xticks = ax.get_xticks()
    yticks = ax.get_yticks()

    # Shift tick positions to center
    ax.set_xticks(xticks + 0.5)
    ax.set_yticks(yticks + 0.5) 

    figure_name = title + '_MF_enrich_hm.png'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches = 'tight', pad_inches = 0.1)

def goea_bar_chart(terms_pvals_tuples, title):
    title_split = title.split('__')

    fig = plt.figure(dpi = 300)
    fig.tight_layout()
    ax = plt.subplot()

    ax.grid(axis = 'x', zorder = 0)

    bar_width = 0.35

    terms = list(zip(*terms_pvals_tuples))[0]
    pvals_1 = [-math.log10(pval) for pval in list(zip(*terms_pvals_tuples))[1]]
    pvals_2 = [-math.log10(pval) for pval in list(zip(*terms_pvals_tuples))[2]]

    x = np.arange(0, len(terms))
    ax.barh(x - bar_width / 2, pvals_1[::-1], color = exper_source_to_color[title_split[0]],
            height = bar_width, zorder = 3, edgecolor = 'black', linewidth = 0.75)
    
    ax.barh(x + bar_width / 2, pvals_2[::-1], color = exper_source_to_color[title_split[1]],
            height = bar_width, zorder = 3, edgecolor = 'black', linewidth = 0.75) #Flipping order so top-enriched are at the top
    
    ax.tick_params(axis = 'x', width = 1)
    ax.tick_params(axis = 'y', width = 1)
    #ax.get_yaxis().set_ticks([])
    plt.yticks(x, labels = [])
    plt.xticks(np.arange(0, 80, 10), labels = [])

    for axis in ['bottom']:
        ax.spines[axis].set_linewidth(1)
    for axis in ['top','right', 'left']:
        ax.spines[axis].set_linewidth(0) 
    
    figure_name = title + '_MF_enrich_bars.png'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches = 'tight', pad_inches = 0.1)


def main():

    global figure_dir
    figure_dir = '../../results/figures/target_info'


    dd_list = functools.partial(defaultdict, list)
    dd_set = functools.partial(defaultdict, set)

    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'
    exper_source_all = ['DREAM', 'NCI_ALMANAC', 'Nair', 'Prostate_Organoid']


    exper_source_to_cl = {'DREAM': ['BT-20', 'BT-474', 'CAMA-1', 'MCF7_DREAM', 'MDA-MB-157'],
                        'NCI_ALMANAC': ['BT-549', 'HS 578T', 'MCF7_NCI', 'MDA-MB-231/ATCC', 'MDA-MB-468', 'T-47D'],
                        'Nair': ['A549', 'EKVX', 'HOP-62', 'NCI-H23', 'NCI-H322M', 'NCI-H460', 'NCI-H522'], 
                        'NCI_ALMANAC_prostate': ['DU-145', 'PC-3'],
                        'Prostate_Organoid': ['22RV1', 'MDVR']}
    
    global exper_source_to_color
    exper_source_to_color = {'DREAM': 'dodgerblue',
                            'NCI_ALMANAC': 'indianred', #'navy',
                            'NCI_ALMANAC_prostate': 'indianred',
                            'Prostate_Organoid': 'firebrick'}
    
    nbhd_dir = '../../results/neighborhoods'
    
    #INITIALIZE GOEA
    fin_gene2go = download_ncbi_associations() #download associations
    obodag = GODag("../../resources/go-basic.obo") #load ontologies

    # Read NCBI's gene2go. Store annotations in a list of namedtuples
    objanno = Gene2GoReader(fin_gene2go, taxids=[9606]) #homo sapien taxid = 9606
    ns2assoc = objanno.get_ns2assc() #Dictionary with MF, CC, BP keys and gene -> GO term annotations ex: 109617024: {'GO:0006396'}
    
    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated homo sapien genes". format(NS = nspc, N = len(id2gos)))

    #Load background genes (PathFX 0.5 threshold interactome)
    sym2id_df = pd.read_csv('../../data/0.5thr_PathFX_int_GeneSymbol2GeneID.csv', sep = ',')

    #Create dictionary that maps gene ID to symbol in order to convert target neighborhoods to gene IDs
    global id2sym_dict
    id2sym_dict = dict()
    for i, row in sym2id_df.iterrows():
        if row['converted_alias'] != 'None':
            id = int(row['converted_alias'])
            sym = row['initial_alias']
            id2sym_dict[id] = sym
    
    global goeaobj
    goeaobj = GOEnrichmentStudyNS(
        id2sym_dict.keys(), # List of PathFX protein-coding genes with Gene IDs
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    
    targset_across_screens_dict = dict()
    for exper_source in exper_source_all:

        if exper_source == 'Nair':
            filename = 'PRN_nbhd_opt_alpha_with_db_targs.pkl'
        else:
            filename = 'PRN_nbhd_opt_alpha.pkl'
        

        nbhd_path = os.path.join(nbhd_dir, interactome_aname, exper_source, filename)
        nbhd_dict = pickle.load(open(nbhd_path, 'rb'))

        targs_per_screen = [targ for targ, nbhd in nbhd_dict.items() if len(nbhd) > 0]
        targset_across_screens_dict[exper_source] = targs_per_screen
    
    #NCI_ALMANAC_prostate are a subset of NCI_ALMANAC
    exper_source = 'NCI_ALMANAC_prostate'
    data_dir = os.path.join('../../data', 'NCI_ALMANAC_inputs')
    filename = 'prostate_cl_DB00648_DB00675_combos_dbid_to_targs_dict_db_022825.pkl'
    data_path = os.path.join(data_dir, filename)
    dbid_to_targs_dict = pickle.load(open(data_path, 'rb'))
    targset = set()
    for targset_list in dbid_to_targs_dict.values():
        for targ in targset_list:
            if targ in targset_across_screens_dict['NCI_ALMANAC']:
                targset.add(targ)
    
    
    targset_across_screens_dict['NCI_ALMANAC_prostate'] = list(targset)
    
    """
    #GOATools enrichment
    global goea_dict
    goea_dict = defaultdict(list)
    mf_dict = defaultdict(list)
    depth_cutoff = 7
    for exper_source in ['Prostate_Organoid', 'NCI_ALMANAC_prostate']: #targset_across_screens_dict.keys():
        aname = exper_source + '_targset'
        print(aname, len(targset_across_screens_dict[exper_source]))
        run_goea(aname, targset_across_screens_dict[exper_source], 0)
        print(goea_dict[aname])

        mf_dict[aname] = list(zip(*goea_dict[aname]))[1]
    
    print([mf for mf in mf_dict['NCI_ALMANAC_prostate_targset'] if mf not in mf_dict['Prostate_Organoid_targset']])
    print([mf for mf in mf_dict['Prostate_Organoid_targset'] if mf not in mf_dict['NCI_ALMANAC_prostate_targset']])

    for exper_source in ['Prostate_Organoid']:
        print(exper_source)
        for targ in targset_across_screens_dict[exper_source]:
            print(targ)
    """

    #Enrichr enrichment (NCI ALMANAC and prostate organoid)
    exper_source_to_goea_dict = dict()
    top_terms_per_screen = dict()
    for exper_source in ['NCI_ALMANAC', 'Prostate_Organoid']:
        data_dir = os.path.join('../../data', exper_source + '_inputs')
        filename = exper_source + '_GO_MF_2025_table.txt'
        data_path = os.path.join(data_dir, filename)
        goea_df = pd.read_csv(data_path, sep = '\t')
        exper_source_to_goea_dict[exper_source] = {term: pval for (term, pval) in zip(goea_df['Term'].to_list(), goea_df['Adjusted P-value'].to_list())}

        top_terms_per_screen[exper_source] = sorted(exper_source_to_goea_dict[exper_source].items(), key = lambda x:x[1])[:15]
    

    nci_top_terms = list(zip(*(top_terms_per_screen['NCI_ALMANAC'])))[0]
    org_top_terms = list(zip(*(top_terms_per_screen['Prostate_Organoid'])))[0]

    top_terms_to_plot = list(set(nci_top_terms) | set(org_top_terms))
    terms_pvals_tuples = list()
    for term in top_terms_to_plot:
        term_list = list()
        term_list.append(term)
        for exper_source in ['NCI_ALMANAC', 'Prostate_Organoid']:
            if term in exper_source_to_goea_dict[exper_source].keys():
                term_list.append(exper_source_to_goea_dict[exper_source][term])
            else:
                term_list.append(1)
        
        terms_pvals_tuples.append(tuple(term_list))
    
    terms_pvals_tuples_sorted = sorted(terms_pvals_tuples, key = lambda x:x[2])
    #goea_hm(terms_pvals_tuples_sorted, '__'.join(['NCI_ALMANAC', 'Prostate_Organoid']))
    goea_bar_chart(terms_pvals_tuples_sorted, '__'.join(['NCI_ALMANAC', 'Prostate_Organoid']))

    for (term, pval_1, pval_2) in terms_pvals_tuples_sorted:
        term_split = term.split(' ')
        go_term_p = term_split[-1].replace('(', '')
        go_term = go_term_p.replace(')', '')
        print(term, pval_1, pval_2)




if __name__ == "__main__":
    import pickle, os, functools, math
    from collections import defaultdict
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import networkx as nx
    from scipy.stats import ttest_ind, pearsonr, spearmanr, rankdata

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