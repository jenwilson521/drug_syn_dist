"""
07/16/2025 - takes in protein class text file from PANTHER to analyze target protein classes across all drug combo screens


"""
def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)


def protein_class_pie_chart(class_to_gene_dict, exper_source):

    gene_count_list = [len(gene_list) for gene_list in class_to_gene_dict.values()]
    protein_classes = list(class_to_gene_dict.keys())

    fig = plt.figure(dpi = 300)
    fig.tight_layout()
    ax = plt.subplot()

    plt.pie(gene_count_list, labels = protein_classes)

    figure_name = exper_source + '_PANTHER_protein_classes_pie_char.png'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches = 'tight', pad_inches = 0.1)

def parent_protein_class_pie_chart(parent_class_prop_df, exper_source):

    #parent_class_prop_df_sorted = parent_class_prop_df.sort_values(by = 'Proportion_of_targets', ascending = False)

    parent_protein_classes = list(parent_class_prop_df.index)
    target_prop = parent_class_prop_df['Proportion_of_targets'].to_list() 

    for (prop, p_class) in zip(target_prop, parent_protein_classes):
        print(p_class, prop / np.sum(target_prop))

    fig = plt.figure(dpi = 300)
    fig.tight_layout()
    ax = plt.subplot()

    plt.pie(target_prop, colors = parent_class_colors,
                    wedgeprops={"linewidth": 1, "edgecolor": "black"})

    figure_name = exper_source + '_PANTHER_parent_protein_classes_pie_chart.png'
    plt.savefig(os.path.join(figure_dir, figure_name), bbox_inches = 'tight', pad_inches = 0.1) 


def main():


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
    
    exper_source_to_color = {'DREAM': 'dodgerblue',
                            'NCI_ALMANAC': 'indianred', #'navy',
                            'NCI_ALMANAC_prostate': 'indianred',
                            'Prostate_Organoid': 'firebrick'}
    
    data_dir = '../../data'
    nbhd_dir = '../../results/neighborhoods'

    global figure_dir
    figure_dir = '../../results/figures/target_info'
    check_dir(figure_dir)
    
    parent_class_prop_df_dict = dict()
    for exper_source in ['Prostate_Organoid', 'NCI_ALMANAC']:

        if exper_source == 'Nair':
            filename = 'PRN_nbhd_opt_alpha_with_db_targs.pkl'
        else:
            filename = 'PRN_nbhd_opt_alpha.pkl'

        nbhd_path = os.path.join(nbhd_dir, interactome_aname, exper_source, filename)
        nbhd_dict = pickle.load(open(nbhd_path, 'rb'))

        targs_per_screen = [targ for targ, nbhd in nbhd_dict.items() if len(nbhd) > 0]

        filename = exper_source + '_panther_protein_classes.txt'
        data_path = os.path.join(data_dir, exper_source + '_inputs', filename)

        protein_class_df = pd.read_csv(data_path, header = None, sep = '\t')
        protein_class_df.columns = ['Gene_ID', 'Input_gene_symbol', 'Gene_name', 'PANTHER_family/subfamily', 'PANTHER_protein_class', 'Species']
        protein_class_df['Mapped_gene_symbol'] = protein_class_df.apply(lambda row: row.Gene_name.split(';')[1], axis = 1)
        protein_class_df_filtered = protein_class_df[protein_class_df['Mapped_gene_symbol'] == protein_class_df['Input_gene_symbol']]
        protein_class_df_filtered_filled = protein_class_df_filtered.fillna('No_mapping')

        class_to_gene_dict = protein_class_df_filtered_filled.groupby('PANTHER_protein_class')['Input_gene_symbol'].apply(list).to_dict()
        print(class_to_gene_dict)

        c = 0
        targs_from_panther = set()
        for protein_class, gene_list in class_to_gene_dict.items():
            print(protein_class, len(gene_list))

            c += len(gene_list)

            for gene in gene_list:
                targs_from_panther.add(gene)
        
        #for targ in targs_per_screen:
        #    print(targ)

        print('Total number of genes:', c)
        print('Missing genes:', [targ for targ in targs_per_screen if targ not in targs_from_panther])

        filename = exper_source + '_panther_parent_protein_class_pie_chart_data.txt'
        data_path = os.path.join(data_dir, exper_source + '_inputs', filename)

        parent_class_prop_df = pd.read_csv(data_path, header = None, sep = '\t')
        parent_class_prop_df.columns = ['Parent_protein_class', 'Proportion_of_targets']
        parent_class_prop_df_dict[exper_source] = parent_class_prop_df.set_index('Parent_protein_class')

        
        #protein_class_pie_chart(class_to_gene_dict, exper_source)


    num_classes = max([len(parent_class_prop_df_dict[exper_source]) for exper_source in parent_class_prop_df_dict.keys()])
    print(num_classes)
    base_color_rgb = mcolors.to_rgb('brown') #mcolors.to_rgb(exper_source_to_color[exper_source])
    base_color_hsl = mcolors.rgb_to_hsv(base_color_rgb)
    global parent_class_colors
    parent_class_colors = list()
    for i in range(num_classes):
        lightness_factor = i / (num_classes - 1) * 1
        current_color = [base_color_rgb[j] + (1 - base_color_rgb[j]) * lightness_factor for j in range(3)]
        #lightness_factor = (i * (1/num_classes))
        #new_hsl = ((base_color_hsl[0] + lightness_factor) % 1.0, base_color_hsl[1], base_color_hsl[2])
        #parent_class_colors.append(mcolors.hsv_to_rgb(new_hsl))
        parent_class_colors.append(current_color)
    
    #Sort terms according to organoid percentages
    parent_class_prop_df_sorted = parent_class_prop_df_dict['Prostate_Organoid'].sort_values(by = 'Proportion_of_targets', ascending = False)
    parent_protein_classes = list(parent_class_prop_df_sorted.index)
    parent_class_prop_df_dict['Prostate_Organoid'] = parent_class_prop_df_dict['Prostate_Organoid'].reindex(parent_protein_classes)

    parent_class_prop_df_sorted = parent_class_prop_df_dict['NCI_ALMANAC'].sort_values(by = 'Proportion_of_targets', ascending = False)
    for pc in list(parent_class_prop_df_sorted.index):
        if pc not in parent_protein_classes:
            parent_protein_classes.append(pc)
    parent_class_prop_df_dict['NCI_ALMANAC'] = parent_class_prop_df_dict['NCI_ALMANAC'].reindex(parent_protein_classes)
    

    for exper_source in ['NCI_ALMANAC', 'Prostate_Organoid']:
        parent_protein_class_pie_chart(parent_class_prop_df_dict[exper_source], exper_source)




if __name__ == "__main__":
    import pickle, os, functools, math
    from collections import defaultdict, Counter
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import seaborn as sns
    import pandas as pd
    import networkx as nx
    from scipy.stats import ttest_ind, pearsonr, spearmanr, rankdata
    main()