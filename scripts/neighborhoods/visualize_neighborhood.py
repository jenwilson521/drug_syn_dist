def main():

    #Load interactome
    interactome_source = 'PathFX'
    threshold = 0.50
    exper_source = 'DREAM'
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'

    interactome_path = '../../data'
    interactome_file = str(threshold) + 'thr_filtered_' + interactome_source + '_scored_interactome.txt'
    interactome_dir = os.path.join(interactome_path, interactome_file)
    interactome_df = pd.read_csv(interactome_dir, sep = '\t', na_filter = False)
    interactome = nx.from_pandas_edgelist(interactome_df, 'protein1', 'protein2', edge_attr = 'cost')

    nhbd_list = ['HDAC3', 'NCOA1', 'MED13', 'KAT2B'] #Example target neighborhood

    #Degree for all neighborhood proteins
    for n in nhbd_list:
        print(interactome.degree(n))
    
    #Neighborhood subgraph, saving it as a text file
    nhbd_subgraph = interactome.subgraph(nhbd_list)
    subgraph_df = nx.to_pandas_edgelist(nhbd_subgraph)
    print(subgraph_df)
    subgraph_df.to_csv('../../results/HDAC3_nhbd.txt', index = False, sep = ' ')

    #All first degree neighbors of a target, using list comprehension
    print([n for n in interactome.neighbors('HDAC3')])


if __name__ == "__main__":
    import os, pickle
    import networkx as nx
    import pandas as pd
    main()