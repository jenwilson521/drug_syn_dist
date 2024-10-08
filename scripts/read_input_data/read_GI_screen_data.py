"""
    Reads and gene interaction screen data from Han et al (NAt Biotechnol, 2017)
    Reconfigures data structure so it is compatible with the rest of the pipeline

"""
def check_dir(dir_name):
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)

def main():

    #Load GI screen data
    GI_file = '../../data/Han_et_al_GI_Scores.xlsx'
    GI_df = pd.read_excel(GI_file)

    GT_scores = dict([x for x in zip(GI_df['Gene_Pair'].to_list(), GI_df['GeneAB GIT score'].to_list())]) #['Gene1__Gene2'] = GI score
    growth_phen = dict([x for x in zip(GI_df['Gene_Pair'].to_list(), GI_df['GeneAB Observed\nγ Phenotype'].to_list())]) #['Gene1__Gene2'] = growth phenotype

    #Save results
    outpath = '../../data/GI_screen_inputs'
    check_dir(outpath)

    file_name = 'GT_scores_dict.pkl'
    pickle.dump(GT_scores, open(os.path.join(outpath, file_name), 'wb'))

    file_name = 'growth_phen_dict.pkl'
    pickle.dump(growth_phen, open(os.path.join(outpath, file_name), 'wb'))

if __name__ == "__main__":
    import os, pickle
    import pandas as pd
    main()