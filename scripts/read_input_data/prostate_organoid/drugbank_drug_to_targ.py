"""
02/04/2025 - determining which drugs have target in drug bank

"""

def main():

    #Load DrugBank
    dint = pickle.load(open('../../../resources/Pfx050120_dint.pkl','rb'))
    drugbank_file = '../../../resources/Drugbank050120.xlsx'
    drugbank_df = pd.read_excel(drugbank_file, sheet_name = 'DrugsToTargets')


    #Load data
    cl_file = '../../../data/Bliss_score_list_92_drugs_22RV1.xlsx'
    cl_df = pd.read_excel(cl_file, sheet_name = '40%O2_bliss_22RV1')

    drug_list = cl_df['Drug'].tolist()

    d2t_dict = defaultdict(list)
    for drug in drug_list:
        if '(' in drug:
            drug_split = drug.split(' ')
            drug = drug_split[0]
        
        if drug in dint:
            for n in dint[drug]:
                if 'target' in drugbank_df[(drugbank_df['DrugName'] == drug) & 
                                                       (drugbank_df['geneSym'] == n)].category.tolist():
                    
                    d2t_dict[drug].append(n)
    
    print(len(d2t_dict.keys())) #23/90 drugs in DrugBank (2020 version...)

    
        



if __name__ == "__main__":
    import os, pickle
    import pandas as pd
    import numpy as np
    from collections import defaultdict
    main()