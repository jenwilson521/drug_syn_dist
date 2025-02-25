"""

02/19/2025 - plotting cross-validation results from multi-target curve fitting

"""
def train_vs_test_scatter(corr_dict_train, corr_dict_test, corr_metric, cl):
    fig = plt.figure(dpi = 300)
    ax = plt.subplot()
    fig.tight_layout() 

    fold_colors = ['red', 'goldenrod', 'forestgreen', 'dodgerblue', 'darkviolet']
    x = np.arange(0, 38)
    print(len(x))

    loc = 1
    for metric in corr_dict_train.keys():
        for alphaD in corr_dict_train[metric].keys():
            corr_train_list = list(zip(*corr_dict_train[metric][alphaD]))[0] #p-vals are also saved here...
            corr_test_list = list(zip(*corr_dict_test[metric][alphaD]))[0] #p-vals are also saved here...


            for i, (corr_train, corr_test) in enumerate(zip(corr_train_list, corr_test_list)):
                x_rand = np.random.normal(x[loc], 0.05, size = 2)
                plt.scatter(x_rand[0], corr_train, c = fold_colors[i], marker = 'o', alpha = 0.5, s = 10)
                plt.scatter(x_rand[1], corr_test, c = fold_colors[i], marker = '*', alpha = 0.5, s = 10)
            
            x_rand = np.random.normal(x[loc], 0.05, size = 2)
            plt.scatter(x_rand[0], np.mean(corr_train_list), c = 'black', marker = 'o', alpha = 0.5, s = 10)
            plt.scatter(x_rand[1], np.mean(corr_test_list), c = 'black', marker = '*', alpha = 0.5, s = 10)

            loc += 1
    
    plt.ylabel(corr_metric)

    for axis in ['bottom','left']:
        ax.spines[axis].set_linewidth(1.5)
    for axis in ['top','right']:
        ax.spines[axis].set_linewidth(0)  
    
    ax.tick_params(axis = 'y', width = 1.5)
    ax.tick_params(axis = 'x', width = 0)
    plt.xticks(x, labels = [])

    figure_name = cl + '_cv_train_vs_test_' + corr_metric + '.png'
    plt.savefig(os.path.join(figure_path, figure_name), bbox_inches = 'tight', pad_inches = 0.1)






def main():

    global figure_path
    figure_path = '../../results/figures'

    interactome_source = 'PathFX'
    threshold = 0.50
    interactome_aname = interactome_source + '_' + str(threshold) + 'thr_int'
    exper_source = 'DREAM'

    model_name = 'posyn'

    outpath = os.path.join('../../results/Multi_Targ_Curve_Fitting', exper_source, interactome_aname)

    post_fit_pcc_train = pickle.load(open(os.path.join(outpath, 'post_fit_pcc_and_sig_dict_train_' + model_name + '.pkl'), 'rb'))
    post_fit_spearman_train = pickle.load(open(os.path.join(outpath, 'post_fit_spearman_and_sig_dict_train_' + model_name + '.pkl'), 'rb'))
    post_fit_pcc_test = pickle.load(open(os.path.join(outpath, 'post_fit_pcc_and_sig_dict_test_' + model_name + '.pkl'), 'rb'))
    post_fit_spearman_test = pickle.load(open(os.path.join(outpath, 'post_fit_spearman_and_sig_dict_test_' + model_name + '.pkl'), 'rb'))

    for cl in post_fit_pcc_train.keys():
        train_vs_test_scatter(post_fit_pcc_train[cl], post_fit_pcc_test[cl], 'Spearman', cl)


if __name__ == "__main__":
    import os, pickle
    import matplotlib.pyplot as plt
    import numpy as np
    main()