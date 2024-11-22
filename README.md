# drug_syn_dist

This suite of code executes our network model of protein target combinations within the interactome. The model only requires a weighted, undirected interactome text file, which we provide in data/0.5thr_filtered_PathFX_scored_interactome, and combined perturbed protein targets with their corresponding experimental readout out of combination effects. Therefore, while we use two example data sources, our model can be applied to any combination screen as long as protein targets and synergy scores (or a similar readout) are available.

## Input data

We use two example datasets of protein target combinations - a gene interaction (GI) screen and small molecule combination screen. 

The GI screen data are publicly available at https://www.nature.com/articles/nbt.3834, Supplemental Table 4. With access to the GI screen data, all analyses with the exception of those within the multi_targ_analysis folder can be executed. We provide the code to read and save this data in a format compatible with all other scripts in scripts/read_input_data/read_GI_screen_data.py.

The small molecule combination screen comes from a DREAM challenge and is available upon application at the portal: https://www.synapse.org/Synapse:syn18496666. Once given access, the code provided in scripts/read_input_data/read_DREAM_data.py reads and saves this data in a format compatible with the rest of the model. Select small molecules in this dataset had an identifiable name in DrugBank. We provide the necessary files to extract their corresponding protein targets at resources/Drugbank050120.xlsx and resources/Pfx050120_dint.pkl.

## Target neighborhoods

To define protein target neighborhoods (i.e., clusters of neighboring proteins surrounding a target), we provide code that runs the local community detection algorithm, PageRank-Nibble (PRN), at scripts/neighborhoods/run_PageRank_Nibble.py. The script calls the function myGraph (scripts/neighborhoods/myGraph.py) which we pulled from the publicly available repository, https://github.com/Jx-Liu-229/PageRank-Nibble. Additionally, within run_PageRank_Nibble.py, the restart probability hyperparameter(s) can be changed in Line 94; the script saves neighborhood results for all hyperparameters in a pickled dictionary.

To determine the hyperparameter that yields an optimal target neighborhood (we provide examples for searching for those that are large and compact), we used the script script/neighborhoods/find_opt_prn_alpha.py. This script resaves neighborhood results using each target's optimal available (pickled dictionary) and also generates figures to compare these results to those for neighborhoods defined by standard means (i.e., taking all first-degree target-neighbor interactions).

To perform Gene Ontology Enrichment Analysis (GOEA), was used the python-based library GOATOOLs. This analysis is executed for both optimal hyperparameter PRN and first-neighbors neighborhoods in scripts/targ_nbhd_goea.py. Much of the code used in this script was lifted from the GOATOOLs repository, https://github.com/tanghaibao/goatools/tree/main. However, users can easily perform this analysis for non-optimal PRN neighborhoods or opt out of running it for first-neighbors neighborhoods by commenting out Line 197. Notably, to be compatible with GOATOOLs, we converted our interactome's gene symbols (i.e., background list used in GOEA) to gene IDs using https://biit.cs.ut.ee/gprofiler/convert. We provide these gene IDs at data/0.5thr_PathFX_int_GeneSymbol2GeneID.csv.

## Interactome distances

To generate the main inputs for our interactome distances, we generated diffusion profiles at scripts/distances/compute_diffusion_profiles.py. Diffusion profiles – representations of single target effects in the interactome – are computed as the steady-state distriubtion of visitation probabilities from the random walk with restart (RWR) algorithm. We provide examples for running the RWR algorithm relative to varying numbers of starting nodes as well as at 6 different diffusion depth parameters. We save all such results as dictionaries in .npz files.

The diffusion profiles are then used as inputs to run scripts/distances/calc_distance_metrics.py. For each target pair (or target set pair for multi-target combinations), we provide functions to compute a standard average shortest paths metric as well as for diffusion-based distances that leverage various configurations of the diffusion profiles for all diffusion depth parameters. We save the results as separate pickled dictionaries, however, we also provide code to condense all distance results into a single, nested dictionary at scripts/distances/rearrange_distance_metric_dicts.py.

## Multi-target analysis

We provide analyses for two representations of multi-target combinations (i.e., those provided in the DREAM challenge data). The first representation involves merging all target neighborhoods within a multi-target drug (scripts/multi_targ_analysis/merge_multi_targ_nbhds.py) and computing interactome distances relative to these target sets and merged neighborhoods, as described above.

The second representation involves considering all target pair distandes for a given combination of target sets and analyzing which, if any, are better predictors of their experimental readout. This analysis is executed in scripts/multi_targ_analysis/multi_targ_curve_fitting.py. The framework is currently configured for first-degree polynomial fitting, however, different curves (e.g., higher-order polynomials, exponential functions, etc) can easily implemented by changing the function called in line 227. Within this script, we also provide example randomizations to explore the robustness of fit, such as initializing the models with random distances as well as running the framework on shuffled experimental data. 

To evaluate curve fitting performance for all distance metrics, we provide the script script/mutli_targ_analysis/evaluate_curve_fitting.py which computes and plots a normalized performance score and performs enrichment analysis for top performing distance metrics and their encoded network features. Results from this enrichment analysis are saved in a pickled dictionary.

## Network feature sensitivity analysis

Our diffusion-based distances were all encoded by differing assumptions of three network features (diffusion depth, target neighborhood inclusion, and directionality). We provide a sensitivity analysis to determine the importance of each feature in scripts/net_feat_sens_analysis.py. The analysis is split into comparisons of feature sets that only vary in a single parameter and sensitiivty trends are plotted. Here we provide examples for comparing distance-experiment Pearson's correlation coefficients in order to determine feature sensitivities, however, this metric can easily be changed given distance-experiment relationships (i.e., comparisons of mean squared error may be better suited for higher-order relationships). Notably, this analysis can be executed for both cases of single or multi-target combinations.
