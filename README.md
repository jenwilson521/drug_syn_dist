# drug_syn_dist

This suite of code consists executes our network model of protein target combinations within the interactome. 

Can be applied to any combination screen as long as protein targets and experimental readouts (i.e., synergy scores) are available.

## Input data

We use two example datasets of protein target combinations: a gene interaction (GI) screen and small molecule combination screen. 

With access to the GI screen data, all analyses with the exception of those within the multi_targ_analysis folder can be executed. 

We provide the code to read and save this data in a format compatible with all other scripts/read_input_data/read_GI_screen_data.py.

## Target neighborhoods

We used https://biit.cs.ut.ee/gprofiler/convert to convert our interactome's gene symbols to gene IDs in order to be compatible with GOATOOL's enrichment analysis. Further, we lifted code from the GOATOOLs repository, https://github.com/tanghaibao/goatools/tree/main, in order to implement enrichment analysis. 

## Multi-target analysis

The framework is currently configured for first-degree polynomial fitting, however, different curves (e.g., higher-order polynomials, exponential functions, etc) can easily implemented by changing the function called in line 227.
