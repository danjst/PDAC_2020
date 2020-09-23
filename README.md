# Supplementary Code and Scripts

Contains code associated with the manuscript "Intron retention is a robust marker of intertuomral heterogeneity in pancreatic ductal adenocarcinoma" (2020) Tan et al.

1. NMF_code.R : R workflow for running Non-negative matrix factorization (NMF) algorithms on 76 high purity PDAC tumor samples. In addition to clustering samples and extracting top NMF events. 
2. Metrics_Formulas.R: R script that contains formulas for calcuating the following clustering metrics: RMSSTD, R-squared, and SD Index. 
3. rbp_and_event_correlation.R: R script which contains workflow for identifying highly correlated rbp and intron-retention event pairings from PDAC, PRAD, and KIRC data. 
4. IR_spladder_76.sh: Shell script for running quantification of intron-retention events for 76 high-purity PDAC samples using Spladder.
5. IR_spladder_DS.sh: Shell script for running differential splicing (intron retention) analysis of IR-1 vs IR-2 clusters using Spladder. 
6. TF_high_cor_rbps_script.R: R code for identifying significant, highly correlated RBP targets for each TF reported in Chea3 enrichment analysis. 
7. coding_region_match.R: R code for identifying region of IR event occurence. Categorizes IR event as among coding region, 5' UTR, or 3' UTR regions. 
8. nmf_and_survival_analysis_20_cancers.R: R workflow for analyzing Spladder IR quantification data of 20 GDC cancers. Includes NMF clustering and cluster Kaplan-Meier survival analysis.
9.ptc_finder.py: Python script that takes a .fasta file (e.g. sequences of IR events) and reports any detected premature termination codons in three ORFs (open reading frames).
10. IR_DS_filtering_workflow.R: R workflow for filtering Spladder differential intron retention results for events meeting significance and non-null value thresholds. 
