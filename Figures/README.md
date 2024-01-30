# Figures in Topological organization of the inter-chromosomal architecture across cell types

## Figure 1E: K-means clustering
One scripts required: 1E_K-means_clustering.R
### 1E_K-means_clustering.R
Arguments: 
1. Pathway to input (cis and trans interaction frequency)
2. Pathway tO output
## Figure 1H: Validation tickplots (edit by Jordan)
Two scripts required: part 1 (signature_validation_data_processing.R) and part 2 (signature_validation_graphs.R)
### signature_validation_data_processing.R 
Arguments: 
1. Signature zscore 1vsAll output data (ALL interactions)
2. Signature pvalue 1vsAll output data (ALL interactions)
3. Output pathway for RDS object
4. Type of Signature analysis
### signature_validation_graphs.R
Arguments: 
1. RDS output from signature_validation_data_processing.R
2. Text label given that descibes the two interacting loci
3. Pvalue cutoff (numeric)
4. First interacting region bed file (tab separated)
5. Second interacting region bed file (tab separated)
6. A text file (tab separated) formatted like the example below, "Dataset" and "Tissue_Identifier" columns are required, rest are optional
   
    | Gonosomoal_Sex | Dataset               | Tissue_Identifier      | Tissue_Group | Germ_Layer | Group_ID |
    |----------------|-----------------------|------------------------|--------------|------------|----------|
    | F              | Astrocyte_Cerebellum  | b_Astrocyte_Cerebellum | Brain        | Ectoderm   | 1        |
    | M              | Germinal_center_Bcell | bl_GCBC                | Blood        | Mesoderm   | 2        |
    | F              | IMR90_Dixon           | lu_IMR90_D             | Lung         | Endoderm   | 3        |
    
8. Output pathway


## Figure 2A: Genome-wide contact matrix
Five scripts required: pt1_genomewide_matrix_zscoreAll.R, pt2_genomewide_matrix_zscoreSig.R, pt3_genomewide_matrix_zscoreAll_split.R, pt4_genomewide_matrix_zscoreSig_split.R, and pt5_genomewide_matrix_graph.R
### pt1_genomewide_matrix_zscoreAll.R
Arguments:
1. Signature zscore 1vsAll output data (ALL interactions)
2. Output pathway for RDS object
### pt2_genomewide_matrix_zscoreSig.R
Arguments:
1. Signature zscore 1vsAll output data (ALL interactions)
2. Output pathway for RDS object
### pt3_genomewide_matrix_zscoreAll_split.R
Arguments:
1. RDS object from pt1
2. Output pathway for RDS object
3. Text file containing all female (XX) dataset names on a new line
4. Text file containing all male (XY) dataset names on a new line
### pt4_genomewide_matrix_zscoreSig_split.R
Arguments:
1. Signature zscore 1vsAll output data (ALL interactions)
2. Output pathway for RDS object
3. Text file containing all female (XX) dataset names on a new line
4. Text file containing all male (XY) dataset names on a new line
### pt5_genomewide_matrix_graph.R
Arguments:
1. RDS object from pt1 (zscore_df.rds)
2. RDS object from pt3 (zscore_df_X.rds)
3. RDS object from pt3 (zscore_df_Y.rds)
4. RDS object from pt2 (zscore_df_sig.rds)
5. RDS object from pt4 (zscore_df_sig_X.rds)
6. RDS object from pt4 (zscore_df_sig_Y.rds)
7. Output pathway for PDF


## Figure 2B: PQ analysis
Two scripts required: pt1_cooler_data_column_sum.R, and pt2_pq_arms_analysis.R
### pt1_cooler_data_column_sum.R
Arguments:
1. Cooler (1vsAll) output data
2. Output pathway for txt file
### pt2_2D_p_q_arms_analysis.R
Arguments:
1. Input pathway to the txt file from pt1
2. Output pathway to the results (txt file including interaction types, binomial test results, and figure)
## Figure 2E: PQ distribution
One script required
### PQ_Interaction_Distrubtion.R
Arguments:
1. Signature qvalue 1vsAll output data (positive significant interactions)
2. Output pathway for PDF


## Figure 2F/G: PQ alluvials
2 scripts required: prepare_data_pqarm.R and signature_alluvial_pqarm.R
### prepare_data_pqarm.R
Arguments:
1. Signature qvalue 1vsAll output data (positive significant interactions)
2. Signature qvalue 1vsAll output data (ALL interactions)
3. Text file containing all female (XX) dataset names on a new line
4. Text file containing all male (XY) dataset names on a new line
5. Output pathway for RData
6. Supplemental table S4-S8 (i.e. P/Q arm annotation)
### signature_alluvial_pqarm.R
Arguments:
1. Pathway to Rdata output from prepare_data_pqarm.R


## Figure 3C: AB analysis
Three scripts required: pt1_compartment.R, pt2_compartment_parallel.R, and pt3_merge_visualization.R
### pt1_data_prepration.R
Arguments:
1. Pathway to Signature zscore 1vsAll data
2. Output pathway
### pt2_compartment_parallel.R
Arguments:
1. Index of cell types ($1 in the shell script)
2. Output pathway (same pathway as the pt1_data_prepration.R)
3. Pathway to aggregated compartment table for 62 cell types (Table S13)
###  pt3_merge_visualization.R   
Arguments:
1. Pathway to the cell types (txt files)
2. Pathway to the output (merged file)
