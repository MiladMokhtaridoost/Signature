# Figures in "Inter-chromosomal contacts demarcate genome topology along a spatial gradient"
<br/>

## Figure 1F: Validation tickplots
Two scripts required: part 1 (signature_validation_data_processing.R) and part 2 (signature_validation_graphs.R)
### signature_validation_data_processing.R 
Arguments: 
1. Signature zscore output data (ALL interactions)
2. Signature pvalue output data (ALL interactions)
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
    
8. Output pathway for the results
<br/>

## Figure 2A: Genome-wide contact matrix
Five scripts required: pt1_genomewide_matrix_zscoreAll.R, pt2_genomewide_matrix_zscoreSig.R, pt3_genomewide_matrix_zscoreAll_split.R, pt4_genomewide_matrix_zscoreSig_split.R, and pt5_genomewide_matrix_graph.R
### pt1_genomewide_matrix_zscoreAll.R
Arguments:
1. Signature zscore output data (ALL interactions)
2. Output pathway for RDS object
### pt2_genomewide_matrix_zscoreSig.R
Arguments:
1. Signature zscore output data (ALL interactions)
2. Output pathway for RDS object
### pt3_genomewide_matrix_zscoreAll_split.R
Arguments:
1. RDS object from pt1
2. Output pathway for RDS object
3. Text file containing all female (XX) dataset names on a new line
4. Text file containing all male (XY) dataset names on a new line
### pt4_genomewide_matrix_zscoreSig_split.R
Arguments:
1. Signature zscore output data (ALL interactions)
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
<br/>

## Figure 2C: AB compartment analysis
Three scripts required: 2C_pt1_data_prepration.R, 2C_pt2_compartment_parallel.R, and 2C_pt3_merge_visualization.R
### 2C_pt1_data_prepration.R
Arguments:
1. Signature zscore output data (significant interactions)
2. Output pathway for the results
### 2C_pt2_compartment_parallel.R
Arguments:
1. Index of cell types ($1 in the shell schedular script)
2. Output pathway (argument 2 from 2C_pt1_data_prepration.R)
3. Aggregated compartment table for all cell types (Table S4)
###  2C_pt3_merge_visualization.R   
Arguments:
1. Output pathway (argument 2 from 2C_pt1_data_prepration.R)
<br/>

## Figure 2E: PQ arm interaction analysis
One script required: 2E_pq_arms_analysis.R
### 2E_pq_arms_analysis.R
Arguments:
1. Cooler data
2. Output pathway for the results
<br/>

## Figure 2F: PQ arm interaction distribution
One script required: 2F_pq_interaction_distrubtion.R
### 2F_pq_interaction_distrubtion.R
Arguments:
1. Signature qvalue output data (positive significant interactions)
2. Output pathway for PDF
<br/>

## Figure 2I: Telomeric and centromeric analysis
One script required: 2I_centromeric_telomeric_interactions.R (note: the pre-processing of data is same for figure 2J)
### 2I_centromeric_telomeric_interactions.R
Arguments:
1. Signature qvalue output data (positive significant interactions)
<br/>

## Figure 4A: K-means clustering
One scripts required: 4A_Kmeans_clustering.R
### 4A_Kmeans_clustering.R
Arguments: 
1. Tab-separated text file containing columns of i) cell types, ii) summed cooler data in cis, iii) summed cooler data in trans

    | Cell                    | Cis_IntFreq   | Trans_IntFreq   |
    |-------------------------|---------------|-----------------|
    | Adrenal_gland_Schmitt   | 1946.46       | 792.381         |
    | Astrocyte_Cerebellum    | 1867.59       | 765.822         |
    | Astrocyte_Spine         | 1554.89       | 542.854         |
   
3. Output pathway for results
<br/>

## Figure 5A/B: PQ alluvials
2 scripts required: 5AB_pt1_prepare_data_pqarm.R and 5AB_pt2_alluvial_pqarm.R
### 5AB_pt1_prepare_data_pqarm.R
Arguments:
1. Signature qvalue output data (positive significant interactions)
2. Signature qvalue output data (ALL interactions)
3. Text file containing all female (XX) dataset names on a new line
4. Text file containing all male (XY) dataset names on a new line
5. Output pathway for RData
6. P/Q arm annotations of Signature data
### 5AB_pt2_alluvial_pqarm.R
Arguments:
1. Pathway to Rdata output (argument 5 from 5AB_pt1_prepare_data_pqarm.R)

## Extended Figure 3C: Merfish comparison
One script required: Ex.3C_Merfish_data_analysis.R
### Ex.3C_Merfish_data_analysis.R
Arguments:
1. Signature qvalue output data (all interactions)
2. Signature qvalue output data (positive significant interactions)
3. Merfish genomic scale data (PMID: 25858977)

<br/>


## Extended Figure 3D: HiCAN comparison
One script required: Ex.3D_HiCAN_overlap.R
### Ex.3D_HiCAN_overlap.R
Arguments: 
1. HiCAN's top 100 nucleolar loci (Figure 2B PMID: 37013988)
2. HiCAN's top 100 speckle loci (Figure 2B PMID: 370139881)
3. Signature qvalue output data (positive significant interactions) [GM12878]
4. Output pathway for the results
5. Tab-separated text file containing all of Signature's 1Mb bins and their annotated genes (genes comma separated)
<br/>
 
## Extended Figure 3E: SPRITE comparison
One script required: Ex.3E_Sprite_overlap.R
### Ex.3E_Sprite_overlap.R
Arguments: 
1. SPRITE's nucleolar and active hubs (Table S2 PMID: 29887377)
2. Signature qvalue output data (positive significant interactions) [GM12878]
3. Output pathway for the results
4. Tab-separated text file containing all of Signature's 1Mb bins and their annotated genes (genes comma separated)
<br/>

