options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

DATA_PATH <- args[1]
RESULT_PATH <- args[2]

library(data.table)
library(dplyr)
library(tidyr)

cells <- c("Adrenal_gland_Schmitt",
           "Aorta_Leung",
           "Astrocyte_Cerebellum",
           "Astrocyte_Spine",
           "Bladder_Schmitt",
           "Cardiac_mesoderm_cell_day05_Zhang",
           "Cardiac_progenitor_cell_day07_Zhang",
           "Chondrocyte_day01_Maass",
           "Chondrocyte_day03_Maass",
           "Chondrocyte_day07_Maass",
           "Chondrocyte_day14_Maass",
           "Chondrocyte_day21_Maass",
           "Dorsolateral_prefrontal_cortex_Schmitt",
           "Endothelial_cell_Microvasculature",
           "Germinal_center_Bcell",
           "Glia_Rajarajan",
           "H1hESC_Dixon",
           "H1hESC_Oksuz",
           "H9hESC_day00_Zhang",
           "Hippocampus_Schmitt",
           "HMEC_Rao",
           "hTERT_HPNE_Ren",
           "HUVEC_Rao",
           "IMR90_Dixon",
           "IMR90_Rao",
           "Islet_Lawlor",
           #"Keratinocyte_day0_Rubin",
           "Left_ventricle_Leung",
           "Liver_Leung",
           "Lung_PatientF_Schmitt",
           "Lung_PatientM_Schmitt",
           "Memory_Bcell",
           "Mesenchymal_stem_cell_day00_Maass",
           "Mesenchymal_stem_cell_Dixon",
           "Mesendoderm_cell_Dixon",
           "Mesoderm_cell_day02_Zhang",
           "Naive_Bcell",
           "Neuron_Rajarajan",
           "NHEK_dilution_Rao",
           "NHEK_insitu_Rao",
           "NP69_Animesh",
           "NPC_Dixon",
           "NPC_Rajarajan",
           "Ovary_Schmitt",
           "Pancreas_PatientF_Schmitt",
           #"Pancreas_PatientM_Schmitt",
           #"PHH_patient342_NonInfect_Moreau",
           #"PHH_patient399_NonInfect_Moreau",
           "Plasma_Bcell",
           "Primitive_cardiomyocyte_day15_Zhang",
           "Psoas_Schmitt",
           "Right_ventricle_Schmitt",
           "RPE1_Control_Darrow",
           "Small_bowel_Schmitt",
           "Spleen_Schmitt",
           "Tcell_CD4_Kloetgen",
           "TeloHAEC_TNFalphaControl_0h_Lalonde",
           #"TeloHAEC_Wilson",
           "Thymus_Leung",
           "Trophectoderm_cell_Dixon",
           "Ventricular_cardiomyocyte_day80_Zhang",
           "VSMC_day21_Maass")


final_domains <- c()

for (cell in cells) {
  
  cat(sprintf("cell = %s", cell), sep="\n")
  
  ######### significant p-value (interacting regions - positive zscore)
  qval_pos <- fread(sprintf("%s/signature_trans1vsAll_1MB_merged_qvalue_neg.txt", DATA_PATH),
                    select = c("ID" ,sprintf("%s", cell)))
  
  #dd <- read.table("Z:/Signature/results/merged_output/diploid/signature_trans1vsAll_1MB_merged_qvalue_neg.txt")
  qval_pos <- qval_pos[which(qval_pos[,2] != "NA"), ]

  ### changing the format of IDs
  IDs <- data.frame("ID_A","ID_B", "pval")
  for (i in 1:nrow(qval_pos)) {
    IDs[i,1:2] <- unlist(strsplit( gsub("B","B~",qval_pos$ID[i]), "~" ))
    IDs[i,1] <- unlist(gsub("Achr","", IDs[i,1]))
    IDs[i,2] <- unlist(gsub("chr","", IDs[i,2]))
    IDs[i,3] <- qval_pos[i,2]
  }
  
  IDs <- separate(IDs, X.ID_A., into = c("chrA", "stA", "endA"), sep = "\\.")
  IDs <- separate(IDs, X.ID_B., into = c("chrB", "stB", "endB"), sep = "\\.")
  
  IDs <- IDs[,c("chrA", "stA", "chrB", "stB", "X.pval.")]
  colnames(IDs) <- c("chrA", "stA", "chrB", "stB", "pval")
  
  IDs$stA <- as.numeric(IDs$stA)
  IDs$stB <- as.numeric(IDs$stB)
  
  colnames(IDs) <- c("chr", "st", "chr", "st", "pval")
  
  sig_IDs <- rbind(IDs[,c(1,2)], IDs[,c(3,4)])
  sig_IDs$st <- sig_IDs$st/1000000
  
  uniq_sigs <- unique(sig_IDs)
  uniq_chrs <- unique(uniq_sigs$chr)
  
  domain_chr <- c()
  
  for(c in uniq_chrs){
    
    uniq_ch <- uniq_sigs[which(uniq_sigs$chr == sprintf("%s",c)),]
    uniq_ch <- uniq_ch[order(uniq_ch$st),] 
    bins <- uniq_ch[,2]
    
    if(length(bins) > 1){
      
      consecutive_bins <- c()
      current_length <- 1
      
      # Iterate over the vector
      for (j in 2:length(bins)) {
        # Check if the current number is consecutive or has a gap of at most 2
        if (bins[j] == bins[j - 1] + 1 || bins[j] == bins[j - 1] + 2) {
          current_length <- current_length + 1
        } else {
          # Store the length of the previous consecutive sequence
          consecutive_bins <- c(consecutive_bins, bins[j-1] - bins[j-current_length] + 1)
          # Reset the length counter
          current_length <- 1
        }
      }
      
      # Store the length of the last consecutive sequence
      consecutive_bins <- c(consecutive_bins, current_length)
      
    } else {
      consecutive_bins <- 1
    }
    
    domain_chr <- c(domain_chr, consecutive_bins)
    
  }
  
  final_domains <- c(final_domains, consecutive_bins)
  
}


save(final_domains, file = sprintf("%s/result_negative_qval.RData", RESULT_PATH))

