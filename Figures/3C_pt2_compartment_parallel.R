#_____Load required packages_____________________________________________________

library(dplyr)
library(tidyr)
library(parallel) 

#_____Read in arguments_________________________________________________________

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

outpath <- args[2]
compartment_pathway <- args[3]

#_____Input list_________________________________________________________
cells <- c("Adrenal_gland_Schmitt.txt","Aorta_Leung.txt",
"Astrocyte_Cerebellum.txt",
"Astrocyte_Spine.txt",
"Bladder_Schmitt.txt",
"Cardiac_mesoderm_cell_day05_Zhang.txt",
"Cardiac_progenitor_cell_day07_Zhang.txt",
"Chondrocyte_day01_Maass.txt",
"Chondrocyte_day03_Maass.txt",
"Chondrocyte_day07_Maass.txt",
"Chondrocyte_day14_Maass.txt",
"Chondrocyte_day21_Maass.txt",
"Dorsolateral_prefrontal_cortex_Schmitt.txt",
"Endothelial_cell_Microvasculature.txt",
"Germinal_center_Bcell.txt",
"Glia_Rajarajan.txt",
"H1hESC_Dixon.txt",
"H1hESC_Oksuz.txt",
"H9hESC_day00_Zhang.txt",
"Hippocampus_Schmitt.txt",
"HMEC_Rao.txt",
"hTERT_HPNE_Ren.txt",
"HUVEC_Rao.txt",
"IMR90_Dixon.txt",
"IMR90_Rao.txt",
"Islet_Lawlor.txt",
"Keratinocyte_day0_Rubin.txt",
"Left_ventricle_Leung.txt",
"Liver_Leung.txt",
"Lung_PatientF_Schmitt.txt",
"Lung_PatientM_Schmitt.txt",
"Memory_Bcell.txt",
"Mesenchymal_stem_cell_day00_Maass.txt",
"Mesenchymal_stem_cell_Dixon.txt",
"Mesendoderm_cell_Dixon.txt",
"Mesoderm_cell_day02_Zhang.txt",
"Naive_Bcell.txt",
"Neuron_Rajarajan.txt",
"NHEK_dilution_Rao.txt",
"NHEK_insitu_Rao.txt",
"NP69_Animesh.txt",
"NPC_Dixon.txt",
"NPC_Rajarajan.txt",
"Ovary_Schmitt.txt",
"Pancreas_PatientF_Schmitt.txt",
"Pancreas_PatientM_Schmitt.txt",
"PHH_patient342_NonInfect_Moreau.txt",
"PHH_patient399_NonInfect_Moreau.txt",
"Plasma_Bcell.txt",
"Primitive_cardiomyocyte_day15_Zhang.txt",
"Psoas_Schmitt.txt",
"Right_ventricle_Schmitt.txt",
"RPE1_Control_Darrow.txt",
"Small_bowel_Schmitt.txt",
"Spleen_Schmitt.txt",
"Tcell_CD4_Kloetgen.txt",
"TeloHAEC_TNFalphaControl_0h_Lalonde.txt",
"TeloHAEC_Wilson.txt",
"Thymus_Leung.txt",
"Trophectoderm_cell_Dixon.txt",
"Ventricular_cardiomyocyte_day80_Zhang.txt",
"VSMC_day21_Maass.txt")

#______________________Add compartment_A and compartment_B and total compartment columns__________________________ 

cell <- cells[as.numeric(args[[1]])]
cat(sprintf("cell = %s",cell), sep="\n")


  txt_df <- read.table(sprintf("%s/%s",outpath, cell), header = TRUE)
  compartment <- read.csv(sprintf("%s", compartment_pathway), header=T)

  compartment_col_name <- gsub(".txt", "", cell)
  
  txt_df$compartment_A <- NA
  txt_df$compartment_B <- NA
  
  for (i in 1:nrow(txt_df)) {
    
    indexA <- which(compartment$chrom == sub("^A", "", txt_df[i, "chrA"]) & 
                      compartment$start == txt_df[i, "st1"])
    
    if (length(indexA) > 0) {
      txt_df[i, "compartment_A"] <- compartment[indexA, compartment_col_name]
    }
    indexB <- which(compartment$chrom == sub("^B", "", txt_df[i, "chrB"]) & 
                      compartment$start == txt_df[i, "st2"])
    
    if (length(indexB) > 0) {
      txt_df[i, "compartment_B"] <- compartment[indexB, compartment_col_name]
    }
  }
  

  txt_df$final_compartment <- ifelse(txt_df$compartment_A > 0 & txt_df$compartment_B > 0, "AA",
                                    ifelse(txt_df$compartment_A <= 0 & txt_df$compartment_B <= 0, "BB", "AB"))
  
  write.table(txt_df, file = cell, sep = "\t", quote = FALSE, row.names = FALSE)

