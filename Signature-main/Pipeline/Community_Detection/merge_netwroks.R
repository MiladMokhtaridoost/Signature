library(dplyr) 
library(purrr)

DATA_PATH <- [pathway/to/your/netwrok/folder]
RESULT_PATH <- [pathway/to/your/output/folder]

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
           "Keratinocyte_day0_Rubin",
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
           "Pancreas_PatientM_Schmitt",
           "PHH_patient342_NonInfect_Moreau",
           "PHH_patient399_NonInfect_Moreau",
           "Plasma_Bcell",
           "Primitive_cardiomyocyte_day15_Zhang",
           "Psoas_Schmitt",
           "Right_ventricle_Schmitt",
           "RPE1_Control_Darrow",
           "Small_bowel_Schmitt",
           "Spleen_Schmitt",
           "Tcell_CD4_Kloetgen",
           "TeloHAEC_TNFalphaControl_0h_Lalonde",
           "TeloHAEC_Wilson",
           "Thymus_Leung",
           "Trophectoderm_cell_Dixon",
           "Ventricular_cardiomyocyte_day80_Zhang",
           "VSMC_day21_Maass")

read_function <- function(i){
  
  read.table(sprintf("%s/%s_network.txt", DATA_PATH, cells[i]), header = T)
  
}

print("reading in files ...")
l_dfs <- lapply(1:length(cells), read_function)

# merge all batches together
print("merging ...")
mdat <- reduce(l_dfs, full_join, by = c("ID_chrA","ID_chrB"))

for(i in 3:ncol(mdat)){
  mdat[,i] <- as.numeric(as.character(mdat[,i])) 
}
  
colnames(mdat) <- c("ID_chrA", "ID_chrB", cells)

mdat$average <- rowMeans(mdat[, -(1:2)], na.rm = TRUE)

avrg <- mdat[,c("ID_chrA", "ID_chrB", "average")]

write.table(avrg, file = sprintf("%s/average_1MB_network.txt", RESULT_PATH), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

