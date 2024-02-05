library(dplyr)
library(tidyr)

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
DATA_PATH <- args[2]
RESULT_PATH <- args[3]
DATA_PATH <- "Z:/3D-flow/normalized_data_4DNuc_pipeline/human/Test"
RESULT_PATH <- "Z:/Milad/Community_detection/GitHub"

#resolution <- as.numeric(args[[4]])
resolution <- 1000000
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

cells <- c("Astrocyte_Spine",
           "H9hESC_day00_Zhang")

CELL <- cells[as.numeric(args[[1]])]
CELL <- cells[2]
###DATA_PATH <- "Z:/3D-flow/normalized_data_4DNuc_pipeline/human/CHLA9_Maass"

#int_table <- read.table(sprintf("%s/%s/trans.100000_iced.sorted.txt", DATA_PATH, CELL))
#int_table <- read.table(sprintf("%s/%s/trans.100000_iced.sorted.txt", DATA_PATH, CELL))
int_table <- read.delim(sprintf("%s/%s/%s.pairs.res1000000.cool.txt",DATA_PATH,CELL,CELL), header = F)

cat(sprintf("Cell type = %s", CELL), sep="\n")
colnames(int_table) <- c("chrA", "stA", "endA","chrB", "stB", "endB", "rawread", "freq")

int_table <- int_table[!is.na(int_table$freq),]
int_table <- int_table[,c("chrA", "stA", "chrB", "stB", "freq")]

########deleting rows correspond to chrM
int_table <- int_table[!(int_table$chrA == "chrM" | int_table$chrB == "chrM"), ]

int_table$chrA <- sub("chr", "", as.character(int_table$chrA))
int_table$chrB <- sub("chr", "", as.character(int_table$chrB))

resolution <- as.numeric(resolution)

int_table$stA <- as.numeric(int_table$stA)/resolution 
int_table$stB <- as.numeric(int_table$stB)/resolution

int_table$ID_chrA <- paste(int_table$chrA, int_table$stA, sep = "_")
int_table$ID_chrB <- paste(int_table$chrB, int_table$stB, sep = "_")

int_table$freq <- as.numeric(int_table$freq)

int_table <- int_table[,c("ID_chrA", "ID_chrB", "freq")]

#########################################################
#write.csv(int_table, sprintf("%s/%s_network.csv", RESULT_PATH, CELL))
write.table(int_table, sprintf("%s/%s_network.txt", RESULT_PATH, CELL), row.names = F)
