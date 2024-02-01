#_____Load required packages_____________________________________________________
library(tidyr)
library(parallel)
library(dplyr)
library(parallel)
library(data.table)
library(ggplot2)
library(RColorBrewer)

#_____Read in arguments_________________________________________________________
args = commandArgs(trailingOnly = TRUE)

input_pathway <- args[1]
save_output_path <- args[2]

#__________________________Set the theme__________________________________________________
theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            theme(plot.title = element_text(hjust = 0.6, size = 14),
                  panel.background = element_rect(fill = "white", colour = NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  #panel.spacing = unit(0.25, "lines"),
                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                  #change the colour of facet label background
                  #strip.background = element_rect(fill = "#E6E1EA"),
                  #remove space between facets
                  # panel.spacing.x=unit(0, "lines"),
                  #legend.key = element_blank(),
                  legend.background=element_blank(),
                  #legend background
                  legend.key = element_rect(fill = NA),
                  legend.title = element_text(vjust = 0.8),
                  legend.position="right")
)
#_____________________________Input___________________________________________

cells <- c("Adrenal_gland_Schmitt.txt",
           "Aorta_Leung.txt",
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

#_____________________________Merging all files_______________________________________________

merged_compartment <- data.frame()

for (i in seq_along(cells))
{  
  
  file_name <- cells[i] 
  
  full_path <- file.path(input_pathway, file_name) 
  
  temp_all_compartments <- read.table(full_path, header = FALSE) 
  
  if (i == 1) {   
    
    merged_compartment <- temp_all_compartments  
    
  } else {   
    
    merged_compartment <- rbind(merged_compartment, temp_all_compartments[-1, ])   
  } 
}

print(merged_compartment) 
 

write.table(merged_compartment, sprintf("%s/merged_compartment.txt", save_output_path), sep = "\t", quote = FALSE, row.names = FALSE)



#_____________________________visualization__________________________________________

colnames(merged_compartment) <- merged_compartment[1,]
merged_compartment<-merged_compartment[-1,]
merged_compartment$zscore <- as.numeric(as.character(merged_compartment$zscore))
merged_compartment$zscore<- round(merged_compartment$zscore,digits=1)

#statistics
AA<- merged_compartment[merged_compartment$final_compartment == "AA",]

BB<- merged_compartment[merged_compartment$final_compartment == "BB",]

AB<- merged_compartment[merged_compartment$final_compartment == "AB",]



AA_BB_test<-wilcox.test(AA$zscore,BB$zscore)
AA_AB_test<-wilcox.test(AA$zscore,AB$zscore)
AB_BB_test<-wilcox.test(AB$zscore,BB$zscore)

pvalues_1<- AA_BB_test$p.value
pvalues_2<- AA_AB_test$p.value
pvalues_3<- AB_BB_test$p.value

#all pvalues are less than 0.001


mean_AA<- mean(AA$zscore)
median_AA <- median(AA$zscore)
min_AA <- min(AA$zscore)
min_AA
max_AA <- max(AA$zscore)
max_AA
first_quartile_AA <- quantile(AA$zscore,0.25)
third_quartile_AA <- quantile(AA$zscore,0.75)


mean_BB<- mean(BB$zscore)
median_BB <- median(BB$zscore)
min_BB <- min(BB$zscore)
min_BB
max_BB <- max(BB$zscore)
max_BB
first_quartile_BB <- quantile(BB$zscore,0.25)
first_quartile_BB
third_quartile_BB <- quantile(BB$zscore,0.75)
third_quartile_BB



mean_AB<- mean(AB$zscore)
mean_AB
median_AB <- median(AB$zscore)
median_AB
min_AB <- min(AB$zscore)
min_AB
max_AB <- max(AB$zscore)
max_AB
first_quartile_AB <- quantile(AB$zscore,0.25)
first_quartile_AB
third_quartile_AB <- quantile(AB$zscore,0.75)
third_quartile_AB


final_compartment<- c("AA","AB","BB")
mean<-c(mean_AA,mean_AB,mean_BB)
median<-c(median_AA,median_AB,median_BB)
first<-c(first_quartile_AA,first_quartile_AB,first_quartile_BB)
third<-c(third_quartile_AA,third_quartile_AB,third_quartile_BB)
third
min<-c(min_AA,min_AB,min_BB)
min
max<-c(max_AA,max_AB,max_BB)
max
stat<- data.frame(final_compartment,mean,median,first,third)
stat
data_long <- stat %>%
  pivot_longer(cols = c(mean, median, first, third),
               names_to = "Statistic",
               values_to = "zscore")
data_long

data_long<- data.frame(data_long)
row_names <-rownames(data_long[,1])
row_names

  
count_AA<- sum(merged_compartment$final_compartment== "AA")
count_AA
count_BB<- sum(merged_compartment$final_compartment== "BB")
count_BB

count_AB<- sum(merged_compartment$final_compartment== "AB")
count_AB
 df <- data.frame(   Compartment = rep(c("AA", "AB", "BB"), each = 3),   Statistic = rep(c("median", "first", "third"), times = 3),   ZScore = c(2.3, 2.1, 2.7, -2.3, -2.6, -2.1, -2.2, -2.5, -2) ) 

  colors <- c("AA" = "pink", "AB" = "green", "BB" = "cyan")
 {plot <- ggplot(df, aes(x = Compartment, y = ZScore, fill = Compartment)) +  
     geom_boxplot(coef = 0, outlier.shape = NA) +   labs(x = "Compartment", y = "Significant Z-Score") + 
     scale_fill_manual(values = colors) +   ggtitle("Boxplot of Z-Scores by Compartment") +  
     theme(legend.position = "bottom")
   
   filename <- "compartment_comparison_zscore_box_plot.pdf" 
   pdf(filename, width = 12, height = 16)
   print(plot)
   dev.off() 
   print("done")}

  
  

