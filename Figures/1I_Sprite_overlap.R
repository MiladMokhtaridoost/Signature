#_____Load required packages_____________________________________________________
options(scipen=999)
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(ggVennDiagram)
  library(VennDiagram)
})


#_____Read in arguments_________________________________________________________
args = commandArgs(trailingOnly = TRUE)

sprite_nhub <- read.xlsx(args[1], sheet = "2D_nucleolar_hub_GM12878")
sprite_ahub <- read.xlsx(args[1], sheet = "2E_active_hub_GM12878")
sign_GM <- read.table(args[2], header = T)
out <- args[3]


#_____Prepare data______________________________________________________________

# split up interaction ID information into new columns
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
sign_GM <- sign_GM %>% select(ID)
sign_GM$ID <- sub("B", "\\.B", as.character(sign_GM$ID))
sign_GM <- sign_GM %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
sign_GM$chrA <- gsub("A", "", sign_GM$chrA)
sign_GM$chrB <- gsub("B", "", sign_GM$chrB)

# prepare list of bins (signature GM)
sign_GM$bin1 <- paste(sign_GM$chrA, sign_GM$st1, sign_GM$end1, sep = ".")
sign_GM$bin2 <- paste(sign_GM$chrB, sign_GM$st2, sign_GM$end2, sep = ".")
sign_GM_list <- c(sign_GM$bin1,sign_GM$bin2)
sign_GM_list <- unique(sign_GM_list)

# prepare list of bins (SPRITE GM-a)
sprite_ahub$bin <- paste(sprite_ahub$chr, sprite_ahub$st, sprite_ahub$end, sep = ".")
sprite_ahub_list <- unique(sprite_ahub$bin)

# prepare list of bins (SPRITE GM--n)
sprite_nhub$bin <- paste(sprite_nhub$chr, sprite_nhub$st, sprite_nhub$end, sep = ".")
sprite_nhub_list <- unique(sprite_nhub$bin)



#_____Produce Venn Diagram______________________________________________________

gglist <- list("sign_GM" = sign_GM_list,
               "SPRITE_nucHub" = sprite_ahub_list,
               "SPRITE_actHub" = sprite_nhub_list)

ggVennDiagram(gglist, set_size = 3, label_alpha = 0, label_size = 4) +
              labs(title = "Signature vs SPRITE (q<0.05) [GM12878]", fill = "1MB loci") +
              scale_fill_distiller(palette = "Reds", direction = 1) +
              scale_color_brewer(palette = "Greys") +
              theme(plot.title = element_text(face="bold", colour="black", size=14, hjust=0.5, vjust=3),
                    plot.subtitle = element_text(face="plain", colour="black", size=12, hjust=0.5),
                    legend.title = element_text(face="plain", colour="black", size=10, hjust=0.5),
                    legend.text = element_text(face="plain", colour="black", size=8),
                    legend.position = NULL)

ggsave(paste0(out,"/sprite_VD.pdf"), width = 8, height = 6)


# get list of unique bins to extract genes later
all_overlaps <- calculate.overlap(gglist) 

total <- lengths(all_overlaps)

ggg <- c(rep("all",total[1]),
         rep("nu+si",total[2]),
         rep("ac+si",total[3]),
         rep("ac+nu",total[4]),
         rep("si",total[5]),
         rep("nu",total[6]),
         rep("ac",total[7]))

bbb <- c(all_overlaps[[1]],
         all_overlaps[[2]],
         all_overlaps[[3]],
         all_overlaps[[4]],
         all_overlaps[[5]],
         all_overlaps[[6]],
         all_overlaps[[7]])

genes_df <- data.frame("bins" = bbb,
                       "group" = ggg)



#_____Annotate with genes_______________________________________________________

dat <- read.table(args[4], header = T, sep = "\t")

dat$bin_ID <- gsub(":",".",dat$bin_ID)
dat$bin_ID <- gsub("-",".",dat$bin_ID)

genes_df$gene_IDs <- "NA"

for (b in 1:nrow(genes_df)){
  
  gbi <- which(dat$bin_ID == genes_df$bins[b])
  
  if(length(gbi) > 0){genes_df$gene_IDs[b] <- dat$gene_IDs[gbi]}
  
}

# export
write.table(genes_df, paste0(out,"/overlap_list_annotated_geneIDs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


