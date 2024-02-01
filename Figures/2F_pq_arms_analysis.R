#______Read in arguments________________________________________________________
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

data <- read.table(args[1]
outpath <- args[2]

#______Load required packages_____________________________________________

library(dplyr)
library(tidyr)
library(parallel)
library(data.table)
library(ggplot2)

#______Get the data imported and ready_________________________________________

row_sums <- rowSums(data[, -1], na.rm = TRUE) 
data$Total <- row_sums
print(data)
df_1MB <- data[, c(1, ncol(data))] 

#______Data prepration___________________________________________________

colnames(df_1MB)[1] <- "ID"

df_1MB <- df_1MB[-1, ]
head(df_1MB)
print("reading table 1MB .....")



df_1MB <- data.frame(
  df_1MB$ID,
  do.call(rbind, strsplit(as.character(df_1MB$ID), "[_/ ]")),
  df_1MB$Total,
  stringsAsFactors = FALSE
)
df_1MB<- df_1MB[,-6]

colnames(df_1MB) <- c("ID", "anchor", "chr1", "start", "end", "target", "chr2", "start2", "end2", "Total")

chrInf <- data.frame( chrom = c("chr1","chr2",
                                "chr3","chr4",
                                "chr5","chr6",
                                "chr7","chrX",
                                "chr8","chr9",
                                "chr11","chr10",
                                "chr12","chr13",
                                "chr14","chr15",
                                "chr16","chr17",
                                "chr18","chr20",
                                "chr19","chrY",
                                "chr22","chr21"),
                      centromere = c(123252373.5,93787431.5,
                                     90856062,50074452.5,
                                     48585285.5,60557102.5,
                                     60058972.5,61016889,
                                     45249872,43893383.5,
                                     53454152,39800499.5,
                                     35764400,17692000.5,
                                     17117352,19037747.5,
                                     36878628.5,25067566.5,
                                     18464134,28099979.5,
                                     26161912,10470308,
                                     15520235.5,11917946),
                      chrClass = c("Metacentric","Metacentric",
                                   "Metacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Acrocentric",
                                   "Acrocentric","Acrocentric",
                                   "Metacentric","Submetacentric",
                                   "Submetacentric","Metacentric",
                                   "Metacentric","Acrocentric",
                                   "Acrocentric","Acrocentric"),
                      size = c(248956422,242193529,
                               198295559,190214555,
                               181538259,170805979,
                               159345973,156040895,145138636,
                               138394717,135086622,
                               133797422,133275309,
                               114364328,107043718,
                               101991189,90338345,
                               83257441,80373285,
                               64444167,58617616,57227415,
                               50818468,46709983))

#_________Locate centromere for each chromosome and add interaction type column for anchor and target___________

df_1MB$chr1_centromere <- sapply(df_1MB$chr1, function(x) chrInf[chrInf$chrom == x, "centromere"])
df_1MB$chr2_centromere <- sapply(df_1MB$chr2, function(x) chrInf[chrInf$chrom == x, "centromere"])

for (i in 1:nrow(df_1MB)) {
  if (df_1MB[i,"start"] <= df_1MB[i, "chr1_centromere"]) {
    df_1MB[i, "interaction_type_anchor"] <- "p_arm"
  } else {
    df_1MB[i, "interaction_type_anchor"] <- "q_arm"
  }
  
  if (df_1MB[i,"start2"] <= df_1MB[i, "chr2_centromere"]) {
    df_1MB[i, "interaction_type_target"] <- "p_arm"
  } else {
    df_1MB[i, "interaction_type_target"] <- "q_arm"
  }
}


df_1MB$chr2_centromere <- df_1MB$chr2_centromere[[1]]


write.table(df_1MB, file = sprintf("%s/df_1MB_genome_wide.txt", outpath), sep = "\t", row.names = FALSE)

#________Add the final interaction_type column to each file____________________________________

  check_interaction <- function(interaction_type_anchor, interaction_type_target) {
    if (interaction_type_anchor == "p_arm" & interaction_type_target == "p_arm") {
      return("pp")
    } else if (interaction_type_anchor == "p_arm" & interaction_type_target == "q_arm") {
      return("pq")
    } else if (interaction_type_anchor == "q_arm" & interaction_type_target == "p_arm") {
      return("pq")
    } else if (interaction_type_anchor == "q_arm" & interaction_type_target == "q_arm") {
      return("qq")
    } else {
      return("")
    }
  }
  
df_1MB$interaction_type <- apply(df_1MB[, c("interaction_type_anchor", "interaction_type_target")], 1, function(x)
  check_interaction(x[1], x[2]))


write.table(df_1MB, file = sprintf("%s/df_1MB_genome_wide_final.txt", outpath), sep = "\t", row.names = FALSE)


#_________Apply binomial_test___________________________________________
data <- read.table("df_1MB_genome_wide_final.txt", header = TRUE)

pp_rows <- data$interaction_type == "pp"
total_pp <- sum(data$Total[pp_rows])
total_pp

pq_rows <- data$interaction_type == "pq"
total_pq <- sum(data$Total[pq_rows])
total_pq 
 

qq_rows <- data$interaction_type == "qq"
total_qq <- sum(data$Total[qq_rows])
total_qq 
total_int_freq<-c(total_pp,total_pq,total_qq)


prob_total_qq<-total_qq/(total_pp+total_pq+total_qq)
prob_total_pq<-total_pq/(total_pp+total_pq+total_qq)
prob_total_pp<-total_pp/(total_pp+total_pq+total_qq)
prob_total_int_freq<- c(prob_total_pp,prob_total_pq,prob_total_qq)
chrinf <- read.table("chrInf.txt", header = TRUE)
chrinf
chrinf$pq_new<- chrinf$pq+chrinf$qp
chrinf
#chrinf<-chrinf[,-14]
chrinf
q_arm_sum<- sum(chrinf$q_arm)
p_arm_sum<- sum(chrinf$p_arm)
size_sum<- sum(chrinf$size)

probability_p_arm<-(p_arm_sum/size_sum)
probability_q_arm<-(q_arm_sum/size_sum)
probability_pp<-(probability_p_arm)*(probability_p_arm)
probability_pq<-2*(probability_p_arm)*(probability_q_arm)
probability_qq<-(probability_q_arm)*(probability_q_arm)  

interaction_types <- c("pp", "pq", "qq")
interaction_prob<-c(probability_pp,probability_pq,probability_qq)

df <- data.frame(interaction_types, interaction_prob)

results_df <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(results_df) <- c("pp_pvalue", "pp_direction", "pq_pvalue", "pq_direction", "qq_pvalue", "qq_direction")

# perform binomial test for all three interaction types
for (i in 1:length(interaction_types)) {
  
  samples <- data[data$interaction_type == interaction_types[i], ]
  success <- round(sum(samples$Total))
  trials <- round(sum(data$Total))
  
  prob <- df$interaction_prob[df$interaction_types == interaction_types[i]]
  
  test_result <- binom.test(success, trials, prob)
  
  # save p-value and direction for each interaction type in results_df
  if (interaction_types[i] == "pp") {
    results_df$pp_pvalue <- test_result$p.value
    if (test_result$p.value <= 0.05) {
      if (success / trials > prob) {
        results_df$pp_direction <- "over_represented"
      } else {
        results_df$pp_direction <- "under_represented"
      }
    } else {
      results_df$pp_direction <- "not_significant"
    }
  }
  if (interaction_types[i] == "pq") {
    results_df$pq_pvalue <- test_result$p.value
    if (test_result$p.value <= 0.05) {
      if (success / trials > prob) {
        results_df$pq_direction <- "over_represented"
      } else {
        results_df$pq_direction <- "under_represented"
      }
    } else {
      results_df$pq_direction <- "not_significant"
    }
  }
  if (interaction_types[i] == "qq") {
    results_df$qq_pvalue <- test_result$p.value
    if (test_result$p.value <= 0.05) {
      if (success / trials > prob) {
        results_df$qq_direction <- "over_represented"
      } else {
        results_df$qq_direction <- "under_represented"
      }
    } else {
      results_df$qq_direction <- "not_significant"
    }
  }
}

data.frame(interaction_types,interaction_prob)
interaction_pvalues <- c(results_df$pp_pvalue, results_df$pq_pvalue, results_df$qq_pvalue)
interaction_directions <- c(results_df$pp_direction, results_df$pq_direction, results_df$qq_direction)
interaction_results_genome_wide <- data.frame(interaction_types=interaction_types,prob_total_int_freq=prob_total_int_freq,interaction_prob = interaction_prob,
                                              interaction_pvalues = interaction_pvalues,
                                              interaction_directions = interaction_directions)


write.table(
  interaction_results_genome_wide,
  file = sprintf("%s/interaction_results_genome_wide.txt", outpath),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
#______Bar_plot________________________________
Interaction_type<-c("pp","pp","pq","pq","qq","qq")
Probability<-c(0.2699555, 0.1118366,0.5028044,0.4451664,0.2273771,0.4429970)
Prob_type<- c("observed","expected","observed","expected","observed", "expected")

data<-data.frame(Interaction_type,Probability,Prob_type)
ggplot(data,                                      # Grouped barplot using ggplot2
       aes(x =Interaction_type ,
           y = Probability,
           fill = Prob_type)) +
  geom_bar(stat = "identity",
           position = "dodge")+ theme(panel.background = element_blank())

p_values <- c("p < 1e-10", "p < 1e-10", "p < 1e-10")

ggplot(data, aes(x = Interaction_type, y = Probability, fill = Prob_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  annotate("text", x = 0.82, y = 0.32, label = p_values[1], size = 3, fontface = "bold") +
  annotate("text", x = 2, y = 0.55, label = p_values[2], size = 3, fontface = "bold") +
  annotate("text", x = 3, y = 0.48, label = p_values[3], size = 3, fontface = "bold") +
  geom_segment(aes(x = .5, y = 0.29, xend = 1.2, yend = 0.29), size = 0.8) +
  geom_segment(aes(x = 1.6, y = 0.53, xend = 2.35, yend = 0.53), size = 0.8) +
  geom_segment(aes(x = 2.7, y = 0.46, xend = 3.4, yend = 0.46), size = 0.8) +
  scale_fill_manual(values = c("darkorange1", "grey90")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(sprintf("%s/AllvsAll.pdf", outpath), width = 8, height = 5, dpi = 300)
