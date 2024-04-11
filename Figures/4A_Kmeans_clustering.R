#______Read in arguments________________________________________________________
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
IntFreq_df <- args[1]
outpath <- args[2]

#______Load required packages___________________________________________________
library(factoextra) 
library(ggplot2)

#______Number of clusters___________________________________________________
#df <- read.table("IntFreqSum_df_diploid.txt",header=TRUE)
df <- read.table(IntFreq_df), header= TRUE)
cell<-df$Cell
rownames(df)<-df[,1]
df<-df[,-1]
# Extract the two columns for clustering
df_cluster <- df[, c("Cis_IntFreq", "Trans_IntFreq")]

# Scale the data
df_cluster_scaled <- scale(df_cluster)

# Determine the optimal number of clusters using the elbow method
tp1<-{fviz_nbclust(df_cluster_scaled, kmeans, method = "wss") +
  ggtitle("Elbow Method for Optimal Number of Clusters")}

filename <- paste0(outpath,"/test_cis-trans number of clusters.pdf")
pdf(filename, width = 14, height = 8)
print(tp1)
dev.off()
print("done")
#______Kmeans_4 clusters___________________________________________________

# Perform K-Means clustering on the scaled data
set.seed(123)
best_k <- 4 #  optimal number of clusters
kmeans_fit <- kmeans(df_cluster_scaled, best_k, nstart = 25)
cluster<-kmeans_fit$cluster
cluster<-data.frame(cluster)
df$cluster <- as.factor(kmeans_fit$cluster)

# Plot the K-Means clustering results
tp2<-{fviz_cluster(kmeans_fit, df_cluster_scaled, geom = "point",
             fill = "cluster", ggtheme = theme_classic()) +
  ggtitle("K-Means Clustering Results")}

filename <- paste0(outpath,"/cis-trans kmean clustering (4 clusters).pdf")
pdf(filename, width = 14, height = 8)
print(tp2)
dev.off()
print("done") 
write.csv(df, file = sprintf("%s/Cis_Trans_kmean_clustering(4 clusters).csv", outpath), row.names = TRUE)
