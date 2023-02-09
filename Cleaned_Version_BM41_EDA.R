library(tidyverse)
library(mlr)
#library(dplyr)
library(ggplot2)
library(reshape2)
library(GGally)
library(FactoMineR)
library(data.table)
library(pheatmap)

### BM41 dataset with times points ####

setwd("~/R Data Science/BM41/")

dataset <- read_tsv("k-meansBM41StandartVar.txt")
str(dataset)
dim(dataset)
colnames(dataset)

result <- dataset[-1]
dim(result)
result <- as.data.frame(result[-1])# putting the first column as the row name
dim(result)
colnames(result)
row.names(result) <- dataset$TimeID
result
result <- result[complete.cases(result), ]

str(result)
colnames(result)
dim(result)
str(result)

  

################# Data distribution  ##################

nfDT <- dataset %>%
  pivot_longer(!Group & !TimeID, names_to = "measurement", values_to = "value" )


ggplot(nfDT, aes(x=value)) + geom_density() + 
  facet_wrap(~measurement, scales = "free")





################# Normalization ################

dataset_new <-  dataset
dataset_new [, 6:18] <-  log(dataset[,6:18]+0.001)


dataset_log_scaled <- scale(dataset_new[,3:18])



result_new <-  as.matrix(result)

result_new[, 4:16] <-  log(result_new[,4:16]+0.001)



result_log_scaled <- scale(result_new)

boxplot(result_log_scaled, las=2)



##################### Data distribution normalized #########################

nfDT_log_scaled <-
  dataset_log_scaled %>%
  as_tibble() %>% 
  pivot_longer(everything(), names_to = "measurement", values_to = "value" )



ggplot(nfDT_log_scaled, aes(x=value)) + geom_density() + 
  facet_wrap(~measurement, scales = "free")

#Substracting the mean and dividing by the standard deviation 




##############################PCA ############################

library("ggfortify")

pca_log_scaled <- prcomp(result_log_scaled)

autoplot(pca_log_scaled, data= dataset, colour= "Group")

#screeplot(pca_log_scaled)


pca.var <- pca_log_scaled$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

loading_scores <- pca_log_scaled$rotation[,1]

biomarker_scores <- abs(loading_scores) ## get the magnitudes
biomarker_score_ranked <- sort(biomarker_scores, decreasing=TRUE)
top_10_biomarker <- names(biomarker_score_ranked[1:10])





################################# Correlation- Plot measurements ######


corr_matrix <- round(cor(result_log_scaled, method = "spearman"),digits = 2)


library("Hmisc")
library(lattice)
library(survival)
library(Formula)
library(ggcorrplot)
library(corrplot)

res <- rcorr(as.matrix(result_log_scaled), type = "spearman") 
round(res$P,3 )



colnames(corr_matrix)[1] <- "Inhibition RBL "
colnames(corr_matrix)[2] <- "Inhibition FAB "
colnames(corr_matrix)[3] <- "Competition FAB "
colnames(corr_matrix)[4] <- "Standard FAB "


rownames(corr_matrix)[1] <- "Inhibition RBL "
rownames(corr_matrix)[2] <- "Inhibition FAB "
rownames(corr_matrix)[3] <- "Competition FAB "
rownames(corr_matrix)[4] <- "Standard FAB "


p_mat <- cor_pmat(result_log_scaled, method= "spearman")

corrplot.mixed(corr_matrix, 
               lower = "circle",
               upper = "number", 
               tl.pos = "lt",
               tl.cex = 0.8, 
               tl.col = "black", 
               diag = "n",
               p.mat = p_mat,
               insig = "blank",
        
               title = "Significant Correleations between measurements")


###################### Heat maps #######################



cor_data_las <- cor(t(result_log_scaled))

cor_data_las <- cor(t(result_log_scaled), method="spearman")

pheatmap(cor_data_las, scale = "none", col=colorRampPalette(c("blue", "white", "red"))(100))


# split with colsplit - see cluster with the Dentrogram 

library(ComplexHeatmap)

Heatmap(cor_data_las, col=colorRampPalette(c("blue", "white", "red"))(100), row_split = 3, column_split = 3 )


dev.off()
pdf("CorrelationHeatmap_las.pdf", w=20, h=20)

dev.off()

# i really don´t know why this gives this other plot but it does alothough it is the same code and there are no value above and under 1/-1 


cor_data <- cor(t(result_log_scaled))

cor_data <- cor(result_log_scaled, method="spearman") 

pheatmap(result_log_scaled, scale = "none", col=colorRampPalette(c("blue", "white", "red"))(100))

pdf("CorrelationHeatmap_las_measurements.pdf", w=20, h=20)

Heatmap(result_log_scaled, col=colorRampPalette(c("blue", "white", "red"))(100),row_split = 3, column_split = 3 )

###################### K-means Algorithm #############################


# Run the K means algorithm, for 3 groups in values


set.seed(1234)

result_cluster_las <- dataset_log_scaled


clusters_kmean_k3_las<- kmeans(x = result_log_scaled, centers = 3)

result_cluster_las$Groups_k3 <- as.factor(clusters_kmean_k3_las$cluster)


#Cluster Plot 

autoplot(clusters_kmean_k3_las,result_log_scaled, frame= TRUE)


#Cluster centers 

clusters_kmean_k3_las$centers


# Plot

pheatmap(clusters_kmean_k3_las$centers)



######################### UMAP #################################

library(umap)

las_UMAP = umap(result_log_scaled)

plot(las_UMAP$layout)


umapdataframe_las <- las_UMAP$layout

umapdataframe_las <-  as.data.frame(umapdataframe_las)

result_cluster_las$UMAP1 <- umapdataframe_las$V1

result_cluster_las$UMAP2 <-  umapdataframe_las$V2




ggplot(umapdataframe_las, aes(x= V1, y=V2,shape= result_cluster_las$Groups_k3, color=result_cluster_las$Group)) +
  geom_point() +
  scale_color_manual(result_cluster_las$Groups_k3, name= "Kmeans Cluster", values = c( ""))

  ggsave("UMAP_Plot_las_all_timepoints.pdf")


##################################################Separate############
  
result$rows  <- row.names(result)
  
result<- result %>% separate(rows, c("Groups", "TimePoint"), sep="_")

result$TimePoint <- last(result$TimePoint)

view(result)



######################### New column###################



