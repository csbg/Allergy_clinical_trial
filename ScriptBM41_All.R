library(tidyverse)
library(mlr)
#library(dplyr)
library(ggplot2)
library(reshape2)
library(GGally)
library(FactoMineR)
library(data.table)
library(pheatmap)

### BM41 dataset with times point ####

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

df_result_log_scaled <- as_data_frame(result_log_scaled)

write_csv(result_log_scaled, file='result_log_scaled.csv')





set.seed(1234) # because it is an iterative algorithm, different configurations can generate different results

# Run the K means algorithm, for 1, 2, 3, 4, 5 and 6 groups
intervalo_k <- 1:6 # set the number of groups to run wss <- 0 # initialize an object that will store the wss
wss <- 0
for (k in intervalo_k) {
  km.out <- kmeans(x = result, centers = k) # for each iteration, run the k means to the number "k" and save the result in km.out. Select only columns 1 and 2
  wss[k] <- km.out$tot.withinss # save the wss result of that iteration in the indexed 'wss' object for the respective iteration
}

# View the results in a graph
ggplot(data=data.frame(intervalo_k, wss), aes(x=intervalo_k, y=wss, group=1)) +
  geom_line() +
  geom_point() +
  xlab("Number of clusters") +
  ylab("WSS") +
  scale_x_continuous(breaks=c(1:6))

cluster.mining<- data.frame(intervalo_k, wss)
write.csv(cluster.mining, file = "./cluster_mining.csv")

# Run the K means algorithm for 2 groups 

set.seed(1234)

clusters_kmean_k2_todos <- kmeans(x = result, centers = 2) 

ResultsCluster <- dataset # create new object not to overwrite original data
ResultsCluster$Groups_k2 <- as.factor(clusters_kmean_k2_todos$cluster)

# Plot

ggplot(ResultsCluster, aes(x=Groups_k2, y=Group,  shape=Groups_k2, color=Group)) +
  geom_point()

ggplot(ResultsCluster, aes(x=TimeID, y=Groups_k2,  shape=Groups_k2, color=Group)) +
  geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Run the algorithm K means, for 3 groups in all parameters

set.seed(1234)

clusters_kmean_k3_todos <- kmeans(x = result, centers = 3) 


ResultsCluster$Groups_k3 <- as.factor(clusters_kmean_k3_todos$cluster)

# Plot

ggplot(ResultsCluster, aes(x=TimeID, y=Groups_k3,  shape=Groups_k3, color=Group)) +
  geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Run the K means algorithm, for 4 groups in all parameters

set.seed(1234)

clusters_kmean_k4_todos <- kmeans(x = result, centers = 4) 


ResultsCluster$Groups_k4 <- as.factor(clusters_kmean_k4_todos$cluster)

# Plot

ggplot(ResultsCluster, aes(x=TimeID, y=Groups_k4, shape=Groups_k4, color=Group)) +
  geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Run the K means algorithm, for 5 groups in all genes

set.seed(1234)

clusters_kmean_k5_todos <- kmeans(x = result, centers = 5) 


ResultsCluster$Groups_k5 <- as.factor(clusters_kmean_k5_todos$cluster)

# Plot

ggplot(ResultsCluster, aes(x=TimeID, y=Groups_k5, shape=Groups_k5, color=Group)) +
  geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Re-write the table

write.csv2(ResultsCluster, file = "k-meansBM41STimePoint.csv")

write.csv2(clusters_kmean_k5_todos$centers, file = "k-meansBM41TimePoint.csv")


## Variables vs clusters

pheatmap(clusters_kmean_k5_todos$centers)

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(clusters_kmean_k5_todos$centers, 1, cal_z_score))
pheatmap(data_subset_norm)

#### UMAP #####

library(umap)


bm41.umap = umap(result)

head (bm41.umap$layout)
plot(bm41.umap$layout)


umapdataframe <- bm41.umap$layout
umapdataframe <-as.data.frame(umapdataframe)

ResultsCluster$UMAP1 <- umapdataframe$V1

ResultsCluster$UMAP2 <- umapdataframe$V2


ggplot(umapdataframe, aes(x=V1, y=V2, shape=ResultsCluster$Group, color=ResultsCluster$Groups_k5)) +
  geom_point()

ggplot(umapdataframe, aes(x=V1, y=V2, shape=ResultsCluster$Groups_k5, color=ResultsCluster$Group)) +
  geom_point()

############ PCA ################


install.packages("ggfortify")
library("ggfortify")

pca <- prcomp(result)

head(pca)


pca_res <- prcomp(result, scale. = TRUE)

autoplot(pca_res, data= dataset, colour= 'Group')
 
autoplot(pca_res, data= dataset, colour= 'Group', 
         loadings= TRUE)

autoplot(pca_res, data= dataset, colour= 'Group', 
         loadings= TRUE, 
         loadings.label= TRUE)

#--- PCA Log_adjusted_scaled --------------


pca_log_scaled <- prcomp(result_log_scaled)

autoplot(pca_log_scaled, data= dataset, colour= "Group")

pca.var <- pca_log_scaled$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

loading_scores <- pca_log_scaled$rotation[,1]

biomarker_scores <- abs(loading_scores) ## get the magnitudes
biomarker_score_ranked <- sort(biomarker_scores, decreasing=TRUE)

top_10_biomarker <- names(biomarker_score_ranked[1:10])

view(top_10_biomarker)

############# Heat map #############################

cor_data <- cor(t(result))

cor_data <- cor(t(scale(result)), method="spearman") 

pheatmap(scale(log10(result + 0.0001)))

dev.off()
dev.off()
pdf("CorrelationHeatmap.pdf", w=20, h=20)
pheatmap(cor_data, scale = "none", col=colorRampPalette(c("blue", "white", "red"))(100))
dev.off()


boxplot(result_log_scaled, las=2)



cor_data_las <- cor(t(result_log_scaled))

cor_data_las <- cor(t(result_log_scaled), method="spearman")

pheatmap(cor_data_las, scale = "none", col=colorRampPalette(c("blue", "white", "red"))(100))

dev.off()
pdf("CorrelationHeatmap_las.pdf", w=20, h=20)

dev.off()

library(ComplexHeatmap)
Heatmap(kmeans_model_5_a$centers)

################ Density Plot #######################



colnames(result) <- make.names(colnames(result))
plot(density(log2(result$rBet.v.1sIgG4..µg.L.)))
str(result)

plot(density(log2(result$Serology.Birch)))




plot(density(log2(result$ rBet.v.1.sIgE..kU.L.)))

plot(density(log2(result$rBet.v.2.sIgE..kU.L.)))

plot(density(log2(result$BM41.sIgE..kU.L.)))



plot(density(log2(result$Birch.sIgG4..µg.L.)))

plot(density(log2(result$ rBet.v.2.sIgG4..µg.L.)))

plot(density(log2(result$ BM41.sIgG4..µg.L.)))

     
plot(density(log2(result$ Birch.sIgG..mg.L.)))     

plot(density(log2(result$ rBet.v.1.sIgG..mg.L.)))

plot(density(log2(result$ rBet.v.2.sIgG..mg.L.)))

plot(density(log2(result$ BM41.sIgG..mg.L.)))
     
plot(density(log2(result$ Standard.FAB.aIgE..gated.. )))





result %>%
  ggplot(aes_string(x="rBet.v.1sIgG4..µg.L.")) +
  geom_density()



tidy_dataset <- dataset %>%
  
  pivot_longer(c("Inhibition RBL normalized in Inhibition", 
                 "Inhibition FAB ITC-aIgE+ (%) normalized compared all in% inhibition",
                  "Competition FAB FITC-aIgE+ (%) normalized compare all in% inhibition",
                 "Standard FAB aIgE+ gated %",
                 "Serology Birch",
                 "rBet v 1 sIgE (kU/L)",
                 "rBet v 2 sIgE (kU/L)",
                 "BM41 sIgE (kU/L)",
                 "Birch sIgG4 (µg/L)",
                 "rBet v 1sIgG4 (µg/L)", 
                 "rBet v 2 sIgG4 (µg/L)",
                 "BM41 sIgG4 (µg/L)",
                 "Birch sIgG (mg/L)",
                 "rBet v 1 sIgG (mg/L)",
                 "rBet v 2 sIgG (mg/L)",
                 "BM41 sIgG (mg/L)"),
               
                names_to = "measurement", values_to = "value" )


#shorter expression

nfDT <- dataset %>%
  pivot_longer(!Group & !TimeID, names_to = "measurement", values_to = "value" )



ggplot(nfDT, aes(x=value)) + geom_density() + 
  facet_wrap(~measurement, scales = "free")


ggplot(nfDT, aes(x=log10(value + 0.0001))) + geom_density() + 
  facet_wrap(~measurement, scales = "free")
ggsave("nfDF_density_log10.pdf")

nfDT[which(is.na(log10(nfDT$value))),]$value


ggplot(nfDT, aes(x=time)) + geom_density() + facet_wrap(~measurement, scales = "free")

table(nfDT$time, nfDT$Group)

nfDT$time <- as.numeric(gsub("^.+?_t", "", nfDT$TimeID))



ggplot(nfDT, aes(x=time, y=value, color=Group)) + 
  
  geom_smooth() + facet_wrap(~measurement, scales = "free")


ggplot(nfDT, aes(x=time, y=log10(value), color=Group)) + 
  
  geom_point() + facet_wrap(~measurement, scales = "free")


str(dataset)

ggplot(data= tidy_dataset) + 
  geom_point( mapping= aes_string(x= "TimeID", y= "rBet v 1 sIgE (kU/L)"))+
  facet_wrap(~ Group, nrow = 3)





#################log Transformation ###################

dataset_new <-  dataset
  
dataset_new [, 6:18] <-  log(dataset[,6:18]+0.001)

  
dim(dataset_new)
str(dataset_new)


nfDT_log <- dataset_new %>%
  pivot_longer(!Group & !TimeID, names_to = "measurement", values_to = "value" )


ggplot(nfDT_log, aes(x=value)) + geom_density() + 
  facet_wrap(~measurement, scales = "free")


view(dataset_new)

nfDT_log %>% mutate(scaled = scale(value))

################### Density Plot for log adjusted and scaled Data #######################
  


nfDT_log_scaled <-
  dataset_log_scaled %>%
  as_tibble() %>% 
  pivot_longer(everything(), names_to = "measurement", values_to = "value" )


  
  ggplot(nfDT_log_scaled, aes(x=value)) + geom_density() + 
  facet_wrap(~measurement, scales = "free")


#make a result table for the log transformed data without character 

  
  
result_new <-  as.matrix(result)
str(result_new)

result_new[, 4:16] <-  log(result_new[,4:16]+0.001)
str(result_new)
quantile(result_new)


view(result_new)

result_log_scaled <- scale(result_new)

pca_log <- prcomp(x= result_new, scale= TRUE)

stopifnot(all(dataset$TimeID == row.names(result_new)))

autoplot(pca_log, data= dataset, colour= 'Group')

autoplot(pca_log, data= dataset, colour= 'Group', 
         loadings= TRUE)

autoplot(pca_log, data= dataset, colour= 'Group', 
         loadings= TRUE, 
         loadings.label= TRUE)





ggplot(nfDT_log, aes(x=time, y=value, color=Group)) + 
  
  geom_smooth() + facet_wrap(~measurement, scales = "free")



boxplot(result_new, las=2)
boxplot(scale(result_new), las=2)


cor_data <- cor(t(result_new))


cor_data <- cor(t(scale(result_new)), method="spearman") 

pheatmap(scale(result_new))




# K means - log_adapted data

# Run the K means algorithm, for 5 groups in values


set.seed(1234)

result_cluster_log <- dataset_new

view(result_cluster_log)

clusters_kmean_k5_log <- kmeans(x = result_new, centers = 5) 

#Cluster Plot 

autoplot(clusters_kmean_k5_log,result_new, frame= TRUE)

#Cluster centers 

clusters_kmean_k5_log$centers

# new table with Kmean group 

result_cluster_log$Groups_k5_log <- as.factor(clusters_kmean_k5_log$cluster)


# Plot

ggplot(result_cluster_log, aes(x=TimeID, y=Groups_k5_log, shape=Groups_k5_log, color=Group))+
  geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

pheatmap(clusters_kmean_k5_log$centers)
     
# UMAP 

library(umap)

log_UMAP = umap(result_new)

plot(log_UMAP$layout)

umapdataframe_log <- log_UMAP$layout

umapdataframe_log <-  as.data.frame(umapdataframe_log)

result_cluster_log$UMAP1 <- umapdataframe_log$V1

result_cluster_log$UMAP2 <-  umapdataframe_log$V2

# Plot 

ggplot(umapdataframe_log, aes(x=V1, y=V2,shape=result_cluster_log$Group,color= result_cluster_log$Groups_k5_log)) +
        
        geom_point()


ggplot(umapdataframe_log, aes(x=V1, y=V2,shape=result_cluster_log$Groups_k5_log,color= result_cluster_log$Group)) +
  
  geom_point()


################################# Correlation- Plots######

# correlation for all variables 

corr_matrix <- round(cor(result_log_scaled, method = "spearman"),digits = 2)

install.packages("Hmisc")
library("Hmisc")
library(lattice)
library(survival)
library(Formula)

res <- rcorr(as.matrix(result_log_scaled), type = "spearman") 
round(res$P,3 )

install.packages("corrplot")

library(ggcorrplot)
library(corrplot)



colnames(corr_matrix)[1] <- "Inhibition RBL "
colnames(corr_matrix)[2] <- "Inhibition FAB "
colnames(corr_matrix)[3] <- "Competition FAB "
colnames(corr_matrix)[4] <- "Standard FAB "


rownames(corr_matrix)[1] <- "Inhibition RBL "
rownames(corr_matrix)[2] <- "Inhibition FAB "
rownames(corr_matrix)[3] <- "Competition FAB "
rownames(corr_matrix)[4] <- "Standard FAB "


#rownames(corr_matrix) <- rownames(corr_matrix) %>% str_sub(1, 10)
corrplot(corr_matrix, 
         method = "circle", 
         type = "full", 
         order = "hclust",
         title = "Correlogram measurements",
         cl.pos = "r")
         
corrplot.mixed(corr_matrix, 
               lower = "circle",
                upper = "number", 
               tl.pos = "lt",
               tl.cex = 0.8, 
               tl.col = "black", 
               diag = "n",
               title = "Correlogram measurements")


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


