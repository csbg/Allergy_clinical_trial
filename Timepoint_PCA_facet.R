# Data set first time point 

library(tidyverse)
library(readxl)
#install.packages("openxlsx")
library(openxlsx)
#install.packages("tableone")
library(tableone)
#install.packages("summarytools")
library(summarytools)
library(ggfortify)

dataset_tp <- read_xlsx("C:/Users/Tomic/OneDrive/Dokumente/R Data Science/Tables/dataset_timepoints.xlsx")

dataset_tp_new <-  dataset_tp

dataset_tp_new [, 6:18] <-  log(dataset_tp[,6:18]+0.001)


dataset_tp_log_scaled <- scale(dataset_tp_new[,3:18])

#result_tp <- dataset_tp[-1]

#result_tp_new <-  as.matrix(result_tp)

#result_tp_new[, 4:16] <-  log(result_tp_new[,4:16]+0.001)


#result_tp_log_scaled <- scale(result_tp_new)

###################other try################

df <- dataset_tp

#log transform and add 0.001
df[,6:18] <- log(df[,6:18] + 0.001)

#scale all 3:18
df[,3:18] <- scale(df[,3:18])

#save as excel table
write.xlsx(df, "df.xlsx")

#df is the data frame with the timepoint and the Group

# getting the patient ID to track the samples 

# Subset the data according to Timepoints

df_timepoint1 <- subset(df, Timepoint == '1')
df_timepoint2 <- subset(df, Timepoint == '2')
df_timepoint3 <- subset(df, Timepoint == '3')

##Make a new column in each subset a number starting with 1

df_timepoint1$ID <- 1:nrow(df_timepoint1)
df_timepoint2$ID <- 1:nrow(df_timepoint2)
df_timepoint3$ID <- 1:nrow(df_timepoint3)

# Merging the columns together to one 

df_timepoint1$TimeID <- paste(df_timepoint1$Group, "_ ", df_timepoint1$Timepoint, "_ ", df_timepoint1$ID, sep = "_ ")
df_timepoint2$TimeID <- paste(df_timepoint2$Group, "_ ", df_timepoint2$Timepoint, "_ ", df_timepoint2$ID, sep = "_ ")
df_timepoint3$TimeID <- paste(df_timepoint3$Group, "_ ", df_timepoint3$Timepoint, "_ ", df_timepoint3$ID, sep = "_ ")

#Combine the datasets 

df <- rbind(df_timepoint1, df_timepoint2, df_timepoint3)

View(df)

#save as excel table
write.xlsx(df, "df.xlsx")
#####################summary of the data###############


summarytools::dfSummary(df)

str(df)
df$Timepoint<-as.factor(df$Timepoint)

summary(df)

CreateTableOne( data=df)

tab1<-CreateTableOne(data = df,strata = "cluster_5K_a")

tab1


############PCA###############



library("ggfortify")

dataset_tp_log_scaled <-  as.matrix(dataset_tp_log_scaled)

pca_tp_log_scaled <- prcomp(dataset_tp_log_scaled)


autoplot(pca_tp_log_scaled, data= dataset_tp) + facet_grid( row= "Timepoint") +
  theme_bw()

autoplot(pca_tp_log_scaled, data= dataset_tp, colour= "Group") + facet_grid( row= "Timepoint") +
  theme_bw()


ggsave(filename = "PCA_facet_timepoints.pdf")


###### Timepoint 3##############

# Filter data frame to only include timepoint 3

df_timepoint3 <- df %>%
  filter(Timepoint == 3)

# Perform PCA on columns 3:18

pca_3 <- prcomp(df_timepoint3[,3:18], scale = TRUE)

autoplot(pca_3, data= df, colour= "Group") +theme_bw()


pca.var <- pca_3$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

loading_scores <- pca_3$rotation[,1]

biomarker_scores <- abs(loading_scores) ## get the magnitudes
biomarker_score_ranked <- sort(biomarker_scores, decreasing=TRUE)

top_10_biomarker_PCA3 <- names(biomarker_score_ranked[1:10])
top_5_biomarker_PCA3 <- names(biomarker_score_ranked[1:5])

Biomarker <- data.frame(
  Biomarker = top_10_biomarker_PCA3)

View(Biomarker)

barplot(top_5_biomarker_PCA5)

##############Heatmap Timepoint 3 ###################

set.seed(1234)
kmeans_clusters_timepoint3_5 <- kmeans(df_timepoint3[,3:18], centers = 5)

autoplot(kmeans_clusters_timepoint3_5, frame = TRUE, shape = 'Group', data = df_timepoint3, main= " Kmeans Clustering Timepoint 3")

Heatmap(kmeans_clusters_timepoint3_5$centers)

################################################

grid.arrange(tableGrob(Biomarker),autoplot(kmeans_model_5_a,df, sclae=TRUE, shape = 'Group', main= "Kmeans_5 Clustering splitted by timepoints") + facet_grid( row= "Timepoint") +
               theme_bw()),ncol=2)


##################other try#######################

#view the name of the top 10 biomarkers and the % of variance # fail 

top10 <- names(sort(pca_3$sdev^2, decreasing = TRUE))[1:10]
top10_var <- round(pca.var/sum(pca.var)*100, 1)[1:10]

data.frame(top10, top10_var)

################################################################################################




pheatmap(dataset_tp_log_scaled)

cor_data_las_tp <- cor(t(dataset_tp_log_scaled), method="spearman")

#annotation_col<- dataset_tp("TimePoint") 

pheatmap(cor_data_las_tp, scale = "none", col=colorRampPalette(c("blue", "white", "red"))(100)) 

#########################Clustermining###############

intervalo_k <- 1:6 # set the number of groups to run wss <- 0 # initialize an object that will store the wss
wss <- 0
for (k in intervalo_k) {
  km.out <- kmeans(x = df[,3:18], centers = k) # for each iteration, run the k means to the number "k" and save the result in km.out. Select only columns 1 and 2
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



###################### K-means ###########

#perform Kmeans clustering with all timepoints 

set.seed(1234)
kmeans_model_3_a <- kmeans(df[,3:18], centers = 3)
kmeans_model_5_a <- kmeans(df[,3:18], centers = 5)
kmeans_model_4_a <- kmeans(df[,3:18], centers = 4)

autoplot(kmeans_model_3_a,df, sclae=TRUE, frame= TRUE, shape = 'Group', main= "Kmeans_3 Clustering all timepoints")


autoplot(kmeans_model_3_a,df, sclae=TRUE, shape = 'Group', main= "Kmeans_3 Clustering all timepoints") + facet_grid( row= "Timepoint") +
  theme_bw()+ scale_shape_manual(values=c(5, 16,3))



autoplot(kmeans_model_4_a, df, scale=TRUE, frame=TRUE, shape = 'Group', main= "Kmeans_4 Clustering all timepoints")


autoplot(kmeans_model_4_a,df, sclae=TRUE, shape = 'Group', main= "Kmeans_3 Clustering all timepoints") + facet_grid( row= "Timepoint") +
  theme_bw()

autoplot(kmeans_model_5_a,df, scale=TRUE, frame=TRUE, shape = 'Group', main="Kmeans_5 Clustering all timepoints ")

autoplot(kmeans_model_5_a,df, sclae=TRUE, shape = 'Group', main= "Kmeans_5 Clustering all timepoints") + facet_grid( row= "Timepoint") +
  theme_bw()+scale_shape_manual(values=c(5, 16,3))


#make a new column with the cluster as factor
df$cluster_3K_a <- as.factor(kmeans_model_3_a$cluster)

df$cluster_5K_a <- as.factor(kmeans_model_5_a$cluster)
df$cluster_4K_a <- as.factor(kmeans_model_4_a$cluster)

write.xlsx(df, "df.xlsx")


##################PCA Kmeans3 and Kmeans5 ##############

autoplot(kmeans_model_5,df, frame= TRUE)


mat_log <- as.matrix(df[,3:18])

pca <- prcomp(mat_log, scale = TRUE)

df$PC1 <- pca$x[,1]
df$PC2 <- pca$x[,2]

write.xlsx(df, "df.xlsx")
autoplot(pca, data= df, colour= "cluster_3K",shape= "cluster_5K") + theme_bw()


################K-means timepoint 2 and 3#####################

#Subset data frame to only include timepoint 3 and 2

set.seed(1234)

df_subset <- df %>%
  filter(Timepoint %in% c('2', '3'))

# Perform kmeans clustering on columns 3:18
kmeans_clusters_3 <- kmeans(df_subset[,3:18], centers = 3)
kmeans_clusters_4 <- kmeans(df_subset[,3:18], centers = 4)
kmeans_clusters_5 <- kmeans(df_subset[,3:18], centers = 5)

# Autoplot the clusters
autoplot(kmeans_clusters_3, frame = TRUE, shape = 'Group', data = df_subset, main= "Kmeans_3 Clustering Timepoint 2 &3")
pdf(file="Clustering Kmeans 3 with Timepoint 2 and 3.pdf", width=8, height=8)

autoplot(kmeans_clusters_4, frame = TRUE, shape = 'Group', data = df_subset, main= "Kmeans_4 Clustering Timepoint 2 &3")
pdf(file="Clustering Kmeans 4 with Timepoint 2 and 3.pdf", width=8, height=8)

autoplot(kmeans_clusters_5, frame = TRUE, shape = 'Group', data = df_subset, main= "Kmeans_5 Clustering Timepoint 2 &3")
pdf(file="Clustering Kmeans 5 with Timepoint 2 and 3.pdf", width=8, height=8)

################Hierarchical Clustering############################


# Create a subset of the data containing only the biomarkers
df_biomarkers <- df[,3:18]

# Calculate the Euclidean distance between all samples
dist_matrix <- dist(df_biomarkers, method = "euclidean")

# Plot the distance matrix
library(ggplot2)

ggplot(as.matrix(dist_matrix)) +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_gradient(low = "white", high = "steelblue")

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "ward.D2")

# Plot the dendrogram
library(dendextend)
plot(as.dendrogram(hc))

# Color the dendrogram by group
group_colors <- c("Alutard" = "red", "BM41" = "green", "Placebo" = "blue")
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = 3, col = group_colors[df$group])
plot(dend)

#From the dendrogram, we can see that the Alutard samples form one cluster, 
#the BM41 samples form another cluster, and the Placebo samples form a third cluster. 
#herefore, we can conclude that the samples from each group are different from each other, but similar within the same group.

###########################some inspiration###################

# Generate clusters
clusters <- kmeans(df[,3:18], centers = 3)

# Identify responders and non-responders
responders <- df[clusters$cluster == 1,]
non_responders <- df[clusters$cluster == 2,]

# Calculate overall response and effect of the treatment
alutard_response <- mean(responders$Alutard)
bm41_response <- mean(responders$BM41)
placebo_response <- mean(responders$Placebo)

# Calculate contribution of each measurement to the variance of the treatment
alutard_variance <- var(responders$Alutard)
bm41_variance <- var(responders$BM41)
placebo_variance <- var(responders$Placebo)

# Calculate correlations between measurements and clusters
alutard_corr <- cor(responders$Alutard, clusters$cluster)
bm41_corr <- cor(responders$BM41, clusters$cluster)
placebo_corr <- cor(responders$Placebo, clusters$cluster)


#?aov  Fit an Analysis of Variance 


# Cluster the data
df_cluster <- kmeans(df[,3:18], 3)

# Calculate contribution of each measurement to the variance of the treatment

contrib <- apply(df[,3:18], 2, function(x) cor(x, df_cluster$cluster))


############################ T-test significant Biomarkers############################

# Create a data frame containing only the three groups of interest 
df2 <- df[df$Group %in% c("Alutard", "BM41", "Placebo"), ] 

# Create a list of the biomarker columns 
biomarker_cols <- colnames(df2)[3:18] 

# Run a series of t-tests to compare the effects of the treatments 
res <- lapply(biomarker_cols, function(x) 
  t.test(df2[df2$Group == "Alutard", x], df2[df2$group == "BM41", x], 
         df2[df2$Group == "Placebo", x], 
         paired = FALSE, var.equal = TRUE)) 

# Get the p-values for each of the tests 
p_values <- sapply(res, function(x) x$p.value) 











 



