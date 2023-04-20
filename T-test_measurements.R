
library(dplyr)
library(openxlsx)
library(tidyverse)
library(broom)
library(tidyr)
library(ggplot2)
library(readxl)
library(tableone)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)
library(grDevices)



#Assign response levels based on rBetv sIgG values-----------


# Load data
df_a <- read.xlsx("data_generated/df_scaled.xlsx")

# Subset data for Timepoint 3 and Alutard treatment
alutard_df <- df_a %>% 
  filter(Timepoint == 3, Group == "Alutard") 

# Subset data for Timepoint 3 and BM41 treatment
bm41_df <- df_a %>% 
  filter(Timepoint == 3, Group == "BM41") 

# Assign response levels based on rBetv sIgG values
assign_response_levels <- function(df, top_values, bottom_values) {
  df %>%
    mutate(response_degree= case_when(
      rBetv1_sIgG %in% top_values$rBetv1_sIgG ~ "high",
      rBetv1_sIgG %in% bottom_values$rBetv1_sIgG ~ "low",
      TRUE ~ "intermediate"
    ))
}

# Get the top 5 and bottom 5 values of rBetv sIgG for Alutard treatment
alutard_top5 <- alutard_df %>% 
  arrange(desc(rBetv1_sIgG)) %>% 
  slice(1:5)

alutard_bottom5 <- alutard_df %>% 
  arrange(rBetv1_sIgG) %>% 
  slice(1:5)

# Assign response levels based on rBetv sIgG values for Alutard treatment
alutard_df <- assign_response_levels(alutard_df, alutard_top5, alutard_bottom5)

# Get the top 5 and bottom 5 values of rBetv sIgG for BM41 treatment
bm41_top5 <- bm41_df %>% 
  arrange(desc(rBetv1_sIgG)) %>% 
  slice(1:5)

bm41_bottom5 <- bm41_df %>% 
  arrange(rBetv1_sIgG) %>% 
  slice(1:5)

# Assign response levels based on rBetv sIgG values for BM41 treatment
bm41_df <- assign_response_levels(bm41_df, bm41_top5, bm41_bottom5)

# Combine the two dataframes
combined_df <- bind_rows(alutard_df, bm41_df)

view(combined_df)

# Generate table one
CreateTableOne(data = combined_df, strata = "response_degree")

# Write to excel file
write.xlsx(combined_df, "data_generated/rBet_v_1_sIgG_Response_table.xlsx")



# High vs. Low Responder t-test --------------------------

#Loop through each biomarker and perform a t-test

##################Alutard############################

alutard_high_response <- alutard_top5[3:18]

alutard_low_response <- alutard_bottom5[3:18]

# Initialize vectors to store results
p_values_alutard <- numeric(ncol(alutard_high_response))
t_statistic_alutard <- numeric(ncol(alutard_high_response))
effect_size_alutard <- numeric(ncol(alutard_high_response))



# Loop through each biomarker and perform a t-test
measurements <- colnames(alutard_df)[3:18]


for (i in 1:ncol(alutard_high_response)) {
  biomarker1 <- alutard_high_response[, i]
  biomarker2 <- alutard_low_response[, i]
  ttest_result <- t.test(biomarker1, biomarker2)
  pvalue <- ttest_result$p.value
  t_stat <- ttest_result$statistic
  cohens_d <- abs(mean(biomarker1) - mean(biomarker2)) / sd(c(biomarker1, biomarker2))
  biomarker_name <- colnames(alutard_high_response)[i]
  print(paste("P-value for", biomarker_name, ":", pvalue))
  print(paste("t-statistic for", biomarker_name, ":", t_stat))
  print(paste("Cohen's d for", biomarker_name, ":", cohens_d))
  p_values_alutard[i] <- pvalue
  t_statistic_alutard[i] <- t_stat
  effect_size_alutard[i] <- cohens_d
  
}



# new dataframe with the Biomarker and p-values 

biomarker_alutard <- data.frame(Parameter = names(alutard_df)[3:18], p_value = p_values_alutard, t_statistic = t_statistic_alutard, effect_size = effect_size_alutard)

biomarker_alutard$Treatment <- "Alutard"



#######################BM41############################## 

#Subset data 

BM41_high_response <- bm41_top5[3:18]

BM41_low_response <- bm41_bottom5[3:18]


# Initialize vectors to store results
p_values_BM41 <- numeric(ncol(BM41_high_response))
t_statistic_BM41 <- numeric(ncol(BM41_high_response))
effect_size_BM41 <- numeric(ncol(BM41_high_response))


# Loop through each biomarker and perform a t-test
for (i in 1:ncol(BM41_high_response)) {
  biomarker3 <- BM41_high_response[, i]
  biomarker4 <- BM41_low_response[, i]
  ttest_result <- t.test(biomarker3, biomarker4)
  pvalue <- ttest_result$p.value
  t_stat <- ttest_result$statistic
  cohens_d <- abs(mean(biomarker3) - mean(biomarker4)) / sd(c(biomarker3, biomarker4))
  biomarker_name <- colnames(BM41_high_response)[i]
  print(paste("P-value for", biomarker_name, ":", pvalue))
  print(paste("t-statistic for", biomarker_name, ":", t_stat))
  print(paste("Cohen's d for", biomarker_name, ":", cohens_d))
  p_values_BM41[i] <- pvalue
  t_statistic_BM41[i] <- t_stat
  effect_size_BM41[i] <- cohens_d
  
}



# new dataframe Biomarker and p-values 

biomarker_BM41 <- data.frame(Parameter = names(bm41_df)[3:18],p_value = p_values_BM41, t_statistic = t_statistic_BM41, effect_size = effect_size_BM41)

biomarker_BM41$Treatment <- "BM41"


# DF Biomarker and p-values 

Treatment_pvalues <- bind_rows(biomarker_alutard, biomarker_BM41)


Roc_curve_analysis <- read_excel("~/Clinical_Trial_Master_project/Clinical_Trail_BM41_Analysis/data_generated/Roc_curve_analysis.xlsx")
View(Roc_curve_analysis)


roc_timepoint3 <- Roc_curve_analysis %>% 
  filter(Timepoint == 3) 

roc_timepoint3 <- rename(roc_timepoint3, Parameter = Measurement)


### BIOMARKER DF-------------------- 


Biomarker_df <- merge.data.frame(Treatment_pvalues,roc_timepoint3, by= c("Parameter", "Treatment"), all.x = TRUE )
    
view(Biomarker_df)


# Calculate adjusted p-value using the Benjamini-Hochberg method
Biomarker_df$adj_pvalue <- p.adjust(Biomarker_df$p_value, method = "BH")

# Take the negative logarithm of the adjusted p-value
Biomarker_df$log_adj_pvalue <- -log10(Biomarker_df$adj_pvalue)




##Visulizing differences in high vs. low reponder-----------------

# Create a variable to label significant differences
Biomarker_df$Significant <- ifelse(Biomarker_df$adj_pvalue < 0.05, "*", "")



# Plot with significant differences marked in red
ggplot(Biomarker_df, aes(x = log_adj_pvalue, y = Parameter)) +
  geom_bar(aes(fill = Significant), stat = "identity", position = position_dodge()) +
  facet_wrap(~ Treatment) +
  scale_fill_manual(values = c("grey", "red")) +
  labs(fill = "Significance: adj_p-value < 0.05")+
  theme_bw()+
  labs(title = " High vs.Low response based on proxy ",x = " -log(adj_p-value)", y = " Biomarker")



# create a bar plot with error bars
ggplot(Biomarker_df, aes(x = adj_pvalue, y = Parameter)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Treatment)

ggplot(Biomarker_df, aes(x = adj_pvalue, y = Parameter, fill = Treatment)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("grey60", "grey80")) +
  facet_wrap(~ Treatment, scales = "free_y") +
  labs(x = "Adjusted p-value", y = "") +
  theme_minimal() +
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.y = element_text(size = 8),
        legend.position = "none")





##########Visulizing AUC vs. p_value ############

# Define a threshold for AUC and  p_value

AUC_threshold <- 0.8
p_value_threshold <- 0.05

# Create a new column indicating whether the measurement is high AUC and low p_value
Biomarker_df$Significant <- ifelse(Biomarker_df$AUC >= AUC_threshold & Biomarker_df$p_value <= p_value_threshold, "yes", "no")



# filter and print highlighted parameters
significant_param <- Biomarker_df %>%
  filter(Significant == "yes") %>%
  select(Parameter, AUC,p_value,Treatment)




# Create the scatter plot
ggplot(Biomarker_df, aes(x = AUC, y = p_value, color = Significant)) +
  geom_point() +
  scale_color_manual(values = c("no" = "black", "yes" = "red")) +
  facet_wrap(~ Treatment) + 
  theme_bw()+
  labs(title = "AUC vs. p value Biomarker in high a low response degree group (based on rBet v 1 sIgG)",
       x = "AUC Prediction power Treatment vs. Placebo",
       y = "p value high vs. low response to treatment")




# Create a scatter plot with t-statistic as color

ggplot(Biomarker_df, aes(x = AUC, y = log_adj_pvalue, color= t_statistic)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  facet_wrap(~ Treatment) + 
  theme_bw() +
  xlim(0, 1) +
  labs(title = "AUC vs. adjusted p-value Biomarker in high a low response degree group (based on rBet v 1 sIgG)",
       x = "AUC Prediction power Treatment vs. Placebo",
       y = "-log10(adjusted p-value) high vs. low response to treatment",
       color = "t-statistic") +
  guides(color = guide_colorbar(title = "t-statistic"))+
  geom_text_repel(aes(label = ifelse(AUC > 0.8, Parameter, "")))


####

ggplot(Biomarker_df, aes(x = AUC, y = log_adj_pvalue, color= t_statistic > 2)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red")) +
  facet_wrap(~ Treatment) + 
  theme_bw() +
  xlim(0, 1) +
  labs(title = "AUC vs. adjusted p-value Biomarker in high a low response degree group (based on rBet v 1 sIgG)",
       x = "AUC Prediction power Treatment vs. Placebo",
       y = "-log10(adjusted p-value) high vs. low response to treatment",
       color = "Significance") +
  guides(color = guide_legend(title = NULL))+
  geom_text_repel(aes(label = ifelse(AUC > 0.8, Parameter, "")))



####

ggplot(Biomarker_df, aes(x = AUC, y = log_adj_pvalue, color = sign(t_statistic))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red"), name = "t-statistic", labels = c("Negative", "Positive")) +
  facet_wrap(~ Treatment) + 
  theme_bw() +
  xlim(0, 1) +
  labs(title = "AUC vs. adjusted p-value Biomarker in high a low response degree group (based on rBet v 1 sIgG)",
       x = "AUC Prediction power Treatment vs. Placebo",
       y = "-log10(adjusted p-value) high vs. low response to treatment") +
  geom_text_repel(aes(label = ifelse(AUC > 0.8, Parameter, "")), color = "black") +
  theme(legend.position = "bottom")




# Correlation t=1 and t=3 ------------------------------

view(df_a)

# Select columns for timepoint 1 and 3

tp1_cols <- df_a[df_a$Timepoint == 1, 3:18]

tp3_cols <- df_a[df_a$Timepoint == 3, 3:18]
library(dplyr)

# Rename the column names and add "T3" in front of each column name
colnames(tp3_cols) <- paste0("T3_", colnames(tp3_cols))

# Print the updated column names
colnames(tp3_cols)


# Calculate correlation matrix

cor_matrix <- cor(tp1_cols, tp3_cols)

cor_matrix_spearman <- cor(tp1_cols,tp3_cols, method = "spearman")

# View correlation matrix

cor_matrix


library(corrplot)

corrplot(cor_matrix)

corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", addCoef.col = "black", tl.col = "black")

library(ComplexHeatmap)
library(RColorBrewer)
library(grDevices)

corrplot(cor_matrix_spearman)
col <- colorRampPalette(c("blue", "white", "red"))(n = 101)

# create the plot
corrplot(cor_matrix_spearman, col = col, method = "number", is.corr = TRUE, tl.col = "black",add = TRUE, xlab = " Timepoint 3", ylab = "Timepoint 3", 
         title = "Correlation Matrix")

corrplot(cor_matrix_spearman, col = col, is.corr = TRUE, tl.col = "black")



Heatmap(
  cor_matrix,
  name = "Correlation",
  column_title = "Timepoint 3",
  row_title = "Timepoint 1",
  col = colorRampPalette(c("blue", "white", "red"))(100))

Heatmap(cor_matrix_spearman, name= "Spearman Correlation", column_title = "Timepoint 3",
        row_title = "Timepoint 1",
        col = colorRampPalette(c("blue", "white", "red"))(100))




#################

library(ggplot2)
library(reshape2)

# Create a data frame from the correlation matrix
cor_df <- melt(cor_matrix_spearman)

# Set up the plot
ggplot(data = cor_df, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), na.value = "white") +
  coord_fixed() +
  labs(x = "Timepoint 3", y = "Timepoint 1", title = "Spearman Correlation")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## Corr t=2 t=3 --------------------
# Select columns for timepoint 1 and 3

tp2_cols <- df_a[df_a$Timepoint == 2, 3:18]


# Calculate correlation matrix

cor_matrix <- cor(tp2_cols, tp3_cols)

cor_matrix_spearman2 <- cor(tp2_cols,tp3_cols, method = "spearman")

corrplot(cor_matrix_spearman2, col = col, method = "number", is.corr = TRUE, tl.col = "black")
corrplot(cor_matrix_spearman2, col = col,  is.corr = TRUE, tl.col = "black")


cor_matrix_spearman2_m <- melt(cor_matrix_spearman2)



Heatmap(cor_matrix_spearman2, name= "Spearman Correlation", column_title = "Timepoint 3",
        row_title = "Timepoint 2",
        col = colorRampPalette(c("blue", "white", "red"))(100))

# Set up the plot
ggplot(data = cor_matrix_spearman2_m, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), na.value = "white") +
  coord_fixed() +
  labs(x = "Timepoint 3", y = "Timepoint 2", title = "Spearman Correlation")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





#%Diff t1 and t3---------------------------------

#Trandforming the data

                          
df_wide_Betv1_sIgG <- pivot_wider(df, id_cols = c("Group", "ID"), names_from = Timepoint, values_from = c( "rBetv1_sIgG"))
df_wide_Betv1_sIgE <- pivot_wider(df, id_cols = c("Group", "ID"), names_from = Timepoint, values_from = c( "rBetv1_sIgE"))

colnames(df_wide_Betv1_sIgG)
  

# Create plot object


# Create new columns
df_wide_Betv1_sIgG <- df_wide_Betv1_sIgG %>%
  mutate(start = `1`,
         end = `3`)

df_wide_Betv1_sIgE <- df_wide_Betv1_sIgE %>%
  mutate(start = `1`,
         end = `3`)


# Create plot object

ggplot(df_wide_Betv1_sIgG, aes(y=ID)) +
  geom_point(aes(x=start, color= "T1"))+
  geom_point(aes(x=end, color= "T3"))+
  facet_grid(~ Group)+
  geom_segment(aes(x = `1`, y = ID, xend = `3`, yend = ID),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt",
               color = "black")+
labs(x = "rBet v 1 sIgG", y = "ID", title = "rBet v 1 sIgG during treatment")

#reareange ID order based in T1 

df_wide_Betv1_sIgG <- df_wide_Betv1_sIgG %>%
  arrange(`1`) %>%
  mutate(ID = factor(ID, levels=unique(ID), ordered = TRUE))

df_wide_Betv1_sIgE <- df_wide_Betv1_sIgE %>%
  arrange(`1`) %>%
  mutate(ID = factor(ID, levels=unique(ID), ordered = TRUE))


# Create plot object
ggplot(df_wide_Betv1_sIgG, aes(y=ID)) +
  geom_point(aes(x=start, color= "T1"))+
  geom_point(aes(x=end, color= "T3"))+
  facet_grid( row= "Group", scales = "free_y")+
  geom_segment(aes(x = `1`, y = ID, xend = `3`, yend = ID),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt",
               color = "black")+
  labs(x = "rBet v 1 sIgG", y = "ID", title = "rBet v 1 sIgG during treatment")


# Create plot object
ggplot(df_wide_Betv1_sIgE, aes(y=ID)) +
  geom_point(aes(x=start, color= "T1"))+
  geom_point(aes(x=end, color= "T3"))+
  facet_grid( row= "Group", scales = "free_y")+
  geom_segment(aes(x = `1`, y = ID, xend = `3`, yend = ID),
               arrow = arrow(length = unit(0.2, "inches")),
               lineend = "butt",
               color = "black") +
  labs(x = "rBet v 1 sIgE", y = "ID", title = "rBet v 1 sIgE during treatment")


#Heatmap-------------------


# Convert data frame to matrix
df_corr <- df_a[, 3:18]

# Calculate correlation matrix
corr_mat <- cor(df_corr, method = "spearman")

mat_df_a <- as.matrix(df_a)

# Select columns for correlation matrix

# Define the color scale for the heatmap

my_palette <- brewer.pal(11, "RdBu")


# Create heatmap
Heatmap(corr_mat, 
        show_row_names = TRUE, 
        show_column_names = TRUE,
        clustering_distance_columns = "correlation",
        clustering_distance_rows = "correlation",
        name = "Correlation Heatmap",
        top_annotation = HeatmapAnnotation(df_a, which = c("Group", "Timepoint")),
        row_title = "TimeID",
        col = my_palette)



df_subset <- df_a[,c(3:18)]

#### Group anno

treatment_annotation <- df_a %>% select(Group)  
row.names(treatment_annotation) <- df_a$TimeID

#### Timepoint Anno

timepoint_annotation <- df_a %>% select(ID,Timepoint)  
row.names(timepoint_annotation) <- df_a$TimeID


# view the converted vector
timepoint_annotation
treatment_annotation

# create matrix
mat <- as.matrix(df_subset)



# create Heatmap object
ha <- Heatmap(
  mat,
  cluster_rows = TRUE,
  show_row_names = TRUE
)


# Customize Heatmap object
ha <- Heatmap(
  mat ,
  split = df_a$Timepoint,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  col = my_palette,
  row_title = "Treatment",
  column_title = "Biomarker",
  name = "Value",
  width = unit(12, "cm"),
  height = unit(15, "cm"),
  row_dend_reorder = TRUE,
  column_names_rot = 75,
  column_names_gp = grid::gpar(fontsize = 8),
  row_labels = as.matrix(treatment_annotation))


# Draw heatmap
draw(ha)



dev.off()













































































# Tidy the dataset to use pivot 

combined_df <- combined_df %>%
  mutate(response_degree = recode(response_degree, "high" = 1, "low" = 2, "intermediate" = 3))


ndf <- subset(combined_df, select = -c(Group,ID, Timepoint,cluster_3K_a,cluster_5K_a, cluster_4K_a))

# Set TimeID as row names
rownames(ndf) <- ndf$TimeID 
ndf <- ndf[, -17] 
str(ndf)

# Testing the parameter 

#ALUTARD 

# filter for rows with Alutard in the row name


alutard_data <- ndf %>%
  rownames_to_column(var = "TimeID") %>%
  filter(str_detect(TimeID, "Alutard")) %>%
  select(-c(TimeID)) # remove unnecessary columns



# filter for rows with BM41 in the row name
bm41_data <- ndf %>%
  rownames_to_column(var = "TimeID") %>%
  filter(str_detect(TimeID, "BM41")) %>%
  select(-c(TimeID)) # remove unnecessary columns

# Pivot the data into a longer format

alutard_data_longer <- alutard_data %>% 
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  separate(parameter, into = c("parameter", "response_degree"), sep = "_") %>%
  select(-response_degree) # remove response_degree column


# perform t-test for each parameter between response_degree 1 and 2
alutard_results <- alutard_data %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")
  group_by(Parameter) %>%
  summarize(p_value = t.test(Value[response_degree == 1], Value[response_degree == 2])$p.value)



























































































# Calculate t-test for each parameter between high and low response groups for each treatment
t_test_results <- ndf %>%
  filter(response_degree != "3") %>%
  pivot_longer(cols = -c(Group, response_degree)) %>%
  group_by(group,name, response_degree) %>%
  summarise(p.value = t.test(value ~ response_degree)$p.value,
            diff = diff(t.test(value ~ response_degree)$estimate)) %>%
  pivot_wider(names_from = response_degree, values_from = c(p.value, diff))


# Format table to show parameter name, Alutard and BM41 p-values, Alutard and BM41 parameter differences
table_output <- t_test_results %>%
  pivot_longer(cols = c(p.value_high, diff_high, p.value_low, diff_low), 
               names_to = c(".value", "response_degree"), 
               names_sep = "_") %>%
  pivot_wider(names_from = group, values_from = c(p.value, diff)) %>%
  select(name, 
         Alutard_p_value_high, Alutard_diff_high, 
         Alutard_p_value_low, Alutard_diff_low, 
         BM41_p_value_high, BM41_diff_high, 
         BM41_p_value_low, BM41_diff_low) %>%
  mutate(across(where(is.numeric), round, digits = 3)) 

# Display table
table_output





# Create a plot of the parameter differences for each treatment
ggplot(t_test_results, aes(x = diff, y = reorder(name, diff), color = group)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = diff - conf.low, xmax = diff + conf.high), height = 0.2) +
  labs(title = "Parameter Differences between High and Low Response Groups",
       x = "Difference",
       y = "Parameter",
       color = "Treatment") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top")
















#rBet v 1 sIgG is also significantly different in the clusters 1 and 2 


# Create a subset of the data frame with only the rows that have cluster_5k_a = 1
high_response <- subset(df_a, cluster_5K_a == '1')

# Create a subset of the data frame with only the rows that have cluster_5k_a = 2
low_response <- subset(df_a, cluster_5K_a == '2')

# make a loop to test every measurements in high and low response group 


for (col in names(df_a)[3:18]) {
  ttest_result <- t.test(high_response[[col]], low_response[[col]])
  df_t <- as.data.frame(print(paste0(col, ": t-value=", ttest_result$statistic, ", p-value=", ttest_result$p.value)))
}

# Perform a t-test to compare the r_bet_v_1_s_ig_g levels between the two groups

t.test(rBet.v.1.sIgG ~ cluster_5K_a, data = df_a, subset = cluster_5K_a %in% c('1', '2'))














 




