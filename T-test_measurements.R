
library(dplyr)
library(openxlsx)
library(tidyverse)
library(broom)
library(tidyr)
library(ggplot2)
library(readxl)
library(tableone)




# Load data ---------------------------------------------------------------


df_a <- read.xlsx("data_generated/df_scaled.xlsx")

colnames(df_a)

# Subset data for Timepoint 3 and Alutard treatment
alutard_df <- df_a %>% 
  filter(Timepoint == 3, Group == "Alutard") 



# Subset data for Timepoint 3 and BM41 treatment
bm41_df <- df_a %>% 
  filter(Timepoint == 3, Group == "BM41") 


# Direct Analysis rBetv1 sIgG as proxy for response -----------------------

# Get the top 5 and bottom 5 values of rBetv sIgG for Alutard treatment
alutard_top5 <- alutard_df %>% 
  arrange(desc(rBet.v.1.sIgG)) %>% 
  slice(1:5)


alutard_bottom5 <- alutard_df %>% 
  arrange(rBet.v.1.sIgG) %>% 
  slice(1:5)


# Assign response levels based on rBetv sIgG values for Alutard treatment

alutard_df <- alutard_df %>% 
  mutate(response_degree= case_when(
    rBet.v.1.sIgG %in% alutard_top5$rBet.v.1.sIgG ~ "high",
    rBet.v.1.sIgG %in% alutard_bottom5$rBet.v.1.sIgG ~ "low",
    TRUE ~ "intermediate"
  ))



# Get the top 5 and bottom 5 values of rBetv sIgG for BM41 treatment
bm41_top5 <- bm41_df %>% 
  arrange(desc(rBet.v.1.sIgG)) %>% 
  slice(1:5)


bm41_bottom5 <- bm41_df %>% 
  arrange(rBet.v.1.sIgG) %>% 
  slice(1:5)

# Assign response levels based on rBetv sIgG values for BM41 treatment
bm41_df <- bm41_df %>% 
  mutate(response_degree = case_when(
    rBet.v.1.sIgG %in% bm41_top5$rBet.v.1.sIgG ~ "high",
    rBet.v.1.sIgG%in% bm41_bottom5$rBet.v.1.sIgG ~ "low",
    TRUE ~ "intermediate"
  ))


# Combine the two dataframes
combined_df <- bind_rows(alutard_df, bm41_df)

CreateTableOne( data=combined_df, strata="response_degree")


write.xlsx(combined_df, "C:\\Users\\Tomic\\OneDrive\\Dokumente\\Clinical_Trial_Master_project\\Clinical_Trail_BM41_Analysis\\data_generated\\rBet_v_1_sIgG_Response_table.xlsx") 


########################Alutard##################### 

#Subset data 

alutard_high_response <- alutard_top5[3:18]

alutard_low_response <- alutard_bottom5[3:18]

# Create an empty vector to store the p-values
p_values_alutard <- vector(length = ncol(alutard_high_response))

# Loop through each biomarker and perform a t-test

(measurements <- colnames(alutard_df)[3:18])
treatx <- 'Alutard'

for (i in 1:ncol(alutard_high_response)) {
  biomarker1 <- alutard_high_response[, i]
  biomarker2 <- alutard_low_response[, i]
  ttest_result <- t.test(biomarker1, biomarker2)
  pvalue <- ttest_result$p.value
  biomarker_name <- colnames(alutard_high_response)[i]
  print(paste("P-value for", biomarker_name, ":", pvalue))
  p_values_alutard[i] <- pvalue
  
  
}

# Create empty vectors to store the p-values, t-statistics, and effect sizes
p_values_alutard <- numeric(ncol(alutard_high_response))
t_statistic_alutard <- numeric(ncol(alutard_high_response))
effect_size_alutard <- numeric(ncol(alutard_high_response))


# Loop through each biomarker and perform a t-test
measurements <- colnames(alutard_df)[3:18]
treatx <- 'Alutard'

for (i in 1:ncol(alutard_high_response)) {
  biomarker1 <- alutard_high_response[, i]
  biomarker2 <- alutard_low_response[, i]
  ttest_result <- t.test(biomarker1, biomarker2)
  pvalue <- ttest_result$p.value
  tstatistic <- ttest_result$statistic
  cohens_d <- abs(mean(biomarker1) - mean(biomarker2)) / sd(c(biomarker1, biomarker2))
  biomarker_name <- colnames(alutard_high_response)[i]
  print(paste("P-value for", biomarker_name, ":", pvalue))
  print(paste("t-statistic for", biomarker_name, ":", t_stat))
  print(paste("Cohen's d for", biomarker_name, ":", cohens_d))
  p_values_alutard[i] <- pvalue
  t_statistic_alutard[i] <- tstatistic
  effect_size_alutard[i] <- cohens_d
  
}




# new dataframe with the Biomarker and p-values 

biomarker_alutard <- data.frame(Parameter = names(alutard_df)[3:18], p_value = p_values_alutard, t_statistic = t_statistic_alutard, effect_size = effect_size_alutard)

biomarker_alutard$Treatment <- "Alutard"



#######################BM41############################## 

#Subset data 

BM41_high_response <- bm41_top5[3:18]

BM41_low_response <- bm41_bottom5[3:18]

# Create an empty vector to store the p-values
p_values_BM41 <- vector(length = ncol(BM41_high_response))

# Loop through each biomarker and perform a t-test

for (i in 1:ncol(BM41_high_response)) {
  biomarker3 <- BM41_high_response[, i]
  biomarker4 <- BM41_low_response[, i]
  ttest_result <- t.test(biomarker3, biomarker4)
  pvalue <- ttest_result$p.value
  biomarker_name <- colnames(BM41_high_response)[i]
  print(paste("P-value for", biomarker_name, ":", pvalue))
  p_values_BM41[i] <- pvalue
}



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



# new dataframe with the Biomarker and p-values 

biomarker_BM41 <- data.frame(Parameter = names(bm41_df)[3:18],p_value = p_values_BM41, t_statistic = t_statistic_BM41, effect_size = effect_size_BM41)

biomarker_BM41$Treatment <- "BM41"




# Combine the two dataframes

Treatment_pvalues <- bind_rows(biomarker_alutard, biomarker_BM41)



Roc_curve_analysis <- read_excel("~/Clinical_Trial_Master_project/Clinical_Trail_BM41_Analysis/data_generated/Roc_curve_analysis.xlsx")
View(Roc_curve_analysis)


roc_timepoint3 <- Roc_curve_analysis %>% 
  filter(Timepoint == 3) 

roc_timepoint3 <- rename(roc_timepoint3, Parameter = Measurement)


# Umbenennen der Spaltennamen des zweiten Dataframes


Biomarker_df <- merge.data.frame(Treatment_pvalues,roc_timepoint3, by= c("Parameter", "Treatment"), all.x = TRUE )
    
view(Biomarker_df)




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


# Calculate adjusted p-value using the Benjamini-Hochberg method
Biomarker_df$adj_pvalue <- p.adjust(Biomarker_df$p_value, method = "BH")

# Take the negative logarithm of the adjusted p-value
Biomarker_df$log_adj_pvalue <- -log10(Biomarker_df$adj_pvalue)

# Create a scatter plot with t-statistic as color
ggplot(Biomarker_df, aes(x = AUC, y = log_adj_pvalue)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_wrap(~ Treatment) + 
  theme_bw() +
  xlim(0, 1) +
  labs(title = "AUC vs. adjusted p-value Biomarker in high a low response degree group (based on rBet v 1 sIgG)",
       x = "AUC Prediction power Treatment vs. Placebo",
       y = "-log10(adjusted p-value) high vs. low response to treatment",
       color = "t-statistic") +
  guides(color = guide_colorbar(title = "t-statistic"))


# Correlation t=1 and t=3 ------------------------------

view(df_a)

# Select columns for timepoint 1 and 3
tp1_cols <- df_a[df_a$Timepoint == 1, 3:18]

tp3_cols <- df_a[df_a$Timepoint == 3, 3:18]

# Calculate correlation matrix
cor_matrix <- cor(tp1_cols, tp3_cols)

# View correlation matrix
cor_matrix



library(corrplot)



corrplot(cor_matrix)

corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", addCoef.col = "black", tl.col = "black")

library(ComplexHeatmap)
library(RColorBrewer)
library(grDevices)



Heatmap(
  cor_matrix,
  name = "Correlation",
  column_title = "Timepoint 3",
  row_title = "Timepoint 1",
  col = colorRampPalette(c("blue", "white", "red"))(100),
  
)























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














 




