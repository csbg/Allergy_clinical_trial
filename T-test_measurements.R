
library(dplyr)
library(openxlsx)
library(tidyverse)
library(broom)



df_a <- read.xlsx("C:/Users/Tomic/OneDrive/Dokumente/Clinical_Trial_Master_project/Clinical_Trail_BM41_Analysis/data_generated/df.xlsx")

colnames(df_a)

# Subset data for Timepoint 3 and Alutard treatment
alutard_df <- df_a %>% 
  filter(Timepoint == 3, Group == "Alutard") 



# Subset data for Timepoint 3 and BM41 treatment
bm41_df <- df_a %>% 
  filter(Timepoint == 3, Group == "BM41") 



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


write.xlsx(combined_df, "C:\\Users\\Tomic\\OneDrive\\Dokumente\\Clinical_Trial_Master_project\\Clinical_Trail_BM41_Analysis\\data_generated\\rBet_v_1_sIgG_Response_table.xlsx") 



# Calculate t-test for each parameter between high and low response groups for each treatment
t_test_results <- combined_df %>%
  filter(response_degree != "intermediate") %>%
  pivot_longer(cols = -c(Group, response_degree)) %>%
  group_by(Group,response_degree) %>%
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


##################Direct Analysis rBetv1 sIgG as proxy for response #################


#Create a new column called "Response"

df$Response <- NA
# Select the rows with timepoint 3 and Alutard and BM41
df_subset <- df[df$timepoint == '3' & df$Alutard == 'Alutard' & df$BM41 == 'BM41',]

# Sort the data frame by r_bet_v_1_s_ig_g
df_subset <- df_subset[order(df_subset$r_bet_v_1_s_ig_g),]

# Select the top 5 rows
top_5 <- df_subset[1:5,]

# Assign the value 'High Response' to the top 5 rows
top_5$Response <- 'High Response'

# Select the bottom 5 rows
bottom_5 <- df_subset[(nrow(df_subset)-4):nrow(df_subset),]

# Assign the value 'Low Response' to the bottom 5 rows
bottom_5$Response <- 'Low Response'

# Select the remaining rows
intermediate <- df_subset[6:(nrow(df_subset)-5),]

# Assign the value 'Intermediate Response' to the remaining rows
intermediate$Response <- 'Intermediate Response'

# Combine the three data frames
df_subset <- rbind(top_5, intermediate, bottom_5)

# Replace the original data frame with the new one
df[df$timepoint == '3' & df$Alutard == 'Alutard' & df$BM41 == 'BM41',] <- df_subset




# Get the highest 5 and lowest 5 values from r_bet_v_1_s_ig_g in Alutard and BM41 group 

# Subset the data according to Treatment

df_a_Alutard <- subset(df_a, Group == 'Alutard')

df_a_BM41 <- subset(df_a, Group == 'BM41')

# sort 

#Alutard r_bet_v_1_s_ig_g

high_5_proxy <- sort(df_a_Alutard$rBet.v.1.sIgG, decreasing = TRUE)[1:5]

low_5_proxy <- sort(df_a_Alutard$rBet.v.1.sIgG)[1:5]

# Test if the highest 5 and lowest 5 values are significantly different
t.test(high_5_proxy, low_5_proxy)



#BM41  r_bet_v_1_s_ig_g

high_5_proxy_b <- sort(df_a_BM41$rBet.v.1.sIgG, decreasing = TRUE)[1:5]

low_5_proxy_b <- sort(df_a_BM41$rBet.v.1.sIgG)[1:5]

# Test if the highest 5 and lowest 5 values are significantly different
t.test(high_5_proxy_b, low_5_proxy_b)



# display the differences between this top 5 and bottom 5 r_bet_v_1_s_ig_g in other measurements 




