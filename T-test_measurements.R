



df_a <- read.xlsx("C:/Users/Tomic/OneDrive/Dokumente/Clinical_Trial_Master_project/Clinical_Trail_BM41_Analysis/data_generated/df.xlsx")


# Create a subset of the data frame with only the rows that have cluster_5k_a = 1
high_response <- subset(df_a, cluster_5K_a == '1')

# Create a subset of the data frame with only the rows that have cluster_5k_a = 2
low_response <- subset(df_a, cluster_5K_a == '2')

# Perform a t-test to compare the r_bet_v_1_s_ig_g levels between the two groups

t.test(rBet.v.1.sIgG ~ cluster_5K_a, data = df_a, subset = cluster_5K_a %in% c('1', '2'))


# Get the highest 5 and lowest 5 values from r_bet_v_1_s_ig_g

high_5 <- sort(df_a$rBet.v.1.sIgG, decreasing = TRUE)[1:5]
low_5 <- sort(df_a$rBet.v.1.sIgG)[1:5]

# Test if the highest 5 and lowest 5 values are significantly different
t.test(high_5, low_5)


# make a loop to test every measurements --------------------- bad is not working 

(measurements <- colnames(df_a)[3:18])

mx <- measurements[1]

for(mx in measurements)
  t.test(mx~ cluster_5K_a, data = df_a, subset = cluster_5K_a %in% c('1', '2'))

