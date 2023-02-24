#ROC Analysis


library(pROC)
library(tidyverse)
library(readxl)
library(openxlsx)

library(ggplot2)
library(plotROC)
library(dplyr)
library(tidyr)

library(ROCR)
library(tableone)
library(summarytools)
library(ggfortify)
library(readxl)
library(lubridate)
library(xlsx)


dataset_tp <- read_xlsx("C:/Users/Tomic/OneDrive/Dokumente/R Data Science/Tables/dataset_timepoints.xlsx", col_names = TRUE)
df <- dataset_tp

#log transform and add 0.001
df[,6:18] <- log(df[,6:18] + 0.001)

#scale all 3:18
df[,3:18] <- scale(df[,3:18])

df$Timepoint <- as.factor(df$Timepoint)



#########################Labels######################################

# 
# # Create a data frame with only the Alutard and Placebo groups
# 
# df_subset_AP <- df[df$Group %in% c('Alutard', 'Placebo'),header=TRUE]
# 
# # Create a new column in the data frame to indicate the treatment group
# df_subset_AP$treatment <- ifelse(df_subset_AP$Group == "Alutard", 1, 0) %>%
#                           as.factor()
# 
# df_subset_AP$treatment
# 
# write.csv(df_subset_AP, file = "df_subset_AP_treatment.csv", row.names = FALSE)
# 
# #Subset the data according to Timepoints
# 
# df_subset_AP_timepoint1 <- subset(df_subset_AP, Timepoint == '1', header=TRUE)
# df_subset_AP_timepoint2 <- subset(df_subset_AP, Timepoint == '2', header=TRUE)
# df_subset_AP_timepoint3 <- subset(df_subset_AP, Timepoint == '3',header=TRUE)
# 
# 
# df_1 <- df_subset_AP_timepoint1[,3:19]
# df_2 <- df_subset_AP_timepoint2[,3:19]
# df_3 <- df_subset_AP_timepoint3[,3:19]




################### AUC for each measurement at differnt timepoints###############

(measurements <- colnames(df)[3:18])

res <- list()
res.roc <- list()


tx <- 1
treatx <- 'Alutard'
mx <- measurements[1]
for(tx in unique(df$Timepoint)){
  for(mx in measurements){
    for(treatx in setdiff(unique(df$Group), 'Placebo')){
      x <- df %>% 
        filter(Group %in% c(treatx, 'Placebo')) %>%
        filter(Timepoint == tx)
      pred <- prediction(predictions=x[[mx]], labels = as.numeric(x$Group == treatx))
      perf <- performance(pred, "auc")
      res[[length(res) + 1]] <- data.frame(
        Timepoint = tx,
        Measurement=mx,
        Treatment=treatx,
        AUC=attr(perf, "y.values")[[1]]
        )
      
      roc <- performance(pred, "tpr", "fpr")
      res.roc[[length(res) + 1]] <- data.frame(
        Timepoint = tx,
        Measurement=mx,
        Treatment=treatx,
        tpr=attr(roc, "y.values")[[1]],
        fpr=attr(roc, "x.values")[[1]]
      )
    }
  }
}



require(tidyverse)
bind_rows(res)
bind_rows(res.roc)

roc_df <- bind_rows(res) 

write.xlsx(roc_df, "C:\\Users\\Tomic\\OneDrive\\Dokumente\\Clinical_Trial_Master_project\\Clinical_Trail_BM41_Analysis\\data_generated\\Roc_curve_analysis.xlsx")


Top_variables <- roc_df %>% 
  filter(AUC > 0.8)


roc_res_df <- bind_rows(res.roc)

roc_df$Timepoint <- as.factor(roc_df$Timepoint)



# Create a plot of the AUC values for each measurement, separated by treatment

ggplot(Top_variables, aes(x = Timepoint, y = AUC, color = Treatment)) +
  geom_point() +
  facet_wrap(~ Measurement) +
  theme_bw() +
  labs(title = "AUC values for each measurement, separated by treatment",
       x = "Timepoint",
       y = "AUC")





ggplot(roc_df, aes(x= Timepoint, y=AUC)) + facet_wrap(~Measurement, scales = "free") + theme_bw()

ggplot(data = roc_df, aes(x = Timepoint, y = AUC, fill = Measurement)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Measurement, scales = "free")




###########Correlation Heatmap ###############

library(ComplexHeatmap)

Heatmap(roc_df)


###################ROC Curves################################


Competition_FAB_df <- filter(roc_res_df, Measurement == "Competition FAB FITC" & Treatment == "Alutard")


ggplot(data= Competition_FAB_df, aes(x=fpr,y=tpr)) + geom_line() + facet_grid( row= "Timepoint") + theme_bw() + ggtitle("Competition_FAB ROC")
   





Inhibition_FAB_df <- filter(roc_res_df, Measurement == "Inhibition FAB ITC" & Treatment == "Alutard")


ggplot(data= Inhibition_FAB_df, aes(x=fpr,y=tpr)) + geom_line() + facet_grid( row= "Timepoint") + theme_bw() + ggtitle("Inhibition FAB ROC")





Inhibtion_RBL_df <- filter(roc_res_df, Measurement == "Inhibition RBL" & Treatment == "Alutard")


ggplot(data= Inhibtion_RBL_df, aes(x=fpr,y=tpr)) + geom_line() + facet_grid( row= "Timepoint") + theme_bw() + ggtitle("Inhibition RBL ROC")
  




rBet_v_1_sIgG_df <- filter(roc_res_df, Measurement == "rBet v 1 sIgG" & Treatment == "Alutard")


ggplot(data= rBet_v_1_sIgG_df, aes(x=fpr,y=tpr)) + geom_line() + facet_grid( row= "Timepoint") + theme_bw() + ggtitle("rBet v 1 sIgG ROC")

 






##################### Logistic regression timepoint 1######################


logistic_regression_results_1 <- apply(df_subset_AP_timepoint1[,3:18], 2, function(x) {
  
  logistic_regression <- glm(treatment ~ x, data = df_subset_AP_timepoint1, family = "binomial")
  return(summary(logistic_regression))
  glm.fit=glm(treatment ~ x, data = df_subset_AP_timepoint1, family = "binomial")
  
  
})

# Create a table with the coefficients of each model

coefficients <- lapply(logistic_regression_results_1, coef)
coefficients_table <- as.data.frame(coefficients)


# Lowest p value at timepoint 1 =  "Standard FAB"

min_p_value_1 <- min(sapply(logistic_regression_results_1, function(x) x$coefficients[2,4]))
variable_with_lowest_p_value_1 <- names(which.min(sapply(logistic_regression_results_1, function(x) x$coefficients[2,4])))

variable_with_lowest_p_value_1

+++++++++++++++++++++++++++++++++++++
  
  ####effect size with effectsize package ? #Serology Birch, Standard FAB
  

 

  # Calculate the AUC for each model
  
  auc <- lapply(logistic_regression_results_1, function(x) {
    pred <- predict(x, type = "response")
    perf <- performance(prediction(pred, df_1$treatment), "auc")
    attr(perf, "y.values")[[1]]
  })
  
# Plot a ROC curve for each model
roc_curves <- lapply(models, function(x) {
  pred <- predict(x, type = "response")
  perf <- performance(prediction(pred, df$treatment), "tpr", "fpr")
  plot(perf, col = "blue", lwd = 2, main = paste("ROC Curve for", names(models)[i]))
})

# Show AUC for each variable
auc_table <- as.data.frame(auc)
print(auc_table)#
  
  
  

###################################Logistic regression for each measurement at timepoint 2######### 


logistic_regression_results_2 <- apply(df_subset_AP_timepoint2[,3:18], 2, function(x) {
  
  logistic_regression_ <- glm(treatment ~ x, data = df_subset_AP_timepoint2, family = "binomial")
  return(summary(logistic_regression_))
})

logistic_regression_results_2 

# Lowest p value at timepoint 1 =  "rBet v 1 sIgG"
min_p_value <- min(sapply(logistic_regression_results_2, function(x) x$coefficients[2,4]))
variable_with_lowest_p_value <- names(which.min(sapply(logistic_regression_results_2, function(x) x$coefficients[2,4])))

variable_with_lowest_p_value




######################Logistic regression for each measurement at timepoint 3###### not sure if this is really right 

logistic_regression_results_3 <- apply(df_subset_AP_timepoint3[,3:18], 2, function(x) {
  
  logistic_regression_3 <- glm(treatment ~ x, data = df_subset_AP_timepoint3, family = "binomial")
  return(summary(logistic_regression_3))
})

logistic_regression_results_3  

# Lowest p value at timepoint 3

min_p_value <- min(sapply(logistic_regression_results_3, function(x) x$coefficients[2,4]))
variable_with_lowest_p_value <- names(which.min(sapply(logistic_regression_results_3, function(x) x$coefficients[2,4])))

variable_with_lowest_p_value





#####################All variables logistic model #####################


#first split it in train and test ? 

df1_m <- as.matrix(df_1)



require(glmnet)

glmnet(treatment ~ ., data = df1_m, family = binomial)


fit <- cv.glmnet(keep=TRUE)
fit@fit.preval

prediction(predictions = fit@fit.preval)






model_1 <- glm(treatment ~ ., data = df_1, family = binomial)

summary(model_1)

pvals_1 <- summary(model_1)$coefficients[,4]

# Fit a logistic regression model to predict treatment
model_1 <- glm(treatment ~ ., data = df_1, family = binomial)


# Calculate the AIC for the model
aic <- AIC(model_1)

# Calculate the AUC for the model
auc <- roc(df_1$treatment, model_1$fitted.values)

plot(auc)

roc(df_1$treatment, model_1$fitted.values, plot=TRUE)

# Calculate the accuracy of the model
accuracy <- mean(df_1$treatment == round(model_1$fitted.values))

# Print the results
cat("AIC:", aic, "\n")
cat("AUC:", auc$auc, "\n")
cat("Accuracy:", accuracy, "\n")

