library(ggplot2)
library(dplyr)
library(FSA)
library(rcompanion)
library(vegan)
library(tidyverse)

M1_filelist <- list.files("./Mock1_confmatstats/", pattern = ".csv", full.names = TRUE)
M2_filelist <- list.files("./Mock2_confmatstats/", pattern = ".csv", full.names = TRUE)

Mock1 <- list()
Mock2 <- list() 

for (i in seq_along(M1_filelist)) {
  Mock1[[i]] <- read.csv(M1_filelist[[i]], header = TRUE)
}
for (i in seq_along(M2_filelist)) {
  Mock2[[i]] <- read.csv(M2_filelist[[i]], header = TRUE)
}
Mock1_df <- do.call(rbind, Mock1)
Mock2_df <- do.call(rbind, Mock2)
Mock1_df[is.na(Mock1_df)] <- 0
Mock2_df[is.na(Mock2_df)] <- 0

# First create a version which doesn't contain ML classifiers, for later assessment
Mock1_df_noML <- Mock1_df %>% filter(!Classifier %in% c("vtd", "mycoAICNN", "mycoAIBert"))
Mock2_df_noML <- Mock2_df %>% filter(!Classifier %in% c("vtd", "mycoAICNN", "mycoAIBert"))

# Manually set the simNCBI database to NCBI to prevent a third panel during plot faceting
Mock2_df$Database[Mock2_df$Database == "simNCBI"] <- "NCBI"
Mock1_df$Database[Mock1_df$Database == "simNCBI"] <- "NCBI"

# Try transforming F1 score data to achieve normality
# Logit transformation with residual = 0.001 because values are between 0 and 1
Mock1_df$F1_logit <- log((Mock1_df$F1_score + 0.001)/(1 - Mock1_df$F1_score + 0.001))
Mock1_df$F1_arcsine <- asin(sqrt(Mock1_df$F1_score))

hist(Mock1_df$F1_logit) # Skewed to the right
hist(Mock1_df$F1_arcsine) # Also skewed to the right
qqnorm(Mock1_df$F1_logit)
qqnorm(Mock1_df$F1_arcsine)

# Histograms and QQ plots indicate that transformation did not achieve normal distribution. Split data by database - performance (F1_score value) highly dependant on NCBI vs GSD
Mock1_df_GSD <- filter(Mock1_df, Database == "GSD")
Mock1_df_NCBI <- filter(Mock1_df, Database == "NCBI" | Database == "simNCBI")

hist(Mock1_df_GSD$F1_logit)
hist(Mock1_df_NCBI$F1_logit)
hist(Mock1_df_GSD$F1_arcsine)
hist(Mock1_df_NCBI$F1_arcsine)
qqnorm(Mock1_df_GSD$F1_logit)
qqnorm(Mock1_df_NCBI$F1_logit)
qqnorm(Mock1_df_GSD$F1_arcsine)
qqnorm(Mock1_df_NCBI$F1_arcsine)
# Data is still indicating that transformation did not achieve normal distribution. Proceed with non-parametric testing. 

# Univariate pairwise hypothesis testing
M1_F1_Dunn <- dunnTest(Mock1_df$F1_score, Mock1_df$Classifier, method = "sidak")
print(M1_F1_Dunn)
M2_F1_Dunn <- dunnTest(Mock2_df$F1_score, Mock2_df$Classifier, method = "sidak")
print(M2_F1_Dunn)

# Multivariate hypothesis testing with PERMANOVA on full data
M1_F1_adonis <- adonis2(Mock1_df$F1_score ~ as.factor(Mock1_df$Classifier) + as.factor(Mock1_df$Database) + as.factor(Mock1_df$Species) + as.factor(Mock1_df$MinimumSeqQuality), 
                        method = "euclidean", 
                        by = "margin")
summary(M1_F1_adonis)
M2_F1_adonis <- adonis2(Mock2_df$F1_score ~ as.factor(Mock2_df$Classifier) + as.factor(Mock2_df$Database) + as.factor(Mock2_df$Species) + as.factor(Mock2_df$MinimumSeqQuality), 
                        method = "euclidean", 
                        by = "margin")
summary(M2_F1_adonis)
# PERMANOVA doesn't support interaction/nesting terms when using by = "margin". Fit a random forest model to better test the contribution of factors to overall variance. 
#install.packages(c("randomForest", "caret", "vip"))
library(randomForest)
library(caret)
library(vip)

# Convert the explanatory variable columns to factors
Mock1_df$Classifier <- as.factor(Mock1_df$Classifier)
Mock1_df$Database <- as.factor(Mock1_df$Database)
Mock1_df$MinimumSeqQuality <- as.factor(Mock1_df$MinimumSeqQuality)
Mock1_df$Species <- as.factor(Mock1_df$Species)
set.seed(116)
M1_F1_rf <- randomForest(F1_score ~ Classifier + Database + Species + MinimumSeqQuality, # RF models inherently test interactions/nesting
                         data = Mock1_df,
                         ntree = 500, # Number of trees. More = more robust.
                         mtry = sqrt(4), # Number of variables tried at each split. Default = number of predictors / 3 for regression, but sqrt(p) is common in classification.
                         importance = TRUE) # Compute variable importance.
print(M1_F1_rf) #%VarExp = 73.02%
randomForest::importance(M1_F1_rf)
M1_vip <- vip(M1_F1_rf, title = "Mock1 F1 score Variable Importance") # Visualise the importance of each term in the RF model

# Database is (expectedly) the most important variable in explaining F1 score, followed by Species, then Classifier. Use the Database-split data to eliminate the effect of this term. 
set.seed(116)
M1_F1_GSD_rf <- randomForest(F1_score ~ Classifier + Species + MinimumSeqQuality, 
                             data = Mock1_df_GSD, 
                             ntree = 500, 
                             mtry = sqrt(3), # Different number of predictors, so different p value
                             importance = TRUE)
print(M1_F1_GSD_rf) # 49.3% variation explained
randomForest::importance(M1_F1_GSD_rf)
varImpPlot(M1_F1_GSD_rf, main = "Variable Importance, Mock1, Gold Standard Database")

set.seed(116)
M1_F1_NCBI_rf <- randomForest(F1_score ~ Classifier + Species + MinimumSeqQuality, 
                              data = Mock1_df_NCBI, 
                              ntree = 500, 
                              mtry = sqrt(3), 
                              importance = TRUE)
print(M1_F1_NCBI_rf) # % of variation explained is 55.01%
randomForest::importance(M1_F1_NCBI_rf)
varImpPlot(M1_F1_NCBI_rf, main = "Variable Importance, Mock1, NCBI Database")

# Repeat for Mock2 data
set.seed(116)
M2_F1_rf <- randomForest(F1_score ~ Classifier + Database + Species + MinimumSeqQuality, 
                         data = Mock2_df, 
                         ntree = 500, 
                         mtry = sqrt(4), 
                         importance = TRUE)
print(M2_F1_rf) # % of variation explained = 72.53%
# Evaluate the model with out-of-bag (OOB) prediction
M2_F1_rf_oob <- predict(M2_F1_rf, Mock1_df)
cor(Mock2_df$F1_score, M2_F1_rf_oob)^2 # Rough estimate of R^2 = 0.79
plot(Mock2_df$F1_score, M2_F1_rf_oob, 
     xlab = "Actual F1 Score", ylab = "Predicted F1 Score", 
     main = "OOB Predictions vs Actual")
abline(0, 1, col="red")
# R^2 is high, but the data doesn't seem to have a great correlation to the OOB predictions. This is likely a quirk of non-linear method. 

randomForest::importance(M2_F1_rf)
M2_vip <- varImpPlot(M2_F1_rf, main = "Variable Importance") # Visualise the importance of each term in the RF model
# Interestingly, here, Species and Classifier seem to have roughly the same influence on F1 score. Split the data by Database as above and re-model.

Mock2_df_GSD <- filter(Mock2_df, Database == "GSD")
Mock2_df_NCBI <- filter(Mock2_df, Database == "NCBI" | Database == "simNCBI")

set.seed(116)
M2_F1_GSD_rf <- randomForest(F1_score ~ Classifier + Species + MinimumSeqQuality, 
                             data = Mock2_df_GSD, 
                             ntree = 500, 
                             mtry = sqrt(3), # Different number of predictors, so different p value
                             importance = TRUE)
print(M2_F1_GSD_rf) # % of variation is lower at 53.72%
randomForest::importance(M2_F1_GSD_rf)
varImpPlot(M2_F1_GSD_rf, main = "Variable Importance, Mock2, Gold Standard Database") # Unlike when Database is a factor, Species far outweighs Classifier in terms of importance.

set.seed(116)
M2_F1_NCBI_rf <- randomForest(F1_score ~ Classifier + Species + MinimumSeqQuality, 
                              data = Mock2_df_NCBI, 
                              ntree = 500, 
                              mtry = sqrt(3), 
                              importance = TRUE)
print(M2_F1_NCBI_rf) # % of variation explained is 54.91%
randomForest::importance(M2_F1_NCBI_rf)
varImpPlot(M2_F1_NCBI_rf, main = "Variable Importance, Mock2, NCBI Database")
# In the most realistic situation, which is NCBI database + uneven community, Classifier is the most important variable, while species does contribute some variation. This model does only explain ~55% of overall variance.

# Realistically, we know that Classifier, Species and Database interact with each other. Test the interaction significance using the 
# The vivid package can be used to look at variable importance, and variable interaction simultaneously.
install.packages("vivid")
library(vivid)
vivi_M1_F1_rf <- vivi(fit = M1_F1_rf, data = Mock1_df, response = "F1_score")
viviHeatmap(vivi_M1_F1_rf)

# Plotting F1 score in boxplot with jitter
library(FSA)
library(rcompanion)
KW_M2GSD_F1 <- dunnTest(x = Mock2_df_GSD$F1_score, g = Mock2_df_GSD$Classifier, method = "sidak")
KW_M2GSD_F1
KW_M2GSD_F1.res <- KW_M2GSD_F1$res
Dunncld <- cldList(comparison = KW_M2GSD_F1.res$Comparison, p.value = KW_M2GSD_F1.res$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "Classifier"
GSD_table_M2 <- group_by(Mock2_df_GSD, Classifier) %>% summarise(mean=mean(F1_score, na.rm = TRUE), third_quantile = quantile(F1_score, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
GSD_table_M2 <- left_join(GSD_table_M2, Dunncld, by = "Classifier")
GSD_table_M2$Database <- "GSD"
print(GSD_table_M2)

KW_M2NCBI_F1 <- dunnTest(x = Mock2_df_NCBI$F1_score, g = Mock2_df_NCBI$Classifier, method = "sidak")
KW_M2NCBI_F1
KW_M2NCBI_F1.res <- KW_M2NCBI_F1$res
Dunncld <- cldList(comparison = KW_M2NCBI_F1.res$Comparison, p.value = KW_M2NCBI_F1.res$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "Classifier"
NCBI_table_M2 <- group_by(Mock2_df_NCBI, Classifier) %>% summarise(mean=mean(F1_score, na.rm = TRUE), third_quantile = quantile(F1_score, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
NCBI_table_M2 <- left_join(NCBI_table_M2, Dunncld, by = "Classifier")
NCBI_table_M2$Database <- "NCBI"
print(NCBI_table_M2)

grouped_table_M2 <- rbind(GSD_table_M2, NCBI_table_M2)

F1_dotplot_CLD <- ggplot(Mock2_df, aes(x = factor(Classifier), y = F1_score)) +
  geom_boxplot(aes(fill = Classifier), show.legend = FALSE, outlier.shape = NA) +
  geom_jitter(color = "black", size = 0.5, width = 0.2, alpha = 0.6) +
  facet_grid(. ~ Database) + 
  geom_text(data = grouped_table_M2, aes(x = Classifier, y = third_quantile, label = Letter), size = 3, vjust = -1, hjust = -1, color = "black") +
  labs(title = "F1 Score of Tested Classifiers", subtitle = "Qmin15", x = "Classifier", y = "F1 Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_brewer(palette = "Dark2")
F1_dotplot_CLD

# Because the F1 score for ML classifiers using the simNCBI database is so low, this is likely to skew the importance of Classifier in the RF. 
# Removing the three ML classifiers from M1 and M2 datasets and re-running the RF.
Mock1_df_NCBI_noML <- filter(Mock1_df_noML, Database == "NCBI")
Mock2_df_NCBI_noML <- filter(Mock2_df_noML, Database == "NCBI")

set.seed(116)
M1_F1_NCBI_rf_noML <- randomForest(F1_score ~ Classifier + Species + MinimumSeqQuality, 
                                   data = Mock1_df_NCBI_noML, 
                                   ntree = 500, 
                                   mtry = sqrt(3), 
                                   importance = TRUE)
print(M1_F1_NCBI_rf_noML) # % of variation explained is 42.15%
randomForest::importance(M1_F1_NCBI_rf_noML)
varImpPlot(M1_F1_NCBI_rf_noML, main = "Variable Importance, Mock1, NCBI Database")

set.seed(116)
M2_F1_NCBI_rf_noML <- randomForest(F1_score ~ Classifier + Species + MinimumSeqQuality, 
                                   data = Mock2_df_NCBI_noML, 
                                   ntree = 500, 
                                   mtry = sqrt(3), 
                                   importance = TRUE)
print(M2_F1_NCBI_rf_noML) # % of variation explained is 40.82%
randomForest::importance(M2_F1_NCBI_rf_noML)
varImpPlot(M2_F1_NCBI_rf_noML, main = "Variable Importance, Mock2, NCBI Database")
