library(readxl)
library(ggplot2)
library(dplyr)
library(FSA)
library(rcompanion)
library(AER)
library(MASS)
library(ordinal)
library(MuMIn)
library(nlme)
library(mgcv)

# Read in the data. Round the sequences_classified column because the Emu data isn't an integer. 
Mock1 <- read.csv("./Mock1_NotInNCBI.csv", header = TRUE, sep = ",")
Mock1$sequences_classified <- round(Mock1$sequences_classified)
Mock2 <- read.csv("./Mock2_NotInNCBI.csv", header = TRUE, sep = ",")
Mock2$sequences_classified <- round(Mock2$sequences_classified)
ML <- read.csv('./MLModels_NotInNCBI.csv', header = TRUE, sep = ",")
ML$sequences_classified <- round(ML$sequences_classified)
ML_M1 <- filter(ML, mock_replicate == "Mock1")
ML_M2 <- filter(ML, mock_replicate == "Mock2")
Mock1_all <- rbind(Mock1, ML_M1)
Mock2_all <- rbind(Mock2, ML_M2)

# Because the closeness score is an ordinal variable, by definition it cannot be parametric. 
# Use CLMMs to model closeness score variance for alignment-based (including Kraken2) classifiers.
# Make test variables factors
Mock1$classifier <- as.factor(Mock1$classifier)
Mock2$classifier <- as.factor(Mock2$classifier)
Mock1$mock_species <- as.factor(Mock1$mock_species)
Mock2$mock_species <- as.factor(Mock2$mock_species)
Mock1$quality_filt <- as.factor(Mock1$quality_filt)
Mock2$quality_filt <- as.factor(Mock2$quality_filt)

# Test whether including minimum quality score improves model fitness using full Mock1/Mock2 data. 
Mock1$sequences_classified <- scale(Mock1$sequences_classified)
Mock1_full_CLMM <- clmm(as.factor(closeness_score) ~ sequences_classified + factor(classifier) + (1|mock_species) + factor(quality_filt), data = Mock1)
Mock2$sequences_classified <- scale(Mock2$sequences_classified)
Mock2_full_CLMM <- clmm(as.factor(closeness_score) ~ sequences_classified + factor(classifier) + (1|mock_species) + factor(quality_filt), data = Mock2)
summary(Mock1_full_CLMM) #AIC = 374.64, random effect variance = 7.07, quality_filt not significant, sequences_classified slightly significant. 
summary(Mock2_full_CLMM) #AIC = 382, random effect variance = 5.421, quality_filt not significant or sequences_classified. 

Mock1_full_CLMM_noS <- clmm(as.factor(closeness_score) ~ sequences_classified + factor(classifier) + (1|mock_species), data = Mock1)
Mock2_full_CLMM_noS <- clmm(as.factor(closeness_score) ~ sequences_classified + factor(classifier) + (1|mock_species), data = Mock2)
summary(Mock1_full_CLMM_noS) #AIC = 372.78
summary(Mock2_full_CLMM_noS) #AIC = 380.01
# AIC values are lower without quality_filt considered in the model. Can filter the data to reduce model complexity.

# Separate the data based on minimum Q-score.
Mock1_15 <- Mock1 %>% filter(quality_filt == "Qmin15")
Mock2_15 <- Mock2 %>% filter(quality_filt == "Qmin15")

# Because the negative binomial model indicated that classifier may influence the number of sequences classified, may need to include sequences classified as a predictor of closeness score, in case it has a latent effect. First, check whether there is a significant correlation between sequences classified and classifier, because this will influence the behavior of the model. 
KW <- kruskal.test(sequences_classified ~ classifier, data = Mock1_15)
KW # p = 0.033
Dunn <- dunnTest(Mock1_15$sequences_classified, Mock1_15$classifier, method = "sidak")
Dunn
# These tests show that the number of sequences classified does change with level of classifier, and particularly that minimap2 has a different effect on the number of sequences classified in comparison to kraken2. 
# To better account for this in the model, classifier will be a factor, as well as being a fixed effect. This means the model will test each level of classifier separately, rather than assuming a continuous, linear relationship between classifier and closeness score. 
# Also scale the sequences classified column, because it has such a different scale to the other numeric variables. Species is included as a random effect. Because the negative binomial model indicated no significant effect of quality filtering on the closeness score, model the Qmin15 data alone. 
Mock1_15$sequences_classified <- scale(Mock1_15$sequences_classified)
Mock1_CLMM <- clmm(as.factor(closeness_score) ~ sequences_classified + factor(classifier) + (1|mock_species), data = Mock1_15)
summary (Mock1_CLMM) # AIC = 215.87
# Dnabarcoder was interpreted as the 'reference' level of classifier. The only significant result was that minimap2 would have a positive influence on closeness score, increasing it's value. Minimap2 is therefore likely to classify unknown sequences further away from their true taxonomy. 
# Emu and Kraken2 were not significantly different from Dnabarcoder in their influence on closeness score. Sequences classified did not have a significant effect on closeness score. 
# The random effect for mock_species has a variance of 5.38 and standard deviation of 2.32, which indicates that there is some variability in closeness score across different species. To investigate this further, generate the same model but without the random effect using ordered logistic regression, and compare AIC values to see whether adding the random effect improves model fit: 
Mock1_OLR <- polr(as.factor(closeness_score) ~ sequences_classified + factor(classifier), data = Mock1_15)
AIC(Mock1_CLMM) # = 215.87
AIC(Mock1_OLR) # = 233.94
# AIC is the Akaike Information Criterion. Lower values indicate a better-fitting model. While the difference is only small, including species as a random effect does influence the fitness of the model.

KW2 <- kruskal.test(sequences_classified ~ classifier, data = Mock2_15)
KW2 # p = 0.267
Dunn2 <- dunnTest(Mock2_15$sequences_classified, Mock2_15$classifier, method = "sidak")
Dunn2
# The Kruskal-Wallis test indicated that there is not a significant influence of classifier on the number of sequences classified for Mock 2, there was no significant difference between some classifiers during pairwise comparisons. 
Mock2_15$sequences_classified <- scale(Mock2_15$sequences_classified)
Mock2_CLMM <- clmm(as.factor(closeness_score) ~ sequences_classified + factor(classifier) + (1|mock_species), data = Mock2_15)
summary(Mock2_CLMM) # AIC = 217.7. Random effect variance = 4.62, standard deviation = 2.15. 
Mock2_CLMM3 <- clmm(as.factor(closeness_score) ~ factor(classifier) + (1|mock_species), data = Mock2_15)
summary(Mock2_CLMM3) # AIC = 215.79. Random effect variance = 4.22, standard deviation = 2.06
# Again, minimap2 is the only classifier which differently impacts closeness score in comparison to the reference (dnabarcoder). It increases the value of the closeness score, therefore minimap2 places unknown taxa further away from their true taxonomy. 
# As expected from KW test results, removing number of sequences classified from the model improves fitness.
# Because the random effect has a variance and stdev > 1, again test whether a model without the random effect would fit better: 
Mock2_OLR <- polr(as.factor(closeness_score) ~ factor(classifier), data = Mock2_15)
AIC(Mock2_CLMM) #217.7
AIC(Mock2_CLMM3) #215.79
AIC(Mock2_OLR)#233.51
# As smaller AIC values indicate better fit, including mock_species as a random effect has a small, but positive contribution to model fitness. 

# Because the AIC values in general are quite large, testing whether a non-linear model would fit better, starting with Mock 1 data. 
# Try the NLME:
Mock1_NLME <- nlme(closeness_score ~ a * exp(b * sequences_classified), data = Mock1_15, fixed = a + b ~ 1, random = a + b ~ 1 | mock_species, start = c(a = 1, b = 0.1))
summary(Mock1_NLME)
# AIC value for this model is 284.5, which is worse than the linear model. 

# Use Dunn's test (pairwise Kruskal-Wallis) to test pairwise differences between classifiers.
# Use delta closeness score to enable comparisons between ML and non-ML classifiers.
Mock1_KW <- dunnTest(delta_closeness_score ~ classifier, data = Mock1, method = "sidak", kw = TRUE)
Mock2_KW <- dunnTest(delta_closeness_score ~ classifier, data = Mock2, method = "sidak", kw = TRUE)

# Assess whether there's a significant difference between Q score levels for delta closeness score
M1_Qscore_KW <- kruskal.test(delta_closeness_score ~ factor(quality_filt), data = Mock1)
M1_Qscore_KW # p = 0.82
M2_Qscore_KW <- kruskal.test(delta_closeness_score ~ factor(quality_filt), data = Mock2)
M2_Qscore_KW # p = 0.75

# Filter to only Qmin15
M1_15 <- filter(Mock1_all, quality_filt == "Qmin15")
M2_15 <- filter(Mock2_all, quality_filt == "Qmin15")

# Use KW test to see whether classifier influences closeness score. Use Dunn's test to correct for multiple comparisons.
M1_KW <- kruskal.test(formula = delta_closeness_score ~ factor(classifier), data = M1_15)
M1_D <- dunnTest(x = M1_15$delta_closeness_score, g = M1_15$classifier, method = "sidak")
M1_D
M1_Dres <- M1_D$res # Extract the result of the Dunn test to a data frame
Dunncld <- cldList(comparison = M1_Dres$Comparison, p.value = M1_Dres$P.adj, threshold = 0.05)[1:2] # Generate CLD
names(Dunncld)[1] <- "classifier" # Change the name of the first column to reflect the grouping variable, consistent with the column names in the original data frame. 
grouped_table_M1.15 <- group_by(M1_15, classifier) %>% summarise(mean=mean(delta_closeness_score, na.rm = TRUE), third_quantile=quantile(delta_closeness_score, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_M1.15 <- left_join(grouped_table_M1.15, Dunncld, by = "classifier")

M1_15$classifier <- factor(M1_15$classifier, levels = c("minimap2", "Kraken2", "Dnabarcoder", "Emu", "vtd", "mycoAICNN", "mycoAIBERT"))
M2_15$classifier <- factor(M2_15$classifier, levels = c("minimap2", "Kraken2", "Dnabarcoder", "Emu", "vtd", "mycoAICNN", "mycoAIBERT"))

# Open both the 'grouped_table_NP' and 'Dunncld' objects and double-check that the CLD letters have matched properly during merging
closeness_boxplot_M1.15 <- ggplot(M1_15, aes(classifier, delta_closeness_score)) + 
  geom_boxplot(aes(fill = classifier), show.legend = FALSE, outlier.shape = NA) +
  geom_jitter(color = "black", size = 1, width = 0.2, height = 0, alpha = 0.6) +
  labs(title = "Classification of 'Unknown' taxa", subtitle = "Mock1, Qmin15", x = "", y = "Delta Closeness Score") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red") +
  theme_minimal() +
  geom_text(data = grouped_table_M1.15, aes(x = classifier, y = third_quantile, label = Letter), size = 5, vjust = -1, hjust = -1) +
  scale_fill_manual(values = c("Dnabarcoder" = "#E7298A", 
                               "Emu" = "#66A61E", 
                               "Kraken2" = "#D95F02", 
                               "minimap2" = "#1B9E77", 
                               "vtd" = "#E6AB02", 
                               "mycoAIBERT" = "#A6761D", 
                               "mycoAICNN" = "#666666")) +
  theme(axis.text = element_text(size = 14, color = "black"), 
        plot.title = element_text(size = 16), 
        plot.subtitle = element_text(size = 14), 
        axis.title.y = element_text(size = 15),
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.background = element_rect(fill = "white")) +
  scale_x_discrete(labels = c("minimap2" = "minimap2", 
                              "Kraken2" = "Kraken2", 
                              "Dnabarcoder" = "Dnabarcoder",
                              "Emu" = "Emu", 
                              "vtd" = "vtd", 
                              "mycoAICNN" = "mycoAI CNN", 
                              "mycoAIBERT" = "mycoAI BERT"))
closeness_boxplot_M1.15
ggsave("./ClosenessScoreBoxplot_M1.15.png", plot = closeness_boxplot_M1.15, dpi = 1000, bg = "white")

# Get mean values and check
M1_15_means <- aggregate(M1_15, list(M1_15$classifier), mean)

M2_KW <- kruskal.test(formula = delta_closeness_score ~ factor(classifier), data = M2_15)
M2_D <- dunnTest(x = M2_15$delta_closeness_score, g = M2_15$classifier, method = "sidak")
M2_D
M2_Dres <- M2_D$res
Dunncld <- cldList(comparison = M2_Dres$Comparison, p.value = M2_Dres$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "classifier"
grouped_table_M2.15 <- group_by(M2_15, classifier) %>% summarise(mean=mean(delta_closeness_score, na.rm = TRUE), third_quantile = quantile(delta_closeness_score, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_M2.15 <- left_join(grouped_table_M2.15, Dunncld, by = "classifier")

closeness_boxplot_M2.15 <- ggplot(M2_15, aes(classifier, delta_closeness_score)) + 
  geom_boxplot(aes(fill = classifier), show.legend = FALSE, outlier.shape = NA) + 
  geom_jitter(color = "black", size = 1, width = 0.2, height = 0, alpha = 0.6) +
  labs(title = "Classification of 'Unknown' taxa", subtitle = "Mock2, Qmin15", x = "", y = "Delta Closeness Score") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red") +
  geom_text(data = grouped_table_M2.15, aes(x = classifier, y = third_quantile, label = Letter), size = 5, vjust = -1, hjust = -1) +
  scale_fill_manual(values = c("Dnabarcoder" = "#E7298A", 
                               "Emu" = "#66A61E", 
                               "Kraken2" = "#D95F02", 
                               "minimap2" = "#1B9E77", 
                               "vtd" = "#E6AB02", 
                               "mycoAIBERT" = "#A6761D", 
                               "mycoAICNN" = "#666666")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 14, color = "black"), 
        plot.title = element_text(size = 16), 
        plot.subtitle = element_text(size = 14), 
        axis.title.y = element_text(size = 15),
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(fill = NA)) +
  scale_x_discrete(labels = c("minimap2" = "minimap2", 
                              "Kraken2" = "Kraken2", 
                              "Dnabarcoder" = "Dnabarcoder",
                              "Emu" = "Emu", 
                              "vtd" = "vtd", 
                              "mycoAICNN" = "mycoAI CNN", 
                              "mycoAIBERT" = "mycoAI BERT"))
closeness_boxplot_M2.15
ggsave("./ClosenessScoreBoxplot_M2.15.png", plot = closeness_boxplot_M2.15, dpi = 1000, bg = "white")

M2_15_means <- aggregate(M2_15, list(M2_15$classifier), mean)
