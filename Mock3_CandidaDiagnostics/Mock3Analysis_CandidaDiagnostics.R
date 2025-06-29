library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(vegan)
library(tidyr)
library(stringr)

### ---------- Mock 3 Data Analysis ----------
file_list <- list.files(path = "./Input_ConfusionMatrixStats/", pattern = ".csv", full.names = TRUE)
confmat <- list()
for (i in seq_along(file_list)) {
  file <- file_list[[i]]
  df <- read.csv(file, header = TRUE)
  df_name <- tools::file_path_sans_ext(basename(file))
  confmat[[df_name]] <- df
}

# Extract only the row results for Candida species
subset_rows <- function(df) {
  subset(df, grepl("Candida", df$MockSpecies))
}
confmat_Candida <- lapply(confmat, subset_rows)

# Bind all data together
for (i in seq_along(confmat_Candida)) {
  confmat_Candida[[i]]$Classifier <- str_replace(confmat_Candida[[i]]$Classifier, "vuthuyduong", "vtd")
  confmat_Candida[[i]]$Classifier <- str_replace(confmat_Candida[[i]]$Classifier, "mycoAI_bert", "mycoAIBert")
  confmat_Candida[[i]]$Classifier <- str_replace(confmat_Candida[[i]]$Classifier, "mycoAI_cnn", "mycoAICNN")
}
all <- do.call(rbind, confmat_Candida)
write.csv(all, "./CompleteSummary_CandidaOnly.csv", quote = FALSE)

all_15 <- subset(all, all$MinimumSeqQuality == "Qmin15")

# Set parameters for plotting 
library(ggplot2)

classifier_colours <- c("minimap2" = "#1B9E77", 
                        "kraken2" = "#D95F02", 
                        "aodp" = "#7570B3", 
                        "dnabarcoder" = "#E7298A", 
                        "emu" = "#66A61E", 
                        "vtd" = "#E6AB02", 
                        "mycoAIBert" = "#A6761D", 
                        "mycoAICNN" = "#666666")
classifier_order <- c("aodp", "minimap2", "kraken2", "dnabarcoder", "emu", "vtd", "mycoAICNN", "mycoAIBert")
all_15$Classifier <- factor(all_15$Classifier, levels = classifier_order)

# Determine the distribution of the data
F1_shap_full <- shapiro.test(all$F1)
F1_Bart_full <- bartlett.test(all$F1, all$Classifier)

F1_PERMANOVA <- adonis2(all$F1_Score ~ all$Database/all$Classifier + all$MockSpecies + all$MinimumSeqQuality, by = "margin", method = "manhattan") # The / operator indicates that Classifier is nested within Database

# Test for the influence of classifier without database effects
F1_PERMANOVA_byDB <- list()
for (i in seq_along(confmat_Candida)) {
  F1_PERMANOVA_byDB[[i]] <- adonis2(confmat_Candida[[i]]$F1_Score ~ confmat_Candida[[i]]$Classifier + confmat_Candida[[i]]$MockSpecies + confmat_Candida[[i]]$MinimumSeqQuality, by = "margin", method = "manhattan")
}
names(F1_PERMANOVA_byDB) <- names(confmat_Candida)
print(F1_PERMANOVA_byDB)

# Use a Dunn's test for pairwise comparisons between groups
library(FSA)
confmat_Candida15 <- list()
for (i in seq_along(confmat_Candida)) {
  confmat_Candida15[[i]] <- subset(confmat_Candida[[i]], confmat_Candida[[i]]$MinimumSeqQuality == "Qmin15")
}
names(confmat_Candida15) <- names(confmat_Candida)

F1_Dunn_list <- list()
for (i in seq_along(confmat_Candida15)) {
  F1_Dunn_list[[i]] <- dunnTest(confmat_Candida15[[i]]$F1_Score, confmat_Candida15[[i]]$Classifier, method = "sidak")
}
names(F1_Dunn_list) <- names(confmat_Candida15)

# Extract CLD by database
library(rcompanion)
F1_Dunn_ResGSD <- F1_Dunn_list$`20250324_Mock3_GSD_ConfMatStats`$res
Dunncld_GSD <- cldList(comparison = F1_Dunn_ResGSD$Comparison, p.value = F1_Dunn_ResGSD$P.adj, threshold = 0.05) [1:2]
names(Dunncld_GSD)[1] <- "Classifier"
grouped_table_GSD <- group_by(confmat_Candida$`20250324_Mock3_GSD_ConfMatStats`, Classifier) %>% summarise(mean=mean(F1_Score, na.rm = TRUE), third_quantile = quantile(F1_Score, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_GSD <- left_join(grouped_table_GSD, Dunncld_GSD, by = "Classifier")

F1_Dunn_ResNCBI <- F1_Dunn_list$`20250324_Mock3_NCBI_ConfMatStats`$res
Dunncld_NCBI <- cldList(comparison = F1_Dunn_ResNCBI$Comparison, p.value = F1_Dunn_ResNCBI$P.adj, threshold = 0.05) [1:2]
names(Dunncld_NCBI)[1] <- "Classifier"
grouped_table_NCBI <- group_by(confmat_Candida$`20250324_Mock3_NCBI_ConfMatStats`, Classifier) %>% summarise(mean = mean(F1_Score, na.rm = TRUE), third_quantile = quantile(F1_Score, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_NCBI <- left_join(grouped_table_NCBI, Dunncld_NCBI, by = "Classifier")

F1_Dunn_ResNCBINoCd <- F1_Dunn_list$`20250324_Mock3_NCBINoCd_ConfMatStats`$res
Dunncld_NCBINoCd <- cldList(comparison = F1_Dunn_ResNCBINoCd$Comparison, p.value = F1_Dunn_ResNCBINoCd$P.adj, threshold = 0.05) [1:2]
names(Dunncld_NCBINoCd)[1] <- "Classifier"
grouped_table_NCBINoCd <- group_by(confmat_Candida$`20250324_Mock3_NCBINoCd_ConfMatStats`, Classifier) %>% summarise(mean=mean(F1_Score, na.rm = TRUE), third_quantile = quantile(F1_Score, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_NCBINoCd <- left_join(grouped_table_NCBINoCd, Dunncld_NCBINoCd, by = "Classifier")

F1_Dunn_RessimNCBI <- F1_Dunn_list$`20250324_Mock3_simNCBI_ConfMatStats`$res
Dunncld_simNCBI <- cldList(comparison = F1_Dunn_RessimNCBI$Comparison, p.value = F1_Dunn_RessimNCBI$P.adj, threshold = 0.05) [1:2]
names(Dunncld_simNCBI)[1] <- "Classifier"
grouped_table_simNCBI <- group_by(confmat_Candida$`20250324_Mock3_simNCBI_ConfMatStats`, Classifier) %>% summarise(mean=mean(F1_Score, na.rm = TRUE), third_quantile = quantile(F1_Score, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_simNCBI <- left_join(grouped_table_simNCBI, Dunncld_simNCBI, by = "Classifier")

# Generate plots and add CLD from Dunn's tests
confmat_Candida15$`20250324_Mock3_GSD_ConfMatStats`$Classifier <- factor(confmat_Candida15$`20250324_Mock3_GSD_ConfMatStats`$Classifier, levels = classifier_order)
boxplotCLD_GSD <- ggplot(data = confmat_Candida15$`20250324_Mock3_GSD_ConfMatStats`, mapping = aes(x = Classifier, y = F1_Score, fill = Classifier)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  # Boxplot without showing outliers separately
  geom_jitter(width = 0.2, height = 0, alpha = 0.7, size = 1) +  # Jitter points for each individual species
  labs(subtitle = "Gold Standard Database", x = "Classifier", y = "F1 score") + 
  geom_text(data = grouped_table_GSD, aes(x = Classifier, y = third_quantile, label = Letter), size = 5, vjust = -1, hjust = -1) +
  scale_fill_manual(values = classifier_colours) +
  theme_minimal() +
  scale_x_discrete(labels = c("aodp" = "aodp", 
                              "minimap2" = "minimap2", 
                              "kraken2" = "Kraken2", 
                              "dnabarcoder" = "Dnabarcoder", 
                              "emu" = "Emu", 
                              "vtd" = "vtd", 
                              "mycoAICNN" = "mycoAI CNN", 
                              "mycoAIBert" = "mycoAI BERT")) +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "lightgrey")) +
  ylim(0,1.05)
boxplotCLD_GSD
ggsave("./F1BoxplotCLD_GSD.png", boxplotCLD_GSD, bg = "white", height = 1400, width = 2200, units = "px")

confmat_Candida15$`20250324_Mock3_NCBI_ConfMatStats`$Classifier <- factor(confmat_Candida15$`20250324_Mock3_NCBI_ConfMatStats`$Classifier, levels = classifier_order)
boxplotCLD_NCBI <- ggplot(data = confmat_Candida15$`20250324_Mock3_NCBI_ConfMatStats`, mapping = aes(x = Classifier, y = F1_Score, fill = Classifier)) + 
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.7, size = 1) +
  labs(subtitle = "NCBI RefSeq ITS Database", x = "Classifier", y = "F1 score") +
  geom_text(data = grouped_table_NCBI, aes(x = Classifier, y = third_quantile, label = Letter), size = 5, vjust = -1, hjust = -1) + 
  scale_fill_manual(values = classifier_colours) +
  theme_minimal() + 
  scale_x_discrete(labels = c("aodp" = "aodp", 
                              "minimap2" = "minimap2", 
                              "kraken2" = "Kraken2", 
                              "dnabarcoder" = "Dnabarcoder", 
                              "emu" = "Emu", 
                              "vtd" = "vtd", 
                              "mycoAICNN" = "mycoAI CNN", 
                              "mycoAIBert" = "mycoAI BERT")) +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
        axis.title.y = element_text(size = 15),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "lightgrey")) +
  ylim(0,1.05)
boxplotCLD_NCBI
ggsave("./F1BoxplotCLD_NCBI.png", boxplotCLD_NCBI, bg = "white", height = 1300, width = 2150, units = "px")

confmat_Candida15$`20250324_Mock3_NCBINoCd_ConfMatStats`$Classifier <- factor(confmat_Candida15$`20250324_Mock3_NCBINoCd_ConfMatStats`$Classifier, levels = classifier_order)
boxplotCLD_NCBINoCd <- ggplot(data = confmat_Candida15$`20250324_Mock3_NCBINoCd_ConfMatStats`, mapping = aes(x = Classifier, y = F1_Score, fill = Classifier)) + 
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.7, size = 1) +
  labs(subtitle = expression("NCBI RefSeq ITS Database, no " * italic("C. dubliniensis")), x = "Classifier", y = "F1 score") +
  geom_text(data = grouped_table_NCBINoCd, aes(x = Classifier, y = third_quantile, label = Letter), size = 5, vjust = -1, hjust = -1) +
  scale_fill_manual(values = classifier_colours) +
  theme_minimal() +
  scale_x_discrete(labels = c("aodp" = "aodp", 
                              "minimap2" = "minimap2", 
                              "kraken2" = "Kraken2", 
                              "dnabarcoder" = "Dnabarcoder", 
                              "emu" = "Emu", 
                              "vtd" = "vtd", 
                              "mycoAICNN" = "mycoAI CNN", 
                              "mycoAIBert" = "mycoAI BERT")) +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "lightgrey")) +
  ylim(0,1.05)
boxplotCLD_NCBINoCd
ggsave("./F1BoxplotCLD_NCBINoCd.png", boxplotCLD_NCBINoCd, bg = "white", height = 1250, width = 2100, units = "px")

confmat_Candida15$`20250324_Mock3_simNCBI_ConfMatStats`$Classifier <- factor(confmat_Candida15$`20250324_Mock3_simNCBI_ConfMatStats`$Classifier, levels = classifier_order)
boxplotCLD_simNCBI <- ggplot(data = confmat_Candida15$`20250324_Mock3_simNCBI_ConfMatStats`, mapping = aes(x = Classifier, y = F1_Score, fill = Classifier)) + 
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.05, height = 0, alpha = 0.7, size = 1) + 
  labs(subtitle = "Simulated NCBI RefSeq ITS Database", x = "Classifier", y = "F1 score") +
  geom_text(data = grouped_table_simNCBI, aes(x = Classifier, y = third_quantile, label = Letter), size = 5, vjust = -1, hjust = -1) +
  scale_fill_manual(values = classifier_colours) +
  theme_minimal() + 
  scale_x_discrete(labels = c("aodp" = "aodp", 
                              "minimap2" = "minimap2", 
                              "kraken2" = "Kraken2", 
                              "dnabarcoder" = "Dnabarcoder", 
                              "emu" = "Emu", 
                              "vtd" = "vtd", 
                              "mycoAICNN" = "mycoAI CNN", 
                              "mycoAIBert" = "mycoAI BERT")) +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        plot.subtitle = element_text(hjust = 0.5, size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"), 
        axis.title.y = element_text(size = 13),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(color = "lightgrey")) +
  ylim(0, 1.05)
boxplotCLD_simNCBI
ggsave("./F1BoxplotCLD_simNCBI.png", boxplotCLD_simNCBI, bg = "white", height = 1080, width = 1950, units = "px")















































