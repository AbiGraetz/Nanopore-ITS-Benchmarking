library(readxl)
library(tibble)
library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)
library(purrr)

# ---------- GSD data wrangling ---------- 
columnnames<- c("Species", "AssignedSpecies", "TotalCount")

GSD_Mock1 <- list()
files_GSDM1 <- list.files(path = './Mock1/GSD', pattern = ".csv", full.names = TRUE)
for (i in seq_along(files_GSDM1)) {
  file = files_GSDM1[[i]]
  df <- read.csv(file, header = FALSE)
  colnames(df) <- columnnames
  df_name <- tools::file_path_sans_ext(basename(file)) 
  GSD_Mock1[[df_name]] <- df
}
GSD_Mock2 <- list()
files_GSDM2 <- list.files(path = './Mock2/GSD', pattern = ".csv", full.names = TRUE)
for (i in seq_along(files_GSDM2)) {
  file = files_GSDM2[[i]]
  df <- read.csv(file, header = FALSE)
  colnames(df) <- columnnames
  df_name <- tools::file_path_sans_ext(basename(file))
  GSD_Mock2[[df_name]] <- df
}

# Separate out the AssignedSpecies column to get at species-level assignments
for (i in seq_along(GSD_Mock1)) {
  GSD_Mock1[[i]] <- separate_wider_delim(GSD_Mock1[[i]], AssignedSpecies, delim = ",", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "AssignedSpecies"), too_many = "merge", too_few = "align_start")
  GSD_Mock1[[i]] <- dplyr::select(GSD_Mock1[[i]], c("Species", "AssignedSpecies", "TotalCount"))
}
print(GSD_Mock1[[1]])

# Reformat results to be consistent with Kraken2
reformat_results <- function(df) {
  df %>% 
    group_by(Species) %>%
    summarise(
      CorrectlyMapped = sum(TotalCount[Species == AssignedSpecies]),
      .groups = "drop"
    ) %>%
    rename(MockSpecies = Species) %>%
    left_join(
      df %>%
        group_by(AssignedSpecies) %>%
        summarise(TotalMapped = sum(TotalCount), .groups = "drop"), 
      by = c("MockSpecies" = "AssignedSpecies")
    )
}
GSD_Mock1_reformat <- map(GSD_Mock1, reformat_results)

# Add in the Kraken2 data to the list, as everything should be in the same format now.
K2cols <- c("MockSpecies", "CorrectlyMapped", "TotalMapped")
files_K2M1 <- list.files(path = './Mock1/GSD/Kraken2/', pattern = ".csv", full.names = TRUE)
K2list <- list()
for (i in seq_along(files_K2M1)) {
  file = files_K2M1[[i]]
  df <- read.csv(file, header = FALSE)
  colnames(df) <- K2cols
  df_name <- tools::file_path_sans_ext(basename(file)) 
  K2list[[df_name]] <- df
}
GSD_Mock1_reformat <- c(GSD_Mock1_reformat, K2list)

# Sanity check
check_df <- filter(GSD_Mock1[[21]], AssignedSpecies == "Austropuccinia_psidii")
check_count <- sum(check_df$TotalCount)
reformat_count <- GSD_Mock1_reformat[[21]] %>% filter(MockSpecies == "Austropuccinia_psidii") %>% pull(TotalMapped)

if (check_count == reformat_count) {
  print("Counts are equal, sanity checked!")
}

# Convert to BIOM format
all_species <- unique(unlist(lapply(GSD_Mock1, function(df) df$Species)))
GSD_Mock1_OTUtable <- data.frame(Species = all_species)

for (i in seq_along(GSD_Mock1_reformat)) {
  df <- GSD_Mock1_reformat[[i]]
  temp <- data.frame(Species = all_species)
  temp$TotalCount <- df$TotalMapped[match(all_species, df$MockSpecies)]
  temp$TotalCount <- round(as.numeric(temp$TotalCount))
  GSD_Mock1_OTUtable[[names(GSD_Mock1_reformat)[i]]] <- temp$TotalCount
  GSD_Mock1_reformat[[i]] <- df
}

# Calculate whether there are any 'unclassified' sequences, given that the total number is known. For Mock 1, this is 400 * 54 sequences total in the Mock Community
total_sequences_M1 <- 54 * 400
unclassified_counts <- total_sequences_M1 - colSums(sapply(GSD_Mock1_OTUtable[ , -1], as.numeric))
GSD_Mock1_OTUtable <- rbind(GSD_Mock1_OTUtable, c("Unclassified", unclassified_counts))

# Read in the ground truth and append
Mock1_GT <- read.csv('./Mock1/Mock1Abundance_split.csv', header = TRUE)
Mock1_GT <- dplyr::select(Mock1_GT, c('organism', 'Mock1'))
colnames(Mock1_GT) <- c('Species', 'Mock1_GroundTruth_GroundTruth_Qmin15')
Mock1_GT$Species <- gsub(" ", "_", Mock1_GT$Species)
GSD_Mock1_OTUtable <- left_join(GSD_Mock1_OTUtable, Mock1_GT, by = "Species")
GSD_Mock1_OTUtable[is.na(GSD_Mock1_OTUtable)] <- 0

# Sanity check: check column sums are equal
# Ensure all columns except the first are numeric.
i <- c(2:ncol(GSD_Mock1_OTUtable))
GSD_Mock1_OTUtable[, i] <- apply(GSD_Mock1_OTUtable[, i], 2, function(x) as.numeric(x))
# Calculate sums of all columns
column_sums <- colSums(GSD_Mock1_OTUtable[, -1], na.rm = TRUE)
# Check if all sums are equal
all_equal <- length(unique(column_sums)) == 1
# Print results
cat("Column sums:\n")
print(column_sums)
cat("All columns sum to the same number:", all_equal, "\n")
# If not equal, show which columns differ
if (!all_equal) {
  cat("Unique sum values:", unique(column_sums), "\n")
}

# Calculate Bhattacharyya distance and Shannon and Simpson alpha diversity metrics
GSD_Mock1_OTUmatrix <- as.matrix(GSD_Mock1_OTUtable[, -1])
GSD_Mock1_OTUmatrix <- apply(GSD_Mock1_OTUmatrix, 2, as.numeric)
# Calculate the Bhattacharyya distance using the ground truth species distribution. 
ground_truth <- GSD_Mock1_OTUmatrix[, "Mock1_GroundTruth_GroundTruth_Qmin15"]
norm_matrix <- apply(GSD_Mock1_OTUmatrix, 2, function(x) x / sum(x))
norm_GT <- ground_truth / sum(ground_truth)
bhattacharyya_distance <- function(p, q) {
  -log(sum(sqrt(p * q)))
}
bhattacharyya <- apply(norm_matrix, 2, function(x) bhattacharyya_distance(x, norm_GT))
# Remove the ground truth column to calculate the alpha diversity metrics. Calculate this separately for the GT. 
GSD_Mock1_OTUmatrix <- GSD_Mock1_OTUmatrix[, colnames(GSD_Mock1_OTUmatrix) != "Mock1_GroundTruth_GroundTruth_Qmin15"]
shannon_GSDM1 <- apply(GSD_Mock1_OTUmatrix, 2, function(x) diversity(x, index = "shannon"))
simpson_GSDM1 <- apply(GSD_Mock1_OTUmatrix, 2, function(x) diversity(x, index = "simpson")) 
# Calculate the alpha diversity metrics for the ground truth alone, then append to the end of the appropriate vectors. Otherwise, the diversity table won't form, because the Bhattacharyya distance vector will be one longer than the other vectors (GT-GT comparison).
shannon_M1GT <- diversity(as.matrix(Mock1_GT[, -1]), index = "shannon")
simpson_M1GT <- diversity(as.matrix(Mock1_GT[, -1]), index = "simpson")
shannon_GSDM1 <- append(shannon_GSDM1, shannon_M1GT, after = length(shannon_GSDM1))
simpson_GSDM1 <- append(simpson_GSDM1, simpson_M1GT, after = length(simpson_GSDM1))
GSD_Mock1_diversity <- data.frame(
  Sample = c(colnames(GSD_Mock1_OTUmatrix), "Mock1_GroundTruth_GroundTruth_Qmin15"), 
  Shannon = shannon_GSDM1, 
  Simpson = simpson_GSDM1, 
  Bhattacharyya = bhattacharyya
)
GSD_Mock1_diversity <- separate_wider_delim(GSD_Mock1_diversity, Sample, "_", names = c("MockCommunity", "Classifier", "SubRep", "MinimumSeqQuality"), too_few = "align_start")
GSD_Mock1_diversity$Database <- "Gold Standard Database"

# Calculate delta Shannon and Simpson
GT_shannon <- GSD_Mock1_diversity$Shannon[[nrow(GSD_Mock1_diversity)]] # Use value from final row, because the ground truth is appended to the end.
GT_simpson <- GSD_Mock1_diversity$Simpson[[nrow(GSD_Mock1_diversity)]]
GSD_Mock1_diversity$DeltaShannon <- GSD_Mock1_diversity$Shannon - GT_shannon
GSD_Mock1_diversity$DeltaSimpson <- GSD_Mock1_diversity$Simpson - GT_simpson

# Use technical replicates to calculate means and 95% CIs
aodp15 <- filter(GSD_Mock1_diversity, Classifier == "aodp" & MinimumSeqQuality == "Qmin15")
aodp17 <- filter(GSD_Mock1_diversity, Classifier == "aodp" & MinimumSeqQuality == "Qmin17")
dnabarcoder15 <- filter(GSD_Mock1_diversity, Classifier == "dnabarcoder" & MinimumSeqQuality == "Qmin15")
dnabarcoder17 <- filter(GSD_Mock1_diversity, Classifier == "dnabarcoder" & MinimumSeqQuality == "Qmin17")
emu15 <- filter(GSD_Mock1_diversity, Classifier == "emu" & MinimumSeqQuality == "Qmin15")
emu17 <- filter(GSD_Mock1_diversity, Classifier == "emu" & MinimumSeqQuality == "Qmin17")
kraken15 <- filter(GSD_Mock1_diversity, Classifier == "kraken2" & MinimumSeqQuality == "Qmin15")
kraken17 <- filter(GSD_Mock1_diversity, Classifier == "kraken2" & MinimumSeqQuality == "Qmin17")
minimap15 <- filter(GSD_Mock1_diversity, Classifier == "minimap2" & MinimumSeqQuality == "Qmin15")
minimap17 <- filter(GSD_Mock1_diversity, Classifier == "minimap2" & MinimumSeqQuality == "Qmin17")
mycoAIbert15 <- filter(GSD_Mock1_diversity, Classifier == "mycoAIbert" & MinimumSeqQuality == "Qmin15")
mycoAIbert17 <- filter(GSD_Mock1_diversity, Classifier == "mycoAIbert" & MinimumSeqQuality == "Qmin17")
mycoAICNN15 <- filter(GSD_Mock1_diversity, Classifier == "mycoAICNN" & MinimumSeqQuality == "Qmin15")
mycoAICNN17 <- filter(GSD_Mock1_diversity, Classifier == "mycoAICNN" & MinimumSeqQuality == "Qmin17")
vtd15 <- filter(GSD_Mock1_diversity, Classifier == "vtd" & MinimumSeqQuality == "Qmin15")
vtd17 <- filter(GSD_Mock1_diversity, Classifier == "vtd" & MinimumSeqQuality == "Qmin17")

GSD_M1_splitlist <- list(aodp17, aodp15, dnabarcoder17, dnabarcoder15, emu15, emu17, kraken15, kraken17, minimap15, minimap17, mycoAIbert15, mycoAIbert17, mycoAICNN15, mycoAICNN17, vtd15, vtd17)

calculate_summary <- function(df, metric_cols, metadata_cols) {
  summary_list <- lapply(metric_cols, function(metric) {
    values <- df[[metric]]
    n <- length(values)
    mean_val <- mean(values, na.rm = TRUE)
    stderr <- sd(values, na.rm = TRUE) / sqrt(n)
    ci_lower <- mean_val - qt(0.975, df = n - 1) * stderr
    ci_upper <- mean_val + qt(0.975, df = n - 1) * stderr
    
    data.frame(
      Metric = metric,
      Mean = mean_val,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper,
      df[1, metadata_cols, drop = FALSE]
    )
  })
  do.call(rbind, summary_list)
}
metric_cols <- c("DeltaSimpson", "DeltaShannon", "Bhattacharyya")
metadata_cols <- c("Classifier", "Database", "MinimumSeqQuality")
summary_list <- lapply(GSD_M1_splitlist, calculate_summary, metric_cols = metric_cols, metadata_cols = metadata_cols)
GSD_M1_summary <- do.call(rbind, summary_list)
GSD_M1_summary <- dplyr::select(GSD_M1_summary, c("Classifier", "Database", "MinimumSeqQuality", "Metric", "Mean", "CI_Lower", "CI_Upper"))

# Repeat for Mock 2
for (i in seq_along(GSD_Mock2)) {
  GSD_Mock2[[i]] <- separate_wider_delim(GSD_Mock2[[i]], AssignedSpecies, delim = ",", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "AssignedSpecies"), too_many = "merge", too_few = "align_start")
  GSD_Mock2[[i]] <- dplyr::select(GSD_Mock2[[i]], c("Species", "AssignedSpecies", "TotalCount"))
}
GSD_Mock2_reformat <- map(GSD_Mock2, reformat_results)
print(GSD_Mock2_reformat[[1]])
files_K2M2 <- list.files(path = './Mock2/GSD/Kraken2/', pattern = ".csv", full.names = TRUE)
K2list <- list()
for (i in seq_along(files_K2M2)) {
  file = files_K2M2[[i]]
  df <- read.csv(file, header = FALSE)
  colnames(df) <- K2cols
  df_name <- tools::file_path_sans_ext(basename(file)) 
  K2list[[df_name]] <- df
}
GSD_Mock2_reformat <- c(GSD_Mock2_reformat, K2list)

# Sanity check
check_df <- filter(GSD_Mock2[[21]], AssignedSpecies == "Austropuccinia_psidii")
check_count <- sum(check_df$TotalCount)
reformat_count <- GSD_Mock2_reformat[[21]] %>% filter(MockSpecies == "Austropuccinia_psidii") %>% pull(TotalMapped)
if (check_count == reformat_count) {
  print("Counts are equal, sanity checked!") } else {
    print("Counts are not equal - something is wrong!")
  }

# Create the BIOM table
all_species <- unique(unlist(lapply(GSD_Mock2, function(df) df$Species))) 

GSD_Mock2_OTUtable <- data.frame(Species = all_species)

for (i in seq_along(GSD_Mock2_reformat)) {
  df <- GSD_Mock2_reformat[[i]]
  temp <- data.frame(Species = all_species)
  temp$TotalCount <- df$TotalMapped[match(all_species, df$MockSpecies)]
  temp$TotalCount <- round(as.numeric(temp$TotalCount))
  GSD_Mock2_OTUtable[[names(GSD_Mock2_reformat)[i]]] <- temp$TotalCount
  GSD_Mock2_reformat[[i]] <- df
}

# Calculate whether there are any 'unclassified' sequences, given that the total number is known. For Mock 2, this is 21,779 sequences total in the Mock Community.
total_sequences_M2 <- 21779
unclassified_counts <- total_sequences_M2 - colSums(GSD_Mock2_OTUtable[ , -1])
GSD_Mock2_OTUtable <- rbind(GSD_Mock2_OTUtable, c("Unclassified", unclassified_counts))
Mock2_GT <- read.delim('./Mock2/Mock2Abundance_part.tsv', header = TRUE, sep = ",")
Mock2_GT <- dplyr::select(Mock2_GT, c('organism', 'sequence_abundance'))
colnames(Mock2_GT) <- c('Species', 'Mock2_GroundTruth_GroundTruth_Qmin15')
Mock2_GT$Species <- gsub(" ", "_", Mock2_GT$Species)
print(Mock2_GT)
GSD_Mock2_OTUtable <- left_join(GSD_Mock2_OTUtable, Mock2_GT, by = "Species")
GSD_Mock2_OTUtable[is.na(GSD_Mock2_OTUtable)] <- 0

# Sanity check
# Ensure all columns except the first are numeric.
i <- c(2:ncol(GSD_Mock2_OTUtable))
GSD_Mock2_OTUtable[, i] <- apply(GSD_Mock2_OTUtable[, i], 2, function(x) as.numeric(x))
# Calculate sums of all columns
column_sums <- colSums(GSD_Mock2_OTUtable[, -1], na.rm = TRUE)
# Check if all sums are equal
all_equal <- length(unique(column_sums)) == 1
# Print results
cat("Column sums:\n")
print(column_sums)
cat("All columns sum to the same number:", all_equal, "\n")
# If not equal, show which columns differ
if (!all_equal) {
  cat("Unique sum values:", unique(column_sums), "\n")
}

# Calculate diversity metrics
GSD_Mock2_OTUmatrix <- as.matrix(GSD_Mock2_OTUtable[, -1])
GSD_Mock2_OTUmatrix <- apply(GSD_Mock2_OTUmatrix, 2, as.numeric)
# Calculate the Bhattacharyya distance with the ground truth appended. 
ground_truth <- GSD_Mock2_OTUmatrix[, "Mock2_GroundTruth_GroundTruth_Qmin15"]
norm_matrix <- apply(GSD_Mock2_OTUmatrix, 2, function(x) x / sum(x))
norm_GT <- ground_truth / sum(ground_truth)
bhattacharyya_GSDM2 <- apply(norm_matrix, 2, function(x) bhattacharyya_distance(x, norm_GT))
# Remove the ground truth and calculate the alpha diversity metrics.
GSD_Mock2_OTUmatrix <- GSD_Mock2_OTUmatrix[, colnames(GSD_Mock2_OTUmatrix) != "Mock2_GroundTruth_GroundTruth_Qmin15"]
shannon_GSDM2 <- apply(GSD_Mock2_OTUmatrix, 2, function(x) diversity(x, index = "shannon"))
simpson_GSDM2 <- apply(GSD_Mock2_OTUmatrix, 2, function(x) diversity(x, index = "simpson")) 
# Calculate the alpha diversity metrics for the ground truth alone, then append to the end of the appropriate vectors. Otherwise, the diversity table won't form, because the Bhattacharyya distance vector will be one longer than the other vectors (GT-GT comparison).
shannon_M2GT <- diversity(as.matrix(Mock2_GT[, -1]), index = "shannon")
simpson_M2GT <- diversity(as.matrix(Mock2_GT[, -1]), index = "simpson")
shannon_GSDM2 <- append(shannon_GSDM2, shannon_M2GT, after = length(shannon_GSDM1))
simpson_GSDM2 <- append(simpson_GSDM2, simpson_M2GT, after = length(simpson_GSDM1))
GSD_Mock2_diversity <- data.frame(
  Sample = c(colnames(GSD_Mock2_OTUmatrix), "Mock2_GroundTruth_GroundTruth_Qmin15"), 
  Shannon = shannon_GSDM2, 
  Simpson = simpson_GSDM2, 
  Bhattacharyya = bhattacharyya_GSDM2
)
rownames(GSD_Mock2_diversity) <- NULL
GSD_Mock2_diversity <- separate_wider_delim(GSD_Mock2_diversity, Sample, "_", names = c("MockCommunity", "Classifier", "SubRep", "MinimumSeqQuality"), too_few = "align_start")
GSD_Mock2_diversity$Database <- "Gold Standard Database"

# Calculate delta Shannon and Simpson
GT_shannon <- GSD_Mock2_diversity$Shannon[[nrow(GSD_Mock2_diversity)]] # Use final row value, as ground truth is appended to the end.
GT_simpson <- GSD_Mock2_diversity$Simpson[[nrow(GSD_Mock2_diversity)]]
GSD_Mock2_diversity$DeltaShannon <- GSD_Mock2_diversity$Shannon - GT_shannon
GSD_Mock2_diversity$DeltaSimpson <- GSD_Mock2_diversity$Simpson - GT_simpson

# Use technical replicates to calculate mean and 95% CIs
aodp15 <- filter(GSD_Mock2_diversity, Classifier == "aodp" & MinimumSeqQuality == "Qmin15")
aodp17 <- filter(GSD_Mock2_diversity, Classifier == "aodp" & MinimumSeqQuality == "Qmin17")
dnabarcoder15 <- filter(GSD_Mock2_diversity, Classifier == "dnabarcoder" & MinimumSeqQuality == "Qmin15")
dnabarcoder17 <- filter(GSD_Mock2_diversity, Classifier == "dnabarcoder" & MinimumSeqQuality == "Qmin17")
emu15 <- filter(GSD_Mock2_diversity, Classifier == "emu" & MinimumSeqQuality == "Qmin15")
emu17 <- filter(GSD_Mock2_diversity, Classifier == "emu" & MinimumSeqQuality == "Qmin17")
kraken15 <- filter(GSD_Mock2_diversity, Classifier == "kraken2" & MinimumSeqQuality == "Qmin15")
kraken17 <- filter(GSD_Mock2_diversity, Classifier == "kraken2" & MinimumSeqQuality == "Qmin17")
minimap15 <- filter(GSD_Mock2_diversity, Classifier == "minimap2" & MinimumSeqQuality == "Qmin15")
minimap17 <- filter(GSD_Mock2_diversity, Classifier == "minimap2" & MinimumSeqQuality == "Qmin17")
mycoAIbert15 <- filter(GSD_Mock2_diversity, Classifier == "mycoAIbert" & MinimumSeqQuality == "Qmin15")
mycoAIbert17 <- filter(GSD_Mock2_diversity, Classifier == "mycoAIbert" & MinimumSeqQuality == "Qmin17")
mycoAICNN15 <- filter(GSD_Mock2_diversity, Classifier == "mycoAICNN" & MinimumSeqQuality == "Qmin15")
mycoAICNN17 <- filter(GSD_Mock2_diversity, Classifier == "mycoAICNN" & MinimumSeqQuality == "Qmin17")
vtd15 <- filter(GSD_Mock2_diversity, Classifier == "vtd" & MinimumSeqQuality == "Qmin15")
vtd17 <- filter(GSD_Mock2_diversity, Classifier == "vtd" & MinimumSeqQuality == "Qmin17")

GSD_M2_splitlist <- list(aodp17, aodp15, dnabarcoder17, dnabarcoder15, emu15, emu17, kraken15, kraken17, minimap15, minimap17, mycoAIbert15, mycoAIbert17, mycoAICNN15, mycoAICNN17, vtd15, vtd17)

metric_cols <- c("DeltaSimpson", "DeltaShannon", "Bhattacharyya")
metadata_cols <- c("Classifier", "Database", "MinimumSeqQuality")
summary_list <- lapply(GSD_M2_splitlist, calculate_summary, metric_cols = metric_cols, metadata_cols = metadata_cols)
GSD_M2_summary <- do.call(rbind, summary_list)
GSD_M2_summary <- dplyr::select(GSD_M2_summary, c("Classifier", "Database", "MinimumSeqQuality", "Metric", "Mean", "CI_Lower", "CI_Upper"))

## ---------- NCBI data wrangling -----------
columnnames <- c("Species", "AssignedSpecies", "TotalCount")

NCBI_Mock1 <- list()
files_NCBIM1 <- list.files(path = './Mock1/NCBI/CSV', pattern = ".csv", full.names = TRUE)
for (i in seq_along(files_NCBIM1)) {
  file = files_NCBIM1[[i]]
  df <- read.csv(file, header = FALSE)
  colnames(df) <- columnnames
  df_name <- tools::file_path_sans_ext(basename(file)) 
  NCBI_Mock1[[df_name]] <- df
}
NCBI_Mock2 <- list()
files_NCBIM2 <- list.files(path = './Mock2/NCBI/CSV', pattern = ".csv", full.names = TRUE)
for (i in seq_along(files_NCBIM2)) {
  file = files_NCBIM2[[i]]
  df <- read.csv(file, header = FALSE)
  colnames(df) <- columnnames
  df_name <- tools::file_path_sans_ext(basename(file))
  NCBI_Mock2[[df_name]] <- df
}

# Get species-level assignment
for (i in seq_along(NCBI_Mock1)) {
  df <- NCBI_Mock1[[i]]
  df <- separate_wider_delim(df, AssignedSpecies, ",", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "AssignedSpecies"), too_few = "align_start", too_many = "merge")
  df <- dplyr::select(df, c("Species", "AssignedSpecies", "TotalCount"))
  NCBI_Mock1[[i]] <- df
}
for (i in seq_along(NCBI_Mock2)) {
  df <- NCBI_Mock2[[i]]
  df <- separate_wider_delim(df, AssignedSpecies, ",", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "AssignedSpecies"), too_few = "align_start", too_many = "merge")
  df <- dplyr::select(df, c("Species", "AssignedSpecies", "TotalCount"))
  NCBI_Mock2[[i]] <- df
}

# Because of synonymous nomenclature, rename some species so the matching function works
patterns <- c("Candida_boleticola" = "[Candida]_boleticola", 
              "Candida_caryicola" = "[Candida]_caryicola", 
              "Candida_catenulata" = "Diutina_catenulata", 
              "Candida_glabrata" = "Nakaseomyces_glabratus", 
              "Candida_zeylanoides" = "[Candida]_zeylanoides", 
              "Candida_orthopsilosis" = "Candida_metapsilosis", 
              "Leptosphaeria_maculans" = "Plenodomus_lingam", 
              "Cryptococcus_neoformans_VN_IIII" = "Cryptococcus_neoformans", 
              "Cryptococcus_gattii_VG_I" = "Cryptococcus_gattii")
replace_strings <- function(df, patterns) {
  df[] <- lapply(df, function(col) {
    if (is.character(col)) {
      stringr::str_replace_all(col, patterns) 
    } else {
      col
    }
  })
  return(df)
}

# Need a specific function for T. brumale, because the name appears in both columns and ran into problems trying to change it in the above function. 
species_replacement <- c("Tuber_brumale" = "Tuber_brumale_var._gorii")
replace_tuber <- function(df) {
  if ("Species" %in% names(df)) {
    df$Species <- stringr::str_replace_all(df$Species, species_replacement)
  }
  return(df)
}
NCBI_Mock1 <- lapply(NCBI_Mock1, replace_strings, patterns)
NCBI_Mock1 <- lapply(NCBI_Mock1, replace_tuber)
NCBI_Mock2 <- lapply(NCBI_Mock2, replace_strings, patterns)
NCBI_Mock2 <- lapply(NCBI_Mock2, replace_tuber)
NCBI_Mock1 <- lapply(NCBI_Mock1, function(df) {
  df %>% mutate(AssignedSpecies = ifelse(is.na(AssignedSpecies), "Unclassified", AssignedSpecies))
})
NCBI_Mock2 <- lapply(NCBI_Mock2, function(df) {
  df %>% mutate(AssignedSpecies = ifelse(is.na(AssignedSpecies), "Unclassified", AssignedSpecies))
})

# Calculate the unclassified sequences
# Ensure Mock1_GT is accessible and correct
total_ground_truth <- sum(Mock1_GT$Mock1_GroundTruth_GroundTruth_Qmin15)
cat("Total ground truth sequences:", total_ground_truth, "\n")

# Apply the function to each data frame in NCBI_Mock1
NCBI_Mock1 <- lapply(NCBI_Mock1, function(df) {
  # Calculate current total classified sequences
  current_total <- sum(df$TotalCount)
  cat("Current total sequences:", current_total, "\n")
  # Calculate difference from ground truth
  difference <- total_ground_truth - current_total
  cat("Difference from ground truth:", difference, "\n")
  # Find Unclassified rows
  unclass_rows <- which(df$AssignedSpecies == "Unclassified")
  if (length(unclass_rows) > 0) {
    # Handle existing Unclassified rows
    if (length(unclass_rows) > 1) {
      warning("Multiple Unclassified rows detected! Consolidating counts.")
      # Sum existing Unclassified counts
      current_unclass_count <- sum(df$TotalCount[unclass_rows])
      # Remove extra Unclassified rows, keep first
      df <- df[-unclass_rows[-1], ]
      unclass_rows <- unclass_rows[1]
    } else {
      current_unclass_count <- df$TotalCount[unclass_rows]
    }
    # Update Unclassified count
    new_unclass_count <- current_unclass_count + difference
    cat("Current Unclassified count:", current_unclass_count, "\n")
    cat("New Unclassified count:", new_unclass_count, "\n")
    # Ensure non-negative
    if (new_unclass_count < 0) {
      warning("Negative Unclassified count detected! Setting to 0. Check data consistency.")
      new_unclass_count <- 0
    }
    df$TotalCount[unclass_rows] <- new_unclass_count
  } else {
    # No existing Unclassified row
    new_unclass_count <- difference
    cat("New Unclassified count:", new_unclass_count, "\n")
    # Ensure non-negative
    if (new_unclass_count < 0) {
      warning("Negative Unclassified count detected! Setting to 0. Check data consistency.")
      new_unclass_count <- 0
    }
    # Create new Unclassified row
    new_row <- data.frame(
      Species = "Unclassified",
      AssignedSpecies = "Unclassified",
      TotalCount = new_unclass_count
    )
    df <- rbind(df, new_row)
  }
  # Verify final total
  final_total <- sum(df$TotalCount)
  cat("Final total sequences:", final_total, "\n\n")
  return(df)
})

# Repeat for Mock 2
# Ensure Mock2_GT is accessible and correct
total_ground_truth2 <- sum(Mock2_GT$Mock2_GroundTruth_GroundTruth_Qmin15)
cat("Total ground truth sequences:", total_ground_truth2, "\n")

# Apply the function to each data frame in NCBI_Mock2
NCBI_Mock2 <- lapply(NCBI_Mock2, function(df) {
  # Calculate current total classified sequences
  current_total <- sum(df$TotalCount)
  cat("Current total sequences:", current_total, "\n")
  # Calculate difference from ground truth
  difference <- total_ground_truth2 - current_total
  cat("Difference from ground truth:", difference, "\n")
  # Find Unclassified rows
  unclass_rows <- which(df$AssignedSpecies == "Unclassified")
  if (length(unclass_rows) > 0) {
    # Handle existing Unclassified rows
    if (length(unclass_rows) > 1) {
      warning("Multiple Unclassified rows detected! Consolidating counts.")
      # Sum existing Unclassified counts
      current_unclass_count <- sum(df$TotalCount[unclass_rows])
      # Remove extra Unclassified rows, keep first
      df <- df[-unclass_rows[-1], ]
      unclass_rows <- unclass_rows[1]
    } else {
      current_unclass_count <- df$TotalCount[unclass_rows]
    }
    # Update Unclassified count
    new_unclass_count <- current_unclass_count + difference
    cat("Current Unclassified count:", current_unclass_count, "\n")
    cat("New Unclassified count:", new_unclass_count, "\n")
    # Ensure non-negative
    if (new_unclass_count < 0) {
      warning("Negative Unclassified count detected! Setting to 0. Check data consistency.")
      new_unclass_count <- 0
    }
    df$TotalCount[unclass_rows] <- new_unclass_count
  } else {
    # No existing Unclassified row
    new_unclass_count <- difference
    cat("New Unclassified count:", new_unclass_count, "\n")
    # Ensure non-negative
    if (new_unclass_count < 0) {
      warning("Negative Unclassified count detected! Setting to 0. Check data consistency.")
      new_unclass_count <- 0
    }
    # Create new Unclassified row
    new_row <- data.frame(
      Species = "Unclassified",
      AssignedSpecies = "Unclassified",
      TotalCount = new_unclass_count
    )
    df <- rbind(df, new_row)
  }
  # Verify final total
  final_total <- sum(df$TotalCount)
  cat("Final total sequences:", final_total, "\n\n")
  return(df)
})

# Sanity check!
# Check each data frame in NCBI_Mock1
lapply(seq_along(NCBI_Mock1), function(i) {
  df <- NCBI_Mock1[[i]]
  total_in_df <- sum(df$TotalCount)
  diff <- total_ground_truth - total_in_df
  cat(sprintf("Classifier %d: Total sequences = %d, Ground truth = %d, Difference = %d\n", 
              i, total_in_df, total_ground_truth, diff))
  if (abs(diff) > 0.001) {  # Allowing for small floating-point errors
    warning(sprintf("Classifier %d: Sum of TotalCount (%d) does not match ground truth (%d)", 
                    i, total_in_df, total_ground_truth))
  }
})

# Check each data frame in NCBI_Mock2
lapply(seq_along(NCBI_Mock2), function(i) {
  df <- NCBI_Mock2[[i]]
  total_in_df <- sum(df$TotalCount)
  diff <- total_ground_truth2 - total_in_df
  cat(sprintf("Classifier %d: Total sequences = %d, Ground truth = %d, Difference = %d\n", 
              i, total_in_df, total_ground_truth2, diff))
  if (abs(diff) > 0.001) {  # Allowing for small floating-point errors
    warning(sprintf("Classifier %d: Sum of TotalCount (%d) does not match ground truth (%d)", 
                    i, total_in_df, total_ground_truth2))
  }
})

NCBIspecies <- read_xlsx("./SpeciesInNCBI_Taxonomy.xlsx", col_names = TRUE)
NCBIspecies <- dplyr::select(NCBIspecies, "organism") 
NCBIspecies <- replace_strings(NCBIspecies, patterns)
NCBIspecies <- replace_tuber(NCBIspecies)
colnames(NCBIspecies) <- "Species"

Mock1_notinNCBI <- lapply(NCBI_Mock1, function(df) {
  df %>%
    filter(!(Species %in% NCBIspecies$Species))
})
Mock2_notinNCBI <- lapply(NCBI_Mock2, function(df) {
  df %>% 
    filter(!(Species %in% NCBIspecies$Species))
})

reformat_results_NCBI <- function(df) {
  # Ensure all AssignedSpecies appear in MockSpecies
  all_species <- unique(df$AssignedSpecies)
  df %>%
    group_by(AssignedSpecies) %>%
    summarise(
      CorrectlyMapped = sum(TotalCount[Species == AssignedSpecies], na.rm = TRUE),
      TotalMapped = sum(TotalCount), 
      .groups = "drop"
    ) %>%
    rename(MockSpecies = AssignedSpecies) %>%
    # Ensure all AssignedSpecies are included, even those without matches in Species
    complete(MockSpecies = all_species, fill = list(CorrectlyMapped = 0, TotalMapped = 0))
}
NCBI_Mock1_reformat <- map(NCBI_Mock1, reformat_results_NCBI)
# Using this method, the 'Unclassified' CorrectlyMapped count will be >0. Change back to 0. 
NCBI_Mock1_reformat <- lapply(NCBI_Mock1_reformat, function(df) {
  df[df$MockSpecies == "Unclassified", "CorrectlyMapped"] <- 0
  return(df)
})
NCBI_Mock2_reformat <- map(NCBI_Mock2, reformat_results_NCBI)
NCBI_Mock2_reformat <- lapply(NCBI_Mock2_reformat, function(df) {
  df[df$MockSpecies == "Unclassified", "CorrectlyMapped"] <- 0
  return(df)
})

# Create BIOM tables
all_species <- unique(unlist(lapply(NCBI_Mock1, function(df) df$AssignedSpecies))) 
NCBI_Mock1_OTUtable <- data.frame(Species = all_species)
for (i in seq_along(NCBI_Mock1_reformat)) {
  df <- NCBI_Mock1_reformat[[i]]
  temp <- data.frame(Species = all_species)
  temp$TotalCount <- df$TotalMapped[match(all_species, df$MockSpecies)]
  temp$TotalCount <- round(as.numeric(temp$TotalCount))
  NCBI_Mock1_OTUtable[[names(NCBI_Mock1_reformat)[i]]] <- temp$TotalCount
  NCBI_Mock1_reformat[[i]] <- df
}
colnames(Mock1_GT) <- c("Species", "Mock1_GroundTruth_ncbi_GroundTruth_Qmin15") # Change the column name so it'll split correctly along with the other data
NCBI_Mock1_OTUtable <- left_join(NCBI_Mock1_OTUtable, Mock1_GT, by = "Species")
NCBI_Mock1_OTUtable[is.na(NCBI_Mock1_OTUtable)] <- 0
NCBI_Mock1_OTUtable[, 2:ncol(NCBI_Mock1_OTUtable)] <- lapply(NCBI_Mock1_OTUtable[, 2:ncol(NCBI_Mock1_OTUtable)], function(x) round(as.numeric(x)))

# Calculate diversity metrics
NCBI_Mock1_OTUmatrix <- as.matrix(NCBI_Mock1_OTUtable[, -1])
# Calculate the Bhattacharyya distance while the ground truth column is a part of the matrix
ground_truth <- NCBI_Mock1_OTUtable[, "Mock1_GroundTruth_ncbi_GroundTruth_Qmin15"]
norm_matrix <- apply(NCBI_Mock1_OTUmatrix, 2, function(x) x / sum(x))
norm_GT <- ground_truth / sum(ground_truth)
bhattacharyya_NCBIM1 <- apply(norm_matrix, 2, function(x) bhattacharyya_distance(x, norm_GT))
# Remove the ground truth and calculate the alpha diversity metrics separately. Zero-abundance species will change the diversity metric calculation for the GT, so use the value calculated from the GSD dataset chunk above. 
NCBI_Mock1_OTUmatrix <- NCBI_Mock1_OTUmatrix[, colnames(NCBI_Mock1_OTUmatrix) != "Mock1_GroundTruth_ncbi_GroundTruth_Qmin15"]
shannon_NCBIM1 <- apply(NCBI_Mock1_OTUmatrix, 2, function(x) diversity(x, index = "shannon"))
simpson_NCBIM1 <- apply(NCBI_Mock1_OTUmatrix, 2, function(x) diversity(x, index = "simpson")) 
# Append the ground truth alpha diversity calculations to the vectors. No need to re-calculate.
shannon_NCBIM1 <- append(shannon_NCBIM1, shannon_M1GT, after = length(shannon_NCBIM1))
simpson_NCBIM1 <- append(simpson_NCBIM1, simpson_M1GT, after = length(simpson_NCBIM1))
NCBI_Mock1_diversity <- data.frame(
  Sample = c(colnames(NCBI_Mock1_OTUmatrix), "Mock1_GroundTruth_ncbi_GroundTruth_Qmin15"), 
  Shannon = shannon_NCBIM1, 
  Simpson = simpson_NCBIM1, 
  Bhattacharyya = bhattacharyya_NCBIM1
)
rownames(NCBI_Mock1_diversity) <- NULL
NCBI_Mock1_diversity <- separate_wider_delim(NCBI_Mock1_diversity, Sample, "_", names = c("MockCommunity", "Classifier", "Database", "SubRep", "MinimumSeqQuality"), too_few = "align_start")

# Repeat for Mock 2
all_species <- unique(unlist(lapply(NCBI_Mock2, function(df) df$AssignedSpecies))) 
NCBI_Mock2_OTUtable <- data.frame(Species = all_species)
for (i in seq_along(NCBI_Mock2_reformat)) {
  df <- NCBI_Mock2_reformat[[i]]
  temp <- data.frame(Species = all_species)
  temp$TotalCount <- df$TotalMapped[match(all_species, df$MockSpecies)]
  temp$TotalCount <- round(as.numeric(temp$TotalCount))
  NCBI_Mock2_OTUtable[[names(NCBI_Mock2_reformat)[i]]] <- temp$TotalCount
  NCBI_Mock2_reformat[[i]] <- df
}
colnames(Mock2_GT) <- c("Species", "Mock2_GroundTruth_ncbi_GroundTruth_Qmin15") # Change the column name so it'll split correctly along with the other data
NCBI_Mock2_OTUtable <- left_join(NCBI_Mock2_OTUtable, Mock2_GT, by = "Species")
# Calculate whether there are any 'unclassified' sequences, given that the total number is known. For Mock 2, this is 21,779 sequences total in the Mock Community
total_sequences_M2 <- 21779
unclassified_counts <- total_sequences_M2 - colSums(NCBI_Mock2_OTUtable[ , -1])
NCBI_Mock2_OTUtable <- rbind(NCBI_Mock2_OTUtable, c("Unclassified", unclassified_counts))
NCBI_Mock2_OTUtable[is.na(NCBI_Mock2_OTUtable)] <- 0
NCBI_Mock2_OTUtable[, 2:ncol(NCBI_Mock2_OTUtable)] <- lapply(NCBI_Mock2_OTUtable[, 2:ncol(NCBI_Mock2_OTUtable)], function(x) round(as.numeric(x)))

# Calculate diversity metrics for Mock 2
NCBI_Mock2_OTUmatrix <- as.matrix(NCBI_Mock2_OTUtable[, -1])

# Calculate Bhattacharyya distance with the ground truth data still in the matrix.
ground_truth <- NCBI_Mock2_OTUtable[, "Mock2_GroundTruth_ncbi_GroundTruth_Qmin15"]
norm_matrix <- apply(NCBI_Mock2_OTUmatrix, 2, function(x) x / sum(x))
norm_GT <- ground_truth / sum(ground_truth)
bhattacharyya_NCBIM2 <- apply(norm_matrix, 2, function(x) bhattacharyya_distance(x, norm_GT))
# Remove the ground truth and calculate the alpha diversity metrics
NCBI_Mock2_OTUmatrix <- NCBI_Mock2_OTUmatrix[, colnames(NCBI_Mock2_OTUmatrix) != "Mock2_GroundTruth_ncbi_GroundTruth_Qmin15"]
shannon_NCBIM2 <- apply(NCBI_Mock2_OTUmatrix, 2, function(x) diversity(x, index = "shannon"))
simpson_NCBIM2 <- apply(NCBI_Mock2_OTUmatrix, 2, function(x) diversity(x, index = "simpson")) 
# Append the ground truth alpha diversity metrics for Mock 2. No need to re-calculate. 
shannon_NCBIM2 <- append(shannon_NCBIM2, shannon_M2GT, after = length(shannon_NCBIM2))
simpson_NCBIM2 <- append(simpson_NCBIM2, simpson_M2GT, after = length(simpson_NCBIM2))
NCBI_Mock2_diversity <- data.frame(
  Sample = c(colnames(NCBI_Mock2_OTUmatrix), "Mock2_GroundTruth_ncbi_GroundTruth_Qmin15"), 
  Shannon = shannon_NCBIM2, 
  Simpson = simpson_NCBIM2, 
  Bhattacharyya = bhattacharyya_NCBIM2
)
rownames(NCBI_Mock2_diversity) <- NULL
NCBI_Mock2_diversity <- separate_wider_delim(NCBI_Mock2_diversity, Sample, "_", names = c("MockCommunity", "Classifier", "Database", "SubRep", "MinimumSeqQuality"), too_few = "align_start")

# Calculate delta Shannon and delta Simpson
GT_shannon <- NCBI_Mock1_diversity$Shannon[[nrow(NCBI_Mock1_diversity)]]
GT_simpson <- NCBI_Mock1_diversity$Simpson[[nrow(NCBI_Mock1_diversity)]]
NCBI_Mock1_diversity$DeltaShannon <- NCBI_Mock1_diversity$Shannon - GT_shannon
NCBI_Mock1_diversity$DeltaSimpson <- NCBI_Mock1_diversity$Simpson - GT_simpson

GT_shannon <- NCBI_Mock2_diversity$Shannon[[nrow(NCBI_Mock2_diversity)]]
GT_simpson <- NCBI_Mock2_diversity$Simpson[[nrow(NCBI_Mock2_diversity)]]
NCBI_Mock2_diversity$DeltaShannon <- NCBI_Mock2_diversity$Shannon - GT_shannon
NCBI_Mock2_diversity$DeltaSimpson <- NCBI_Mock2_diversity$Simpson - GT_simpson

# Use technical replicates to calculate mean and 95% CI's
dnabarcoder15 <- filter(NCBI_Mock1_diversity, Classifier == "dnabarcoder" & MinimumSeqQuality == "Qmin15")
dnabarcoder17 <- filter(NCBI_Mock1_diversity, Classifier == "dnabarcoder" & MinimumSeqQuality == "Qmin17")
emu15 <- filter(NCBI_Mock1_diversity, Classifier == "emu" & MinimumSeqQuality == "Qmin15")
emu17 <- filter(NCBI_Mock1_diversity, Classifier == "emu" & MinimumSeqQuality == "Qmin17")
kraken15 <- filter(NCBI_Mock1_diversity, Classifier == "kraken2" & MinimumSeqQuality == "Qmin15")
kraken17 <- filter(NCBI_Mock1_diversity, Classifier == "kraken2" & MinimumSeqQuality == "Qmin17")
minimap15 <- filter(NCBI_Mock1_diversity, Classifier == "minimap2" & MinimumSeqQuality == "Qmin15")
minimap17 <- filter(NCBI_Mock1_diversity, Classifier == "minimap2" & MinimumSeqQuality == "Qmin17")
mycoAIbert15 <- filter(NCBI_Mock1_diversity, Classifier == "mycoAIbert" & MinimumSeqQuality == "Qmin15")
mycoAIbert17 <- filter(NCBI_Mock1_diversity, Classifier == "mycoAIbert" & MinimumSeqQuality == "Qmin17")
mycoAICNN15 <- filter(NCBI_Mock1_diversity, Classifier == "mycoAICNN" & MinimumSeqQuality == "Qmin15")
mycoAICNN17 <- filter(NCBI_Mock1_diversity, Classifier == "mycoAICNN" & MinimumSeqQuality == "Qmin17")
vtd15 <- filter(NCBI_Mock1_diversity, Classifier == "vtd" & MinimumSeqQuality == "Qmin15")
vtd17 <- filter(NCBI_Mock1_diversity, Classifier == "vtd" & MinimumSeqQuality == "Qmin17")
NCBI_M1_splitlist <- list(dnabarcoder17, dnabarcoder15, emu15, emu17, kraken15, kraken17, minimap15, minimap17, mycoAIbert15, mycoAIbert17, mycoAICNN15, mycoAICNN17, vtd15, vtd17)
metric_cols <- c("DeltaSimpson", "DeltaShannon", "Bhattacharyya")
metadata_cols <- c("Classifier", "Database", "MinimumSeqQuality")
summary_list <- lapply(NCBI_M1_splitlist, calculate_summary, metric_cols = metric_cols, metadata_cols = metadata_cols)
NCBI_M1_summary <- do.call(rbind, summary_list)
NCBI_M1_summary <- dplyr::select(NCBI_M1_summary, c("Classifier", "Database", "MinimumSeqQuality", "Metric", "Mean", "CI_Lower", "CI_Upper"))

dnabarcoder15 <- filter(NCBI_Mock2_diversity, Classifier == "dnabarcoder" & MinimumSeqQuality == "Qmin15")
dnabarcoder17 <- filter(NCBI_Mock2_diversity, Classifier == "dnabarcoder" & MinimumSeqQuality == "Qmin17")
emu15 <- filter(NCBI_Mock2_diversity, Classifier == "emu" & MinimumSeqQuality == "Qmin15")
emu17 <- filter(NCBI_Mock2_diversity, Classifier == "emu" & MinimumSeqQuality == "Qmin17")
kraken15 <- filter(NCBI_Mock2_diversity, Classifier == "kraken2" & MinimumSeqQuality == "Qmin15")
kraken17 <- filter(NCBI_Mock2_diversity, Classifier == "kraken2" & MinimumSeqQuality == "Qmin17")
minimap15 <- filter(NCBI_Mock2_diversity, Classifier == "minimap2" & MinimumSeqQuality == "Qmin15")
minimap17 <- filter(NCBI_Mock2_diversity, Classifier == "minimap2" & MinimumSeqQuality == "Qmin17")
mycoAIbert15 <- filter(NCBI_Mock2_diversity, Classifier == "mycoAIbert" & MinimumSeqQuality == "Qmin15")
mycoAIbert17 <- filter(NCBI_Mock2_diversity, Classifier == "mycoAIbert" & MinimumSeqQuality == "Qmin17")
mycoAICNN15 <- filter(NCBI_Mock2_diversity, Classifier == "mycoAICNN" & MinimumSeqQuality == "Qmin15")
mycoAICNN17 <- filter(NCBI_Mock2_diversity, Classifier == "mycoAICNN" & MinimumSeqQuality == "Qmin17")
vtd15 <- filter(NCBI_Mock2_diversity, Classifier == "vtd" & MinimumSeqQuality == "Qmin15")
vtd17 <- filter(NCBI_Mock2_diversity, Classifier == "vtd" & MinimumSeqQuality == "Qmin17")
NCBI_M2_splitlist <- list(dnabarcoder17, dnabarcoder15, emu15, emu17, kraken15, kraken17, minimap15, minimap17, mycoAIbert15, mycoAIbert17, mycoAICNN15, mycoAICNN17, vtd15, vtd17)
summary_list <- lapply(NCBI_M2_splitlist, calculate_summary, metric_cols = metric_cols, metadata_cols = metadata_cols)
NCBI_M2_summary <- do.call(rbind, summary_list)
NCBI_M2_summary <- dplyr::select(NCBI_M2_summary, c("Classifier", "Database", "MinimumSeqQuality", "Metric", "Mean", "CI_Lower", "CI_Upper"))

## ---------- Plotting means and CIs -----------
NCBI_M1_summary$MockCommunity <- "Mock1"
NCBI_M2_summary$MockCommunity <- "Mock2"
GSD_M1_summary$MockCommunity <- "Mock1"
GSD_M2_summary$MockCommunity <- "Mock2"
GSD_full_data <- rbind(GSD_M1_summary, GSD_M2_summary)
NCBI_full_data <- rbind(NCBI_M1_summary, NCBI_M2_summary)

# Filter data frames by metric and sequence quality. Use Qmin15 for now. Come back later after statistical testing if there's a significant difference between Qmin15 and Qmin17.
shannon_GSD <- filter(GSD_full_data, Metric == "DeltaShannon" & MinimumSeqQuality == "Qmin15")
simpson_GSD <- filter(GSD_full_data, Metric == "DeltaSimpson" & MinimumSeqQuality == "Qmin15")
bhattacharyya_GSD <- filter(GSD_full_data, Metric == "Bhattacharyya" & MinimumSeqQuality == "Qmin15")
shannon_NCBI <- filter(NCBI_full_data, Metric == "DeltaShannon" & MinimumSeqQuality == "Qmin15")
simpson_NCBI <- filter(NCBI_full_data, Metric == "DeltaSimpson" & MinimumSeqQuality == "Qmin15")
bhattacharyya_NCBI <- filter(NCBI_full_data, Metric == "Bhattacharyya" & MinimumSeqQuality == "Qmin15")

library(ggforce)
library(scales)
library(RColorBrewer)

# Define which colour in the Dark2 palette from RColorBrewer corresponds to which classifier, so this is kept consistent throughout plots.
classifier_colours <- c("minimap2" = "#1B9E77",
                        "kraken2" = "#D95F02", 
                        "aodp" = "#7570B3", 
                        "dnabarcoder" = "#E7298A", 
                        "emu" = "#66A61E", 
                        "vtd" = "#E6AB02", 
                        "mycoAIbert" = "#A6761D", 
                        "mycoAICNN" = "#666666")

# Define the order for the classifiers to appear in
classifier_order <- c("aodp", "minimap2", "kraken2", "dnabarcoder", "emu", "vtd", "mycoAIbert", "mycoAICNN")

plot_metric <- function(data, facet_var, output_plot, plot_subtitle, plot_title, y_label, plot_palette = "Dark2") {
  # Make sure factors appear in the right order
  data$Classifier <- factor(data$Classifier, levels = classifier_order)
  
  metric_dotplot <- ggplot(data, aes(x = factor(Classifier), y = Mean, color = Classifier)) +
    geom_dotplot(binaxis = "y", stackdir = "center") +
    facet_grid(. ~ get(facet_var)) + # Use the variable set dynamically in the function. 
    theme(panel.background = element_rect(fill = "white"), panel.grid = element_line(color = "lightgrey")) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, linewidth = 0.7) +
    labs(title = plot_title, subtitle = plot_subtitle, x = "Classifier", y = y_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 14), 
          plot.title = element_text(size = 16), 
          axis.text.y = element_text(size = 12, color = "black"), 
          axis.title.x = element_blank()) +
    scale_y_continuous(breaks = pretty_breaks(n = 8)) + 
    scale_color_manual(values = classifier_colours) +
    theme(legend.position = "none") 
  return(metric_dotplot)
}

shannon_GSD_dotplot <- plot_metric(shannon_GSD, facet_var = "MockCommunity", output_plot = "20250514_GSD_ShannonDotplot_AlignOnly.png", plot_subtitle = "Gold Standard Database", plot_title = "Delta Shannon Index for Assessed Classifiers", y_label = "Delta Shannon Index") 
simpson_GSD_dotplot <- plot_metric(simpson_GSD, facet_var = "MockCommunity", output_plot = "20250514_GSD_SimpsonDotplot_AlignOnly.png", plot_subtitle = "Gold Standard Database", plot_title = "Delta Simpson Index for Assessed Classifiers", y_label = "Delta Simpson Index")
bhattacharyya_GSD_dotplot <- plot_metric(bhattacharyya_GSD, facet_var = "MockCommunity", output_plot = "20250514_GSD_BhattacharyyaDotplot_AlignOnly.png", plot_subtitle = "Gold Standard Database", plot_title = "Bhattacharyya Distance for Assessed Classifiers", y_label = "Bhattacharyya Distance")
shannon_NCBI_dotplot <- plot_metric(shannon_NCBI, facet_var = "MockCommunity", output_plot = "20250514_NCBI_ShannonDotplot_AlignOnly.png", plot_subtitle = "NCBI RefSeq ITS Database", plot_title = "Delta Shannon Index for Assessed Classifiers", y_label = "Delta Shannon Index")
simpson_NCBI_dotplot <- plot_metric(simpson_NCBI, facet_var = "MockCommunity", output_plot = "20250514_NCBI_SimpsonDotplot_AlignOnly.png", plot_subtitle = "NCBI RefSeq ITS Database", plot_title = "Delta Simpson Index for Assessed Classifiers", y_label = "Delta Simpson Index")
bhattacharyya_NCBI_dotplot <- plot_metric(bhattacharyya_NCBI, facet_var = "MockCommunity", output_plot = "20250214_NCBI_BhattacharyyaDotplot_AlignOnly.png", plot_subtitle = "NCBI RefSeq ITS Database", plot_title = "Bhattacharyya Distance for Assessed Classifiers", y_label = "Bhattacharyya Distance")

## ----------- Statistical Analysis -----------
fullresults <- bind_rows(GSD_Mock1_diversity, GSD_Mock2_diversity, NCBI_Mock1_diversity, NCBI_Mock2_diversity)
print(fullresults)

M1div <- bind_rows(GSD_Mock1_diversity, NCBI_Mock1_diversity)
M2div <- bind_rows(GSD_Mock2_diversity, NCBI_Mock2_diversity)

data_list <- list(GSD_Mock1_diversity, GSD_Mock2_diversity, NCBI_Mock1_diversity, NCBI_Mock2_diversity)

# Iterate over each distribution separately to determine normality of data distribution.
shapiro.bhattacharyya <- list()
shapiro.deltashannon <- list()

for (i in seq_along(data_list)) {
  shapiro.bhattacharyya[[i]] <- shapiro.test(data_list[[i]]$Bhattacharyya)
  shapiro.deltashannon[[i]] <- shapiro.test(data_list[[i]]$DeltaShannon) #Result for DeltaShannon and Shannon will be the same, as inherently the data will follow the same distribution, even if the values are different. 
}

adonis.M1.shannon <- adonis2(M1div$Shannon ~ M1div$Classifier + M1div$Database + M1div$MinimumSeqQuality, M1div, method = "bray", by = "margin")
adonis.M2.shannon <- adonis2(M2div$Shannon ~ M2div$Classifier + M2div$Database + M2div$MinimumSeqQuality, M2div, method = "bray", by = "margin")
adonis.M1.bhattacharyya <- adonis2(M1div$Bhattacharyya ~ M1div$Classifier + M1div$Database + M1div$MinimumSeqQuality, M1div, method = "euclidean", by = "margin")
adonis.M2.bhattacharyya <- adonis2(M2div$Bhattacharyya ~ M2div$Classifier + M2div$Database + M2div$MinimumSeqQuality, M2div, method = "euclidean", by = "margin")

# For DeltaShannon, use Euclidean distance, because some values are negative. Bray-Curtis won't accurately handle these.
adonis.M1.deltashannon <- adonis2(M1div$DeltaShannon ~ M1div$Classifier + M1div$Database + M1div$MinimumSeqQuality, M1div, method = "euclidean", by = "margin")
adonis.M2.deltashannon <- adonis2(M2div$DeltaShannon ~ M2div$Classifier + M2div$Database + M2div$MinimumSeqQuality, M2div, method = "euclidean", by = "margin")

library(DescTools)
GSD_M1 <- filter(fullresults, MockCommunity == "Mock1") %>% filter(Database == "Gold Standard Database") %>% filter(MinimumSeqQuality == "Qmin15")
NCBI_M1 <- filter(fullresults, MockCommunity == "Mock1") %>% filter(Database == "ncbi") %>% filter(MinimumSeqQuality == "Qmin15")
GSD_M2 <- filter(fullresults, MockCommunity == "Mock2") %>% filter(Database == "Gold Standard Database") %>% filter(MinimumSeqQuality == "Qmin15")
NCBI_M2 <- filter(fullresults, MockCommunity == "Mock2") %>% filter(Database == "ncbi") %>% filter(MinimumSeqQuality == "Qmin15")

dunn.M1GSD.deltashannon <- DunnTest(GSD_M1$DeltaShannon, GSD_M1$Classifier, method = "BH")
dunn.M1GSD.bhattacharyya <- DunnTest(GSD_M1$Bhattacharyya, GSD_M1$Classifier, method = "BH")
print(dunn.M1GSD.bhattacharyya)
print(dunn.M1GSD.deltashannon)

dunn.M2GSD.deltashannon <- DunnTest(GSD_M2$DeltaShannon, GSD_M2$Classifier, method = "BH")
dunn.M2GSD.bhattacharyya <- DunnTest(GSD_M2$Bhattacharyya, GSD_M2$Classifier, method = "BH")
print(dunn.M2GSD.bhattacharyya)
print(dunn.M2GSD.deltashannon)

dunn.M1NCBI.deltashannon <- DunnTest(NCBI_M1$DeltaShannon, NCBI_M1$Classifier, method = "BH")
dunn.M1NCBI.bhattacharyya <- DunnTest(NCBI_M1$Bhattacharyya, NCBI_M1$Classifier, method = "BH")
print(dunn.M1NCBI.bhattacharyya)
print(dunn.M1NCBI.deltashannon)

dunn.M2NCBI.deltashannon <- DunnTest(NCBI_M2$DeltaShannon, NCBI_M2$Classifier, method = "BH")
dunn.M2NCBI.bhattacharyya <- DunnTest(NCBI_M2$Bhattacharyya, NCBI_M2$Classifier, method = "BH")
print(dunn.M2NCBI.bhattacharyya)
print(dunn.M2NCBI.deltashannon)

# Generate CLD and add to plots
library(rcompanion)
library(FSA)

# Remove the ground truth for the purpose of generating the CLD. This test won't be represented in the plot.
GSD_M1 <- filter(GSD_M1, Classifier != "GroundTruth")
GSD_M2 <- filter(GSD_M2, Classifier != "GroundTruth")

KW_GSDM1_DS <- kruskal.test(formula = DeltaShannon ~ factor(Classifier), data = GSD_M1)
Dunn_GSDM1_DS <- dunnTest(x = GSD_M1$DeltaShannon, g = GSD_M1$Classifier, method = "sidak")
Dunn_GSDM1_DS
Dun_GSDM1_DS.res <- Dunn_GSDM1_DS$res
Dunncld <- cldList(comparison = Dun_GSDM1_DS.res$Comparison, p.value = Dun_GSDM1_DS.res$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "Classifier"
grouped_table_GSDM1 <- group_by(GSD_M1, Classifier) %>% summarise(mean=mean(DeltaShannon, na.rm = TRUE), third_quantile = quantile(DeltaShannon, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_GSDM1 <- left_join(grouped_table_GSDM1, Dunncld, by = "Classifier")
grouped_table_GSDM1$MockCommunity <- "Mock1"
print(grouped_table_GSDM1)

KW_GSDM2_DS <- kruskal.test(formula = DeltaShannon ~ factor(Classifier), data = GSD_M2)
Dunn_GSDM2_DS <- dunnTest(x = GSD_M2$DeltaShannon, g = GSD_M2$Classifier, method = "sidak")
Dunn_GSDM2_DS
Dun_GSDM2_DS.res <- Dunn_GSDM2_DS$res
Dunncld <- cldList(comparison = Dun_GSDM2_DS.res$Comparison, p.value = Dun_GSDM2_DS.res$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "Classifier"
grouped_table_GSDM2 <- group_by(GSD_M2, Classifier) %>% summarise(mean=mean(DeltaShannon, na.rm = TRUE), third_quantile = quantile(DeltaShannon, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_GSDM2 <- left_join(grouped_table_GSDM2, Dunncld, by = "Classifier")
grouped_table_GSDM2$MockCommunity <- "Mock2"
print(grouped_table_GSDM2)

shannon_GSD_min15 <- filter(shannon_GSD, MinimumSeqQuality == "Qmin15")
# Make sure classifiers are in the right order for plotting
shannon_GSD_min15$Classifier <- factor(shannon_GSD_min15$Classifier, levels = classifier_order)

grouped_table_GSD <- rbind(grouped_table_GSDM1, grouped_table_GSDM2)
# Add position for CLD letters based on the position of the error bar for each classifier, plus a buffer of 5%
grouped_table_GSD <- grouped_table_GSD %>%
  left_join(shannon_GSD_min15 %>% dplyr::select(Classifier, MockCommunity, CI_Upper), 
            by = c("Classifier", "MockCommunity"))
y_range <- range(shannon_GSD_min15$CI_Upper, na.rm = TRUE)
buffer <- diff(y_range) * 0.05
grouped_table_GSD <- grouped_table_GSD %>%
  mutate(label_y = CI_Upper + buffer)
grouped_table_GSD$Classifier <- factor(grouped_table_GSD$Classifier, levels = classifier_order)

library(ggh4x)
deltashannon_dotplot_CLD <- ggplot(shannon_GSD_min15, aes(x = factor(Classifier), y = Mean, color = Classifier)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  facet_grid2(. ~ MockCommunity, axes = "all") +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, linewidth = 0.7) +
  geom_text(data = grouped_table_GSD, 
            aes(x = Classifier, y = label_y, label = Letter), 
            size = 3.8, 
            vjust = 0, 
            color = "black", 
            inherit.aes = FALSE) +
  labs(subtitle = "Gold Standard Database", 
       y = "Delta Shannon Index") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"), 
        #panel.grid = element_line(color = "lightgrey"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"), 
        legend.position = "none", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.title.y = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(), 
        axis.ticks = element_line(color = "black")) +
  scale_y_continuous(breaks = pretty_breaks(n = 8), 
                     expand = expansion(mult = c(0.05,0.2))) +
  scale_color_manual(values = classifier_colours) +
  scale_x_discrete(labels = c(
    "aodp" = "aodp", 
    "minimap2" = "minimap2",
    "kraken2" = "Kraken2",
    "emu" = "Emu",
    "dnabarcoder" = "Dnabarcoder",
    "vtd" = "vtd", 
    "mycoAICNN" = "mycoAI CNN",
    "mycoAIbert" = "mycoAI BERT"))
deltashannon_dotplot_CLD
ggsave("./DeltaShannonDotplot_GSD_CLD.png", deltashannon_dotplot_CLD, bg = "white", width = 7.29, height = 4.51, units = "in")

# Remove the ground truth
NCBI_M1 <- filter(NCBI_M1, Classifier != "GroundTruth")
NCBI_M2 <- filter(NCBI_M2, Classifier != "GroundTruth")

KW_NCBIM1_DS <- kruskal.test(formula = DeltaShannon ~ factor(Classifier), data = NCBI_M1)
Dunn_NCBIM1_DS <- dunnTest(x = NCBI_M1$DeltaShannon, g = NCBI_M1$Classifier, method = "sidak")
Dunn_NCBIM1_DS
Dun_NCBIM1_DS.res <- Dunn_NCBIM1_DS$res
Dunncld <- cldList(comparison = Dun_NCBIM1_DS.res$Comparison, p.value = Dun_NCBIM1_DS.res$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "Classifier"
grouped_table_NCBIM1 <- group_by(NCBI_M1, Classifier) %>% summarise(mean=mean(DeltaShannon, na.rm = TRUE), third_quantile = quantile(DeltaShannon, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_NCBIM1 <- left_join(grouped_table_NCBIM1, Dunncld, by = "Classifier")
grouped_table_NCBIM1$MockCommunity <- "Mock1"
print(grouped_table_NCBIM1)

KW_NCBIM2_DS <- kruskal.test(formula = DeltaShannon ~ factor(Classifier), data = NCBI_M2)
Dunn_NCBIM2_DS <- dunnTest(x = NCBI_M2$DeltaShannon, g = NCBI_M2$Classifier, method = "sidak")
Dunn_NCBIM2_DS
Dun_NCBIM2_DS.res <- Dunn_NCBIM2_DS$res
Dunncld <- cldList(comparison = Dun_NCBIM2_DS.res$Comparison, p.value = Dun_NCBIM2_DS.res$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "Classifier"
grouped_table_NCBIM2 <- group_by(NCBI_M2, Classifier) %>% summarise(mean=mean(DeltaShannon, na.rm = TRUE), third_quantile = quantile(DeltaShannon, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_NCBIM2 <- left_join(grouped_table_NCBIM2, Dunncld, by = "Classifier")
grouped_table_NCBIM2$MockCommunity <- "Mock2"
print(grouped_table_NCBIM2)

shannon_NCBI_min15 <- filter(shannon_NCBI, MinimumSeqQuality == "Qmin15")

grouped_table_NCBI <- rbind(grouped_table_NCBIM1, grouped_table_NCBIM2)
# Add position for CLD letters based on the position of the error bar for each classifier
grouped_table_NCBI <- grouped_table_NCBI %>%
  left_join(shannon_NCBI_min15 %>% dplyr::select(Classifier, MockCommunity, CI_Upper), 
            by = c("Classifier", "MockCommunity")) %>%
  mutate(label_y = CI_Upper + 0.005)
grouped_table_NCBI$Classifier <- factor(grouped_table_NCBI$Classifier, levels = classifier_order)

# Make sure classifiers are in the right order for plotting
shannon_NCBI_min15$Classifier <- factor(shannon_NCBI_min15$Classifier, levels = classifier_order)

deltashannon_dotplot_CLD2 <- ggplot(shannon_NCBI_min15, aes(x = factor(Classifier), y = Mean, color = Classifier)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  facet_grid2(. ~ MockCommunity, axes = "all") + 
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, linewidth = 0.7) +
  geom_text(data = grouped_table_NCBI, 
            aes(x = Classifier, y = label_y, label = Letter), 
            size = 3.8, 
            vjust = -1, 
            color = "black", 
            inherit.aes = FALSE) +
  labs(subtitle = "NCBI RefSeq ITS Database", 
       y = "Delta Shannon Index") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"), 
        #panel.grid = element_line(color = "lightgrey"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"), 
        legend.position = "none", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.title.y = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(), 
        axis.ticks = element_line(color = "black")) +
  scale_y_continuous(breaks = pretty_breaks(n = 8), 
                     expand = expansion(mult = c(0.05,0.2))) +
  scale_color_manual(values = classifier_colours) +
  scale_x_discrete(labels = c(
    "minimap2" = "minimap2", 
    "kraken2" = "Kraken2", 
    "emu" = "Emu",
    "dnabarcoder" = "Dnabarcoder", 
    "vtd" = "vtd", 
    "mycoAICNN" = "mycoAI CNN", 
    "mycoAIbert" = "mycoAI BERT"))
deltashannon_dotplot_CLD2
ggsave("./DeltaShannonDotplot_NCBI_CLD.png", deltashannon_dotplot_CLD2, bg = "white")

# Repeat for Bhattacharyya distance
#GSD
# Ground Truth is already removed
KW_GSDM1_B <- kruskal.test(formula = Bhattacharyya ~ factor(Classifier), data = GSD_M1)
Dunn_GSDM1_B <- dunnTest(x = GSD_M1$Bhattacharyya, g = GSD_M1$Classifier, method = "sidak")
Dunn_GSDM1_B
Dun_GSDM1_B.res <- Dunn_GSDM1_B$res
Dunncld <- cldList(comparison = Dun_GSDM1_B.res$Comparison, p.value = Dun_GSDM1_B.res$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "Classifier"
grouped_table_GSDM1 <- group_by(GSD_M1, Classifier) %>% summarise(mean=mean(Bhattacharyya, na.rm = TRUE), third_quantile = quantile(Bhattacharyya, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_GSDM1 <- left_join(grouped_table_GSDM1, Dunncld, by = "Classifier")
grouped_table_GSDM1$MockCommunity <- "Mock1"
print(grouped_table_GSDM1)

KW_GSDM2_B <- kruskal.test(formula = Bhattacharyya ~ factor(Classifier), data = GSD_M2)
Dunn_GSDM2_B <- dunnTest(x = GSD_M2$Bhattacharyya, g = GSD_M2$Classifier, method = "sidak")
Dunn_GSDM2_B
Dun_GSDM2_B.res <- Dunn_GSDM2_B$res
Dunncld <- cldList(comparison = Dun_GSDM2_B.res$Comparison, p.value = Dun_GSDM2_B.res$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "Classifier"
grouped_table_GSDM2 <- group_by(GSD_M2, Classifier) %>% summarise(mean=mean(Bhattacharyya, na.rm = TRUE), third_quantile = quantile(Bhattacharyya, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_GSDM2 <- left_join(grouped_table_GSDM2, Dunncld, by = "Classifier")
grouped_table_GSDM2$MockCommunity <- "Mock2"
print(grouped_table_GSDM2)

bhattacharyya_GSD_min15 <- filter(bhattacharyya_GSD, MinimumSeqQuality == "Qmin15")

grouped_table_GSD <- rbind(grouped_table_GSDM1, grouped_table_GSDM2)
# Add position for CLD letters based on the position of the error bar for each classifier
grouped_table_GSD <- grouped_table_GSD %>%
  left_join(bhattacharyya_GSD_min15 %>% dplyr::select(Classifier, MockCommunity, CI_Upper), 
            by = c("Classifier", "MockCommunity")) %>%
  mutate(label_y = CI_Upper + 0.003)

# Make sure classifiers are in the right order for plotting
bhattacharyya_GSD_min15$Classifier <- factor(bhattacharyya_GSD_min15$Classifier, levels = classifier_order)
grouped_table_GSD$Classifier <- factor(grouped_table_GSD$Classifier, levels = classifier_order)

bhattacharyya_dotplot_CLD <- ggplot(bhattacharyya_GSD_min15, aes(x = factor(Classifier), y = Mean, color = Classifier)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  facet_grid2(. ~ MockCommunity, axes = "all") +   
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, linewidth = 0.7) +
  geom_text(data = grouped_table_GSD, 
            aes(x = Classifier, y = label_y, label = Letter), 
            size = 3.8, 
            vjust = -1, 
            color = "black", 
            inherit.aes = FALSE) +
  labs(subtitle = "Gold Standard Database", 
       y = "Bhattacharyya Distance") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"), 
        #panel.grid = element_line(color = "lightgrey"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"), 
        legend.position = "none", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.title.y = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(), 
        axis.ticks = element_line(color = "black")) +
  scale_y_continuous(breaks = pretty_breaks(n = 8), 
                     expand = expansion(mult = c(0.05,0.2))) + 
  scale_color_manual(values = classifier_colours) +
  scale_x_discrete(labels = c(
    "aodp" = "aodp", 
    "minimap2" = "minimap2", 
    "kraken2" = "Kraken2", 
    "emu" = "Emu", 
    "dnabarcoder" = "Dnabarcoder", 
    "vtd" = "vtd", 
    "mycoAICNN" = "mycoAI CNN", 
    "mycoAIbert" = "mycoAI BERT"))
bhattacharyya_dotplot_CLD
ggsave("./BhattacharyyaDotplot_GSD_CLD.png", bhattacharyya_dotplot_CLD, bg = "white", width = 7.29, height = 4.51, units = "in")

#NCBI
KW_NCBIM1_B <- kruskal.test(formula = Bhattacharyya ~ factor(Classifier), data = NCBI_M1)
Dunn_NCBIM1_B <- dunnTest(x = NCBI_M1$Bhattacharyya, g = NCBI_M1$Classifier, method = "sidak")
Dunn_NCBIM1_B
Dun_NCBIM1_B.res <- Dunn_NCBIM1_DS$res
Dunncld <- cldList(comparison = Dun_NCBIM1_B.res$Comparison, p.value = Dun_NCBIM1_B.res$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "Classifier"
grouped_table_NCBIM1 <- group_by(NCBI_M1, Classifier) %>% summarise(mean=mean(Bhattacharyya, na.rm = TRUE), third_quantile = quantile(Bhattacharyya, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_NCBIM1 <- left_join(grouped_table_NCBIM1, Dunncld, by = "Classifier")
grouped_table_NCBIM1$MockCommunity <- "Mock1"
print(grouped_table_NCBIM1)

KW_NCBIM2_B <- kruskal.test(formula = Bhattacharyya ~ factor(Classifier), data = NCBI_M2)
Dunn_NCBIM2_B <- dunnTest(x = NCBI_M2$Bhattacharyya, g = NCBI_M2$Classifier, method = "sidak")
Dunn_NCBIM2_B
Dun_NCBIM2_B.res <- Dunn_NCBIM2_DS$res
Dunncld <- cldList(comparison = Dun_NCBIM2_B.res$Comparison, p.value = Dun_NCBIM2_B.res$P.adj, threshold = 0.05) [1:2]
names(Dunncld)[1] <- "Classifier"
grouped_table_NCBIM2 <- group_by(NCBI_M2, Classifier) %>% summarise(mean=mean(Bhattacharyya, na.rm = TRUE), third_quantile = quantile(Bhattacharyya, 0.75, na.rm = TRUE)) %>% arrange(desc(mean))
grouped_table_NCBIM2 <- left_join(grouped_table_NCBIM2, Dunncld, by = "Classifier")
grouped_table_NCBIM2$MockCommunity <- "Mock2"
print(grouped_table_NCBIM2)

bhattacharyya_NCBI_min15 <- filter(bhattacharyya_NCBI, MinimumSeqQuality == "Qmin15")

grouped_table_NCBI <- rbind(grouped_table_NCBIM1, grouped_table_NCBIM2)
grouped_table_NCBI <- filter(grouped_table_NCBI, Classifier != "GroundTruth")
# Add position for CLD letters based on the position of the error bar for each classifier
grouped_table_NCBI <- grouped_table_NCBI %>%
  left_join(bhattacharyya_NCBI_min15 %>% dplyr::select(Classifier, MockCommunity, CI_Upper), 
            by = c("Classifier", "MockCommunity")) %>%
  mutate(label_y = CI_Upper + 0.005)

# Make sure classifiers are in the right order for plotting
bhattacharyya_NCBI_min15$Classifier <- factor(bhattacharyya_NCBI_min15$Classifier, levels = classifier_order)
grouped_table_NCBI$Classifier <- factor(grouped_table_NCBI$Classifier, levels = classifier_order)

bhattacharyya_dotplot_CLD2 <- ggplot(bhattacharyya_NCBI_min15, aes(x = factor(Classifier), y = Mean, color = Classifier)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  facet_grid2(. ~ MockCommunity, axes = "all") + 
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, linewidth = 0.7) +
  geom_text(data = grouped_table_NCBI, 
            aes(x = Classifier, y = label_y, label = Letter), 
            size = 3.8, 
            vjust = -1, 
            color = "black") +
  labs(subtitle = "NCBI RefSeq ITS Database",
       y = "Bhattacharyya Distance") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"), 
        #panel.grid = element_line(color = "lightgrey"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"), 
        legend.position = "none", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.title.y = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(), 
        axis.ticks = element_line(color = "black")) +
  scale_y_continuous(breaks = pretty_breaks(n = 8), 
                     expand = expansion(mult = c(0.05,0.2))) + 
  scale_color_manual(values = classifier_colours) +
  scale_x_discrete(labels = c(
    "aodp" = "aodp", 
    "minimap2" = "minimap2", 
    "kraken2" = "Kraken2", 
    "emu" = "Emu", 
    "dnabarcoder" = "Dnabarcoder", 
    "vtd" = "vtd", 
    "mycoAICNN" = "mycoAI CNN", 
    "mycoAIbert" = "mycoAI BERT"))
bhattacharyya_dotplot_CLD2
ggsave("./BhattacharyyaDotplot_NCBI_CLD.png", bhattacharyya_dotplot_CLD2, bg = "white", width = 7.29, height = 4.51, units = "in")

# Save results for species not in NCBI for later analyses
# Get current date in YYYYMMDD format
date_str <- format(Sys.Date(), "%Y%m%d")

# Loop through the list and save each data frame
for (name in names(Mock1_notinNCBI)) {
  filename <- paste0(date_str, "_", name, ".csv")
  write.csv(Mock1_notinNCBI[[name]], file = paste0("./UnseenSpecies/", filename), row.names = FALSE)
}

for (name in names(Mock2_notinNCBI)) {
  filename <- paste0(date_str, "_", name, ".csv")
  write.csv(Mock2_notinNCBI[[name]], file = paste0("./UnseenSpecies/", filename), row.names = FALSE)
}


