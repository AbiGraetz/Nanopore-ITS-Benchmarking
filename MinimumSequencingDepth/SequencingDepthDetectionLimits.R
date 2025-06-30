library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)

# Read in the mock community classification results for the Gold Standard Database(GSD)
M1_list <- list.files("./Mock1_Results/GSD/", pattern = ".xlsx", full.names = TRUE)
M2_list <- list.files("./Mock2_Results/GSD/", pattern = ".xlsx", full.names = TRUE)

M1_GSD <- list()
M2_GSD <- list()

for (i in seq_along(M1_list)) {
  file = M1_list[[i]]
  df <- read_xlsx(file, col_names = TRUE)
  df_name <- tools::file_path_sans_ext(basename(file))
  M1_GSD[[df_name]] <- df
}

for (i in seq_along(M2_list)) {
  file = M2_list[[i]]
  df <- read_xlsx(file, col_names = TRUE)
  df$`mapped to this species` <- round(df$`mapped to this species`)
  df_name <- tools::file_path_sans_ext(basename(file))
  M2_GSD[[df_name]] <- df
}

# Read in the ground truth distribution for each Mock community, for comparison
GSD_list <- list.files('./GroundTruth/', pattern = ".csv", full.names = TRUE)
GT <- list()
for (i in seq_along(GSD_list)) {
  GT[[i]] <- read.csv(GSD_list[[i]], header = TRUE)
  GT[[i]] <- dplyr::select(GT[[i]], c("organism", "sequence_abundance"))
}
names(GT) <- c("Mock1_GT", "Mock2_GT")

K2_list <- list()
K2_list[["Mock1_GSD_kraken2_Qmin15"]] <- M1_GSD[["Mock1_GSD_kraken2_Qmin15"]]
K2_list[["Mock2_GSD_kraken2_Qmin15"]] <- M2_GSD[["Mock2_GSD_kraken2_Qmin15"]]
# Because Kraken2 only identified 53 species, running an adjusted version of the function
# Remove Botrytis cinerea, for which there should be zero counts. Check this first.
botrytis_counts <- sapply(K2_list, function(df) {
  val <- df %>%
    filter(species == "Botrytis_cinerea") %>%
    pull(`mapped to this species`)
  if (length(val) == 0) NA else val
})
print(botrytis_counts)
if (all(botrytis_counts == 0, na.rm = TRUE)) {
  K2_list <- lapply(K2_list, function(df) {
    df[df$species != "Botrytis_cinerea", ]
  })
}
# All data frames in K2_list should now have 53 rows. 

# Remove existing Kraken2 data from M1_GSD and M2_GSD, it's unneccessary
M1_GSD <- M1_GSD[names(M1_GSD) != "Mock1_GSD_kraken2_Qmin15"]
M2_GSD <- M2_GSD[names(M2_GSD) != "Mock2_GSD_kraken2_Qmin15"]

# To use the function, input must be a named vector. Initialise empty lists to store these. 
abundances_M1 <- list()
for (i in seq_along(M1_GSD)) {
  df <- M1_GSD[[i]]
  abundances <- setNames(df$`mapped to this species`, df$species)
  
  #Calculate number of unclassified sequences
  unclassified_seqs <- sum(GT$Mock1_GT$sequence_abundance) - sum(df$`mapped to this species`)
  
  # Check for negative unclassified sequences 
  if (unclassified_seqs < 0) {
    warning(paste("Classifier", i, ": Total input sequences less than sum of classified sequences. Setting unclassified to 0."))
    unclassified_seqs <- 0
  }
  
  # Append unclassified sequences to named vector
  abundances <- c(abundances, setNames(unclassified_seqs, "unclassified"))
  
  # Store the result
  abundances_M1[[i]] <- abundances
}
names(abundances_M1) <- names(M1_GSD)

abundances_M2 <- list()
for (i in seq_along(M2_GSD)) {
  df <- M2_GSD[[i]]
  abundances <- setNames(df$`mapped to this species`, df$species)
  unclassified_seqs <- sum(GT$Mock2_GT$sequence_abundance) - sum(df$`mapped to this species`)
  if (unclassified_seqs < 0) {
    warning(paste("Classifier", i, ": Total input sequences less than sum of classified sequences. Setting unclassified to 0."))
    unclassified_seqs <- 0
  }
  abundances <- c(abundances, setNames(unclassified_seqs, "unclassified"))
  abundances_M2[[i]] <- abundances
}
names(abundances_M2) <- names(M2_GSD)

abundances_K2 <- list()
abundances <- setNames(K2_list[[1]]$`mapped to this species`, K2_list[[1]]$species)
unclassified_seqs <- sum(GT$Mock1_GT$sequence_abundance) - sum(K2_list[[1]]$`mapped to this species`)
abundances <- c(abundances, unclassified = unclassified_seqs)
abundances_K2[[1]] <- setNames(as.integer(abundances), names(abundances)) # Coerce to integer, default is double

abundances <- setNames(K2_list[[2]]$`mapped to this species`, K2_list[[2]]$species)
unclassified_seqs <- sum(GT$Mock2_GT$sequence_abundance) - sum(K2_list[[2]]$`mapped to this species`)
abundances <- c(abundances, unclassified = unclassified_seqs)
abundances_K2[[2]] <- setNames(as.integer(abundances), names(abundances)) # Coerce to integer, default is double
names(abundances_K2) <- names(K2_list) # Name list elements

abundances_GT <- list()
for (i in seq_along(GT)) {
  abundances_GT[[i]] <- setNames(GT[[i]]$sequence_abundance, GT[[i]]$organism)
}
names(abundances_GT) <- names(GT)
abundances_GT <- abundances_GT[names(abundances_GT) != "Mock3_GT"]

# Function to find minimum sequencing depth for a given classifier's output
min_seq_depth <- function(abundances, prob_threshold = 0.95, num_sim = 500, num_iter = 100, base_seed = 123) {
  # abundances: named vector of sequence counts with species names and optionally "unclassified"
  # prob_threshold: desired probability (e.g., 0.95)
  # num_sim: number of simulations per n
  # num_iter: number of iterations to average over (default 100)
  # base_seed: base seed for reproducibility (default 123)
  
  if (!interactive()) {
    stop("This function requires an interactive session. Please try again in an R or RStudio session.")
  }
  
  # Check if abundances is a named vector
  if (is.null(names(abundances)) || any(names(abundances) == "")) {
    stop("Abundances must be a named vector with species names.")
  }
  
  # Identify species to detect (exclude "unclassified")
  species_names <- names(abundances)[names(abundances) != "unclassified"]
  num_species <- length(species_names)
  
  # Interactive check for number of species
  cat("Detected", num_species, "species (excluding 'unclassified').\n")
  answer <- tolower(readline(prompt = "Is this the expected number of species to detect? (y/n): "))
  if (answer == "n") {
    cat("Exiting function. Please check your input data.\n")
    return(invisible(NULL))
  } else if (answer != "y") {
    cat("Invalid input. Please enter 'y' or 'n'.\n")
    return(invisible(NULL))
  }
  
  cat("Great! Continuing with function...\n")
  
  # Total sequences
  N <- sum(abundances)
  
  # Vector of species labels (including "unclassified" if present)
  species <- rep(names(abundances), times = abundances)
  
  # Function to estimate probability for a given n
  estimate_prob <- function(n, condition) {
    results <- replicate(num_sim, {
      subsample <- sample(species, size = n, replace = FALSE)
      tbl <- table(subsample)
      # Initialize counts for all species (including "unclassified")
      counts <- numeric(length(abundances))
      names(counts) <- names(abundances)
      counts[names(tbl)] <- tbl
      # Subset to species of interest (exclude "unclassified")
      species_counts <- counts[species_names]
      if (condition == "detect_all") {
        return(length(species_counts[species_counts > 0]) == num_species)  # True if all species detected
      } else if (condition == "min_abund") {
        min_count <- min(species_counts)
        return(min_count >= 109)  # True if all species have >= 109 sequences
      }
    })
    mean(results)  # Estimated probability
  }
  
  # Binary search for minimum n
  binary_search <- function(condition, n_low = 1, n_high = N) {
    while (n_low < n_high) {
      n_mid <- floor((n_low + n_high) / 2)
      prob <- estimate_prob(n_mid, condition)
      if (prob >= prob_threshold) {
        n_high <- n_mid
      } else {
        n_low <- n_mid + 1
      }
    }
    return(n_low)
  }
  
  # Run iterations with different seeds
  results_detect_all <- numeric(num_iter)
  results_min_abund <- numeric(num_iter)
  
  for (i in 1:num_iter) {
    cat("Iteration", i, "of", num_iter, "\n")
    set.seed(base_seed + i - 1)
    results_detect_all[i] <- binary_search("detect_all")
    results_min_abund[i] <- binary_search("min_abund")
  }
  
  # Return mean results across iterations
  return(list(
    depth_detect_all = mean(results_detect_all),
    depth_min_abund = mean(results_min_abund)
  ))
}

# Apply the function
result_list <- list()
for (i in seq_along(abundances_M1)) {
  result_list[[i]] <- min_seq_depth(abundances_M1[[i]])
}
names(result_list) <- names(abundances_M1)

result_list2 <- list()
for (i in seq_along(abundances_M2)) {
  result_list2[[i]] <- min_seq_depth(abundances_M2[[i]])
}
names(result_list2) <- names(abundances_M2)

result_listK2 <- list() 
for (i in seq_along(abundances_K2)) {
  result_listK2[[i]] <- min_seq_depth(abundances_K2[[i]])
}
names(result_listK2) <- names(K2_list)

result_listGT <- list()
for (i in seq_along(abundances_GT)) {
  result_listGT[[i]] <- min_seq_depth(abundances_GT[[i]])
}
names(result_listGT) <- names(abundances_GT)

# Get the data out of the result list and into a data frame
data_M1 <- result_list %>% 
  enframe(name = "Classifier", value = "Depths") %>%
  unnest_wider(Depths) %>%
  rename(Detection = depth_detect_all, MinimumAbundance = depth_min_abund)
data_M1 <- separate_wider_delim(data_M1, cols = Classifier, delim = "_", names = c("MockCommunity", "Database", "Classifier", "MinimumSeqQuality")) %>%
  dplyr::select(c("Classifier", "MockCommunity", "Database", "MinimumSeqQuality", "Detection", "MinimumAbundance")) %>%
  filter(MinimumSeqQuality == "Qmin15")

data_M2 <- result_list2 %>% 
  enframe(name = "Classifier", value = "Depths") %>%
  unnest_wider(Depths) %>%
  rename(Detection = depth_detect_all, MinimumAbundance = depth_min_abund)
data_M2 <- separate_wider_delim(data_M2, cols = Classifier, delim = "_", names = c("MockCommunity", "Database", "Classifier", "MinimumSeqQuality")) %>%
  dplyr::select(c("Classifier", "MockCommunity", "Database", "MinimumSeqQuality", "Detection", "MinimumAbundance")) %>%
  filter(MinimumSeqQuality == "Qmin15")

data_GT <- result_listGT %>% 
  enframe(name = "Classifier", value = "Depths") %>%
  unnest_wider(Depths) %>%
  rename(Detection = depth_detect_all, MinimumAbundance = depth_min_abund)
# Add columns to data_GT so that this data matches the experimental data, and can be row bound to the appropriate results data frame.
data_GT$Classifier <- "Ideal Classifier"
data_GT$MockCommunity <- c("Mock1", "Mock2")
data_GT$Database <- "GSD"
data_GT$MinimumSeqQuality <- "Qmin15"
data_GT <- dplyr::select(data_GT, c("Classifier", "MockCommunity", "Database", "MinimumSeqQuality", "Detection", "MinimumAbundance"))

data_K2 <- result_listK2 %>%
  enframe(name = "Classifier", value = "Depths") %>%
  unnest_wider(Depths) %>%
  rename(Detection = depth_detect_all, MinimumAbundance = depth_min_abund)
data_K2 <- separate_wider_delim(data_K2, cols = Classifier, delim = "_", names = c("MockCommunity", "Database", "Classifier", "MinimumSeqQuality")) %>%
  dplyr::select(c("Classifier", "MockCommunity", "Database", "MinimumSeqQuality", "Detection", "MinimumAbundance"))

# Bind the appropriate GT data to the mock community data. Remove Kraken2 data calculated with 54 species and bind data calculated with 53 species
data_M1 <- rbind(data_M1, data_GT[1,]) %>% rbind(data_K2[1,])
data_M2 <- rbind(data_M2, data_GT[2,]) %>% rbind(data_K2[2,])

library(ggforce)
library(colorspace)

M1_long <- data_M1 %>%
  pivot_longer(
    cols = c(Detection, MinimumAbundance), 
    names_to = "Threshold",
    values_to = "Depth"
  ) %>%
  mutate(
    Panel = case_when(
      Depth <= 500 ~ "Low",
      Threshold == "MinimumAbundance" ~ "High"
    )
  )
# Calculate threshold means to add to the plot 
means_M1 <- M1_long %>% 
  group_by(Threshold) %>%
  summarise(mean_depth = mean(Depth, na.rm = TRUE)) %>%
  mutate(Threshold = factor(Threshold, levels = c("Detection", "MinimumAbundance")))

#Define the colours for points and mean lines
point_colours <- c("Detection" = "#009688", "MinimumAbundance" = "#a020f0") # Second hex code is for R default "purple", used col2rgb() to find RGB values and converted to hex.
line_colours <- c("Detection" = lighten("#009688", amount = 0.4), 
                  "MinimumAbundance" = lighten("#a020f0"), amount = 0.4)
all_colours <- c(point_colours, line_colours)

library(ggbreak)

# Define the order of classifiers and re-name so that everything is consistent
classifier_order <- c("Ideal Classifier", "MycoAIBert", "MycoAICNN", "vtd", "Dnabarcoder", "emu", "kraken2", "minimap2", "aodp")
M1_long$Classifier <- factor(M1_long$Classifier, levels = classifier_order)

M1_dumbbell <- ggplot(M1_long, aes(y = Classifier, x = Depth)) +
  # Plot points for each threshold
  geom_point(aes(color = Threshold), size = 6) + 
  # Plot segments to connect points and make the dumbbell
  geom_segment(data = data_M1, 
               aes(x = Detection, xend = MinimumAbundance, 
                   y = Classifier, yend = Classifier), 
               color = "#b2b2b2", size = 1) +
  #scale_x_continuous(minor_breaks = seq(0, 550, by = 25)) +
  #scale_x_break(breaks = c(1000, 5000), # Create broken x-axis
  #              scale = c(0.8, 0.2), # Scale each side of the break so one isn't squished
  #              space = 0.1) +
  # Add vertical line for means, colour by threshold
  geom_vline(data = means_M1, 
             aes(xintercept = mean_depth, color = paste0(Threshold, "_mean")), 
             linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
  scale_color_manual(values = c("Detection" = "#009688", 
                                "MinimumAbundance" = "#a020f0",
                                "Detection_mean" = lighten("#009688", amount = 0.4), 
                                "MinimumAbundance_mean" = lighten("#a020f0", amount = 0.4)), 
                     breaks = c("Detection", "MinimumAbundance"), 
                     name = "Threshold") +
  labs(title = "Mock1", 
       x = "Sequencing Depth", 
       y = "Classifier") + 
  theme_minimal() + 
  scale_y_discrete(labels = c(
    "aodp" = "aodp", 
    "minimap2" = "minimap2", 
    "kraken2" = "Kraken2", 
    "emu" = "Emu", 
    "Dnabarcoder" = "Dnabarcoder", 
    "vtd" = "vtd", 
    "MycoAICNN" = "mycoAI CNN", 
    "MycoAIBert" = "mycoAI BERT")) +
  scale_x_continuous(minor_breaks = seq(0, 13000, by = 500)) +
  theme(axis.text.y = element_text(size = 16, color = "black"), 
        axis.text.x = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = "bottom", 
        legend.margin = margin(t = 0, b = 0), 
        legend.text = element_text(size = 16, color = "black"),
        axis.text.y.right = element_blank(),
        legend.title = element_blank(), 
        axis.ticks.x.bottom = element_line(color = "black"), 
        axis.ticks.y.left = element_line(color = "black"))
M1_dumbbell #+ xlim(200, 15000)

# Break up the x-axis for Mock2, because the minimum abundance data is much higher value than the detection data
M2_long <- data_M2 %>%
  filter(Classifier != "emumix") %>%
  pivot_longer(
    cols = c(Detection, MinimumAbundance), 
    names_to = "Threshold",
    values_to = "Depth"
  ) %>%
  mutate(
    Panel = case_when(
      Depth <= 700 ~ "Low", # Manually inspect the data to find the right value to create the break at
      TRUE ~ "High"
    )
  )
M2_long$Classifier <- factor(M2_long$Classifier, levels = classifier_order)
# Calculate threshold means to add to the plot 
means_M2 <- M2_long %>% 
  group_by(Threshold) %>%
  summarise(mean_depth = mean(Depth, na.rm = TRUE)) %>%
  mutate(Threshold = factor(Threshold, levels = c("Detection", "MinimumAbundance")))

M2_dumbbell <- ggplot(M2_long, aes(y = Classifier, x = Depth)) +
  # Plot points for each threshold
  geom_point(aes(color = Threshold), size = 6) + 
  # Plot secments to connnect points and make the dumbbell
  geom_segment(data = data_M2, 
               aes(x = Detection, xend = MinimumAbundance, 
                   y = Classifier, yend = Classifier), 
               color = "#b2b2b2", size = 1) +
  # Add vertical line for means, colour by threshold
  geom_vline(data = means_M2, 
             aes(xintercept = mean_depth, color = paste0(Threshold, "_mean")), 
             linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
  # Apply broken x-axis, and scale the segment widths
  #scale_x_break(breaks = c(1000, 10000), # Create broken x-axis
  #              scale = c(0.8, 0.2), # Scale each side of the break so one isn't squished
  #              space = 0.05) + # Increase the amount of space the break takes, so x-axis text doesn't overlap
  scale_color_manual(values = c("Detection" = "#009688", 
                                "MinimumAbundance" = "#a020f0",
                                "Detection_mean" = lighten("#009688", amount = 0.4), 
                                "MinimumAbundance_mean" = lighten("#a020f0", amount = 0.4)), 
                     breaks = c("Detection", "MinimumAbundance"), 
                     name = "Threshold") +
  labs(title = "Mock2", 
       x = "Sequencing Depth", 
       y = "Classifier") + 
  theme_minimal() + 
  scale_y_discrete(labels = c(
    "aodp" = "aodp", 
    "minimap2" = "minimap2", 
    "kraken2" = "Kraken2", 
    "emu" = "Emu", 
    "Dnabarcoder" = "Dnabarcoder", 
    "vtd" = "vtd", 
    "MycoAICNN" = "mycoAI CNN", 
    "MycoAIBert" = "mycoAI BERT")) +
  scale_x_continuous(minor_breaks = seq(0, 21000, by = 1000)) +
  theme(axis.text.y = element_text(size = 16, color = "black"), 
        axis.text.x = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 18, color = "black"),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = "bottom", 
        legend.margin = margin(t = 0, b = 0), 
        legend.text = element_text(size = 16, color = "black"),
        axis.text.y.right = element_blank(),
        legend.title = element_blank(), 
        axis.ticks.x.bottom = element_line(color = "black"), 
        axis.ticks.y.left = element_line(color = "black")) # Remove redundant y-axis labels on the right hand side
M2_dumbbell

