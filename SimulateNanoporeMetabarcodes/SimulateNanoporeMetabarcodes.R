library(dplyr)
library(data.table)
library(tidyr)
library(Biostrings)
library(ggplot2)
library(dyngen)
library(seqinr)
library(readxl)

mocktaxonomy <- read_xlsx("./SpeciesInNCBI_Taxonomy.xlsx", col_names = TRUE)

#Reading in the sequences also as a data frame, because this is much easier to manipulate in R than a DNAStringSet object.
seqs <- readLines("./ConsensusAlignments_NCBIRefSeq_LongITS.fasta")
headers <- character()
sequences <- character()

for (line in seqs) {
  if (startsWith(line, ">")) {
    #This line assigns the headers based on the ">" character at the start of the heading of all FASTA format files
    headers <- c(headers, substr(line, 2, nchar(line)))
    sequences <- c(sequences, "")
  }
  else {
    #This line assigns the sequences
    sequences[length(sequences)] <- paste0(sequences[length(sequences)], line)
  }
}

seqsDF <- data.frame(accession = headers, sequence = sequences, stringsAsFactors = FALSE)

seqsDF <- separate_wider_delim(seqsDF, accession, delim = "_", names = c("accession", "header"), too_many = "merge") %>% select(c("accession", "sequence"))
seqsDF$sequence <- gsub("\\?", "N", seqsDF$sequence) # Remove the ? characters introduced by Geneious prime to bridge assembly gaps.

trainingset <- list()

replicate_rows <- function(df, n) {
  df_rep <- df[rep(seq_len(nrow(df)), each = n), ]
  return(df_rep)
}
trainingset <- replicate_rows(seqsDF, 2200)

accession_counts <- trainingset %>% group_by(accession) %>% summarise(count = n())
all_expected_counts <- all(accession_counts$count == 2200)
unique_accession_counts <- trainingset %>% distinct(accession) %>% nrow()

if(all_expected_counts) {
  print("all accessions have the expected count of 2200")
} else {
  print("some accessions do not have the expected count of 2200")
}

if(unique_accession_counts == 35) {
  print("there are 35 unique accessions") 
} else {
  print("there are not 35 unique accessions")
}

# Split the training set by accession
trainingset_mut <- split(trainingset, trainingset$accession)

for (n in seq_along(trainingset_mut)) {
  trainingset_mut[[n]]$seqlen <- nchar(trainingset_mut[[n]]$sequence) #get the length of each sequence and assign value to a column. Useful for checking that indel mutations have been produced.
}

# Define a function to introduce mutations
introduce_indels <- function(df, seed) {
  set.seed(seed)  # Set a unique seed for reproducibility per data frame
  
  # Randomly select 1% of rows
  sample_rows <- sample(nrow(df), size = round(0.01 * nrow(df)))
  
  # Generate consistent insertion and deletion positions for the entire data frame
  set.seed(seed)  # Use the same seed for reproducibility
  insertion_positions <- sample(1:(nchar(df$sequence[1]) - 1), 2)
  
  set.seed(seed + 1)  # Change the seed slightly to generate new randomness
  insertion_letters <- sample(c('A', 'C', 'G', 'T'), 2, replace = TRUE)
  
  set.seed(seed + 2)  # Change the seed again for deletion positions
  deletion_positions <- sample(1:(nchar(df$sequence[1])), 4)
  
  # Add mutated_sequence column
  df$mutated_sequence <- df$sequence
  
  # Apply the consistent InDels to all selected rows
  for (i in sample_rows) {
    # Apply the same insertions
    for (j in seq_along(insertion_positions)) {
      df$mutated_sequence[i] <- paste0(
        substr(df$mutated_sequence[i], 1, insertion_positions[j]),
        insertion_letters[j],
        substr(df$mutated_sequence[i], insertion_positions[j] + 1, nchar(df$mutated_sequence[i]))
      )
    }
    
    # Apply the same deletions
    for (j in seq_along(deletion_positions)) {
      df$mutated_sequence[i] <- paste0(
        substr(df$mutated_sequence[i], 1, deletion_positions[j] - 1),
        substr(df$mutated_sequence[i], deletion_positions[j] + 1, nchar(df$mutated_sequence[i]))
      )
    }
  }
  
  # Add a column for mutated sequence length
  df$mutseqlen <- nchar(df$mutated_sequence)
  
  return(df)
}

# Iterate over each data frame in the list and apply the indel function
for (i in seq_along(trainingset_mut)) {
  trainingset_mut[[i]] <- introduce_indels(trainingset_mut[[i]], seed = 123 + i)  # Unique seed per data frame
}

# Re-bind data frames if needed
trainingset_mutated <- do.call(rbind, trainingset_mut)

# Automated sanity check function
sanity_check <- function(data_list) {
  for (i in seq_along(data_list)) {
    cat("\nChecking data frame:", i, "\n")
    
    df <- data_list[[i]]
    
    # Check if mutated rows are consistent with the expected changes
    mutated_rows <- df$sequence != df$mutated_sequence
    actual_length_change <- df$mutseqlen - nchar(df$sequence)
    expected_length_change <- 2 - 4
    
    # Check length changes
    if (any(actual_length_change[mutated_rows] != expected_length_change)) {
      warning(paste("Length change mismatch detected in mutated rows for data frame", i))
    } else {
      print("All length changes are as expected.")
    }
    
    # Check if all mutated sequences are identical within the same data frame
    mutated_sequences <- unique(df$mutated_sequence[mutated_rows])
    mutated_lengths <- unique(df$mutseqlen[mutated_rows])
    
    if (length(mutated_sequences) > 1) {
      warning(paste("Inconsistent mutated sequences detected in data frame", i))
      cat("Mutated sequences found:\n")
      print(mutated_sequences)
    } else {
      print("All mutated sequences are identical within this data frame.")
    }
    
    if (length(mutated_lengths) > 1) {
      warning(paste("Inconsistent mutated sequence lengths detected in data frame", i))
      cat("Mutated sequence lengths found:\n")
      print(mutated_lengths)
    } else {
      print("All mutated sequence lengths are identical within this data frame.")
    }
  }
  
  cat("\nSanity check complete!\n")
}

# Apply sanity check to the list of data frames
sanity_check(trainingset_mut)

seqerror_values <- seq(from = 0, to = 3.16, length.out = 2200) # Define the values to be randomly sampled.  

set.seed(123) # Set seed for reproducibility
seqerror_dist <- sample(seqerror_values, 2200, replace = TRUE) # Sample values randomly with replacement to create the distribution. Length = number of sequences/species.
seqerror_dist <- rnorm_bounded(seqerror_dist, mean = 1.58, sd = 2, min = 0, max = 3.16) # Redefine the distribution of values to a bounded normal distribution. max = maximum allowed error. 
seqerror_dist <- as.data.frame(seqerror_dist)
colnames(seqerror_dist) <- "percerror_value"

# Plot the distribution of error values to visualise. All values are percentages of sequence length.
seqerror_densityplot <- ggplot(data = seqerror_dist, aes(x = percerror_value)) +
  geom_density() +
  labs(title = "Distribution of sequencing error values for NCBI Training Set", x = "Sequencing error (%)", y = "Frequency") +
  theme_minimal()
seqerror_densityplot

#Iterate over each data frame, and calculate the length of each sequence. Append the length and percentage error columns to the sequences, and use these values to calculate how many bases should be altered in each sequence. 
trainingset_fin <- split(trainingset_mutated, trainingset_mutated$accession)
for (i in seq_along(trainingset_fin)) {
  df <- trainingset_fin[[i]]
  df$percerror <- seqerror_dist$percerror_value
  df$basestomutate <- round((df$seqlen * (df$percerror/100)), digits = 0) #use 'round()' to make sure the value is an integer
  trainingset_fin[[i]] <- df
}
print(trainingset_fin[[1]])

#Use the newly calculated 'basestomutate' value to mutate the sequences.
for (i in seq_along(trainingset_fin)) {
  df <- trainingset_fin[[i]]
  df$final_sequence = NA
  df$final_seqlen = NA
  
  #define a function to mutate a given DNA sequence
  mutate_sequence <- function(sequence, positions) {
    iupac_codes <- c("A", "C", "G", "T") # Define the pool of bases from which R can draw the mutant base
    set.seed(234)  # Set seed for reproducibility
    new_bases <- sample(iupac_codes, length(positions), replace = TRUE)  #randomly select new bases for mutation
    
    sequence_mutated <- strsplit(sequence, "")[[1]] #split the string into distinct letters, so each position can be accessed
    sequence_mutated[positions] <- new_bases #replace bases at selected positions with mutant bases
    
    final_sequence <- paste(sequence_mutated, collapse = "") #join the mutated bases back into a string, collapse white space
    return(final_sequence)
  }
  
  #iterate over each row in the current data frame df and mutate the sequence
  for (j in 1:nrow(df)) {
    # Extract information from the current row
    sequence <- df$mutated_sequence[j]
    basestomutate <- df$basestomutate[j]
    set.seed(234)
    positions_to_mutate <- sample(1:nchar(sequence), basestomutate) # Randomly select positions to mutate, size = number of mutations
    
    final_sequence <- mutate_sequence(sequence, positions_to_mutate) # Mutate the sequence using the defined function
    
    df$final_sequence[j] <- final_sequence # Add the mutated sequence to the new column created in the first steps of the loop
    df$final_seqlen[j] <- nchar(df$final_sequence[j])
  }
  
  # Print the updated data frame, put the i'th dataframe back into position
  print(df)
  trainingset_fin[[i]] <- df
}

# Re-bind accession-split data frames back together
trainingset_final <- do.call(rbind, trainingset_fin)

trainingset_final$fasta_header <- paste(trainingset_final$accession, "_NCBIRefSeq_AssembledConsensus")

write.fasta(as.list(trainingset_final$final_sequence), names = trainingset_final$fasta_header, file.out = "./simNCBI_TrainingSet.fasta", open = "w")



