# Clear environment and load necessary libraries
rm(list = ls())
library(seqinr)
library(tidyverse)

# Load pruned FASTA files
H2A_Pruned_wSQ <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A/PRUNED/Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_OG2.fa", seqtype = "AA")
H2A_Pruned_woSQ <- read.fasta("/Users/perrinm/Documents/AHocher_Project/TREE/H2A_Minus_SQ_Minus_H2AZ_2.fasta", seqtype = "AA")


getName(H2A_Pruned_wSQ)[getName(H2A_Pruned_wSQ) %in% getName(H2A_Pruned_woSQ)]

table(getName(H2A_Pruned_wSQ) %in% getName(H2A_Pruned_woSQ))

# Identify sequences that appear in H2A_Pruned_wSQ but not in H2A_Pruned_woSQ
tree_diff <- getName(H2A_Pruned_wSQ)[!getName(H2A_Pruned_wSQ) %in% getName(H2A_Pruned_woSQ)]
tree_diff <- data.frame(species_ID = tree_diff)

# Extract clean species ID for tree_diff
tree_diff$species_ID_Clean <- gsub(".*\\$(.*?)@.*", "\\1", tree_diff$species_ID)

# Extract clean species ID for H2A_Pruned_wSQ and H2A_Pruned_woSQ
H2A_Pruned_wSQ_df <- data.frame(species_ID = getName(H2A_Pruned_wSQ))
H2A_Pruned_woSQ_df <- data.frame(species_ID = getName(H2A_Pruned_woSQ))
H2A_Pruned_wSQ_df$species_ID_Clean <- gsub(".*\\$(.*?)@.*", "\\1", H2A_Pruned_wSQ_df$species_ID)
H2A_Pruned_woSQ_df$species_ID_Clean <- gsub(".*\\$(.*?)@.*", "\\1", H2A_Pruned_woSQ_df$species_ID)

# Define the "SQ" motif pattern
sq_pattern <- "sq(e|d)|SQ(E|D)"

# Create a function to check if a sequence contains the "SQ" motif
contains_sq <- function(seq) {
  grepl(sq_pattern, paste(seq, collapse = ""))
}

# Identify sequences in H2A_Pruned_wSQ that contain the SQ motif
H2A_Pruned_wSQ_df$contains_SQ <- sapply(H2A_Pruned_wSQ, contains_sq)
H2A_Pruned_wSQ_df$contains_SQ <- as.integer(H2A_Pruned_wSQ_df$contains_SQ)  # Convert to 1 (TRUE) / 0 (FALSE)

# Identify sequences in H2A_Pruned_woSQ that contain the SQ motif
H2A_Pruned_woSQ_df$contains_SQ <- sapply(H2A_Pruned_woSQ, contains_sq)
H2A_Pruned_woSQ_df$contains_SQ <- as.integer(H2A_Pruned_woSQ_df$contains_SQ)  # Convert to 1 (TRUE) / 0 (FALSE)

# Merge the SQ motif presence info into tree_diff based on species_ID_Clean
tree_diff <- tree_diff %>%
  left_join(H2A_Pruned_wSQ_df %>% select(species_ID_Clean, contains_SQ_wSQ = contains_SQ), by = "species_ID_Clean") %>%
  left_join(H2A_Pruned_woSQ_df %>% select(species_ID_Clean, contains_SQ_woSQ = contains_SQ), by = "species_ID_Clean")

# Display the final table with SQ motif information
print(tree_diff)

tree_diff <- drop_na(tree_diff)
