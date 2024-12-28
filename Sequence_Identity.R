#This script is to calculate the 100%, between all the reconstruced Seqs
rm(list = ls())
library(seqinr)
library(tidyverse)
#Set WD and import fasta
setwd("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/RECONSTRUCTION")

H2ASequences <- read.fasta("Eukaryota_Reconstructions.fasta", seqtype = "AA")

seq_list <- lapply(H2ASequences, as.character)

# Function to calculate sequence identity between two sequences
calculate_identity <- function(seq1, seq2) {
  if (length(seq1) != length(seq2)) {
    stop("Sequences must be of the same length (aligned).")
  }
  matches <- sum(seq1 == seq2)
  total_positions <- length(seq1)
  identity <- (matches / total_positions) * 100
  return(identity)
}

# Get the number of sequences
num_sequences <- length(seq_list)

# Initialize a matrix to store sequence identities
identity_matrix <- matrix(NA, nrow = num_sequences, ncol = num_sequences,
                          dimnames = list(names(H2ASequences), names(H2ASequences)))

# Loop through all pairs of sequences
for (i in 1:num_sequences) {
  for (j in 1:num_sequences) {
    if (i <= j) {  # Only compute upper triangle (or diagonal)
      identity_matrix[i, j] <- calculate_identity(seq_list[[i]], seq_list[[j]])
    } else {  # Copy values to lower triangle for symmetry
      identity_matrix[i, j] <- identity_matrix[j, i]
    }
  }
}
#Save file as .csv
setwd("/Users/perrinm/Documents/AHocher_Project/ANALYSIS")
write.csv(identity_matrix ,file = "Sequence_Identity_Matrix.csv")

#generate heatmap for 
library(pheatmap)
pheatmap(identity_matrix,
         cluster_rows = FALSE,  # Disable clustering for rows
         cluster_cols = FALSE,  # Disable clustering for columns
         display_numbers = TRUE,  # Show sequence identity values on the heatmap
         color = colorRampPalette(c("grey80", "steelblue1"))(50),  # Customize color gradient
         main = "Sequence Identity Heatmap")  # Add a title
