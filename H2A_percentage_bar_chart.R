#This is a script to produce a bar graph of the percentage of H2AX containing H2As in each species list
#Clear environment and load libraries
rm(list = ls())

# Set working directory
setwd("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/FINAL")

# Load required libraries
library(seqinr)   # For read.fasta
library(tidyverse)    # For data manipulation
library(ggpubr)
# List all .fasta files in the directory
fasta_files <- list.files(pattern = "\\.fasta$")

# Initialize an empty data frame to store counts
results <- data.frame(File = character(), Count = integer(), stringsAsFactors = FALSE)

# Loop through each .fasta file to count sequences containing SQ patterns
for (file in fasta_files) {
  
  # Read each .fasta file
  fasta_data <- read.fasta(file, seqtype = "AA")
  
  # Count the sequences containing 'SQ' patterns
  contains_SQ <- sapply(fasta_data, function(seq) grepl("sq(e|d)|SQ(E|D)", paste(seq, collapse = "")))
  count_SQ <- sum(contains_SQ) # Number of sequences containing 'SQ'
  
  # Append results to the data frame
  results <- rbind(results, data.frame(File = file, Count = count_SQ))
}

# View the table with counts for each file
print(results)

results$percentage <- round(((results$Count/184)*100), digits = 0)

results$File <- gsub("H2A_Pruned_H2A_Forced_2.fasta", "H2A Forced", results$File)
results$File <- gsub("H2A_Pruned_H2AX_Forced_2.fasta", "H2AX Forced", results$File)
results$File <- gsub("H2A_Truncated_2.fasta", "H2A Truncated", results$File)
results$File <- gsub("H2A_Unbiased_2.fasta", "H2A Unbiased", results$File)

# Set the custom order for the 'File' column
results$File <- factor(results$File, levels = c("H2A Unbiased", "H2A Truncated", "H2A Forced", "H2AX Forced"))

# Plot the data as a bar chart
ggplot(results, aes(x = File, y = percentage)) +
  geom_bar(stat = "identity", fill = "steelblue1") +
  labs(title = "Percentage of Species Where H2AX is Chosen After Pruning",
       x = "H2A Pruned",
       y = "Percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme_pubclean()
  scale_fill_viridis_d()  # Optional: use a color scale from the 'viridis' package for better visibility

ggsave(filename = "/Users/perrinm/Documents/AHocher_Project/PLOTS/Percentage_H2AX_Chosen.pdf")
ggsave()
