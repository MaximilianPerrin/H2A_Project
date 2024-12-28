#This script is to produce a table with protein domains and correlated them with precence or absence of H2AX

rm(list = ls())
library(seqinr)

#Set the working directory
setwd("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/")

#Import sequences
AllH2A=read.fasta("Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_OG2.fa",seqtype = "AA")

getSequence(AllH2A)
getName(AllH2A)


# Step 2: Check for the "SQ" pattern in each sequence
contains_SQ <- sapply(AllH2A, function(seq) grepl("sq(e|d)|SQ(E|D)", paste(seq, collapse = "")))

# Step 2: Create a data frame with the names of the proteins
H2A_table <- data.frame(OriginalName = names(AllH2A))

# Step 3: Create another column with only the name of the species
H2A_table$Specie <- gsub(".*\\$(.*?)@.*", "\\1", H2A_table$OriginalName)

# Step 4: Count H2As with and without SQ
H2A_table$Contains_SQ <- contains_SQ


# Step 5: Summarize counts
library(dplyr)

# Step 1: Summarize species by H2A copy number and SQ presence
H2ASummary <- H2A_table %>%
  group_by(Specie) %>%
  summarize(
    Total_H2A = n(),                           # Total H2A copies for each species
    Contains_SQ = as.integer(sum(Contains_SQ) > 0)  # Binary indicator for "SQ" motif presence
  )

setwd("/Users/perrinm/Documents/AHocher_Project/PROTEIN_DOMAINS/")

#Read the text file containing the domains
pfam_dom <- read.table("Pfam_vs_EukprotCLS.txt", header = TRUE, stringsAsFactors = FALSE)

# Merge the Contains_SQ column from H2ASummary into pfam_dom by Specie_Name
pfam_dom <- merge(pfam_dom, H2ASummary, by = "Specie", all.x = TRUE)

# Convert all numerical data into binary format (1 if >= 1, else 0)
pfam_dom[sapply(pfam_dom, is.numeric)] <- lapply(pfam_dom[sapply(pfam_dom, is.numeric)], function(x) ifelse(x >= 1, 1, 0))

library(tidyverse)

#Remove NA 
pfam_cor <- drop_na(pfam_dom, Total_H2A, Contains_SQ)

# Initialize a list to store correlation and p-value results
cor_results <- list()
p_values <- list()  # New list to store p-values

# Identify numeric columns in the data frame
numeric_columns <- sapply(pfam_cor, is.numeric)

# Loop through each numeric column
for (col_name in names(pfam_cor)[numeric_columns]) {
  if (col_name != "Contains_SQ") {  # Skip Contains_SQ itself
    # Calculate standard deviation to check for zero variance
    sd_Contains_SQ <- sd(pfam_cor$Contains_SQ, na.rm = TRUE)
    sd_other <- sd(pfam_cor[[col_name]], na.rm = TRUE)
    
    # Check if both columns have variation
    if (sd_Contains_SQ > 0 && sd_other > 0) {
      # Perform correlation test to obtain both correlation and p-value
      test_result <- cor.test(pfam_cor$Contains_SQ, pfam_cor[[col_name]], method = "pearson", use = "complete.obs")
      
      # Store correlation and p-value in separate lists
      cor_results[[col_name]] <- test_result$estimate  # Correlation coefficient
      p_values[[col_name]] <- test_result$p.value      # P-value
    } else {
      # Handle case where one of the columns has no variation
      cor_results[[col_name]] <- NA  # Or any other indicator you prefer
      p_values[[col_name]] <- NA     # Same for p-value
    }
  }
}

# Convert the results to a data frame for better readability
cor_results_df <- data.frame(
  Variable = names(cor_results), 
  Correlation = unlist(cor_results),
  P_Value = unlist(p_values)  # Add p-values as a new column
)

# Display the results
print(cor_results_df)

#Adjust PValue to account for multiple tests

cor_results_df$Adjusted_PValue = p.adjust(cor_results_df$P_Value,method = "fdr")

#filter out significant correlations
Sig_Cor_Pfam <- filter(cor_results_df, Adjusted_PValue <= 0.05)

sig_pfam <- select(pfam_dom, Specie, Sig_Cor_Pfam$Variable)
