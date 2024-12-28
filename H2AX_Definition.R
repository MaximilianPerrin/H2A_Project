#The aim of this script is to define what we select as H2AX
#The tail of the histone will be read and the SQ closest to the end will be chosen
#Subsequerntly it will be determined which AA follows SQ
# Clear the environment
rm(list = ls())

# Load required libraries
library(seqinr)
library(tidyverse)

# Load in fasta file
H2A_tail <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_OG2_Tail/H2A_OG2_Tail_HSapiens_121.fasta", seqtype = "AA")

# Look for SQ motif and get the following AAs
All_SQ <- lapply(getSequence(H2A_tail), function(seq) {
  seq_collapsed <- paste(seq, collapse = "")
  
  matches <- gregexpr("sq|SQ", seq_collapsed, ignore.case = TRUE)[[1]]
  
  if (matches[1] == -1) {
    # No match found
    return(list(Closest_SQ_Position = NA, Next_AA = NA, AA_After_Next = NA))
  } else {
    # Get the closest match to the end of the sequence
    closest_to_end <- max(matches)
    
    # Get the AA following the SQ motif
    next_pos <- closest_to_end + 2  # Add 2 because "SQ" has length 2
    next_aa <- ifelse(next_pos <= nchar(seq_collapsed), substr(seq_collapsed, next_pos, next_pos), NA)
    
    # Get the AA following the next AA
    after_next_pos <- next_pos + 1
    after_next_aa <- ifelse(after_next_pos <= nchar(seq_collapsed), substr(seq_collapsed, after_next_pos, after_next_pos), NA)
    
    return(list(Closest_SQ_Position = closest_to_end, Next_AA = next_aa, AA_After_Next = after_next_aa))
  }
})

# Extract results into a data frame
All_SQ_df <- do.call(rbind, lapply(All_SQ, function(x) data.frame(t(unlist(x)), stringsAsFactors = FALSE)))

# Add sequence IDs
All_SQ_df$ID <- getName(H2A_tail)

# Reorder columns
All_SQ_df <- All_SQ_df[, c("ID", "Closest_SQ_Position", "Next_AA", "AA_After_Next")]

# View the final data frame
print(All_SQ_df)

All_SQ_df$Closest_SQ_Position[!is.na(All_SQ_df$Closest_SQ_Position)] = "SQ"

#Filter for SQ

filtered_data <- All_SQ_df %>% 
  filter(!is.na(Closest_SQ_Position))

#Filter E or D
E_or_D <- filtered_data %>%
  filter(Next_AA == "E" | Next_AA == "D")

#Filter for A
A <- filtered_data %>%
  filter(Next_AA == "A")
#Filter for hydrophobic followed by E

Hydrophobe_E <- E %>%
  filter(AA_After_Next == "A" | AA_After_Next == "I" | AA_After_Next == "L"| AA_After_Next == "M"| AA_After_Next == "P"| AA_After_Next == "V" | 
           AA_After_Next == "F" | AA_After_Next == "W"| AA_After_Next == "C"| AA_After_Next == "Y")

#Filter for D followed by hydrophobic residue
Hydrophobe_D <- D %>%
  filter(AA_After_Next == "A" | AA_After_Next == "I" | AA_After_Next == "L"| AA_After_Next == "M"| AA_After_Next == "P"| AA_After_Next == "V" | 
           AA_After_Next == "F" | AA_After_Next == "W"| AA_After_Next == "C"| AA_After_Next == "Y")

#Filter fore A followed by hydrophobic residue
Hydrophobe_A <- A %>%
  filter(AA_After_Next == "A" | AA_After_Next == "I" | AA_After_Next == "L"| AA_After_Next == "M"| AA_After_Next == "P"| AA_After_Next == "V" | 
           AA_After_Next == "F" | AA_After_Next == "W"| AA_After_Next == "C"| AA_After_Next == "Y")

E <- filtered_data %>%
  filter(Next_AA == "E")
D <- filtered_data %>%
  filter(Next_AA == "D")

