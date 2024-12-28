#The goal of this script is to compare the FASTA files from two different 
#phylogentic trees to determine if the tree pruning is the same indpendant of SQ motif

rm(list = ls())
library(seqinr)
library(tidyverse)

#Load in pruned FASTA Files

H2A_Pruned_wSQ <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A/PRUNED/Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_OG2.fa", seqtype = "AA")
H2A_Pruned_woSQ <- read.fasta("/Users/perrinm/Documents/AHocher_Project/TREE/H2A_Minus_SQ_Minus_H2AZ_2.fasta", seqtype = "AA")

#Identify sequences that appear in tree 1 (unbiased) one but not on tree 2 (with sq removed)

tree_diff <- getName(H2A_Pruned_wSQ)[!getName(H2A_Pruned_wSQ) %in% getName(H2A_Pruned_woSQ)]
print(tree_diff)  # Elements in T1 but not in T2

#Build a table with list of species that vary. 3 possible scenaros: Not_SQ -> SQ; SQ -> Not_SQ; Not_SQ -> Not_SQ

# Create a data frame with a character column
tree_diff <- data.frame(species_ID = tree_diff)
H2A_Pruned_wSQ_df <- data.frame(species_ID = getName(H2A_Pruned_wSQ))
H2A_Pruned_woSQ_df <- data.frame(species_ID = getName(H2A_Pruned_woSQ))


# create species_ID_Clean
tree_diff$species_ID_Clean = gsub(".*\\$(.*?)@.*", "\\1", tree_diff$species_ID)
H2A_Pruned_wSQ_df$species_ID_Clean <- gsub(".*\\$(.*?)@.*", "\\1", H2A_Pruned_wSQ_df$species_ID)
H2A_Pruned_woSQ_df$species_ID_Clean <- gsub(".*\\$(.*?)@.*", "\\1", H2A_Pruned_woSQ_df$species_ID)

