# The aim of this script is to produce a venn diagram of the overlap 

rm(list = ls())
library(seqinr)
library(tidyverse)
library(VennDiagram)


# Load pruned FASTA files
H2A_unbiased <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A/PRUNED/Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_OG2.fa", seqtype = "AA")
H2A_truncated <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/FINAL/H2A_Truncated_2.fasta", seqtype = "AA")

#Create a list of SQ containing species for each list

H2A_unbiasedSQ <- data.frame(species_ID = getName(H2A_unbiased))
H2A_unbiasedSQ$Species <- gsub(".*\\$(.*?)@.*", "\\1", H2A_unbiasedSQ$species_ID)
H2A_truncatedSQ <- data.frame(species_ID = getName(H2A_truncated))
H2A_truncatedSQ$Species <- gsub(".*\\$(.*?)@.*", "\\1", H2A_truncatedSQ$species_ID)

# Define the "SQ" motif pattern
sq_pattern <- "sq(e|d)|SQ(E|D)"

# Create a function to check if a sequence contains the "SQ" motif
contains_sq <- function(seq) {
  grepl(sq_pattern, paste(seq, collapse = ""))
}

# Identify sequences in H2A_Pruned_wSQ that contain the SQ motif
H2A_unbiasedSQ$contains_SQ <- sapply(getSequence(H2A_unbiased), contains_sq)
H2A_unbiasedSQ$contains_SQ <- as.integer(H2A_unbiasedSQ$contains_SQ)  # Convert to 1 (TRUE) / 0 (FALSE)

# Identify sequences in H2A_Pruned_woSQ that contain the SQ motif
H2A_truncatedSQ$contains_SQ <- sapply(getSequence(H2A_truncated), contains_sq)
H2A_truncatedSQ$contains_SQ <- as.integer(H2A_truncatedSQ$contains_SQ)  # Convert to 1 (TRUE) / 0 (FALSE)


# Join data frames by the Species column
merged_df <- full_join(H2A_unbiasedSQ %>% select(Species, contains_SQ) %>% rename(contains_SQ_unbiased = contains_SQ),
                       H2A_truncatedSQ %>% select(Species, contains_SQ) %>% rename(contains_SQ_truncated = contains_SQ),
                       by = "Species")

# Replace NA values with 0 for the Venn Diagram
merged_df$contains_SQ_unbiased[is.na(merged_df$contains_SQ_unbiased)] <- 0
merged_df$contains_SQ_truncated[is.na(merged_df$contains_SQ_truncated)] <- 0

# Create lists of species for each condition
species_unbiased <- merged_df$Species[merged_df$contains_SQ_unbiased == 1]
species_truncated <- merged_df$Species[merged_df$contains_SQ_truncated == 1]


# Draw Venn Diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(species_unbiased),
  area2 = length(species_truncated),
  cross.area = length(intersect(species_unbiased, species_truncated)),
  category = c("Unbiased", "Truncated"),
  fill = c("grey80", "steelblue1"),
  lty = "blank"
)

#Save the venn diagram as a .pdf

pdf(file = "/Users/perrinm/Documents/AHocher_Project/PLOTS/h2ax_venn_diagram.pdf")

# Draw the Venn diagram and close the device
grid.draw(venn.plot)
dev.off()


