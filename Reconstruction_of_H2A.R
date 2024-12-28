#This script is going to be used to reconstruct the ancestral sequence of H2A
rm(list = ls())
library(tidyverse)
library(seqinr)
library(ape)
library(adephylo)

#Import Fasta and add archea sequence
fasta_file <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/RECONSTRUCTION_3/H2A_Unbiased_3.fasta", seqtype = "AA")
archea <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/FINAL/HmfA.fa", seqtype = "AA")

#Read the tree to get names to filter
tree <- read.tree("/Users/perrinm/Documents/AHocher_Project/TREE/RECONSTRUCTION/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")

#Change the names in the fasta file to be the species only
names(fasta_file) <- gsub(".*\\$(.*?)@.*","\\1",getName(fasta_file))
tree$tip.label

#Filter the fasta for complete histone set
fasta_file <- fasta_file[names(fasta_file) %in% tree$tip.label]

#Align the Files using MAFFT

system(paste0("mafft-linsi --op 0.75 ",fasta_file,"> ",gsub(".fa","_linsi_op075.fa",i)))
