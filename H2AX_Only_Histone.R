#The aim of this code is to produce a list of species for which H2AX is their only H2A

#Clear workspace and load librarys
rm(list = ls())
library(seqinr)
library(tidyverse)

#Set the working directory
setwd("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/")

#Import sequences
AllH2A=read.fasta("Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_OG2.fa",seqtype = "AA")

getSequence(AllH2A)
getName(AllH2A)

#Create a table with the names of the proteins
H2A_table=data.frame(OriginalName=getName(AllH2A))

head(H2A_table)

#Create another column with only the name of the specie
H2A_table$Specie_Name=gsub(".*\\$(.*?)@.*","\\1",H2A_table$OriginalName)

#Create a summary table (using the table() function) of the specie, as a way to count H2As
H2APerSpecie=as.data.frame(table(H2A_table$Specie_Name))

#rename the columns
names(H2APerSpecie)=c("Specie","Nb_H2A")

head(H2APerSpecie)

# Step 2: Check for the "SQ" pattern in each sequence
contains_SQ <- sapply(AllH2A, function(seq) grepl("sq(e|d)|SQ(E|D)", paste(seq, collapse = "")))

# Step 2: Create a data frame with the names of the proteins
H2A_table <- data.frame(OriginalName = names(AllH2A))

# Step 3: Create another column with only the name of the species
H2A_table$Specie_Name <- gsub(".*\\$(.*?)@.*", "\\1", H2A_table$OriginalName)

# Step 4: Count H2As with and without SQ
H2A_table$Contains_SQ <- contains_SQ

#Add a column for if species contains SQ

SQ_Species <- H2A_table %>%
  group_by(Specie_Name) %>%
  filter(any(Contains_SQ == TRUE)) %>%
  summarise()

SQ_Multiple <- H2A_table %>%
  group_by(Specie_Name) %>%
  filter(!any(Contains_SQ == FALSE)) %>%
  summarise()

#New column for containing sq
H2APerSpecie$H2AX_Only = 0

#Look into H2APerSpecie and if Species is found in SQ list then change conatins SQ to 1 
#Edit code to include or found in list and all values are true 
H2APerSpecie[which(H2APerSpecie$Specie %in% SQ_Multiple$Specie_Name), ]$H2AX_Only = 1

#Filter table to just binary H2AX containing
H2AX_Only <- H2APerSpecie %>%
  select(Specie, H2AX_Only)


