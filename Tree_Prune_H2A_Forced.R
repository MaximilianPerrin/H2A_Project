#Aim: keep a single protein per ortholog group of interest

rm(list = ls())

# Load the required library
# install.packages("ape")
# install.packages("adephylo")
# install.packages("seqinr")
library(ape)
library(adephylo)
library(seqinr)
library(tidyverse)
# Load your tree
HistoneTreePath <- "/Users/perrinm/Documents/AHocher_Project/TREE/H2A_TREE/Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_Curated_wHmfAB_mafftaligned.fa.clipkit.Rooted.treefile.txt"

  HistoneType=gsub("Individual_Euk_HF_hits_PostProcessed_","",basename(HistoneTreePath))
  HistoneType=strsplit(HistoneType,"_HMM_Subset")[[1]][1]
  
  tree <- read.tree(HistoneTreePath)
  
  Names=tree$tip.label
  tree$tip.label=sub("_","$",tree$tip.label)
  tree$tip.label = sub("_(?=[^_]*$)", "@", tree$tip.label, perl=TRUE)
  tree$tip.label=gsub("-","_",tree$tip.label)
  tree$tip.label <- gsub("\\$MAST_04C_sp_MAST_4C_sp1@", "\\$MAST-04C_sp_MAST-4C-sp1@", tree$tip.label)
  tree$tip.label <- gsub("\\$Colponemidia_sp_Colp_10@", "\\$Colponemidia_sp_Colp-10@", tree$tip.label)
  tree$tip.label <- gsub("\\$Telonema_sp_P_2@", "\\$Telonema_sp_P-2@", tree$tip.label)
  tree$tip.label <- gsub("\\$Hematodinium_sp_SG_2012@", "\\$Hematodinium_sp_SG-2012@", tree$tip.label)
  tree$tip.label <- gsub("\\$Blastocystis_sp_subtype4_isolateWR1@", "\\$Blastocystis_sp_subtype4-isolateWR1@", tree$tip.label)
  tree$tip.label <- gsub("\\$Apusomonadida_sp_AF_17@", "\\$Apusomonadida_sp_AF-17@", tree$tip.label)
  tree$tip.label <- gsub("\\$Colponemidia_sp_Colp_15@", "\\$Colponemidia_sp_Colp-15@", tree$tip.label)
  tree$tip.label <- gsub("\\$Choanocystis_sp_HF_7@", "\\$Choanocystis_sp_HF-7@", tree$tip.label)
  tree$tip.label <- gsub("\\$Acanthocystis_sp_HF_20@", "\\$Acanthocystis_sp_HF-20@", tree$tip.label)

  # Function to calculate branch length to the root for a given node
  Distances=adephylo::distRoot(tree,method="patristic")
  
  DistancesDt=data.frame(ID=names(Distances),Dist=as.numeric(Distances))
  DistancesDt$specie=gsub(".*\\$(.*?)@.*","\\1", DistancesDt$ID)
  DistancesDt$ToPrune=0
  
  #TO do remove h2aZ
  #Import H2AZ Subset
  H2AZ <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A/Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_OG3.fa", seqtype = "AA")
  og4 <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A/Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_OG4.fa")
  #Create list of all names 
  
  H2AZ_names <- getName(H2AZ)
  og4_name <- getName(og4)
  
  
  #Filter H2AZ from the list
  #DistancesDt$ToPrune[DistancesDt$ID %in% H2AZ_names] <- 1
  DistancesDt <- DistancesDt[which(!DistancesDt$ID %in% H2AZ_names), ]
  DistancesDt <- DistancesDt[which(!DistancesDt$ID %in% og4_name), ]
  
  #Read The fasta file to get sequences
  H2AFastaSeq <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A/Individual_Euk_HF_hits_PostProcessed_H2A_HMM_Subset_OG2.fa", seqtype = "AA")

  #Check for the "SQ" pattern in each sequence
  H2AX <- data.frame(sapply(H2AFastaSeq, function(seq) grepl("sq(e|d)|SQ(E|D)", paste(seq, collapse = ""))))
  H2AX_to_prune <- data.frame(getName(H2AFastaSeq), H2AX$sapply.H2AFastaSeq..function.seq..grepl..sq.e.d..SQ.E.D....paste.seq..)
  colnames(H2AX_to_prune)[1] = "Species_ID"
  colnames(H2AX_to_prune)[2] = "contains_sq"
  H2AX_to_prune$Species = gsub(".*\\$(.*?)@.*","\\1", H2AX_to_prune$Species_ID)
  H2AX_to_prune <- H2AX_to_prune[which(H2AX_to_prune$Species_ID %in% DistancesDt$ID), ]
  #Identify duplicate species and keep only species with more than one H2A
  Duplicalte_Species <- H2AX_to_prune %>%
    count(Species) %>%
    filter(n > 1)
  H2AX_to_prune <- H2AX_to_prune[which(H2AX_to_prune$Species %in% Duplicalte_Species$Species), ]
  
  #Group by species, keep only species that contain at least one sq, Remove species that only contain sq, keep the remaining sq
  H2AX_to_prune <- H2AX_to_prune %>%
    group_by(Species) %>%
    filter(any(contains_sq == TRUE)) %>%
    filter(any(contains_sq == FALSE)) %>%
    filter(contains_sq == TRUE)
  
  #set Distances_dt toprune = 1
  DistancesDt[which(DistancesDt$ID%in%H2AX_to_prune$Species_ID),]$ToPrune=1
  DistancesDt <- DistancesDt[which(!DistancesDt$ID %in% H2AX_to_prune$Species_ID), ]
  
  Lost_species <- DistancesDt %>%
    group_by(specie) %>%
    filter(all(ToPrune == 1))
    
    SpecieCount=table(DistancesDt$specie)
    
      #isolating the name of the species
      DupSpecies=unique(names(SpecieCount[which(SpecieCount>1)]))
      
      #Now this is the important bit, for each specie we select the histone with the shortest distance to the root
      #In pratice we flag all the seqeunce but the one to keep and we prune those from the tree later onn
      for (S in DupSpecies) {
        Subselecta=DistancesDt[which(DistancesDt$specie==S),]
        Subselecta=Subselecta[-which.min(Subselecta$Dist),]
        DistancesDt[which(DistancesDt$ID%in%Subselecta$ID),]$ToPrune=1}
      
  
  #Here we remove all the Flagged sequences
  PrunedTree=drop.tip(tree,H2AZ_names)
  PrunedTree=drop.tip(PrunedTree,H2AX_to_prune$Species_ID)
  PrunedTree=drop.tip(PrunedTree,DistancesDt[which(DistancesDt$ToPrune==1),]$ID)
  PrunedTree$tip.label=unlist(lapply(PrunedTree$tip.label,function(x) strsplit(x,"\\|")[[1]][1]))
  write.tree(PrunedTree,paste0("/Users/perrinm/Documents/AHocher_Project/TREE/",HistoneType,"_H2A_Forced_3.tree"))
  
  #Now output the sequences of that tree by ortholog group
  
  setwd("/Users/perrinm/Documents/AHocher_Project/TREE/")
  
  # Replace the text between $ and @ with "new_text"
  PrunedTree$tip.label <- gsub("\\$MAST_04C_sp_MAST_4C_sp1@", "\\$MAST-04C_sp_MAST-4C-sp1@", PrunedTree$tip.label)
  PrunedTree$tip.label <- gsub("\\$Colponemidia_sp_Colp_10@", "\\$Colponemidia_sp_Colp-10@", PrunedTree$tip.label)
  PrunedTree$tip.label <- gsub("\\$Telonema_sp_P_2@", "\\$Telonema_sp_P-2@", PrunedTree$tip.label)
  PrunedTree$tip.label <- gsub("\\$Hematodinium_sp_SG_2012@", "\\$Hematodinium_sp_SG-2012@", PrunedTree$tip.label)
  PrunedTree$tip.label <- gsub("\\$Blastocystis_sp_subtype4_isolateWR1@", "\\$Blastocystis_sp_subtype4-isolateWR1@", PrunedTree$tip.label)
  PrunedTree$tip.label <- gsub("\\$Apusomonadida_sp_AF_17@", "\\$Apusomonadida_sp_AF-17@", PrunedTree$tip.label)
  PrunedTree$tip.label <- gsub("\\$Colponemidia_sp_Colp_15@", "\\$Colponemidia_sp_Colp-15@", PrunedTree$tip.label)
  PrunedTree$tip.label <- gsub("\\$Choanocystis_sp_HF_7@", "\\$Choanocystis_sp_HF-7@", PrunedTree$tip.label)
  PrunedTree$tip.label <- gsub("\\$Acanthocystis_sp_HF_20@", "\\$Acanthocystis_sp_HF-20@", PrunedTree$tip.label)
  
  #Remove all sequences not found in pruned tree
  H2AFastaSeq <- H2AFastaSeq[names(H2AFastaSeq) %in% PrunedTree$tip.label]
  
  write.fasta(sequences = H2AFastaSeq, 
              names = names(H2AFastaSeq), 
              file.out = "H2A_H2A_Forced_3.fasta")
  
  