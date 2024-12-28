#This script is going to be used to reconstruct the ancestral sequence of H2A
rm(list = ls())
library(tidyverse)
library(seqinr)
library(ape)
library(adephylo)

#Import Fasta and add archea sequence
fasta_file <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/R5/A_Hocher_Seqs.fasta", seqtype = "AA")
archea <- read.fasta("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/FINAL/HmfA.fa", seqtype = "AA")

#Read the tree to get names to filter
tree <- read.tree("/Users/perrinm/Documents/AHocher_Project/TREE/RECONSTRUCTION/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")

#Change the names in the fasta file to be the species only
names(fasta_file) <- gsub(".*\\$(.*?)@.*","\\1",getName(fasta_file))
tree$tip.label

#Filter the fasta for complete histone set
fasta_file <- fasta_file[names(fasta_file) %in% tree$tip.label]
#Write Fasta file 
setwd("/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/R5/Aligned")
write.fasta(sequences = fasta_file, 
            names = names(fasta_file), 
            file.out = "A_Hocher_Seqs_H2A_Unbiased_4_Species_Only.fasta")


#Align the Files using MAFFT


# Define the directory where the histone sequence files are located.
FastaDir="/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/R5/Aligned"

# List all files in the FastaDir that match the pattern "_Complete_subset.fa" to identify histone sequence files.
HistonesFiles=list.files(FastaDir)

# Set the working directory to FastaDir where the fasta files are stored.
setwd(FastaDir)

# Loop through each histone file for alignment and renaming based on species.
for(i in HistonesFiles){
  # Align sequences within each file using MAFFT with an open penalty of 0.75, and output the aligned sequences.
  system(paste0("mafft-linsi --op 0.75 ",i,"> ",gsub(".fa","_linsi_op075.fa",i)))
  
  # Replace the file extension of the input file to match the output file from the alignment.
  FastaPath=gsub(".fa","_linsi_op075.fa",i)
  
  # Read the aligned fasta file.
  FastaIn=read.fasta(FastaPath)
  
  write.fasta(getSequence(FastaIn),getName(FastaIn),FastaPath)
}

HmfFilepath="/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/FINAL/HmfA.fa"
AlignmentPath="/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/R4/Aligned/H2A_Unbiased_4_Species_Only_linsi_op075.fasta"
AlignmentPathNew="/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/R4/Aligned/H2A_Unbiased_Species_Only_linsi_op075_Warchea.fasta"
Command2=paste("mafft --addfragments ",HmfFilepath," ",AlignmentPath," > ",AlignmentPathNew,sep="")
system(Command2)


#Use tree to reconstruct the ancestral sequence

setwd("/Users/perrinm/Documents/AHocher_Project/TREE/RECONSTRUCTION/H2A_Unbiased")
#Specify File
AlignedFile = "/Users/perrinm/Documents/AHocher_Project/SEQUENCES/H2A_MP/H2A_Sequence_Reconstruction/RECONSTRUCTION/H2A_H2AX_Forced_linsi_op075.fasta"

# For each aligned file, perform ancestral sequence reconstruction using IQ-TREE.
  # Specify the rooted tree for use with IQ-TREE.
  RootedTree="/Users/perrinm/Documents/AHocher_Project/TREE/RECONSTRUCTION/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree"
  # Construct the command for IQ-TREE, specifying input files and options for ancestral state reconstruction.
  CommandLine=paste0("/Applications/iqtree-2.2.5-MacOSX/bin/iqtree2 -nt AUTO --redo -asr --keep-ident -s ",AlignedFile," -te ",RootedTree)
  # Execute the command.
  system(CommandLine)

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Analysis section: processes the results of ASR to analyze and write out the reconstructed sequences and their states.
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  setwd("/Users/perrinm/Documents/AHocher_Project/TREE/RECONSTRUCTION/H2A_H2AX_Forced")
  
    State = read.table("/Users/perrinm/Documents/AHocher_Project/TREE/RECONSTRUCTION/H2A_H2AX_Forced/H2A_H2AX_Forced_linsi_op075.fasta.state", skip=1,header=T)
    # Extract and process unique node names and their associated sequences.
    NodeNames=unique(State$Node)
    Sequences=data.frame(Node=NodeNames,Sequence=unlist(lapply(NodeNames,function(x) paste(State[which(State$Node==x),]$State,collapse = ""))))
    write.fasta(as.list(Sequences$Sequence),Sequences$Node,file.out = "H2A_H2AX_Forced_linsi_op075_States.fasta")
  
  
  

  
  # Further analysis includes linking ancestral state probabilities with gap percentages in the original sequences
  #To get all tips from a given node
  SpeciesTree=ape::read.tree(RootedTree)
  NodesCorresp=adephylo::listTips(SpeciesTree)
  NodesCorresp[1]
  
  #Building a table with nodenames and corresponding species
  DTNodes=data.frame(NodeID=rep("Eukaryota",length(SpeciesTree$tip.label)),SpeciesName=SpeciesTree$tip.label)
  for(k in 2:length(names(NodesCorresp))){
    TheNode=names(NodesCorresp)[k]
    Species=names(NodesCorresp[TheNode][[1]])
    TmpDT=data.frame(NodeID=rep(TheNode,length(Species)),SpeciesName=Species)
    DTNodes=rbind(DTNodes,TmpDT)}
  
  
  
    # Read the original sequences to calculate gap percentages.
    Sequences=seqinr::read.fasta(AlignedFile,seqtype = "AA")
    
    # Calculate gap percentages for each sequence.
    seqMatrix <- sapply(Sequences, function(x) unlist(strsplit(x, "")))
    gapPercentages <- rowMeans(seqMatrix == "-") * 100
    
    # Read the state reconstruction results again.
    States = read.table("/Users/perrinm/Documents/AHocher_Project/TREE/RECONSTRUCTION/H2A_H2AX_Forced/H2A_H2AX_Forced_linsi_op075.fasta.state", skip=1,header=T)
    
    # Update the data with gap percentages.
    NodeNames=unique(States$Node)
    States$GapPercent=rep(gapPercentages,length(NodeNames))
    
    States$GapPercentPerNode=NA
    #Now for each node compute the %of gap per position (not just global)
    for (UNode in NodeNames){
      NodeSequenceNames=DTNodes[which(DTNodes$NodeID==UNode),]$SpeciesName
      if(length(NodeSequenceNames)>0){
        seqMatrix <- sapply(Sequences[which(getName(Sequences)%in%NodeSequenceNames)], function(x) unlist(strsplit(x, "")))
        gapPercentages <- rowMeans(seqMatrix == "-") * 100
        States[which(States$Node==UNode),]$GapPercentPerNode=gapPercentages}
    }
    
    
    
    
    SequencesDT=data.frame(Node=NodeNames,Sequence=unlist(lapply(NodeNames,function(x) paste(States[which(States$Node==x),]$State,collapse = ""))))
    
    
    
    
    # Identify the maximum probability for each state.
    States$MaxP=apply((States[,4:23]),1,max)
    
    #Count the number of states to reach >0.95 certainty
    # Apply the function to each row of the dataframe
    States$NumStates9 = apply(States[,4:23], 1, function(row_data) {
      # Calculate the cumulative sum
      cum_sum <- cumsum(sort(row_data, decreasing = TRUE))
      # Find where the cumulative sum exceeds 0.9
      return(which(cum_sum > 0.9)[1])
    })
    States$NumStates95 = apply(States[,4:23], 1, function(row_data) {
      # Calculate the cumulative sum
      cum_sum <- cumsum(sort(row_data, decreasing = TRUE))
      # Find where the cumulative sum exceeds 0.9
      return(which(cum_sum > 0.95)[1])
    })
   
    #Now filtering on the percent of Gaps per node
    States50perNode=States[which(States$GapPercentPerNode<50),]
    States50perNode$TotalNumberOfSeq09=NA
    for(u in unique(States50perNode$Node)){
      States50perNode[which(States50perNode$Node==u),]$TotalNumberOfSeq09=prod(States50perNode[which(States50perNode$Node==u),]$NumStates9)}
    write.table(States50perNode,file = "H2A_H2AX_Forced_linsi_op075.fasta.statestatewithoutGapsPerNode.txt",sep="\t",quote=F,row.names = F)
    
    
    #Export sequences with length adjusted per node 
    SequencesDTPerNode=data.frame(Node=NodeNames,Sequence=unlist(lapply(NodeNames,function(x) paste(States50perNode[which(States50perNode$Node==x),]$State,collapse = ""))))
    
    #Order by node 
    SeqReordered=SequencesDTPerNode[match(unique(DTNodes$NodeID),SequencesDTPerNode$Node),]
    
    write.fasta(as.list(SeqReordered$Sequence),SeqReordered$Node, file.out =  "H2A_H2AX_Forced_linsi_op075_nsequences_GapsPerNodes.fa")
    
  
  
  # Filtering steps, plotting, and saving analysis results involve selective sequence analysis, plotting alignment, and gap percentage visualization, concluding with saving the plots.
  
  

    
   
    Sequences=read.fasta(AlignedFile,seqtype = "AA")
    
    seqMatrix <- sapply(Sequences, function(x) unlist(strsplit(x, "")))
    gapPercentages <- rowMeans(seqMatrix == "-") * 100
    
    #Analyse results
    States = read.table("/Users/perrinm/Documents/AHocher_Project/TREE/RECONSTRUCTION/H2A_H2AX_Forced/H2A_H2AX_Forced_linsi_op075.fasta.state", skip=1,header=T)
    States$MaxP=apply((States[,4:23]),1,max)
    
    NodeNames=unique(States$Node)
    States$GapPercent=rep(gapPercentages,length(NodeNames))
    
    StatesNoGap=States[which(States$GapPercent<50),]
    
    #Added to make sure every node is in the right order when plotting with gmsa
    StatesNoGap=StatesNoGap[order(StatesNoGap$Node),]
    
    SequencesDTNoGap=data.frame(Node=NodeNames,Sequence=unlist(lapply(NodeNames,function(x) paste(StatesNoGap[which(StatesNoGap$Node==x),]$State,collapse = ""))))
    
    NodeToPlot=rev(c("Node_1","Node_6","Eukaryota"))
    
    write.fasta(as.list(SequencesDTNoGap$Sequence[which(NodeNames%in%NodeToPlot)]),SequencesDTNoGap$Node[which(NodeNames%in%NodeToPlot)],file.out = "H2A_H2AX_Forced_linsi_op075_sequencesGap50.fa")
    
    #Export all nodes
    write.fasta(as.list(SequencesDTNoGap$Sequence),SequencesDTNoGap$Node,file.out = "H2A_H2AX_Forced_linsi_op075_AllNodes.sequencesGap50.fa")
    

    library(ggmsa);library(ggplot2);library(ggpubr)
    
    # Visualization of the sequences and their attributes using ggmsa and ggplot2 for a comprehensive analysis.
    # ggmsa is used for visualizing the multiple sequence alignment of the selected sequences.
    SeqGap50 = "/Users/perrinm/Documents/AHocher_Project/TREE/RECONSTRUCTION/H2A_H2AX_Forced/H2A_H2AX_Forced_linsi_op075_sequencesGap50.fa"
    PLOTA=ggmsa(SeqGap50, start =1, end = 200, char_width = 0.5, seq_name = T,color = "Shapely_AA")
    
    # ggplot is used to create a bar chart of gap percentages for a specific node, providing insights into sequence conservation.
    SeqLength=dim(StatesNoGap[which(StatesNoGap$Node=="Eukaryota"),])[1]
    BarGap=ggplot(data = StatesNoGap[which(StatesNoGap$Node=="Eukaryota"),],aes(x=as.factor(1:SeqLength),y=GapPercent))+geom_bar(stat="identity")+theme_pubclean()+theme(aspect.ratio=1/16)+xlab("")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
    # Subset data for selected nodes and prepare for probability visualization.
    SubSelecta=StatesNoGap[which(StatesNoGap$Node%in%NodeToPlot),c("Site","Node","MaxP")]
    SubSelecta$Site=rep(1:length(which(StatesNoGap$Node=="Eukaryota")),length(NodeToPlot))
    
    SubSelecta$Node=factor(SubSelecta$Node,levels = NodeToPlot)
    
    #Discretize the probabilities
    SubSelecta$DiscreteMaxP=NA
    if(length(which(SubSelecta$MaxP<0.8))>0){SubSelecta[which(SubSelecta$MaxP<0.8),]$DiscreteMaxP="inf08"}
    if(length(which(SubSelecta$MaxP>=0.8 & SubSelecta$MaxP<0.9))>0){SubSelecta[which(SubSelecta$MaxP>=0.8 & SubSelecta$MaxP<0.9),]$DiscreteMaxP="08_09"}
    if(length(which(SubSelecta$MaxP>=0.9 & SubSelecta$MaxP<0.95))>0){SubSelecta[which(SubSelecta$MaxP>=0.9 & SubSelecta$MaxP<0.95),]$DiscreteMaxP="09_095"}
    if(length(which(SubSelecta$MaxP>=0.95 & SubSelecta$MaxP<1))>0){SubSelecta[which(SubSelecta$MaxP>=0.95 & SubSelecta$MaxP<1),]$DiscreteMaxP="095_less1"}
    
    SubSelecta[which(SubSelecta$MaxP==1),]$DiscreteMaxP="1"
    
    Mycolor=c("orange","lightblue","lightgreen","#2dba4e","#FF7F7F")
    names(Mycolor)=c("08_09","09_095","095_less1","1","inf08")
    # Use ggplot to visualize the maximum probability of ancestral states across selected sites and nodes.
    Probs=ggplot(data=SubSelecta,aes(x=Site,y=Node,fill=DiscreteMaxP))+geom_tile(color = "white") +theme_msa()+scale_fill_manual(values=Mycolor)+theme(legend.position = "top")
    
    # Use ggpubr to arrange and save the plots into a single file for publication or further analysis.
    library(ggpubr)
    
    
    # Use Sys.Date() to get the current date in YYYY-MM-DD format and convert it to a string.
    currentDate <- as.character(Sys.Date())
    
    # Construct the filename by including the current date in the desired location within the string.
    filenameWithDate <- "/Users/perrinm/Documents/AHocher_Project/TREE/RECONSTRUCTION/H2A_H2AX_Forced/H2A_H2AX_Forced_linsi_op075_sequencesGap50.pdf"
    
    # Use ggsave with the filename that now includes the date.
    ggsave(plot = ggarrange(PLOTA, Probs, BarGap, ncol=1, align = "v", nrow=3), filename = filenameWithDate,height=8,width=8)

  
  

  
  
  
  
  
  #Manually export yeast and human sequences to compare with the consensus
  H2A=read.fasta("H2A_Complete_subset_linsi_op075_SpecieNameOnly.fa")
  
  H2A.subset=H2A[c(grep("sapiens",getName(H2A)),grep("cerevisiae",getName(H2A)))]
  write.fasta(getSequence(H2A.subset),getName(H2A.subset),"H2A_Human_yeast.fa")
  system("mafft-linsi H2A_Human_yeast.fa > H2A_Human_yeast_ali.fa")
  YHPLot=ggmsa("H2A_Human_yeast_ali.fa", start =1, end = 200, char_width = 0.5, seq_name = T,color = "Shapely_AA")
  ggsave(plot = YHPLot,"/Users/ahocher/Dropbox/Laboratory/ASGARD/PLOTS/EUK_HISTONES_PARALOGS/ANCESTRAL/Human_yeast_H2A.pdf")
  
  H2B=read.fasta("H2B_Complete_subset_linsi_op075_SpecieNameOnly.fa")
  H2B.subset=H2B[c(grep("sapiens",getName(H2B)),grep("cerevisiae",getName(H2B)))]
  write.fasta(getSequence(H2B.subset),getName(H2B.subset),"H2B_Human_yeast.fa")
  system("mafft-linsi H2B_Human_yeast.fa > H2B_Human_yeast_ali.fa")
  YHPLot=ggmsa("H2B_Human_yeast_ali.fa", start =1, end = 200, char_width = 0.5, seq_name = T,color = "Shapely_AA")
  ggsave(plot = YHPLot,"/Users/ahocher/Dropbox/Laboratory/ASGARD/PLOTS/EUK_HISTONES_PARALOGS/ANCESTRAL/Human_yeast_H2B.pdf")
  
  H3=read.fasta("H3_r_Complete_subset_linsi_op075_SpecieNameOnly.fa")
  H3.subset=H3[c(grep("sapiens",getName(H3)),grep("cerevisiae",getName(H3)))]
  write.fasta(getSequence(H3.subset),getName(H3.subset),"H3_Human_yeast.fa")
  system("mafft-linsi H3_Human_yeast.fa > H3_Human_yeast_ali.fa")
  YHPLot=ggmsa("H3_Human_yeast_ali.fa", start =1, end = 200, char_width = 0.5, seq_name = T,color = "Shapely_AA")
  ggsave(plot = YHPLot,"/Users/ahocher/Dropbox/Laboratory/ASGARD/PLOTS/EUK_HISTONES_PARALOGS/ANCESTRAL/Human_yeast_H3.pdf")
  
  H4=read.fasta("H4_Complete_subset_linsi_op075_SpecieNameOnly.fa")
  H4.subset=H4[c(grep("sapiens",getName(H4)),grep("cerevisiae",getName(H4)))]
  write.fasta(getSequence(H4.subset),getName(H4.subset),"H4_Human_yeast.fa")
  system("mafft-linsi H4_Human_yeast.fa > H4_Human_yeast_ali.fa")
  YHPLot=ggmsa("H4_Human_yeast_ali.fa", start =1, end = 200, char_width = 0.5, seq_name = T,color = "Shapely_AA")
  ggsave(plot = YHPLot,"/Users/ahocher/Dropbox/Laboratory/ASGARD/PLOTS/EUK_HISTONES_PARALOGS/ANCESTRAL/Human_yeast_H4.pdf")
  
  
  PLOTA=ggmsa(gsub(".fa",".sequencesGap50.fa",A), start =1, end = 200, char_width = 0.5, seq_name = T,color = "Shapely_AA")
  
  
  
  
  
  
  #To visualise aligned sequences side by side with the tree, we clip them
  clipkit /Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/H2A_Complete_subset_linsi_op075_SpecieNameOnly.fa --mode gappy
  clipkit /Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/H2A_Z_Complete_subset_linsi_op075_SpecieNameOnly.fa --mode gappy 
  clipkit /Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/H2B_Complete_subset_linsi_op075_SpecieNameOnly.fa --mode gappy
  clipkit /Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/H3_r_Complete_subset_linsi_op075_SpecieNameOnly.fa --mode gappy
  clipkit /Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/H4_Complete_subset_linsi_op075_SpecieNameOnly.fa --mode gappy
  
  
  
  
  
  ###############################################
  #To visualise the Nb of ancestral possiblities
  # On the ancestral tree
  ##############################################
  
  SpecieTree=ape::read.tree("/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")
  
  #H2A
  NodesCorresp=names(adephylo::listTips(SpecieTree))
  
  StatesFile="/Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/AROGER_TREE/A_H2A_Complete_subset_linsi_op075_SpecieNameOnly.fa.statewithoutGapsPerNode"
  StatesH2A=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  StatesH2A.nd=StatesH2A[which(duplicated(StatesH2A$Node)==F),]
  NbOfSeauences=log10(StatesH2A.nd[match(NodesCorresp,StatesH2A.nd$Node),]$TotalNumberOfSeq09)
  NbOfSeauences[which(NbOfSeauences>6)]=6
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H2A_NumberOfSeqslog10.tree")
  
  
  # Focusing on AA changes located on H2A core (defined based on S.cerevisiae histone sequence and conservation)
  #Residues positions 275 to 371
  StatesH2A=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  
  StatesH2ACore=StatesH2A[which(StatesH2A$Site> 275 & StatesH2A$Site< 371),]
  
  for(u in unique(StatesH2ACore$Node)){
    StatesH2ACore[which(StatesH2ACore$Node==u),]$TotalNumberOfSeq09=prod(StatesH2ACore[which(StatesH2ACore$Node==u),]$NumStates9)}
  
  StatesH2A.nd=StatesH2ACore[which(duplicated(StatesH2ACore$Node)==F),]
  NbOfSeauences=log10(StatesH2A.nd[match(NodesCorresp,StatesH2A.nd$Node),]$TotalNumberOfSeq09)
  if(length(which(NbOfSeauences>6))>0){
    NbOfSeauences[which(NbOfSeauences>6)]=6}
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H2A_CORE_NumberOfSeqslog10.tree")
  
  
  
  
  #H2B
  
  
  SpecieTree=ape::read.tree("/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")
  
  NodesCorresp=names(adephylo::listTips(SpecieTree))
  
  StatesFile="/Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/AROGER_TREE/A_H2B_Complete_subset_linsi_op075_SpecieNameOnly.fa.statewithoutGapsPerNode"
  StatesH2B=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  StatesH2B.nd=StatesH2B[which(duplicated(StatesH2B$Node)==F),]
  NbOfSeauences=log10(StatesH2B.nd[match(NodesCorresp,StatesH2B.nd$Node),]$TotalNumberOfSeq09)
  NbOfSeauences[which(NbOfSeauences>6)]=6
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H2B_NumberOfSeqslog10.tree")
  
  
  
  # Focusing on AA changes located on H2B core (defined based on S.cerevisiae histone sequence and conservation)
  #Residues positions 275 to 371
  StatesH2B=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  
  StatesH2BCore=StatesH2B[which(StatesH2B$Site> 265 & StatesH2B$Site< 355),]
  
  for(u in unique(StatesH2BCore$Node)){
    StatesH2BCore[which(StatesH2BCore$Node==u),]$TotalNumberOfSeq09=prod(StatesH2BCore[which(StatesH2BCore$Node==u),]$NumStates9)}
  
  StatesH2B.nd=StatesH2BCore[which(duplicated(StatesH2BCore$Node)==F),]
  NbOfSeauences=log10(StatesH2B.nd[match(NodesCorresp,StatesH2B.nd$Node),]$TotalNumberOfSeq09)
  if(length(which(NbOfSeauences>6))>0){
    NbOfSeauences[which(NbOfSeauences>6)]=6}
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H2B_CORE_NumberOfSeqslog10.tree")
  
  
  #H3_r
  
  
  SpecieTree=ape::read.tree("/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")
  
  NodesCorresp=names(adephylo::listTips(SpecieTree))
  
  StatesFile="/Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/AROGER_TREE/A_H3_r_Complete_subset_linsi_op075_SpecieNameOnly.fa.statewithoutGapsPerNode"
  StatesH3_r=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  StatesH3_r.nd=StatesH3_r[which(duplicated(StatesH3_r$Node)==F),]
  NbOfSeauences=log10(StatesH3_r.nd[match(NodesCorresp,StatesH3_r.nd$Node),]$TotalNumberOfSeq09)
  NbOfSeauences[which(NbOfSeauences>6)]=6
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H3_r_NumberOfSeqslog10.tree")
  
  #H4
  
  
  
  SpecieTree=ape::read.tree("/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")
  
  NodesCorresp=names(adephylo::listTips(SpecieTree))
  
  StatesFile="/Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/AROGER_TREE/A_H4_Complete_subset_linsi_op075_SpecieNameOnly.fa.statewithoutGapsPerNode"
  StatesH4=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  StatesH4.nd=StatesH4[which(duplicated(StatesH4$Node)==F),]
  NbOfSeauences=log10(StatesH4.nd[match(NodesCorresp,StatesH4.nd$Node),]$TotalNumberOfSeq09)
  NbOfSeauences[which(NbOfSeauences>6)]=6
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H4_NumberOfSeqslog10.tree")
  
  
  
  
  
  
  
  
  
  ###############################################
  #To visualise the Nb of ancestral possiblities
  # On the ancestral tree
  #Alternative approach, based on number of states with probability > 0.25
  ##############################################
  
  SpecieTree=ape::read.tree("/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")
  
  #H2A
  NodesCorresp=names(adephylo::listTips(SpecieTree))
  NodesCorresp=NodesCorresp[is.na(NodesCorresp)==F]
  StatesFile="/Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/AROGER_TREE/A_H2A_Complete_subset_linsi_op075_SpecieNameOnly.fa.statewithoutGapsPerNode"
  StatesH2A=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  
  #Compute the number of states > 0.15 for each position
  StatesH2A$States_Sup02=NA
  for(i in 1:dim(StatesH2A)[1]){
    Sel=as.list(StatesH2A[i,grep("p_",names(StatesH2A))])
    
    #Here, we want to count the number of states to take into account, arbitrarly chose as > 0.2 . 
    #If there are none, we choose 5 states (assuming that this is a lot of states to look at for just one position)
    StatesH2A[i,]$States_Sup02=ifelse(length(which(Sel>0.25))>0,length(which(Sel>0.25)),5)
  }
  
  #Now looping over each node, we multiple the number of states over threshold to compute how much sequences are need to have a clearer picture 
  StatesH2A$TotalNumberSeq02=NA
  for(u in unique(StatesH2A$Node)){
    StatesH2A[which(StatesH2A$Node==u),]$TotalNumberSeq02=prod(StatesH2A[which(StatesH2A$Node==u),]$States_Sup02)}
  
  StatesH2A.nd=StatesH2A[which(duplicated(StatesH2A$Node)==F),]
  NbOfSeauences=log10(StatesH2A.nd[match(NodesCorresp,StatesH2A.nd$Node),]$TotalNumberSeq02)
  
  #Putting a cap at 1M sequences
  if(length(which(NbOfSeauences>6))>0){
    NbOfSeauences[which(NbOfSeauences>6)]=6}
  
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H2A_Threshold025_approachlog10.tree")
  
  
  
  ################
  #H2A_Z
  ################
  
  
  SpecieTree=ape::read.tree("/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")
  NodesCorresp=names(adephylo::listTips(SpecieTree))
  NodesCorresp=NodesCorresp[is.na(NodesCorresp)==F]
  StatesFile="/Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/AROGER_TREE/A_H2A_Z_Complete_subset_linsi_op075_SpecieNameOnly.fa.statewithoutGapsPerNode"
  StatesH2A_Z=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  
  #Compute the number of states > 0.15 for each position
  StatesH2A_Z$States_Sup02=NA
  for(i in 1:dim(StatesH2A_Z)[1]){
    Sel=as.list(StatesH2A_Z[i,grep("p_",names(StatesH2A_Z))])
    
    #Here, we want to count the number of states to take into account, arbitrarly chose as > 0.2 . 
    #If there are none, we choose 5 states (assuming that this is a lot of states to look at for just one position)
    StatesH2A_Z[i,]$States_Sup02=ifelse(length(which(Sel>0.25))>0,length(which(Sel>0.25)),5)
  }
  
  #Now looping over each node, we multiple the number of states over threshold to compute how much sequences are need to have a clearer picture 
  StatesH2A_Z$TotalNumberSeq02=NA
  for(u in unique(StatesH2A_Z$Node)){
    StatesH2A_Z[which(StatesH2A_Z$Node==u),]$TotalNumberSeq02=prod(StatesH2A_Z[which(StatesH2A_Z$Node==u),]$States_Sup02)}
  
  StatesH2A_Z.nd=StatesH2A_Z[which(duplicated(StatesH2A_Z$Node)==F),]
  NbOfSeauences=log10(StatesH2A_Z.nd[match(NodesCorresp,StatesH2A_Z.nd$Node),]$TotalNumberSeq02)
  
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H2A_Z_Threshold025_approachlog10.tree")
  
  
  ##############
  #H2B
  ##############
  SpecieTree=ape::read.tree("/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")
  NodesCorresp=names(adephylo::listTips(SpecieTree))
  NodesCorresp=NodesCorresp[is.na(NodesCorresp)==F]
  StatesFile="/Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/AROGER_TREE/A_H2B_Complete_subset_linsi_op075_SpecieNameOnly.fa.statewithoutGapsPerNode"
  StatesH2B=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  
  #Compute the number of states > 0.15 for each position
  StatesH2B$States_Sup02=NA
  for(i in 1:dim(StatesH2B)[1]){
    Sel=as.list(StatesH2B[i,grep("p_",names(StatesH2B))])
    
    #Here, we want to count the number of states to take into account, arbitrarly chose as > 0.2 . 
    #If there are none, we choose 5 states (assuming that this is a lot of states to look at for just one position)
    StatesH2B[i,]$States_Sup02=ifelse(length(which(Sel>0.25))>0,length(which(Sel>0.25)),5)
  }
  
  #Now looping over each node, we multiple the number of states over threshold to compute how much sequences are need to have a clearer picture 
  StatesH2B$TotalNumberSeq02=NA
  for(u in unique(StatesH2B$Node)){
    StatesH2B[which(StatesH2B$Node==u),]$TotalNumberSeq02=prod(StatesH2B[which(StatesH2B$Node==u),]$States_Sup02)}
  
  StatesH2B.nd=StatesH2B[which(duplicated(StatesH2B$Node)==F),]
  NbOfSeauences=log10(StatesH2B.nd[match(NodesCorresp,StatesH2B.nd$Node),]$TotalNumberSeq02)
  
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H2B_Threshold025_approachlog10.tree")
  
  
  
  #H3_r
  
  
  ##############
  SpecieTree=ape::read.tree("/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")
  NodesCorresp=names(adephylo::listTips(SpecieTree))
  NodesCorresp=NodesCorresp[is.na(NodesCorresp)==F]
  StatesFile="/Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/AROGER_TREE/A_H3_r_Complete_subset_linsi_op075_SpecieNameOnly.fa.statewithoutGapsPerNode"
  StatesH3_r=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  
  #Compute the number of states > 0.15 for each position
  StatesH3_r$States_Sup02=NA
  for(i in 1:dim(StatesH3_r)[1]){
    Sel=as.list(StatesH3_r[i,grep("p_",names(StatesH3_r))])
    
    #Here, we want to count the number of states to take into account, arbitrarly chose as > 0.2 . 
    #If there are none, we choose 5 states (assuming that this is a lot of states to look at for just one position)
    StatesH3_r[i,]$States_Sup02=ifelse(length(which(Sel>0.25))>0,length(which(Sel>0.25)),5)
  }
  
  #Now looping over each node, we multiple the number of states over threshold to compute how much sequences are need to have a clearer picture 
  StatesH3_r$TotalNumberSeq02=NA
  for(u in unique(StatesH3_r$Node)){
    StatesH3_r[which(StatesH3_r$Node==u),]$TotalNumberSeq02=prod(StatesH3_r[which(StatesH3_r$Node==u),]$States_Sup02)}
  
  StatesH3_r.nd=StatesH3_r[which(duplicated(StatesH3_r$Node)==F),]
  NbOfSeauences=log10(StatesH3_r.nd[match(NodesCorresp,StatesH3_r.nd$Node),]$TotalNumberSeq02)
  
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H3_r_Threshold025_approachlog10.tree")
  
  
  
  #H4
  
  
  
  ##############
  SpecieTree=ape::read.tree("/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea.tree")
  NodesCorresp=names(adephylo::listTips(SpecieTree))
  NodesCorresp=NodesCorresp[is.na(NodesCorresp)==F]
  StatesFile="/Users/ahocher/Dropbox/Laboratory/ASGARD/HMMSEARCH/PARALOGS/FINAL/AROGER_TREE/A_H4_Complete_subset_linsi_op075_SpecieNameOnly.fa.statewithoutGapsPerNode"
  StatesH4=read.table(StatesFile,header=T,sep="\t",stringsAsFactors = F,quote="")
  
  #Compute the number of states > 0.15 for each position
  StatesH4$States_Sup02=NA
  for(i in 1:dim(StatesH4)[1]){
    Sel=as.list(StatesH4[i,grep("p_",names(StatesH4))])
    
    #Here, we want to count the number of states to take into account, arbitrarly chose as > 0.2 . 
    #If there are none, we choose 5 states (assuming that this is a lot of states to look at for just one position)
    StatesH4[i,]$States_Sup02=ifelse(length(which(Sel>0.25))>0,length(which(Sel>0.25)),5)
  }
  
  #Now looping over each node, we multiple the number of states over threshold to compute how much sequences are need to have a clearer picture 
  StatesH4$TotalNumberSeq02=NA
  for(u in unique(StatesH4$Node)){
    StatesH4[which(StatesH4$Node==u),]$TotalNumberSeq02=prod(StatesH4[which(StatesH4$Node==u),]$States_Sup02)}
  
  StatesH4.nd=StatesH4[which(duplicated(StatesH4$Node)==F),]
  NbOfSeauences=log10(StatesH4.nd[match(NodesCorresp,StatesH4.nd$Node),]$TotalNumberSeq02)
  
  SpecieTree$node.label=round(NbOfSeauences,2)
  
  ape::write.tree(SpecieTree,"/Users/ahocher/Dropbox/Laboratory/CAMBRIDGE/ANCESTRAL_NUC/TREE/ARogers_tree_Grafted_Eukprot_Specieonly_5HistoneSet.wArchaea_H4_Threshold025_approachlog10.tree")
  
  
  
  
