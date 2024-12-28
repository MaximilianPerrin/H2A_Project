#Aim : import data and export it as itol metadata, realted to Eukprot main tree



#import data

H2ASummary=read.table("/Users/perrinm/Documents/AHocher_Project/ANALYSIS/H2APerSp_SQ.txt",header=T,sep="\t")


#
H2ASummary$TID=H2ASummary$Specie_Name


H2ASummary$binary_SQ=H2ASummary$SQ_Present

#binarise values
H2ASummary$binary_SQ[which(H2ASummary$binary_SQ>0)]=1

OutDataPath="/Users/perrinm/Documents/AHocher_Project/ANALYSIS/H2APerSp_SQ_export_itol.txt"

write.table(H2ASummary,"/Users/perrinm/Documents/AHocher_Project/ANALYSIS/H2APerSp_SQ_export_itol.txt",row.names = F,sep="\t",quote=F)


setwd("/Users/perrinm/Documents/AHocher_Project/ITOL")
#import the function to export metadata
source("/Users/perrinm/Documents/SCRIPTS/FROM_OTHERS/table2itol-master/table2itol.R")
create_itol_files(OutDataPath,identifier = "TID",separator = "\t")
