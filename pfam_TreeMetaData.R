#Aim : import data and export it as itol metadata, realted to Eukprot main tree



#import data

sig_pfam$TID=sig_pfam$Specie


OutDataPath="/Users/perrinm/Documents/AHocher_Project/ANALYSIS/pfam_correlations_itol.txt"

write.table(sig_pfam,"/Users/perrinm/Documents/AHocher_Project/ANALYSIS/pfam_correlations_itol.txt",row.names = F,sep="\t",quote=F)


setwd("/Users/perrinm/Documents/AHocher_Project/ITOL")
#import the function to export metadata
source("/Users/perrinm/Documents/SCRIPTS/FROM_OTHERS/table2itol-master/table2itol.R")
create_itol_files(OutDataPath,identifier = "TID",separator = "\t")
