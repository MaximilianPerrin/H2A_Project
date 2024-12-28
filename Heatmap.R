#Clear environmment
rm(list = ls())
#Install packages
library(seqinr)
library(pheatmap)
library(tidyverse)

#Set working directory
setwd("/Users/perrinm/Documents/AHocher_Project/PROTEIN_DOMAINS/")

#Read CSV
pfam_dom <- read.table("pfam_dom_No_NA.csv", header = TRUE, row.names = 1, sep = ",")
pfam_dom <- select(pfam_dom, !OriginalID)
pfam_dom <- pfam_dom %>% select(any_of(names(sig_pfam)), H2A_SQ)

#calcualte correlation matrix
cor_mat <- cor(pfam_dom)
print(cor_mat)

#Visulise in heatmap
pheatmap(cor_mat, fontsize_row = 5, fontsize_col = 5, main = "PFAM Corelations" )

#Save the heatmap diagram as a .pdf

my_pheatmap <- pheatmap(cor_mat, fontsize_row = 5, fontsize_col = 5, main = "PFAM Corelations" )[[4]]
my_plot <- grid.arrange(my_pheatmap, nrow=1, ncol=1)
ggsave(filename="my_pheatmap.svg", plot=my_plot)