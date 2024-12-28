#Clear environmment
rm(list = ls())
#Install packages

library(seqinr)
library(pheatmap)
library(tidyverse)
#The aim of this code is to produce a heatmap of protein domain corelation
#Code must be run with Protein_Domain_Correlation.R

#Set working directory
setwd("/Users/perrinm/Documents/AHocher_Project/PROTEIN_DOMAINS/")

#Read CSV
pfam_dom <- read.table("pfam_dom_No_NA.csv", header = TRUE, row.names = 1, sep = ",")
pfam_dom <- select(pfam_dom, !OriginalID)
pfam_dom <- pfam_dom %>% select(any_of(names(sig_pfam)), H2A_SQ, GPR_Gpa2_C, Las1)

#calcualte correlation matrix
cor_mat <- cor(pfam_dom)
cor_mat <- data.frame(cor_mat) %>% select(any_of(names(sig_pfam)), H2A_SQ)
  
#Visulise in heatmap
pheatmap(cor_mat, fontsize_row = 5, fontsize_col = 5, main = "PFAM Correlations" )

#Save the heatmap diagram as a .pdf

my_pheatmap <- pheatmap(cor_mat, fontsize_row = 5, fontsize_col = 5, main = "PFAM Corelations" )[[4]]
my_plot <- grid.arrange(my_pheatmap, nrow=1, ncol=1)
ggsave(filename="my_pheatmap.svg", plot=my_plot)