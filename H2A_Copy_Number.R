install.packages("seqinr")
library(seqinr)

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

#Install ggplot2 to make a chart
install.packages("ggplot2")
library(ggplot2)


ggplot(H2APerSpecie, aes(x = Nb_H2A)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "H2A Copy number per species",
       x = "H2A Copy Number",
       y = "Species")


# Step 2: Check for the "SQ" pattern in each sequence
contains_SQ <- sapply(AllH2A, function(seq) grepl("SQ", paste(seq, collapse = "")))

# Step 3: Split sequences into two groups
sequences_with_SQ <- AllH2A[contains_SQ]     # Sequences that contain "SQ"
sequences_without_SQ <- AllH2A[!contains_SQ] # Sequences that do not contain "SQ"

# Print the results: Sequence names in each group
cat("Sequences that contain 'SQ':\n")
print(names(sequences_with_SQ))

cat("\nSequences that do not contain 'SQ':\n")
print(names(sequences_without_SQ))
