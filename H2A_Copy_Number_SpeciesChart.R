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
contains_SQ <- sapply(AllH2A, function(seq) grepl("sq(e|d)|SQ(E|D)", paste(seq, collapse = "")))

# Step 2: Create a data frame with the names of the proteins
H2A_table <- data.frame(OriginalName = names(AllH2A))

# Step 3: Create another column with only the name of the species
H2A_table$Specie_Name <- gsub(".*\\$(.*?)@.*", "\\1", H2A_table$OriginalName)

# Step 4: Count H2As with and without SQ
H2A_table$Contains_SQ <- contains_SQ


      
# Step 5: Summarize counts
library(dplyr)

H2APerSp_SQ <- H2A_table %>%
  group_by(Specie_Name) %>%
  summarize(
    Total_H2A_Sequences = n(),                      # Total sequences per species
    Contains_SQ = ifelse(sum(Contains_SQ) > 0, 1, 0)  # 1 if any sequence contains SQ, otherwise 0
  )

print(H2APerSp_SQ)

# Step 1: Count H2As with and without the SQ motif per species
H2APerSp_SQ <- H2A_table %>%
  group_by(Specie_Name) %>%
  summarize(
    Total_H2A_Sequences = n(),                 # Total sequences per species
    SQ_Present = sum(Contains_SQ),             # Number of sequences containing SQ per species
    SQ_Absent = n() - sum(Contains_SQ)         # Number of sequences not containing SQ
  )

# Step 2: Reshape the data for ggplot (melt it to long format)
library(tidyr)
H2A_long <- H2APerSp_SQ %>%
  pivot_longer(cols = c(SQ_Present, SQ_Absent), 
               names_to = "SQ_Status", 
               values_to = "Count")

# Step 3: Create a stacked bar chart
ggplot(H2A_long, aes(x = Specie_Name, y = Count, fill = SQ_Status)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "H2A Copy Number per Species with SQ Distinction",
       x = "Species",
       y = "Number of H2A Copies",
       fill = "SQ Motif") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

