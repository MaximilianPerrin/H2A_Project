# Clear the workspace
rm(list = ls())

# Set the directory containing the files
input_directory <- "/Users/perrinm/Documents/AHocher_Project/ITOL/ITOL_pfam_correlations/"
output_directory <- "/Users/perrinm/Documents/AHocher_Project/ITOL/ITOL_pfam_correlations/ITOL_Fixed/"

# Ensure the output directory exists
if (!dir.exists(output_directory)) {
  dir.create(output_directory)
}

# List all files in the directory
files <- list.files(path = input_directory, full.names = TRUE)

# Loop through each file, read it, apply the transformations, and save it
for (file in files) {
  # Load the data, assuming tab-separated format
  data <- read.delim(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Identify lines where the last element is 0, and replace it with -1
  data[data[, ncol(data)] == 0, ncol(data)] <- -1
  
  # Apply the specific replacements for species names
  data <- as.data.frame(lapply(data, function(x) gsub("Acanthocystis_sp_HF_20", "Acanthocystis", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("Apusomonadida_sp_AF_17", "Apusomonadida_sp_AF-17", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("Blastocystis_sp_subtype4_isolateWR1", "Blastocystis_sp_subtype4-isolateWR1", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("Colponemidia_sp_Colp_10", "Colponemidia_sp_Colp-10", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("BColponemidia_sp_Colp_15", "Colponemidia_sp_Colp-15", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("MAST_03A_sp_MAST_3A_sp1", "MAST-03A_sp_MAST-3A-sp1", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("MAST_04A_sp_MAST_4A0_sp1", "MAST-04A_sp_MAST-4A0-sp1", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("MAST_04C_sp_MAST_4C_sp1", "MAST-04C_sp_MAST-4C-sp1", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("Choanocystis_sp_HF_7", "Choanocystis", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("Telonema_sp_P_2", "Telonema_sp_P-2", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("Anaeramoeba_flamelloides", "AnaeramoebaFlam", x)))
  data <- as.data.frame(lapply(data, function(x) gsub("Hematodinium_sp_SG_2012", "Hematodinium_sp_SG-2012", x)))
  
  # Define the output file name
  output_file <- paste0(output_directory, basename(tools::file_path_sans_ext(file)), "_Fixed.txt")
  
  # Save the updated data to the output file
  write.table(data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
