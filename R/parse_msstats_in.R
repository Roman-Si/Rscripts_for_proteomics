# Required Libraries
library(dplyr)
library(stringr)
library(Biostrings)
source("scripts/generic_functions.R")

# Function to extract sequence, number and indexes of quantified peptides for each protein 
# Maybe too long and needs splitting
#' @param msstats_input_path Path to msstats input style file from quantms (probably works for other programs as well?)
#' @param fasta_path Path to fasta file with protein sequences, can be gzipped
#' @param contaminant_prefix The prefix for contaminants, default is CONTAMINANT_
#' @param delimiter The delimiter in case of peptides shared between different proteins, default is ;
#' @param output_file Name of the output file (optional). better to merge it with the output of parse_mzTab::extract_protein_stats first and save together
#' @return Data frame with peptide info (without PTMs) for each detected protein with following columns
# - peptide sequences
# - count of unmodified peptides
# - min and max index of identified peptides in protein sequence, useful to check which part of the protein is detected

# Define the function
extract_peptides_per_protein <- function(msstats_input_path, fasta_path, contaminant_prefix = "CONTAMINANT", delimiter = ";", output_file = NULL) {
  
  # 1. Read msstats_input, remove the PTM info from peptides
  msstats_input <- read.csv(msstats_input_path, header = TRUE, sep = ',')
  msstats_input <- msstats_input %>%
    filter(!grepl(contaminant_prefix, ProteinName)) %>%
    mutate(PeptideSequence = str_replace_all(PeptideSequence, "\\(.*?\\)", "") %>% # Remove modifications
            str_replace_all("\\.", "")) # Remove dots
  
  # 2. Read fasta file and keep only detected proteins
  protein_ids <- unlist(strsplit(msstats_input$ProteinName, ",")) %>% unique()
  fasta_file <- readAAStringSet(fasta_path, format = "fasta")
  fasta_file <- fasta_file[names(fasta_file) %in% protein_ids]

  
  # 3. create the peptide_df and iterate over each protein
  peptide_df <- data.frame(accession = protein_ids, peptide_seqs = NA, nr_unmodified_peptides = NA,  peptide_index = NA, stringsAsFactors = FALSE)

  for (i in seq_along(fasta_file)) {
    protein_id <- names(fasta_file)[i]
    protein_seq <- as.character(fasta_file[[i]])
    start_index <- +Inf
    end_index <- 0
    
    # Filter for peptides of the current protein
    peptides <- msstats_input %>%
      filter(grepl(protein_id, ProteinName)) %>%
      pull(PeptideSequence) %>%
      unique()
    
    # Iterate peptides and update indices
    for (pep in peptides) {
      pep_indices <- str_locate(protein_seq, pep)
      
      if (!is.na(pep_indices[1])) { # checks if peptide is found, unnecessary check probably
        start_index <- min(start_index, pep_indices[1])
        end_index <- max(end_index, pep_indices[2])
      }
    }
    
    # If no valid peptides are found, set indices to NA
    if (start_index == +Inf || end_index == 0) {
      start_index <- NA
      end_index <- NA
    }
    
    # Update the peptide_df
    peptide_df <- peptide_df %>%
      mutate(
        peptide_seqs = ifelse(accession == protein_id, paste(peptides, collapse = delimiter), peptide_seqs),
        nr_unmodified_peptides = ifelse(accession == protein_id, length(peptides), nr_unmodified_peptides),
        peptide_index = ifelse(accession == protein_id, paste0(start_index, "-", end_index), peptide_index)
      )
  }
  
  # if save to csv
  if (!is.null(output_file)) {
    save_gzipped_csv(peptide_df, output_file)
  }
  
  return(peptide_df)
}