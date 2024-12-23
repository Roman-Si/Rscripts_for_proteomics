#' @import dplyr
#' @importFrom stringr str_replace_all str_locate
#' @importFrom Biostrings readAAStringSet
NULL

#' Extract peptide sequences, count and indices per protein
#'
#' Processes the msstats_in (formatted as input for msstats) file to extract peptide sequences, number and indexes of quantified peptides for each protein. Maybe too long and needs splitting
#'
#' @param msstats_input_path Path to msstats input style file from quantms (probably works for other programs as well?)
#' @param fasta_path Path to fasta file with protein sequences (can be gzipped)
#' @param contaminant_prefix Prefix for contaminants (default is CONTAMINANT)
#' @param delimiter Delimiter for peptides shared between different proteins (default is ;)
#' @param output_file Optional. Name of the output file (better to merge it with the output of parse_mzTab::extract_protein_stats first and save together)
#' @return Data frame with peptide info (without PTMs) for each detected protein including:
# - `peptide_seqs`: ; delimited peptide sequences
# - `nr_unmodified_peptides`: count of unmodified peptides
# - `peptide_index`: min and max index of identified peptides in protein sequence, useful to check which part of the protein is detected
#' @details This function performs the following steps:
#'   - Reads the msstats input file and removes PTM information from peptide sequences
#'   - Filters the FASTA file to include only detected proteins from msstats_in
#'   - Extracts peptide sequences, counts, and positional indices for each protein.
#' @export
extract_peptides_per_protein <- function(msstats_input_path, fasta_path, contaminant_prefix = "CONTAMINANT", delimiter = ";", output_file = NULL) {
  
  # Helper function to preprocess msstats input
  preprocess_msstats_input <- function(msstats_input_path, contaminant_prefix) {
    msstats_input <- read.csv(msstats_input_path, header = TRUE, sep = ',')
    msstats_input <- msstats_input %>%
      filter(!grepl(contaminant_prefix, ProteinName)) %>%
      mutate(PeptideSequence = stringr::str_replace_all(PeptideSequence, "\\(.*?\\)", "") %>% # Remove modifications
               stringr::str_replace_all("\\.", "")) # Remove dots
    return(msstats_input)
  }

  # Helper function to filter the FASTA file for detected proteins (move it in general functions?)
  filter_fasta_file <- function(fasta_path, protein_ids) {
    fasta_file <- readAAStringSet(fasta_path, format = "fasta")
    fasta_file <- fasta_file[names(fasta_file) %in% protein_ids]
    return(fasta_file)
  }
  
  # Read and preprocess inputs
  msstats_input <- preprocess_msstats_input(msstats_input_path, contaminant_prefix)
  protein_ids <- unlist(strsplit(msstats_input$ProteinName, delimiter)) %>% unique()
  fasta_file <- filter_fasta_file(fasta_path, protein_ids)

  # Initialize peptide_df
  peptide_df <- data.frame(accession = protein_ids, nr_unmodified_peptides = NA,  peptide_index = NA, peptide_seqs = NA, stringsAsFactors = FALSE)

  # Iterate over proteins in fasta
  for (i in seq_along(fasta_file)) {
    protein_id <- names(fasta_file)[i]
    protein_seq <- as.character(fasta_file[[i]])
    start_index <- NA
    end_index <- NA
    peptides <- msstats_input %>%
      filter(grepl(protein_id, ProteinName)) %>%
      pull(PeptideSequence) %>%
      unique()
    

    # Find peptide indices
    for (pep in peptides) {
      pep_indices <- stringr::str_locate(protein_seq, pep)
      start_index <- min(start_index, pep_indices[1], na.rm = TRUE)
      end_index <- max(end_index, pep_indices[2], na.rm = TRUE)
      }

    # Update peptide_df
    peptide_df[peptide_df$accession == protein_id, ] <- list(
      accession = protein_id,
      nr_unmodified_peptides = length(peptides),
      peptide_index = paste0(start_index, "-", end_index),
      peptide_seqs = paste(peptides, collapse = delimiter)
      
    )
  }
  
  
  # if save to csv
  if (!is.null(output_file)) {
    save_gzipped_csv(peptide_df, output_file)
  }
  
  return(peptide_df)
}