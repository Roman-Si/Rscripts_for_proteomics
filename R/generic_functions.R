#' @import dplyr
#' @importFrom stringr str_replace_all str_locate
#' @importFrom Biostrings readAAStringSet
NULL


#' Save gzipped csv file
#'
#' Saves a dataframe to a gzipped csv file, optionally including row names. . If the specified directory does not exist, it will be created.
#'
#' @param df The dataframe
#' @param output_file The name of the output fle, if given in a directory that does not exist it will be created
#' @param row_names Logical; whether to include row names in the output (default is FALSE)
#' @return None
#' @export
save_gzipped_csv <- function(df, output_file, row_names = FALSE) {
  output_dir <- dirname(output_file)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(df, gzfile(output_file), row.names = row_names)
}



#' Get Valid Protein IDs
#'
#' Filter vector of IDs, retaining only main IDs  (as determined from proteomicsLFQ mzTab) in case of proteinGroups
#'
#' @param all_ids A vector with one or multiple IDs per item
#' @param valid_ids A vector with the IDs to keep
#' @param delimiter The delimiter in case of multiple IDs per row (default is ;)
#' @return The filtered vector with one leading ID per group
#' @export
get_valid_ids <- function(all_ids, valid_ids, delimiter = ";") {
  ids_split <- strsplit(all_ids, delimiter)[[1]]
  valid_ids <- ids_split[ids_split %in% valid_ids]
  paste(valid_ids, collapse = ";")
}


#' Filter FASTA with IDs
#'
#' Filter a FASTA file for proteins present in a vector
#'
#' @param fasta_path Path to the FASTA file
#' @param protein_ids A vector with the IDs to keep
#' @param uniprot_fasta_header If FASTA headers are in UniProt format parses the accession (deafult FALSE)
#' @return The filtered Biostrings XStringSet object
#' @export
filter_fasta_file_with_ids <- function(fasta_path, protein_ids, uniprot_fasta_header = FALSE) {
    fasta_file <- readAAStringSet(fasta_path, format = "fasta")

    # Remove everything after first space
    names(fasta_file) <- sub("(.*?) .*", "\\1", names(fasta_file))

    if (uniprot_fasta_header) {
    names(fasta_file) <- sub(".*\\|(.*)\\|.*", "\\1", names(fasta_file))
    }

    fasta_file <- fasta_file[names(fasta_file) %in% protein_ids]
    return(fasta_file)
  }