#' @importFrom Biostrings readAAStringSet AAStringSet
#' @importFrom cleaver cleave
NULL

#' Read protein FASTA and parse headers
#'
#' Reads a FASTA file as an AAStringSet and simplifies the sequence names:
#'   - Removes everything after the first space.
#'   - If UniProt-style headers, extracts the accession.
#'
#' @param fasta_path Path to the FASTA file
#' @param uniprot_fasta_header Logical. If TRUE, parse UniProt-style headers.
#'
#' @return Biostrings AAStringSet with parsed names
#' @export
read_and_parse_fasta <- function(fasta_path, uniprot_fasta_header = FALSE) {
  fasta_file <- Biostrings::readAAStringSet(fasta_path, format = "fasta")
  
  # Remove everything after first space
  names(fasta_file) <- sub("(.*?) .*", "\\1", names(fasta_file))
  
  # If UniProt headers, extract the accession
  if (uniprot_fasta_header) {
    names(fasta_file) <- sub(".*\\|(.*)\\|.*", "\\1", names(fasta_file))
  }
  
  fasta_file
}

#' Filter FASTA with IDs
#'
#' @param fasta_path Path to the FASTA file
#' @param protein_ids Vector of IDs to keep
#' @param uniprot_fasta_header Logical. If TRUE, parse UniProt-style headers.
#'
#' @return Filtered AAStringSet
#' @export
filter_fasta_file_with_ids <- function(fasta_path, protein_ids, uniprot_fasta_header = FALSE) {
  fasta_file <- read_and_parse_fasta(fasta_path, uniprot_fasta_header)
  fasta_file[names(fasta_file) %in% protein_ids]
}

#' Generate pseudo-reversed decoy sequences as Biostrings AAStringSet
#'
#' Keeps the first residue fixed (methionine) and reverses remaining sequence.
#'
#' @param fasta_path Path to the FASTA file
#' @param uniprot_fasta_header Logical. If TRUE, parse UniProt-style headers.
#'
#' @return AAStringSet of decoy sequences prefixed by "DECOY_"
#' @export
make_decoy_AAStringSet <- function(fasta_path, uniprot_fasta_header = FALSE) {
  fasta_file <- read_and_parse_fasta(fasta_path, uniprot_fasta_header)
  
  decoy_seqs <- vapply(as.character(fasta_file), function(seq) {
    n <- nchar(seq)
    if (n <= 1) {
      # if for some mysterious reason there are sequence of length 1 or empty: just return it
      return(seq)
    }
    first <- substr(seq, 1, 1)
    rest  <- substr(seq, 2, n)
    rest_rev <- paste(rev(strsplit(rest, "", fixed = TRUE)[[1]]), collapse = "")
    paste0(first, rest_rev)
  }, character(1))
  
  decoys <- Biostrings::AAStringSet(decoy_seqs)
  names(decoys) <- paste0("DECOY_", names(fasta_file))
  decoys
}


#' Cleave an AAStringSet into peptides
#'
#' @param aa_set Biostrings AAStringSet of protein sequences
#' @param enzyme_cleaver The enzyme rule used for digestion with cleaver:cleave (default "trypsin-simple" allowing cuts before proline)
#' @param mc Maximum number of missed cleavages (default 2)
#' @param return_peptide_list_only Logical. If TRUE, return a character vector
#'   of all peptides (unlisted, protein names dropped). If FALSE (default), return the
#'   list returned by cleaver::cleave.
#'
#' @return Either a list (default) of peptides per sequence, or a flat character
#'   vector with peptide sequences if `return_peptide_list_only = TRUE`.
#' @export
cleave_AAStringSet <- function(aa_set, enzyme_cleaver = "trypsin-simple", mc = 2, return_peptide_list_only = FALSE) {
  pep_list <- cleaver::cleave(as.character(aa_set),
                              enzym = enzyme_cleaver,
                              missedCleavages = 0:mc)
  
  if (return_peptide_list_only) {
    return(unlist(pep_list, use.names = FALSE))
  } else {
    return(pep_list)
  }
}

