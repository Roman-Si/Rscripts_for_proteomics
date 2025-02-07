#' @import dplyr
#' @importFrom arrow read_parquet
NULL

#' Convert diann parquet to msstats_in format
#'
#' Converts the diann parquet to msstats_in format after filtering for contaminants and qvalue of precursors, then renames columns to prepare dataframe for functions like the extract_peptides_per_protein
#'
#' @param diann_parquet Path to diann parquet main report
#' @param qvalue FDR for precursors (default 0.01)
#' @param contaminant_prefix Prefix for contaminants (default is Cont)
#' @return Data frame with reliable precursors and without contaminants and the following columns for msstats_input format
# - `ProteinName`: derived from the Protein.Group column
# - `PeptideSequence`: derived from the Stripped.Sequence column
#' @export
convert_diann_to_msstats_in <- function(diann_parquet, qvalue = 0.01, contaminant_prefix = "Cont") {
    diann_df <- arrow::read_parquet(diann_parquet)
    diann_df <- diann_df %>%
        filter(!grepl(contaminant_prefix, Protein.Ids)) %>%
        filter(Q.Value <= 0.01)
    diann_df$ProteinName <- diann_df$Protein.Group
    diann_df$PeptideSequence <- diann_df$Stripped.Sequence
  return(diann_df)
}


#' Prepare diann for prolfqua
#'
#' Converts the diann parquet to a protein-level dataframe that will be converted to prolfqua object with create_lfqdata function
#'
#' @param diann_parquet Path to diann parquet main report
#' @param contaminant_prefix Prefix for contaminants (default is Cont)
#' @param response_level Specifies the level to use:
#'   - `"protein"` (default): Will use the proteinGroup level MaxLFQ
#'   - `"gene"`: Will use the gene level MaxLFQ
#' @param qvalue FDR for precursors and response level (default 0.01)
#' @return Data frame with reliable precursors and without contaminants and the following columns for msstats_input format
# - `accession`: derived from the Protein.Group column
# - `PeptideSequence`: derived from the Stripped.Sequence column
#' @export
prepare_diann_for_prolfqua <- function(diann_parquet, qvalue = 0.01, response_level = "protein", contaminant_prefix = "Cont") {
    diann_df <- arrow::read_parquet(diann_parquet)

    # Filter for contaminants and precursors 1% FDR
    diann_df <- diann_df %>%
        filter(!grepl(contaminant_prefix, Protein.Ids)) %>%
        filter(Q.Value <= 0.01)
        
    # check the response level is correct
    if (!response_level %in% c("protein", "gene")) {
        stop("Invalid response_level. Choose either 'protein' or 'gene'.")
        }

    if (response_level == "protein") {
        quant_df <- diann_df %>% 
                        filter(PG.Q.Value <= 0.01) %>% 
                        select(Run, Protein.Group, PG.MaxLFQ) %>%
                        distinct()
        }
    
    if (response_level == "gene") {
        quant_df <- diann_df %>% 
                        filter(GG.Q.Value <= 0.01) %>% 
                        select(Run, Genes, Genes.MaxLFQ) %>% 
                        distinct()
        }
    colnames(quant_df) <- c("Reference", "ProteinName", "Intensity")
    
  return(quant_df)
}
