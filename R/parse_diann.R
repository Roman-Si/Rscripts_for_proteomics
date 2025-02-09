#' @import dplyr
#' @importFrom arrow read_parquet
#' @importFrom tidyr unnest
NULL

#' Filter diann parquet for 1% precursor FDR and contaminants
#'
#' @param diann_parquet Path to diann parquet main report
#' @param qvalue FDR for precursors (default 0.01)
#' @param contaminant_prefix Prefix for contaminants (default is Cont)
#' @return Data frame with reliable precursors and without contaminants
#' @export
filter_diann_parquet <- function(diann_parquet, qvalue = 0.01, contaminant_prefix = "Cont") {
    diann_df <- arrow::read_parquet(diann_parquet)
    diann_df <- diann_df %>%
        filter(!grepl(contaminant_prefix, Protein.Ids)) %>%
        filter(Q.Value <= 0.01)
  return(diann_df)
}

#' Count PSMs per proteinGroup or gene in diann_dataframe
#'
#' Extracts PSM count for each proteinGroup or gene in the diann_dataframe. PENDING to add peptidoform count, not urgent since I do not include variable PTMs for now
#'
#' @param diann_dataframe The diann dataframe, prefiltered with filter_diann_parquet
#' @param response_level Specifies the level to count:
#'   - `"protein"` (default): Will count peptidoforms and PSMs per proteinGroup
#'   - `"gene"`: Will count peptidoforms and PSMs per gene
#' @return Data frame with the following columns:
#'   - `Protein.Group or Gene`: Protein or gene accession, one per row
#'   - `PSMs`: Number of PSMs.
#' @export
count_psms_diann <- function(diann_dataframe,  response_level = "protein") {
    # check the response level is correct
    if (!response_level %in% c("protein", "gene")) {
        stop("Invalid response_level. Choose either 'protein' or 'gene'.")
        }

    if (response_level == "protein") {
        psms_count <- diann_dataframe %>% 
                        mutate(Protein.Group = strsplit(as.character(Protein.Group), ";")) %>% # ; separated if multiple proteins per group
                        unnest(Protein.Group) %>%
                        distinct() %>%
                        group_by(Protein.Group) %>%
                        summarize(PSMs = n())
        }
    
    if (response_level == "gene") {
        psms_count <- diann_dataframe %>% 
                        mutate(Genes = strsplit(as.character(Genes), ";")) %>% # ; separated if multiple genes per group
                        unnest(Genes) %>%
                        distinct() %>%
                        group_by(Genes) %>%
                        summarize(PSMs = n())
        }
 
  return(psms_count)
  
}




#' Rename columns of diann dataframe to match msstats_in names
#'
#' Converts the diann parquet to msstats_in format by renaming columns to prepare dataframe for functions like the extract_peptides_per_protein
#'
#' @param diann_dataframe The diann dataframe, prefiltered with filter_diann_parquet
#' @return Data frame with reliable precursors and without contaminants and the following columns for msstats_input format
# - `ProteinName`: derived from the Protein.Group column
# - `PeptideSequence`: derived from the Stripped.Sequence column
#' @export
rename_diann_cols_like_msstats_in <- function(diann_dataframe) {
    diann_dataframe$ProteinName <- diann_dataframe$Protein.Group
    diann_dataframe$PeptideSequence <- diann_dataframe$Stripped.Sequence
  return(diann_dataframe)
}


#' Prepare diann dataframe for prolfqua
#'
#' Converts the diann dataframe to a protein-level dataframe that will be converted to prolfqua object with create_lfqdata function
#'
#' @param diann_dataframe The diann dataframe, prefiltered with filter_diann_parquet
#' @param response_level Specifies the level to use:
#'   - `"protein"` (default): Will use the proteinGroup level MaxLFQ
#'   - `"gene"`: Will use the gene level MaxLFQ
#' @param qvalue FDR for response level (default 0.01)
#' @return Data frame with reliable precursors and without contaminants and the following columns for msstats_input format
# - `accession`: derived from the Protein.Group column
# - `PeptideSequence`: derived from the Stripped.Sequence column
#' @export
prepare_diann_for_prolfqua <- function(diann_dataframe, qvalue = 0.01, response_level = "protein") {        
    # check the response level is correct
    if (!response_level %in% c("protein", "gene")) {
        stop("Invalid response_level. Choose either 'protein' or 'gene'.")
        }

    if (response_level == "protein") {
        quant_df <- diann_dataframe %>% 
                        filter(PG.Q.Value <= 0.01) %>% 
                        select(Run, Protein.Group, PG.MaxLFQ) %>%
                        distinct()
        }
    
    if (response_level == "gene") {
        quant_df <- diann_dataframe %>% 
                        filter(GG.Q.Value <= 0.01) %>% 
                        select(Run, Genes, Genes.MaxLFQ) %>% 
                        distinct()
        }
    colnames(quant_df) <- c("Reference", "ProteinName", "Intensity")
    
  return(quant_df)
}
