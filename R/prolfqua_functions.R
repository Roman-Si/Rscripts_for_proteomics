#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr separate_rows
#' @import prolfqua
NULL


#' Create an LFQData object for prolfqua from msstats_in
#'
#' Converts a data frame of msstats input format or a custom data frame to LFQData object for the prolfqua package. Sample names in must be in Reference column, response in Intensity column
#'
#' @param data Dataframe with peptide at msstats_in format or protein data
#' @param contaminant_prefix Prefix for contaminants, default is CONTAMINANT_
#' @param response_level Specifies the level to use:
#'   - `"peptide"` (default): Aggregates peptide-level intensities.
#'   - `"protein"`: Uses protein-level intensities.
#' @param proteinId_column Name of column with protein IDs (deafult is "ProteinName")
#' @param factor_column Name of column with the experimental/biological factor (default is "Condition")
#' @param extra_factor Optional. Name of additional factor for protein level data, batches for example.
#' @return The LFQData object
#' @export
create_lfqdata <- function(data, contaminant_prefix = 'CONTAMINANT_', response_level = "peptide", proteinId_column = "ProteinName",factor_column = "Condition", extra_factor = NULL) {

    # check the response level is correct
    if (!response_level %in% c("peptide", "protein")) {
        stop("Invalid response_level. Choose either 'peptide' or 'protein'.")
        }

    if (response_level == "peptide") {
        data <- data |> dplyr::group_by(Reference, Run, Condition, BioReplicate, PeptideSequence, ProteinName) |> 
        dplyr::summarize(nprecursor = dplyr::n(),Intensity = sum(Intensity, na.rm = TRUE), .groups = "drop")
        }

    atable <- prolfqua::AnalysisTableAnnotation$new()
    atable$fileName = "Reference"
    atable$hierarchy[["protein_Id"]] <- c(proteinId_column)
    if (response_level == "peptide") {
        atable$hierarchy[["peptide_Id"]] <- c("PeptideSequence")
        }    
    atable$hierarchyDepth <- 1
    atable$set_response("Intensity")
    atable$factors[["condition_"]] = factor_column
    atable$factors[["run"]] = "Run"
    if(!is.null(extra_factor)) {
        atable$factors[["extrafactor_"]] = extra_factor
}
    atable$factorDepth <- 1
    config <- prolfqua::AnalysisConfiguration$new(atable)
    
    adata <- prolfqua::setup_analysis(data, config)
    lfqdata <- prolfqua::LFQData$new(adata, config)
    lfqdata$data <- lfqdata$data %>% filter(!grepl("^DECOY_", protein_Id)) %>% filter(!grepl(contaminant_prefix, protein_Id))
    lfqdata$remove_small_intensities()
    
    return(lfqdata)
  
}


#' Aggregate peptide intensities to protein level
#'
#' Function for protein aggregation of peptide normalized data
#' - Use Tukey's median polish protein aggregation
#' - Remove artifacts with same value in all samples 
#' - Replace the artifact intensities with top3 method
#' - Merge the two methods
#'
#' @param lfqdataPeptideNorm The prolfqua object after peptide normalization
#' @return The a dataframe with protein intensities
#' @export
protein_aggregation <- function(lfqdataPeptideNorm) {
    
    # Helper function to check if all values in a row are the same
    all_same <- function(x) {
        unique_values <- unique(na.omit(x))
        length(unique_values) == 1}

    # TMP
    agg <- lfqdataPeptideNorm$get_Aggregator()
    agg$medpolish()
    mp_df_tmp <- data.frame(agg$lfq_agg$to_wide()$data)
    mp_df_tmp$isotopeLabel <- NULL
    # Identify and filter proteins with the same values across all samples
    same_values_proteins <- mp_df_tmp %>%
        dplyr::mutate(same = apply(.[, -1], 1, all_same)) %>%
        dplyr::filter(same) %>%
        dplyr::pull(protein_Id)
    mp_df_tmp_filtered <- mp_df_tmp %>%
        dplyr::filter(!protein_Id %in% same_values_proteins)

    # top3
    agg <- lfqdataPeptideNorm$get_Aggregator()
    agg$mean_topN(3)
    mp_df_top3 <- data.frame(agg$lfq_agg$to_wide()$data)
    mp_df_top3$isotopeLabel <- NULL

    mp_df_top3_filtered <- mp_df_top3 %>%
        dplyr::filter(protein_Id %in% same_values_proteins)
  
    # Merge both filtered datasets
    mp_df <- rbind(mp_df_tmp_filtered, mp_df_top3_filtered)
  
    return(mp_df)
}


#' Filter protein groups to retain leading protein IDs
#'
#' For protein groups with multiple IDs, filters the data frame to retain only the leading protein ID based on a given list.
#'
#' @param mp_df A dataframe with proteinGroups in `protein_Id` column
#' @param proteomicslfq_ids_to_keep A vector containing the proteomicsLFQ leading proteins
#' @param delimiter Delimiter for multiple proteins in a protein group (default is ;)
#' @return A dataframe with only leading protein in `protein_Id` column
#' @export
filter_for_leading_protein <- function(mp_df, proteomicslfq_ids_to_keep, protein_separator = ";") {
    proteins_in_groups <- mp_df %>%
        dplyr::filter(grepl(protein_separator, protein_Id)) %>%
        tidyr::separate_rows(protein_Id, sep = protein_separator) %>%
        dplyr::pull(protein_Id)
        
    if (length(proteins_in_groups) > 0) { # check there are proteinGroups
        mp_df_multiple_ids <- mp_df %>%
        dplyr::filter(grepl(protein_separator, protein_Id))
        # no mapply iterates over each element of mp_df_multiple_ids$protein_Id and applies to it get_valid_ids with proteomicslfq_ids_to_keep as fixed argument
        mp_df_multiple_ids$protein_Id <- mapply(get_valid_ids, mp_df_multiple_ids$protein_Id, MoreArgs = list(proteomicslfq_ids_to_keep))
    
        # Merge the filtered rows with single IDs and those with previously multiple IDs that were processed
        mp_df <- dplyr::bind_rows(
            mp_df %>% dplyr::filter(!grepl(';', protein_Id)), 
            mp_df_multiple_ids)
        }
  
    return(mp_df)
}


#' Filter proteins by replicate counts per condition
#'
#' Retains proteins that are quantified in at least a specified number of replicates for one condition.
#'
#' @param mp_df The dataframe with protein quantifications
#' @param conditions The vector of the conditions to search the mp_df colnames
#' @param min_replicates Minimum number of min replicates needed
#' @return The dataframe with proteins passing the filtering
#' @export
filter_proteins_by_replicates <- function(mp_df, conditions, min_replicates) {
    quantified_proteins <- c()
  
    for (condition in conditions) {
        condition_cols <- grep(condition, colnames(mp_df), value = TRUE)
        max_nan <- length(condition_cols) - min_replicates
        mp_df[[paste0("count_na_", condition)]] <- rowSums(is.na(mp_df[, condition_cols]))

        passing_proteins <- rownames(mp_df[mp_df[[paste0("count_na_", condition)]] <= max_nan, ])
        quantified_proteins <- unique(c(quantified_proteins, passing_proteins))
        }
    
    mp_df_quant <- mp_df[rownames(mp_df) %in% quantified_proteins, ]
    mp_df_quant <- mp_df_quant %>% select(-contains('count_na_'))
    
    return(mp_df_quant)
}
