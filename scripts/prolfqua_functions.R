library(tidyverse)
library(ggplot2)
library(prolfqua)
source("scripts/generic_functions.R")

# Function to turn the msstats input dataframe or a cutom dataframe to lfqdata object for prolfqua 
# Sample names in Reference column, response in Intensity column
#' @param data Dataframe with peptide at msstats_in format or protein data
#' @param contaminant_prefix The prefix for contaminants, default is CONTAMINANT_
#' @param response_level Default is "peptide" and assumes msstats_in column names, set to "protein" to use protein level intensities
#' @param proteinId_column Default is "ProteinName"
#' @param factor_column Default is "Condition"
#' @return The LFQData object

create_lfqdata <- function(data, contaminant_prefix = 'CONTAMINANT_', response_level = "peptide", proteinId_column = "ProteinName",factor_column = "Condition") {

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
    atable$factorDepth <- 1
    config <- prolfqua::AnalysisConfiguration$new(atable)
    
    adata <- prolfqua::setup_analysis(data, config)
    lfqdata <- prolfqua::LFQData$new(adata, config)
    lfqdata$data <- lfqdata$data %>% filter(!grepl("^DECOY_", protein_Id)) %>% filter(!grepl(contaminant_prefix, protein_Id))
    lfqdata$remove_small_intensities()
    
    return(lfqdata)
  
}



# Function for protein aggregation of peptide normalized data
# - Use Tukey's median polish protein aggregation
# - Remove artifacts with same value in all samples 
# - Replace the artifact intensities with top3 method
# - Merge the two methods
#' @param lfqdataPeptideNorm The prolfqua object after peptide normalization
#' @return The a dataframe with protein intensities

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


# Function to keep leading protein of proteinGroups
#' @param mp_df The dataframe with the proteinGroups in protein_Id column
#' @param proteomicslfq_ids_to_keep The vector of the proteomicsLFQ leading proteins
#' @param protein_separator The separator used for multiple protein IDs, default ";"
#' @return The dataframe with only leading protein in protein_Id column
filter_for_leading_protein <- function(mp_df, proteomicslfq_ids_to_keep, protein_separator = ";") {
    proteins_in_groups <- mp_df %>%
        dplyr::filter(grepl(protein_separator, protein_Id)) %>%
        tidyr::separate_rows(protein_Id, sep = protein_separator) %>%
        dplyr::pull(protein_Id)
        
    if (length(proteins_in_groups) > 0) {
        mp_df_multiple_ids <- mp_df %>%
        dplyr::filter(grepl(protein_separator, protein_Id))
    
        mp_df_multiple_ids$protein_Id <- mapply(get_valid_ids, mp_df_multiple_ids$protein_Id, MoreArgs = list(proteomicslfq_ids_to_keep))
    
        # Merge the filtered rows with single IDs and those with multiple IDs that were processed
        mp_df <- dplyr::bind_rows(
            mp_df %>% dplyr::filter(!grepl(';', protein_Id)), 
            mp_df_multiple_ids)
        }
  
    return(mp_df)
}


# Function to filter dataframe for proteins quantified in at least n replicates of one condition
#' @param mp_df The dataframe with protein quantifications
#' @param conditions The vector of the conditions to search the mp_df colnames
#' @param min_replicates The number of min replicates needed
#' @return The dataframe with proteins passing the filtering
filter_proteins_by_replicates <- function(mp_df, conditions, min_replicates) {
    quantified_proteins <- c()
  
    for (condition in conditions) {
        condition_cols <- grep(condition, colnames(mp_df_quant), value = TRUE)
        max_nan <- length(condition_cols) - min_replicates
        mp_df[[paste0("count_na_", condition)]] <- rowSums(is.na(mp_df[, condition_cols]))

        passing_proteins <- rownames(mp_df[mp_df[[paste0("count_na_", condition)]] <= max_nan, ])
        quantified_proteins <- unique(c(quantified_proteins, passing_proteins))
        }
    
    mp_df_quant <- mp_df[rownames(mp_df) %in% quantified_proteins, ]
    mp_df_quant <- mp_df_quant %>% select(-contains('count_na_'))
    
    return(mp_df_quant)
}
