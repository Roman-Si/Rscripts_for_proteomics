# Required Libraries
library(tidyverse)
library(ggplot2)
library(MSnbase)
source("scripts/generic_functions.R")

# Function to extract protein coverage, peptide and PSM count from mzTab
#' @param mzt Path to mzTab file
#' @param contaminant_prefix The prefix for contaminants, default is CONTAMINANT_
#' @param output_file Name of the output file (optional)
#' @return Data frame with filtered protein information

extract_protein_stats <- function(mzt, contaminant_prefix = 'CONTAMINANT_', output_file = NULL) {
  prot_info <- data.frame(proteins(MzTab(mzt))) %>% 
    filter(!grepl(contaminant_prefix, accession), opt_global_result_type == "protein_details") %>% 
    select(accession, protein_coverage, opt_global_nr_found_peptides)
 
  psms <- data.frame(psms(MzTab(mzt))) %>% 
    filter(!grepl(contaminant_prefix, accession)) %>% 
    mutate(accession = strsplit(as.character(accession), ",")) %>% # in case of shared PSM proteins are comma separated
    unnest(accession) %>%
    distinct() %>%
    group_by(accession) %>%
    summarize(PSMs = n_distinct(PSM_ID))

  # merge the two dataframes
  prot_info <- merge(prot_info, psms, by = "accession", all = TRUE)

  if (!is.null(output_file)) {
    save_gzipped_csv(prot_info, output_file)
  }

  return(prot_info)
  
}


# Function that returns a vector with leading proteinId per proteinGroup based on proteomicsLFQ mzTab
#' @param mzt Path to mzTab file
#' @param contaminant_prefix The prefix for contaminants, default is CONTAMINANT_
#' @return Vector with one protein per proteinGroug

get_proteinIds_from_proteomicslfq <- function(mzt, contaminant_prefix = 'CONTAMINANT_') {
  proteomicslfq <- data.frame(proteins(MzTab(mzt))) 
  proteomicslfq <- proteomicslfq %>% filter(!grepl(contaminant_prefix, accession), opt_global_result_type != "protein_details")
  proteomicslfq_ids_to_keep <- proteomicslfq$accession
  return(proteomicslfq_ids_to_keep)
}



# Function to plot PTM counts (counting each modification once per spectrum)
#' @param mzt Path to mzTab file
#' @param contaminant_prefix The prefix for contaminants, default is CONTAMINANT_
#' @param plot_type Specifies the type of plot to generate: "separate" for plotting the frequency of each modification individually (recommended for more than 3 PTMs), or "mixed" for plotting the frequency of all possible combinations of modifications.
#' @param output_jpg Name of the output jpg image (optional)
#' @param jpg_width Image width, default 10
#' @param jpg_height Image heigth, default 8
#' @return The plot with PTM frequencies

count_peptide_modifications <- function(mzt, contaminant_prefix = 'CONTAMINANT_', plot_type = "mixed",  output_jpg = NULL, jpg_width = 10, jpg_height = 8) {

  if (!plot_type %in% c("mixed", "separate")) {
    stop("Invalid plot_type. Choose either 'mixed' or 'separate'.")
  }

  # Mapping UniMod numbers to names
  replacement_mapping <- c(
    "UNIMOD:4" = "Carbamidomethyl",
    "UNIMOD:35" = "Oxidation",
    "UNIMOD:1" = "Acetyl",
    "UNIMOD:7" = "Deamidation",
    "UNIMOD:28" = "Pyroglu",
    "UNIMOD:385" = "Ammonia-loss")

  psms <- data.frame(psms(MzTab(mzt))) %>% 
    filter(!grepl(contaminant_prefix, accession))

  psms$modifications <- gsub("[0-9]+-UNIMOD", "UNIMOD", psms$modifications)
  # Replace Unimod accessions with names
  psms$modifications <- sapply(psms$modifications, function(mod_string) {
    for (pattern in names(replacement_mapping)) {
      mod_string <- gsub(pattern, replacement_mapping[pattern], mod_string)
    }
    return(mod_string)
  })
  # Sort modifications for each PSM
  psms$modifications <- sapply(psms$modifications, function(mod_string) {
    mods <- unlist(strsplit(mod_string, ","))
    mods <- unique(mods)  # Keep each modification once per PSM
    mods <- sort(mods)    # Sort them alphabetically
    return(paste(mods, collapse = ","))
  })
  
  term_counts <- table(psms$modifications)
  term_counts_df <- data.frame(term = names(term_counts), frequency = term_counts)
  
  # Plot frequency of all PTM combinations
  plot_mixed <- ggplot(term_counts_df %>% filter(term != ""), aes(x = term, y = frequency.Freq)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = frequency.Freq), vjust = -0.5, size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.2)) +
    labs(title = "Frequency of PTMs", subtitle = paste("Total number of PSMs is:", nrow(psms), "\nPlotting each PTM once per PSM")) +
    xlab("Modification") +
    ylab("Frequency")
  

  # Split terms and sum counts for each separate modification, suggested plot in case of many PTMs
  separate_modifications <- unique(unlist(strsplit(term_counts_df$term, ",")))
  modification_sums <- sapply(separate_modifications, function(mod) {
    sum(term_counts_df$frequency.Freq[grep(mod, term_counts_df$term)])
  })
  modification_sums_df <- data.frame(modification = names(modification_sums), total_count = modification_sums)
  
  # Plot total counts for each separate modification
  plot_separate <- ggplot(modification_sums_df, aes(x = modification, y = total_count)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = total_count), vjust = -0.5, size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.2)) +
    labs(title = "Frequency of PTMs", subtitle = paste("Total number of PSMs is:", nrow(psms), "\nPlotting each PTM once per PSM")) +
    xlab("Modification") +
    ylab("Frequency")
  
  if (!is.null(output_jpg)) {
    if (plot_type == "mixed") {
      ggsave(file.path(output_jpg), plot = plot_mixed, width = jpg_width, height = jpg_height)
    } else {
      ggsave(file.path(output_jpg), plot = plot_separate, width = jpg_width, height = jpg_height)
    }
  }
  
  if (plot_type == "mixed") {
    return(plot_mixed)
  } else {
    return(plot_separate)
  }
}


