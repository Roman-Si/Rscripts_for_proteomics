# Function to save gzipped dataframe in csv excluding the rownames
#' @param df The dataframe
#' @param output_file The name of the output fle, if given in a directory that does not exist it will be created
#' @param row_names Logical; whether to include row names in the output (default is FALSE)
#' @return None
save_gzipped_csv <- function(df, output_file, row_names = FALSE) {
  output_dir <- dirname(output_file)
  
  # Create the directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(df, gzfile(output_file), row.names = row_names)
}



### Function to keep the main proteinId (as determined from proteomicsLFQ mzTab) in case of proteinGroups
get_valid_ids <- function(all_ids, valid_ids, separator = ";") {
#' @param all_ids The vector with one or multiple IDs per item
#' @param valid_ids The vector with the leading IDs to keep
#' @param separator The separator in case of multiple IDs per row (default is ;)
#' @return The filtered vector with one leading ID per group

  # Takes as input 
  # - the vector ids where multiple proteinIds are separated with ;
  # - the vector valid_ids that contains the main proteinId for each group
  
  ids_split <- strsplit(all_ids, ";")[[1]]
  valid_ids <- ids_split[ids_split %in% valid_ids]
  paste(valid_ids, collapse = ";")
}
