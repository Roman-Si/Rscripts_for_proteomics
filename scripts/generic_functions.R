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