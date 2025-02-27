import_salmon_data <- function(star_salmon_path) {
  # Read tx2gene file
  require(stringr)
  require(tximport)
  tx2gene <- read.delim(paste0(star_salmon_path, "/tx2gene.tsv"), header = FALSE)
  
  # Get list of quant.sf files
  files <- list.files(path = star_salmon_path,
                      pattern = "quant.sf",
                      recursive = TRUE,
                      full.names = TRUE)
  
  # Create names for the files
  toremove <- c(star_salmon_path, "/", "quant.sf")
  names(files) <- str_remove_all(files, paste(toremove, collapse = "|"))
  
  # Import quantification data
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  
  return(txi)
}

