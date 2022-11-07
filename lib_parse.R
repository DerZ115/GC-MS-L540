split.vec <- function(vec, sep = "") {
  is.sep <- vec == sep
  split(vec[!is.sep], cumsum(is.sep)[!is.sep])
}

parse_libsearch <- function(filename) {
  require(stringr)
  # Empty vectors for data from match
  rts <- c()
  qualities <- c()
  cas_n <- c()
  compounds <- c()
  
  search_data <- readLines(filename)
  search_data <- head(tail(search_data, -18), -3) # Remove header & date
  search_data <- trimws(search_data, whitespace="[ \t\r\n.]")
  search_data <- split.vec(search_data, sep = "") # Split at empty lines
  
  for (entry in search_data) {
    line1 <- entry[1] # Take most likely match
    rt <- as.double(unlist(strsplit(line1, "\\s+"))[2]) # Retention time
    
    # Extract compound name (can be multiple lines)
    cmp_starts <- grep(" \\d{1,6} \\d{1,7}-\\d{2}-\\d ", entry)
    cmp_starts <- c(cmp_starts, length(entry)+1)
    
    rest <- entry[cmp_starts[1]:(cmp_starts[2]-1)]
    tmp <- str_match(rest[1], "\\s+\\d{1,6} (\\d{1,7}-\\d{2}-\\d)\\s+(\\d+)")
    
    cas <- tmp[1,2]
    quality <- tmp[1,3]
    cmpd_name <- str_remove(rest, "\\s+\\d{1,6} (\\d{1,7}-\\d{2}-\\d)\\s+(\\d+)")
    
    rts <- c(rts, rt)
    qualities <- c(qualities, quality)
    cas_n <- c(cas_n, cas)
    compounds <- c(compounds, paste0(cmpd_name, collapse = ""))
    
    
  }
  
  result <- data.frame(rts, cas_n, compounds, qualities)
  colnames(result) <- c("Ret.Time", "CAS", "Compound", "Quality")
  return(result)
}