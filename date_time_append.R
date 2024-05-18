date.time.append <- function(str, sep = '_', date.format ="%y%m%d%H%M%S") {
  stopifnot(is.character(str))
  return(paste(str, format(Sys.time(), date.format), sep = sep))  
}