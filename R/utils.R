#' @noRd
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Utils
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Remove Sub-cluster Suffix
#' 
#' @param x character.\cr Sub-cluster names with number suffix
#' @param correct.cycling logical.\cr 
#' 
#' @return Return a cosine distance matrix.
#' 
remove.subc.suffix <- function(x, correct.cycling = F){
  x <- gsub("\\:[0-9]*$", "", x)
  if(correct.cycling){
    x[x %in% "Late cycling pro-B"] <- "Late pro-B"
    x[x %in% "Early cycling pro-B"] <- "Early pro-B"
    x[grepl("cycling NK", x)] <- "NK"
  }
  return(x)
}


#' Print Current Time and Message
#' 
#' @param ... ANY.\cr Message.
#' 
#' @return Print string with time.
#' 
message.time <- function (...){
  do.call(message, args = c(list(sprintf("[INFO] %s: ", Sys.time())), ...))
  invisible()
}


#' Determine whether a Data.frame has Row Names
#' 
#' @param data data.frame.\cr A data frame.
#' 
#' @return Return TRUE or FALSE for whether a data.frame has row names.
#' 
has.rownames <- function(data){
  .row_names_info(data) > 0L && !is.na(.row_names_info(data, 0L)[[1L]])
}
