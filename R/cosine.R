#' @importFrom Matrix tcrossprod
#' @importFrom Matrix as.matrix
#' @importFrom Matrix rowSums
#' @importFrom Matrix Diagonal
#' @importFrom Matrix crossprod
#' @noRd
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Cosine Distance
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Compute Cosine Similarity
#' 
#' @param mat matrix.\cr A matrix combined candidate and reference sub-cluster 
#'   expression matrix.
#' 
#' @return Return a cosine similarity matrix.
#' 
computeCosineSim <- function(mat) {
  mat <- t(mat)
  return(as.matrix(tcrossprod( rowNormalize(mat) )))
}


#' Compute Cosine Distance
#' 
#' @param mat matrix.\cr Cosine similarity matrix
#' 
#' @return Return a cosine distance matrix.
#' 
computeCosineDist <- function(mat) {
  cosine.dist <- as.matrix(as.dist(1 - mat))
  # cosine.dist <- Matrix::Matrix(as.dist(1 - mat))
  return(cosine.dist)
}



#' Normalize Matrix by Row
#' 
#' @param m matrix.\cr A matrix for row normalization.
#' 
#' @return Return a normalized matrix.
#' 
rowNormalize <- function(m) {
  d <- Diagonal(x = 1 / sqrt(rowSums(m^2)))
  return(t(crossprod(m, d)))
}
