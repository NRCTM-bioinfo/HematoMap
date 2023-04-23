#' @include zzz.R
#' @importFrom methods setClass
#' @importClassesFrom Matrix dgCMatrix
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Classes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title HematoSubCluster-class
#' 
#' @description A class for sub-cluster to store information, including expression matrix, 
#'   distance matrix, meta-info and other.
#' 
#' @slot data matrix.\cr Sub-cluster expression.
#' @slot dist list.\cr Distance informations of sub-cluster.
#' @slot meta.data data.frame.\cr Meta data information of sub-cluster.
#' @slot misc list.\cr Utility slot for storing additional data.
#' 
HematoSubCluster <- setClass(
  Class = 'HematoSubCluster',
  slots = c(
    data = 'matrix',
    dist  = 'list',
    meta.data = 'data.frame',
    misc = 'list'
  )
)



#' @title CosineDist-class
#' 
#' @description A class for storing Cosine-Distance information.
#' 
#' @slot sim matrix.\cr Cosine similarity matrix.
#' @slot dist matrix.\cr Cosine distance matrix.
#' @slot features.enroll character.\cr Common features between input and reference data.
#' @slot subc.ref character.\cr The name of reference sub-cluster. It can distinguish 
#'   between reference sub-cluster and candidate sub-cluster in cosine similarity matrix and 
#'   cosine distance matrix.
#' 
CosineDist <- setClass(
  Class = 'CosineDist',
  slots = c(
    sim = "matrix",
    dist = "matrix", 
    features.enroll = "character",
    subc.ref = "character"
  )
)

