#' @include AllClasses.R
#'
NULL


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Create Sub-cluster Object of HematoSubCluster
#'
#' @param object Seurat, matrix.\cr Seurat object or sub-cluster expression matrix.
#' @param ... ANY.\cr Arguments passed to other methods.
#'
#' @return Return a object of \code{\link{HematoSubCluster}}, containing sub-cluster
#'   expression matrix, and cosine similarity matrix.
#'
#' @export
#'
CreateSubclusterObject <- function(object, ...){
  UseMethod(generic = "CreateSubclusterObject", object = object)
}


#' Run MinMax-Kmean Cluster for Sub-cluster
#'
#' @param object ANY.\cr Input object.
#' @param ... ANY.\cr Arguments passed to other methods.
#'
#' @return Return a vector containing sub-cluster info.
#'
runSubcluster <- function(object, ...) UseMethod(generic = "runSubcluster", object = object)



#' Run Cosine Similarity Analysis
#'
#' @param object ANY.\cr Input object.
#' @param ... ANY.\cr Arguments passed to other methods.
#'
#' @export
#'
runCosineSimilarity <- function(object, ...) UseMethod(generic = "runCosineSimilarity", object = object)





#' @noRd
setMethod("show", signature = signature("HematoSubCluster"),
          function(object){
            cat("An object of class HematoSubCluster\n")
            if(length(object@dist)){
              cat(ncol(object@data), "sub-clusters with", names(object@dist), "distance")
            }else{
              cat(ncol(object@data), "sub-clusters without dist")
            }
          })



#' @noRd
setMethod("show", signature = signature("CosineDist"),
          function(object){
            cat("An object of class CosineDist\n")
            cat(sum(!colnames(object@sim) %in% object@subc.ref), "sub-clusters with", length(object@features.enroll), "features enroll")
          })


