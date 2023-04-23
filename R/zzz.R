#' @importFrom methods new setClass setClassUnion setMethod setOldClass
#' @importClassesFrom Matrix dgCMatrix
#' @noRd
NULL


# Class definitions
setClassUnion(name = 'OptionalCharacter', members = c('NULL', 'character'))
setClassUnion(name = 'OptionalList', members = c('NULL', 'list'))

setOldClass(Classes = 'package_version')



# Hooks
.onAttach <- function(libname, pkgname) {
  return(invisible(x = NULL))
}

.onLoad <- function(libname, pkgname) {
  return(invisible(x = NULL))
}




