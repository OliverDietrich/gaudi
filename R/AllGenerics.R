########################################
# Accessors for genomeCollection

#' @export
setGeneric("index", function(x) standardGeneric("index"))

#' @export
setGeneric("path", function(x) standardGeneric("path"))

#' @export
setGeneric("metadata", function(x) standardGeneric("metadata"))

#' @export 
setGeneric("metadata<-", function(x, value) standardGeneric("metadata<-"))

#' @export
setGeneric("Reads", function(x, name) standardGeneric("Reads"))

#' @export 
setGeneric("Reads<-", function(x, name, value) standardGeneric("Reads<-"))

#' @export 
setGeneric("ReadsNames", function(x) standardGeneric("ReadsNames"))

########################################
# Getter/setters for Reads.

           
########################################
# Accessors for interactionSet