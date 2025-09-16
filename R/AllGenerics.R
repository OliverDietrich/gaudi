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

#' @export
setGeneric("Assembly", function(x, name) standardGeneric("Assembly"))

#' @export 
setGeneric("Assembly<-", function(x, name, value) standardGeneric("Assembly<-"))

#' @export 
setGeneric("Assemblies", function(x) standardGeneric("Assemblies"))

#' @export
setGeneric("Annotation", function(x, name) standardGeneric("Annotation"))

#' @export 
setGeneric("Annotation<-", function(x, name, value) standardGeneric("Annotation<-"))

#' @export 
setGeneric("Annotations", function(x) standardGeneric("Annotations"))

#' @export
setGeneric("Distance", function(x, name) standardGeneric("Distance"))

#' @export 
setGeneric("Distance<-", function(x, name, value) standardGeneric("Distance<-"))

#' @export 
setGeneric("Distances", function(x) standardGeneric("Distances"))

########################################
# Getter/setters for Reads.

#' @export 
setGeneric("paired", function(x) standardGeneric("paired"))

#' @export 
setGeneric("unpaired", function(x) standardGeneric("unpaired"))

#' @export 
setGeneric("long", function(x) standardGeneric("long"))

########################################
# Getter/setters for Annotation.

#' @export 
setGeneric("Genes", function(x) standardGeneric("Genes"))

#' @export 
setGeneric("Proteins", function(x) standardGeneric("Proteins"))

#' @export 
setGeneric("CDS", function(x) standardGeneric("CDS"))

#' @export 
setGeneric("Log", function(x) standardGeneric("Log"))

########################################
# Accessors for interactionSet