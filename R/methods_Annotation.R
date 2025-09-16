### Methods for the Annotation class

#-------------------------------------------------------------------------------
# validity check for Annotation class object
#-------------------------------------------------------------------------------

.valid.Annotation <- function(object) {
    msg <- NULL

    # ...

    if (length(msg)) { return(msg) }
    return(TRUE)
}

methods::setValidity("Annotation", .valid.Annotation)

#-------------------------------------------------------------------------------
# show
#-------------------------------------------------------------------------------

.show.Annotation <- function(object) {

    # Summarize
    n <- length(object$index)
    print.index <- if (n > 5) c(head(object$index,3), "...", tail(object$index,3)) else object$index
    file.summary <- list(
        'genes' = file.exists(object$genes),
        'cds' = file.exists(object$cds),
        'proteins' = file.exists(object$proteins),
        'log files' = file.exists(object$log)
    )
    file.summary <- sapply(file.summary, sum)
    file.summary <- paste(names(file.summary), file.summary, sep=': ')
    file.summary <- paste(file.summary, collapse='   ')
    
    # Print
    cat(
        is(object),"\n","containing", length(object$index), "entries:", print.index, 
        "\n",
        paste0("annotated from object[['", object$genome, "']]", " using ", object$tool, " (",object$type,")", collapse = ''),
        "\n",
        "Files exist?", file.summary, "\n"
    )
}

#' @export 
setMethod("show", "Annotation", .show.Annotation)

#-------------------------------------------------------------------------------
# Accessors
#-------------------------------------------------------------------------------

#' Dollar-sign autocompletion
#'
#' @importFrom utils .DollarNames
#'
#' @export
.DollarNames.Annotation <- function(x, pattern = "") {
    grep(pattern, slotNames(x), value=TRUE)
}

#' Data
#'
#' Get Annotation data
#'
#' @param x An Annotation object
#' @param name Name of Annotation slot
#'
setMethod("$", "Annotation", function(x, name) slot(x, name))

#' Accessors for the components of an Annotation object.
#'
#' @author Oliver Dietrich
#' @export
#'
setMethod("Genes", "Annotation", function(x) setNames(x@genes, x@index))

setMethod("Proteins", "Annotation", function(x) setNames(x@proteins, x@index))

setMethod("CDS", "Annotation", function(x) setNames(x@CDS, x@index))

setMethod("Log", "Annotation", function(x) setNames(x@log, x@index))