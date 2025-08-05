### Methods for the Reads class

#-------------------------------------------------------------------------------
# constructor function for Reads class
#-------------------------------------------------------------------------------


#' Add a \code{Reads} object to a \code{\link{genomeCollection}} object
#'
#' Create a \code{Reads} object from a data.frame with defined components
#' and add it to the reads slot of a genomeCollection object.
#' 
#' @param object A \code{\link{genomeCollection}} object
#' @param data Data.frame containing the paths to reads files
#' @param name Slot name of the Reads object in the genomeCollection
#' @param samples Name of the column containing sample names
#' @param R1 Name of the column containing file paths with forward paired-end reads
#' @param R2 Name of the column containing file paths with reverse paired-end reads
#' @param S Name of the column containing file paths with unpaired reads
#' @param L Name of the column containing file paths with long reads
#'
#' @return A \code{\link{genomeCollection}} object
#'
#' @rdname Reads
#' @export
#'
AddReads <- function(object = NULL, data=NULL, name=NULL, samples = 'index', R1 = 'R1', R2 = 'R2', S = 'S', L = 'L') {

    # Minimal check
    all_cols <- c(samples, R1, R2, S, L)
    cols_names <- c('index','R1','R2','S','L')
    stopifnot(
        class(object) == 'genomeCollection',
        !is.null(name),
        class(data) == 'data.frame',
        !any(duplicated(all_cols)),
        samples %in% names(data),
        any(c(R1, S, L) %in% names(data))
    )

    # Check output
    if (name %in% ReadsNames(object)) {
        msg <- paste('Slot',name,'already exists in ReadsNames(object). No changes made.')
        warning(msg)
        return(object)
    }

    # Subset data
    ind_cols <- all_cols %in% names(data)
    data <- subset(data, select=all_cols[ind_cols])
    names(data) <- cols_names[ind_cols]

    # Comparing IDs to index(object)
    ind.match <- data[['index']] %in% index(object)
    ind.miss <- which(!ind.match)
    if (length(ind.miss)) {
        non_matching_ids <- paste0(data[['index']][ind.miss], collapse=', ')
        msg <- paste0('Some samples (',non_matching_ids,') are not registered with the genomeCollection and will be removed.')
        warning(msg)
    }
    data <- data[ind.match, ]

    # Remove samples that do not have associated files
    df_remove <- data.frame(row.names = 1:nrow(data))
    for (i in names(data)) {
        if (i == 'index') {
            df_remove[[i]] <- !is.na(data[[i]])
        } else {
            df_remove[[i]] <- file.exists(data[[i]])
        }
    }
    row_index <- apply(df_remove, 1, all)
    data <- data[row_index, ]
    
    # Format data
    vectors <- list()
    for (i in c('R1', 'R2', 'S', 'L')) {
        if (i %in% names(data)) {
            vectors[[i]] <- as.character(data[[i]])
        } else {
            vectors[[i]] <- character(length = 0)
        }
    }

    ## Concatenate duplicated reads of the same type
    bool_duplicates <- duplicated(data[['index']])
    name_duplicates <- data[['index']][bool_duplicates]
    if (length(name_duplicates)) {
        msg <- paste0('Some samples (', paste(name_duplicates, collapse=', '), ') are duplicated. Creating alias...')
        warning(msg)
        vectors[['alias']] <- make.unique(data[['index']], '-')
    }

    # Infer type
    v_len <- sapply(vectors, length)
    type <- character(length = 0)
    if (v_len[['L']]) {type <- c(type, 'long')}
    if (v_len[['R1']] & v_len[['R2']]) {type <- c(type, 'paired')}
    if (v_len[['S']]) {type <- c(type, 'unpaired')}

    # Read summary
    smry <- data.frame()
    
    # Create object
    Reads(object, name) <- methods::new("Reads", 
                                         "index" = as.character(data[['index']]),
                                        "alias" = as.character(vectors[['alias']]),
                                         "R1" = vectors[['R1']],
                                         "R2" = vectors[['R2']],
                                         "S" = vectors[['S']],        
                                         "L" = vectors[['L']],                
                                         "counts" = smry,
                                         "quality" = smry,
                                         "type" = type
                                        )
    return(object)
}

#-------------------------------------------------------------------------------
# validity check for genomeCollection class object
#-------------------------------------------------------------------------------

.valid.Reads <- function(object) {
    msg <- NULL

    if (length(object$index) < 1) {
      msg <- c(msg, "Too few samples registered for this Reads object.")
    }

    if (any(duplicated(object$alias))) {
      msg <- c(msg, "Duplicated alias detected. Something went wrong.")
    }

    if (length(msg)) { return(msg) }
    return(TRUE)
}

methods::setValidity("Reads", .valid.Reads)

#-------------------------------------------------------------------------------
# show
#-------------------------------------------------------------------------------

.show.Reads <- function(object) {

    # Print
    cat(
        str(object, max.level = 2)
    )
}

#' @export 
setMethod("show", "Reads", .show.Reads)

#-------------------------------------------------------------------------------
# Accessors
#-------------------------------------------------------------------------------

#' Dollar-sign autocompletion
#'
#' @export
.DollarNames.Reads <- function(x, pattern = "") {
    grep(pattern, slotNames(x), value=TRUE)
}

#' Data
#'
#' Get Reads data
#'
#' @param x A Reads object
#' @param name Name of Reads slot
#'
setMethod("$", "Reads", function(x, name) slot(x, name))