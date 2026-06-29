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
AddReads <- function(object, data, name, samples = 'index', R1 = 'R1', R2 = 'R2', S = 'S', L = 'L') {

    # Variables
    all_cols <- c(samples, R1, R2, S, L)
    cols_names <- c('index','R1','R2','S','L')
    
    # Minimal check
    stopifnot(
        is(object) == 'genomeCollection',
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

    # Remove samples that do not have ANY associated files
    # BUG: when multiple reads present (R1, R2, L) and only some files are hybrid, any NA removes the whole sample
    df_remove <- data.frame(row.names = 1:nrow(data))
    for (i in names(data)) {
        if (i == 'index') {
            df_remove[[i]] <- !is.na(data[[i]])
        } else {
            df_remove[[i]] <- file.exists(data[[i]])
        }
    }
    row_index <- apply(df_remove, 1, all)
    #data <- data[row_index, ] # BUG - currently there is NO FILTERING !!!
    
    # Format data
    vectors <- list()
    for (i in c('R1', 'R2', 'S', 'L')) {
        if (i %in% names(data)) {
            vectors[[i]] <- as.character(data[[i]])
        } else {
            vectors[[i]] <- character(length = 0)
        }
    }

    # Concatenate duplicated reads of the same type
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
    # This is super weird. Re-phrase ASAP...

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
# validity check for Reads class object
#-------------------------------------------------------------------------------

.valid.Reads <- function(object) {
    msg <- NULL

    if (length(object@index) < 1) {
        msg <- c(msg, "Too few samples registered for this Reads object.")
    }

    if (length(object@R1) != length(object@R2)) {
        msg <- c(msg, "Forward (R1) and reverse (R2) reads do not match.")
    }

    if (any(duplicated(object@alias))) {
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

    # Summarize
    n <- length(unique(object@index))
    m <- length(object@alias)
    m <- if (length(m)) m else n
    print.index <- if (n > 5) c(head(object@index,3), "...", tail(object@index,3)) else object@index
    file.summary <- list(
        'paired' = file.exists(object$R1) & file.exists(object$R2),
        'unpaired' = file.exists(object$S),
        'long' = file.exists(object$L)
    )
    file.summary <- sapply(file.summary, sum)
    file.summary <- paste(names(file.summary), file.summary, sep=': ')

    # Print
    cat(
        is(object),"\n","containing", n, "samples with", m, "FASTQ files:", print.index,
        "\n",
        "Files exist?", file.summary, "\n"
    )
}

#' @export 
setMethod("show", "Reads", .show.Reads)

#-------------------------------------------------------------------------------
# Accessors
#-------------------------------------------------------------------------------

#' Dollar-sign autocompletion
#'
#' @importFrom utils .DollarNames
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

#' Accessors for the FASTQ files of a Reads object.
#'
#' @author Oliver Dietrich
#' @export
#' 
setMethod("paired", "Reads", function(x) {
    if (!length(x@R1)) return(NULL)
    forward <- setNames(x@R1, x@index)
    reverse <- setNames(x@R2, x@index)
    ind <- file.exists(forward) & file.exists(reverse)
    short <- list(
        'R1' = forward[ind],
        'R2' = reverse[ind]
    )
    return(short)
})

setMethod("unpaired", "Reads", function(x) {
    if (!length(x@S)) return(NULL)
    short <- setNames(x@S, x@index)
    ind <- file.exists(short)
    return(short[ind])
})

setMethod("long", "Reads", function(x) {
    if (!length(x@L)) return(NULL)
    long <- setNames(x@L, x@index)
    ind <- file.exists(long)
    return(long[ind])
})

#' Accessors for the FASTQ quality summaries of a Reads object.
#'
#' @export
#' 
setMethod("Counts", "Reads", function(x) {
    if (!length(x@counts)) return(NULL)
    return(x@counts)
})

setMethod("Counts<-", "Reads", function(x, value) {
  x@counts <- value
  validObject(x)
  x
})