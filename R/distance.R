#' Compute jaccard similarity/distance
#' 
#' Compute jaccard similarity of two character vectors. 
#'
#' @param a Character vector
#' @param b Character vector
#' @param type Character, what type should be reported (similarity, distance)
#'
#' @returns Numeric
#'
#' @export
#' 
jaccard <- function(a=NULL, b=NULL, type='similarity') {

    # Minimal check
    stopifnot(
        !is.null(a),
        !is.null(b),
        class(a) == 'character',
        class(b) == 'character',
        type %in% c('similarity','distance')
    )

    # Compute sample sets
    intersection <- length(intersect(a,b))
    union <- length(a) + length(b) - intersection

    # Similarity
    sim <- intersection / union

    # Exit 0
    if (type == 'similarity') {
        return(sim)
    } else if (type == 'distance') {
        return(1 - sim)
    }
}

#' CLI wrapper for Mash
#'
#' Compute Mash distance for genomic sequences
#'
#' @param input Character, vector of input FASTA sequences. If names are present, will be used for output IDs.
#' @param output Path, output prefix (first input file used if unspecified). The suffix '.msh' will be appended.
#' @param k Integer, k-mer size
#' @param s Integer, sketch size
#' @param sketch.archive File name, destination of the sketch archive
#' @param threads Integer, this many threads will be spawned for processing.
#' @param recompute.sketch Boolean, whether to re-compute the sketch file
#' @param recompute.dist Boolean, whether to re-compute the Mash distances
#'
#' @export
#'
mash_distance <- function(input=NULL, output.dir=NULL, 
                          k=21, s=1000,
                          threads=n_proc(),
                          recompute.sketch=FALSE,
                          recompute.dist=FALSE
                         ) {

    # Minimal check
    stopifnot(
        !is.null(input),
        !is.null(output.dir),
        file.exists(output.dir)
    )
    check_installed('mash', silent=TRUE)
    check_version('mash')

    # Variables
    if (!endsWith(output.dir,'/')) {
        output.dir <- paste0(output.dir,'/')
    }
    seq.list <- paste0(output.dir,'seq_list.txt')
    sketch.archive <- paste0(output.dir,'reference.msh')
    dist.tsv <- paste0(output.dir,'distances.tsv')
    result <- NULL
    
    # Check input
    files_present <- file.exists(input)
    if (!all(files_present)) {
        msg <- paste('Files',paste(input[!files_present],collapse=', '),'missing. Will be skipped')
        input <- input[files_present]
        warning(msg)
    }
    writeLines(input, seq.list)

    # Check output
    if (file.exists(dist.tsv) & !recompute.dist) {
        msg <- 'Output already present. Reading...'
        message(msg)
        result <- readr::read_tsv(dist.tsv, col_names = c('reference','query','distance','p_value','matches'))
    }

    # Create sketch
    if (!is.null(result)) {
        # Skip 
    } else if (file.exists(sketch.archive) & !recompute.sketch) {
        msg <- paste('Sketch archive',sketch.archive,'already exists.')
        message(msg)
    } else {
        cmd <- paste('mash sketch','-k',k,'-s',s,'-p',threads,'-o',sketch.archive,'-l',seq.list)
        cmd <- paste(cmd,'2>&1')
        stdout <- system(cmd, intern=TRUE)
        stdout <- paste(stdout, collapse='\n')
        cat(stdout)
    }
    
    # Run Mash
    if (is.null(result)) {
        cmd <- paste('mash dist','-p',threads,sketch.archive,sketch.archive,'>',dist.tsv)
        cmd <- paste(cmd,'2>&1')
        stdout <- system(cmd, intern=TRUE)
        stdout <- paste(stdout, collapse='\n')
        cat(stdout)

        # Read output
        result <- readr::read_tsv(dist.tsv, col_names = c('reference','query','distance','p_value','matches'))
    }

    # Change file path to name
    input.names <- names(input)
    if (!is.null(input.names) & !any(duplicated(input.names))) {
        lv <- setNames(names(input), input)
        result$reference <- lv[result$reference]
        result$query <- lv[result$query]        
    }

    # Exit
    return(result)
}