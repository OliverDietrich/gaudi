#' Check if program is installed on the system
#'
#' Helper function to determine if programs are installed (available) on the system.
#' Missing programs lead to an error that are intended to break higher-level functions.
#' 
#' @param program Name of the program
#' @param silent Whether to report what is found
#' 
#' @export
check_installed <- function(program=NULL, silent=FALSE) {

    # Check input
    stopifnot(
        !is.null(program),
        class(program) == 'character'
    )

    # Detect missing programs
    PATH <- Sys.which(program)
    index <- PATH == ''
    if (any(index)) {
        missing <-names(PATH[which(index)])
        msg <- paste('Programs missing. Not found:', paste(missing, collapse=', '))
        stop(msg)
    } else if (silent) {
    } else {
        msg <- paste('All programs installed. Found:',paste(names(PATH), collapse=', '))
        message(msg)
    }
}

#' Check program version
#' 
#' @param program Name of the program
#' @param command Name of the version command (default: --version)
#' @param return.version Whether to return the version number (default FALSE will only print message)
#'
#' @export
check_version <- function(program=NULL, command='--version', return.version = FALSE) {
    
    # Minimal check
    stopifnot(
        !is.null(program),
        class(program) == 'character',
        length(program) == 1
    )

    # Exceptions
    exceptions <- list(
        'dataformat' = 'dataformat version',
        'defense-finder' = 'defense-finder version'
    )

    # Version command
    cmd <- paste(program,command)
    if (program %in% names(exceptions)) {
        cmd <- exceptions[[program]]
    }
    version <- system3(cmd, include.errors = FALSE, return.stdout = TRUE)

    # Exit 2
    if (length(version) == 0) {
        stop('No output to version call.')
    }
    if (length(version) > 1) {
        version <- paste(version, collapse='\n')
    }
    
    # Message
    if (stringr::str_detect(version, program)) {
        msg <- version
    } else {
        msg <- paste(program,version)
    }
    message(msg)

    # Exit
    if (return.version) return(version)
}

#' Detect the number of processing units available
#' 
#' @export
n_proc <- function() {

    # Check program
    check_installed('nproc', silent=TRUE)

    # Call
    n <- system('nproc', intern=TRUE)

    # Format
    n <- as.numeric(n)
    
    return(n)    
}

#' Check if object refers to a single file
#'
#' @param file_name Character, file path
#' @param suffix Character, file extension (e.g. '.csv','.fasta')
#' @param silent Boolean, whether to print warnings
#'
#' @export
is_file <- function(file_name=NULL, suffix=NULL, silent=FALSE) {

    # NULL
    if (is.null(file_name)) {
        msg <- 'File name not specified. Returning FALSE!'
        if (!silent) warning(msg)
        return(FALSE)
    }

    # Ensure single file name
    if (length(file_name) < 1) {
        msg <- 'Incorrect file specification: object of length 0.'
        if (!silent) warning(msg)
        return(FALSE)
    }
    if (class(file_name) != 'character') {
        msg <- 'Incorrect file specification: object not a character vector.'
        if (!silent) warning(msg)
        return(FALSE)
    }
    if (length(file_name) > 1) {
        msg <- 'Incorrect file specification: object length > 1.'
        if (!silent) warning(msg)
        return(FALSE)
    }

    # File suffix
    if (!is.null(suffix)) {
        index <- endsWith(file_name, suffix)
        if (!any(index)) {
            msg <- paste('File',file_name,'does not carry the extension',paste(suffix, collapse=', '))
            if (!silent) warning(msg)
            return(FALSE)
        }
    }

    # Return TRUE
    if (file.exists(file_name)) return(TRUE) else return(FALSE)
}

#' Check FASTA formatting
#'
#' @param file Path to FASTA file
#' @param type Define whether 'DNA', 'RNA', or 'protein' sequences are expected.
#'
is_fasta <- function(file, type = 'any') {

    # Check input
    missing <- !file.exists(file)
    msg <- paste(sum(missing), 'out of', length(missing), 'files not found.')
    if (any(missing)) stop(msg)

    # Try reading
    FUN <- NULL
    FUN <- if (type == 'protein') Biostrings::readAAStringSet else FUN;
    FUN <- if (type == 'DNA') Biostrings::readDNAStringSet else FUN;
    FUN <- if (type == 'RNA') Biostrings::readRNAStringSet else FUN;
    n_found <- if (is.null(FUN)) length(file) else sum(unlist(lapply(lapply(file, FUN, nrec = 1), length)));
    if (n_found < length(file)) {
        warning('Not able to read all FASTA files. Aborting...')
        return(FALSE)
    }

    return(TRUE)
}

#' Check FASTQ formatting
#'
#' @param file File name
#' @param suffix List of allowed file suffixes
#'
#' @export
is_valid_fastq <- function(file=NULL, 
                           suffix=c('.fastq','.fastq.gz'),
                           multiple=FALSE
                          ) {

    stopifnot(
        !is.null(file)        
    )

    # Checks
    N <- length(file)
    if (N > 1) {
        msg <- paste(N,'files supplied.')
        if (multiple) {
            warning(msg)
        } else (
            stop(msg)
        )
    }

    ## In-memory FASTQ (class ShortReadQ)

    ## File presence
    index <- !file.exists(file)
    if (any(index)) {
        index <- which(!file.exists(file))
        msg <- paste('File(s)',paste(file[index],sep=', '),'do(es) not exist.')
        warning(msg)
        return(FALSE)
    }

    ## File name
    index <- endsWith(file, suffix)
    if (!any(index)) {
        msg <- paste('File(s) do(es) not end with an accepted suffix.')
        warning(msg)
        return(FALSE)
    }

    ## Add more sophisticated checks ...

    # Exit 0
    return(TRUE)
}