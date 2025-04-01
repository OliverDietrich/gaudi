#' Jellyfish k-mer counting
#'
#' CLI wrapper for Jellyfish to count k-mers
#'
#' @param input.sequence File path
#' @param output.tsv File path
#' @param k Integer, length of k-mer
#' @param hash.size Integer, Initial hash size
#' @param threads Integer, Number of threads
#' 
#' @export
jellyfish_kmer_count <- function(input.sequence=NULL,
                                 out.dir=NULL,
                                 k = 32,
                                 hash.size = 64,
                                 threads = n_proc(),
                                 remove.intermediates = TRUE
                                ) {

    stopifnot(
        !is.null(input.sequence),
        !is.null(output.txt)
    )
    check_installed('jellyfish', silent=TRUE)
    check_version('jellyfish')

    # Check input

    # Check output

    # Run jellyfish

    # Exit 0
    return(counts)
}

#' Count nucleotide k-mers using jellyfish
#' 
#' @param seq Nucleotide sequence(s) to count k-mers. Input must be character vector or list.
#' @param k Size of k-mer
#' @param threads Number of cores to use (default: all)
#'
#' @export 
jellyfish_count_kmers <- function(
    infile=NULL,
    sequence=NULL, 
    outfile=NULL,
    k=31, 
    threads=NULL,
    return.counts=TRUE
) {

    # Minimal checks
    stopifnot(
        !is.null(sequence) | !is.null(infile)
    )

    ## Check input
    if (is.null(outfile)) {
        outfile <- tempfile(pattern = 'count.tsv')        
    }

    if (file.exists(infile)) {
        fa_file <- sequence
    } else {
        fa_file <- tempfile(pattern = 'nucleotide.fa')
        sequence <- rev(unlist(lapply(sequence, c, '> Any nucleotide sequence')))
        writeLines(sequence, fa_file)
    }

    if (is.null(threads)) {
        threads <- system('nproc', intern = TRUE)
    }

    # Variables
    jf_file <- tempfile(pattern = 'count.jf')

    # Nucleotide alphabet
    # ...

    # Count k-mers
    cmd <- paste0("jellyfish count -m ",k," -s 100M -t ",threads," ", infile, " -o ", jf_file, ' 2>&1')
    system(cmd, intern=TRUE)

    # Convert jellyfish to tsv
    cmd <- paste0("jellyfish dump -c ",jf_file, ' 2>&1') # ," > ",outfile
    system(cmd, intern=TRUE)

    # Read counts
    counts <- data.table::fread(jf_tsv) # , col.names=c('kmer','freq')

    # Exit
    if (return.counts) {
        return(counts)
    } else {
        msg <- paste0('Writing output file:', outfile)
        message(msg)
        
        return(NULL)
    }
}

#' Count k-mers using mercat2
#' 
#' @param sequence
#' 
mercat_count_kmers <- function(
    sequence = NULL
) {
    
}