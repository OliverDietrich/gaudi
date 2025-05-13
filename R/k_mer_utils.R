#' Jellyfish k-mer counting
#'
#' CLI wrapper for Jellyfish to count k-mers
#'
#' @param input.sequence File path
#' @param output.tsv File path
#' @param k Integer, length of k-mer
#' @param hash.size Character, Initial hash size (e.g. '100M')
#' @param threads Integer, Number of threads
#' 
#' @export
jellyfish_count_kmers <- function(input.sequence=NULL,
                                  output.tsv=NULL,
                                  temp.dir='tmp/jellyfish_count/',
                                  k = 32,
                                  hash.size = '100M',
                                  threads = n_proc(),
                                  remove.intermediates = TRUE
                                 ) {

    # Minimal check
    stopifnot(
        !is.null(input.sequence)
    )
    check_installed('jellyfish', silent=TRUE)
    check_version('jellyfish')

    # Set variables
    input_file <- paste0(temp.dir,'input_sequence.txt')
    jf_file <- paste0(temp.dir,'mer_count.jf')
    if (!is.null(output.tsv)) {
        output_file <- output.tsv
        if (file.exists(output.tsv)) {
            msg <- paste('Output file',output.tsv,'already exists.')
            warning(msg)
            counts <- vroom::vroom(output.tsv, show_col_types=FALSE, col_names = c('kmer','count'))
            counts <- setNames(counts$count, counts$kmer)
            return(counts)
        }
    } else {
        output_file <- paste0(temp.dir,'mer_count.tsv')
    }

    # Check input
    if (is_file(input.sequence)) {
        input_file <- input.sequence
    } else {
        writeLines(input.sequence, input_file)
    }
    if (!is_file(input_file)) {
        msg <- paste('Input file',input_file,'does not exist')
        warning(msg)
        return(NULL)
    } else if (dir.exists(temp.dir)) {
        msg <- paste('Temporary directory',temp.dir,'already exists. Consider cleaning up your working directory or specify a custom directory.')
        warning(msg)
    } else {
        dir.create(temp.dir, recursive=TRUE)
    }

    # Check output

    # Count k-mers
    cmd <- paste('jellyfish','count','-m',k,'-s',hash.size,'-t',threads,' ',input_file,'-o',jf_file,'2>&1')
    system(cmd, intern=TRUE)

    # Dump counts to file
    cmd <- paste('jellyfish','dump','-c',jf_file,'>',output_file,'2>&1')
    stdout <- system(cmd, intern=TRUE)

    # Read file
    counts <- vroom::vroom(output_file, show_col_types=FALSE, col_names = c('kmer','count'))
    counts <- setNames(counts$count, counts$kmer)

    # Remove intermediates
    if (remove.intermediates) {
        unlink(temp.dir, recursive=TRUE)
    }

    # Exit 0
    return(counts)
}

#' Count k-mers using mercat2
#' 
#' @param sequence
#' 
mercat_count_kmers <- function(
    sequence = NULL
) {
    
}