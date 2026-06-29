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

#' Compute Mash distance for genomic sequences
#'
#' Estimate sequence identity (average nucleotide identity, ANI) and containment
#' (alignment fraction, AF) for biological sequences. Input can be either Biostrings
#' objects (DNAStringSet, AAStringSet) or character vector of FASTA file paths. If
#' no reference sequences are supplied, query sequences will be compared to each other.
#' Sequence names will be derived from ...
#'
#' @param query XStringSet or vector of file paths.
#' @param reference XStringSet or vector of file paths.
#' @param output.dir Path to output directory
#' @param k Integer, k-mer size
#' @param s Integer, sketch size
#' @param individual.seqs Boolean, whether to sketch individual sequences rather than whole files
#' @param alphabet Character, use 'nucleotide' (ACGT) or 'aminoacid' (A-Z, except BJOUXZ) alphabet
#' @param threads Integer, this many threads will be spawned for processing
#' @param recompute.sketch Boolean, whether to re-compute the sketch file
#' @param recompute.dist Boolean, whether to re-compute the Mash distances
#'
#' @importFrom Biostrings writeXStringSet
#' @importFrom readr read_tsv
#'
#' @export
#'
mash_distance <- function(query, reference=NULL, out.dir, 
                          k=21, 
                          s=1e4,
                          alphabet = 'nucleotide',
                          threads=n_proc(),
                          recompute = FALSE
                         ) {

    # Minimal check
    check_installed('mash', silent=TRUE)
    check_version('mash')

    # Variables
    cols <- list(
        'query' = col_character(),
        'reference' = col_character(),
        'distance' = col_number(),
        'p-value' = col_number(),
        'shared-hashes' = col_character()
    )
    out.dir <- if (endsWith(out.dir,'/')) out.dir else paste0(out.dir,'/');
    if (!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE);
    distances.tsv <- paste0(out.dir,'mash_dist','_k',k,'_s',s,'.tsv')

    # Timestamp
    t0 <- Sys.time()

    # Output
    if (recompute) unlink(distances.tsv)
    existing <- if (file.exists(distances.tsv)) {
        readr::read_tsv(distances.tsv, col_names = names(cols), col_types = cols)
    } else {
        data.frame('query'=character(), 'reference'=character())
    }
    
    # Input
    symmetric <- if (is.null(reference)) TRUE else FALSE;
    reference <- if (is.null(reference)) query else reference;
    if (any(duplicated(query))) stop('Duplicated entries in Query. Check your input!');
    if (any(duplicated(reference))) stop('Duplicated entries in Reference. Check your input!');
    
    # Handling missing files
    qind <- file.exists(query)
    rind <- file.exists(reference)
    if (any(!qind) | any(!rind)) warning('Some files do not exists and will be removed...')
    query <- query[qind]
    reference <- reference[rind]

    # Decompose matrix
    ind <- query %in% existing$query
    query.old <- query[ind]
    query.new <- query[!ind]
    ind <- reference %in% existing$reference
    ref.old <- reference[ind]
    ref.new <- reference[!ind]

    X <- length(query) * length(reference)
    A <- length(query.old) * length(ref.old)
    B <- length(query.old) * length(ref.new)
    C <- length(query.new) * length(ref.old)
    D <- length(query.new) * length(ref.new)
    msg <- paste0('Found ',A,' distances out of ',X,'. Computing ', X-A)
    message(msg)
    cat('','[',A,'|',B,']','\n','[',C,'|',D,']','\n')

    # Exit 1
    if (A & !B & !C & !D) {

        # Format
        result <- format_mash_distances(existing, query, reference)

        # Timestamp
        t1 <- Sys.time()
        print(t1 - t0)
        
        # Exit
        return(result)
    }

    # Compute sketches
    if (A | B) {
        Q.old <- mash_sketch(query.old, 
                             out.dir = out.dir, 
                             k = k, 
                             s = s, 
                             alphabet = alphabet,
                             prefix = 'oldquery', 
                             recompute.sketch = recompute, 
                             verbose=FALSE
                            )
    }
    if (A | C) {
        R.old <- if (symmetric) Q.old else mash_sketch(ref.old, 
                                                       out.dir = out.dir, 
                                                       k = k, 
                                                       s = s, 
                                                       alphabet = alphabet,
                                                       prefix = 'oldref', 
                                                       recompute.sketch = recompute, 
                                                       verbose=FALSE
                                                      );
    }
    if (C | D) {
        Q.new <- mash_sketch(query.new, 
                             out.dir = out.dir, 
                             k = k, 
                             s = s, 
                             alphabet = alphabet,
                             prefix = 'query', 
                             recompute.sketch = recompute, 
                             verbose = FALSE
                            )
    }
    if (B | D) {
        R.new <- if (symmetric) Q.new else mash_sketch(ref.new, 
                                                       out.dir = out.dir, 
                                                       k = k, 
                                                       s = s, 
                                                       alphabet = alphabet,
                                                       prefix = 'reference', 
                                                       recompute.sketch = recompute, 
                                                       verbose = FALSE);
    }
    
    # # Compute distances
    if (A) {
        cmd.A <- paste('mash dist','-p',threads,Q.old,R.old,'>>',distances.tsv)
        system3(cmd.A, include.errors = FALSE)
    }
    if (B) {
        cmd.B <- paste('mash dist','-p',threads,Q.old,R.new,'>>',distances.tsv)
        system3(cmd.B, include.errors = FALSE)
    }
    if (C) {
        cmd.C <- paste('mash dist','-p',threads,Q.new,R.old,'>>',distances.tsv)
        system3(cmd.C, include.errors = FALSE)
    }
    if (D) {
        cmd.D <- paste('mash dist','-p',threads,Q.new,R.new,'>>',distances.tsv)
        system3(cmd.D, include.errors = FALSE)
    }
    
    # Read output & format
    result <- readr::read_tsv(distances.tsv, col_names = names(cols), col_types = cols)
    result <- format_mash_distances(result, query, reference)

    # Timestamp
    t1 <- Sys.time()
    print(t1 - t0)

    # Exit 0
    return(result)
}

#' Format Mash distances
#'
format_mash_distances <- function(distances, query, reference) {
    
    # Minimal check
    stopifnot(
        is.data.frame(distances),
        is.character(query),
        is.character(reference)
    )

    # Variables
    cols <- list(
        'query' = col_character(),
        'reference' = col_character(),
        'distance' = col_number(),
        'p-value' = col_number(),
        'shared-hashes' = col_character()
    )

    # Input
    check <- names(distances) == names(cols)
    if (!all(check)) stop('Formatting error in mash distances: Mismatch in column names.')

    # Compute ANI
    distances$ANI <- 1 - distances$distance

    # Lookup names
    if (!is.null(names(query))) {
        lookup <- setNames(names(query), query)
        distances$query <- lookup[distances[['query']]]
        if (is.null(reference)) {
            distances$reference <- lookup[distances[['reference']]]
        }
    }
    if (!is.null(names(reference))) {
        lookup <- setNames(names(reference), reference)
        distances$reference <- lookup[distances[['reference']]]
    }

    # Compute jaccard sim/dist from matches
    matches <- str_split(distances[['shared-hashes']], '\\/', simplify=TRUE)    
    distances$jaccard_similarity <- as.numeric(matches[,1]) / as.numeric(matches[,2])
    distances$jaccard_similarity <- distances$jaccard_similarity
    distances$jaccard_distance <- 1 - distances$jaccard_similarity

    # Exit
    return(distances)
}

#' Mash screen
#'
#' Screen the containment of a query in a set of reference sequences.
#'
#' @importFrom readr read_tsv
#' 
mash_screen <- function(query, reference, out.dir, 
                        k = 21, s = 1e3, 
                        threads = n_proc(), 
                        recompute = FALSE
                       ) {

    # Minimal check
    stopifnot(
        is_file(query),
        is.character(reference)
    )

    # Variables
    cols <- list(
        'identity' = col_number(),
        'shared-hashes' = col_character(),
        'median-multiplicity' = col_number(),
        'p-value' = col_number(),
        'query-ID' = col_character(),
        'query-comment' = col_character()
    )
    out.dir <- if (endsWith(out.dir, '/')) out.dir else paste0(out.dir, '/')
    dir.create(out.dir, recursive = TRUE)
    distances.tsv <- paste0(out.dir,'mash_screen','_k',k,'_s',s,'.tsv')

    # Input
    query.sketch <- mash_sketch(query, 
                                out.dir = out.dir, 
                                k = k, 
                                s = s, 
                                prefix = 'query', 
                                recompute.sketch = recompute
                               )
    reference <- if (length(reference) > 1) paste(reference, collapse = ' ') else reference

    # Output
    if (recompute) unlink(distances.tsv)

    # Main
    cmd <- paste('mash screen','-p',threads,'-i -1',query.sketch,reference,'>',distances.tsv)
    system3(cmd, include.errors = FALSE)

    # Read
    result <- readr::read_tsv(distances.tsv, col_names = names(cols), col_types = cols)

    # Exit
    return(result)
}

#' Create mash sketch from sequences
#'
#' @param input XStringSet or vector of FASTA file paths
#' @param prefix Name added as file prefix (e.g. query, reference)
#' @param output.dir Path to output directory
#' @param k Integer, k-mer size
#' @param s Integer, sketch size
#'
#' @importFrom readr read_tsv
#' 
#' @export
#'
mash_sketch <- function(input, out.dir, prefix, 
                        k=21, 
                        s=1e3,
                        alphabet = 'nucleotide',
                        threads=n_proc(),
                        recompute.sketch = FALSE,
                        verbose = TRUE
                       ) {

    # Minimal check
    check_installed('mash', silent=TRUE)
    check_version('mash')
    stopifnot(
        is.character(input),
        length(input) > 0
    )

    # Variables
    out.dir <- if (endsWith(out.dir,'/')) out.dir else paste0(out.dir,'/')
    input.fasta <- paste0(out.dir,prefix,'.fasta')
    input.list <- paste0(out.dir,prefix,'_sequence-list.txt')
    input.names <- paste0(out.dir,prefix,'_sequence-names.txt')
    output.sketch <- paste0(out.dir,prefix,'_k',k,'_s',s,'.msh')
    sketch.info <- paste0(out.dir,prefix,'_k',k,'_s',s,'.info')
    aa_flag <- if (alphabet == 'aminoacid') '-a' else NULL;

    # Input
    ind <- duplicated(input)
    if (sum(ind)) stop('Duplicated sequence file(s) in query. Aborting...')

    # Set names
    names(input) <- if (is.null(names(input))) input else names(input)

    # Check file
    ind <- file.exists(input)
    if (!all(ind)) {
        missing <- paste(input[which(!ind)], collapse=', ')
        input <- input[which(ind)]
        msg <- paste('Some input file(s) do not exist and have been removed:', missing)
        warning(msg)
    }

    # Output
    if (recompute.sketch) unlink(output.sketch)
    if (file.exists(output.sketch)) {
        
        # Fetch sketch info
        cmd <- paste('mash info','-t', output.sketch,'>',sketch.info)
        system3(cmd)
        # Read
        info <- readr::read_tsv(sketch.info, 
                         col_names = c('Hashes','Length','ID','Comment'), 
                         col_types = list(col_number(), col_number(), col_character(), col_character()), 
                         comment = '#'
                        )
        # Check
        recompute.sketch <- if (!all(info$Hashes == s)) TRUE else recompute.sketch
        recompute.sketch <- if (!all(info$ID %in% names(input))) TRUE else recompute.sketch
    } else {
        recompute.sketch <- TRUE
    }

    # Create sequence list
    if (!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
    writeLines(input, input.list)

    # Main
    cmd <- paste('mash sketch',aa_flag,'-k',k,'-s',s,'-p',threads,'-o',output.sketch,'-l',input.list)
    if (recompute.sketch) system3(cmd)

    # View sketch header
    cmd <- paste('mash info','-H', output.sketch)
    if (verbose) system3(cmd)

    # Read sketch info
    cmd <- paste('mash info','-t', output.sketch,'>',sketch.info)
    system3(cmd)
    info <- readr::read_tsv(sketch.info, 
                            col_names = c('Hashes','Length','ID','Comment'), 
                            col_types = list(col_number(), col_number(), col_character(), col_character()), 
                            comment = '#'
                           )
    
    # Exit 0
    return(output.sketch)
}

#' Mash info
#'
#' Display information about sketch files
#'
mash_info <- function(x) {

    # Variables
    cols <- list(
        'Hashes' = col_number(), 
        'Length' = col_number(), 
        'ID' = col_character(), 
        'Comment' = col_character()
    )
    info.tsv <- str_replace(x, '.msh', '.info')

    # Input
    if (!endsWith(x, '.msh')) stop('Sketch file (.msh) required for mash info. Aborting...')

    # Main
    cmd <- paste('mash info','-t', x, '>', info.tsv)
    stdout <- system3(cmd)
    df <- readr::read_tsv(info.tsv, 
                   col_names = names(cols), 
                   col_types = cols, 
                   comment = '#'
                  )
    
    # Exit
    return(df)
}

#' Compute Mash distance between pre-computed sketches
#' 
#' @param query Character, path to query sketch
#' @param reference Character, path to reference sketch
#' @param out.file Character, file name to output.tsv
#' @param symmetric Whether to also compare reference to query, in contrast to only query to ref
#'
#' @importFrom readr read_tsv col_character col_number
#'
#' @export
#'
mash_sketch_distance <- function(query, reference=NULL, out.file,
                                 symmetric = FALSE,
                                 threads = n_proc()
                                ) {

    # Variables
    reference <- if (is.null(reference)) query else reference
    output_cols <- list(
        'reference' = col_character(),
        'query' = col_character(),
        'distance' = col_number(),
        'p_value' = col_number(),
        'matches' = col_character()
    )

    # Check input
    if (!endsWith(query, '.msh')) stop('Query must be a sketch file (.msh). Aborting...')
    if (!endsWith(reference, '.msh')) stop('Reference must be a sketch file (.msh). Aborting...')
    query.sketch <- mash_info(query)
    reference.sketch <- mash_info(reference)

    # Check output
    if (endsWith(query, '.tsv')) stop('Output must be a TSV file. Aborting...')
    if (file.exists(out.file)) {
        msg <- paste('File',out.file,'already exists. Reading...')
        warning(msg)
        result <- readr::read_tsv(out.file, col_names = names(output_cols), col_types = output_cols)
        check <- all(result$query %in% query.sketch$ID) & all(result$reference %in% reference.sketch$ID)
        if (check) return(result)
    }

    # Run Mash
    cmd <- paste('mash dist','-p',threads,reference,query,'>',out.file)
    system3(cmd)

    # Run Mash
    cmd <- paste('mash dist','-p',threads,query,reference,'>>',out.file)
    if (symmetric) system3(cmd)

    # Read output
    result <- readr::read_tsv(out.file, col_names = names(output_cols), col_types = output_cols)

    # Exit
    return(result)
}

#' Skani sketch
#' 
#' Sketch (index) genomes for fast, robust ANI calculation
#'
#' @param input Path(s) to FASTA files
#' @param out.dir Path to output directory
#' @param individual.seqs Whether to create sketches for individual sequences instead of entire files
#' @param threads Number of threads
#'
#' @export
#'
skani_sketch <- function(input, out.dir, 
                         individual.seqs = FALSE,
                         threads = n_proc(),
                         recompute = FALSE
                        ) {

    # Minimal check
    check_installed('skani', silent=TRUE)
    check_version('skani')

    # Variables
    out.dir <- if (endsWith(out.dir,'/')) out.dir else paste0(out.dir,'/')
    i_flag <- if (individual.seqs) '-i' else NULL
    input.lines <- 'sequence.list'
    out.sketch <- paste0(out.dir,'sketches.db')

    # Check output
    if (recompute) unlink(out.dir, recursive=TRUE)
    if (file.exists(out.sketch)) {
        msg <- paste('Output file',out.sketch,'already exists.')
        message(msg)
        return()
    }
    if (dir.exists(out.dir)) unlink(out.dir, recursive=TRUE)

    # Check input
    if (!is.character(input)) stop('Input must be a vector of FASTA files.')
    ind <- file.exists(input)
    input <- if (sum(ind)) input[which(ind)] else input
    if (sum(ind)) warning('Some input files were not present and have been removed.')
    writeLines(input, input.lines)

    # Main
    cmd <- paste('skani sketch','-t',threads,i_flag,'-o',out.dir,'-l',input.lines)
    system3(cmd)

    # Exit 0
    return()
}