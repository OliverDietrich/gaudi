#' Compute read quality 
summarize_read_quality <- function(fastq=NULL, read=NULL, min.base.quality=15, per.base.quality=FALSE) {

    # Check
    stopifnot(
        !is.null(fastq),
        !is.null(read),
        class(fastq) == 'ShortReadQ',
        length(read) == 1
    )

    # Extract quality scores
    qual <- as.character(as.matrix(fastq@quality@quality[[read]]))
    result <- data.frame(
        'quality' = encoding(quality(fastq))[qual],
        'position' = 1:length(qual),
        'read' = read
    )

    # EXIT 1
    if (per.base.quality) {
        msg <- 'Returning per-base quality scores. No summary across read will be computed.'
        warning(msg)
        return(result)
    }

    # Summarize quality
    result <- result %>% mutate(qualified = quality >= min.base.quality) %>% group_by(read)
    result <- result %>% summarize(avg_quality = mean(quality), length = max(position), percent_qualified = sum(qualified)/length*100)

    # EXIT 0
    return(result)
}

#' Compute FASTQ summary
#'
#' @param in_file File name.
#' @param min.base.quality PHRED score of individual bases to be qualified (impacts percent qualified).
#' @param max.reads Integer, number of reads to sample
#'
#' @export
summarize_fastq_quality <- function(in_file=NULL, min.base.quality=15, max.reads=5000) {

    # Check
    stopifnot(
        !is.null(in_file),
        is_valid_fastq(in_file)
    )

    # Long vs. short
    # Paired end reads
    
    # Read file
    # TODO: accept either file or in-memory FASTQ (reduce read time)
    all_reads <- ShortRead::readFastq(in_file) # bottleneck

    # Index sample reads
    n_reads <- length(all_reads)
    n_reads <- ifelse(n_reads < max.reads, n_reads, max.reads)
    some_reads <- sample(1:length(all_reads), n_reads)
    names(some_reads) <- some_reads

    # Compute read quality
    some_reads <- lapply(some_reads, summarize_read_quality, fastq=all_reads, min.base.quality=min.base.quality)
    some_reads <- dplyr::bind_rows(some_reads)

    return(some_reads)
}

#' Fastplong filtering of long reads
#'
#' CLI wrapper for fastplong to filter and trim long reads
#'
#' @param file_in Input file
#' @param file_out Output file
#' @param length_required Integer, reads shorter than length_required will be discarded, default is 15.
#' @param mean_qual Integer, if one read's mean_qual quality score < mean_qual, then this read is discarded. Default 0 means no requirement.
#' @param unqualified_percent_limit Integer, how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%.
#' @param qualified_quality_phred Integer, the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.
#' @param threads Numeric, worker thread number, default is all (see n_proc for details)
#' @param overwrite 
#'
#' @export
fastplong <- function(file_in=NULL, file_out=NULL, 
                      length_required=100L, mean_qual=20L,
                      unqualified_percent_limit=40L, qualified_quality_phred = 15L,
                      threads=n_proc(), overwrite=FALSE
                     ) {

    # Minimal check
    stopifnot(
        !is.null(file_in),
        is_valid_fastq(file_in),
        !is.null(file_out)
    )
    check_installed('fastplong', silent=TRUE)
    check_version('fastplong')

    # Check output
    if (file.exists(file_out) & overwrite) {
        msg <- paste('Output file',file_out,'exists. Will be replaced.')
        warning(msg)
    } else if (file.exists(file_out) & !overwrite) {
        msg <- paste('Output file',file_out,'exists. Aborting...')
        return(msg)
    } else if (!endsWith(file_out, '.fastq')) {
        msg <- paste0('Output file "',file_out,'" is not formatted well. Please add .fastq.')
        stop(msg)
    }

    # TODO: function writes fastplong.html and fastplong.json to working directory. Should be stored in log/ or similar...
    
    # Main
    cmd <- paste0('fastplong -i ',file_in,' -o ',file_out)
    cmd <- paste0(cmd,' --length_required ',length_required,' --mean_qual ',mean_qual,' --unqualified_percent_limit ',unqualified_percent_limit,
                  ' --qualified_quality_phred ',qualified_quality_phred,' --thread ',threads)
    cmd <- paste(cmd,'2>&1')
    system(cmd, intern=TRUE)

    # EXIT 0
    msg <- 'fastplong executed successfully.'
    message(msg)
}

#' FastP filtering of short reads
#'
#' CLI wrapper for FastP to filter and trim short reads
#'
#' @param input_read_1 Character, Input file name.
#' @param output_read_1 Character, Output file name.
#' @param input_read_2 Character, Input file name.
#' @param output_read_2 Character, Output file name.
#' @param qualified_quality_phred Numeric, the quality value that a base is qualified.
#' @param unqualified_percent_limit Numeric, how many percents of bases are allowed to be unqualified (0~100).
#' @param length_required Numeric, shorter reads will be discarded.
#' @param average_quality Numeric, if one read's average quality score is lower, then this read/pair is discarded.
#' @param threads Integer, number of worker threads
#' @param overwrite Boolean, whether to overwrite output.
#'
#' @export
fastp <- function(input_read_1 = NULL,
                  output_read_1 = NULL,
                  input_read_2 = NULL,
                  output_read_2 = NULL,
                  qualified_quality_phred = 15,
                  unqualified_percent_limit = 40,
                  length_required = 15,
                  average_quality = 0,
                  threads = n_proc(),
                  overwrite = FALSE
                 ) {

    # Minimal check
    stopifnot(
        !is.null(input_read_1),
        !is.null(output_read_1)
    )
    check_installed('fastp', silent=TRUE)
    check_version('fastp')

    # Check input
    paired <- !is.null(input_read_2)
    if (paired & is.null(output_read_2)) {
        msg <- 'Output file name for read 2 missing.'
        stop(msg)
    }

    # Check output
    if (file.exists(output_read_1)) {
        msg <- paste('Output file',output_read_1,'already exists.')
        if (overwrite) {
            msg <- paste(msg, 'Will be overwritten...')
            warning(msg)
        } else {
            return('Exit 2: Output exists.')
        }
    }
    output.dir <- dirname(output_read_1)
    output.json <- paste0(output.dir,'/',str_replace(basename(output_read_1),'.fastq.gz','.json'))
    output.html <- paste0(output.dir,'/',str_replace(basename(output_read_1),'.fastq.gz','.html'))

    # Run FastP
    if (paired) {
        cmd <- paste('fastp','-i',input_read_1,'-I',input_read_2,'-o',output_read_1,'-O',output_read_2)
    } else {
        cmd <- paste('fastp','-i',input_read_1,' -o ',output_read_1)
    }
    cmd <- paste(cmd,
                 '--qualified_quality_phred',qualified_quality_phred,
                 '--unqualified_percent_limit',unqualified_percent_limit,
                 '--length_required',length_required,
                 '--average_qual',average_quality,
                 '--json',output.json,
                 '--html',output.html
                )
    cmd <- paste(cmd,'2>&1')
    stdout <- system(cmd, intern=TRUE)
    cat(paste(stdout, collapse='\n'))
}