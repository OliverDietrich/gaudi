#' ReadsCounts
#'
#' Count FASTQ files for a genomeCollection
#'
#' @param object genomeCollection
#' @param name Slot name of Reads object in genomeCollection
#' @param threads Integer, number of threads to use
#'
#' @export
#'
ReadsCounts <- function(object, name, recompute=FALSE, multicore=TRUE, threads=n_proc()) {

    # Minimal check
    stopifnot(
        name %in% ReadsNames(object)
    )

    # Variables & data
    data <- Reads(object, name)
    csv <- paste0(path(object), 'reads/', name, '_counts.csv')
    if (!dir.exists(dirname(csv))) {
        dir.create(dirname(csv))
    }
    id <- if (length(data$alias)) 'alias' else 'index'
    if (id == 'alias') {
        lookup <- setNames(data$index, data$alias)
    }

    # Check output
    df <- data$counts
    if (recompute) {
        compute <- TRUE
    } else
    if (!length(df) & file.exists(csv)) {
        msg <- paste('File',csv,'found in path(object)')
        message(msg)
        df <- readr::read_csv(csv)
        compute <- FALSE
    } else
    if (length(df) & all(data$index == df$index)) {
        msg <- 'Read counts have been detected for all samples.'
        message(msg)
        compute <- FALSE
    } else {
        compute <- TRUE
    }

    # Run countFastq
    if (compute) {
        # Timestamp
        time_start <- Sys.time()
        
        msg <- paste0('Counting FASTQ entries for Reads(object, "', name, '").')
        message(msg)

        df <- list()
        for (i in c('R1','R2','S','L')) {
            files <- slot(data, i)
            ind <- file.exists(files)
            files <- files[ind]
            if (!length(files)) next

            # Set names
            names(files) <- slot(data, id)[ind]
            
            # Estimate runtime
            single.time.start <- Sys.time()
            discard <- countFastq(files[[1]])
            single.time.stop <- Sys.time()
            time.diff <- single.time.stop-single.time.start
            if (multicore) {
                time.est <- time.diff*length(files)/threads
            } else {
                time.est <- time.diff*length(files)
            }
            msg <- paste('Estimated runtime for slot', i, 'is', time.est, units(time.est))
            message(msg)

            # Apply to all
            if (multicore) {
                df[[i]] <- parallel::mclapply(files, ShortRead::countFastq, mc.cores = threads)
            } else {
                df[[i]] <- lapply(files, ShortRead::countFastq)
            }
            df[[i]] <- dplyr::bind_rows(df[[i]], .id = id)
        }

        # Combine across read types
        df <- dplyr::bind_rows(df, .id = 'type')
        if (length(df$alias)) {
            df$index <- lookup[df$alias]
        }
        
        # Save counts
        readr::write_csv(df, csv)

        ## Timestamp
        time_stop <- Sys.time()
        time.diff <- time_stop-time_start
        msg <- paste('Time elapsed:', time.diff, units(time.diff), 'on', threads, 'cores.')
        message(msg)
    }

    # Replace
    slot(Reads(object, name), 'counts') <- as.data.frame(df)

    return(object)
}

#' PlotReadsCounts
#'
#' Visualize FASTQ counts
#' 
#' @param object genomeCollection
#' @param name Slot name of Reads object in genomeCollection
#'
#' @importFrom ggplot2 ggplot aes theme labs geom_point theme_classic
#'
#' @export
#' 
PlotReadsCounts <- function(object, name, base.size=20, label.size = 5, grid.cols=1,
                            top.point.size = 4, top.point.stroke = 1, jitter=FALSE,
                            bottom.point.size = 3, bottom.point.stroke = 1, bottom.text.size = 10
                           ) {

    # Collect data
    df <- slot(Reads(object, name), 'counts')

    # Minimal check
    stopifnot(
        !is.null(df)
    )

    # Translate type
    lookup <- c('R1' = 'Short reads (forward)', 'R2' = 'Short reads (reverse)', 'S' = 'Short reads (unpaired)', 'L' = 'Long reads')
    df$type <- factor(lookup[df$type], lookup)

    # Set colors
    cols <- list()
    cols$type <- RColorBrewer::brewer.pal(4, 'Set1')
    names(cols$type) <- lookup

    # Modify scales
    scale_jitter <- if (jitter) ggplot2::position_jitter(width=0.25) else ggplot2::position_identity()

    # Select label
    df$label <- if (length(df$alias)) df[['alias']] else df[['index']]

    # Plot
    p1 <- ggplot(data = df, mapping = aes(x=records, y=nucleotides, label=label)) +
                ggrepel::geom_text_repel(force_pull = 0, force = 25, alpha=.8, size=label.size) +
                geom_point(aes(col = type), shape=21, size=top.point.size, stroke=top.point.stroke) +
                ggplot2::scale_color_manual(values = cols$type) +
                ggplot2::scale_y_continuous(trans='log10') +
                ggplot2::scale_x_continuous(trans='log10') +
                theme_classic(base.size) +
                theme(
                  panel.grid.major = ggplot2::element_line(color = 'grey', linewidth = .25),
                  panel.grid.minor = ggplot2::element_line(color = 'grey', linewidth = .1)
                ) +
                labs(x='Number of reads', y='Nucleotides')

    # Re-order
    df <- df[order(df$nucleotides, decreasing=TRUE), ]
    df$index <- factor(df$index, unique(df$index))
    
    # Plot
    p2 <- ggplot(data = df, mapping = aes(x = index, y = nucleotides, label=label)) +
                ggrepel::geom_text_repel(force_pull = 0, force = 1, alpha=.8, size=label.size) +
                geom_point(aes(col = type), shape = 21, size = bottom.point.size, stroke=bottom.point.stroke, position=scale_jitter) +
                ggplot2::scale_color_manual(values = cols$type) +
                ggplot2::scale_y_continuous(trans='log10') +
                theme_classic(base.size) +
                theme(
                    axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1, size=bottom.text.size),
                    panel.grid.major.y = ggplot2::element_line()
                ) +
                labs(x='Sample', y='Nucleotides')

    plot <- cowplot::plot_grid(p1, p2, ncol=grid.cols)

    return(plot)
}

#' Compute read quality
#'
#' Compute average PHRED score for a single read within a FASTQ file
#'
#' @param fastq Character, path to FASTQ file
#' @param read Integer, which read to summarize
#' @param min.base.quality Integer, PHRED score of individual bases to be qualified
#' @param per.base.quality Boolean, whether to return per-base quality scores (TRUE) or summary scores (FALSE)
#'
#' @export
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
#' Compute average PHRED score across a subset of reads within a FASTQ file
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

#' ReadsQuality
#'
#' Summarize quality of FASTQ files in a genomeCollection
#'
#' @param object genomeCollection
#' @param read.col Column in metadata(object) containing FASTQ file paths
#' @param threads Integer, number of threads to use
#'
ReadsQuality <- function(object, name, max.reads=5000, min.base.quality=15, 
                         recompute=FALSE, threads=n_proc()
                        ) {

    # Minimal check
    stopifnot(
        name %in% ReadsNames(object)
    )

    # Variables & data
    data <- Reads(object, name)
    csv <- paste0(path(object), 'reads/', name, '_quality.csv')
    if (!dir.exists(dirname(csv))) {
        dir.create(dirname(csv))
    }
    id <- if (length(data$alias)) 'alias' else 'index'
    if (id == 'alias') {
        lookup <- setNames(data$index, data$alias)
    }

    # Check output
    df <- data$quality
    if (recompute) {
        compute <- TRUE
    } else
    if (!length(df) & file.exists(csv)) {
        msg <- paste('File',csv,'found in path(object)')
        message(msg)
        df <- readr::read_csv(csv)
        compute <- FALSE
    } else
    if (length(df) & all(data$index %in% df$index)) {
        msg <- 'Read quality summary has been detected for all samples.'
        message(msg)
        compute <- FALSE
    } else {
        compute <- TRUE
    }


    # Run countFastq
    if (compute) {
        # Timestamp
        time_start <- Sys.time()
        
        msg <- paste0('Summarizing FASTQ quality for Reads(object, "', name, '").')
        message(msg)

        df <- list()
        for (i in c('R1','R2','S','L')) {
            files <- slot(data, i)
            ind <- file.exists(files)
            files <- files[ind]
            if (!length(files)) next

            # Set names
            names(files) <- slot(data, id)[ind]
            
            # Estimate runtime
            single.time.start <- Sys.time()
            discard <- summarize_fastq_quality(files[[1]], min.base.quality=min.base.quality, max.reads=max.reads)
            single.time.stop <- Sys.time()
            time.diff <- single.time.stop-single.time.start
            time.est <- time.diff*length(files)/threads
            msg <- paste('Estimated runtime for slot', i, 'is', time.est, units(time.est), 'on', threads, 'cores.')
            message(msg)

            # Apply for all
            df[[i]] <- parallel::mclapply(files, summarize_fastq_quality, min.base.quality=min.base.quality, max.reads=max.reads, mc.cores = threads)
            df[[i]] <- dplyr::bind_rows(df[[i]], .id = id)
        }

        # Combine across read types
        df <- dplyr::bind_rows(df, .id = 'type')
        if (length(df$alias)) {
            df$index <- lookup[df$alias]
        }

        # Save counts
        readr::write_csv(df, csv)

        ## Timestamp
        time_stop <- Sys.time()
        time.diff <- time_stop-time_start
        msg <- paste('Time elapsed:', time.diff, units(time.diff), 'on', threads, 'cores.')
        message(msg)
    }

    # Replace
    slot(Reads(object, name), 'quality') <- as.data.frame(df)

    return(object)
}

#' PlotReadsQuality
#'
#' View FASTQ quality
#' 
#' @param object genomeCollection
#' @param read.col Column in metadata(object) containing FASTQ paths
#'
#' @importFrom ggplot2 ggplot aes theme labs geom_point theme_classic
#'
#' @export
#' 
PlotReadsQuality <- function(object, name, x='length', y='avg_quality', col='type', wrap='index', 
                             L.lenght_required = 100L, L.average_quality = 20L, L.percent_qualified = 50L,
                             S.lenght_required = 50L, S.average_quality = 30L, S.percent_qualified = 50L,
                             base.size=20, pt.size=.5, pt.shape=21, pt.stroke=.5,
                             wrap.rows=NULL, wrap.cols=NULL
                            ) {

    # Collect data
    df <- slot(Reads(object, name), 'quality')

    # Minimal check
    stopifnot(
        !is.null(df)
    )

    # Translate type
    lookup <- c('R1' = 'Short reads (forward)', 'R2' = 'Short reads (reverse)', 'S' = 'Short reads (unpaired)', 'L' = 'Long reads')
    df$type <- factor(lookup[df$type], lookup)

    # Set index levels
    df$index <- factor(df$index, index(object))
    
    # Set axes
    df$x <- df[[x]]
    df$y <- df[[y]]
    df$col <- df[[col]]
    df$wrap <- df[[wrap]]
    df$wrap <- factor(df$wrap, unique(df$wrap[order(df$index)]))

    # Set colors
    cols <- list()
    cols$type <- RColorBrewer::brewer.pal(4, 'Set1')
    names(cols$type) <- lookup

    # Scales
    x_scale <- if (x %in% c('length','avg_quality')) ggplot2::scale_x_continuous(trans='log10') else NULL
    y_scale <- if (y %in% c('length','avg_quality')) ggplot2::scale_y_continuous(trans='log10') else NULL

    # Lines
    thresh <- data.frame(
        type = c('R1','R2','S','L'),
        size = c('short','short','short','long'),
        length = c(S.lenght_required, S.lenght_required, S.lenght_required, L.lenght_required),
        avg_quality = c(S.average_quality, S.average_quality, S.average_quality, L.average_quality),
        percent_qualified = c(S.percent_qualified, S.percent_qualified, S.percent_qualified, L.percent_qualified)
    )
    thresh$type <- factor(lookup[thresh$type], lookup)
    thresh <- thresh[thresh$type %in% df$type, ]
    thresh$wrap <- thresh[[wrap]]
    thresh$quality <- thresh[[y]]

    # Plot
    plot <- ggplot(df, aes(x = x, y = y, col = col)) +
                ggplot2::geom_point(size = pt.size, shape=pt.shape, stroke=pt.stroke) +
                ggplot2::geom_vline(data = thresh, aes(xintercept = length, linetype=size)) + ggplot2::geom_hline(data = thresh, aes(yintercept = quality, linetype=size)) +
                ggplot2::facet_wrap(~wrap, nrow=wrap.rows, ncol=wrap.cols) +
                x_scale + y_scale +
                ggplot2::scale_color_manual(values = cols$type) +
                ggplot2::theme_classic(base.size) +
                ggplot2::theme(
                    legend.position = 'top',
                    axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1),
                    panel.grid.major.x = element_line(),
                    panel.grid.major.y = element_line(),
                ) +
                ggplot2::guides(
                    col = ggplot2::guide_legend(override.aes = list(size = 5, stroke=1))
                ) +
                ggplot2::labs(x = x, y = y, col = col, linetype = 'filter')

    return(plot)
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
                      length_required=100L, 
                      mean_qual=20L,
                      unqualified_percent_limit=40L, 
                      qualified_quality_phred = 15L,
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
    cmd <- paste0(cmd,
                  ' --length_required ',length_required,
                  ' --mean_qual ',mean_qual,
                  ' --unqualified_percent_limit ',unqualified_percent_limit,
                  ' --qualified_quality_phred ',qualified_quality_phred,
                  ' --thread ',threads
                 )
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
                  length_required = 40,
                  average_quality = 20,
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

#' FilterReads
#'
#' Filter FASTQ reads using fastp/fastplong
#'
#' @param object genomeCollection
#' @param name.from Name of the Reads object storing unfiltered FASTQ files
#' @param name.to Name of the Reads object to create for filtered FASTQ files
#' @param fastp.length_required Numeric, see ?fastp()
#' @param fastp.average_quality Numeric, see ?fastp()
#' @param fastp.unqualified_percent_limit Numeric, see ?fastp()
#' @param fastp.qualified_quality_phred Numeric, see ?fastp()
#' @param fastplong.length_required Numeric, see ?fastplong()
#' @param fastplong.mean_qual Numeric, see ?fastplong()
#' @param fastplong.unqualified_percent_limit Numeric, see ?fastplong()
#' @param fastplong.qualified_quality_phred Numeric, see ?fastplong()
#'
#' @export
#'
FilterReads <- function(object, name.from = 'raw', name.to = 'filtered',
                        fastp.length_required = 40L, 
                        fastp.average_quality = 20L, 
                        fastp.unqualified_percent_limit = 40L,
                        fastp.qualified_quality_phred = 15L, 
                        fastplong.length_required = 100L, 
                        fastplong.mean_qual = 20L, 
                        fastplong.unqualified_percent_limit = 40L, 
                        fastplong.qualified_quality_phred = 15L
                       ) {

    # Minimal check
    if (name.to %in% ReadsNames(object)) {
        msg <- paste0('Reads(object, "',name.to,'") already exists. Aborting...')
        warning(msg)
        return(object)
    }
    
    # Collect data
    data <- Reads(object, name.from)
    types <- c('R1','R2','S','L')

    # Concatenate multiple files of same type
    if (length(data$alias)) {
        msg <- paste0('Alias present in Reads(object, "', name.from, '"). Multiple read files of the same sample will be concatenated...')
        warning(msg)

        samples <- data$index[which(data$alias != data$index)]
        ind.keep <- which(data$alias == data$index)
        
        for (type in types) {
            vector <- slot(data, type)
            if (!length(vector)) next # Skip empty vectors
            for (sample in samples) {
                cat(sample, type, '\n')
                ind <- which(data$index == sample)
                old.files <- slot(data, type)[ind]
                suffix <- unique(stringr::str_split(old.files, '\\.', simplify=TRUE)[, 2])
                if (length(suffix) > 1) {
                    msg <- paste('Trying to concatenate files of different suffix:', paste(suffix, collapse=', '))
                    stop(msg)
                }
                old.files <- paste0(old.files, collapse=' ')
                new.file <- paste0(path(object),'assemblies/',sample,'/reads/raw/',type,'.',suffix)
                dir.create(dirname(new.file), recursive=TRUE, showWarnings = FALSE)
                cmd <- paste('cat',old.files,'>',new.file)
                system3(cmd)
                ind.replace <- which(data$alias == sample)
                vector[ind.replace] <- new.file
            }
            # Replace components
            slot(data, type) <- vector[ind.keep]
        }
        # Replace components
        slot(data, 'index') <- data$index[ind.keep]
        slot(data, 'alias') <- character()
        Reads(object, name.from) <- data
    }

    # Create PATHS for filtered FASTQ files
    transfer.ind <- sapply(lapply(types, slot, object = data), length) > 0
    transfer.types <- types[transfer.ind]
    df <- data.frame('index' = data$index)
    for (type in transfer.types) {
        df[[type]] <- paste0(path(object),'assemblies/', df$index, '/reads/filtered/',type,'.fastq')
    }

    # Run FastP(long)
    if (length(df$R1) & length(df$R2)) {
        message('Running fastp for paired, short reads (R1 + R2).')
        for (n in 1:length(df$index)) {
            from.1 <- data$R1[[n]]
            from.2 <- data$R2[[n]]
            to.1 <- df$R1[[n]]
            to.2 <- df$R2[[n]]
            if (!file.exists(fromfile)) { # Needed to not break hybrid assemblies for some samples ...
                df$R1[[n]] <- NA
                df$R2[[n]] <- NA
                next
            }
            if (file.exists(to.1)) next
            print(paste0(data$index[[n]],' (',n,'/',length(data$index),')'))
            # Run fastp
            dir.create(dirname(to.1), recursive=TRUE, showWarnings=FALSE)
            fastp(from.1, to.1, from.2, to.2,
                  length_required = fastp.length_required, average_quality = fastp.average_quality, 
                  unqualified_percent_limit = fastp.unqualified_percent_limit, 
                  qualified_quality_phred = fastp.qualified_quality_phred
                 )
        }
    }
    
    if (length(df$S)) {
        stop('Running fastp for unpaired, short reads (S) IS NOT IMPLEMENTED YET!!!')
    }
    
    if (length(df$L)) {
        message('Running fastplong for long reads (L).')
        for (n in 1:length(df$index)) {
            fromfile <- data$L[[n]]
            tofile <- df$L[[n]]
            if (!file.exists(fromfile)) { # Needed to not break hybrid assemblies for some samples ...
                df$L[[n]] <- NA
                next
            }
            if (file.exists(tofile)) next
            print(paste0(data$index[[n]],' (',n,'/',length(data$index),')'))
            # Run fastplong
            dir.create(dirname(tofile), recursive=TRUE, showWarnings=FALSE)
            fastplong(fromfile, tofile, 
                      length_required = fastplong.length_required, mean_qual = fastplong.mean_qual, 
                      unqualified_percent_limit = fastplong.unqualified_percent_limit, 
                      qualified_quality_phred = fastplong.qualified_quality_phred
                     )
        }
    }

    # Create Reads object
    df <- as.data.frame(df)
    object <- AddReads(object, data = df, name = name.to)

    return(object)
}