#' BLAST
#'
blast <- function() {

    cmd <- paste(
        'blastn',
        '-task','blastn',
        '-query',input['ps.009'],
        '-strand','both',
        '-subject',input['ps.021'],
        '-outfmt','7',
        '-out','tmp/blastn.aln'
    )
    system3(cmd)
}

#' Minimap2
#'
#' Minimap2 alignment of FASTQ reads to FASTA reference.
#'
#' @param reference FASTA file
#' @param out.dir Path to output directory
#' @param short.1 FASTQ file of first short reads in each pair (R1)
#' @param short.2 FASTQ file of second short reads in each pair (R2)
#' @param unpaired FASTQ file of unpaired short reads (S)
#' @param long FASTQ or FASTA file of long reads
#' @param tool Name of alignment tool (minimap2 or bowtie2)
#' @param preset Preset for tuning the minimap 2 params. One of 'map-ont', 'map-pb', or 'map-ilcr' for long reads. 
#' Will default to 'sr' for short reads.
#' @param threads Number of threads
#'
#' @export
#'
map_reads_to_reference <- function(reference, out.dir,
                                   short.1=NULL, short.2=NULL, unpaired=NULL, long=NULL,
                                   tool = 'minimap2',
                                   preset = 'map-ont',
                                   file.prefix = NULL,
                                   remove.intermediates = TRUE,
                                   recompute = FALSE,
                                   threads = n_proc()
                                  ) {

    # Input
    dir.create(out.dir, recursive=TRUE, showWarnings=FALSE)
    read_1_exists <- if (is.null(short.1)) FALSE else file.exists(short.1);
    read_2_exists <- if (is.null(short.2)) FALSE else file.exists(short.2);
    paired <- if (read_1_exists & read_2_exists) TRUE else FALSE;
    read_S_exists <- if (is.null(unpaired)) FALSE else file.exists(unpaired);
    read_L_exists <- if (is.null(long)) FALSE else file.exists(long);
    index <- c('paired'=paired, 'unpaired'=read_S_exists,'long'=read_L_exists)
    if (sum(index) == 0) {
        stop('No reads (pairs) found. Aborting...')
    } else if (sum(index) == 1) {
        prefix <- names(index)[index]
        msg <- paste('Found',prefix,'reads. Mapping to',reference,'...')
        message(msg)
    } else if (sum(index) > 1) {
        prefix <- paste(names(index)[index], collapse=' and ')
        msg <- paste('Found', prefix, 'reads. Please specify a single type.')
        stop(msg)
    }

    # Variables
    preset <- if (prefix %in% c('paired','unpaired')) 'sr' else preset;
    prefix <- if (!is.null(file.prefix)) file.prefix else prefix;
    out.dir <- if(endsWith(out.dir,'/')) out.dir else paste0(out.dir,'/');
    base.name <- paste0(out.dir,prefix,'-',tool)
    out.sam <- paste0(base.name,'.','sam')
    out.bam <- paste0(base.name,'.','bam')
    out_sorted.bam <- paste0(base.name,'_sorted','.','bam')
    out.bai <- paste0(base.name,'_indexed','.','bai')
    out_mapped.bam <- paste0(base.name,'_mapped','.','bam')
    out.report <- paste0(base.name,'_report','.','txt')
    out.bcf <- paste0(base.name,'_variants','.','bcf')
    out.vcf <- paste0(base.name,'_variants','.','vcf')
    out.coverage <- paste0(base.name,'_coverage','.','tsv')

    # Alignment
    if (tool == 'minimap2') {
        fastq <- if (paired) paste(short.1, short.2) else if (unpaired) unpaired else if (long) long;
        cmd <- paste('minimap2','-a','-x',preset,reference,fastq,'-o',out.sam)
    } else if (tool == 'bowtie2') {
        fastq <- if (paired) {
            paste('-1',short.1,'-2',short.2)
        } else if (unpaired) {
            unpaired
        } else if (long) {
            stop('Bowtie for long reads is not supported.')
        }
        cmd <- paste('bowtie2','-x',out.dir,fastq)
        stop('No implemented yet.')
    } else {
        stop('Unknown alignment tool selected. Please specify minimap2 or bowtie2.')
    }
    if (!file.exists(out.sam) | recompute) {
        message('Running alignment...')
        system3(cmd)
    } else {
        message('Alignment file (SAM) already exists.')
        cat(cmd)
    }

    # Convert to BAM
    # samtools view -bT reference.fa test.sam > test.bam #(if header absent)
    # samtools view -bS test.sam > test.bam #(if header present)
    cmd <- paste('samtools','view','-bT',reference,out.sam,'>',out.bam)
    if (!file.exists(out.bam) | recompute) {
        message('Converting SAM to BAM...')
        system3(cmd)
    } else {
        message('Alignment file (BAM) already exists.')
        cat(cmd)
    }

    # Sort BAM
    # samtools sort test.bam test_sorted
    cmd <- paste('samtools','sort',out.bam,'-o',out_sorted.bam)
    if (!file.exists(out_sorted.bam) | recompute) {
        message('Sorting BAM file...')
        system3(cmd)
    } else {
        message('BAM file has already been sorted.')
        cat(cmd)
    }

    # Index BAM
    # samtools index test_sorted.bam test_sorted.bai
    cmd <- paste('samtools','index',out_sorted.bam,'-o',out.bai)
    if (!file.exists(out.bai) | recompute) {
        message('Indexing BAM file...')
        system3(cmd)
    } else {
        message('BAM file has already been indexed.')
        cat(cmd)
    }

    # Filter BAM
    # samtools view -h -F 4 blah.bam > blah_only_mapped.sam
    cmd <- paste('samtools','view','-h','-F','4',out_sorted.bam,'-o',out_mapped.bam)
    if (!file.exists(out_mapped.bam) | recompute) {
        system3(cmd)
    } else {
        message('BAM file has already been filtered.')
        cat(cmd)
    }    

    # Create report
    cmd <- paste('samtools','stats',out_sorted.bam,'--reference',reference,'--threads',threads,'>',out.report)
    if (!file.exists(out.report) | recompute) {
        system3(cmd, include.errors = FALSE)
    } else {
        message('Alignment report already exists.')
        cat(cmd)
    }
    
    # Call variants
    # bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
    cmd <- paste('bcftools','mpileup','-Ou','-f',reference,out_mapped.bam,'|','bcftools','call','-mv','-Ob','-o',out.bcf)
    if (!file.exists(out.bcf) | recompute) {
        message('Calling variants...')
        system3(cmd)
    } else {
        message('Variant calling file (VCF) present.')
        cat(cmd)
    }

    # Convert BCF -> VCF
    cmd <- paste('bcftools','view',out.bcf,'>',out.vcf)
    if (!file.exists(out.vcf) | recompute) {
        message('Converting BCF to VCF...')
        system3(cmd, include.errors = FALSE)
    } else {
        message('VCF file present.')
        cat(cmd)
    }

    # Compute coverage
    cmd <- paste('samtools','depth',out_sorted.bam,'-o',out.coverage)
    if (!file.exists(out.coverage) | recompute) {
        message('Computing coverage...')
        system3(cmd)
    } else {
        message('Coverage already computed.')
        cat(cmd)
    }

    # Cleanup
    remove.intermediates <- FALSE # During DEVELOPMENT
    if (remove.intermediates) {
        message('Removing (large) intermediates...')
        unlink(out.sam)
        unlink(out.bam)
        unlink(out_sorted.bam)
        unlink(out.bai)
    }

    # Exit
    message('Done.')
}

#' Create MMseqs database
#'
#' @param files Path to (multiple) file name(s)
#' @param output.file Path to output database
#' @param tmp.dir Path to temporary directory. Important for some HPC systems for I/O load to set to /tmp subdirectories.
#' @param create.index Whether to create an index file. Takes time, only recommended if database is used multiple times.
#' @param overwrite Whether to remove and overwrite output.file
#' @param threads Number of threads to use
#'
mmseqs_create_database <- function(files, output.file, tmp.dir = tempdir(), 
                                   create.index = FALSE,
                                   overwrite = FALSE, 
                                   threads = n_proc(),
                                   verbosity = 3,
                                   clean = TRUE
                                  ) {

    # Program
    check_installed('mmseqs')

    # Check input
    if (file.exists(output.file) & overwrite) file.remove(output.file)
    if (endsWith(output.file, '/')) stop(paste('Output file', output.file, 'looks like a directory.'))
    if (file.exists(output.file)) {
        message(paste('Output file', output.file, 'already exists.'))
        return()
    }
    print_command <- TRUE
    if (length(files) > 5) {
        print_command <- FALSE
        message('Multiple input files detected. MMseqs call will not be shown...')
    }

    # Handle large argument lists    
    if (length(files) > 1000) {
        message('Found > 1000 files. Concatenating in chunks of 1000.')
        files.fasta <- paste0(tmp.dir, 'input_files.fasta')
        n_groups <- ceiling(length(files) / 1000)
        ind <- cut(seq_along(files), n_groups, labels = FALSE)
        files_list <- split(files, ind)
        file.create(files.fasta) # Create empty file
        for (n in 1:length(files_list)) { # Fill
            cmd <- paste(
                'cat',paste0(files_list[[n]], collapse = ' '),'>>',files.fasta
            )
            system3(cmd, verbose = FALSE)
        }
        files <- files.fasta # Re-assign input file
    }

    # Main
    cmd <- paste(
        'mmseqs','createdb',
        paste(files, collapse = ' '),
        output.file,
        #'--threads',threads, # BUG, results in mmseqs runtime error!
        '-v',verbosity
    )
    system3(cmd, verbose = print_command)

    # Index
    cmd <- paste(
        'mmseqs','createindex',output.file,tmp.dir,
        '--threads',threads,
        '-v',verbosity
    )
    if (create.index) {
        system3(cmd)
    }

    # Add another file to the database 

    # Check output
    if (!is_mmseqs_database(output.file, verbose = TRUE)) stop('Database formatting issues detected. Aborting...')

    # Cleanup
    if (clean) {
        unlink(tmp.dir, recursive = TRUE)
        dir.create(tmp.dir) # If the temp.dir() is missing, other stuff breaks!
    }
    
    # Exit 0
    message('Done.')
}

#' Check whether file is a mmseqs database
#'
#' @param file
#'
is_mmseqs_database <- function(file, verbose = FALSE) {

    # Input checks
    if (length(file) > 1) {
        if (verbose) warning('Single file expected, not multiple...')
        return(FALSE)
    }
    if (!file.exists(file)) {
        if (verbose) warning('File not found. Aborting...')
        return(FALSE)
    }

    # Variables
    DBTYPE <- paste0(file, '.dbtype')
    HEADER <- paste0(file, '_h')
    HEADER_INDEX <- paste0(file, '_h.index')
    INDEX <- paste0(file, '.index')
    LOOKUP <- paste0(file, '.lookup')
    SOURCE <- paste0(file, '.source')

    DB_INDEX <- paste0(file, '.idx')
    DB_INDEX_DBTYPE <- paste0(file, '.idx.dbtype')
    DB_INDEX_INDEX <- paste0(file, '.idx.index')

    # Check regular file presence
    files_exist <- file.exists(c(DBTYPE, HEADER, HEADER_INDEX, INDEX, LOOKUP, SOURCE))
    if (!all(files_exist)) {
        if (verbose) warning('Incomplete mmseqs database. Aborting...')
        return(FALSE)
    }

    return(TRUE)
}

#' MMseqs search
#'
#' Sensitive homology search of protein sequences using MMseqs2
#'
#' @param query Input query sequences (.faa, .fasta) or mmseqs database (.db). File suffix must be exact match.
#' @param reference Input query sequences (.faa, .fasta) or mmseqs database (.db). File suffix must be exact match.
#'
mmseqs_search <- function(query, reference, output.file, 
                          tmp.dir = tempdir(), 
                          overwrite = FALSE, 
                          sensitivity = 5.7,
                          max.seqs = 300,
                          start.sens = 4,
                          sens.steps = 1,
                          alignment.mode = 0,
                          comp.bias.corr = 1,
                          filter.hits = 0,
                          threads = n_proc(),
                          verbosity = 3,
                          clean = FALSE
                         ) {
    
    # Program
    check_installed('mmseqs')

    # Variables
    cn <- c('QuerySetId','TargetSetId','Query','Target','Identity')
    output.dir <- dirname(output.file)
    output.dir <- if (endsWith(output.dir, '/')) output.dir else paste0(output.dir, '/');
    output.dir_notrail <- stringr::str_replace(output.dir, '/$', '')
    tmp.dir <- if (endsWith(tmp.dir, '/')) tmp.dir else paste0(tmp.dir, '/');
    tmp.dir_notrail <- stringr::str_replace(tmp.dir, '/$', '')
    seqDB <- paste0(tmp.dir, 'seqDB')
    query_names <- names(query)
    reference_names <- names(reference)

    # Check input
    if (file.exists(output.file) & overwrite) file.remove(output.file)
    if (file.exists(output.file)) {
        warning(paste('Output file', output.file, 'already exists. Returning previous results...'))
        result <- readr::read_tsv(output.file, col_names = cn)
        return(result)
    }

    # Set up directories
    dir.create(output.dir, recursive = TRUE)
    
    # Query    
    if (is_mmseqs_database(query)) {
        queryDB <- query
    } else if (is_fasta(query)) {
        queryDB <- paste0(tmp.dir, 'queryDB')
        mmseqs_create_database(query, queryDB)
    } else {
        stop('Query is neither protein FASTA nor mmseqs database. Aborting...')
    }
    
    # Reference
    if (is_mmseqs_database(reference)) {
        targetDB <- reference
    } else if (is_fasta(reference)) {
        targetDB <- paste0(tmp.dir, 'targetDB')
        mmseqs_create_database(reference, targetDB)
    } else {
        stop('Reference is neither protein FASTA nor mmseqs database. Aborting...')
    }

    # Simplify data input by writing a function that checks mmseqs input formatting...
    # Must have tmp.dir passed to it though.

    # Search
    cmd <- paste(
        'mmseqs','search',queryDB,targetDB,seqDB,tmp.dir_notrail,
        '-s',sensitivity,
        '--max-seqs',max.seqs,
        '--start-sens',start.sens,
        '--sens-steps',sens.steps,
        '--alignment-mode',alignment.mode,
        '--comp-bias-corr',comp.bias.corr,
        '--filter-hits',filter.hits,
        '--threads',threads,
        '-v',verbosity
    )
    system3(cmd, verbose = FALSE)

    # Convert output
    cmd <- paste(
        'mmseqs','convertalis',
        queryDB,targetDB,seqDB,output.file,
        '--format-output','qsetid,tsetid,query,target,fident'
    )
    system3(cmd, verbose = FALSE)

    # Read
    result <- readr::read_tsv(output.file, col_names = cn)

    # Format
    # I want to add the name (e.g. phage) for each input file. 

    # For mmseqs databases, I would need to write an additional file (e.g. sample_header)
    # so I can read that and index the query.sample_header by QuerySetId
    
    # This code works for FASTA input (with or without names)
    result$QuerySet <- names(query)[result$QuerySetId+1]
    result$TargetSet <- names(reference)[result$TargetSetId+1]

    # Cleanup
    if (clean) {
        unlink(tmp.dir, recursive = TRUE)
        dir.create(tmp.dir) # If the temp.dir() is missing, other stuff breaks!
    }

    # Exit
    return(result)
    message('Done.')
}

#' MMseqs clustering of protein sequences
#'
#' @param input Path to input FASTA file(s)
#' @param output.file Path to output file (TSV)
#' @param sensitivity Numeric, 1.0 (faster); 4.0 (fast); 7.5 (sensitive)
#'
mmseqs_cluster <- function(input, output.file, 
                           cluster.mode = 0,
                           sensitivity = 4,
                           min.seq.id = 0,
                           cov.mode = 0,
                           tmp.dir = tempdir(), 
                           overwrite = FALSE,
                           clean = TRUE
                          ) {

    # Program
    check_installed('mmseqs')

    # Variables
    cn <- c('QuerySetId','TargetSetId','Query','Target','Identity')
    output.dir <- dirname(output.file)
    output.dir <- if (endsWith(output.dir, '/')) output.dir else paste0(output.dir, '/');
    output.dir_notrail <- stringr::str_replace(output.dir, '/$', '')
    tmp.dir <- if (endsWith(tmp.dir, '/')) tmp.dir else paste0(tmp.dir, '/');
    tmp.dir_notrail <- stringr::str_replace(tmp.dir, '/$', '')
    DB_clu <- paste0(output.dir, 'DB_clu')

    # Check input
    if (overwrite) file.remove(output.file)
    if (file.exists(output.file)) {
        warning(paste('Output file', output.file, 'already exists. Returning previous results...'))
        result <- readr::read_tsv(output.file, col_names = cn)
        return(result)
    }

    # Set up directories
    dir.create(output.dir, recursive = TRUE)

    # Query    
    if (is_mmseqs_database(input)) {
        seqDB <- query
    } else if (is_fasta(input)) {
        seqDB <- paste0(tmp.dir, 'seqDB')
        mmseqs_create_database(input, seqDB)
    } else {
        stop('Input is neither protein FASTA nor mmseqs database. Aborting...')
    }

    # Main
    cmd <- paste(
        'mmseqs','cluster',seqDB,DB_clu,tmp.dir_notrail,
        'cluster-mode',cluster.mode,
        '-s',sensitivity,
        '--min-seq-id',min.seq.id,
        '--cov-mode',cov.mode
    )
    system3(cmd)

    # Reformat
    cmd <- paste(
        'mmseqs','createtsv',seqDB,seqDB,DB_clu,output.file
    )
    system3(cmd)

    # Read
    result <- readr::read_tsv(output.file)

    # Cleanup
    if (clean) {
        unlink(tmp.dir, recursive = TRUE)
        dir.create(tmp.dir) 
    }

    # Exit
    return(result)
    message('Done.')
}

#' MUMmer alignment
#' 
#' @param query FASTA file
#' @param reference FASTA file
#' @param out.dir Path to output directory
#' @param program Which mummer program to use (nucmer or promer)
#' @param breaklen Distance an alignment extension will attempt to extend poor scoring regions before giving up
#' @param mincluster Minimum length of a cluster of matches
#' @param minmatch Minimum length of a single match
#' @param maxgap Maximum gap between two adjacent matches in a cluster
#' @param prefix Character prefix to temporary files
#' @param threads Number of threads
#' @param recompute Wether to delete the out.dir and recompute
#'
#' @export
#'
mummer_alignment <- function(query, reference, out.dir, 
                             program = 'nucmer',
                             breaklen = if (program == 'nucmer') 200 else 60,
                             mincluster = if (program == 'nucmer') 65 else 20,
                             minmatch = if (program == 'nucmer') 20 else 6,
                             maxgap = if (program == 'nucmer') 90 else 30,
                             prefix = program,
                             recompute = FALSE
                            ) {

    # Minimal check
    stopifnot(
        is_file(reference),
        is_file(query)
    )
    check_installed("mummer", silent = TRUE)
    if (!program %in% c('nucmer','promer')) stop('Program must be nucmer or promer.')

    # Variables
    out.dir <- if (endsWith(out.dir, '/')) out.dir else paste0(out.dir, '/')
    if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
    wd <- getwd()
    query.abspath <- paste0(wd,'/',query)
    ref.abspath <- paste0(wd,'/',reference)
    temp.delta <- paste0(prefix,'.delta')
    out.delta <- paste0(out.dir,'nucmer.delta')
    log.file <- paste0(out.dir,'nucmer.log')
    out.coords <- paste0(out.dir,'nucmer.coords')

    # Timestamp
    t0 <- Sys.time()

    # Input
    ref.seqs <- readDNAStringSet(ref.abspath)
    ref.names <- str_split(names(ref.seqs), ' ', simplify=TRUE)[,1]
    qry.seqs <- readDNAStringSet(query.abspath)
    qry.names <- str_split(names(qry.seqs), ' ', simplify=TRUE)[,1]

    # Output
    if (file.exists('temp.delta')) warning('Temporary delta file present. Indicates problems, will be overwritten & removed...')
    if (recompute) {
        unlink(out.delta)
        unlink(log.file)
        unlink(out.coords)
    }
    if (file.exists(out.delta)) {
        msg <- paste('Output files found in ',out.dir,'...')
        message(msg)
        delta <- readLines(out.delta, n=1)
        delta <- str_split(delta, ' ')[[1]]
        check <- delta[1] == ref.abspath & delta[2] == query.abspath
        if (check) {
            # Read
            coords <- read_mummer_coords(out.coords)
            # Modify
            coords$query <- query
            coords$query_name <- names(query)
            coords$ref <- reference
            coords$ref_name <- names(reference)
            coords$ref_id <- factor(coords$ref_id, ref.names)
            # Exit
            return(coords)
        }
    }

    # nucmer
    cmd <- paste(program,reference,query,'-p',prefix,'-b',breaklen,'-c',mincluster,'-l',minmatch,'-g',maxgap)
    system3(cmd, log.file = log.file)
    file.copy(temp.delta,out.delta, overwrite = TRUE)
    file.remove(temp.delta)

    # show-coords
    cmd <- paste('show-coords','-H',out.delta,'>',out.coords)
    system3(cmd, include.errors = FALSE)

    # Read
    coords <- read_mummer_coords(out.coords)

    # Modify
    coords$query <- query
    coords$query_name <- names(query)
    coords$ref <- reference
    coords$ref_name <- names(reference)
    coords$ref_id <- factor(coords$ref_id, ref.names)

    # Timestamp
    t1 <- Sys.time()
    print(t1 - t0)

    # Exit
    return(coords)
}

#' Plot MUMmer coordinates
#'
#' @param x Data.frame with MUMmer coordinates
#'
#' @export
#'
plot_mummer_alignment <- function(x, max.contigs = 10) {

    # Minimal check
    stopifnot(
        is.data.frame(x),
        c('ref_start','ref_end','query_start','query_end','identity','mum','ref_id','query_id') %in% names(x)
    )

    # Faceting
    facet_x <- if (length(unique(x$ref_id)) <= max.contigs) TRUE else FALSE
    facet_y <- if (length(unique(x$query_id)) <= max.contigs) TRUE else FALSE
    facet <- if (facet_x & facet_y) facet_grid(query_id ~ ref_id, scales = 'free') else
            if (facet_x & !facet_y) facet_wrap(. ~ ref_id, scales = 'free') else
            if (!facet_x & facet_y) facet_wrap(query_id ~ ., scales = 'free') else
            if (!facet_x & !facet_y) NULL
            
    
    
    # Labels
    xname <- unique(x[['ref_name']])
    yname <- unique(x[['query_name']])
    
    # Plot
    plot <- ggplot(x) +
      geom_segment(aes(x = ref_start, xend = ref_end, y = query_start, yend = query_end, col = identity, group = mum), linewidth = 1) +
      facet +
      scale_color_distiller(palette = 'RdYlBu') +
      scale_x_continuous(limits = c(0,NA)) +
      scale_y_continuous(limits = c(0,NA)) +
      theme_light(15) +
      theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text.y = element_text(angle = 0)
      ) +
      labs(x = xname, y = yname)

    return(plot)
}

#' Read MUMmer coordinates
#'
#' @param x MUMmer coordinate file
#' 
read_mummer_coords <- function(file) {

    # Minimal check
    stopifnot(
        is_file(file)
    )

    # Variables
    coords_columns <- c(
        'sep_0',
        'ref_start',
        'ref_end',
        'sep_1',
        'query_start',
        'query_end',
        'sep_2',
        'ref_length',
        'query_length',
        'sep_3',
        'identity',
        'sep_4',
        'ref_id',
        'query_id'
    )
    
    # Read
    x <- readLines(file)
    x <- str_split(x, pattern = '[:space:]+', simplify=TRUE)
    x <- as.data.frame(x)

    # Name
    names(x) <- coords_columns

    # Clean
    ind <- which(!str_detect(names(x), 'sep'))
    x <- x[, ind]

    # Modify
    x$mum <- 1:nrow(x)
    for (i in c('identity','ref_start','ref_end','query_start','query_end','ref_length','query_length')) {
        x[[i]] <- as.numeric(x[[i]])
    }
    
    # Exit
    return(x)
}