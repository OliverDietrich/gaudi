# Wrappers for the NCBI datasets CLI

#' Summary of genomes present on NCBI
#' 
#' Download assembly summary for GenBank and RefSeq from NCBI FTP.
#' 
#' @param dirname Name of directory to save assembly summaries (default: data/)
#' @param refseq Name of refseq assembly summary
#' @param genbank Name of genbank assembly summary
#' @param url_refseq URL to NCBI FTP for refseq summary
#' @param url_genbank URL to NCBI FTP for genbank summary
#' @param max.download.minutes Numeric, maximum download time in minutes
#' @param max.download.minutes Numeric, maximum download time in seconds
#' 
#' @export
ncbi_genome_summary <- function(dirname='data/',
                                summary = 'data/ncbi_genome_summary.tsv',
                                refseq='refseq_assembly_summary.tsv',
                                genbank='genbank_assembly_summary.tsv',
                                url_refseq='ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt',
                                url_genbank='ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt',
                                max.download.minutes = 30,
                                max.download.seconds = max.download.minutes*60
) {

    # Checks
    if(!dir.exists(dirname)) {
        msg <- paste('Directory',dirname,'does not exist.')
        stop(msg)
    }
    if (!endsWith(dirname,'/')) {
        dirname <- paste0(dirname,'/')
    }
    stopifnot(
        is.numeric(max.download.seconds)
    )

    # Exit 1
    if (file.exists(summary)) {
        object <- vroom::vroom(summary, show_col_types=FALSE)
        return(object)
    }

    # Modify file names
    genbank <- paste0(dirname,genbank)
    refseq <- paste0(dirname,refseq)
    
    # Download files
    options(timeout = max.download.seconds)
    if (!file.exists(genbank)) {
        msg <- paste('Downloading',url_genbank,'to',genbank)
        message(msg)
        download.file(url_genbank, genbank)
    }
    if (!file.exists(refseq)) {
        msg <- paste('Downloading',url_refseq,'to',refseq)
        message(msg)
        download.file(url_refseq, refseq)
    }

    # Read summaries
    object <- list()
    
    ## Refseq
    msg <- paste('Reading',refseq)
    message(msg)
    object$refseq <- vroom::vroom(refseq, skip=1, show_col_types = FALSE)

    ## Genbank
    msg <- paste('Reading',genbank)
    message(msg)
    object$genbank <- vroom::vroom(genbank, skip=1, show_col_types=FALSE)

    # Add source, harmonize vector classes
    for (i in names(object)) {
        object[[i]]$source <- i
        for (j in c('total_gene_count','protein_coding_gene_count','non_coding_gene_count')) {
            object[[i]][[j]] <- suppressWarnings(as.numeric(object[[i]][[j]]))
        }
    }

    # Combine
    object <- dplyr::bind_rows(object)

    # Re-name columns
    object$Assembly.version <- object$`#assembly_accession`

    # Write cleaned summary
    vroom::vroom_write(object, summary)
    
    # Return object & summary
    print(str(object, max.level = 0))
    return(object)
}

#' Download genome sequence data from NCBI using Assembly.version IDs
#'
#' Wrapper for the NCBI datasets CLI
#' Based on the call 'datasets download genome accession' it retrieves sequence data
#' for Assembly.version numbers. Checks determine completeness of the downloaded data.
#' 
#' @param accession Character vector with Assembly.version numbers
#' @param inputfile File with accession numbers (one per line)
#' 
#' @export
ncbi_datasets_download_genome <- function(Assembly.version = NULL,
                                          out.dir = NULL,
                                          debug = FALSE,
                                          debug.n = 10
                                         ) {

    # Minimal check
    stopifnot(
        !is.null(Assembly.version),
        !is.null(out.dir)
    )
    check_installed('datasets', silent=TRUE)
    check_version('datasets')

    # Create file names/paths
    download.log <- paste0(out.dir,'download.log')
    out.acc <- paste0(out.dir,'all_files.txt')
    out.missing <- paste0(out.dir,'missing_files.txt')
    out.zip <- paste0(out.dir,'ncbi_dataset.zip')
    out.data <- paste0(out.dir,'ncbi_dataset/data/')
    in.checksums <- paste0(out.dir,'md5sum.txt')
    out.final <- paste0(out.dir,'assemblies/')
    out.checksums <- paste0(out.final,'md5sum.tsv')
    report.jsonl <- paste0(out.final,'assembly_data_report.jsonl')
    report.tsv <- paste0(out.final,'assembly_data_report.tsv')
    categories_genome <- c('genome','rna','protein','cds','gff3','gtf','gbff','seq-report') # Sequence types to include
    
    # Check input
    if (length(Assembly.version) == 1 & all(file.exists(Assembly.version))) {
        msg <- paste('Input via file:', Assembly.version)
        message(msg)
        out.acc <- Assembly.version
        Assembly.version <- readLines(out.acc)
    } else {
        writeLines(Assembly.version, out.acc)
    }

    dups <- duplicated(Assembly.version)
    if (any(dups)) {
        Assembly.version <- Assembly.version[!duplicated(Assembly.version)]
    }

    # Check output
    if (!endsWith(out.dir,'/')) {
        out.dir <- paste0(out.dir,'/')
    }
    output_present <- list.files(out.final)
    index <- Assembly.version %in% output_present
    missing <- Assembly.version[!index]
    missing <- if (debug) head(missing, debug.n) else missing # DEBUGGING
    msg <- paste('Output present for',sum(index),'out of',length(index),'genomes.','Downloading', length(missing),'...')
    message(msg)
    writeLines(missing, out.missing)

    # Download dataset
    if (!file.exists(out.zip)) {
        cmd <- paste('datasets download genome accession','--inputfile',out.missing,'--include',paste(categories_genome,collapse=','),'--filename',out.zip,'2>&1')
        message(cmd)
        stdout <- system(cmd, intern=TRUE)
        stdout <- paste(stdout,'\n')
        writeLines(stdout, download.log)
    } else {
        msg <- paste('NCBI archive', out.zip, 'exists. Skipping download...')
        message(msg)
    }

    # Check successful download
    if (!file.exists(out.zip)) {
        msg <- paste('NCBI archive', out.zip, 'was not retrieved successfully. Aborting.')
        stop(msg)
    }

    # Extract archive
    if (!file.exists(out.data)) {
        msg <- paste('Extracting data from',out.zip)
        message(msg)
        unzip(out.zip, exdir =out.dir)
    } else {
        msg <- paste('NCBI archive', out.data, 'exists. Skipping unzip...')
        message(msg)
    }

    # Check archive extraction
    if (!file.exists(out.data)) {
        msg <- paste('NCBI archive',out.zip,'could not be extracted successfully. Removing archive...')
        file.remove(out.zip)
        stop(msg)
    }

    # Check md5sums and transfer files
    if (file.exists(in.checksums)) {
        checksums <- read.table(in.checksums, col.names = c('ncbi_checksum','file'))
        checksums$file <- paste0(out.dir, checksums$file)        
        checksums$local_checksum <- sapply(checksums$file, md5sums)
        checksums$newfile <- paste0(out.final, str_remove(checksums$file, out.data))
        index <- checksums$local_checksum == checksums$ncbi_checksum
        msg <- paste('Checksums for',sum(index),'of',length(index),'files match. Moving files to keep...')
        checksums <- checksums[index, ]
        dirs_to_create <- unique(dirname(checksums$newfile))
        sapply(dirs_to_create, dir.create)
        Map(file.rename, checksums$file, checksums$newfile)
    
        # Clean working directory    
        file.remove(in.checksums)
        unlink(dirname(out.data), recursive = TRUE)
        unlink(out.zip)
    }
    
    # Re-format files
    if (file.exists(report.jsonl) & !file.exists(report.tsv)) {
        cmd <- paste('dataformat tsv genome','--inputfile',report.jsonl,'>',report.tsv,'2>&1')
        message(cmd)
        system(cmd)
    }
    
    ## Check completeness
    index <- Assembly.version %in% list.files(out.final)
    tbl <- table('Genomes present:'=index)
    
    return(tbl)
}

#' Summary of genomes present on NCBI
#' 
#' Wrapper for NCBI datasets CLI based on dataset summary (virus) genome.
#' Downloads genome reports and returns formatted output.
#' 
#' @param accession Character vector or file name containing accession.version or assembly.version numbers
#' @param out.tsv Name of TSV file (temporary)
#' @param remove.intermediates Boolean, whether to remove TSV file
#'
#' @export
ncbi_datasets_summary_genome <- function(accession = NULL,
                                         out.acc = 'datasets_summary.ids',
                                         out.tsv = 'datasets_summary.tsv',
                                         remove.intermediates = TRUE
                                        ) {

    stopifnot(
        !is.null(accession)
    )

    # Check programs
    check_installed(c('datasets','dataformat'), silent=TRUE)
    check_version('datasets')
    
    # Check input
    # ...

    # Write file
    writeLines(accession, out.acc)
    
    # Get dataset summary
    cmd <- paste('datasets summary genome accession --inputfile',out.acc,'--as-json-lines | dataformat tsv genome','>',out.tsv,'2>&1')
    system(cmd, intern=TRUE)

    # Create output object
    object <- readr::read_tsv(out.tsv)

    # Remove intermediates
    if (remove.intermediates) {
        file.remove(out.tsv)
        file.remove(out.acc)
    }

    # Exit 0
    return(object)
}

#' List NCBI genome dataset files
#'
#' List files of NCBI genomes by category (CDS, genome, gbff, gff, gtf, proteins, Other)
#' 
#' @param path Path to genome dataset (e.g. .../ncbi_dataset/data/GCA_000000000.1/)
#' 
#' @export
ncbi_genome_files <- function(path = NULL) {

    # Minimal check
    stopifnot(
        !is.null(path)
    )

    # Check input
    if (!endsWith(path,'/')) {
        path <- paste0(path,'/')
    }

    # List and categorize files
    object <- list.files(path)
    object <- data.frame(
        'file' = object,
        'type' = dplyr::case_when(
            object == 'cds_from_genomic.fna' ~ "CDS",
            endsWith(object, '_genomic.fna') ~ "genome",
            endsWith(object, '.gbff') ~ "gbff",
            endsWith(object, 'gff') ~ "gff",
            endsWith(object, '.gtf') ~ "gtf",
            object == 'protein.faa' ~ "proteins",
            .default = 'Other'
        )
    )
    object$file <- paste0(path,object$file)

    # Pivot wider
    object <- tidyr::pivot_wider(object, names_from=type, values_from=file, values_fn = ~paste(.x, collapse=','))
    
    # Exit 0
    return(object)
}