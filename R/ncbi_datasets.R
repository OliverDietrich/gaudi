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
                                refseq='refseq_assembly_summary.txt',
                                genbank='genbank_assembly_summary.txt',
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
    object$refseq <- readr::read_tsv(refseq, skip=1, show_col_types = FALSE)

    ## Genbank
    msg <- paste('Reading',genbank)
    message(msg)
    object$genbank <- readr::read_tsv(genbank, skip=1, show_col_types=FALSE)

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

    # Add Accession.version
    # ... 
    # id_file <- tempfile()
    # ids <- object$Assembly.version
    # writeLines(ids,id_file)
    #  datasets summary genome accession GCF_000001405.40 --report sequence --as-json-lines | dataformat tsv genome-seq

    # Return object & summary
    print(str(object, max.level = 0))
    return(object)
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
                                          overwrite=FALSE
                                         ) {

    # Minimal check
    stopifnot(
        !is.null(Assembly.version),
        !is.null(out.dir)
    )
    check_installed('datasets', silent=TRUE)
    check_version('datasets')

    # Create file names/paths
    out.acc <- paste0(out.dir,'AssemblyVersion.txt')
    out.zip <- paste0(out.dir,'ncbi_dataset.zip')
    out.data <- paste0(out.dir,'ncbi_dataset/data/')
    report.jsonl <- paste0(out.data,'assembly_data_report.jsonl')
    report.tsv <- paste0(out.data,'assembly_data_report.tsv')
    
    # Check input
    if (length(Assembly.version) == 1 & all(file.exists(Assembly.version))) {
        msg <- paste('Input via file:', Assembly.version)
        message(msg)
        out.acc <- Assembly.version
    } else {
        writeLines(Assembly.version, out.acc)
    }

    # Check output
    if (!endsWith(out.dir,'/')) {
        out.dir <- paste0(out.dir,'/')
    }
    if (dir.exists(out.data) & !overwrite) {
        msg <- paste('Output directory',out.data,'exists.')
        if (overwrite) {
            msg <- paste(msg, 'Overwriting.')
            warning(msg)
        } else {
            msg <- paste(msg, 'Aborting.')
            warning(msg)
            return('Exit 1: Output present.')
        }
    }

    ## Sequence types to include
    #categories_virus_genome <- c('genome','cds','protein','annotation','biosample') # ,'none'
    categories_genome <- c('genome','rna','protein','cds','gff3','gtf','gbff','seq-report') # ,'none'

    # Download dataset
    cmd <- paste('datasets download genome accession',paste(Assembly.version,collapse=' '),'--include',paste(categories_genome,collapse=','),'--filename',out.zip,'2>&1')
    cmd <- paste('datasets download genome accession','--inputfile',out.acc,'--include',paste(categories_genome,collapse=','),'--filename',out.zip,'2>&1')
    message(cmd)
    system(cmd, intern=TRUE)

    # Extract archive
    msg <- paste('Extracting data from',out.zip)
    message(msg)
    unzip(out.zip, exdir=out.dir)

    # Re-format files
    cmd <- paste('dataformat tsv genome','--inputfile',report.jsonl,'>',report.tsv,'2>&1')
    message(cmd)
    system(cmd, intern=TRUE)
}