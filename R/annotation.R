#' Detect MGEs with geNomad
#'
#' CLI wrapper for geNomad to detect mobile genetic elements (MGE)
#'
#' @param input.fasta Character, path to input FASTA
#' @param output.dir Character, path to output directory
#' @param db Character, path to database
#' @param force Boolean, whether to overwrite output.dir
#'
#' @export
geNomad <- function(input.fasta, output.dir, db, force=FALSE) {

    # Minimal check
    stopifnot(
        is_file(input.fasta, suffix=c('.fasta','.fna')),
        dir.exists(db)
    )
    check_installed('genomad', silent=TRUE)
    check_version('genomad')

    # Input
    output.dir <- if(endsWith(output.dir,'/')) output.dir else paste0(output.dir,'/')

    # Variables
    log.file <- paste0(output.dir,'runtime.log')
    classification.tsv <- paste0(output.dir,'genome_aggregated_classification/genome_aggregated_classification.tsv')
    prophage.tsv <- paste0(output.dir,'genome_find_proviruses/genome_provirus.tsv')
    out.gff <- paste0(output.dir,'regions.gff')

    # Output
    if (force) unlink(output.dir, recursive=TRUE)
    if (file.exists(out.gff)) {
        msg <- paste('Output file',out.gff,'already exists. Aborting...')
        warning(msg)
        return()
    }

    # DB
    if (endsWith(db,'/')) {
        db <- stringr::str_sub(db, start = 0, end = -2)
    }
    if (!endsWith(db,'genomad_db')) {
        msg <- 'geNomad DB must be named "genomad_db"!'
        stop(msg)
    }
    if (!dir.exists(db)) {
        msg <- paste('geNomad database not present. Downloading to',db)
        warning(msg)
        cmd <- paste('genomad download-database',dirname(db))
        system3(cmd)
    }

    # Timestamp
    t.start <- Sys.time()

    # Main
    cmd <- paste('genomad end-to-end',input.fasta,output.dir,db)
    if (!file.exists(classification.tsv)) {
        system3(cmd, log.file = log.file)
    }

    # Read output files
    class <- readr::read_tsv(classification.tsv, show_col_types = FALSE)
    proph <- readr::read_tsv(prophage.tsv, show_col_types = FALSE)
    seqs <- Biostrings::readDNAStringSet(input.fasta)

    # Checkpoint
    # I though about checking if fasta headers match the seq_name in the genomad classification TSV
    # However, I would have to incorporate this check in RunGenomad to be effective...
    # Currently, the best option is to wait for an error (caused by assigning width(seqs) to gff$end) 
    # and continue by setting 'recompute = TRUE' in RunGenomad().

    # Assign contig types as 'Name'
    assign <- pivot_longer(class, cols = c('chromosome_score','plasmid_score','virus_score'), names_to = 'tag', values_to = 'value')
    assign <- as.data.frame(assign)
    assign$label <- str_split(assign$tag, '_', simplify=TRUE)[,1]
    assign <- group_by(assign, seq_name)
    assign <- mutate(assign, max = max(value))
    assign <- subset(assign, value == max)
    ind <- match(assign$seq_name, class$seq_name)
    class$Name <- assign$label[ind]
    class$score <- assign$value[ind]

    # Format as GFF3
    gff <- class
    names(gff)[1] <- 'seqid'
    gff$source <- 'geNomad'
    gff$type <- 'region'
    gff$start <- 1
    gff$end <- width(seqs)
    gff$strand <- '+'
    gff$phase <- '.'
    gff <- format_gff3(gff)

    # Keep seqid order
    lvls <- as.character(unique(gff$seqid))

    # Add prophage regions
    if (nrow(proph) > 0) {
        proph <- data.frame(
            'seqid' = proph[['source_seq']],
            'source' = 'geNomad',
            'type' = 'region',
            'start' = proph[['start']],
            'end' = proph[['end']],
            'score' = proph[['v_vs_c_score']],
            'strand' = '+',
            'phase' = '.',
            'attributes' = 'Name=prophage'
        )
        proph <- format_gff3(proph)
    } else {
        proph <- NULL
    }

    # Combine
    gff <- rbind(gff, proph)

    # Order entries
    gff$seqid <- factor(gff$seqid, lvls)
    gff <- gff[order(gff$seqid, gff$start), ]

    # FASTA
    attr(gff, 'sequence') <- seqs

    # Write GFF
    write_gff3(gff, out.gff)
    
    # Timestamp
    t.stop <- Sys.time()
    print(t.stop - t.start)
}

#' Run GeNomad for a genomeCollection
#'
#' @param object A genomeCollection object
#'
#' @export
#'
RunGenomad <- function(object, 
                       recompute = FALSE,
                       recompute.sample = NULL,
                       debug = FALSE, debug.n = 1
                      ) {

    # Minimal check
    stopifnot(
        is(object) == 'genomeCollection'
    )

    # Variables
    ind <- file.exists(object$genome)
    out.dir <- paste0(object$path[ind],'annotation/genomad/')
    genomes <- object$genome[ind]

    # Output
    data <- Annotation(object, 'geNomad')
    if (is.null(data)) {
        data <- methods::new("Annotation", 
                             index = index(object)[ind],
                             genes = paste0(out.dir,'regions.gff'),
                             cds = character(),
                             proteins = character(),
                             log = paste0(out.dir,'runtime.log'),
                             genome = 'genome',
                             tool = 'geNomad',
                             type = "default"
                            )
    }

    # Input
    ind <- file.exists(data$genes)
    missing <- data$index[!ind]
    missing <- if (recompute) data$index else missing
    missing <- if (debug) head(missing, debug.n) else missing
    ind.rc <- recompute.sample %in% data$index
    missing <- if (any(ind.rc)) recompute.sample[ind.rc] else missing
    msg <- paste('Found genes for',sum(ind),'out of',length(data$index),'genomes.','Running geNomad for',length(missing),'...')
    message(msg)

    # Main
    for (sample in missing) {
        i <- which(data$index == sample)
        dir.create(out.dir[[i]], recursive = TRUE)
        geNomad(input.fasta = genomes[[i]], output.dir = out.dir[[i]], db = '../databases/genomad_db/', force = TRUE)
    }

    # Assign object
    Annotation(object, 'geNomad') <- data

    # Exit
    return(object)
}

#' Bacterial ORF detection by Prodigal
#'
#' CLI wrapper for open reading frame (ORF) detection by Prodigal
#'
#' @param genome Path to genome FASTA file
#' @param output.dir Path to output directory
#' @param output.format Select output format (gbk, gff, sco)
#' @param procedure Select procedure (single or meta)
#' @param translation.table Specify a translation table to use (1-25, default: auto)
#' @param write.protein.translations Whether to write protein translations
#' @param closed.ends Whether genes can run off the edges of a contig
#' @param write.cds Whether to write coding nucleotide sequences (CDS)
#' @param mask.N Whether to treat runs of N as masked sequence; don't build genes across them
#' @param write.potential.genes Whether to write all potential genes (with scores) to the selected file
#'
#' @export
prodigal <- function(genome, output.dir, output.format = 'gff',
                     procedure = 'single',
                     translation.table = 'auto',
                     write.protein.translations = TRUE, 
                     closed.ends = FALSE,
                     write.cds = TRUE,
                     mask.N = TRUE,
                     write.potential.genes = FALSE,
                     recompute = FALSE
                    ) {

    # Minimal check
    stopifnot(
        length(genome) == 1,
        file.exists(genome),
        length(output.dir) == 1
    )
    check_installed('prodigal', silent=TRUE)
    check_version('prodigal', command = '-v')

    # Variables
    if (!endsWith(output.dir,'/')) {
        output.dir <- paste0(output.dir,'/')
    }
    out.file <- paste0(output.dir,'genes','.',output.format)
    out.proteins <- paste0(output.dir,'proteins.faa')
    out.cds <- paste0(output.dir,'cds.fna')
    out.log <- paste0(output.dir,'runtime.log')
    out.potential <- paste0(output.dir,'potential-genes.gff')

    # Output
    if (file.exists(out.file) & !recompute) {
        msg <- paste('Output file', out.file, 'already exists')
        warning(msg)
        return()
    }
    dir.create(output.dir, recursive = TRUE)

    # Timestamp
    t.start <- Sys.time()

    # Main
    cmd <- paste('prodigal','-i',genome,'-f',output.format,'-o',out.file,'-p',procedure)
    cmd <- if (translation.table %in% 1:25) paste(cmd,'-g',translation.table) else cmd
    cmd <- if(write.protein.translations) paste(cmd,'-a',out.proteins) else cmd
    cmd <- if(closed.ends) paste(cmd,'-c') else cmd
    cmd <- if(write.cds) paste(cmd,'-d',out.cds) else cmd
    cmd <- if(mask.N) paste(cmd,'-m') else cmd
    cmd <- if(write.potential.genes) paste(cmd,'-s',out.potential) else cmd
    system3(cmd, log.file = out.log)

    # Timestamp
    t.stop <- Sys.time()
    print(t.stop - t.start)
}

#' Run Prodigal for a genomeCollection
#'
#' @param object A genomeCollection object
#' @param recompute Whether to re-run Prodigal even though an output GFF exists
#' @param debug Whether to run in debug mode (run only for debug.n samples)
#' @param debug.n Number of samples to run while debugging
#'
#' @export
#'
RunProdigal <- function(object, 
                        recompute = FALSE,
                        debug = FALSE, debug.n = 1
                       ) {

    # Minimal check
    stopifnot(
        is(object) == 'genomeCollection'
    )

    # Variables
    ind <- file.exists(object$genome)
    out.dir <- paste0(object$path[ind],'annotation/prodigal/')
    genomes <- object$genome[ind]

    # Output
    data <- Annotation(object, 'Prodigal')
    if (is.null(data)) {
        data <- methods::new("Annotation", 
                             index = index(object)[ind],
                             genes = paste0(out.dir,'genes.gff'),
                             cds = paste0(out.dir,'cds.fna'),
                             proteins = paste0(out.dir,'proteins.faa'),
                             log = paste0(out.dir,'runtime.log'),
                             genome = 'genome',
                             tool = 'Prodigal',
                             type = "default"
                            )
    }

    # Input
    ind <- file.exists(data$genes)
    missing <- data$index[!ind]
    missing <- if (recompute) data$index else missing
    missing <- if (debug) head(missing, debug.n) else missing
    msg <- paste('Found genes for',sum(ind),'out of',length(data$index),'genomes.','Running Prodigal for',length(missing),'...')
    message(msg)

    # Main
    for (sample in missing) {
        i <- which(data$index == sample)
        prodigal(genome = genomes[[i]], output.dir = out.dir[[i]], recompute = TRUE)
    }

    # Assign object
    Annotation(object, 'Prodigal') <- data

    # Exit
    return(object)
}


#' PHANOTATE: A phage genome annotator
#'
#' CLI wrapper for open reading frame (ORF) detection by PHANOTATE
#'
#' @param in.file Input file in FASTA format
#' @param out.dir Path to output directory
#' @param format Output the features in the specified format (genbank, gff, gff3, tabular, fasta, fna, faa)
#' @param start.codons Comma separated list of start codons and frequency [atg:0.85,gtg:0.10,ttg:0.05]
#' @param stop.codons Comma separated list of stop codons [tag,tga,taa]
#' @param min.orf.length Minimum ORF length
#'
#' @export
phanotate <- function(in.file, out.dir, 
                      format = 'gff', 
                      start.codons = c('atg'=0.85, 'gtg'=0.10, 'ttg'=0.05),
                      stop.codons = c('tag','tga','taa'),
                      min.ORF.length = NULL,
                      recompute = FALSE
                     ) {

    # Minimal check
    stopifnot(
        is_file(in.file, suffix = c('fasta','fna'))
    )
    check_installed('phanotate.py', silent=TRUE)
    version <- check_version('phanotate.py', return.version = TRUE)

    # Variables
    out.dir <- if (endsWith(out.dir,'/')) out.dir else paste0(out.dir,'/')
    out.file <- paste0(out.dir,'genes','.',format)
    log.file <- paste0(out.dir,'runtime.log')
    start.codons <- vector_to_python_list(start.codons)
    stop.codons <- vector_to_python_list(stop.codons)

    # Timestamp
    t.start <- Sys.time()

    # Output
    if (recompute) unlink(out.dir, recursive=TRUE)
    if (file.exists(out.file)) {
        msg <- paste('Output file',out.file,'already exists. Skipping PHANOTATE...')
        warning(msg)
        return()
    }
    dir.create(out.dir, recursive=TRUE)

    # Main
    cmd <- paste('phanotate.py',in.file,'-o',out.file,'-f',format) # ,'-s',start.codons,'-e',stop.codons
    cmd <- if (length(min.ORF.length)) paste(cmd,'-l',min.ORF.length) else cmd
    system3(cmd, log.file = log.file)

    # Fix GFF3 formatting
    if (format %in% c('gff', 'gff3')) {
        gff <- read_gff3(out.file)
        gff$source <- paste0('PHANOTATE_v',version)
        gff$score <- as.numeric(str_split(gff$note, ':', simplify=TRUE)[,2])
        write_gff3(gff, out.file)
    }

    # Timestamp
    t.stop <- Sys.time()
    print(t.stop - t.start)
}

#' Run PHANOTATE for a genomeCollection
#'
#' @param object A genomeCollection object
#' @param recompute Whether to re-run PHANOTATE even though an output GFF exists
#' @param debug Whether to run in debug mode (run only for debug.n samples)
#' @param debug.n Number of samples to run while debugging
#'
#' @export
#'
RunPHANOTATE <- function(object, 
                        recompute = FALSE,
                        recompute.sample = NULL,
                        debug = FALSE, debug.n = 1
                        ) {

    # Minimal check
    stopifnot(
        is(object) == 'genomeCollection'
    )

    # Variables
    ind <- file.exists(object$genome)
    out.dir <- paste0(object$path[ind],'annotation/phanotate/')
    genomes <- object$genome[ind]

    # Output
    data <- Annotation(object, 'PHANOTATE')
    if (is.null(data)) {
        data <- methods::new("Annotation", 
                             index = index(object)[ind],
                             genes = paste0(out.dir,'genes.gff'),
                             cds = character(),
                             proteins = character(),
                             log = paste0(out.dir,'runtime.log'),
                             genome = 'genome',
                             tool = 'PHANOTATE',
                             type = "default"
                            )
    }

    # Input
    ind <- file.exists(data$genes)
    missing <- data$index[!ind]
    missing <- if (recompute) data$index else missing
    missing <- if (debug) head(missing, debug.n) else missing
    ind.rc <- recompute.sample %in% data$index
    missing <- if (any(ind.rc)) recompute.sample[ind.rc] else missing
    msg <- paste('Found genes for',sum(ind),'out of',length(data$index),'genomes.','Running PHANOTATE for',length(missing),'...')
    message(msg)

    # Main
    for (sample in missing) {
        i <- which(data$index == sample)
        phanotate(in.file = genomes[[i]], out.dir = out.dir[[i]], recompute = TRUE)
    }
    
    # Assign object
    Annotation(object, 'PHANOTATE') <- data
    
    # Exit
    return(object)
}

#' Combine annotations for a genomeCollection
#'
#' Combine the genes from multiple Annotation objects in a genomeCollection.
#' 
#' @param object A genomeCollection object
#' @param name Name of the slot in metadata(object)
#' @param annotations A character vector of annotations present in Annotations(object)
#' @param overwrite Whether to overwrite existing output files
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#'
CombineAnnotations <- function(object, name, annotations, overwrite = FALSE) {

    # Minimal check
    stopifnot(
        is(object) == 'genomeCollection'
    )
    if (is.null(object[[name]])) {
        msg <- paste('No slot', name, 'found in object. Will be added to metadata(object).')
    } else
    if (name %in% names(object[[]])) {
        msg <- paste0('Column "', name, '" already exists in metadata(object).')
        message(msg)
    } else {
        msg <- paste0('Slot "', name, '" already exists in object. Aborting...')
        stop(msg)
    }

    # Variables
    object[[name]] <- paste0(object$path, name,'.gff3')
    ind.genome <- file.exists(object$genome)
    ind.regions <- file.exists(object$regions)
    ind.anns <- annotations %in% Annotations(object)
    missing <- if (overwrite) index(object)[ind.genome] else index(object)[ind.genome & !ind.regions]

    # Input
    if (length(ind)) {
        msg <- paste('Found', sum(ind.anns), 'out of', length(ind.anns), 'Annotation objects:', paste(annotations[ind.anns], collapse=', '))
        message(msg)
    } else {
        msg <- paste('Did not found any Annotation objects in genomeCollection. Please check what has been computed!')
        stop(msg)
    }
    
    # Output
    msg <- paste('Found', sum(ind.regions), 'out of', sum(ind.genome), 'output files. Combining annotations for', length(missing), 'genomes.')
    message(msg)

    # Main
    for (sample in missing) {
        index <- which(index(object) == sample)
        
        # Get GFF path
        x <- list()
        for (i in annotations) {
            data <- object[[i]]
            ind <- data$index == sample
            x[[i]] <- object[[i]]$genes[ind]   
        }

        # Read GFF
        for (i in names(x)) {
            x[[i]] <- suppressMessages(read_gff3(x[[i]], extract.attributes = FALSE))
            
        }
        x <- dplyr::bind_rows(x)

        # Set sequence attribute
        attr(x, 'sequence') <- readDNAStringSet(object$genome[[index]])

        # Write GFF
        write_gff3(x, object[[name]][[index]])
    }

    # Exit
    return(object)
}

#' Bakta annotation
#' 
#' CLI wrapper for genome annotation using Bakta
#'
#' @param input.fasta Path to input FASTA file
#' @param output.dir Path to output directory
#' @param output.prefix Prefix of file names
#' @param regions Path to GFF3 file with pre-computed CDS regions
#' @param db Character, path to the bakta database
#' @param download.db Boolean, whether to download the database if missing
#' @param threads Integer, number of threads to use (default: all)
#' @param force Boolean, whether to force overwriting existing output folder
#'
#' @importFrom stringr str_detect
#'
#' @export
bakta <- function(input.fasta, output.dir,
                  db='../databases/bakta',
                  prefix='bakta',
                  keep.contig.headers = TRUE,
                  regions=NULL,
                  download.db=FALSE, 
                  threads=n_proc(), 
                  force=FALSE
                 ) {

    # Minimal check
    stopifnot(
        is_file(input.fasta, suffix=c('.fasta','.fna')),
        dir.exists(db)
    )
    check_installed('bakta', silent=TRUE)
    check_version('bakta')

    # Variables
    db <- if (endsWith(db,'/')) db else paste0(db,'/')
    output.dir <- if (endsWith(output.dir,'/')) output.dir else paste0(output.dir,'/')
    download.log <- paste0(db,'download.log')
    bakta.db <- paste0(db,'db/')
    bakta.db.version <- paste0(bakta.db,'/version.json')
    int.regions <- paste0(output.dir,'regions.gff')
    out.gff <- paste0(output.dir,prefix,'.gff3')
    log.file <- paste0(output.dir,'runtime.log')

    # Check output
    if (force) unlink(output.dir, recursive=TRUE)
    if (file.exists(out.gff)) {
        msg <- paste('Output file',out.gff,'already exists. Aborting...')
        warning(msg)
        return()
    }
    dir.create(output.dir)

    # Check DB
    if (is.null(db)) {
        msg <- 'Please supply a path for the Bakta database!'
        stop(msg)
    } else if (download.db) {
        msg <- paste('Downloading bakta database to',bakta.db)
        message(msg)
        cmd <- paste('bakta_db download','--output',db,'--type','full')
        system3(cmd, log.file=download.log)
    } else if (file.exists(bakta.db.version)) {
        bakta.db.version <- jsonlite::read_json(bakta.db.version)
        msg <- paste0('Bakta database (',bakta.db.version$type,'), version ',bakta.db.version$major,'.',bakta.db.version$minor,' (',bakta.db.version$date,')')
        message(msg)
    } else {
        msg <- paste('Bakta database',db,'does not exists.','Consider passing download.db=TRUE.')
        stop(msg)
    }

    # Check log for errors
    if (file.exists(log.file)) {
        log <- readLines(log.file)
        error_regions <- 'ERROR: User-provided regions/features file GFF3 format not valid'
        regions <- if (any(stringr::str_detect(x, error_regions))) NULL else regions
    }

    # Format regions
    if (is_file(regions)) {
        # Read
        gff <- read_gff3(regions)

        # Remove ANY regions identical to those called by Prodigal
        gff$remove <- stringr::str_detect(gff$source, 'Prodigal')
        gff <- dplyr::group_by(gff, seqid, start, end) %>% dplyr::mutate(remove = any(remove))
        index <- !gff$remove
        gff$remove <- NULL
        gff <- gff[index, ]

        # Write
        write_gff3(gff, int.regions)
    }
    
    # Run Bakta
    cmd <- paste('bakta','--db',bakta.db,'--verbose','--force','--output',output.dir,'--threads',threads)
    cmd <- if (length(prefix)) paste(cmd,'-p',prefix) else cmd 
    cmd <- if (length(regions)) paste(cmd,'--regions',int.regions) else cmd 
    cmd <- if (keep.contig.headers) paste(cmd,'--keep-contig-headers') else cmd
    cmd <- paste(cmd,input.fasta)
    stdout <- system3(cmd, log.file=log.file)
}

#' Run Bakta for a genomeCollection
#'
#' @param object A genomeCollection object
#' @param recompute Whether to re-run Bakta even though an output GFF exists
#' @param debug Whether to run in debug mode (run only for debug.n samples)
#' @param debug.n Number of samples to run while debugging
#'
#' @export
#'
RunBakta <- function(object, 
                     regions = 'regions',
                     db = '../databases/bakta/',
                     recompute = FALSE,
                     recompute.sample = NULL,
                     debug = FALSE, debug.n = 1
                    ) {

    # Minimal check
    stopifnot(
        is(object) == 'genomeCollection'
    )

    # Variables
    ind <- file.exists(object$genome)
    out.dir <- paste0(object$path[ind],'annotation/bakta/')
    genomes <- object$genome[ind]
    regions <- if (regions %in% names(object[[]])) object[[regions]] else NULL
    force <- if (recompute) TRUE else FALSE

    # Output
    data <- Annotation(object, 'Bakta')
    if (is.null(data)) {
        prefix <- index(object)[ind]
        data <- methods::new("Annotation", 
                             index = index(object)[ind],
                             genes = paste0(out.dir,prefix,'.gff3'),
                             cds = paste0(out.dir,prefix,'.ffn'),
                             proteins = paste0(out.dir,prefix,'.faa'),
                             log = paste0(out.dir,'runtime.log'),
                             genome = 'genome',
                             tool = 'Bakta',
                             type = "default"
                            )
    }

    # Input
    ind <- file.exists(data$genes)
    missing <- data$index[!ind]
    missing <- if (recompute) data$index else missing
    missing <- if (debug) head(missing, debug.n) else missing
    ind.rc <- recompute.sample %in% data$index
    missing <- if (any(ind.rc)) recompute.sample[ind.rc] else missing
    msg <- paste('Found genes for',sum(ind),'out of',length(data$index),'genomes.','Running Bakta for',length(missing),'...')
    message(msg)
    
    # Main
    for (sample in missing) {
        print(sample)
        i <- which(data$index == sample)
        bakta(input.fasta = genomes[[i]], 
              regions = regions[[i]],
              output.dir = out.dir[[i]],
              prefix = sample,
              db = db,
              force = force
             )
    }
    
    # Assign object
    Annotation(object, 'Bakta') <- data
    
    # Exit
    return(object)
}


#' PhANNs
#'
#' Annotate phage structural proteins
#'
#' @param genome.fasta Path to genome FASTA
#'
#' @export
PhANNs <- function(genome.fasta=NULL, output.dir=NULL, 
                            conda.env='../envs/phanns',
                            repo.dir='../PhANNs',
                            repo.url='https://github.com/Adrian-Cantu/PhANNs.git',
                            model.url='https://edwards.sdsu.edu/phanns/download/model.tar'
                           ) {

    # Minimal check
    stopifnot(
        is_file(genome.fasta),
        !is.null(output.dir)
    )

    # Check program

    return(NULL)
}

#' Pharokka
#'
#' CLI wrapper for phage genome annotation using Pharokka.
#'
#' @param genome.fasta Path to genome FASTA
#'
#' @export
pharokka <- function(genome.fasta=NULL) {

    # Minimal check
    stopifnot(
        is_file(genome.fasta)
    )

    return(NULL)
}


#' PADLOC: Prokaryotic Antiviral Defence LOCator
#' 
#' CLI wrapper for defense gene annotation using PADLOC. Uses either a genome FASTA or
#' a combination of protein FASTA with genes GFF. Database will be downloaded if not present.
#'
#' @param protein.faa Character, path to protein FASTA (only in combination with genes.gff)
#' @param genes.gff Character, path to genomic GFF (only in combination with protein.faa)
#' @param genome.fna Character, path to genomic FASTA
#' @param crispr.gff Character, path to CRISPRfinder GFF.
#' @param db.path Character, path to padloc database.
#' @param overwrite Boolean, whether to overwrite existing results.
#' @param threads Numeric, number of threads. Defaults to nproc.
#' 
#' @export
padloc <- function(protein.faa=NULL, genes.gff=NULL, genome.fna=NULL, 
                   crispr.gff=NULL, 
                   out.dir, 
                   db.path,
                   overwrite=FALSE,
                   threads=n_proc()
                  ) {

    # Minimal check
    check_installed('padloc', silent=TRUE)

    # Variables
    out.dir <- if (endsWith(out.dir,'/')) out.dir else paste0(out.dir,'/')
    log.file <- paste0(out.dir,'runtime.log')
    all_files <- list.files(out.dir)

    # Output
    if (overwrite) unlink(out.dir, recursive = TRUE)
    index <- str_detect(all_files, 'padloc.csv')
    if (any(index)) {
        output_file <- paste0(out.dir,all_files[index])
        msg <- paste('Output',output_file,'already exists. Aborting...')
        message(msg)
        return()
    }
    dir.create(out.dir)

    # Check PADLOC version
    ## Program
    cmd <- paste('padloc','--data',db.path,'--version')
    padloc_version <- try(system(cmd, intern=TRUE))
    padloc_version_number <- str_remove(padloc_version,'padloc v')
    ## Database
    cmd <- paste('padloc','--data',db.path,'--db-version')
    db_version <- try(system(cmd, intern=TRUE))
    db_version_number <- str_remove(db_version,'padloc-db v')
    ## Checks
    if (db_version == 'ERROR: Database version information not found') {
        update_db <- TRUE
        msg <- paste('No PADLOC-DB has been found in.',db.path,'Trying to install it...')
        warning(msg)
    } else if (padloc_version_number > db_version_number) {
        update_db <- TRUE
        msg <- paste('PADLOC-DB version below PADLOC version.')
        warning(msg)
    } else {
        update_db <- FALSE
    }
    msg <- paste('PADLOC version\n',padloc_version,'\n',db_version)
    message(msg)

    if (padloc_version_number < db_version_number) {
        update_db <- FALSE
        msg <- paste('PADLOC-DB is newer than PADLOC. You should fix your PADLOC install.')
        stop(msg)
    }
    
    # Update PADLOC database
    if (update_db) {
        msg <- paste('Updating PADLOC database in',db.path)
        message(msg)
        cmd <- paste('padloc','--db-update','--data',db.path)
        system(cmd, intern=TRUE) %>% paste0(collapse = '\n') %>% cat
    }

    # Check PADLOC dependencies
    cmd <- paste('padloc','--check-deps')
    msg <- try(system(cmd, intern=TRUE)) %>% paste0(collapse = '\n')
    message(msg)
    
    # Base command
    cmd <- paste('padloc','--outdir',out.dir)

    ## Protein
    if (is_file(protein.faa, silent=TRUE)) {
        cmd <- paste(cmd,'--faa',protein.faa)
    }

    ## Genes
    if (is_file(genes.gff, silent=TRUE)) {
        cmd <- paste(cmd,'--gff',genes.gff,'--fix-prodigal')
    }

    ## Genome
    if (is_file(genome.fna, silent=TRUE)) {
        cmd <- paste(cmd,'--fna',genome.fna)
    }

    ## CRISPR array
    if (is_file(crispr.gff, silent=TRUE)) {
        cmd <- paste(cmd,'--crispr',crispr.gff)
    }

    # Main
    cmd <- paste(cmd,'--cpu',threads,'--data',db.path,'--force')
    system3(cmd, log.file = log.file)
}

#' DefenseFinder
#'
#' CLI wrapper for defense gene annotation using DefenseFinder
#'
#' @param file Character, path to the input file.
#' @param out.dir Character, path to target directory.
#' @param db.path Character, path to defense-system-model database (mentioned in
#' https://github.com/mdmparis/defense-finder/issues/8).
#' @param threads Numeric, number of threads.
#' @param update.db Boolean, whether to update the DefenseFinder database
#'
#' @export
defense_finder <- function(file=NULL, out.dir=NULL, db.path=NULL, 
                           anti.defense=TRUE, anti.defense.only=FALSE,
                           threads = n_proc(), update.db=FALSE
                          ) {

    # Minimal check
    stopifnot(
        is_file(file),
        !is.null(out.dir),
        !is.null(db.path),
        length(out.dir) == 1,
        length(file) == 1 # Might be changed later ...
    )
    check_installed('defense-finder', silent=TRUE)
    check_version('defense-finder')    
    
    # Status message
    message(paste('File:',file))

    # Check output 
    out_files <- list.files(out.dir)
    output_present <- str_detect(out_files, 'systems.tsv')
    if (any(output_present)) {
        msg <- paste('Output file',out_files[output_present],'already exists. Aborting...')
        warning(msg)
        return()
    }

    # DB version
    if (!endsWith(db.path,'/')) {
        db.path <- paste0(db.path,'/')
    }
    cmd <- paste0('grep vers: ',db.path,'CasFinder/metadata.yml')
    cat('CasFinder ')
    system3(cmd)
    cat('\n')
    cmd <- paste0('grep vers: ',db.path,'defense-finder-models/metadata.yml')
    cat('DefenseFinder models ')
    system3(cmd)
    cat('\n')

    # Update DB
    db_metadata <- paste0(db.path,'defense-finder-models/metadata.yml')
    if (!file.exists(db_metadata) | update.db) {
        cmd <- paste('defense-finder update','--models-dir',db.path,'2>&1')
        cat(cmd,'\n\n')
        stdout <- system(cmd, intern=TRUE)
        stdout <- paste(stdout, collapse='\n')
        cat(stdout,'\n\n')
    }

    # Antidefense
    if (anti.defense.only) {
        anti_defense <- '-A'
    } else if (anti.defense) {
        anti_defense <- '-a'
    } else {
        anti_defense <- NULL
    }
    
    # Run DefenseFinder
    cmd <- paste('defense-finder run',file,anti_defense,'--out-dir',out.dir,'--models-dir',db.path,'--workers',threads,'2>&1')
    cat(cmd,'\n\n')
    stdout <- system(cmd, intern=TRUE)
    stdout <- paste(stdout, collapse='\n')
    cat(stdout)
}

#' ECTyper
#'
#' ...
ec_typer <- function() {}

#' PSTyper
#'
#' Serotyping Pseudomonas syringae genomes
ps_typer <- function() {}

#' CapsuleFinder
#'
#' implement from https://gitlab.pasteur.fr/gem/capsuledb/-/tree/master?ref_type=heads
capsule_finder <- function() {}

#' Kaptive
#'
#' version 3

#' PhageRBPdetect
#'
#' Detect receptor binding proteins (RBP) in phage genomes
#'
#' @param genome FASTA
#'
#' @export
phage_rbp_detect <- function(genome=NULL) {

    # Minimal check
    stopifnot(
        !is.null(genome)
    )

    # Check output

    # ...
    return(NULL)
}