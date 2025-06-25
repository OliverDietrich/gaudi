#' Bakta annotation
#' 
#' CLI wrapper for genome annotation using Bakta
#'
#' @param input.fasta Character, path to FASTA file
#' @param output.dir Character, path to output dir
#' @param output.prefix Character, 
#' @param db Character, path to the bakta database
#' @param download.db Boolean, whether to download the database if missing
#' @param threads Integer, number of threads to use (default: all)
#' @param force Boolean, whether to force overwriting existing output folder
#'
#' @export
bakta <- function(input.fasta=NULL, output.dir=NULL, output.prefix=NULL, db=NULL, download.db=FALSE, threads=n_proc(), force=FALSE) {

    # Minimal check
    stopifnot(
        is_file(input.fasta, suffix=c('.fasta','.fna')),
        !is.null(output.dir)
    )
    check_installed('bakta', silent=TRUE)
    check_version('bakta')

    # Variables
    download.log <- paste0(db,'download.log')
    if (!endsWith(db,'/')) {
        db <- paste0(db,'/')        
    }
    bakta.db <- paste0(db,'db/')
    bakta.db.version <- paste0(bakta.db,'/version.json')

    # Check output
    if (dir.exists(output.dir) & !force) {
        msg <- paste('Output directory',output.dir,'already exists. Aborting...')
        cat(msg)
        return()
    }

    # Check DB
    if (is.null(db)) {
        msg <- 'Please supply a path for the Bakta database!'
        stop(msg)
    } else if (file.exists(bakta.db.version)) {
        bakta.db.version <- jsonlite::read_json(bakta.db.version)
        msg <- paste0('Bakta database (',bakta.db.version$type,'), version ',bakta.db.version$major,'.',bakta.db.version$minor,' (',bakta.db.version$date,')')
        message(msg)
    } else if (!file.exists(bakta.db) & download.db) {
        msg <- paste('Downloading bakta database to',db)
        message(msg)
        cmd <- paste('bakta_db download','--output',db,'--type','full','2>&1')
        message(cmd)
        system2(cmd, stdout=download.log, stderr=download.log)
    } else {
        msg <- paste('Bakta database',db,'does not exists.','Consider passing download.db=TRUE.')
        stop(msg)
    }

    # Run Bakta
    cmd <- paste('bakta','--db',bakta.db,'--verbose','--force','--output',output.dir,'--threads',threads,input.fasta,'2>&1')
    stdout <- system(cmd, intern=TRUE)
    stdout <- paste(stdout, collapse='\n')
    cat(stdout)
}

#' Prokka
#'
#' ...


#' PhANNs
#'
#' Annotate phage structural proteins
#'
#' @param genome.fasta Path to genome FASTA
#'
#' @export
annotate_PhANNs <- function(genome.fasta=NULL, output.dir=NULL, 
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
                   out.dir=NULL, 
                   db.path=NULL,
                   overwrite=FALSE,
                   threads=n_proc()
                  ) {

    # Minimal check
    stopifnot(
        !is.null(protein.faa) | !is.null(genes.gff) | !is.null(genome.fna),
        !is.null(out.dir),
        !is.null(db.path)
    )
    check_installed('padloc', silent=TRUE)

    # Variables
    if (!endsWith(out.dir,'/')) {
        out.dir <- paste0(out.dir,'/')
    }
    log.file <- paste0(out.dir,'runtime.log')
    all_files <- list.files(out.dir)

    # Output
    index <- str_detect(all_files, 'padloc.csv')
    if (any(index) & !overwrite) {
        output_file <- paste0(out.dir,all_files[index])
        msg <- paste('Output',output_file,'already exists. Aborting...')
        warning(msg)
        return()
    } else if (file.exists(log.file) & !overwrite) {
        log <- readLines(log.file)
        log <- paste(log, collapse='\n')
        msg <- paste('Log file',log.file,'found. Aborting...')
        warning(msg)
        cat(log)
        return()
    } else {
        dir.create(out.dir, showWarnings = FALSE)
    }

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

    # Run PADLOC
    cmd <- paste(cmd,'--cpu',threads,'--data',db.path,'--force','2>&1')
    cat(cmd)
    stdout <- system(cmd, intern=TRUE)
    outsize <- length(stdout)
    stdout <- paste(stdout, collapse='\n')
    writeLines(stdout, log.file)

    # Print statement
    if (outsize < 25) {
        cat(stdout)
    } else {
        msg <- paste('Long stdout, check',log.file)
        warning(msg)
    }
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

#'
#' @export
phANNs_phage_annotation <- function(genome.fasta=NULL) {

    # Minimal check
    stopifnot(
        !is.null(genome)
    )

    return(NULL)
}

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

#' geNomad
#'
#' Plasmid detection
geNomad <- function() {}

#' Virsorter2
#'
#' Prophage detection
virsorter <- function() {}

#' MacSyFinder
#'
#' Macromolecular systems finder
mac_sy_finder <- function() {}