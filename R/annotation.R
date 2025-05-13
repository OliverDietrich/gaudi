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