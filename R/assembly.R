#' SKESA assembly
#'
#' CLI wrapper for SKESA assembly from short reads
#'
#' @param input.1 Character, file name of forward paired-end reads
#' @param input.2 Character, file name of reverse paired-end reads
#' @param input.s Character, file name of unpaired reads
#' @param threads Integer, number of threads to use
#' @param overwrite Boolean, whether to overwrite output
#'
#' @export
skesa_assembly <- function() {

    # Minimal check
    stopifnot()
    check_installed('skesa', silent=TRUE)
    check_version('skesa')

    # Check input

    # Check output

    # Run SKESA
    # skesa --reads short_1.fastq.gz,short_2.fastq.gz --cores 16 --memory 64 > skesa_assembly.fasta
    cmd <- 'skesa'
    cmd <- paste(cmd,'2>&1')
    stdout <- system(cmd, intern=TRUE)
    stdout <- paste(stdout, collapse='\n')
    cat(stdout)
}

#' SPAdes assembly
#'
#' CLI wrapper for SPAdes assembly from short reads
#'
#' @param input.1 File name, file with forward paired-end reads
#' @param input.2 File name, file with reverse paired-end reads
#' @param input.s File name, file with unpaired reads
#' @param threads Integer, number of threads to use
#' @param overwrite Boolean, whether to overwrite output
#'
#' @export
spades_assembly <- function(input.1=NULL, input.2=NULL, input.s=NULL, 
                            output=NULL,
                            threads=n_proc(), overwrite=FALSE
                           ) {

    # Minimal check
    stopifnot(
        !is.null(output)
    )
    check_installed('spades.py', silent=TRUE)
    check_version('spades.py')

    # Check output
    if (file.exists(output)) {
        msg <- paste('Output file',output,'already exists. Will be skipped...')
        warning(msg)
        return(1)
    }

    # Check input
    if (!is.null(input.s)) {
        if (is.null(input.1) & is.null(input.2)) {
            input <- paste('-s',input.s)
            msg <- 'Using unpaired reads...'
            message(msg)
        } else {
            msg <- 'Paired-end input supplied together with unpaired input. Aborting.'
            stop(msg)
        }
    }
    if (!is.null(input.1) & !is.null(input.2)) {
        input <- paste('-1',input.1,'-2',input.2)
        msg <- 'Using paired-end reads...'
        message(msg)
    } else {
        msg <- 'For paired-end reads both input.1 AND input.2 must be supplied. Aborting.'
        stop(msg)
    }

    # Run SPAdes
    cmd <- 'spades.py'
    cmd <- paste(cmd, input,'-o', output,'--threads',threads)
    cmd <- paste(cmd,'2>&1')
    stdout <- system(cmd, intern=TRUE)
    stdout <- paste(stdout, collapse='\n')
    cat(stdout)
    return(0)
}

#' Unicycler assembly
#'
#' CLI wrapper for Unicycler genome assembly.
#' Ideally hybrid (long + short) but accepts both short and long reads only.
#' 
#' @param short.1 FASTQ file of first short reads in each pair (R1)
#' @param short.2 FASTQ file of second short reads in each pair (R2)
#' @param unpaired FASTQ file of unpaired short reads (S)
#' @param long FASTQ or FASTA file of long reads
#' @param out.dir Output directory
#' @param keep Level of file retention 0 = only keep final files: assembly (FASTA, GFA and log), 
#' 1 = also save graphs at main checkpoints, 2 = also keep SAM (enables fast rerun in different mode),3 = keep all temp files and save all graphs (for debugging)
#' @param min_fasta_length Exclude contigs from the FASTA file which are shorter than this length (default: 100)
#' @param mode Bridging mode: conservative, normal, bold (default: normal)
#' @param linear_seqs The expected number of linear (i.e. non-circular) sequences in the underlying sequence (default: 0)
#' @param threads Number of threads used
#'
#' @export
unicycler_assembly <- function(short.1=NULL, short.2=NULL, unpaired=NULL, 
                               long=NULL, out.dir, 
                               keep = 1, mode = 'normal',
                               min_fasta_length=100,
                               linear_seqs = 0,
                               threads = n_proc()
                              ) {

    # Minimal check

    # Variables

    # Check output

    # Check input

    # Main
    stdout <- system3('unicycler')
    
    # Format output
    cat(stdout)
}

#' Raven assembly
#' 
#' CLI wrapper for Raven to assemble genomes from long reads
#' 
#' @param input.fastq Fastq files to assemble
#' @param output.fasta Fasta file for assembly
#' @param threads Number of CPU cores to use
#'
#' @export
raven_assembly <- function(input.fastq, output.fasta,
                           threads=n_proc()
                          ) {

    # Minimal check
    stopifnot(
        !is.null(input.fastq),
        !is.null(output.fasta)
        #is_valid_fastq(input.fastq)
    )
    check_installed('raven', silent=TRUE)
    check_version('raven')

    # Check output
    ## Exit 1
    if (file.exists(output.fasta)) {
        msg <- paste('Assembly',output.fasta,'exists.')
        warning(msg)
        return('Exit status 1: Output exists.')
    }
    
    # Run Raven
    cmd <- paste('raven','--threads',threads,input.fastq,'>',output.fasta,'2>&1')
    message(cmd)
    stdout <- system(cmd, intern=TRUE)
    stdout <- paste(stdout, collapse='\n')
    cat(stdout)
}

#' Flye assembly
#'
#' CLI wrapper for Flye to assemble genomes from long reads
#'
#' @param pacbio.raw Path, PacBio regular CLR reads (<20% error)
#' @param pacbio.corr Path, PacBio reads that were corrected with other methods (<3% error)
#' @param pacbio.hifi Path, PacBio HiFi reads (<1% error)
#' @param nano.raw Path, ONT regular reads, pre-Guppy5 (<20% error)
#' @param nano.corr Path, ONT reads that were corrected with other methods (<3% error)
#' @param nano.hq Path, ONT high-quality reads: Guppy5+ SUP or Q20 (<5% error)
#' @param genome.size Character, estimated genome size (for example, 5m or 2.6g)
#' @param out.dir Path, Output directory
#' @param threads Integer,  number of parallel threads
#'
#' @export
flye_assembly <- function(pacbio.raw=NA, pacbio.corr=NA, pacbio.hifi=NA, nano.raw=NA, nano.corr=NA, nano.hq=NA, # Input formats
                          genome.size=NULL, 
                          out.dir=NULL, 
                          threads=n_proc(), 
                          help=FALSE
                         ) {

    # Variables
    input_files <- setNames(
        c(pacbio.raw,pacbio.corr,pacbio.hifi,nano.raw,nano.corr,nano.hq),
        c('pacbio.raw','pacbio.corr','pacbio.hifi','nano.raw','nano.corr','nano.hq')
    )
    if (!endsWith(out.dir, '/')) {
        out.dir <- paste0(out.dir,'/')
    }
    out_final <- paste0(out.dir, 'assembly.fasta')
    
    # Minimal check
    stopifnot(
        !is.null(out.dir)
    )
    check_installed('flye', silent=TRUE)
    check_version('flye')

    # Check input
    input_present <- !is.na(input_files)
    if (all(!input_present)) {
        msg <- paste('Please provide an input file:', paste(names(input_files), collapse=', '))
        stop(msg)
    } else if (sum(input_present) > 1) {
        msg <- paste('Multiple inputs present. Select one of:', paste(names(input_files[input_present]), collapse=', '))
        stop(msg)
    }
    input_file <- input_files[input_present]
    input_type <- str_replace(names(input_files[input_present]), '\\.', '-')

    # Check output
    if (!dir.exists(out.dir)) {
        msg <- paste('Directory', out.dir, 'does not exits. Will be created...')
        dir.create(out.dir, recursive=FALSE)
        warning(msg)
    }

    # Exit 1
    if (file.exists(out_final)) {
        msg <- paste('Assembly exists as',out_final,'and will be skipped.')
        warning(msg)
        return('Exit status 1: Output exists.')
    }

    # Run Flye
    cmd <- paste0('flye',' --',input_type,' ',input_file,' --out-dir ',out.dir,' --threads ',threads)
    if (!is.null(genome.size)) {
        cmd <- paste0(cmd,' --genome-size ',genome.size)
    }
    cmd <- paste0(cmd,' 2>&1')
    message(cmd)
    stdout <- system(cmd, intern=TRUE)
    stdout <- paste(stdout, collapse='\n')
    cat(stdout)
}

#' Medaka polishing
#'
#' CLI wrapper for Medaka polishing of genome assemblies.
#' Medaka is optimized to work with the Flye assembler.
#' 
#' @param input.fastx FASTx input basecalls
#' @param input.assembly FASTA input assembly
#' @param out.dir Output folder
#' @param fill_gaps Boolean, whether to fill gaps in consensus with draft sequence
#' @param model Medaka model (see "medaka_consensus -h" for choices)
#' @param threads Number of threads with which to create features
#' @param batch_size Numeric, select batch size to control memory use
#'
#' @export
medaka_polish <- function(input.fastx=NULL, input.assembly=NULL, out.dir=NULL, 
                          fill_gaps=TRUE, model = 'r1041_e82_400bps_sup_v5.0.0', 
                          threads=1
                         ) {

    # Minimal check

    # Variables

    # Check output

    # Check input

    # Main
    stdout <- system3('medaka_consensus')

    # Format output
    cat(stdout)
}