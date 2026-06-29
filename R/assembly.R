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
#' @param out.dir Path to output directory
#' @param threads Integer, number of threads to use
#' @param overwrite Boolean, whether to overwrite output
#'
#' @export
spades_assembly <- function(input.1=NULL, input.2=NULL, input.s=NULL, 
                            out.dir=NULL,
                            isolate = FALSE,
                            single.cell = FALSE,
                            metagenomic = FALSE,
                            biosynthetic = FALSE,
                            sewage = FALSE,
                            corona = FALSE,
                            rna = FALSE,
                            plasmid = FALSE,
                            metaviral = FALSE,
                            metaplasmid = FALSE,
                            rnaviral = FALSE,
                            log.file = NULL,
                            threads=n_proc(), 
                            overwrite=FALSE
                           ) {

    # Minimal check
    stopifnot(
        !is.null(out.dir)
    )
    check_installed('spades.py', silent=TRUE)
    check_version('spades.py')

    # Variables
    out.dir <- if(endsWith(out.dir, '/')) out.dir else paste0(out.dir, '/');
    log.file <- if (is.null(log.file)) paste0(out.dir, 'runtime.log') else log.file;
    contigs.fasta <- paste0(out.dir, 'contigs.fasta')

    # Check output
    if (file.exists(contigs.fasta) & !overwrite) {
        msg <- paste('Output file',contigs.fasta,'already exists. Will be skipped...')
        warning(msg)
        return()
    } else if (overwrite) {
        unlink(out.dir, recursive = TRUE)
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

    # Add flags
    flags <- list()
    flags$isolate <- if(isolate) '--isolate' else NULL;
    flags$sc <- if(single.cell) '--sc' else NULL;
    flags$meta <- if(metagenomic) '--meta' else NULL;
    flags$bio <- if(biosynthetic) '--bio' else NULL;
    flags$sewage <- if(sewage) '--sewage' else NULL;
    flags$corona <- if(corona) '--corona' else NULL;
    flags$rna <- if(rna) '--rna' else NULL;
    flags$plasmid <- if(plasmid) '--plasmid' else NULL;
    flags$metaviral <- if(metaviral) '--metaviral' else NULL;
    flags$metaplasmid <- if(metaplasmid) '--metaplasmid' else NULL;
    flags$rnaviral <- if(rnaviral) '--rnaviral' else NULL;
    flags <- paste(unlist(flags), collapse=' ')
    flags <- if (flags == '') NULL else flags;

    # Main
    cmd <- paste('spades.py', input,'-o', out.dir,'--threads',threads, flags)
    system3(cmd, log.file = log.file)
}

#' RunSPAdes
#'
#' Run SPAdes assembly for a genomeCollection
#' 
#' @param object A genomeCollection object
#' @param name Slot to use for Assemblies(object)
#' @param type Type of reads to use for assembly. One of paired, unpaired.
#' @param reads.name Slot in ReadsNames(object
#' @param threads Number of cores to use
#'
#' @export
#'
RunSPAdes <- function(object, name = 'SPAdes',
                      type = 'paired',
                      reads.name='filtered',
                      recompute = FALSE,
                      recompute.sample = NULL,
                      ...
                     ) {
    
    # Minimal check
    stopifnot(
        reads.name %in% ReadsNames(object)
    )

    # Fetch data
    reads <- Reads(object, reads.name)
    
    # Check input
    if (type == 'paired') {
        ind.reads <- file.exists(reads$R1) & file.exists(reads$R2)
        input.1 <- reads$R1[ind.reads]
        input.2 <- reads$R2[ind.reads]
        input.s <- NULL
    } else 
    if (type == 'unpaired') {
        ind.reads <- file.exists(reads$S)
        input.1 <- NULL
        input.2 <- NULL
        input.s <- reads$S[ind.reads]
    } else {
        msg <- paste('Type', type, 'not found. Options: paired, unpaired.')
        stop(msg)
    }
    msg <- paste('Found', sum(ind.reads), 'samples with', type, 'reads.')
    message(msg)

    # Variables
    samples <- reads$index[ind.reads]
    path <- object$path[match(samples, index(object))]
    output <- paste0(path, 'assembly/SPAdes/')
    contigs <- paste0(output, 'contigs.fasta')
    graphs <- paste0(output, 'assembly_graph_after_simplification.gfa')
    log <- paste0(output, 'runtime.log')

    # Create object
    data <- Assembly(object, name)
    if (is.null(data)) {
        data <- methods::new("Assembly", 
                             "index" = samples,
                             "contig" = contigs,
                             "graph" = graphs,
                             "log" = log,
                             "reads" = reads.name,
                             "tool" = 'SPAdes',
                             "type" = type
                            )
    }
    check <- all(data$index == samples) & all(data$contig == contigs) & all(data$log == log) & data$reads == reads.name
    if (!check) {
        msg <- paste0('Some slots in Assembly(object, "', name, '") do not contain the correct data. Aborting...')
        stop(msg)
    }

    # Check output
    ind <- file.exists(contigs)
    missing <- if (recompute) samples else samples[!ind]
    missing <- if (length(recompute.sample)) recompute.sample else missing # Debug mode
    msg <- paste('Assembly found for', sum(ind), 'samples. Running SPAdes for', length(missing), 'samples...')
    message(msg)
    
    # Main
    for (sample in samples) {
        if (!sample %in% missing) next
        cat(sample)
        i <- which(samples == sample)
        spades_assembly(input.1 = input.1[[i]], input.2 = input.2[[i]], input.s = input.s[[i]], 
                        out.dir = output[[i]], 
                        log.file = log[[i]], 
                        overwrite = TRUE, 
                        ...
                       )
    }

    # Replace assembly
    Assembly(object, name) <- data
    return(object)
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
unicycler_assembly <- function(short.1=NULL, short.2=NULL, unpaired=NULL, long=NULL, 
                               out.dir=NULL, 
                               log.file=NULL,
                               keep = 1, 
                               mode = 'normal',
                               min_fasta_length=100,
                               linear_seqs = 0,
                               threads = n_proc(), 
                               overwrite = FALSE
                              ) {

    # Minimal check
    stopifnot(
        !all(is.null(c(short.1, short.2, unpaired, long))),
        !is.null(out.dir)
    )

    # Variables
    out.dir <- if (endsWith(out.dir, '/')) out.dir else paste0(out.dir, '/')
    out.fasta <- paste0(out.dir, 'assembly.fasta')
    out.gfa <- paste0(out.dir, 'assembly.gfa')

    # Check output
    if (overwrite) {
        unlink(out.dir, recursive=TRUE)
    }
    if (!dir.exists(out.dir)) dir.create(out.dir, recursive=TRUE)

    # Exit 1
    if (file.exists(out.fasta)) {
        msg <- paste('Output file', out.fasta, 'already exists.')
        warning(msg)
        return()
    }

    # Exit 2

    # Check input
    short.1 <- if (is.na(short.1)) NULL else short.1
    short.2 <- if (is.na(short.2)) NULL else short.2
    unpaired <- if (is.na(unpaired)) NULL else unpaired
    long <- if (is.na(long)) NULL else long
    if (length(short.1) & length(short.2) & length(unpaired)) {
        msg <- 'Found both paired and unpaired reads. This might fail...'
        warning(msg)
    }
    paired <- if (length(short.1) & length(short.2)) paste('--short1',short.1,'--short2',short.2) else NULL
    unpaired <- if (length(unpaired)) paste('--unpaired',unpaired) else NULL
    long <- if (length(long)) paste('--long',long) else NULL

    # Main
    cmd <- paste('unicycler', paired, unpaired, long, '--out', out.dir, '--keep', keep, '--mode', mode, '--linear_seqs', linear_seqs, '-t', threads)
    system3(cmd, log.file = log.file)

    # Exit 0
    return()
}

#' Run Unicycler
#'
#' Run Unicycler assembly for a genomeCollection
#' 
#' @param object A genomeCollection object
#' @param name Slot to use for Assemblies(object)
#' @param type Type of reads to use for assembly. One of paired, unpaired.
#' @param reads.name Slot in ReadsNames(object
#' @param threads Number of cores to use
#'
#' @export
#'
RunUnicycler <- function(object, name = 'Unicycler',
                         reads.name = 'filtered',
                         recompute = FALSE,
                         recompute.sample = NULL,
                         samples.skip = NULL,
                         debug = FALSE,
                         debug.n = 1,
                         ...
                        ) {

    # Minimal check
    stopifnot(
        reads.name %in% ReadsNames(object)
    )

    # Fetch data
    reads <- Reads(object, reads.name)

    # Variables
    samples <- reads$index
    path <- object$path[match(samples, index(object))]
    output <- paste0(path, 'assembly/',name,'/')
    contigs <- paste0(output, 'assembly.fasta')
    graphs <- paste0(output, 'assembly.gfa')
    log <- paste0(output, 'runtime.log')

    # Create object
    data <- Assembly(object, name)
    if (is.null(data)) {
        data <- methods::new("Assembly", 
                             "index" = samples,
                             "contig" = contigs,
                             "graph" = graphs,
                             "log" = log,
                             "reads" = reads.name,
                             "tool" = 'Unicycler',
                             "type" = 'Will auto-detect type for each sample...'
                            )
    }
    check <- all(data$index == samples) & all(data$contig == contigs) & all(data$log == log) & data$reads == reads.name
    if (!check) {
        msg <- paste0('Some slots in Assembly(object, "', name, '") do not contain the correct data. Aborting...')
        stop(msg)
    }

    # Check output
    ind <- file.exists(contigs)
    missing <- if (recompute) samples else samples[!ind]
    missing <- if (length(recompute.sample)) recompute.sample else missing
    missing <- if (debug) head(missing, debug.n) else missing
    msg <- paste('Assembly found for', sum(ind), 'samples. Running Unicycler for', length(missing), 'samples...')
    message(msg)

    # Main
    for (sample in samples) {
        if (!sample %in% missing) next
        if (sample %in% samples.skip) next
        cat(sample,'\n')
        i <- which(samples == sample)        
        unicycler_assembly(short.1 = reads$R1[i], short.2 = reads$R2[i], unpaired = reads$S[i], long = reads$L[i],
                        out.dir = output[[i]],
                        log.file = log[[i]], 
                        overwrite = TRUE, 
                        ...
                       )
    }

    # Replace assembly
    Assembly(object, name) <- data
    return(object)
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
                          log.file = NULL,
                          threads=n_proc(), 
                          overwrite=FALSE
                         ) {

    # Minimal check
    stopifnot(
        !is.null(out.dir)
    )
    check_installed('flye', silent=TRUE)
    check_version('flye')

    # Variables
    input_files <- setNames(
        c(pacbio.raw,pacbio.corr,pacbio.hifi,nano.raw,nano.corr,nano.hq),
        c('pacbio.raw','pacbio.corr','pacbio.hifi','nano.raw','nano.corr','nano.hq')
    )
    out.dir <- if (endsWith(out.dir, '/')) out.dir else paste0(out.dir,'/')
    out_final <- paste0(out.dir, 'assembly.fasta')
    log.file <- if(is.null(log.file)) paste0(out.dir, 'runtime.log') else log.file

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
    if (overwrite) unlink(out.dir, recursive=TRUE)
    if (!dir.exists(out.dir)) {
        msg <- paste('Directory', out.dir, 'does not exist. Will be created...')
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
    cmd <- if (!is.null(genome.size)) paste0(cmd,' --genome-size ',genome.size) else cmd
    system3(cmd, log.file = log.file)
}

#' Run Flye
#'
#' Run Flye assembly for a genomeCollection
#' 
#' @param object A genomeCollection object
#' @param name Slot to use for Assemblies(object)
#' @param type Type of reads used for assembly. One of pacbio.raw, pacbio.corr, pacbio.hifi, nano.raw, nano.corr, nano.hq.
#' @param reads.name Slot in ReadsNames(object
#' @param threads Number of cores to use
#' @param ... Other arguments to flye_assembly()
#'
#' @export
#'
RunFlye <- function(object, 
                    name = 'Flye',
                    reads.name = 'filtered',
                    type,
                    recompute = FALSE,
                    recompute.sample = NULL,
                    samples.skip = NULL,
                    debug = FALSE,
                    debug.n = 1,
                    ...
                   ) {

    # Minimal check
    stopifnot(
        reads.name %in% ReadsNames(object),
        type %in% c('pacbio.raw', 'pacbio.corr', 'pacbio.hifi', 'nano.raw', 'nano.corr', 'nano.hq')
    )

    # Fetch data
    reads <- Reads(object, reads.name)

    # Variables
    samples <- reads$index
    path <- object$path[match(samples, index(object))]
    output <- paste0(path, 'assembly/Flye/')
    contigs <- paste0(output, 'assembly.fasta')
    graphs <- paste0(output, 'assembly_graph.gfa')
    log <- paste0(output, 'runtime.log')

    # Check input
    if (!length(reads$L)) {
        msg <- paste0('No long reads found in Reads(object, "', reads.name, '") but required for Flye assembly. Aborting...')
        stop(msg)
    }

    # Create object
    data <- Assembly(object, name)
    if (is.null(data)) {
        data <- methods::new("Assembly", 
                             "index" = samples,
                             "contig" = contigs,
                             "graph" = graphs,
                             "log" = log,
                             "reads" = reads.name,
                             "tool" = 'Flye',
                             "type" = type
                            )
    }
    check <- all(data$index == samples) & all(data$contig == contigs) & all(data$log == log) & data$reads == reads.name
    if (!check) {
        msg <- paste0('Some slots in Assembly(object, "', name, '") do not contain the correct data. Aborting...')
        stop(msg)
    }

    # Check output
    ind <- file.exists(contigs)
    missing <- if (recompute) samples else samples[!ind]
    missing <- if (length(recompute.sample)) recompute.sample else missing # Debug mode
    missing <- if (debug) head(missing, debug.n) else missing
    msg <- paste('Assembly found for', sum(ind), 'samples. Running Flye for', length(missing), 'samples...')
    message(msg)

    # Main
    for (sample in samples) {
        if (!sample %in% missing) next
        if (sample %in% samples.skip) next
        cat(sample)
        i <- which(samples == sample)
        if (type == 'pacbio.raw') flye_assembly(pacbio.raw = reads$L[[i]], out.dir = output[[i]], log.file = log[[i]], overwrite = TRUE, ...)
        if (type == 'pacbio.corr') flye_assembly(pacbio.corr = reads$L[[i]], out.dir = output[[i]], log.file = log[[i]], overwrite = TRUE, ...)
        if (type == 'pacbio.hifi') flye_assembly(pacbio.hifi = reads$L[[i]], out.dir = output[[i]], log.file = log[[i]], overwrite = TRUE, ...)
        if (type == 'nano.raw') flye_assembly(nano.raw = reads$L[[i]], out.dir = output[[i]], log.file = log[[i]], overwrite = TRUE, ...)
        if (type == 'nano.corr') flye_assembly(nano.corr = reads$L[[i]], out.dir = output[[i]], log.file = log[[i]], overwrite = TRUE, ...)
        if (type == 'nano.hq') flye_assembly(nano.hq = reads$L[[i]], out.dir = output[[i]], log.file = log[[i]], overwrite = TRUE, ...)
    }

    # Replace assembly
    Assembly(object, name) <- data
    return(object)
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

#' Polypolish
#'
#' CLI wrapper for Polypolish for short-read polishing of long-read assemblies.
#' 
#' @param assembly FASTA input assembly
#' @param out.dir Output folder
#' @param short.1 FASTQ file of first short reads in each pair (R1)
#' @param short.2 FASTQ file of second short reads in each pair (R2)
#' @param unpaired FASTQ file of unpaired short reads (S)
#' @param long FASTQ or FASTA file of long reads
#'
#' @export
polypolish <- function(assembly, out.dir,
                       short.1=NULL, short.2=NULL, unpaired=NULL,
                       threads = n_proc(),
                       overwrite = FALSE
                      ) {

    # Minimal check
    stopifnot(
        is_file(assembly)
    )
    check_installed('polypolish', silent=TRUE)
    check_version('polypolish')

    # Variables
    out.dir <- if (endsWith(out.dir, '/')) out.dir else paste0(out.dir,'/')
    dir.create(out.dir, recursive=TRUE)
    log.align.unpaired <- paste0(out.dir,'minimap_unpaired.log')
    log.align.paired <- paste0(out.dir,'minimap_paired.log')
    out.file <- paste0(out.dir,'polished.fasta')
    paired.alignment.sam <- NULL
    unpaired.alignment.sam <- NULL

    # Output
    if (file.exists(out.file) & !overwrite) {
        msg <- paste('Output file', out.file, 'already exists.')
        warning(msg)
        return()
    }

    # Timestamp
    t.start <- Sys.time()

    # Alignment
    ## None
    if (is.null(short.1) & is.null(short.2) & is.null(unpaired)) stop('Please supply any short reads.')
    # Unpaired
    if (length(unpaired) == 1) {
        msg <- paste('Unpaired read file',unpaired,'does not exist.')
        if (!file.exists(unpaired)) stop(msg)
        unpaired.alignment.sam <- paste0(out.dir,'unpaired.sam')
        message('Running minimap2 for unpaired (S) reads...')
        cmd <- paste('minimap2','-t',threads,'-a','-x','sr',assembly,unpaired,'-o',unpaired.alignment.sam)
        system3(cmd) # , log.file = log.align.unpaired
    }
    ## Paired
    if (length(short.1) | length(short.2)) {
        if (is.null(short.1) | is.null(short.2)) stop('Please supply both forward and reverse reads for paired reads.')
        unpaired <- c(short.1, short.2)
        ind <- file.exists(unpaired)
        msg <- paste('Paired read file(s)',paste(unpaired[ind], collapse=', '),'does not exist.')
        if (!all(ind)) stop(msg)
        paired.alignment.sam <- paste0(out.dir,'paired.sam')
        message('Running minimap2 for paired (R1+R2) reads...')
        cmd <- paste('minimap2','-t',threads,'-a','-x','sr',assembly,short.1,short.2,'-o',paired.alignment.sam)
        system3(cmd) # , log.file = log.align.paired
    }    

    # Main
    message('Running polypolish...')
    cmd <- paste('polypolish','polish',assembly,paired.alignment.sam,unpaired.alignment.sam,'>',out.file)
    system3(cmd, include.errors = FALSE)

    # Cleanup
    unlink(paired.alignment.sam)
    unlink(unpaired.alignment.sam)

    # Timestamp
    t.stop <- Sys.time()
    print(t.stop - t.start)
}

#' Assembly summary
#' 
#' Create a summary of the contigs contained in an Assembly object.
#'
#' @param object genomeCollection
#' @param name Slot name of Assembly object in genomeCollection
#' @param filter Data.frame with filters to apply. Columns: sample, length.min, length.max, cov.min. (Not checked, experimental...)
#' @param samples Names of samples to write genomes for. Used to decide for different assemblies in different samples.
#' @param min.size Minimum size of contigs, either single value or threshold for each sample
#' @param min.coverage Minimum coverage of contigs, either single value or threshold for each sample
#' @param max.size Maximum size of contigs, either single value or threshold for each sample
#' @param return.filter Whether to return the filter data.frame. Used to adjust thresholds for individual samples and pass to 'filter' argument in second iteration.
#' @param scaffold Whether to create scaffolds from groups of contigs 
#' param debug Whether to run in debug mode. Returns a list containing the contigs, filter, and sequences.
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_line geom_point facet_wrap theme theme_classic guides element_line guide_legend
#'
#' @export
#'
WriteAssemblyToGenome <- function(object, name,
                                  filter = NULL,
                                  samples = NULL,
                                  min.size = 0, 
                                  min.coverage = 0, 
                                  max.size = Inf,
                                  return.filter = FALSE,
                                  scaffold = FALSE,
                                  scaffold.length = 50,
                                  debug = FALSE
                                 ) {

    # Minimal check
    stopifnot(
        class(object) == 'genomeCollection',
        name %in% Assemblies(object)
    )

    # Variables
    sep.scaffold <- paste(rep('N', scaffold.length), collapse='')
    slot <- if (scaffold) 'scaffold' else 'genome';

    # Input
    data <- Assembly(object, name)
    samples <- if (is.null(samples)) data$index else samples;
    if (!all(samples %in% data$index)) stop('Not all samples are part of Assembly(object, name). Aborting...')

    # Output
    if (scaffold & is.null(object[['scaffold']])) stop('The slot object$scaffold is empty. Please provide scaffold file paths.')

    # Subset samples
    ind <- match(samples, data$index)
    contigs <- setNames(data$contig, data$index)[ind]
    graphs <- setNames(data$graph, data$index)[ind]

    # Require graph
    ind <- file.exists(graphs) # NO SUPPORT FOR ASSEMBLIES WITHOUT GRAPH !!! (yet)
    msg <- paste('Samples', paste(samples[!ind], collapse=', '), 'have no graph and will be removed.')
    if (sum(!ind) > 0) warning(msg)
    sample <- samples[ind]
    contigs <- contigs[ind]
    graphs <- graphs[ind]

    # Fetch components
    graph <- lapply(graphs, read_assembly_graph, remove.seqs = FALSE)
    #seqs <- purrr::map(graph, 'sequence')
    #seqs <- Biostrings::DNAStringSet(do.call('c', lapply(seqs, as.character)))
    contigs <- purrr::map(graph, 'segments')
    contigs <- dplyr::bind_rows(contigs, .id = 'Sample')
    contigs <- as.data.frame(contigs)
    contigs$sample_group <- paste0(contigs$Sample,'_',contigs$Group)

    # Filter
    if (!length(filter)) {
        filter <- data.frame(
            'sample' = data$index,
            'length.min' = min.size,
            'length.max' = max.size,
            'cov.min' = min.coverage
        )
    }
    contigs <- merge(contigs, filter, by.x = 'Sample', by.y = 'sample')
    contigs$filter <- 'Keep'
    contigs$filter[contigs$Group_Length < contigs$length.min] <- 'Discard'
    contigs$filter[contigs$Group_Length > contigs$length.max] <- 'Discard'
    contigs$filter[contigs$Group_Coverage < contigs$cov.min] <- 'Discard'
    contigs$filter[contigs$Group_Coverage > contigs$cov.max] <- 'Discard'

    # Order
    contigs$Sample <- factor(contigs$Sample, index(object))
    
    # Plot
    plot <- ggplot(contigs, aes(Group_Length, Coverage, fill = filter, group = Group, shape = Circular)) +
      geom_vline(xintercept = 5e3, linetype = 'dashed') +
      geom_vline(xintercept = 2e5, linetype = 'solid', linewidth = 2, col = 'grey') +
      #geom_vline(xintercept = 3e5) +
      geom_vline(xintercept = 6e6, linetype = 'dotted') +
      geom_point(aes(size = Length), stroke=.1) +
      geom_line() +
      facet_wrap(~Sample, ncol = 5) +
      ggplot2::scale_size(limits = c(1, NA), range = c(.25,5), trans = 'log10') +
      #ggplot2::scale_size_area(max_size = 5, trans = 'log10') +
      ggplot2::scale_shape_manual(values = c('TRUE'=21, 'FALSE'=24)) +
      ggplot2::scale_y_continuous(trans = 'log10') +
      ggplot2::scale_x_continuous(trans = 'log10') +
      ggplot2::scale_fill_manual(values = c('Keep'='cyan3','Discard'='darkorange')) +
      theme_classic(20) +
      theme(
          panel.grid.major = ggplot2::element_line(color = 'grey', linewidth = .25),
          panel.grid.minor = ggplot2::element_line(color = 'grey', linewidth = .1)
      ) +
      guides(
          fill = guide_legend(override.aes = list(size = 5, shape = 21)),
          shape = guide_legend(override.aes = list(size = 5))
      )
    suppressMessages(print(plot))

    # Substitute numbers for letters
    # First set
    all_groups <- LETTERS
    
    # Second set
    x <- expand.grid(LETTERS, LETTERS)
    x <- x[, rev(names(x))]
    all_groups <- c(all_groups, apply(x, 1, paste, collapse=''))
    
    # Third set
    x <- expand.grid(LETTERS, LETTERS, LETTERS)
    x <- x[, rev(names(x))]
    all_groups <- c(all_groups, apply(x, 1, paste, collapse=''))

    # Apply filter
    index <- contigs$filter == 'Keep'
    contigs <- contigs[index, ]
    #seqs <- seqs[index]

    # Set group names
    groups <- dplyr::group_by(contigs, sample_group, Sample, Group)
    groups <- dplyr::summarize(groups, 
                        Coverage = sum(Length * Coverage) / sum(Length), 
                        Length = sum(Length),
                        Contigs = sum(N)
                       )
    groups <- dplyr::arrange(groups, Sample, dplyr::desc(Length))
    groups <- dplyr::mutate(dplyr::group_by(groups, Sample),
                            group = all_groups[1:length(Group)]
                           )
    if (any(is.na(groups$group))) stop('NA introduced in groups. Over 18278 contigs in one sample, please check input!')
    lookup <- setNames(groups$group, groups$sample_group)
    contigs$newgroup <- lookup[contigs$sample_group]

    # Order
    index <- order(contigs$Sample, contigs$newgroup, contigs$Length, 
                   decreasing = c(FALSE, FALSE, TRUE), 
                   method = 'radix')
    # contigs <- contigs[index, ] # DEPRECATED, moved ordering to read_assembly_graph.

    # Re-name contigs
    contigs <- mutate(dplyr::group_by(contigs, Sample, newgroup),
       newcontig = 1:length(newgroup)
      )

    # Set sequence names
    seqs <- Biostrings::DNAStringSet(contigs$Sequence)
    names(seqs) <- paste0(contigs$Sample,'_',contigs$newgroup,contigs$newcontig,
                          ' ','sample=',contigs$Sample,
                          ' ','group=',contigs$newgroup,
                          ' ','contig=',contigs$newcontig,
                          ' ','length=',contigs$Length,
                          ' ','cov=',contigs$Coverage,
                          ' ','circular=',contigs$Circular
                         )
    contigs$Sequence <- NULL

    # Check
    test <- contigs$Length == sapply(seqs, length)
    if (!all(test)) {
        msg <- 'Not all sequences have the same length as stated in their header. Aborting...'
        stop(msg)
    }

    # Scaffold groups
    if (scaffold) {
        message('Scaffolding...')

        # Concatenate sequences
        key <- paste(contigs$Sample, contigs$newgroup, sep = '__')
        seq.scaffold <- as.list(split(seqs, key))
        seq.names <- names(seq.scaffold)
        for (i in seq.names) {
            s <- lapply(seq.scaffold[[i]], as.character)
            s <- paste(s, collapse = sep.scaffold)
            names(s) <- i
            seq.scaffold[[i]] <- s
        }
        seq.scaffold <- unlist(seq.scaffold)
        names(seq.scaffold) <- seq.names
        seqs <- DNAStringSet(seq.scaffold)

        # Summarize contigs
        contigs <- group_by(contigs, Sample, newgroup)
        contigs <- summarize(contigs, 
                             Name = unique(newgroup),
                             Coverage = unique(Group_Coverage), 
                             Contigs = unique(Group_Length),
                             Length = sum(Length)
                            )
    }
    
    # Write FASTA
    msg <- 'Writing FASTA files...'
    message(msg)
    for (i in samples) {
        ind <- contigs$Sample == i
        fn <- object[[slot]][index(object) == i]
        if (length(seqs[ind]) > 0) Biostrings::writeXStringSet(seqs[ind], fn) else unlink(fn)
    }

    # Exit
    if (debug) {
        output <- list(
            'filter' = filter,
            'contigs' = contigs,
            'seqs' = seqs
        )
        return(output)
    }
    if (return.filter) {
        return(filter)
    }
}

#' Reorient genomes by Dnaapler
#'
#' @param file Path to input FASTA
#' @param output.file Path to output FASTA
#' @param tmp.dir Path to temporary directory for intermediate files. MMseqs requires permission to execute files (u=rwx) which on some HPC systems can be restricted to /tmp.
#'
dnaapler_reorient <- function(file, output.file, tmp.dir = tempdir()) {

    # Check output
    if (file.exists(output.file)) {
        message('Output file already exists.')
        return()
    }

    # Program
    check_installed('dnaapler')

    # Check input
    if (!file.exists(file)) stop('Input file not found.')

    # Variables
    dir.create(tmp.dir)    
    OUTDIR <- paste0(tmp.dir,'/dnaapler')

    # Main
    cmd <- paste(
        'dnaapler','all',
        '-i',file,
        '-o',OUTDIR,
        '-t',n_proc(),
        '-f'
    )
    system3(cmd)

    # Copy output from temporary directory
    file.copy(paste0(OUTDIR,'/','dnaapler_reoriented.fasta'), output.file, overwrite = TRUE)

    # Exit 0
    message('Done.')
}

#' Run Dnaapler
#'
#' Run dnnaapler_reorient for a genomeCollection
#'
#' @param object genomeCollection
#' @param name.from meta.data slot to use as input FASTA
#' @param name.to meta.data slot to store output FASTA
#' @param recompute Whether to remove output files and re-compute
#' @param recompute.sample Names of samples (must match index(object)) to recompute output for
#'
RunDnaapler <- function(object, name.from = 'genome', name.to = 'genome_reoriented',
                        recompute = FALSE,
                        recompute.sample = NULL,
                        samples.skip = NULL,
                        debug = FALSE,
                        debug.n = 2
                       ) {

    # Check input
    if (class(object) != 'genomeCollection') stop('Object must be a genomeCollection.')
    if (is.null(object[[name.from]])) stop(paste('Slot', name.from, 'not found in metadata(object).'))
    if (is.null(object[[name.to]])) {
        object[[name.to]] <- paste0(object$path, name.to, '.fasta')
        msg <- paste0('Creating slot "', name.to, '"')
        message(msg)
    }
    ind <- file.exists(object[[name.from]])
    msg <- paste('Genome found for', sum(ind), 'out of', length(ind), 'samples.')
    message(msg)

    # Variables
    samples <- index(object)[ind]
    input.fasta <- object[[name.from]][ind]
    output.fasta <- object[[name.to]][ind]
    
    # Check output
    ind <- file.exists(output.fasta)
    missing <- if (recompute) samples else samples[!ind]
    missing <- if (length(recompute.sample)) recompute.sample else missing # Debug mode
    missing <- if (debug) head(missing, debug.n) else missing
    msg <- paste('Output found for', sum(ind), 'samples. Running Dnaapler for', length(missing), 'samples...')
    message(msg)
    
    # Main
    for (sample in samples) {
        i <- which(samples == sample)
        if (!sample %in% missing) next
        if (sample %in% samples.skip) next
        cat(sample)
        dnaapler_reorient(file = input.fasta[[i]], output.file = output.fasta[[i]])
    }
    
    # Exit
    return(object)
}