#' The genomeCollection Class
#'
#' The genomeCollection object contains genome sequence data. Formatted based 
#' on NCBI datasets.
#'
#' @slot index Character, unique identifiers for each genome
#' @slot path Character, path to the directory where progress is tracked
#' @slot meta.data Data.frame containing genome metadata
#' @slot assembly A list of assembly objects
#' @slot annotation A list of annotation objects 
#' @slot distances A list of distance objects
#' @slot misc A list of miscellaneous information
#' 
#' @exportClass genomeCollection
#'
methods::setClass(
  Class = "genomeCollection",
  slots = c(
    index = "character",
    path = "character",
    meta.data = "data.frame",
    reads = "list",
    assembly = "list",
    annotation = "list",
    distances = "list",
    misc = "list"
  )
)

#' The Reads class
#'
#' The Read object stores file paths for raw and filtered FASTQ reads. Reads 
#' can be either paired-end short reads (R1 + R2), single short reads (S) or
#' long reads (L). While samples can have multiple FASTQ files (e.g. hybrid
#' assemblies with L and R1+R2) each row corresponds to a sample. For multiple
#' sequencing runs of the same type (2x long reads) the FASTQ files must be 
#' concatenated before being assigned to the Reads object or registered in
#' two distinct Reads object before being merged into a third Reads object.
#'
#' @slot index Character vector matching to index(object)
#' @slot alias Character vector of unique aliases for duplicated index(object)
#' @slot R1 Character, file paths with forward paired-end reads
#' @slot R2 Character, file paths with reverse paired-end reads
#' @slot S Character, file paths with unpaired reads
#' @slot L Character, file paths with long reads
#' @slot counts Data.frame with per sample summary of read counts
#' @slot quality Data.frame with per FASTQ quality summary (default 5000 reads per FASTQ file)
#' @slot type Summary what type of reads are present in the object
#'
methods::setClass(
    Class = "Reads",
    slots = c(
        index = "character",
        alias = "character",
        R1 = "character",
        R2 = "character",
        S = "character",
        L = "character",
        counts = "data.frame",
        quality = "data.frame",
        type = "character"
    )
)

#' The Assembly class
#'
#' The Assembly object stores file paths for genome assemblies. The core components are
#' the genome FASTA file (necessary) and the assembly graph (optional).
#'
#' @slot index A character vector sample indices matching index(genomeCollection)
#' @slot genome A character vector with paths to the genome FASTA.
#' @slot graph A character vector with paths to the genome assembly graph (GFA).
#' @slot log A character vector with paths to the assembly log file.
#' @slot reads A character vector indicating the slot name in ReadsNames(object).
#' @slot tool A character vector indicating the assembly tool (e.g. SPAdes, Flye, Unicycler, ...)
#' @slot type A character vector indicating the assembly type
#'
methods::setClass(
    Class = "Assembly",
    slots = c(
        index = "character",
        contig = "character",
        graph = "character",
        log = "character",
        reads = "character",
        tool = "character",
        type = "character"
    )
)

#' The Annotation class
#'
#' The Annotation object stores file paths for genome annotation The core components are
#' the genes GFF file, the coding sequences (CDS) FASTA, the protein FASTA.
#'
#' @slot index A character vector sample indices matching index(genomeCollection)
#' @slot genes A character vector with paths to the genes GFF.
#' @slot cds A character vector with paths to the coding sequences FASTA.
#' @slot proteins A character vector with paths to the protein FASTA.
#' @slot log A character vector with paths to the annotation log file.
#' @slot genome A character vector indicating the column in metadata(genomeCollection) storing genome locations.
#' @slot tool A character vector indicating the annotation tool (e.g. Bakta, Padloc, DefenseFinder, ...)
#' @slot type A character vector indicating the annotation type
#'
methods::setClass(
    Class = "Annotation",
    slots = c(
        index = "character",
        genes = "character",
        cds = "character",
        proteins = "character",
        log = "character",
        genome = "character",
        tool = "character",
        type = "character"
    )
)

#' The Distance class
#'
#' The Distance object stores file paths for genome distances to infer phylogeny.
#' The core components are ...
#'
#' @slot index A character vector sample indices matching index(genomeCollection)
#' @slot distance ...
#' @slot sketches ...
#' @slot log A character vector with paths to the annotation log file.
#' @slot tool A character vector indicating the distance tool (e.g. Mash, FastANI, skani, MMseqs, MUMer, Foldseek, ...)
#' @slot type A character vector indicating the distance type (e.g. genome, CDS, protein, ...)
#'
methods::setClass(
    Class = "Distance",
    slots = c(
        index = "character",
        distance = "character",
        sketches = "character",
        log = "character",
        sequences = "character",
        tool = "character",
        type = "character"
    )
)

#' The interactionSet Class
#'
#' The interactionSet object contains an interaction matrix from two genomeCollections.
#'
#' @slot interactions Matrix of interactions
#' @slot phages genomeCollection of bacteriophages
#' @slot bacteria genomeCollection of bacteria
#'
methods::setClass(
    Class = "interactionSet",
    slots = c(
        interactions = "matrix",
        phages = "genomeCollection",
        bacteria = "genomeCollection"
    )
)