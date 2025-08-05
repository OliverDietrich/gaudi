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
#' @slot genome A character vector with paths to the genome FASTA.
#' @slot graph A character vector with paths to the genome assembly graph (GFA).
#' @slot log A character vector with paths to the assembly log file.
#'
methods::setClass(
    Class = "Assembly",
    slots = c(
        genome = "character",
        genome_graph = "character",
        log = "character"
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