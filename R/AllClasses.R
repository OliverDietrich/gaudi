#' The virusCollection Class
#'
#' The virusCollection object contains viral sequence data. Formatted based on NCBI datasets.
#'
#' @slot genomes Data.frame containing genome metadata
#' @slot features Data.frame containing genomic features (based on gff3)
#' @slot key Character vector indicating the column to use as key between layout and data
#' @slot misc A list of miscellaneous information
#' 
#' @exportClass virusCollection
#'
methods::setClass(
  Class = "virusCollection",
  slots = c(
    layout = "data.frame",
    data = "data.frame",
    key = "character",
    misc = "list"
  )
)