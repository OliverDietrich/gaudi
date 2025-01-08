#' Find NCBI database entries for a topic per year
#' 
#' @param db NCBI database to query
#' @param year Year to search entries for
#' @param term Search term
#' 
#' @export 
entrez_search_count_by_year <- function(db='nuccore', year=NULL, term=NULL){

    stopifnot(
        !is.null(year),
        !is.null(term)
    )
    
    # Run NCBI query
    query <- paste(term, "AND (", year, "[PDAT])")
    result <- entrez_search(db=db, term=query, retmax=0)$count

    return(result)
}