#' Find NCBI database entries for a topic per year
entrez_search_count_by_year <- function(db='pubmed', year=NULL, term=NULL){

    stopifnot(
        !is.null(year),
        !is.null(term)
    )
    
    # Run NCBI query
    query <- paste(term, "AND (", year, "[PDAT])")
    result <- entrez_search(db=db, term=query, retmax=0)$count

    return(result)
}