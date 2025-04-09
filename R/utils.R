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
    result <- rentrez::entrez_search(db=db, term=query, retmax=0)$count

    return(result)
}

#' Calculate plaque forming units (PFU) from spot assay
#'
#' Calculate PFU from plaque counts observed for a dilution of phage stock. Requires the spot volume.
#'
#' @param dil.factor Dilution factor of the droplet
#' @param plaque.count Number of plaques counted from the droplet
#' @param droplet.volume Volume of the deposited droplet (mL)
#' 
#' @export
calc_pfu <- function(
    dil.factor = NULL,
    plaque.count = 10,
    droplet.volume = 2.5e-3
) {
    stopifnot(
        !is.null(dil.factor),
        is.numeric(dil.factor),
        is.numeric(plaque.count),
        is.numeric(droplet.volume)
    )
    
    # Combine data
    object <- data.frame(
        'dilution' = dil.factor,
        'droplet_volume' = droplet.volume,
        'plaque_count' = plaque.count
    )
    
    # Correct errors
    
    ## Dilution 0 (should be 1)
    ind <- which(object$dilution == 0)
    object$dilution[ind] <- 1
    
    ## NAs
    ind <- which(is.na(object$dilution))
    object$dilution[ind] <- 1
    
    ## Dilution 1 (implies no plaques)
    ind <- which(object$dilution == 1)
    object$plaque_count[ind] <- 0
    
    # Calculate PFU
    object <- object %>% mutate(PFU = (plaque_count/droplet_volume) * (1/dilution))
    
    # Exit 0
    return(object[['PFU']])
}

#' Calculate efficiency of plating/plaquing (EOP) from plaque forming units (PFU)
#'
#' @param test.pfu Numeric, vector of PFU units
#' @param iso.pfu Numeric, vector of PFU units
#'
#' @export
calc_eop <- function(
    test.pfu = NULL,
    iso.pfu = NULL
) {
    
    # Minimal check
    stopifnot(
        !is.null(test.pfu),
        !is.null(iso.pfu)
    )
    
    # Check iso.pfu
    if (any(iso.pfu == 0)) {
        msg <- 'Isolation PFU must not be 0, aborting...'
        stop(msg)
    }
    
    # Combine data
    object <- data.frame(
        'iso' = iso.pfu,
        'test' = test.pfu
    )
    
    # Calculate EOP
    object <- object %>% mutate(eop = log10(test/iso))
    
    # Correct technical artifacts
    ind <- which(object$test == 0)
    if (length(ind) > 0) {
        msg <- 'Some EOP values have been assumed from isolation PFU.'
        warning(msg)
        object$eop[ind] <- -log10(object$iso[ind] * 0.003)
    }
    
    # Exit 0
    return(object[['eop']])
}

