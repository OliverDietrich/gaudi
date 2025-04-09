#' Create diverging colorscale
#' 
#' Creates a diverging colorscale around a midpoint
#' based on RColorBrewer palettes
#' 
#' @param x Numeric vector or matrix
#' @param midpoint Integer within data range
#' @param breaks_per_side Integer specifying the number of color values on either side of the midpoint
#' @param palette Name of the RColorBrewer palette
#' @param direction Direction of the color values (1, -1)
#' 
#' @export
colorscale_diverging <- function(x = NULL, 
                                 midpoint=0, 
                                 breaks_per_side = 50, 
                                 palette='RdYlBu', 
                                 direction=-1
                                ) {
    
    stopifnot(
        !is.null(x),
        is.numeric(x),
        direction %in% c(1,-1)
    )
    
    # Remove NAs
    if (any(is.na(x))) {
        x <- na.omit(x)
    }
    
    # Errors
    ## Data range
    if (x < min(x) || x > max(x)) {
        stop('Midpoint is out of data range')
    }
    ## Color palette
    if (!palette %in% row.names(RColorBrewer::brewer.pal.info)) {
        stop('Palette not available in RColorBrewer')
    }
    
    # Create breaks
    breaks <- unique(c(
        seq(min(x), midpoint, length.out = breaks_per_side),
        seq(midpoint, max(x), length.out = breaks_per_side)
    ))
    
    # Color values
    color_values <- RColorBrewer::brewer.pal(n = 9, name = palette)
    if (direction == -1) {
        color_values <- rev(color_values)
    }
    
    # Create colorscale
    colorscale <- colorRampPalette(color_values)(length(breaks))
    
    # Add breaks
    names(colorscale) <- breaks
    
    return(colorscale)
}