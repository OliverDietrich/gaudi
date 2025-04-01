#' Minmax scaling
#' 
#' Rescales numeric vector to range a-b (default: 0-1)
#' Based on the formula {\displaystyle x'=a+{\frac {(x-{\text{min}}(x))(b-a)}{{\text{max}}(x)-{\text{min}}(x)}}}
#' 
#' @param x Numeric vector
#' @param a Lower limit
#' @param b Upper limit
#' 
#' @export
minmax <- function(x, a=0, b=1) {

  stopifnot(is.numeric(x))

  y <- a + ( (x-min(x)*(b-a) ) / (max(x) - min(x)) )

  return(y)
}