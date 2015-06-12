#' Square Function
#' 
#' This will squares the value of the numeric input
#' @author Ali Amiryousefi
#' @description This the others.
#' @param x is a numeric value
#' @examples 
#' x<-3
#' square(x) #returns the square value
#' @return The numeric \code{value} that is the result of the square input
#  manip misc # for manipulation and miscelineous
#' @seealso plot.apg
#' @export
#' 
square <- function(x){
  return(x^2)
}
#' Cube a number
#' @param x A numeric value to be cubed.
#' @return A numeric value that is the cube of \code{x}.
cube<- function(x) {
  return(x^3)
}