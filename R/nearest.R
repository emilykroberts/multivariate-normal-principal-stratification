#' nearest
#'
#' @description nearest function from GenKern package (no longer available)
#'
#' @param x vector of values
#' @param xval value to find the nearest value in 'x' to
#'
#' @return index of 'x' with the value nearest to 'xval'
#'
#' @examples 
#' example(nearest(1:10, 3))
     
nearest = function(x, xval) # nearest function from Genkern
{
 
 na.rm=FALSE
 outside=FALSE
 
 # Do the NA handling
 # put it all together into a data frame or na.omit doesn't work
 x1 = x
 z = data.frame(x, x)
 # if NAs not allowed fail the function
 if(na.rm == FALSE){na.fail(z)}
 # get rid of NA cases
 if(na.rm == TRUE){z = na.omit(z)}
 # reassign the vectors with NAs removed
 x = z$x
 
 
 # if the value is outside the range of the vector and it isn't acceptable
 # then issue an error and stop
 if(outside == FALSE){if((max(x) <= xval) || (min(x) >= xval)) {stop("value outside vector range")}}
 
 # if the value is outside the range of the vector and this is acceptable then
 # merely assign one of the index with one of the extreme values in it
 
 if(outside == TRUE)
 {
 if((max(x) <= xval) || (min(x) >= xval))
 {
 sorx = sort(x)
 if(abs(sorx[1] - xval) < abs(sorx[length(sorx)]- xval))
 {index = 1; vally = sorx[index]; index = which(x1 == vally); return(index)}
 
 if(abs(sorx[1] - xval) > abs(sorx[length(sorx)]- xval))
 {index = length(sorx); vally = sorx[index]; index = which(x1 == vally); return(index)}
 }
 }
 
 
 # for most cases in which the value falls within the vector find the nearest
 # value and assign that index to the return value
 sorx = sort(x)
 upp = which(sorx >= xval)[1]
 low = which(sorx <= xval); low = low[length(low)]
 upp = sorx[upp]; low = sorx[low]
 
 if(upp == low) {index = which(x == upp)}
 if((abs(upp - xval)) >= (abs(low - xval))) {index = which(x1 == low)}
 if((abs(upp - xval)) < (abs(low - xval))) {index = which(x1 == upp)}
 
 return(index)
}
