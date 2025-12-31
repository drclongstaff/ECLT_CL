#Smoothing function to generate smoothed data
#xx is set of time that is 5x length of original data
#crit value is a way of adjusting the fit
fun_splsmooth <- function(as, aPlate, crit){
  Len <- length(aPlate[1,])
  #increase number of points by 5x
  xx <- seq(min(aPlate[[1]], na.rm = TRUE), max(aPlate[[1]], na.rm = TRUE), length.out=Len*5)
  absSpl <- smooth.spline(aPlate[[1]], as, df = max(xx)/crit)
  yy <- predict(absSpl, xx)
}