#Smoothing function to generate smoothed data
#xx is set of time that is 5x length of original data
#crit value is a way of adjusting the fit
fun_splsmooth <- function(as, PLATE, crit){
  Len <- length(PLATE[1,])
  xx <- seq(min(PLATE[[1]], na.rm = TRUE), max(PLATE[[1]], na.rm = TRUE), length.out=Len*5)
  absSpl <- smooth.spline(PLATE[[1]], as, df = max(xx)/crit)
  yy <- predict(absSpl, xx)
}