#zero all the curves
fun_baseline <- function(n, off) {
  n <- n-min(n)-off
}
#Calculate the time to percent abs change
fun_Downy <- function(Y, PLATE, PerC){
  Time <- PLATE[[1]]
  maxY <- max(Y, na.rm = TRUE)
  minY <- min(Y, na.rm = TRUE)
  pcChange<-0.01*PerC*(maxY-minY)+minY
  decayPoint<-which(abs(Y-pcChange)==min(abs(Y-pcChange)))[1]#Find closest point
  #Check if decaypoint is close to the end
  ifelse(decayPoint>=length(Y)-3, decayPoint <- 1, decayPoint<-decayPoint)
  #Check if final abs is > pcChange
  ifelse(Y[length(Y)]>=pcChange, decayTime <- 0, decayTime<-Time[decayPoint])
}
