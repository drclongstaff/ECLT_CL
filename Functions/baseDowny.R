#zero all the curves
fun_baseline <- function(n, off) {
  n <- n-min(n)-off
}
#Calculate the time to percent abs change
fun_Downy <- function(Y, aPlate, PerC){
  Time <- aPlate[[1]]
  maxY <- max(Y, na.rm = TRUE)
  minY <- min(Y, na.rm = TRUE)
  pointmax <- which.max(Y)
  pcChange<-0.01*PerC*(maxY-minY)+minY
  downTime <- Time[-c(1:pointmax)]
  downAbs <- Y[-c(1:pointmax)]
  #This deals with wiggly late points
  TC <- which.min(downAbs)
  #downTime <- downTime[1:TC]
  #downAbs <- downAbs[1:TC]
  ifelse(TC<=1, downTime <- downTime, downTime <- downTime[1:TC])
  ifelse(TC<=1, downAbs <- downAbs, downAbs <- downAbs[1:TC])
  #Try replacement calculations
  #decayPoint<-which(abs(Y-pcChange)==min(abs(Y-pcChange)))[1]#Find closest point
  decayPoint<-which(abs(downAbs-pcChange)==min(abs(downAbs-pcChange)))[1] 
  #Check if decaypoint is close to the end
  ifelse(decayPoint>=length(Y)-3, decayPoint <- 1, decayPoint<-decayPoint)
  #Check if final abs is > pcChange
  ifelse(Y[length(Y)]>=pcChange, decayTime <- 0, decayTime<-Time[decayPoint+pointmax])
}
