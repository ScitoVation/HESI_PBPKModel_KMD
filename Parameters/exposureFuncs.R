#Create a perdose function to enable forcing
PerDose <- function(mag,Period,start,ExpDuration,times){
  #Number of replications of the event
  Nrep <- ceiling(max(times) / Period)
  #Find start and end times
  times <- rep(c(start, ExpDuration), Nrep) + rep(Period * (0:(Nrep - 1)), rep(2, Nrep))
  # Repeat start magnitude and 0 alternatingly to create a forcing vector
  y <- rep(c(mag,0), Nrep)
  cbind(times, y)
}