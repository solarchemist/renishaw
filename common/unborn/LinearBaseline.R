##################################################
################ LinearBaseline ##################
##################################################
LinearBaseline <- function (potential, current, iplim) {
   ## Arguments:
   ##   potential: full half-cycle, potentials
   ##     current: full half-cycle, currents
   ##       iplim: interpolation limits along x (potential)
   ## Value:
   ##   A dataframe with two columns, potential and current
   ##   (of the calculated baseline)
   # Construct potential-current dataframe
   sweep <- data.frame(potential = potential, current = current)
   #
   sweep.iplim <- subset(subset(sweep, potential > iplim[1]), potential < iplim[2])
   sweep.baseline <- data.frame(potential = approxExtrap(sweep.iplim$potential,
      sweep.iplim$current, xout = sweep$potential, method = "linear")$x,
      current = approxExtrap(sweep.iplim$potential, sweep.iplim$current, 
      xout = sweep$potential, method = "linear")$y)
   sweep.data <- data.frame(potential = sweep.baseline$potential,
      current = sweep.baseline$current)
   return(sweep.data)
}