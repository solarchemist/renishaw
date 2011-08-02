##################################################
################### It2charge ####################
##################################################
It2charge <- function (time, current) {
   ## **** STOP USING THIS FUNCTION *** CAUSED WEIRD, UNREPRODUCIBLE ERRORS /110304
   ## Description:
   ##    Calculates cumulative charge, differentials, etc. from
   ##    amperometric data (current and time).
   ##    __Intended to be used from within other functions (CHI.R)__
   ## Usage:
   ##    It2charge(time, current)
   ## Arguments:
   ##       time: a vector with time data.
   ##    current: a vector of the same length as time, with current data.
   ##             May be either currents or current densities, no matter.
   ## Value:
   ##    Returns a dataframe with columns:
   ##    timediff, dIdt, charge, sumcharge
   #
   # Calculate the time vector difference
   timediff <- c(time[1], diff(time))
   # timediff times the current gives the charge,
   # since the time vector can be considered as
   # the cumulative time, while we need to multiply
   # the current with the time elapsed since the last 
   # current measurement (the timediff).
   charge <- current * timediff
   dIdt <- current / time
   # Return value
   ff <- data.frame(timediff = timediff, 
      dIdt = dIdt, charge = charge, 
      # perhaps it is more correct to calculate cumsum of the absolute charge?
      sumcharge = cumsum(charge))
   return(ff)
}