##################################################
############### Kelvin2Celsius ###################
##################################################
Kelvin2Celsius <- function(Kelvin) {
   # Converts temperature from Kelvin to Celsius
   #
   # Check and correct for negative values
   if (Kelvin < 0) {
      # If Kelvin is less than zero, set it to zero
      Kelvin <- 0
   }
   Celsius <- Kelvin - 273.15
   return(Celsius)
}
