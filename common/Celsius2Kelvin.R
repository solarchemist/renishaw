##################################################
############### Celsius2Kelvin ###################
##################################################
Celsius2Kelvin <- function(Celsius) {
   # Converts temperature from Celsius to Kelvin
   #
   # Check and correct for values below -273.15
   if (Celsius < -273.15) {
      # If Celsis is less than absolute zero, set it to absolute zero
      Celsius <- -273.15
   }
   Kelvin <- Celsius + 273.15
   return(Kelvin)
}