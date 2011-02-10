# common.R
# General-purpose functions
# Taha Ahmed, Jan 2011

# CONTENTS
# >>>> Celsius2Kelvin
# >>>> Kelvin2Celsius
# >>>> as.radians
# >>>> as.degrees
# >>>> molarity2mass




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
}


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
}


##################################################
################# as.radians #####################
##################################################
as.radians <- function(degrees) {
   # Converts from degrees to radians
   radians <- degrees * (pi / 180)
}


##################################################
################# as.degrees #####################
##################################################
as.degrees <- function(radians) {
   # Converts from radians to degrees
   radians <- radians * (180 / pi)
}


##################################################
############### molarity2mass ####################
##################################################
molarity2mass <- function(formulamass, volume, molarity) {
   # Calculates the required mass of
   # the substance to be dissolved.
   # ARGS: formulamass - formula mass of the substance (in gram per mole)
   #       volume      - volume of the final solution (in liters)
   #       molarity    - molarity (in moles per liter)
   # VALUE: mass of substance (in grams)
   #
   mass <- formulamass * volume * molarity 
   # Unit check:
   # [g * mol-1] * [liter] * [mole * liter-1] = [g]
   return(mass)
}
