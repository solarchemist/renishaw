# common.R
# General-purpose functions
# Taha Ahmed, Jan 2011

# CONTENTS
# >>>> ConvertRefPot
# >>>> Celsius2Kelvin
# >>>> Kelvin2Celsius
# >>>> as.radians
# >>>> as.degrees
# >>>> molarity2mass



##################################################
################# ConvertRefPot ##################
##################################################
ConvertRefPot <- function(argpotential, argrefscale, valuerefscale) {
   # Converts from some reference potential scale into another
   # SHE:     standard hydrogen electrode scale
   # Ag/AgCl: silver silver-chloride electrode scale
   # SCE:     standard calomel scale
   #
   ##### Add more reference electrodes here >>
   refpotatSHEzero <- c(     0,     -0.21,  -0.24,        3)
   refrownames     <- c( "SHE", "Ag/AgCl",  "SCE", "Li/Li+")
   refcolnames     <- c("SHE0",   "AgCl0", "SCE0",    "Li0")
   ##### Add more reference electrodes here <<
   #
   SHE0 <- data.frame(matrix(refpotatSHEzero, ncol=length(refpotatSHEzero), byrow=T))
   refpotmtx <- matrix(NA, length(SHE0), length(SHE0))
   refpotmtx[,1] <- matrix(as.matrix(SHE0), ncol=1, byrow=T)
   for (c in 2:length(SHE0)) {
      # loop over columns (except the first)
      for (r in 1:length(SHE0)) {
         # loop over rows
         refpotmtx[r, c] <- refpotmtx[r, 1] - refpotmtx[c, 1]
      }
   }
   refpotdf <- as.data.frame(refpotmtx)
   names(refpotdf) <- refcolnames
   row.names(refpotdf) <- refrownames
   ## So far we have made a matrix of all the possible combinations,
   ## given the vector refpotatSHEzero. The matrix is not strictly necessary,
   ## but it may prove useful later. It does.
   #
   # Match argrefscale to the refrownames
   argmatch <- match(argrefscale, refrownames, nomatch = 0)
   # Match valuerefscale to the refrownames
   valuematch <- match(valuerefscale, refrownames, nomatch = 0)
   # We simply assume that the match was well-behaved
   valuepotential <- argpotential + refpotdf[valuematch, argmatch]
   # Check that arg and value electrodes are within bounds for a match
   if (argmatch == 0 || valuematch == 0) {
      # No match
      # Perform suitable action
      message("Arg out of bounds in call to ConvertRefPot")
      valuepotential <- NA
   }
   return(valuepotential)
}



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
