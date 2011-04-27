# common.R
# General-purpose functions
# Taha Ahmed, Jan 2011

# CONTENTS              Status      Depends on
# --------              ------      ----------
# >>>> LinearBaseline   deprecated  ?
# >>>> int2padstr                   ?
# >>>> It2charge        deprecated  ?
# >>>> ProvideSampleId              ?
# >>>> Celsius2Kelvin               ?
# >>>> Kelvin2Celsius               ?
# >>>> as.radians                   ?
# >>>> as.degrees                   ?
# >>>> molarity2mass                ?
# >>>> AVS2SHE                      -
# >>>> SHE2AVS                      -
# >>>> ConvertRefPotEC              -
# >>>> ConvertRefPot                ConvertRefPotEC, SHE2AVS, AVS2SHE



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




##################################################
################## int2padstr ####################
##################################################
int2padstr <- function (ii, pchr, w) {
   ## Description:
   ##   Converts an integer or a vector of integers to
   ##   a string padded with characters. 
   ## Usage:
   ##   int2padstr(ii, pchr, w)
   ## Arguments:
   ##     ii: integer or vector of integers
   ##   pchr: a padding character (e.g., "0")
   ##      w: width of the return string (an integer)
   ##         Make sure to set the width longer than
   ##         or equal to the length of the biggest integer.
   ##         For example, if the integers (ii) are
   ##         in the range 1 - 100, set w to at least 3.
   ## Value:
   ##   A character string or a vector of character strings
   gsub(" ", pchr, formatC(ii, format="s", mode="character", width = w))
}


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



##################################################
################ ProvideSampleId #################
##################################################
ProvideSampleId <- function (fullpathwithfilename) {
     ### OBS! Only very rudimentary error-checking.
     ### If the filename is formatted as \w*-\w*-\w*, we use the middle segment, 
     ### otherwise we use the whole string (excluding the extension)
     # Extract the name of the parent directory of the datafilename argument
     substrateid <- basename(dirname(fullpathwithfilename))
     # Extract the name of the method from the filename-part
     # First split the filename over all hyphens
     nameparts <- strsplit(basename(fullpathwithfilename), "-")[[1]]
     # If the number of nameparts exceed 3, save the whole filename as methodid, otherwise use the middle part
     if (length(nameparts) > 3) {
        # We need to lose the file extension from the last namepart
        nameparts[length(nameparts)] <- strsplit(nameparts[length(nameparts)], "\\.")[[1]][1]
        methodid <- paste(nameparts, collapse = "-")
     } else {
        methodid <- nameparts[2]
     }
     # Make an informative sampleid
     sampleid <- paste(substrateid, methodid, sep = "-")
     #
     return(sampleid)
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
   return(Kelvin)
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
   return(Celsius)
}


##################################################
################# as.radians #####################
##################################################
as.radians <- function(degrees) {
   # Converts from degrees to radians
   radians <- degrees * (pi / 180)
   return(radians)
}



##################################################
################# as.degrees #####################
##################################################
as.degrees <- function(radians) {
   # Converts from radians to degrees
   degrees <- radians * (180 / pi)
   return(degrees)
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



##################################################
#################### AVS2SHE #####################
##################################################
AVS2SHE <- function(avs) {
   # Converts from absolute vacuum scale (AVS) to SHE scale
   she <- -(4.5 + avs)
   return(she)
}



##################################################
#################### SHE2AVS #####################
##################################################
SHE2AVS <- function(she) {
   # Converts from SHE scale to absolute vacuum (AVS) scale
   avs <- -(4.5 + she)
   return(avs)
}



##################################################
############### ConvertRefPotEC ##################
##################################################
ConvertRefPotEC <- function(argpotential, argrefscale, valuerefscale) {
   # Converts from an electrochemical reference potential scale into another
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
################# ConvertRefPot ##################
##################################################
ConvertRefPot <- function(argpotential, argrefscale, valuerefscale) {
   # Check that argpotential is valid numeric 
   
   #  IDEA: make a matrix out of these (scale names and flags)
   
   # Valid scales
   scale.names <- list()
   scale.names[["SHE"]] <- c("SHE", "NHE", "she", "nhe")
   scale.names[["AgCl"]] <- c("Ag/AgCl", "AgCl", "ag/agcl", "agcl")
   scale.names[["SCE"]] <- c("SCE", "sce")
   scale.names[["Li"]] <- c("Li/Li+", "Li", "Li+", "li", "li+", "li/li+")
   scale.names[["AVS"]] <- c("AVS", "avs")
   
   # Set flags
   bool.flags <- as.data.frame(matrix(0, nrow = length(scale.names), ncol = 2))
   names(bool.flags) <- c("argref", "valueref")
   row.names(bool.flags) <- names(scale.names)
   
   # argrefscale
   # Check that argrefscale is valid character mode
   # ...
   
   for (j in 1:length(row.names(bool.flags))) {
      if (any(scale.names[[row.names(bool.flags)[j]]] == argrefscale)) {
         bool.flags[row.names(bool.flags)[j], "argref"] <- j
      }
   }
   
   
   # valuerefscale
   # Check that valuerefscale is valid character mode
   # ...
   
   for (k in 1:length(row.names(bool.flags))) {
      if (any(scale.names[[row.names(bool.flags)[k]]] == valuerefscale)) {
         bool.flags[row.names(bool.flags)[k], "valueref"] <- k
      }
   }
   
   # Depending on which flags are set, call the corresponding function
   
   decision.vector <- colSums(bool.flags)
   
   # Check if both scales are the same (no conversion needed). If so, abort gracefully.
   # ...
   
   if (decision.vector["argref"] == 5 || decision.vector["valueref"] == 5) {
      # AVS is requested, deal with it it
      if (decision.vector["argref"] == 5) {
         # Conversion _from_ AVS
         rnpotential <- ConvertRefPotEC(AVS2SHE(argpotential), 
            "SHE", 
            scale.names[[decision.vector["valueref"]]][1])
      } 
      if (decision.vector["valueref"] == 5) {
         # Conversion _to_ AVS
         rnpotential <- SHE2AVS(ConvertRefPotEC(argpotential, 
            scale.names[[decision.vector["argref"]]][1], 
            "SHE"))
      }
   } else {
      rnpotential <- ConvertRefPotEC(argpotential, 
         scale.names[[decision.vector["argref"]]][1], 
         scale.names[[decision.vector["valueref"]]][1])
   }
   return(rnpotential)
}