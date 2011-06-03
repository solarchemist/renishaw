# CHI.R
# Functions to read and manipulate data from the CHI760 potentiostat/galvanostat
# Taha Ahmed, Jan 2011 - Feb 2011

# CONTENTS
source("/home/taha/chepec/chetex/common/R/common.R")
# >>>> mps2df
# >>>> ocp2df
# >>>> chronocm2df
# >>>> chronoamp2df
# >>>> amperometry2df
# >>>> cv2df 
# >>>> lsv2df



##################################################
################### mps2df #######################
##################################################
mps2df <- function(datafilename, wearea = 1) {
   ## Description:
   ##   Reads time vs current data (from CHI 760 potentiostat)
   ##   and returns a dataframe with the data and 
   ##   the data attributes (experimental conditions).
   ## Usage:
   ##   mps2df(datafilename, wearea)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ##         wearea: working electrode area in square centimeters (optional)
   ## Value:
   ##   Dataframe with the following columns (and no extra attributes):
   ##   $ sampleid        : chr
   ##   $ time            : num
   ##   $ current         : num
   ##   $ currentdensity  : num
   ##   $ timediff        : num
   ##   $ dIdt            : num
   ##   $ didt            : num
   ##   $ charge          : num
   ##   $ chargedensity   : num
   ##   $ PotentialSteps  : num
   ##   $ TimeSteps       : num
   ##   $ Cycle           : num
   ##   $ SampleInterval  : num
   ##   $ QuietTime       : num
   ##   $ Sensitivity     : num
   ## Note:
   ##   The CH Instruments 760 potentiostat records all data 
   ##   using standard SI units, therefore this function
   ##   assumes all potential values to be in volts, 
   ##   currents to be in amperes, charges in Coulombs, 
   ##   time in seconds, and so on.
   #
   datafile <- file(datafilename, "r")
   chifile <- readLines(datafile, n = -1) #read all lines of input file
   close(datafile)
   #
   sampleid <- ProvideSampleId(datafilename)
   #
   rgxp.number <- "^\\-?\\d+\\.\\d+[e,]"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers with or without a negative sign (hyphen), 
   # followed by one digit before the decimal, a decimal point,
   # and an arbitrary number of digits after the decimal point,
   # immediately followed by either the letter 'e' or a comma.
   # Note that backslashes are escaped due to the way R handles strings.
   #
   numrow.idx <- regexpr(rgxp.number, chifile)
   # Save the match length attribute to another variable,
   numrow.len <- attr(numrow.idx, "match.length")
   # then scrap the attribute of the original variable.
   attr(numrow.idx, "match.length") <- NULL
   #
   i <- seq(1, length(numrow.idx) - 1, 1)
   j <- seq(2, length(numrow.idx), 1)
   # Start indices of data ranges
   starts <- which(numrow.idx[i] != 1 & numrow.idx[j] == 1) + 1
   # End indices, except for the last
   ends <- which(numrow.idx[i] == 1 & numrow.idx[j] != 1)
   # Fix the last index of end indices
   ends <- c(ends, length(numrow.idx))
   #
   ff <- data.frame(NULL)
   for (s in 1:length(starts)) {
      zz <- textConnection(chifile[starts[s]:ends[s]], "r")
      ff <- rbind(ff,
               data.frame(stringsAsFactors = FALSE,
               sampleid, matrix(scan(zz, what = numeric(), sep = ","),
               ncol = 2, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "time", "current")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff, currentdensity = currentdensity)
   # Calculate time diff and current diff
   timediff <- c(ff$time[1], diff(ff$time))
   currentdiff <- c(ff$current[1], diff(ff$current))
   currentdensitydiff <- c(ff$currentdensity[1], diff(ff$currentdensity))
   # Calculate differential of current and current density
   dIdt <- currentdiff / timediff
   didt <- currentdensitydiff / timediff
   # Calculate charge and charge density
   charge <- cumsum(ff$current)
   chargedensity <- cumsum(ff$currentdensity)
   # Update ff dataframe
   ff <- cbind(ff, 
      timediff = timediff, 
      dIdt = dIdt, 
      didt = didt,
      charge = charge, 
      chargedensity = chargedensity)
   #
   ### Collect attributes of this experiment
   # Potential steps
   PotentialSteps <- ""
   positions.PotentialSteps <- regexpr("^E\\d+\\s\\(V\\)", chifile)
   PotentialSteps <- paste("\\num{", strsplit(chifile[which(positions.PotentialSteps == 1)], 
      "\\s=\\s")[[1]][2], "}", sep = "")
   if (length(which(positions.PotentialSteps == 1)) > 1) {
      for (i in 2:length(which(positions.PotentialSteps == 1))) {
      PotentialSteps <- 
         paste(PotentialSteps, 
            paste("\\num{", strsplit(chifile[which(positions.PotentialSteps == 1)], 
            "\\s=\\s")[[i]][2], "}", sep = ""), sep = ", ")
      }
   }
   ff$PotentialSteps <- PotentialSteps
   # Time steps
   TimeSteps <- ""
   positions.TimeSteps <- regexpr("^T\\d+\\s\\(s\\)", chifile)
   TimeSteps <- paste("\\num{", strsplit(chifile[which(positions.TimeSteps == 1)], 
      "\\s=\\s")[[1]][2], "}", sep = "")
   if (length(which(positions.TimeSteps == 1)) > 1) {
      for (i in 2:length(which(positions.TimeSteps == 1))) {
      TimeSteps <- 
         paste(TimeSteps, 
            paste("\\num{", strsplit(chifile[which(positions.TimeSteps == 1)], 
            "\\s=\\s")[[i]][2], "}", sep = ""), sep = ", ")
      }
   }
   ff$TimeSteps <- TimeSteps
   # Cycles
   position.Cycles <- regexpr("^Cycle", chifile)
   Cycles <- as.numeric(strsplit(chifile[which(position.Cycles == 1)], "\\s=\\s")[[1]][2])
   ff$Cycles <- Cycles
   # Sample interval (s)
   position.SampleInterval <- regexpr("^Sample\\sInterval\\s\\(s\\)", chifile)
   SampleInterval <- as.numeric(strsplit(chifile[which(position.SampleInterval == 1)], "\\s=\\s")[[1]][2])
   ff$SampleInterval <- SampleInterval
   # Quiet Time (s)
   position.QuietTime <- regexpr("^Quiet\\sTime\\s\\(sec\\)", chifile)
   QuietTime <- as.numeric(strsplit(chifile[which(position.QuietTime == 1)], "\\s=\\s")[[1]][2])
   ff$QuietTime <- QuietTime
   # Sensitivity (A/V)
   position.Sensitivity <- regexpr("^Sensitivity\\s\\(A/V\\)", chifile)
   Sensitivity <- as.numeric(strsplit(chifile[which(position.Sensitivity == 1)], "\\s=\\s")[[1]][2])
   ff$Sensitivity <- Sensitivity
   #
   return(ff)
}




##################################################
################### ocp2df #######################
##################################################
ocp2df <- function(datafilename) {
   ## Description:
   ##   Reads time vs potential data (from CHI 760 potentiostat)
   ##   and returns a dataframe with the data and 
   ##   the data attributes (experimental conditions).
   ## Usage:
   ##   ocp2df(datafilename)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ## Value:
   ##   Dataframe with the following columns (and no extra attributes):
   ##   $ sampleid        : chr
   ##   $ time            : num
   ##   $ current         : num
   ##   $ currentdensity  : num
   ##   $ timediff        : num
   ##   $ dIdt            : num
   ##   $ didt            : num
   ##   $ charge          : num
   ##   $ chargedensity   : num
   ##   $ InitE           : num
   ##   $ SampleInterval  : num
   ##   $ RunTime         : num
   ##   $ QuietTime       : num
   ##   $ Sensitivity     : num
   ## Note:
   ##   The CH Instruments 760 potentiostat records all data 
   ##   using standard SI units, therefore this function
   ##   assumes all potential values to be in volts, 
   ##   currents to be in amperes, charges in Coulombs, 
   ##   time in seconds, and so on.
   #
   datafile <- file(datafilename, "r")
   chifile <- readLines(datafile, n = -1) #read all lines of input file
   close(datafile)
   #
   sampleid <- ProvideSampleId(datafilename)
   #
   rgxp.number <- "^\\-?\\d\\.\\d+[e,]"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers with or without a negative sign (hyphen), 
   # followed by one digit before the decimal, a decimal point,
   # and an arbitrary number of digits after the decimal point,
   # immediately followed by either the letter 'e' or a comma.
   # Note that backslashes are escaped due to the way R handles strings.
   #
   numrow.idx <- regexpr(rgxp.number, chifile)
   # Save the match length attribute to another variable,
   numrow.len <- attr(numrow.idx, "match.length")
   # then scrap the attribute of the original variable.
   attr(numrow.idx, "match.length") <- NULL
   #
   i <- seq(1, length(numrow.idx) - 1, 1)
   j <- seq(2, length(numrow.idx), 1)
   # Start indices of data ranges
   starts <- which(numrow.idx[i] != 1 & numrow.idx[j] == 1) + 1
   # End indices, except for the last
   ends <- which(numrow.idx[i] == 1 & numrow.idx[j] != 1)
   # Fix the last index of end indices
   ends <- c(ends, length(numrow.idx))
   #
   ff <- data.frame(NULL)
   for (s in 1:length(starts)) {
      zz <- textConnection(chifile[starts[s]:ends[s]], "r")
      ff <- rbind(ff,
               data.frame(stringsAsFactors = FALSE,
               sampleid, matrix(scan(zz, what = numeric(), sep = ","),
               ncol = 2, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "time", "potential")
   #
   ### Collect attributes of this experiment
   # RunTime (sec)
   position.RunTime <- regexpr("^Run\\sTime\\s\\(sec\\)", chifile)
   RunTime <- as.numeric(strsplit(chifile[which(position.RunTime == 1)], "\\s=\\s")[[1]][2])
   ff$RunTime <- RunTime
   #
   return(ff)
}



##################################################
################# chronocm2df ####################
##################################################
chronocm2df <- function(datafilename) {
   # Function description: chronocoulometry data
   # CH Instruments potentiostat records all data using standard SI units,
   # so all potential values are in volts, currents are in amperes,
   # charges in Coulombs, time in seconds, etc.
   #
   datafile <- file(datafilename, "r")
   chifile <- readLines(datafile, n = -1) #read all lines of input file
   close(datafile)
   #
   sampleid <- ProvideSampleId(datafilename)
   #
   rgxp.number <- "^\\-?\\d\\.\\d+[e,]"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers with or without a negative sign (hyphen), 
   # followed by one digit before the decimal, a decimal point,
   # and an arbitrary number of digits after the decimal point,
   # immediately followed by either the letter 'e' or a comma.
   # Note that backslashes are escaped.
   #
   numrow.idx <- regexpr(rgxp.number, chifile)
   # Save the match length attribute to another variable,
   numrow.len <- attr(numrow.idx, "match.length")
   # then scrap the attribute of the original variable.
   attr(numrow.idx, "match.length") <- NULL
   #
   i <- seq(1, length(numrow.idx) - 1, 1)
   j <- seq(2, length(numrow.idx), 1)
   # Start indices of data ranges
   starts <- which(numrow.idx[i] != 1 & numrow.idx[j] == 1) + 1
   # End indices, except for the last
   ends <- which(numrow.idx[i] == 1 & numrow.idx[j] != 1)
   # Fix the last index of end indices
   ends <- c(ends, length(numrow.idx))
   #
   ff <- data.frame(NULL)
   for (s in 1:length(starts)) {
      zz <- textConnection(chifile[starts[s]:ends[s]], "r")
      ff <- rbind(ff,
               data.frame(sampleid, step = factor(s),
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 2, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "step", "time", "charge")
   #
   ### Collect attributes of this experiment
   # These attributes are specific for each kind of experiment,
   # be careful when adapting to other electrochemical data
   rgxp.attr <- c("^Init\\sE\\s\\(V\\)",
                  "^Final\\sE\\s\\(V\\)",
                  "^Step\\s",
                  "^Pulse\\sWidth\\s\\(sec\\)",
                  "^Sample\\sInterval\\s\\(s\\)",
                  "^Quiet\\sTime\\s\\(sec\\)",
                  "^Sensitivity\\s\\(A/V\\)")
   names.attr <- c("InitE",
                   "FinalE",
                   "Steps",
                   "PulseWidth",
                   "SamplingInterval",
                   "QuietTime",
                   "Sensitivity")
   for (n in 1:length(rgxp.attr)) {
      attrow.idx <- regexpr(rgxp.attr[n], chifile)
      attrow.len <- attr(attrow.idx, "match.length")
      attr(attrow.idx, "match.length") <- NULL
      attr(ff, names.attr[n]) <- strsplit(chifile[which(attrow.idx == 1)],
         "\\s=\\s")[[1]][2]
   }
   #
   return(ff)
}



##################################################
################# chronoamp2df ###################
##################################################
chronoamp2df <- function(datafilename, wearea = 1) {
   ## Description:
   ##   Reads current-time data (from CHI 760 potentiostat)
   ##   and returns a dataframe with the data and the 
   ##   data attributes (experimental conditions).
   ## Usage:
   ##   chronoamp2df(datafilename, wearea)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ##         wearea: (optional) area of working electrode (in square centimeters)
   ## Value:
   ##   Dataframe with the following columns (and no extra attributes):
   ##   $ sampleid        : chr
   ##   $ step            : num
   ##   $ time            : num
   ##   $ current         : num
   ##   $ currentdensity  : num
   ##   $ InitE           : num
   ##   $ HighE           : num
   ##   $ LowE            : num
   ##   $ InitPN          : chr
   ##   $ Step            : num
   ##   $ Pulse width     : num
   ##   $ Sample interval : num
   ##   $ Quiet Time      : num
   ##   $ Sensitivity     : num
   ## Note:
   ##   The CH Instruments 760 potentiostat records all data 
   ##   using standard SI units, therefore this function
   ##   assumes all potential values to be in volts, 
   ##   currents to be in amperes, charges in Coulombs, 
   ##   time in seconds, and so on.
   #
   datafile <- file(datafilename, "r")
   chifile <- readLines(datafile, n = -1) #read all lines of input file
   close(datafile)
   #
   sampleid <- ProvideSampleId(datafilename)
   #
   rgxp.number <- "^\\-?\\d\\.\\d+[e,]"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers with or without a negative sign (hyphen), 
   # followed by one digit before the decimal, a decimal point,
   # and an arbitrary number of digits after the decimal point,
   # immediately followed by either the letter 'e' or a comma.
   # Note that backslashes are escaped.
   #
   numrow.idx <- regexpr(rgxp.number, chifile)
   # Save the match length attribute to another variable,
   numrow.len <- attr(numrow.idx, "match.length")
   # then scrap the attribute of the original variable.
   attr(numrow.idx, "match.length") <- NULL
   #
   i <- seq(1, length(numrow.idx) - 1, 1)
   j <- seq(2, length(numrow.idx), 1)
   # Start indices of data ranges
   starts <- which(numrow.idx[i] != 1 & numrow.idx[j] == 1) + 1
   # End indices, except for the last
   ends <- which(numrow.idx[i] == 1 & numrow.idx[j] != 1)
   # Fix the last index of end indices
   ends <- c(ends, length(numrow.idx))
   #
   ff <- data.frame(NULL)
   for (s in 1:length(starts)) {
      zz <- textConnection(chifile[starts[s]:ends[s]], "r")
      ff <- rbind(ff,
               data.frame(sampleid, step.current = s,
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 2, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "step.current", "time", "current")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff, currentdensity = currentdensity)
   #
   
   ### Collect attributes of this experiment
   # InitE (volt)
   position.InitE <- regexpr("^Init\\sE\\s\\(V\\)", chifile)
   InitE <- as.numeric(strsplit(chifile[which(position.InitE == 1)], "\\s=\\s")[[1]][2])
   ff$InitE <- InitE
   # HighE (volt)
   position.HighE <- regexpr("^High\\sE\\s\\(V\\)", chifile)
   HighE <- as.numeric(strsplit(chifile[which(position.HighE == 1)], "\\s=\\s")[[1]][2])
   ff$HighE <- HighE
   # LowE (volt)
   position.LowE <- regexpr("^Low\\sE\\s\\(V\\)", chifile)
   LowE <- as.numeric(strsplit(chifile[which(position.LowE == 1)], "\\s=\\s")[[1]][2])
   ff$LowE <- LowE
   # InitPN
   position.InitPN <- regexpr("^Init\\sP/N", chifile)
   InitPN <- strsplit(chifile[which(position.InitPN == 1)], "\\s=\\s")[[1]][2]
   ff$InitPN <- InitPN
   # Steps (total number of steps)
   position.Steps <- regexpr("^Step\\s", chifile)
   Steps <- as.numeric(strsplit(chifile[which(position.Steps == 1)], "\\s=\\s")[[1]][2])
   ff$Steps <- Steps
   # Pulse width (s)
   position.PulseWidth <- regexpr("^Pulse\\sWidth\\s\\(sec\\)", chifile)
   PulseWidth <- as.numeric(strsplit(chifile[which(position.PulseWidth == 1)], "\\s=\\s")[[1]][2])
   ff$PulseWidth <- PulseWidth
   # Sample interval (s)
   position.SampleInterval <- regexpr("^Sample\\sInterval\\s\\(s\\)", chifile)
   SampleInterval <- as.numeric(strsplit(chifile[which(position.SampleInterval == 1)], "\\s=\\s")[[1]][2])
   ff$SampleInterval <- SampleInterval
   # Quiet Time (s)
   position.QuietTime <- regexpr("^Quiet\\sTime\\s\\(sec\\)", chifile)
   QuietTime <- as.numeric(strsplit(chifile[which(position.QuietTime == 1)], "\\s=\\s")[[1]][2])
   ff$QuietTime <- QuietTime
   # Sensitivity (A/V)
   position.Sensitivity <- regexpr("^Sensitivity\\s\\(A/V\\)", chifile)
   Sensitivity <- as.numeric(strsplit(chifile[which(position.Sensitivity == 1)], "\\s=\\s")[[1]][2])
   ff$Sensitivity <- Sensitivity
   #
   return(ff)
}



##################################################
############### amperometry2df ###################
##################################################
amperometry2df <- function(datafilename, wearea = 1) {
   ## Description:
   ##   Reads current-time data (from CHI 760 potentiostat)
   ##   and returns a dataframe with the data, 
   ##   the data attributes (experimental conditions),
   ##   and some calculated parameters (charge, didt, etc.) 
   ## Usage:
   ##   amperometry2df(datafilename, wearea)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ##         wearea: (optional) area of working electrode (in square centimeter)
   ## Value:
   ##   Dataframe with the following columns (and no extra attributes):
   ##   $ sampleid        : chr
   ##   $ time            : num
   ##   $ current         : num
   ##   $ currentdensity  : num
   ##   $ timediff        : num
   ##   $ dIdt            : num
   ##   $ didt            : num
   ##   $ charge          : num
   ##   $ chargedensity   : num
   ##   $ InitE           : num
   ##   $ SampleInterval  : num
   ##   $ RunTime         : num
   ##   $ QuietTime       : num
   ##   $ Sensitivity     : num
   ## Note:
   ##   The CH Instruments 760 potentiostat records all data 
   ##   using standard SI units, therefore this function
   ##   assumes all potential values to be in volts, 
   ##   currents to be in amperes, charges in Coulombs, 
   ##   time in seconds, and so on.
   #
   datafile <- file(datafilename, "r")
   chifile <- readLines(datafile, n = -1) #read all lines of input file
   close(datafile)
   #
   sampleid <- ProvideSampleId(datafilename)
   #
   rgxp.number <- "^\\-?\\d\\.\\d+[e,]"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers with or without a negative sign (hyphen), 
   # followed by one digit before the decimal, a decimal point,
   # and an arbitrary number of digits after the decimal point,
   # immediately followed by either the letter 'e' or a comma.
   # Note that backslashes are escaped due to the way R handles strings.
   #
   numrow.idx <- regexpr(rgxp.number, chifile)
   # Save the match length attribute to another variable,
   numrow.len <- attr(numrow.idx, "match.length")
   # then scrap the attribute of the original variable.
   attr(numrow.idx, "match.length") <- NULL
   #
   i <- seq(1, length(numrow.idx) - 1, 1)
   j <- seq(2, length(numrow.idx), 1)
   # Start indices of data ranges
   starts <- which(numrow.idx[i] != 1 & numrow.idx[j] == 1) + 1
   # End indices, except for the last
   ends <- which(numrow.idx[i] == 1 & numrow.idx[j] != 1)
   # Fix the last index of end indices
   ends <- c(ends, length(numrow.idx))
   #
   ff <- data.frame(NULL)
   for (s in 1:length(starts)) {
      zz <- textConnection(chifile[starts[s]:ends[s]], "r")
      ff <- rbind(ff,
            data.frame(stringsAsFactors = FALSE,
               sampleid, matrix(scan(zz, what = numeric(), sep = ","),
               ncol = 2, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "time", "current")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff, currentdensity = currentdensity)
   # Calculate time and current diffs
   timediff <- c(ff$time[1], diff(ff$time))
   currentdiff <- c(ff$current[1], diff(ff$current))
   currentdensitydiff <- c(ff$currentdensity[1], diff(ff$currentdensity))
   # Calculate differential of current and current density
   dIdt <- currentdiff / timediff
   didt <- currentdensitydiff / timediff
   # Calculate charge and charge density
   charge <- cumsum(ff$current)
   chargedensity <- cumsum(ff$currentdensity)
   # Update ff dataframe
   ff <- cbind(ff, 
            timediff = timediff, 
            dIdt = dIdt, 
            didt = didt,
            charge = charge, 
            chargedensity = chargedensity)
   #
   ### Collect attributes of this experiment
   # InitE (volt)
   position.InitE <- regexpr("^Init\\sE\\s\\(V\\)", chifile)
   InitE <- as.numeric(strsplit(chifile[which(position.InitE == 1)], "\\s=\\s")[[1]][2])
   ff$InitE <- InitE
   # SampleInterval (volt)
   position.SampleInterval <- regexpr("^Sample\\sInterval\\s\\(s\\)", chifile)
   SampleInterval <- as.numeric(strsplit(chifile[which(position.SampleInterval == 1)], "\\s=\\s")[[1]][2])
   ff$SampleInterval <- SampleInterval
   # Run time (seconds)
   position.RunTime <- regexpr("^Run\\sTime\\s\\(sec\\)", chifile)
   RunTime <- as.numeric(strsplit(chifile[which(position.RunTime == 1)], "\\s=\\s")[[1]][2])
   ff$RunTime <- RunTime
   # Quiet time (seconds)
   position.QuietTime <- regexpr("^Quiet\\sTime\\s\\(sec\\)", chifile)
   QuietTime <- as.numeric(strsplit(chifile[which(position.QuietTime == 1)], "\\s=\\s")[[1]][2])
   ff$QuietTime <- QuietTime
   # Sensitivity (ampere per volt)
   position.Sensitivity <- regexpr("^Sensitivity\\s\\(A/V\\)", chifile)
   Sensitivity <- as.numeric(strsplit(chifile[which(position.Sensitivity == 1)], "\\s=\\s")[[1]][2])
   ff$Sensitivity <- Sensitivity
   #
   return(ff)
}



##################################################
#################### cv2df #######################
##################################################
cv2df <- function(datafilename, wearea = 1) {
   # Function description: 
   # CH Instruments potentiostat records all data using standard SI units,
   # so all potential values are in volts, currents are in amperes,
   # charges in Coulombs, time in seconds, etc.
   #
   cvfile <- file(datafilename, "r")
   chifile <- readLines(cvfile, n = -1) #read all lines of input file
   close(cvfile)
   #
   sampleid <- ProvideSampleId(datafilename)
   #
   rgxp.number <- "^\\-?\\d\\.\\d+,"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers with or without a negative sign (hyphen), 
   # followed by one digit before the decimal, a decimal point,
   # and an arbitrary number of digits after the decimal point,
   # immediately followed by a comma.
   # Note that backslashes are escaped.
   #
   numrow.idx <- regexpr(rgxp.number, chifile)
   # Save the match length attribute to another variable,
   numrow.len <- attr(numrow.idx, "match.length")
   # then scrap the attribute of the original variable.
   attr(numrow.idx, "match.length") <- NULL
   #
   i <- seq(1, length(numrow.idx) - 1, 1)
   j <- seq(2, length(numrow.idx), 1)
   # Start indices of data ranges
   starts <- which(numrow.idx[i] != 1 & numrow.idx[j] == 1) + 1
   # End indices, except for the last
   ends <- which(numrow.idx[i] == 1 & numrow.idx[j] != 1)
   # Fix the last index of end indices
   ends <- c(ends, length(numrow.idx))
   #
   ff <- data.frame(NULL)
   for (s in 1:length(starts)) {
      zz <- textConnection(chifile[starts[s]:ends[s]], "r")
      ff <- rbind(ff,
               data.frame(sampleid, cycle = as.integer(ceiling(s/2)), segment = s,
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 3, byrow = T)))
      close(zz)
   }
   # Column names after initial assignment
   names(ff) <- c("sampleid", "cycle", "segment", "potential", "current", "charge")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff[, 1:5], currentdensity = currentdensity, charge = ff[, 6])
   ## Collect attributes of this experiment
   # InitE (volt)
   position.InitE <- regexpr("^Init\\sE\\s\\(V\\)", chifile)
   InitE <- as.numeric(strsplit(chifile[which(position.InitE == 1)], "\\s=\\s")[[1]][2])
   ff$InitE <- InitE
   # HighE (volt)
   position.HighE <- regexpr("^High\\sE\\s\\(V\\)", chifile)
   HighE <- as.numeric(strsplit(chifile[which(position.HighE == 1)], "\\s=\\s")[[1]][2])
   ff$HighE <- HighE
   # LowE (volt)
   position.LowE <- regexpr("^Low\\sE\\s\\(V\\)", chifile)
   LowE <- as.numeric(strsplit(chifile[which(position.LowE == 1)], "\\s=\\s")[[1]][2])
   ff$LowE <- LowE
   # InitPN (positive or negative)
   position.InitPN <- regexpr("^Init\\sP/N", chifile)
   InitPN <- strsplit(chifile[which(position.InitPN == 1)], "\\s=\\s")[[1]][2]
   ff$InitPN <- InitPN
   # ScanRate (volt per second)
   position.ScanRate <- regexpr("^Scan\\sRate\\s\\(V/s\\)", chifile)
   ScanRate <- as.numeric(strsplit(chifile[which(position.ScanRate == 1)], "\\s=\\s")[[1]][2])
   ff$ScanRate <- ScanRate
   # Segments, number of
   position.Segments <- regexpr("^Segment\\s=", chifile)
   Segments <- as.numeric(strsplit(chifile[which(position.Segments == 1)], "\\s=\\s")[[1]][2])
   ff$Segments <- Segments
   # SampleInterval (volt)
   position.SampleInterval <- regexpr("^Sample\\sInterval\\s\\(V\\)", chifile)
   SampleInterval <- as.numeric(strsplit(chifile[which(position.SampleInterval == 1)], "\\s=\\s")[[1]][2])
   ff$SampleInterval <- SampleInterval
   # Quiet time (seconds)
   position.QuietTime <- regexpr("^Quiet\\sTime\\s\\(sec\\)", chifile)
   QuietTime <- as.numeric(strsplit(chifile[which(position.QuietTime == 1)], "\\s=\\s")[[1]][2])
   ff$QuietTime <- QuietTime
   # Sensitivity (ampere per volt)
   position.Sensitivity <- regexpr("^Sensitivity\\s\\(A/V\\)", chifile)
   Sensitivity <- as.numeric(strsplit(chifile[which(position.Sensitivity == 1)], "\\s=\\s")[[1]][2])
   ff$Sensitivity <- Sensitivity
   #
   return(ff)
}




##################################################
################### lsv2df #######################
##################################################
lsv2df <- function(datafilename, wearea = 1) {
   ## Description:
   ##   Reads LSV datafiles from CHI 760 potentiostat
   ##   (potential, current, and charge)
   ##   and returns a dataframe with the data, 
   ##   the data attributes (experimental conditions),
   ##   and calculated current density and charge density.
   ## Usage:
   ##   lsv2df(datafilename, wearea)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ##         wearea: (optional) area of working electrode (in square centimeter)
   ## Value:
   ##   Dataframe with the following columns (and no extra attributes):
   ##   $ sampleid        : chr (id)
   ##   $ segment         : num (id)
   ##   $ potential       : num (measure)
   ##   $ current         : num (measure)
   ##   $ charge          : num (measure)
   ##   $ currentdensity  : num (measure)
   ##   $ chargedensity   : num (measure)
   ##   $ InitE           : num (id)
   ##   $ FinalE          : num (id)
   ##   $ ScanRate        : num (id)
   ##   $ SampleInterval  : num (id)
   ##   $ QuietTime       : num (id)
   ##   $ Sensitivity     : num (id)
   ## Note:
   ##   The CH Instruments 760 potentiostat records all data 
   ##   using standard SI units, therefore this function
   ##   assumes all potential values to be in volts, 
   ##   currents to be in amperes, charges in Coulombs, 
   ##   time in seconds, and so on.
   #
   lsvfile <- file(datafilename, "r")
   chifile <- readLines(lsvfile, n = -1) #read all lines of input file
   close(lsvfile)
   #
   sampleid <- ProvideSampleId(datafilename)
   #
   rgxp.number <- "^\\-?\\d\\.\\d+,"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers with or without a negative sign (hyphen), 
   # followed by one digit before the decimal, a decimal point,
   # and an arbitrary number of digits after the decimal point,
   # immediately followed by a comma.
   # Note that backslashes are escaped.
   #
   numrow.idx <- regexpr(rgxp.number, chifile)
   # Save the match length attribute to another variable,
   numrow.len <- attr(numrow.idx, "match.length")
   # then scrap the attribute of the original variable.
   attr(numrow.idx, "match.length") <- NULL
   #
   i <- seq(1, length(numrow.idx) - 1, 1)
   j <- seq(2, length(numrow.idx), 1)
   # Start indices of data ranges
   starts <- which(numrow.idx[i] != 1 & numrow.idx[j] == 1) + 1
   # End indices, except for the last
   ends <- which(numrow.idx[i] == 1 & numrow.idx[j] != 1)
   # Fix the last index of end indices
   ends <- c(ends, length(numrow.idx))
   #
   ff <- data.frame(NULL)
   for (s in 1:length(starts)) {
      zz <- textConnection(chifile[starts[s]:ends[s]], "r")
      ff <- rbind(ff,
               data.frame(stringsAsFactors = FALSE,
               sampleid, segment = s,
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 3, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "segment", "potential", "current", "charge")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff, currentdensity = currentdensity)
   # Calculate charge density
   chargedensity <- ff$charge / wearea
   ff <- cbind(ff, chargedensity = chargedensity)
   #
   ### Collect attributes of this experiment
   # InitE (volt)
   position.InitE <- regexpr("^Init\\sE\\s\\(V\\)", chifile)
   InitE <- as.numeric(strsplit(chifile[which(position.InitE == 1)], "\\s=\\s")[[1]][2])
   ff$InitE <- InitE
   # FinalE (volt)
   position.FinalE <- regexpr("^Final\\sE\\s\\(V\\)", chifile)
   FinalE <- as.numeric(strsplit(chifile[which(position.FinalE == 1)], "\\s=\\s")[[1]][2])
   ff$FinalE <- FinalE
   # ScanRate (volt per second)
   position.ScanRate <- regexpr("^Scan\\sRate\\s\\(V/s\\)", chifile)
   ScanRate <- as.numeric(strsplit(chifile[which(position.ScanRate == 1)], "\\s=\\s")[[1]][2])
   ff$ScanRate <- ScanRate
   # SampleInterval (volt)
   position.SampleInterval <- regexpr("^Sample\\sInterval\\s\\(V\\)", chifile)
   SampleInterval <- as.numeric(strsplit(chifile[which(position.SampleInterval == 1)], "\\s=\\s")[[1]][2])
   ff$SampleInterval <- SampleInterval
   # Quiet time (seconds)
   position.QuietTime <- regexpr("^Quiet\\sTime\\s\\(sec\\)", chifile)
   QuietTime <- as.numeric(strsplit(chifile[which(position.QuietTime == 1)], "\\s=\\s")[[1]][2])
   ff$QuietTime <- QuietTime
   # Sensitivity (ampere per volt)
   position.Sensitivity <- regexpr("^Sensitivity\\s\\(A/V\\)", chifile)
   Sensitivity <- as.numeric(strsplit(chifile[which(position.Sensitivity == 1)], "\\s=\\s")[[1]][2])
   ff$Sensitivity <- Sensitivity
   #
   return(ff)
}
