# CHI.R
# Functions to read and manipulate data from the CHI760 potentiostat/galvanostat
# Taha Ahmed, Jan 2011 - Feb 2011

# CONTENTS
source("/home/taha/chepec/chetex/common/R/common.R")
# >>>> ocp2df
# >>>> chronocm2df
# >>>> chronoamp2df
# >>>> amperometry2df
# >>>> cv2df 
# >>>> lsv2df




##################################################
################### ocp2df #######################
##################################################
ocp2df <- function(datafilename) {
   # Function description: for recorded amperometric i-T curves
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
               data.frame(sampleid, matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 2, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "time", "potential")
   #
   ### Collect attributes of this experiment
   # These attributes are specific for each kind of experiment,
   # be careful when adapting to other electrochemical data
   rgxp.attr <- c("^Run\\sTime\\s\\(sec\\)")
   names.attr <- c("RunTime")
   for (n in 1:length(rgxp.attr)) {
      attrow.idx <- regexpr(rgxp.attr[n], chifile)
      attrow.len <- attr(attrow.idx, "match.length")
      attr(attrow.idx, "match.length") <- NULL
      # attrow.idx should now contain only one matching row
      attr(ff, names.attr[n]) <- strsplit(chifile[which(attrow.idx == 1)],
         "\\s=\\s")[[1]][2]
   }
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
   # Function description: chronoamperometry data
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
   names(ff) <- c("sampleid", "step", "time", "current")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff, currentdensity = currentdensity)
   #
   ### Collect attributes of this experiment
   # These attributes are specific for each kind of experiment,
   # be careful when adapting to other electrochemical data
   rgxp.attr <- c("^Init\\sE\\s\\(V\\)",
                  "^High\\sE\\s\\(V\\)",
                  "^Low\\sE\\s\\(V\\)",
                  "^Init\\sP/N",
                  "^Step\\s",
                  "^Pulse\\sWidth\\s\\(sec\\)",
                  "^Sample\\sInterval\\s\\(s\\)",
                  "^Quiet\\sTime\\s\\(sec\\)",
                  "^Sensitivity\\s\\(A/V\\)")
   names.attr <- c("InitE",
                   "HighE",
                   "LowE",
                   "InitPN",
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
############### amperometry2df ###################
##################################################
amperometry2df <- function(datafilename, wearea = 1) {
   # Function description: for recorded amperometric i-T curves
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
               data.frame(sampleid, matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 2, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "time", "current")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff, currentdensity = currentdensity)
   #
   ### Collect attributes of this experiment
   # These attributes are specific for each kind of experiment,
   # be careful when adapting to other electrochemical data
   rgxp.attr <- c("^Init\\sE\\s\\(V\\)",
                  "^Sample\\sInterval\\s\\(s\\)",
                  "^Run\\sTime\\s\\(sec\\)",
                  "^Quiet\\sTime\\s\\(sec\\)",
                  "^Sensitivity\\s\\(A/V\\)")
   names.attr <- c("InitE",
                   "SamplingInterval",
                   "RunTime",
                   "QuietTime",
                   "Sensitivity")
   for (n in 1:length(rgxp.attr)) {
      attrow.idx <- regexpr(rgxp.attr[n], chifile)
      attrow.len <- attr(attrow.idx, "match.length")
      attr(attrow.idx, "match.length") <- NULL
      # attrow.idx should now contain only one matching row
      attr(ff, names.attr[n]) <- strsplit(chifile[which(attrow.idx == 1)],
         "\\s=\\s")[[1]][2]         
   }
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
               data.frame(sampleid, segment = factor(s), cycle = factor(ceiling(s/2)),
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 3, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "segment", "cycle", "potential", "current", "charge")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff[, 1:5], currentdensity = currentdensity, charge = ff[, 6])
   #
   ### Collect attributes of this experiment
   # These attributes are specific for each kind of experiment,
   # be careful when adapting to other electrochemical data
   rgxp.attr <- c("^Init\\sE\\s\\(V\\)",
                  "^High\\sE\\s\\(V\\)",
                  "^Low\\sE\\s\\(V\\)",
                  "^Init\\sP/N",
                  "^Scan\\sRate\\s\\(V/s\\)",
                  "^Segment\\s=",
                  "^Sample\\sInterval\\s\\(V\\)",
                  "^Quiet\\sTime\\s\\(sec\\)",
                  "^Sensitivity\\s\\(A/V\\)")
   names.attr <- c("InitE",
                   "HighE",
                   "LowE",
                   "InitPN",
                   "ScanRate",
                   "Segments",
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
################### lsv2df #######################
##################################################
lsv2df <- function(datafilename, wearea = 1) {
   # Function description: 
   # CH Instruments potentiostat records all data using standard SI units,
   # so all potential values are in volts, currents are in amperes,
   # charges in Coulombs, time in seconds, etc.
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
               data.frame(sampleid, segment = factor(s),
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 3, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "segment", "potential", "current", "charge")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff[, 1:4], currentdensity = currentdensity, charge = ff[, 5])
   #
   ### Collect attributes of this experiment
   # These attributes are specific for each kind of experiment,
   # be careful when adapting to other electrochemical data
   rgxp.attr <- c("^Init\\sE\\s\\(V\\)",
                  "^Final\\sE\\s\\(V\\)",
                  "^Scan\\sRate\\s\\(V/s\\)",
                  "^Sample\\sInterval\\s\\(V\\)",
                  "^Quiet\\sTime\\s\\(sec\\)",
                  "^Sensitivity\\s\\(A/V\\)")
   names.attr <- c("InitE",
                   "FinalE",
                   "ScanRate",
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
