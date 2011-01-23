# CHI.R
# Functions to read and manipulate data from the CHI760 potentiostat/galvanostat
# Taha Ahmed, Jan 2011

# CONTENTS
# >>>> cv2df 
# >>>> lsv2df




##################################################
#################### cv2df #######################
##################################################
cv2df <- function(cvfilename) {
   # Function description: 
   # CH Instruments potentiostat records all data using standard SI units,
   # so all potential values are in volts, currents are in amperes,
   # charges in Coulombs, time in seconds, etc.
   #
   cvfile <- file(cvfilename, "r")
   chifile <- readLines(cvfile, n = -1) #read all lines of input file
   close(cvfile)
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
               data.frame(segment = factor(s),
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 3, byrow = T)))
      close(zz)
   }
   names(ff) <- c("segment", "potential", "current", "charge")
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
lsv2df <- function(lsvfilename) {
   # Function description: 
   # CH Instruments potentiostat records all data using standard SI units,
   # so all potential values are in volts, currents are in amperes,
   # charges in Coulombs, time in seconds, etc.
   #
   lsvfile <- file(lsvfilename, "r")
   chifile <- readLines(lsvfile, n = -1) #read all lines of input file
   close(lsvfile)
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
               data.frame(segment = factor(s),
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 3, byrow = T)))
      close(zz)
   }
   names(ff) <- c("segment", "potential", "current", "charge")
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
