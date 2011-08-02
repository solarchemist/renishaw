source("/home/taha/chepec/chetex/common/R/common.R")

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