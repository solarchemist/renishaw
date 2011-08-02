source("/home/taha/chepec/chetex/common/R/common.R")

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