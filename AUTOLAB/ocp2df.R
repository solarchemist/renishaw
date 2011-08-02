source("/home/taha/chepec/chetex/common/R/common.R")

##################################################
#################### ocp2df ######################
##################################################
ocp2df <- function(datafilename) {
   ## Description:
   ##   Reads voltage-time data (from AUTOLAB potentiostat)
   ##   and returns a dataframe with the data.
   ## Usage:
   ##   ocp2df(datafilename)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ## Value:
   ##   Dataframe with the following columns:
   ##   $ sampleid        : chr
   ##   $ time            : num
   ##   $ potential       : num
   #
   datafile <- file(datafilename, "r")
   chifile <- readLines(datafile, n = -1) #read all lines of input file
   close(datafile)
   #
   sampleid <- ProvideSampleId(datafilename)
   #
   rgxp.number <- "^\\d+\\.\\d+"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers _without_ a negative sign (hyphen), 
   # followed by one or more digits before the decimal, a decimal point,
   # and an one or more digits after the decimal point.
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
               data.frame(sampleid,
                  stringsAsFactors = FALSE,                  
                  matrix(scan(zz, what = numeric(), sep = ""),
                     ncol = 2, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "time", "potential")
   #
   return(ff)
}