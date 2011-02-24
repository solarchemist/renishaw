# Renishaw.R
# Functions to read data from the Renishaw Raman spectrometer
# Taha Ahmed, Feb 2011

# CONTENTS
# >>>> Raman2df





##################################################
################### Raman2df #######################
##################################################
Raman2df <- function(datafilename) {
   # Function description: for reading Raman spectrum into dataframe
   #
   datafile <- file(datafilename, "r")
   chifile <- readLines(datafile, n = -1) #read all lines of input file
   close(datafile)
   #
   #####
   # A nice algorithm that extracts the filename from the datafilename argument
   # and uses that as a sampleid in the returned dataframe
   #####
   rgxp.sampleid <- "[^/]*(?=\\.\\w*)" ## THIS REQUIRES perl=TRUE
   # Regular expression that extracts the filename out of a full path.
   # Matches and extracts everything from the last forward slash (assuming Unix slashes)
   # up until a dot folllowed by an arbitrary number of alphanumeric characters.
   sampleidmtch <- regexpr(rgxp.sampleid, datafilename, perl=TRUE)
   # Check that there was a match
   if (sampleidmtch < 0) {
      # -1 means no match
      sampleid <- datafilename
      # If match was unsuccessful we use the argument as passed to this function as sampleid
   }
   sampleid <- substr(datafilename, sampleidmtch, (sampleidmtch + attr(sampleidmtch, "match.length") - 1))
   #
   ff <- data.frame(NULL)
   zz <- textConnection(chifile, "r")
   ff <- rbind(ff, data.frame(sampleid,
            matrix(scan(zz, what = numeric(), sep = "\t"),
            ncol = 2, byrow = T)))
   close(zz)
   names(ff) <- c("sampleid", "shift", "counts")
   #
   return(ff)
}
