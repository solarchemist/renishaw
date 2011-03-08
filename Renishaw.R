# Renishaw.R
# Functions to read data from the Renishaw Raman spectrometer
# Taha Ahmed, Feb 2011

# CONTENTS
source("/home/taha/chepec/chetex/common/R/common.R")
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
   sampleid <- ProvideSampleId(datafilename)
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
