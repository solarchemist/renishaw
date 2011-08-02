source("/home/taha/chepec/chetex/common/R/common.R")

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
   ff <- rbind(ff, data.frame(stringsAsFactors = FALSE,
            sampleid,
            matrix(scan(zz, what = numeric(), sep = "\t"),
            ncol = 2, byrow = T)))
   close(zz)
   names(ff) <- c("sampleid", "shift", "counts")
   # Re-order by increasing shift
   ff <- ff[order(ff$shift), ]
   # And fix the row.names
   row.names(ff) <- seq(1, dim(ff)[1])
   # Do not re-calculate the spectrum with evenly spaced points here!
   # You must first remove cosmic peaks, and as long as that is done
   # manually, re-calculation to evenly spaced shifts must also be
   # done manually.
   ##
   return(ff)
}