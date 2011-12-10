source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")

Raman2df <- function(datafilename) {
   ## Description:
   ##   Reads Raman data in ASCII format 
   ##   (wavenumber, counts)
   ##   and returns a dataframe with the original data, 
   ##   as well as interpolated wavenumber and counts values
   ##   (interpolated data is evenly spaced along x-axis)
   ## Usage:
   ##   Raman2df(datafilename)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ## Value:
   ##   Dataframe with the following columns (and no extra attributes):
   ##   $ sampleid        : chr (id)
   ##   $ wavenum         : num (measure)
   ##   $ counts          : num (measure)
   ##   $ wnum.interp     : num (measure)
   ##   $ cts.interp      : num (measure)
   ## Note:
   ##   
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
   names(ff) <- c("sampleid", "wavenum", "counts")
   # Sort dataframe by increasing wavenumbers
   ff <- ff[order(ff$wavenum), ]
   # ... sort the rownames as well
   row.names(ff) <- seq(1, dim(ff)[1])
   # Add interpolated, evenly spaced data to dataframe
   ff <- cbind(ff, 
               wnum.interp = approx(x = ff$wavenum, 
                                    y = ff$counts, 
                                    method = "linear", 
                                    n = length(ff$wavenum))$x,
               cts.interp = approx(x = ff$wavenum, 
                                    y = ff$counts, 
                                    method = "linear", 
                                    n = length(ff$wavenum))$y)
   ##
   return(ff)
}