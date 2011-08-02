##################################################
################### muxd2mtx #####################
##################################################
muxd2mtx <- function(uxdfile, range) {
   # Function that reads an UXD file which contains several ranges
   # (created in a programmed run, for example)
   # Arguments
   # :: uxdfile (filename with extension)
   # :: range (an integer, the number of ranges in the uxdfile)
   # Returns: matrix with 2 columns
   
   cchar <- "[;_]" #regexpr matching the comment characters used in Bruker's UXD
   cdata <- "[0-9]" #regexpr matching one character of any digit
   # Create filenames for the output # no longer used, return dataframe instead
   #datafile <- paste(uxdfile,"-",range.descriptor,".data",sep="")
   
   # Read the input multirange file
   ufile <- file(uxdfile, "r")
   f <- readLines(ufile, n=-1) #read _all_ lines from UXD file
   close(ufile)
   
   # We identify header rows. We will assume that all other rows are data
   wh <- regexpr(cchar, f)
   # length(wh) equals length(f)
   # wh has either 1 or -1 for each element
   # value = 1 means the first character of that row is cchar (row is header)
   # value =-1 means absence of cchar on that row (row is data)
   
   # Since wh contains some attributes (given by regexpr function), we strip out everything
   # but the index vector and assign it to the new vector mh
   mh <- wh[1:length(wh)] # this gives you the corresponding index vector
   
   #length(mh[mh == -1]) #total number of datarows in uxdfile
   #mh[mh > 1 | mh < 0] <- 0 #set all header-rows to zero (just to make things easier)
   
   # Set counters i and j used in assignment below
   i <- seq(1, length(mh) - 1, 1)
   j <- seq(2, length(mh), 1)
   
   starts <- which(mh[i] == 1 & mh[j] != 1) + 1 #start indices
   ends   <- which(mh[i] != 1 & mh[j] == 1) #end indices, except the last
   ends   <- c(ends, length(mh)) #fixes the last index of ends
   # note that length of starts (or ends) gives the number of ranges
   
   ff <- matrix(NA, length(f), 2)
   for (s in 1:range) {
      zz <- textConnection(f[starts[s]:ends[s]], "r")
      ff <- rbind(ff, matrix(scan(zz, what = numeric()), ncol=2, byrow=T))
      close(zz)
   }
   
   # Clean up matrix: remove extra rows
   ff[apply(ff,1,function(x)any(!is.na(x))), ]

   # Return matrix
   ff
}
