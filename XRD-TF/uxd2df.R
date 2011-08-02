##################################################
#################### uxd2df ######################
##################################################
uxd2df <- function(uxdfile) {
   # Function for reading UXD files  # Assumptions: data in two columns
   # Args: uxdfile (filename with extension)
   # Returns: dataframe with three columns
   
   cchar <- "[;_]" #regexpr matching the comment characters used in Bruker's UXD
   cdata <- "[0-9]" #regexpr matching one character of any digit
   
   # A new file (datafile) containing only data will be created,
   # extension ".data" appended to uxdfile
   #datafile <- paste(uxdfile,".data",sep="")
   
   ufile <- file(uxdfile, "r")
   f <- readLines(ufile, n=-1) #read _all_ lines from UXD file
   close(ufile)
   
   # This way we identify data rows by looking for numeric characters.
   #wh <- regexpr("[0-9]", f)
   # This way we identify header rows
   # We assume that all other rows are data
   wh <- regexpr(cchar, f)
   
   mh <- wh[1:length(wh)] # this gives you the corresponding index vector
   # the value of each element corresponds to the position of the regexp match.
   # value = 1 means the first character of the row is cchar (row is header)
   # value =-1 means no cchar occur on the row (row is data)
   
   i <- seq(1, length(mh) - 1, 1)
   j <- seq(2, length(mh), 1)
   
   starts <- which(mh[i] == 1 & mh[j] != 1) + 1
   ends   <- length(mh)
   f <- f[starts:ends]
   
   rgxp.sampleid <- "[^/]*(?=\\.\\w*)" ## THIS REQUIRES perl=TRUE
   # Regular expression that extracts the filename out of a full path.
   # Matches and extracts everything from the last forward slash (assuming Unix slashes)
   # up until a dot folllowed by an arbitrary number of alphanumeric characters.
   sampleidmtch <- regexpr(rgxp.sampleid, uxdfile, perl=TRUE)
   # Check that there was a match
   if (sampleidmtch < 0) {
      # -1 means no match
      sampleid <- uxdfile
      # If match was unsuccessful we use the argument as passed to this function as sampleid
   }
   sampleid <- substr(uxdfile, sampleidmtch, (sampleidmtch + attr(sampleidmtch, "match.length") - 1))
   
   zz <- textConnection(f, "r")
   ff <- data.frame(sampleid, matrix(scan(zz,
         what = numeric()), ncol=2, byrow=T))
   names(ff) <- c("sampleid", "angle", "intensity")
   close(zz)
      
   #zz <- file(datafile, "w") #open connection to datafile
   #write.table(ff, file=datafile, row.names=F, sep=",")
   #close(zz)

   # Return dataframe
   ff
}
