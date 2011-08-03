source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")


##################################################
################### muxd2df ######################
##################################################
muxd2df <- function(uxdfile) {
   ## Description:
   ##   Reads UXD files with multiple ranges (converted using XCH v1.0)
   ##   Extracts both data (thth, intensity) and parameters
   ##   Also automatically calculates cps is counts are present, and vice versa
   ##   (note that this depends on specific strings in the UXD format).
   ## Usage:
   ##   muxd2df(uxdfile)
   ## Arguments:
   ##   uxdfile: text string with full path to UXD file
   ## Value:
   ##   Dataframe with the following columns:
   ##   $ sampleid        : chr
   ##   $ thth            : num
   ##   $ counts (or cps) : num
   ##   $ steptime        : num
   ##   $ stepsize        : num
   ##   $ theta           : num
   ##   $ khi             : num
   ##   $ phi             : num
   ##   $ x               : num
   ##   $ y               : num
   ##   $ z               : num
   ##   $ divergence      : num
   ##   $ antiscatter     : num
   ##   $ cps (or counts) : num
   #
   range.header.start.rexp <- "^; \\(Data for Range" #regexp
   range.header.end.rexp <- "^_2THETA[^=]" #regexp
   
   # Read the input multirange file
   ufile <- file(uxdfile, "r")
   # Note that readLines apparently completely skips empty lines. 
   # In that case line numbers do not match between source and f.
   f <- readLines(ufile, n=-1) #read _all_ lines from UXD file
   close(ufile)
   
   # Fetch a sampleid for the current job
   sampleid <- ProvideSampleId(uxdfile)
   
   # Look for header start rows
   range.header.start.rows <- which(regexpr(range.header.start.rexp, f) == 1)
   # Look for header end rows
   range.header.end.rows <- which(regexpr(range.header.end.rexp, f) == 1)
   
   # Calculate number of ranges
   ranges.total <- ifelse(length(range.header.start.rows) == length(range.header.end.rows), length(range.header.start.rows), NA)
   if (is.na(ranges.total)) {
      # Obviously something bad happened.
      # Do something about it. echo an error message perhaps.
      
   }
         
   # Determine whether we have COUNTS or COUNTS PER SECOND in current UXD-file
   # Assuming it is the same for all ranges in this job (a safe assumption).
   if (f[range.header.end.rows][1] == "_2THETACOUNTS") {
      # we got counts
      counts.flag <- TRUE
      cps.flag <- FALSE
   }
   if (f[range.header.end.rows][1] == "_2THETACPS") {
      # we got counts per second
      counts.flag <-FALSE
      cps.flag <- TRUE
   }
   
   # Extract headers (as-is) and put them in a list (by range)
   headers.raw <- list()
   for (range in 1:ranges.total) {
      headers.raw[[range]] <- f[range.header.start.rows[range]:range.header.end.rows[range]]
   }

   # Data always start on the row after header end
   range.data.start.rows <- range.header.end.rows + 1
   # Data end rows precedes header with one row, except for the first range
   range.data.end.rows <- c(range.header.start.rows[2:length(range.header.start.rows)] - 1, length(f))
   
   # Extract data (as-is) and put it an list (by range)
   data.raw <- list()
   for (range in 1:ranges.total) {
      data.raw[[range]] <- f[range.data.start.rows[range]:range.data.end.rows[range]]
   }
   
   # Specify header parameters to include in dataframe
   header.param.rexp <- c(steptime = "^_STEPTIME=", 
                          stepsize = "^_STEPSIZE=", 
                          theta = "^_THETA=",
                          khi = "^_KHI=",
                          phi = "^_PHI=",
                          x = "^_X=",
                          y = "^_Y=",
                          z = "^_Z=",
                          divergence = "^_DIVERGENCE=",
                          antiscatter = "^_ANTISCATTER=")
   
   # Collect data and header parameters in dataframes, by range in a list
   data <- list()
   for (range in 1:ranges.total) {
      zz <- textConnection(data.raw[[range]], "r")
      data[[range]] <- data.frame(stringsAsFactors = F,
                                  sampleid,
                                  matrix(scan(zz, what = numeric()), ncol = 2, byrow = T))
      close(zz)
      # Collect header parameters
      for (param in 1:length(header.param.rexp)) {
         data[[range]] <- cbind(data[[range]], 
            as.numeric(strsplit(headers.raw[[range]][which(regexpr(unname(header.param.rexp[param]), 
               headers.raw[[range]]) == 1)], "=")[[1]][2]))
      }
      names(data[[range]]) <- c("sampleid", "thth", ifelse(counts.flag, "counts", "cps"), names(header.param.rexp))
   }
                                                                   
   # Calculate the other of the pair counts <-> cps
   if (counts.flag) {
      for (range in 1:ranges.total) {
         data[[range]] <- cbind(data[[range]], cps = data[[range]]$counts / data[[range]]$steptime)
      }
   }
   if (cps.flag) {
      for (range in 1:ranges.total) {
         data[[range]] <- cbind(data[[range]], counts = data[[range]]$cps * data[[range]]$steptime)
      }
   }
      
   # Return a unified dataframe
   data.df <- data[[1]]
   for (range in 2:ranges.total) {
      data.df <- rbind(data.df, data[[range]])
   }
   
   return(data.df)
}





#### OLD VERSION - DEPRECATE
##################################################
################### muxd2df ######################
##################################################
muxd2df.old <- function(uxdfile, range.descriptor) {
   # Function that reads an UXD file which contains several ranges
   # (created in a programmed run, for example)
   # Arguments
   # :: uxdfile (filename with extension)
   # :: range.descriptor (an array with as many elements as
   #    there are ranges in the uxdfile)
   # Returns: dataframe with 3 columns
   
   cchar <- "[;_]" #regexpr matching the comment characters used in Bruker's UXD
   cdata <- "[0-9]" #regexpr matching one character of any digit
   # Create filenames for the output # no longer used, return dataframe instead
   #datafile <- paste(uxdfile,"-",range.descriptor,".data",sep="")
   
   # Read the input multirange file
   ufile <- file(uxdfile, "r")
   f <- readLines(ufile, n=-1) #read _all_ lines from UXD file
   close(ufile)
   
   # This way we identify data rows by looking for numeric characters.
   #wh <- regexpr("[0-9]", f)
   # This way we identify header rows
   # Later we will assume that all other rows are data
   wh <- regexpr(cchar, f)
   
   mh <- wh[1:length(wh)] # this gives you the corresponding index vector
   # the value of each element corresponds to the position of the regexp match.
   # value = 1 means the first character of the row is cchar (row is header)
   # value =-1 means no cchar occur on the row (row is data)
   
   #length(mh[mh == -1]) #total number of datarows in uxdfile
   #mh[mh > 1 | mh < 0] <- 0 #set all header-rows to zero (just to make things easier)
   
   i <- seq(1, length(mh) - 1, 1)
   j <- seq(2, length(mh), 1)
   starts <- which(mh[i] == 1 & mh[j] != 1) + 1 #start indices
   ends   <- which(mh[i] != 1 & mh[j] == 1) #end indices, except the last
   ends   <- c(ends, length(mh)) #fixed the last index of ends   
   
   ff <- data.frame(NULL)
   for (s in 1:length(range.descriptor)) {
      zz <- textConnection(f[starts[s]:ends[s]], "r")
      ff <- rbind(ff, data.frame(range.descriptor[s],
            matrix(scan(zz, what = numeric()), ncol=2, byrow=T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "angle", "intensity")
   
   # Return dataframe
   ff
}
