source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")
source("/home/taha/chepec/chetex/common/R/common/int2padstr.R")


##################################################
#!!!!!!!!!!!!!!!!!# mraw2df #!!!!!!!!!!!!!!!!!!!!#
##################################################
mraw2df <- function(txtfile) {
   ## Description:
   ##   Reads TXT files with one or multiple ranges (XRD commander "save as txt")
   ##   Extracts both data (thth, intensity) and parameters
   ## Usage:
   ##   mraw2df(txtfile)
   ## Arguments:
   ##   txtfile: text string with full path to TXT file, which may
   ##            containing single or multiple data ranges
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
   ranges.total <- 
      ifelse(length(range.header.start.rows) == length(range.header.end.rows), 
             length(range.header.start.rows),
             NA) #why would they not be equal?
   if (is.na(ranges.total)) {
      # Obviously something bad happened.
      # Do something about it. echo an error message perhaps.
      # But why would they not be equal?
      
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
   # But only if data contained more than one range, obviously. Let's make the code check for that
   if (ranges.total > 1) {
      range.data.end.rows <- c(range.header.start.rows[2:length(range.header.start.rows)] - 1, length(f))
   } else {
      # Data in fact only contains one range
      range.data.end.rows <- length(f)
   }
   
   ####
   
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
                                  int2padstr(range, "0", 3),
                                  matrix(scan(zz, what = numeric()), ncol = 2, byrow = T))
      close(zz)
      # Collect header parameters
      for (param in 1:length(header.param.rexp)) {
         data[[range]] <- cbind(data[[range]], 
            as.numeric(strsplit(headers.raw[[range]][which(regexpr(unname(header.param.rexp[param]), 
               headers.raw[[range]]) == 1)], "=")[[1]][2]))
      }
      names(data[[range]]) <- 
         c("sampleid", "range", "thth", ifelse(counts.flag, "counts", "cps"), names(header.param.rexp))
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
   if (ranges.total > 1) {
      for (range in 2:ranges.total) {
         data.df <- rbind(data.df, data[[range]])
      }
   }
   
   return(data.df)
}
