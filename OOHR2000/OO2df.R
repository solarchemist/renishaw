source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")
source("/home/taha/chepec/chetex/common/R/common/int2padstr.R")


##################################################
#################### OO2df #######################
##################################################
OO2df <- function(datafile) {
   ## Description:
   ##   
   ##   
   ##   
   ##   
   ## Usage:
   ##   OO2df(datafile)
   ## Arguments:
   ##   datafile: text string with full path to TXT file
   ##             containing single or multiple data ranges
   ## Value:
   ##   Dataframe with the following columns:
   ##   $ sampleid        : chr
   ##   $ wavelength      : num
   ##   $ counts          : num
   ##   $ 
   ##   $ 
   ##   $ 
   ##   $ 
   ##   $ 
   ##   $ 
   ##   $ 
   ##   $ 
   ##   $ 
   ##   $ 
   ##   $ cps ?
   #
#    range.header.start.rexp <- "^; \\(Data for Range" #regexp
#    range.header.end.rexp <- "^_2THETA[^=]" #regexp
   
   range.data.start.rexp <- ">+Begin[\\s\\w]*<+"
   range.data.end.rexp <- ">+End[\\s\\w]*<+"
   
   # Read the input file
   dfile <- file(datafile, "r")
   # Note that readLines apparently completely skips empty lines. 
   # That causes line numbers to not match between source and f vector.
   f <- readLines(dfile, n=-1) # read _all_ lines from data file
   close(dfile)
   
   # Fetch a sampleid for the current job
   sampleid <- ProvideSampleId(datafile)
   
#    # Look for header start rows
#    range.header.start.rows <- which(regexpr(range.header.start.rexp, f) == 1)
#    # Look for header end rows
#    range.header.end.rows <- which(regexpr(range.header.end.rexp, f) == 1)
   
   # Look for data start marker line
   range.data.start.rows <- which(regexpr(range.data.start.rexp, f, perl = TRUE) == 1) + 1
   # Look for data end marker line
   range.data.end.rows <- which(regexpr(range.data.end.rexp, f, perl = TRUE) == 1) - 1
   
   # Calculate number of ranges
   ranges.total <- 
      ifelse(length(range.data.start.rows) == length(range.data.end.rows), 
             length(range.data.start.rows),
             NA) #why would they not be equal?
   if (is.na(ranges.total)) {
      # Obviously something bad happened.
      # Do something about it. echo an error message perhaps.
      # But why would they not be equal?
      
   }
   
   
   
   # Header always precedes start of data
   range.header.end.rows <- range.data.start.rows - 2
   if (ranges.total > 1) {
      range.header.start.rows <- c(1, range.data.end.rows[2:length(range.data.end.rows)])
   } else {
      # Data in fact only contains one range
      range.header.start.rows <- 1
   }
   
   # Extract headers (as-is) and put them in a list (by range)
   headers.raw <- list()
   for (range in 1:ranges.total) {
      headers.raw[[range]] <- f[range.header.start.rows[range]:range.header.end.rows[range]]
   }

   ####
   
   # Extract data (as-is) and put it an list (by range)
   data.raw <- list()
   for (range in 1:ranges.total) {
      data.raw[[range]] <- f[range.data.start.rows[range]:range.data.end.rows[range]]
      # Replace commas by dots
      data.raw[[range]] <- gsub(",", "\\.", data.raw[[range]])
   }   
   
   
   # Specify header parameters to include in dataframe
   header.param.rexp <- 
      c(DateTime           = "^Date:", 
        IntegrationTime    = "^Integration Time \\(usec\\):", 
        n_Averaged         = "^Spectra Averaged:",
        Boxcar             = "^Boxcar Smoothing:",
        CorrElectricDark   = "^Correct for Electrical Dark:",
        StrobeLampEnabled  = "^Strobe/Lamp Enabled:",
        CorrDetectorNonLin = "^Correct for Detector Non-linearity:",
        CorrStrayLight     = "^Correct for Stray Light:",
        n_Pixels           = "^Number of Pixels")
   
   # Collect data and header parameters in dataframes, by range in a list
   data <- list()
   for (range in 1:ranges.total) {
      zz <- textConnection(data.raw[[range]], "r")
      data[[range]] <- data.frame(stringsAsFactors = FALSE,
                                  sampleid,
                                  int2padstr(range, "0", 3),
                                  matrix(scan(zz, what = character()), ncol = 2, byrow = T))
      close(zz)
      # Collect header parameters
      for (param in 1:length(header.param.rexp)) {
         data[[range]] <- 
            cbind(stringsAsFactors = FALSE, 
                  data[[range]], 
                  # Matches any word, digit, plus, or minus character 
                  # surrounded by parentheses at the end of the string
                  sub("\\s\\([\\w\\d\\+\\-]+\\)$", "",
                      strsplit(headers.raw[[range]][which(regexpr(unname(header.param.rexp[param]), 
                               headers.raw[[range]]) == 1)], ": ")[[1]][2], 
                      perl = TRUE))
      }
      names(data[[range]]) <- 
         c("sampleid", "range", "wavelength", "intensity", names(header.param.rexp))
   }
                                                                   
   # Create a unified dataframe
   data.df <- data[[1]]
   if (ranges.total > 1) {
      for (range in 2:ranges.total) {
         data.df <- rbind(data.df, data[[range]])
      }
   }

   # Convert the DateTime column to more legibly (and compact) format
   data.df$DateTime <- 
      format(as.POSIXct(gsub("\\s[A-Z]{4}\\s", " ", data.df$Date), 
                        format = "%a %b %d %H:%M:%S %Y"), 
             format = "%Y-%m-%d %H:%M:%S")
   # Convert wavelength and intensity to numeric format
   mode(data.df$wavelength) <- "numeric"
   mode(data.df$intensity) <- "numeric"
   
   return(data.df)
}