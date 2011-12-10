##################################################
#################### SP2df #######################
##################################################
SP2df <- function(datafile) {
   ## Description:
   ##   For now just extracting the bare minimum (the data itself) from SP (ASCII) spectra files.
   ## Usage:
   ##   SP2df(datafile)
   ## Arguments:
   ##   datafile: text string with full path to TXT file
   ##             containing single or multiple data ranges
   ## Value:
   ##   Dataframe with the following columns:
   ##   $ sampleid        : chr
   ##   $ wavelength      : num
   ##   $ intensity       : num
   #
   range.data.start.rexp <- "\\#DATA"
   #range.data.end.rexp <- ">+End[\\s\\w]*<+"
   
   # Read the input file
   dfile <- file(datafile, "r")
   # Note that readLines apparently completely skips empty lines. 
   # That causes line numbers to not match between source and f vector.
   f <- readLines(dfile, n=-1) # read _all_ lines from data file
   close(dfile)
   
   # Create a sampleid for the current job (use the folder name)
   sampleid <- basename(dirname(datafile))
   
   # Look for data start marker line
   range.data.start.row <- grep(range.data.start.rexp, f, perl = TRUE) + 1
   # Data ends one line before EOF
   range.data.end.row <- length(f) - 1
   
   # Extract data (as-is)
   data.raw <- f[range.data.start.row:range.data.end.row]

   
   # Collect data into dataframe
   zz <- textConnection(data.raw, "r")
   data <- data.frame(stringsAsFactors = FALSE,
                      sampleid,
                      matrix(scan(zz, what = numeric()), ncol = 2, byrow = T))
   close(zz)
   names(data) <- c("sampleid", "wavelength", "intensity")
                                                                   
   return(data)
}
