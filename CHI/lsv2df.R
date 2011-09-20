source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")

##################################################
################### lsv2df #######################
##################################################
lsv2df <- function(datafilename, wearea = 1) {
   ## Description:
   ##   Reads LSV datafiles from CHI 760 potentiostat
   ##   (potential, current, and charge)
   ##   and returns a dataframe with the data, 
   ##   the data attributes (experimental conditions),
   ##   and calculated current density and charge density.
   ## Usage:
   ##   lsv2df(datafilename, wearea)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ##         wearea: (optional) area of working electrode (in square centimeter)
   ## Value:
   ##   Dataframe with the following columns (and no extra attributes):
   ##   $ sampleid        : chr (id)
   ##   $ segment         : num (id)
   ##   $ potential       : num (measure)
   ##   $ current         : num (measure)
   ##   $ charge          : num (measure)
   ##   $ currentdensity  : num (measure)
   ##   $ chargedensity   : num (measure)
   ##   $ InitE           : num (id)
   ##   $ FinalE          : num (id)
   ##   $ ScanRate        : num (id)
   ##   $ SampleInterval  : num (id)
   ##   $ QuietTime       : num (id)
   ##   $ Sensitivity     : num (id)
   ## Note:
   ##   The CH Instruments 760 potentiostat records all data 
   ##   using standard SI units, therefore this function
   ##   assumes all potential values to be in volts, 
   ##   currents to be in amperes, charges in Coulombs, 
   ##   time in seconds, and so on.
   #
   lsvfile <- file(datafilename, "r")
   chifile <- readLines(lsvfile, n = -1) #read all lines of input file
   close(lsvfile)
   #
   sampleid <- ProvideSampleId(datafilename)
   #
   rgxp.number <- "^\\-?\\d\\.\\d+,"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers with or without a negative sign (hyphen), 
   # followed by one digit before the decimal, a decimal point,
   # and an arbitrary number of digits after the decimal point,
   # immediately followed by a comma.
   # Note that backslashes are escaped.
   #
   numrow.idx <- regexpr(rgxp.number, chifile)
   # Save the match length attribute to another variable,
   numrow.len <- attr(numrow.idx, "match.length")
   # then scrap the attribute of the original variable.
   attr(numrow.idx, "match.length") <- NULL
   #
   i <- seq(1, length(numrow.idx) - 1, 1)
   j <- seq(2, length(numrow.idx), 1)
   # Start indices of data ranges
   starts <- which(numrow.idx[i] != 1 & numrow.idx[j] == 1) + 1
   # End indices, except for the last
   ends <- which(numrow.idx[i] == 1 & numrow.idx[j] != 1)
   # Fix the last index of end indices
   ends <- c(ends, length(numrow.idx))
   #
   ff <- data.frame(NULL)
   for (s in 1:length(starts)) {
      zz <- textConnection(chifile[starts[s]:ends[s]], "r")
      ff <- rbind(ff,
               data.frame(stringsAsFactors = FALSE,
               sampleid, segment = s,
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 3, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "segment", "potential", "current", "charge")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff, currentdensity = currentdensity)
   # Calculate charge density
   chargedensity <- ff$charge / wearea
   ff <- cbind(ff, chargedensity = chargedensity)
   #
   ### Collect attributes of this experiment
   # InitE (volt)
   position.InitE <- regexpr("^Init\\sE\\s\\(V\\)", chifile)
   InitE <- as.numeric(strsplit(chifile[which(position.InitE == 1)], "\\s=\\s")[[1]][2])
   ff$InitE <- InitE
   # FinalE (volt)
   position.FinalE <- regexpr("^Final\\sE\\s\\(V\\)", chifile)
   FinalE <- as.numeric(strsplit(chifile[which(position.FinalE == 1)], "\\s=\\s")[[1]][2])
   ff$FinalE <- FinalE
   # ScanRate (volt per second)
   position.ScanRate <- regexpr("^Scan\\sRate\\s\\(V/s\\)", chifile)
   ScanRate <- as.numeric(strsplit(chifile[which(position.ScanRate == 1)], "\\s=\\s")[[1]][2])
   ff$ScanRate <- ScanRate
   # SampleInterval (volt)
   position.SampleInterval <- regexpr("^Sample\\sInterval\\s\\(V\\)", chifile)
   SampleInterval <- as.numeric(strsplit(chifile[which(position.SampleInterval == 1)], "\\s=\\s")[[1]][2])
   ff$SampleInterval <- SampleInterval
   # Quiet time (seconds)
   position.QuietTime <- regexpr("^Quiet\\sTime\\s\\(sec\\)", chifile)
   QuietTime <- as.numeric(strsplit(chifile[which(position.QuietTime == 1)], "\\s=\\s")[[1]][2])
   ff$QuietTime <- QuietTime
   # Sensitivity (ampere per volt)
   position.Sensitivity <- regexpr("^Sensitivity\\s\\(A/V\\)", chifile)
   Sensitivity <- as.numeric(strsplit(chifile[which(position.Sensitivity == 1)], "\\s=\\s")[[1]][2])
   ff$Sensitivity <- Sensitivity
   #
   return(ff)
}