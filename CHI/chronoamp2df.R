source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")

##################################################
################# chronoamp2df ###################
##################################################
chronoamp2df <- function(datafilename, wearea = 1) {
   ## Description:
   ##   Reads current-time data (from CHI 760 potentiostat)
   ##   and returns a dataframe with the data and the 
   ##   data attributes (experimental conditions).
   ## Usage:
   ##   chronoamp2df(datafilename, wearea)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ##         wearea: (optional) area of working electrode (in square centimeters)
   ## Value:
   ##   Dataframe with the following columns (and no extra attributes):
   ##   $ sampleid        : chr
   ##   $ step            : num
   ##   $ time            : num
   ##   $ current         : num
   ##   $ currentdensity  : num
   ##   $ InitE           : num
   ##   $ HighE           : num
   ##   $ LowE            : num
   ##   $ InitPN          : chr
   ##   $ Step            : num
   ##   $ Pulse width     : num
   ##   $ Sample interval : num
   ##   $ Quiet Time      : num
   ##   $ Sensitivity     : num
   ## Note:
   ##   The CH Instruments 760 potentiostat records all data 
   ##   using standard SI units, therefore this function
   ##   assumes all potential values to be in volts, 
   ##   currents to be in amperes, charges in Coulombs, 
   ##   time in seconds, and so on.
   #
   datafile <- file(datafilename, "r")
   chifile <- readLines(datafile, n = -1) #read all lines of input file
   close(datafile)
   #
   sampleid <- ProvideSampleId(datafilename)
   #
   rgxp.number <- "^\\-?\\d\\.\\d+[e,]"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers with or without a negative sign (hyphen), 
   # followed by one digit before the decimal, a decimal point,
   # and an arbitrary number of digits after the decimal point,
   # immediately followed by either the letter 'e' or a comma.
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
               data.frame(sampleid, step.current = s,
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 2, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "step.current", "time", "current")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff, currentdensity = currentdensity)
   #
   
   ### Collect attributes of this experiment
   # InitE (volt)
   position.InitE <- regexpr("^Init\\sE\\s\\(V\\)", chifile)
   InitE <- as.numeric(strsplit(chifile[which(position.InitE == 1)], "\\s=\\s")[[1]][2])
   ff$InitE <- InitE
   # HighE (volt)
   position.HighE <- regexpr("^High\\sE\\s\\(V\\)", chifile)
   HighE <- as.numeric(strsplit(chifile[which(position.HighE == 1)], "\\s=\\s")[[1]][2])
   ff$HighE <- HighE
   # LowE (volt)
   position.LowE <- regexpr("^Low\\sE\\s\\(V\\)", chifile)
   LowE <- as.numeric(strsplit(chifile[which(position.LowE == 1)], "\\s=\\s")[[1]][2])
   ff$LowE <- LowE
   # InitPN
   position.InitPN <- regexpr("^Init\\sP/N", chifile)
   InitPN <- strsplit(chifile[which(position.InitPN == 1)], "\\s=\\s")[[1]][2]
   ff$InitPN <- InitPN
   # Steps (total number of steps)
   position.Steps <- regexpr("^Step\\s", chifile)
   Steps <- as.numeric(strsplit(chifile[which(position.Steps == 1)], "\\s=\\s")[[1]][2])
   ff$Steps <- Steps
   # Pulse width (s)
   position.PulseWidth <- regexpr("^Pulse\\sWidth\\s\\(sec\\)", chifile)
   PulseWidth <- as.numeric(strsplit(chifile[which(position.PulseWidth == 1)], "\\s=\\s")[[1]][2])
   ff$PulseWidth <- PulseWidth
   # Sample interval (s)
   position.SampleInterval <- regexpr("^Sample\\sInterval\\s\\(s\\)", chifile)
   SampleInterval <- as.numeric(strsplit(chifile[which(position.SampleInterval == 1)], "\\s=\\s")[[1]][2])
   ff$SampleInterval <- SampleInterval
   # Quiet Time (s)
   position.QuietTime <- regexpr("^Quiet\\sTime\\s\\(sec\\)", chifile)
   QuietTime <- as.numeric(strsplit(chifile[which(position.QuietTime == 1)], "\\s=\\s")[[1]][2])
   ff$QuietTime <- QuietTime
   # Sensitivity (A/V)
   position.Sensitivity <- regexpr("^Sensitivity\\s\\(A/V\\)", chifile)
   Sensitivity <- as.numeric(strsplit(chifile[which(position.Sensitivity == 1)], "\\s=\\s")[[1]][2])
   ff$Sensitivity <- Sensitivity
   #
   return(ff)
}