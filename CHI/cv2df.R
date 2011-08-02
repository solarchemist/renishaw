source("/home/taha/chepec/chetex/common/R/common.R")

##################################################
#################### cv2df #######################
##################################################
cv2df <- function(datafilename, wearea = 1) {
   # Function description: 
   # CH Instruments potentiostat records all data using standard SI units,
   # so all potential values are in volts, currents are in amperes,
   # charges in Coulombs, time in seconds, etc.
   #
   cvfile <- file(datafilename, "r")
   chifile <- readLines(cvfile, n = -1) #read all lines of input file
   close(cvfile)
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
               data.frame(sampleid, cycle = as.integer(ceiling(s/2)), segment = s,
               matrix(scan(zz, what = numeric(), sep = ","),
                  ncol = 3, byrow = T)))
      close(zz)
   }
   # Column names after initial assignment
   names(ff) <- c("sampleid", "cycle", "segment", "potential", "current", "charge")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff[, 1:5], currentdensity = currentdensity, charge = ff[, 6])
   ## Collect attributes of this experiment
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
   # InitPN (positive or negative)
   position.InitPN <- regexpr("^Init\\sP/N", chifile)
   InitPN <- strsplit(chifile[which(position.InitPN == 1)], "\\s=\\s")[[1]][2]
   ff$InitPN <- InitPN
   # ScanRate (volt per second)
   position.ScanRate <- regexpr("^Scan\\sRate\\s\\(V/s\\)", chifile)
   ScanRate <- as.numeric(strsplit(chifile[which(position.ScanRate == 1)], "\\s=\\s")[[1]][2])
   ff$ScanRate <- ScanRate
   # Segments, number of
   position.Segments <- regexpr("^Segment\\s=", chifile)
   Segments <- as.numeric(strsplit(chifile[which(position.Segments == 1)], "\\s=\\s")[[1]][2])
   ff$Segments <- Segments
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