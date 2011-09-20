source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")

##################################################
################### mps2df #######################
##################################################
mps2df <- function(datafilename, wearea = 1) {
   ## Description:
   ##   Reads time vs current data (from CHI 760 potentiostat)
   ##   and returns a dataframe with the data and 
   ##   the data attributes (experimental conditions).
   ## Usage:
   ##   mps2df(datafilename, wearea)
   ## Arguments:
   ##   datafilename: text string with full path to experimental file
   ##         wearea: working electrode area in square centimeters (optional)
   ## Value:
   ##   Dataframe with the following columns (and no extra attributes):
   ##   $ sampleid        : chr
   ##   $ time            : num
   ##   $ current         : num
   ##   $ currentdensity  : num
   ##   $ timediff        : num
   ##   $ dIdt            : num
   ##   $ didt            : num
   ##   $ charge          : num
   ##   $ chargedensity   : num
   ##   $ PotentialSteps  : num
   ##   $ TimeSteps       : num
   ##   $ Cycle           : num
   ##   $ SampleInterval  : num
   ##   $ QuietTime       : num
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
   rgxp.number <- "^\\-?\\d+\\.\\d+[e,]"
   # regexp that matches a decimal number at the beginning of the line.
   # Matches numbers with or without a negative sign (hyphen), 
   # followed by one digit before the decimal, a decimal point,
   # and an arbitrary number of digits after the decimal point,
   # immediately followed by either the letter 'e' or a comma.
   # Note that backslashes are escaped due to the way R handles strings.
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
               sampleid, matrix(scan(zz, what = numeric(), sep = ","),
               ncol = 2, byrow = T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "time", "current")
   # Calculate current density
   currentdensity <- ff$current / wearea
   ff <- cbind(ff, currentdensity = currentdensity)
   # Calculate time diff and current diff
   timediff <- c(ff$time[1], diff(ff$time))
   currentdiff <- c(ff$current[1], diff(ff$current))
   currentdensitydiff <- c(ff$currentdensity[1], diff(ff$currentdensity))
   # Calculate differential of current and current density
   dIdt <- currentdiff / timediff
   didt <- currentdensitydiff / timediff
   # Calculate charge and charge density
   charge <- cumsum(ff$current)
   chargedensity <- cumsum(ff$currentdensity)
   # Update ff dataframe
   ff <- cbind(ff, 
      timediff = timediff, 
      dIdt = dIdt, 
      didt = didt,
      charge = charge, 
      chargedensity = chargedensity)
   #
   ### Collect attributes of this experiment
   # Potential steps
   PotentialSteps <- ""
   positions.PotentialSteps <- regexpr("^E\\d+\\s\\(V\\)", chifile)
   PotentialSteps <- paste("\\num{", strsplit(chifile[which(positions.PotentialSteps == 1)], 
      "\\s=\\s")[[1]][2], "}", sep = "")
   if (length(which(positions.PotentialSteps == 1)) > 1) {
      for (i in 2:length(which(positions.PotentialSteps == 1))) {
      PotentialSteps <- 
         paste(PotentialSteps, 
            paste("\\num{", strsplit(chifile[which(positions.PotentialSteps == 1)], 
            "\\s=\\s")[[i]][2], "}", sep = ""), sep = ", ")
      }
   }
   ff$PotentialSteps <- PotentialSteps
   # Time steps
   TimeSteps <- ""
   positions.TimeSteps <- regexpr("^T\\d+\\s\\(s\\)", chifile)
   TimeSteps <- paste("\\num{", strsplit(chifile[which(positions.TimeSteps == 1)], 
      "\\s=\\s")[[1]][2], "}", sep = "")
   if (length(which(positions.TimeSteps == 1)) > 1) {
      for (i in 2:length(which(positions.TimeSteps == 1))) {
      TimeSteps <- 
         paste(TimeSteps, 
            paste("\\num{", strsplit(chifile[which(positions.TimeSteps == 1)], 
            "\\s=\\s")[[i]][2], "}", sep = ""), sep = ", ")
      }
   }
   ff$TimeSteps <- TimeSteps
   # Cycles
   position.Cycles <- regexpr("^Cycle", chifile)
   Cycles <- as.numeric(strsplit(chifile[which(position.Cycles == 1)], "\\s=\\s")[[1]][2])
   ff$Cycles <- Cycles
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