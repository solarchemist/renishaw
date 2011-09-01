source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")

xrfspectro2df <- function(smpfile) {
   ## Description:
   ##   Reads XRF textfile from XLAB SPECTRO XRF.
   ##   Stores data in data frame and parameters in an attributed dataframe.
   ## Usage:
   ##   xrfspectro2df(smpfile)
   ## Arguments:
   ##    smpfile: character string, the full filename
   ##                (with path) to one SMP file (ASCII).
   ## Value:
   ##   A dataframe with attributed dataframe
   #
   filecon <- file(smpfile, "r")
   smpcontents <- readLines(filecon, n = -1) #read all lines of input file
   close(filecon)
   #
   sampleid <- ProvideSampleId(smpfile)
   #
   rgxp.data <- "^Kanal\\s[\\d]+:"
   #
   numrow.idx <- regexpr(rgxp.data, smpcontents, perl = TRUE)
   # scrap the match.length attribute
   attr(numrow.idx, "match.length") <- NULL
   #
   # Determine how many columns the data contains
   smpdata.cols <- length(strsplit(smpcontents[which(numrow.idx == 1)][1], "\t")[[1]]) - 1
   # While we are at it, save row count to a variable as well
   smpdata.rows <- length(smpcontents[which(numrow.idx == 1)])
   # strip prefix off each data row
   #smpdata <- matrix(NA, ncol = smpdata.cols, nrow = smpdata.rows)
   smpdata.txt <- vector(length = smpdata.rows)
   for (i in 1:smpdata.rows) {
      smpdata.txt[i] <- strsplit(smpcontents[which(numrow.idx == 1)][i], ":")[[1]][2]
   }
   smpdata.txt.clean <- gsub("^\\s", "", 
                             gsub("\\t", " ", smpdata.txt))
   
   zz <- textConnection(smpdata.txt.clean, "r")
   ff <- data.frame(stringsAsFactors = FALSE, 
                    sampleid = sampleid,
                    channel = seq(1, smpdata.rows),
                    matrix(scan(zz, what = numeric(), sep = " "), 
                           ncol = smpdata.cols, byrow = TRUE))
   close(zz)
   names(ff) <- c("sampleid", "channel", paste("Y", seq(1, smpdata.cols), sep = ""))
   

   #
   ### Collect attributes of this experiment
   SMPattrEdit <- matrix(c("Voltage",      "^Voltage:",
                           "Current",      "^Current:",
                           "Target",       "^Target:",
                           "Duration",     "^Meas\\.\\sDuration:",
                           "Impulse",      "^Imp\\.\\sRate:",
                           "DeadTime",     "^Rel\\.\\sDead\\sTime:",
                           "FirstChannel", "^First\\sChannel:",
                           "LastChannel",  "^Last\\sChannel:",
                           "PeakTime",     "^Peak\\sTime:",
                           "Gain",         "^Gain:",
                           "ZeroPeak",     "^Zero\\sPeak:"),
                         ncol = 2, byrow = T)
                         
   SMPattr <- matrix(NA, nrow = smpdata.cols + 1, ncol = dim(SMPattrEdit)[1])
   
   for (c in 1:dim(SMPattrEdit)[1]) {
      SMPattr[1, c] <- SMPattrEdit[c, 1]
      SMPattr[2:dim(SMPattr)[1], c] <- 
         matrix(strsplit(gsub("^\\t", "", 
            strsplit(smpcontents[which(regexpr(SMPattrEdit[c, 2], smpcontents) == 1)],
               ":")[[1]][2]), "\\t")[[1]], ncol = smpdata.cols)
   }
   SMPdf <- data.frame(stringsAsFactors = FALSE, 
                       SMPattr[2:dim(SMPattr)[1], ])
   colnames(SMPdf) <- SMPattr[1, ]
                           
                           
   ### Now calculate the energy (keV) scale (convert from channels to energy)
   ff$X <- ff$channel * (as.numeric(SMPdf$Gain[1]) / as.numeric(SMPdf$LastChannel[1]))

   # Attach parameters to returned dataframe
   attr(ff, "parameters") <- SMPdf
   #
   return(ff)
}