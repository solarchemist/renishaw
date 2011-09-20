source("/home/taha/chepec/chetex/common/R/common/ProvideSampleId.R")

xrfspectro2df <- function(smpfile) {
   ## Description:
   ##   Total remake of xrfspectro2df(). Idea is to accomodate all 6 possible 
   ##   datasets of each SMP file, plus the attributes.
   ##   Reads XRF textfile from XLAB SPECTRO XRF.
   ##   Stores data in data frame and parameters in an attributed dataframe.
   ## Usage:
   ##   xrfspectro2df(smpfile)
   ## Arguments:
   ##    smpfile: character string, the full filename
   ##                (with path) to one SMP file (ASCII).
   ## Value:
   ##   A dataframe with attributed dataframe
   
   #### ONLY BOTHER WITH THE FIRST MEASUREMENT IN THE SMP-FILE.
   
   filecon <- file(smpfile, "r")
   smpcontents <- readLines(filecon, n = -1) #read all lines of input file
   close(filecon)
   
   # Parameter table
   # Those are the parameter we may access later in this function
   xrf.param <- data.frame(stringsAsFactors = FALSE,
      matrix(c("Method",       "^Method:",
               "Job",          "^Job:",
               "Status",       "^Status:",
               "Description",  "^Description:",
               "Date",         "^Date\\sof\\sMeasurement:",
               "Measurements", "^Measurements:",
               "Voltage",      "^Voltage:",
               "Current",      "^Current:",
               "Target",       "^Target:",
               "Duration",     "^Meas\\.\\sDuration:",
               "Impulse",      "^Imp\\.\\sRate:",
               "DeadTime",     "^Rel\\.\\sDead\\sTime:",
               "FirstChannel", "^First\\sChannel:",
               "LastChannel",  "^Last\\sChannel:",
               "PeakTime",     "^Peak\\sTime:",
               "Gain",         "^Gain:",
               "ZeroPeak",     "^Zero\\sPeak:",
               "Data",         "^Kanal\\s[\\d]+:"),
             ncol = 2, byrow = T))
   names(xrf.param) <- c("parameter", "regexp")
   
   
   # Data table
   # Contains the regexp used for identifiying rows containing data
   xrf.data <- data.frame(stringsAsFactors = FALSE,
      matrix(c("Data", "^Kanal\\s[\\d]+:"), ncol = 2, byrow = T))
   names(xrf.data) <- c("parameter", "regexp")
   
   
   
   # Find out how many measurements we have in this 
   # file by accessing the Measurements field
   n_measurements <- as.numeric(strsplit(gsub("^\\t", "", 
      strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "Measurements", select = "regexp")$regexp, 
         smpcontents) == 1)], ":")[[1]][2]), "\\t")[[1]])
   # If more than one measurement, issue warning
   if (n_measurements > 1) {
      warning(paste(paste(n_measurements, " measurements detected in ", 
                          basename(smpfile), sep = ""), 
                    "Only the first measurement will be recorded", 
                    sep = "\n", collapse = ""))
   }
   
   
   # How many rows of data?
   n_rowsdata <- 
      length(which(regexpr(subset(xrf.data, parameter == "Data", 
         select = "regexp")$regexp, smpcontents, perl = TRUE) == 1))
   
   # Build an empty matrix big enough to hold all data
   # (i.e., ncol = 3, and nrow = n_rowsdata * n_measurements)
   data.long <- data.frame(matrix(NA, ncol = 5, nrow = 6 * n_rowsdata))
   names(data.long) <- c("sampleid", "measurement", "channel", "energy", "counts")
   
   
   data.mtx <- matrix(NA, ncol = 6, nrow = n_rowsdata)
   for (j in 1:n_rowsdata) {
      data.mtx[j, ] <- as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.data, parameter == "Data", 
         select = "regexp")$regexp, smpcontents, perl = TRUE) == 1)], ":\\t")[[j]][2], "\\t")[[1]])
   }
                                                                          
   
   # Sampleid to column 1
   data.long[, 1] <- rep(ProvideSampleId(smpfile), n_rowsdata)
   # Channel to column 3
   data.long[, 3] <- rep(seq(1, n_rowsdata), dim(data.mtx)[2])
   for (c in 1:6) {
      # Measurement no. in column 2
      data.long[((c * n_rowsdata) - n_rowsdata + 1):(((c + 1) * n_rowsdata) - n_rowsdata), 2] <- rep(c, n_rowsdata)
      # Counts in column 5
      data.long[((c * n_rowsdata) - n_rowsdata + 1):(((c + 1) * n_rowsdata) - n_rowsdata), 5] <- data.mtx[, c]
   }
      
   # Drop all rows with measurement-number not equal to 1
   data.long <- subset(data.long, measurement == 1)


   # Fetch the measurement parameters
   data.long[, subset(xrf.param, parameter == "Date")$parameter] <-
   rep(substr(strsplit(gsub("^\\t", "", strsplit(smpcontents[which(regexpr(subset(xrf.param, 
      parameter == "Date")$regexp, smpcontents) == 1)], ":")[[1]][2]), 
         "\\t")[[1]][1], 1, 8), n_rowsdata)
   
   data.long[, subset(xrf.param, parameter == "Method")$parameter] <-
   rep(strsplit(gsub("^\\t", "", strsplit(smpcontents[which(regexpr(subset(xrf.param, 
      parameter == "Method")$regexp, smpcontents) == 1)], ":")[[1]][2]), 
         "\\t")[[1]][1], n_rowsdata)
   
   data.long[, subset(xrf.param, parameter == "Job")$parameter] <-
   rep(strsplit(gsub("^\\t", "", strsplit(smpcontents[which(regexpr(subset(xrf.param, 
      parameter == "Job")$regexp, smpcontents) == 1)], ":")[[1]][2]), 
         "\\t")[[1]][1], n_rowsdata)
   
   data.long[, subset(xrf.param, parameter == "Status")$parameter] <-
   rep(strsplit(gsub("^\\t", "", strsplit(smpcontents[which(regexpr(subset(xrf.param, 
      parameter == "Status")$regexp, smpcontents) == 1)], ":")[[1]][2]), 
         "\\t")[[1]][1], n_rowsdata)
   
   data.long[, subset(xrf.param, parameter == "Description")$parameter] <-
   rep(gsub("^\\t", "", strsplit(smpcontents[which(regexpr(subset(xrf.param, 
      parameter == "Description")$regexp, smpcontents) == 1)], ":")[[1]][2]), 
      n_rowsdata)
   
   data.long[, subset(xrf.param, parameter == "Voltage")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "Voltage")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)

   data.long[, subset(xrf.param, parameter == "Current")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "Current")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)
   
   data.long[, subset(xrf.param, parameter == "Target")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "Target")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)

   data.long[, subset(xrf.param, parameter == "Duration")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "Duration")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)

   data.long[, subset(xrf.param, parameter == "Impulse")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "Impulse")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)

   data.long[, subset(xrf.param, parameter == "DeadTime")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "DeadTime")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)

   data.long[, subset(xrf.param, parameter == "FirstChannel")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "FirstChannel")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)

   data.long[, subset(xrf.param, parameter == "LastChannel")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "LastChannel")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)

   data.long[, subset(xrf.param, parameter == "PeakTime")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "PeakTime")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)

   data.long[, subset(xrf.param, parameter == "Gain")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "Gain")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)

   data.long[, subset(xrf.param, parameter == "ZeroPeak")$parameter] <-
      rep(as.numeric(strsplit(strsplit(smpcontents[which(regexpr(subset(xrf.param, 
         parameter == "ZeroPeak")$regexp, smpcontents, perl = TRUE) == 1)], 
            ":\\t")[[1]][2], "\\t")[[1]])[1], n_rowsdata)

   
   # Convert channel into energy scale
   # Using the following assumptions:
   #    1. Zero peak is always the strongest (highest) peak in the spectrum
   # The channel with maximum counts should correspond to 0 keV
   # This gives a one-channel deviation from what the instrument shows
   # for a 12.5 keV range measurement using 1024 channels (so far)
   # This is good enough for our purposes, since the peak energies for most
   # ions do not match with reference values without a correction term anyway.
   max.channel <- which(data.long$counts == max(data.long$counts))
   data.long$energy <- (data.long$channel * (data.long$Gain / data.long$LastChannel)) -
                       ((max.channel / data.long$LastChannel) * data.long$Gain)
   # Save the maxchannel to the returned dataframe
   data.long$ZeroChannel <- rep(max.channel, n_rowsdata)
   
   # Calculate energy from channel # this is no longer viable
   #data.long$energy <- (data.long$channel * (data.long$Gain / data.long$LastChannel)) -
   #   ((24 / data.long$LastChannel) * data.long$Gain)
                                                                        
   return(data.long)
}
