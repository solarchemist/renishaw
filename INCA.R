# INCA.R
# Functions to read and manipulate data from the Oxford INCA EDS
# Taha Ahmed, April 2011

# CONTENTS
source("/home/taha/chepec/chetex/common/R/common.R")
# >>>> eds2df
# >>>> edspk




##################################################
################### eds2df #######################
##################################################
eds2df <- function(edstxtfile) {
   ## Description:
   ##   Reads EDS textfile from INCA EDS.
   ##   Stores data in data frame.
   ## Usage:
   ##   eds2df(edstxtfile)
   ## Arguments:
   ##    edstxtfile: character string, the full filename
   ##                (with path) to one txt file.
   ## Value:
   ##   A dataframe 
   #
   incatxt <- file(edstxtfile, "r")
   edsfile <- readLines(incatxt, n = -1) #read all lines of input file
   close(incatxt)
   #
   sampleid <- ProvideSampleId(edstxtfile)
   #
   rgxp.comment <- "^\\#"
   #
   numrow.idx <- regexpr(rgxp.comment, edsfile)
   # scrap the match.length attribute
   attr(numrow.idx, "match.length") <- NULL
   #
   i <- seq(1, length(numrow.idx) - 1, 1)
   j <- seq(2, length(numrow.idx), 1)
   # Start index of data range
   start.idx <- which(numrow.idx[i] == 1 & numrow.idx[j] != 1) + 1
   # End index of the data range
   end.idx <- which(numrow.idx[i] != 1 & numrow.idx[j] == 1)
   #
   zz <- textConnection(edsfile[start.idx:end.idx], "r")
   #
   ff <- data.frame()
   ff <- data.frame(stringsAsFactors = FALSE,
            sampleid = sampleid, 
            matrix(scan(zz, what = numeric(), sep = ","), ncol = 2, byrow = T))
   close(zz)
   names(ff) <- c("sampleid", "energy", "counts")
   #
   ### Collect attributes of this experiment
   # XUnit
   position.XUnit <- regexpr("^\\#XUNITS", edsfile)
   XUnit <- as.character(strsplit(edsfile[which(position.XUnit == 1)], ":\\s")[[1]][2])
   ff$XUnit <- XUnit
   # YUnit
   position.YUnit <- regexpr("^\\#YUNITS", edsfile)
   YUnit <- as.character(strsplit(edsfile[which(position.YUnit == 1)], ":\\s")[[1]][2])
   ff$YUnit <- YUnit
   # Date
   position.Date <- regexpr("^\\#DATE", edsfile)
   Date <- strsplit(edsfile[which(position.Date == 1)], ":\\s")[[1]][2]
   ff$Date <- Date
   # Time
   position.Time <- regexpr("^\\#TIME", edsfile)
   Time <- strsplit(edsfile[which(position.Time == 1)], ":\\s")[[1]][2]
   ff$Time <- Time
   # XPerChannel
   position.XPerChannel <- regexpr("^\\#XPERCHAN", edsfile)
   XPerChannel <- as.numeric(strsplit(edsfile[which(position.XPerChannel == 1)], ":\\s")[[1]][2])
   ff$XPerChannel <- XPerChannel
   # Offset
   position.Offset <- regexpr("^\\#OFFSET", edsfile)
   Offset <- as.numeric(strsplit(edsfile[which(position.Offset == 1)], ":\\s")[[1]][2])
   ff$Offset <- Offset
   # ChOffset
   position.ChOffset <- regexpr("^\\#CHOFFSET", edsfile)
   ChOffset <- as.numeric(strsplit(edsfile[which(position.ChOffset == 1)], ":\\s")[[1]][2])
   ff$ChOffset <- ChOffset
   # LiveTime
   position.LiveTime <- regexpr("^\\#LIVETIME", edsfile)
   LiveTime <- as.numeric(strsplit(edsfile[which(position.LiveTime == 1)], ":\\s")[[1]][2])
   ff$LiveTime <- LiveTime
   # DeadTime is calculated from: REALTIME - LIVETIME
   position.RealTime <- regexpr("^\\#REALTIME", edsfile)
   RealTime <- as.numeric(strsplit(edsfile[which(position.RealTime == 1)], ":\\s")[[1]][2])
   DeadTime <- RealTime - LiveTime
   ff$DeadTime <- DeadTime
   # BeamEnergy
   position.BeamEnergy <- regexpr("^\\#BEAMKV", edsfile)
   BeamEnergy <- as.numeric(strsplit(edsfile[which(position.BeamEnergy == 1)], ":\\s")[[1]][2])
   ff$BeamEnergy <- BeamEnergy   
   #
   return(ff)
}




##################################################
#################### edspk #######################
##################################################
edspk <- function(eds.exp, kerpk = 1, fitmaxiter = 50) {
   
   eds.base <- baselinefit(eds.exp, tau=2.0, gam=1.0, scl.factor=3.0, maxwdth=0.20)
   
   # This loop deals with the output from baselinefit()
   # It makes a "melted" dataframe in long form for each 
   # separated peak for some baseline parameters
   eds.pks      <- data.frame()
   eds.pks.basl <- data.frame()
   eds.pks.pmg  <- data.frame()
   eds.pks.spl  <- data.frame()
   peaks <- 1:length(eds.base$npks)
   for (s in peaks) {
      # recorded data in long form by separated peak
      eds.pks <- rbind(eds.pks, # column names assigned after loop
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            eds.exp[eds.base$indlsep[s]:eds.base$indrsep[s], ]))
      # the calculated baseline in long form by separated peak
      eds.pks.basl <- rbind(eds.pks.basl, 
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = eds.exp[eds.base$indlsep[s]:eds.base$indrsep[s]],
            y = eds.base$baseline$basisl[eds.base$indlsep[s]:eds.base$indrsep[s]]))
      # the taut string estimation in long form by separated peak
      eds.pks.pmg <- rbind(eds.pks.pmg,
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = eds.exp[eds.base$indlsep[s]:eds.base$indrsep[s]],
            y = eds.base$pmg$fn[eds.base$indlsep[s]:eds.base$indrsep[s]]))
      # the weighted smoothed spline in long form by separated peak
      eds.pks.spl <- rbind(eds.pks.spl,
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = eds.exp[eds.base$indlsep[s]:eds.base$indrsep[s]],
            y = eds.base$spl$reg[eds.base$indlsep[s]:eds.base$indrsep[s]]))
   }
   # Column names assigned to d.pks
   names(eds.pks) <- c("peak", "kernel", "x", "y")
   
   
   # This loop calls pkdecompint() on each peak separately
   # It makes a "melted" dataframe in long form for:
   eds.fit       <- list()       # holds pkdecompint output
   eds.fit.fitpk <- data.frame() # contains fitting curves
   eds.fit.parpk <- data.frame() # physical parameters by peak and kernel
   eds.nobasl    <- data.frame() # data with baseline removed
   peaks <- 1:length(eds.base$npks)
   for (s in peaks) {
      ######## PKDECOMPINT ########
      if (length(kerpk) > 1) {
         # set number of kernels per peak manually
         eds.fit[[s]] <- pkdecompint(eds.base, intnum = s,
            k = kerpk[s], maxiter = fitmaxiter)
      } else {
         # use number of kernels determined by baselinefit()
         eds.fit[[s]] <- pkdecompint(eds.base, intnum = s,
            k = eds.base$npks[s], maxiter = fitmaxiter)
      }
      # Setup the dataframe that makes up the peak table
      for (kernel in 1:eds.fit[[s]]$num.ker) {
         eds.fit.parpk <- rbind(eds.fit.parpk,
            data.frame(peak   = factor(eds.fit[[s]]$intnr),
                       kernel = factor(kernel),
                       x      = eds.fit[[s]]$parpks[kernel, "loc"],
                       height = eds.fit[[s]]$parpks[kernel, "height"],
                       area   = eds.fit[[s]]$parpks[kernel, "intens"],
                       fwhm   = eds.fit[[s]]$parpks[kernel, "FWHM"],
                       m      = eds.fit[[s]]$parpks[kernel, "m"],
                       accept = eds.fit[[s]]$accept))
         eds.fit.fitpk <- rbind(eds.fit.fitpk, 
            data.frame(peak = factor(peaks[s]),
                       kernel = factor(kernel),
                       x = eds.fit[[s]]$x,
                       y = eds.fit[[s]]$fitpk[kernel, ]))
      }
      eds.nobasl <- rbind(eds.nobasl,
         data.frame(peak = factor(peaks[s]),
                    x = eds.fit[[s]]$x,
                    y = eds.fit[[s]]$y))
   }
   
   
   
   return(list(eds.base = eds.base, 
      eds.peaks = eds.pks, 
      eds.fit.parpk = eds.fit.parpk, 
      eds.fit.fitpk = eds.fit.fitpk,
      eds.nobasl = eds.nobasl))
  
}














