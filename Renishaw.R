# Renishaw.R
# Functions to read data from the Renishaw Raman spectrometer
# Taha Ahmed, Feb 2011

# CONTENTS
source("/home/taha/chepec/chetex/common/R/common.R")
# >>>> Raman2df
# >>>> Ramanpk





##################################################
################### Raman2df #######################
##################################################
Raman2df <- function(datafilename) {
   # Function description: for reading Raman spectrum into dataframe
   #
   datafile <- file(datafilename, "r")
   chifile <- readLines(datafile, n = -1) #read all lines of input file
   close(datafile)
   #
   #####
   sampleid <- ProvideSampleId(datafilename)
   #
   ff <- data.frame(NULL)
   zz <- textConnection(chifile, "r")
   ff <- rbind(ff, data.frame(stringsAsFactors = FALSE,
            sampleid,
            matrix(scan(zz, what = numeric(), sep = "\t"),
            ncol = 2, byrow = T)))
   close(zz)
   names(ff) <- c("sampleid", "shift", "counts")
   # Re-order by increasing shift
   ff <- ff[order(ff$shift), ]
   # And fix the row.names
   row.names(ff) <- seq(1, dim(ff)[1])
   # Re-calculate the spectrum with evenly spaced points
   # (so that the peak-finding algorithms of diffractometry package may be used)
   ff$cshift <- approx(x = ff$shift, y = ff$counts, 
      method = "linear", n = length(ff$shift))$x
   ff$ccounts <- approx(x = ff$shift, y = ff$counts, 
      method = "linear", n = length(ff$shift))$y 
   ##
   return(ff)
}


##################################################
################## Ramanpk #######################
##################################################
Ramanpk <- function(Raman.exp, kerpk = 1, fitmaxiter = 50) {
   
   Raman.base <- baselinefit(Raman.exp, tau=2.0, gam=1.0, scl.factor=1.2, maxwdth=400)
   
   # This loop deals with the output from baselinefit()
   # It makes a "melted" dataframe in long form for each 
   # separated peak for some baseline parameters
   Raman.pks      <- data.frame()
   Raman.pks.basl <- data.frame()
   Raman.pks.pmg  <- data.frame()
   Raman.pks.spl  <- data.frame()
   peaks <- 1:length(Raman.base$npks)
   for (s in peaks) {
      # recorded data in long form by separated peak
      Raman.pks <- rbind(Raman.pks, # column names assigned after loop
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            Raman.exp[Raman.base$indlsep[s]:Raman.base$indrsep[s], ]))
      # the calculated baseline in long form by separated peak
      Raman.pks.basl <- rbind(Raman.pks.basl, 
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = Raman.exp[Raman.base$indlsep[s]:Raman.base$indrsep[s]],
            y = Raman.base$baseline$basisl[Raman.base$indlsep[s]:Raman.base$indrsep[s]]))
      # the taut string estimation in long form by separated peak
      Raman.pks.pmg <- rbind(Raman.pks.pmg,
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = Raman.exp[Raman.base$indlsep[s]:Raman.base$indrsep[s]],
            y = Raman.base$pmg$fn[Raman.base$indlsep[s]:Raman.base$indrsep[s]]))
      # the weighted smoothed spline in long form by separated peak
      Raman.pks.spl <- rbind(Raman.pks.spl,
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = Raman.exp[Raman.base$indlsep[s]:Raman.base$indrsep[s]],
            y = Raman.base$spl$reg[Raman.base$indlsep[s]:Raman.base$indrsep[s]]))
   }
   # Column names assigned to d.pks
   names(Raman.pks) <- c("peak", "kernel", "x", "y")
   
   
   # This loop calls pkdecompint() on each peak separately
   # It makes a "melted" dataframe in long form for:
   Raman.fit       <- list()       # holds pkdecompint output
   Raman.fit.fitpk <- data.frame() # contains fitting curves
   Raman.fit.parpk <- data.frame() # physical parameters by peak and kernel
   Raman.nobasl    <- data.frame() # data with baseline removed
   peaks <- 1:length(Raman.base$npks)
   for (s in peaks) {
      ######## PKDECOMPINT ########
      if (length(kerpk) > 1) {
         # set number of kernels per peak manually
         Raman.fit[[s]] <- pkdecompint(Raman.base, intnum = s,
            k = kerpk[s], maxiter = fitmaxiter)
      } else {
         # use number of kernels determined by baselinefit()
         Raman.fit[[s]] <- pkdecompint(Raman.base, intnum = s,
            k = Raman.base$npks[s], maxiter = fitmaxiter)
      }
      # Setup the dataframe that makes up the peak table
      for (kernel in 1:Raman.fit[[s]]$num.ker) {
         Raman.fit.parpk <- rbind(Raman.fit.parpk,
            data.frame(peak   = factor(Raman.fit[[s]]$intnr),
                       kernel = factor(kernel),
                       x      = Raman.fit[[s]]$parpks[kernel, "loc"],
                       height = Raman.fit[[s]]$parpks[kernel, "height"],
                       area   = Raman.fit[[s]]$parpks[kernel, "intens"],
                       fwhm   = Raman.fit[[s]]$parpks[kernel, "FWHM"],
                       m      = Raman.fit[[s]]$parpks[kernel, "m"],
                       accept = Raman.fit[[s]]$accept))
         Raman.fit.fitpk <- rbind(Raman.fit.fitpk, 
            data.frame(peak = factor(peaks[s]),
                       kernel = factor(kernel),
                       x = Raman.fit[[s]]$x,
                       y = Raman.fit[[s]]$fitpk[kernel, ]))
      }
      Raman.nobasl <- rbind(Raman.nobasl,
         data.frame(peak = factor(peaks[s]),
                    x = Raman.fit[[s]]$x,
                    y = Raman.fit[[s]]$y))
   }
   
   
   
   return(list(Raman.base = Raman.base, 
      Raman.peaks = Raman.pks, 
      Raman.fit.parpk = Raman.fit.parpk, 
      Raman.fit.fitpk = Raman.fit.fitpk,
      Raman.nobasl = Raman.nobasl)) 
}
