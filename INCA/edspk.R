edspk <- 
   function(eds.exp, kerpk = 1, fitmaxiter = 50, gam = 1.0, scl.factor = 0.1, maxwdth=0.20) {
   
   eds.base <- baselinefit(eds.exp, tau=2.0, gam=gam, scl.factor=scl.factor, maxwdth=maxwdth)
   
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
