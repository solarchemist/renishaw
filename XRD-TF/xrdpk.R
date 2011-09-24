xrdpk <- 
   function(xrd.exp, kerpk = 1, fitmaxiter = 50, gam = 1.0, scl.factor = 1.2, maxwdth=5.0) {
   
   xrd.base <- baselinefit(xrd.exp, tau=2.5, gam=gam, scl.factor=scl.factor, maxwdth=maxwdth)
   
   # This loop deals with the output from baselinefit()
   # It makes a "melted" dataframe in long form for each 
   # separated peak for some baseline parameters
   xrd.pks      <- data.frame()
   xrd.pks.basl <- data.frame()
   xrd.pks.pmg  <- data.frame()
   xrd.pks.spl  <- data.frame()
   peaks <- 1:length(xrd.base$npks)
   for (s in peaks) {
      # recorded data in long form by separated peak
      xrd.pks <- rbind(xrd.pks, # column names assigned after loop
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            xrd.exp[xrd.base$indlsep[s]:xrd.base$indrsep[s], ]))
      # the calculated baseline in long form by separated peak
      xrd.pks.basl <- rbind(xrd.pks.basl, 
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = xrd.exp[xrd.base$indlsep[s]:xrd.base$indrsep[s], ],
            y = xrd.base$baseline$basisl[xrd.base$indlsep[s]:xrd.base$indrsep[s]]))
      # the taut string estimation in long form by separated peak
      xrd.pks.pmg <- rbind(xrd.pks.pmg,
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = xrd.exp[xrd.base$indlsep[s]:xrd.base$indrsep[s], ],
            y = xrd.base$pmg$fn[xrd.base$indlsep[s]:xrd.base$indrsep[s]]))
      # the weighted smoothed spline in long form by separated peak
      xrd.pks.spl <- rbind(xrd.pks.spl,
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = xrd.exp[xrd.base$indlsep[s]:xrd.base$indrsep[s], ],
            y = xrd.base$spl$reg[xrd.base$indlsep[s]:xrd.base$indrsep[s]]))
   }
   # Column names assigned to d.pks
   names(xrd.pks) <- c("peak", "kernel", "x", "y")
   
   
   # This loop calls pkdecompint() on each peak separately
   # It makes a "melted" dataframe in long form for:
   xrd.fit       <- list()       # holds pkdecompint output
   xrd.fit.fitpk <- data.frame() # contains fitting curves
   xrd.fit.parpk <- data.frame() # physical parameters by peak and kernel
   xrd.nobasl    <- data.frame() # data with baseline removed
   peaks <- 1:length(xrd.base$npks)
   for (s in peaks) {
      ######## PKDECOMPINT ########
      if (length(kerpk) > 1) {
         # set number of kernels per peak manually
         xrd.fit[[s]] <- pkdecompint(xrd.base, intnum = s,
            k = kerpk[s], maxiter = fitmaxiter)
      } else {
         # use number of kernels determined by baselinefit()
         xrd.fit[[s]] <- pkdecompint(xrd.base, intnum = s,
            k = xrd.base$npks[s], maxiter = fitmaxiter)
      }
      # Setup the dataframe that makes up the peak table
      for (kernel in 1:xrd.fit[[s]]$num.ker) {
         xrd.fit.parpk <- rbind(xrd.fit.parpk,
            data.frame(peak   = factor(xrd.fit[[s]]$intnr),
                       kernel = factor(kernel),
                       x      = xrd.fit[[s]]$parpks[kernel, "loc"],
                       height = xrd.fit[[s]]$parpks[kernel, "height"],
                       area   = xrd.fit[[s]]$parpks[kernel, "intens"],
                       fwhm   = xrd.fit[[s]]$parpks[kernel, "FWHM"],
                       beta   = xrd.fit[[s]]$parpks[kernel, "intens"] / xrd.fit[[s]]$parpks[kernel, "height"],
                       m      = xrd.fit[[s]]$parpks[kernel, "m"],
                       accept = xrd.fit[[s]]$accept))
         xrd.fit.fitpk <- rbind(xrd.fit.fitpk, 
            data.frame(peak = factor(peaks[s]),
                       kernel = factor(kernel),
                       x = xrd.fit[[s]]$x,
                       y = xrd.fit[[s]]$fitpk[kernel, ]))
      }
      xrd.nobasl <- rbind(xrd.nobasl,
         data.frame(peak = factor(peaks[s]),
                    x = xrd.fit[[s]]$x,
                    y = xrd.fit[[s]]$y))
   }
   
   
   
   return(list(xrd.base = xrd.base, 
      xrd.peaks = xrd.pks, 
      xrd.fit.parpk = xrd.fit.parpk, 
      xrd.fit.fitpk = xrd.fit.fitpk,
      xrd.nobasl = xrd.nobasl))
  
}
