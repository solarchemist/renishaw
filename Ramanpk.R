Ramanpk <- 
   function(data.exp, 
            override = FALSE, 
            kerpk = 1, 
            fitmaxiter = 50, 
            gam = 0.6, 
            scl.factor = 0.1,
            tau = 2.0,
            maxwdth = 200) {
   
   print("... Starting baseline fitting")
      
   data.basl <- baselinefit(data.exp, 
                            tau = tau, 
                            gam = gam, 
                            scl.factor = scl.factor, 
                            maxwdth = maxwdth)
   
   print("... Ended baseline fitting")
   
   # This loop deals with the output from baselinefit()
   # It makes a "melted" dataframe in long form for each 
   # separated peak for some baseline parameters
   data.pks      <- data.frame()
   data.pks.basl <- data.frame()
   data.pks.pmg  <- data.frame()
   data.pks.spl  <- data.frame()
   peaks <- 1:length(data.basl$npks)
   
   
   for (s in peaks) {
      # recorded data in long form by separated peak
      data.pks <- rbind(data.pks, # column names assigned after loop
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            data.exp[data.basl$indlsep[s]:data.basl$indrsep[s], ]))
      # the calculated baseline in long form by separated peak
      data.pks.basl <- rbind(data.pks.basl, 
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = data.exp[data.basl$indlsep[s]:data.basl$indrsep[s], 1],
            y = data.basl$baseline$basisl[data.basl$indlsep[s]:data.basl$indrsep[s]]))
      # the taut string estimation in long form by separated peak
      data.pks.pmg <- rbind(data.pks.pmg,
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = data.exp[data.basl$indlsep[s]:data.basl$indrsep[s], 1],
            y = data.basl$pmg$fn[data.basl$indlsep[s]:data.basl$indrsep[s]]))
      # the weighted smoothed spline in long form by separated peak
      data.pks.spl <- rbind(data.pks.spl,
         data.frame(peak = factor(peaks[s]),
            kernel = NA,
            x = data.exp[data.basl$indlsep[s]:data.basl$indrsep[s], 1],
            y = data.basl$spl$reg[data.basl$indlsep[s]:data.basl$indrsep[s]]))
   }
   
   
   
   # Column names assigned to d.pks
   names(data.pks) <- c("peak", "kernel", "x", "y")
   
   
   # This loop calls pkdecompint() on each peak separately
   # It makes a "melted" dataframe in long form for:
   data.fit       <- list()       # holds pkdecompint output
   data.fit.fitpk <- data.frame() # contains fitting curves
   data.fit.parpk <- data.frame() # physical parameters by peak and kernel
   data.fit.basl  <- data.frame() # data with baseline removed
   peaks <- 1:length(data.basl$npks)
   for (s in peaks) {
      ######## PKDECOMPINT ########
      if (length(kerpk) > 1) {
         # set number of kernels per peak manually
         data.fit[[s]] <- pkdecompint(data.basl, intnum = s,
            k = kerpk[s], maxiter = fitmaxiter)
      } else {
         # use number of kernels determined by baselinefit()
         data.fit[[s]] <- pkdecompint(data.basl, intnum = s,
            k = data.basl$npks[s], maxiter = fitmaxiter)
      }
      # Setup the dataframe that makes up the peak table
      for (kernel in 1:data.fit[[s]]$num.ker) {
         data.fit.parpk <- rbind(data.fit.parpk,
            data.frame(peak   = factor(data.fit[[s]]$intnr),
                       kernel = factor(kernel),
                       x      = data.fit[[s]]$parpks[kernel, "loc"],
                       height = data.fit[[s]]$parpks[kernel, "height"],
                       area   = data.fit[[s]]$parpks[kernel, "intens"],
                       fwhm   = data.fit[[s]]$parpks[kernel, "FWHM"],
                       m      = data.fit[[s]]$parpks[kernel, "m"],
                       accept = data.fit[[s]]$accept))
         data.fit.fitpk <- rbind(data.fit.fitpk, 
            data.frame(peak = factor(peaks[s]),
                       kernel = factor(kernel),
                       x = data.fit[[s]]$x,
                       y = data.fit[[s]]$fitpk[kernel, ]))
      }
      data.fit.basl <- rbind(data.fit.basl,
         data.frame(peak = factor(peaks[s]),
                    x = data.fit[[s]]$x,
                    y = data.fit[[s]]$y))
   }
   
   
   
   return(list(data.basl = data.basl, 
      data.peaks = data.pks, 
      data.fit.parpk = data.fit.parpk, 
      data.fit.fitpk = data.fit.fitpk,
      data.fit.basl = data.fit.basl))
}
