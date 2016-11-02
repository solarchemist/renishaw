#' Raman peak fit
#'
#' This is the actual workhorse that does the peak fitting for each spectra.
#' Note: this function only works for a single spectra and does not handle
#'       saving the output and so on. Always use one of the wrapper functions.
#'
#' @param data.exp     dataframe with experimental data (format?)
#' @param override     force re-running of peak analysis
#' @param kerpk        number of kernels per peak (passed to Ramanpk())
#' @param fitmaxiter   number of max iterations while attempting fit (passed to Ramanpk())
#' @param gam          gam (passed to Ramanpk())
#' @param scl.factor   sclerosis factor (passed to Ramanpk())
#' @param tau          tau (passed to Ramanpk())
#' @param maxwdth      peak max width (passed to Ramanpk())
#'
#' @return a list of 4 elements (lists or dataframes)
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

      data.basl <-
         diffractometry::baselinefit(data.exp,
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
            data.fit[[s]] <-
               diffractometry::pkdecompint(data.basl, intnum = s,
                                           k = kerpk[s], maxiter = fitmaxiter)
         } else {
            # use number of kernels determined by baselinefit()
            data.fit[[s]] <-
               diffractometry::pkdecompint(data.basl, intnum = s,
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



#' Use this wrapper function for calling Ramanpk
#'
#' I wanted to run the \pkg{diffractometry} peak-analysis code (which can be both
#' time-consuming and prone to stray errors halting compilation) without the need
#' to re-run the peak analysis every time the document is re-compiled.
#' For that purpose, I created this wrapper function.
#' Its job is to call \code{\link{Ramanpk}} only if necessary (i.e., peak analysis
#' has not already been performed). That is made possible by saving the results
#' of a successful analysis to a \code{raman-peak-data.rda} file in the directory
#' of the report.
#' Also includes an override option to force re-running of peak analysis.
#'
#' @param data.exp     dataframe with experimental data (format?)
#' @param run          index
#' @param override     force re-running of peak analysis
#' @param kerpk        number of kernels per peak (passed to Ramanpk())
#' @param fitmaxiter   number of max iterations while attempting fit (passed to Ramanpk())
#' @param gam          gam (passed to Ramanpk())
#' @param scl.factor   sclerosis factor (passed to Ramanpk())
#' @param tau          tau (passed to Ramanpk())
#' @param maxwdth      peak max width (passed to Ramanpk())
#'
#' @return a Ramanpk dataframe
#' @export
RamanWrapper <-
   function(data.exp,
            run,
            override = FALSE,
            kerpk = 1,
            fitmaxiter = 50,
            gam = 0.6,
            scl.factor = 0.1,
            tau = 2.0,
            maxwdth = 200) {
      # the override flag IS IN USE

      print("... Started RamanWrapper")

      # check if Ramanpk has already completed successfully for the current job
      current.dirname <- getwd()
      print(current.dirname)
      current.filename <- "raman-peak-data.rda"
      ramandatafile <- paste(current.dirname, current.filename, sep = "/")


      # What follows are three if-clauses (containing no else-statements).
      # We can be in one of three states:
      # - previous peak data exists, and override is not requested
      # - - in this case we check the run number vs the length of
      # - - the previously saved data. If the run number is larger
      # - - we run a peak analysis since it won't overwrite anything
      # - previous peak data exists, and override is requested
      # - no previous peak data exists

      if (file.exists(ramandatafile) && !override) {
         # If file does exist AND override flag is FALSE
         print("... Started if-clause 1")

         # File already exists
         # return the data using load() or data()
         # Load the existing data from file
         load(file = ramandatafile)

         # Only run the peak-fitting algorithm if
         # <run> is higher than what the file contains
         if (run > length(ramres)) {
            print("... Started if-clause 1:1")
            ramres[[run]] <- Ramanpk(data.exp,
                                     kerpk = kerpk,
                                     fitmaxiter = fitmaxiter,
                                     gam = gam,
                                     scl.factor = scl.factor,
                                     tau = tau,
                                     maxwdth = maxwdth)
            save(ramres, file = ramandatafile)
            print("... Ended if-clause 1:1")
         }

         print("... Ended if-clause 1")
      }
      if (file.exists(ramandatafile) && override) {
         # If file does exist AND override flag is TRUE
         print("... Started if-clause 2")

         # Load the existing data from file
         load(file = ramandatafile)

         ramres[[run]] <- Ramanpk(data.exp,
                                  kerpk = kerpk,
                                  fitmaxiter = fitmaxiter,
                                  gam = gam,
                                  scl.factor = scl.factor,
                                  tau = tau,
                                  maxwdth = maxwdth)
         save(ramres, file = ramandatafile)
         print("... Ended if-clause 2")
      }
      # If the file does not exist,
      # it doesn't really matter what the override flag says
      if (!file.exists(ramandatafile)) {
         print("... Started if-clause 3")

         ramres <- list()
         print("... ramres list created")

         ramres[[run]] <- Ramanpk(data.exp,
                                  kerpk = kerpk,
                                  fitmaxiter = fitmaxiter,
                                  gam = gam,
                                  scl.factor = scl.factor,
                                  tau = tau,
                                  maxwdth = maxwdth)
         save(ramres, file = ramandatafile)
         print("... Ended if-clause 3")
      }
      return(ramres)
   }





#' Use this wrapper function for silicon spectra
#'
#' This is a more restricted wrapper, specifically intended for
#' use with silicon reference spectra.
#'
#' @param data.exp     dataframe with experimental data (format?)
#' @param run          index
#' @param override     force re-running of peak analysis
#' @param kerpk        number of kernels per peak
#' @param fitmaxiter   number of max iterations while attempting fit (diffractometry option)
#' @param gam          gam (diffractometry option)
#' @param scl.factor   sclerosis factor (diffractometry option)
#'
#' @return a Ramanpk dataframe
#' @export
SiliconWrapper <-
   function(data.exp, run, override = FALSE,
            kerpk = 1, fitmaxiter = 50, gam = 0.6, scl.factor = 0.1) {
      # the override flag is currently not used

      print("... Started SiliconWrapper")

      # check if Ramanpk has already completed successfully for the current job
      current.dirname <- getwd()
      print(current.dirname)
      current.filename <- "silicon-peak-data.rda"
      ramandatafile <- paste(current.dirname, current.filename, sep = "/")



      if (file.exists(ramandatafile) && !override) {
         print("... Started if-clause 1")

         # File already exists
         # return the data using load() or data()

         load(file = ramandatafile)

         if (run > length(ramres)) {

            print("... Started if-clause 1:1")

            # the it does not really exist
            ramres[[run]] <- Ramanpk(data.exp,
                                     kerpk = kerpk,
                                     fitmaxiter = fitmaxiter,
                                     gam = gam,
                                     scl.factor = scl.factor)
            save(ramres, file = ramandatafile)

            print("... Ended if-clause 1:1")
         }

         print("... Ended if-clause 1")

         return(ramres)
      } else {
         # File does not exist, or override requested

         print("... Started else-clause 1")

         if (!exists("ramres")) {
            ramres <- list()
            print("... ramres list created")
         }

         # Need to call Ramanpk() and save its results to file as above
         ramres[[run]] <- Ramanpk(data.exp,
                                  kerpk = kerpk,
                                  fitmaxiter = fitmaxiter,
                                  gam = gam,
                                  scl.factor = scl.factor)
         save(ramres, file = ramandatafile)

         print("... Ended else-clause 1")

         return(ramres)
      }

   }
