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
