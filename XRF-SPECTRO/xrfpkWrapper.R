xrfpkWrapper <- 
   function(data.exp, run, override = FALSE, 
            kerpk = 1, fitmaxiter = 100, gam = 0.64, scl.factor = 0.06, maxwdth = 10) { 
          # kerpk = 1, fitmaxiter = 50, gam = 0.6, scl.factor = 0.1) {
   # the override flag is IN USE
      
   print("... Started xrfpkWrapper")
      
   # check if xrfpk has already completed successfully for the current job
   current.dirname <- getwd()
   print(current.dirname)
   current.filename <- "xrf-peak-data.rda"
   xrfdatafile <- paste(current.dirname, current.filename, sep = "/")
   
   
   
   if (file.exists(xrfdatafile) && !override) {
      print("... Started if-clause 1")
      
      # File already exists
      # return the data using load() or data()
      
      load(file = xrfdatafile)
            
      if (run > length(xrfres)) {
         
         print("... Started if-clause 1:1")
         
         # the it does not really exist
         xrfres[[run]] <- xrfpk(data.exp, 
                                kerpk = kerpk, 
                                fitmaxiter = fitmaxiter, 
                                gam = gam, 
                                scl.factor = scl.factor,
                                maxwdth = maxwdth)
         save(xrfres, file = xrfdatafile)
         
         print("... Ended if-clause 1:1")
      }
      
      print("... Ended if-clause 1")
      
      return(xrfres)
   } else {
      
      print("... Started else-clause 1")
      
      if (!exists("xrfres")) {
         xrfres <- list()
         print("... xrfres list created")
      }   
      
      # Need to call xrfpk() and save its results to file as above
      xrfres[[run]] <- xrfpk(data.exp, 
                             kerpk = kerpk, 
                             fitmaxiter = fitmaxiter, 
                             gam = gam, 
                             scl.factor = scl.factor,
                             maxwdth = maxwdth)
      save(xrfres, file = xrfdatafile)
      
      print("... Ended else-clause 1")
      
      return(xrfres)
   }
      
}