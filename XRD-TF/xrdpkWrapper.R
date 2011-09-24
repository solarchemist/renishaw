xrdpkWrapper <- 
   function(data.exp, run, override = FALSE, 
            kerpk = 1, fitmaxiter = 50, gam = 1.0, scl.factor = 1.2, maxwdth=5.0) { 
      
   print("... Started xrdpkWrapper")
      
   # check if xrdpk has already completed successfully for the current job
   current.dirname <- getwd()
   print(current.dirname)
   current.filename <- "xrd-peak-data.rda"
   xrddatafile <- paste(current.dirname, current.filename, sep = "/")
   
   
   if (file.exists(xrddatafile) && !override) {
      print("... Started if-clause 1")
      
      # File already exists
      # return the data using load() or data()
      
      load(file = xrddatafile)
            
      if (run > length(xrdres)) {
         
         print("... Started if-clause 1:1")
         
         # then it does not really exist
         xrdres[[run]] <- xrdpk(data.exp, 
                                kerpk = kerpk, 
                                fitmaxiter = fitmaxiter, 
                                gam = gam, 
                                scl.factor = scl.factor,
                                maxwdth = maxwdth)
         save(xrdres, file = xrddatafile)
         
         print("... Ended if-clause 1:1")
      }
      
      print("... Ended if-clause 1")
      
      return(xrdres)
   } else {
      # File does not exist
      # OR override is TRUE
      
      print("... Started else-clause 1")
      
      # If file does not exist at all, run all necessary code to re-create it
      if (!file.exists(xrddatafile)) {
         xrdres <- list()
         print("... xrdres list created")
         
         xrdres[[run]] <- 
            xrdpk(data.exp, 
                  kerpk = kerpk, 
                  fitmaxiter = fitmaxiter, 
                  gam = gam, 
                  scl.factor = scl.factor,
                  maxwdth = maxwdth)
         
         save(xrdres, file = xrddatafile)
      } else {
         # File already exists, but override is TRUE
         load(file = xrddatafile)
         
         xrdres[[run]] <- 
            xrdpk(data.exp, 
                  kerpk = kerpk, 
                  fitmaxiter = fitmaxiter, 
                  gam = gam, 
                  scl.factor = scl.factor,
                  maxwdth = maxwdth)
         
         save(xrdres, file = xrddatafile)
      }
      
      print("... Ended else-clause 1")
      
      return(xrdres)
   }     
}
