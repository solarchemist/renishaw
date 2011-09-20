edsWrapper <- 
   function(data.exp, run, override = FALSE, 
            kerpk = 1, fitmaxiter = 50, gam = 0.6, scl.factor = 0.1, maxwdth=0.20) { 
      
   print("... Started edsWrapper")
      
   # check if edspk has already completed successfully for the current job
   current.dirname <- getwd()
   print(current.dirname)
   current.filename <- "eds-peak-data.rda"
   edsdatafile <- paste(current.dirname, current.filename, sep = "/")
   
   
   if (file.exists(edsdatafile) && !override) {
      print("... Started if-clause 1")
      
      # File already exists
      # return the data using load() or data()
      
      load(file = edsdatafile)
            
      if (run > length(edsres)) {
         
         print("... Started if-clause 1:1")
         
         # then it does not really exist
         edsres[[run]] <- edspk(data.exp, 
                                kerpk = kerpk, 
                                fitmaxiter = fitmaxiter, 
                                gam = gam, 
                                scl.factor = scl.factor,
                                maxwdth = maxwdth)
         save(edsres, file = edsdatafile)
         
         print("... Ended if-clause 1:1")
      }
      
      print("... Ended if-clause 1")
      
      return(edsres)
   } else {
      
      print("... Started else-clause 1")
      
      if (!exists("edsres")) {
         edsres <- list()
         print("... edsres list created")
      }   
      
      # Need to call edspk() and save its results to file as above
      edsres[[run]] <- edspk(data.exp, 
                             kerpk = kerpk, 
                             fitmaxiter = fitmaxiter, 
                             gam = gam, 
                             scl.factor = scl.factor,
                             maxwdth = maxwdth)
      save(edsres, file = edsdatafile)
      
      print("... Ended else-clause 1")
      
      return(edsres)
   }     
}
