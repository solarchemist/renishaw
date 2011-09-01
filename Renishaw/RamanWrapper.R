RamanWrapper <- 
   function(data.exp, run, override = FALSE, 
            kerpk = 1, fitmaxiter = 50, gam = 0.6, scl.factor = 0.1) { 
   # the override flag is currently not used
      
   print("... Started RamanWrapper")
      
   # check if Ramanpk has already completed successfully for the current job
   current.dirname <- getwd()
   print(current.dirname)
   current.filename <- "raman-peak-data.rda"
   ramandatafile <- paste(current.dirname, current.filename, sep = "/")
   
   
   
   if (file.exists(ramandatafile)) {
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