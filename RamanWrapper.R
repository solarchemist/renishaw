RamanWrapper <- 
   function(data.exp, run, override = FALSE, 
            kerpk = 1, fitmaxiter = 50, gam = 0.6, scl.factor = 0.1) { 
      
   print("... Started RamanWrapper")
      
   # check if Ramanpk has already completed successfully for the current job
   current.dirname <- getwd()
   print(current.dirname)
   current.filename <- "raman-peak-data.rda"
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
      # File does not exist, OR override is TRUE
      
      print("... Started else-clause 1")

      # If file does not exist at all, run all necessary code to re-create it
      if (!file.exists(ramandatafile)) {
         ramres <- list()
         print("... ramres list created")
         
         ramres[[run]] <- 
            Ramanpk(data.exp, 
                    kerpk = kerpk, 
                    fitmaxiter = fitmaxiter, 
                    gam = gam, 
                    scl.factor = scl.factor)
                     # add maxwdth arg?
         
         save(ramres, file = ramandatafile)
      } else {
         # File already exists, but override IS TRUE
         load(file = ramandatafile)
         
         ramres[[run]] <- 
            Ramanpk(data.exp, 
                    kerpk = kerpk, 
                    fitmaxiter = fitmaxiter, 
                    gam = gam, 
                    scl.factor = scl.factor)
                    # add maxwdth arg?

         save(ramres, file = ramandatafile)
      }
      
      print("... Ended else-clause 1")
      
      return(ramres)
   }
}
