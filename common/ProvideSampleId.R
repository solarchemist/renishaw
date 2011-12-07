ProvideSampleId <- function (pathexpfile) {
   # Returns a "unique" sample ID when supplied 
   # with a path to an experimental file.
   ## Note: the sample ID must derive directly from the file or path.
   
   sampleid <- sub("\\.[\\w]+$", "", basename(pathexpfile), perl = TRUE)
   
   #### The code below is the old ProvideSampleId() function   
   #  ### OBS! Only very rudimentary error-checking.
   #  ### If the filename is formatted as \w*-\w*-\w*, we use the middle segment, 
   #  ### otherwise we use the whole string (excluding the extension)
   #  # Extract the name of the parent directory of the datafilename argument
   #  substrateid <- basename(dirname(fullpathwithfilename))
   #  # Extract the name of the method from the filename-part
   #  # First split the filename over all hyphens
   #  nameparts <- strsplit(basename(fullpathwithfilename), "-")[[1]]
   #  # If the number of nameparts exceed 3, save the whole filename
   #  # as methodid, otherwise use the middle part
   #  if (length(nameparts) > 3) {
   #     # We need to lose the file extension from the last namepart
   #     nameparts[length(nameparts)] <- 
   #        strsplit(nameparts[length(nameparts)], "\\.")[[1]][1]
   #     methodid <- paste(nameparts, collapse = "-")
   #  } else {
   #     methodid <- nameparts[2]
   #  }
   #  # Make an informative sampleid
   #  sampleid <- paste(substrateid, methodid, sep = "-")
   #  #
   #  return(sampleid)
   ####
   
   return(sampleid)
}
