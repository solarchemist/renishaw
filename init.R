# To source a bunch of files in the same directory

sourceDir <- function(path, trace = TRUE) {
   for (nm in list.files(path, pattern = "\\.[Rr]$")) {
      if(trace) {
         cat(nm,":")
      }
      source(file.path(path, nm))
      if(trace) {
         cat("\n")
      }
   }
}