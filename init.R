# To source a bunch of files in the same directory

sourceDir <- function(path, trace = TRUE) {
   lsDir <- list.files(path, pattern = "\\.[Rr]$")
   for (i in lsDir) {
      if(trace) {cat(i, ":")}
      source(file.path(path, i))
      if(trace) {cat("\n")}
   }
}
