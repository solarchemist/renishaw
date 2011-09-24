##################################################
################## matchpdf ######################
##################################################
matchpdf <- function(expcol, pdfrow) {
   # Function for matching each element of two numeric vectors
   # to the closest-lying element (numerically) in the other
   # vector. For example for matching 2th values of PDF to
   # 2th values of recorded diffractogram.
   ## Args: two vectors to be matched (any special considerations on vector lengths?)
   ## Values: A list of 5 numeric vectors. See description inside function.
   # A word of caution. This function creates a matrix with dimensions pdfrow * expcol.
   # Usually both of these vectors are rather short. But watch out if large vectors
   # are submitted --- the looping below could make this function very slow.
#
   diff.thth <- matrix(NA, length(pdfrow), length(expcol))
   diff.indx <- matrix(NA, length(pdfrow), length(expcol))
   
   for (pdf.row in 1:length(pdfrow)) {
      # For every row, calculate differences between that pdf value and all exp values
      diff.thth[pdf.row, ] <- expcol - pdfrow[pdf.row]
   }
#   
   # Now we go on to find the minimum along each row of diff.thth
   # and set all other values to some arbitrary value
   diff.rmin <- diff.thth
   for (pdf.row in 1:dim(diff.thth)[1]) {
      for (exp.col in 1:dim(diff.thth)[2]) {
         if (abs(diff.thth[pdf.row, exp.col]) != min(abs(diff.thth[pdf.row, ]))) {
            diff.rmin[pdf.row, exp.col] <- Inf
         }
      }
   }
#   
   # We likewise find the minimum along each column and set any other,
   # non-Inf values to Inf
   diff.cmin <- diff.rmin
   for (exp.col in 1:dim(diff.rmin)[2]) {
      for (pdf.row in 1:dim(diff.rmin)[1]) {
         if (abs(diff.rmin[pdf.row, exp.col]) != min(abs(diff.rmin[, exp.col]))) {
            diff.cmin[pdf.row, exp.col] <- Inf
         }
      }
   }
#   
   # The matrix now contains at most one non-Inf value per row and column.
   # To make the use of colSums and rowSums possible (next step) all Inf are set zero,
   # and all matches are set to 1.
   for (pdf.row in 1:dim(diff.cmin)[1]) {
      for (exp.col in 1:dim(diff.cmin)[2]) {
         if ((diff.cmin[pdf.row, exp.col]) != Inf) {
            diff.indx[pdf.row, exp.col] <- 1
         } else {
            diff.indx[pdf.row, exp.col] <- 0
         }
      }
   }
#   
   # Attach error message to return list if sum(rowSums(diff.indx)) == sum(colSums(diff.indx))
   matchpdf.err.msg <-"matchpdf() failed to complete. It seems rowSums != colSums."
   mtch <- matchpdf.err.msg
#   
   # Check that sum(rowSums(diff.indx)) == sum(colSums(diff.indx))
   if (sum(rowSums(diff.indx)) == sum(colSums(diff.indx))) {
      # Reset mtch
      mtch <- list()
      mtch <- list(csums = colSums(diff.indx), 
                   rsums = rowSums(diff.indx), 
                   expthth = expcol[colSums(diff.indx) != 0], 
                   pdfthth = pdfrow[rowSums(diff.indx) != 0], 
                   deltathth = expcol[colSums(diff.indx) != 0] - pdfrow[rowSums(diff.indx) != 0])
      # List of 5
      # $ csums     : num - consisting of ones and zeroes. Shows you which positions of expcol matched.
      # $ rsums     : num - consisting of ones and zeroes. Shows you which positions of pdfrow matched.
      # $ expthth   : num - consisting of the matching elements of expcol.
      # $ pdfthth   : num - consisting of the matching elements of pdfrow.
      # $ deltathth : num - element-wise difference of expcol and pdfrow.
   }
#
   return(mtch)
}
