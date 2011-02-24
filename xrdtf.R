# Collection of R functions for analysis of thin-film X-ray diffraction patterns
# Maintainer: Taha Ahmed

# CONTENTS
# >>>> matchpdf
# >>>> scherrer
# >>>> pdf2df
# >>>> uxd2mtx
# >>>> uxd2df
# >>>> muxd2df
# >>>> muxd2mtx
# >>>> muxd2ls
#    - REPAIR SHOP
#    - print.xtable.booktabs
#    - split.muxd
#    - strip.ka2
#    - pearson.beta

## All XRD reading functions need to be re-written to identify the data instead of identifying metadata
## Take inspiration from the functions in CHI.R




##################################################
################ EliminateKa2 ####################
##################################################
EliminateKa2 <- function(xrdata) {
   ##### STILL UNDER CONSTRUCTION ####
   ##### STILL UNDER CONSTRUCTION ####

   # The following lever arm weights are from Dong1999a
   weights <- list()
   # Three-bar weights
   weights[[1]] <- c(0.005134296, 0.491686047, 0.003179657)
   # Five-bar weights   
   weights[[2]] <- c(0.002614410, 0.011928014, 0.480406967, 0.002121807, 0.002928802)
   # Seven-bar weights
   weights[[3]] <- c(0.001580069, 0.003463773, 0.015533472, 0.422601977, 0.053632977,
                     0.001572467, 0.001615265)
   # Nine-bar weights
   weights[[4]] <- c(0.001138001,   0.00195272, 0.004324464, 0.019246541, 0.394175823,
                     0.079159001, -0.003591547, 0.002505604, 0.001089392)
   # 15-bar weights
   weights[[5]] <- c(0.000614225, 0.000810836, 0.001134775, 0.001723265, 0.002968405,
                     0.006433676, 0.02575384,  0.345872599, 0.100578092, 0.014493969,
                    -0.004176171, 0.000678688, 0.001610333, 0.000918077, 0.000585391)
   # 25-bar weights
   weights[[6]] <- c(0.000349669,  0.000408044, 0.000484578,  0.000587457, 0.000730087,
                     0.000935685,  0.001247401, 0.001753233,  0.002657209, 0.004531817,
                     0.009591103,  0.034998436, 0.2876498,    0.074964321, 0.065000871,
                     0.016762729, -0.00306221, -0.002717412, -0.000902322, 0.000915701,
                     0.001036484,  0.000808199, 0.000539899,  0.000398896, 0.000330325)
   
   ##### STILL UNDER CONSTRUCTION ####
}





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
      mtch <- list(csums = colSums(diff.indx), rsums = rowSums(diff.indx), expthth = expcol[colSums(diff.indx) != 0], pdfthth = pdfrow[rowSums(diff.indx) != 0], deltathth = expcol[colSums(diff.indx) != 0] - pdfrow[rowSums(diff.indx) != 0])
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




##################################################
################## scherrer ######################
##################################################
scherrer <- function(integralbreadth, thth, wavelength = 1.54056, shapeconstant = ((4/3)*(pi/6))^(1/3)) {
   # Function for calculating crystallite grain size from reflection data
   # ARGS: integralbreadth - vector with integral breadth of reflections (in degrees)
   #       thth            - vector with 2theta values of reflections (in degrees)
   #       wavelength      - X-ray wavelength used (default 1.54056 Å, Cu Ka)
   #       shapeconstant   - Scherrer constant (default spherical, ~0.9)
   # VALUE: vector with size parameters
   ## REQUIRES: as.radians(), source("/home/taha/chepec/chetex/common/R/common.R")
   D <- (shapeconstant * wavelength) / (as.radians(integralbreadth) * cos(as.radians(thth)))
   # cos() - angles must be in radians, not degrees!
   return(D) #units of angstrom
}




##################################################
################### pdf2df #######################
##################################################
pdf2df <- function(pdffile) {
   # Function for extracting information from ICDD PDF XML-files
   # For example the PDF files produced by the PDF database at Ångström's X-ray lab
   # NOTE: sometimes intensity values are specified as less than some value.
   #       In those cases, this function simply strips the less-than character.
   # ARGS: pdffile (complete path and filename to PDF file)
   # VALUE: dataframe with 9 columns:
   #        thth angles (numeric),
   #        d (numeric),
   #        h index (numeric),
   #        k index (numeric),
   #        l index (numeric),
   #        hkl indices (factor of strings),
   #        hkl.TeX indices formatted for LaTeX (factor of strings),
   #        intensity (numeric),
   #        int.TeX intensity formatted for LaTeX (factor of strings)
   # attr:  This function sets the following attributes:
   #        ApplicationName,
   #        ApplicationVersion,
   #        pdfNumber,
   #        chemicalformula,
   #        wavelength
   #
   require(XML)
   doc <- xmlTreeParse(pdffile)
   pdf <- xmlRoot(doc)
   rmchar <- "[^0-9]"
   #
   angles <- data.frame(NULL)
   for (i in 1:length(pdf[["graphs"]][["stick_series"]])) {
   angles <- rbind(angles, data.frame(#
      thth      = as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["theta"]])),
      d         = as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["da"]])),
      h         = as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["h"]])),
      k         = as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["k"]])),
      l         = as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["l"]])),
      hkl       = paste(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["h"]]),
                     xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["k"]]),
                     xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["l"]]), sep = ""),
      hkl.TeX   = paste("\\mbox{$", ifelse(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["h"]])) < 0,
                     paste("\\bar{", abs(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["h"]]))),
                        "}", sep = ""),
                     xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["h"]])),
                  "\\,", ifelse(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["k"]])) < 0,
                     paste("\\bar{", abs(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["k"]]))),
                        "}", sep = ""),
                     xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["k"]])),
                  "\\,", ifelse(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["l"]])) < 0,
                     paste("\\bar{", abs(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["l"]]))),
                        "}", sep = ""),
                     xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["l"]])),
                  "$}", sep = "", collapse = ""),
      intensity = as.numeric(gsub(rmchar, "", xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["intensity"]]))),
      int.TeX   = paste("{", xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["intensity"]]), "}", sep = "")
      ))
   }
   #
   attr(angles, "ApplicationName") <- xmlAttrs(pdf)[[1]]
   attr(angles, "ApplicationVersion") <- xmlAttrs(pdf)[[2]]
   attr(angles, "pdfNumber") <- xmlValue(pdf[["pdf_data"]][["pdf_number"]])
   attr(angles, "chemicalformula") <- gsub("[ ]", "", xmlValue(pdf[["pdf_data"]][["chemical_formula"]]))
   attr(angles, "wavelength") <- as.numeric(xmlValue(pdf[["graphs"]][["wave_length"]]))
   # Caution: Do not subset. Subsetting causes all attributes to be lost.
   return(angles)
}




##################################################
################### uxd2mtx ######################
##################################################
uxd2mtx <- function(uxdfile) {
   # Function for reading UXD files 
   # Assumptions: data in two columns
   # Args: uxdfile (filename with extension)
   # Return value: matrix with two columns
   
   cchar <- "[;_]" #regexpr matching the comment characters used in Bruker's UXD
   cdata <- "[0-9]" #regexpr matching one character of any digit
   
   # A new file (datafile) containing only data will be created,
   # extension ".data" appended to uxdfile
   #datafile <- paste(uxdfile,".data",sep="")
   
   ufile <- file(uxdfile, "r")
   f <- readLines(ufile, n=-1) #read _all_ lines from UXD file
   close(ufile)
   
   # This way we identify data rows by looking for numeric characters.
   #wh <- regexpr("[0-9]", f)
   # This way we identify header rows
   # We assume that all other rows are data
   wh <- regexpr(cchar, f)
   
   mh <- wh[1:length(wh)] # this gives you the corresponding index vector
   # the value of each element corresponds to the position of the regexp match.
   # value = 1 means the first character of the row is cchar (row is header)
   # value =-1 means no cchar occur on the row (row is data)
   
   i <- seq(1, length(mh) - 1, 1)
   j <- seq(2, length(mh), 1)
   
   starts <- which(mh[i] == 1 & mh[j] != 1) + 1
   ends   <- length(mh)
   f <- f[starts:ends]
   
   zz <- textConnection(f, "r")
   ff <- matrix(scan(zz, what = numeric()), ncol=2, byrow=T)
   close(zz)
      
   #zz <- file(datafile, "w") #open connection to datafile
   #write.table(ff, file=datafile, row.names=F, sep=",")
   #close(zz)

   # Return matrix
   ff
}





##################################################
#################### uxd2df ######################
##################################################
uxd2df <- function(uxdfile) {
   # Function for reading UXD files  # Assumptions: data in two columns
   # Args: uxdfile (filename with extension)
   # Returns: dataframe with three columns
   
   cchar <- "[;_]" #regexpr matching the comment characters used in Bruker's UXD
   cdata <- "[0-9]" #regexpr matching one character of any digit
   
   # A new file (datafile) containing only data will be created,
   # extension ".data" appended to uxdfile
   #datafile <- paste(uxdfile,".data",sep="")
   
   ufile <- file(uxdfile, "r")
   f <- readLines(ufile, n=-1) #read _all_ lines from UXD file
   close(ufile)
   
   # This way we identify data rows by looking for numeric characters.
   #wh <- regexpr("[0-9]", f)
   # This way we identify header rows
   # We assume that all other rows are data
   wh <- regexpr(cchar, f)
   
   mh <- wh[1:length(wh)] # this gives you the corresponding index vector
   # the value of each element corresponds to the position of the regexp match.
   # value = 1 means the first character of the row is cchar (row is header)
   # value =-1 means no cchar occur on the row (row is data)
   
   i <- seq(1, length(mh) - 1, 1)
   j <- seq(2, length(mh), 1)
   
   starts <- which(mh[i] == 1 & mh[j] != 1) + 1
   ends   <- length(mh)
   f <- f[starts:ends]
   
   zz <- textConnection(f, "r")
   ff <- data.frame(uxdfile, matrix(scan(zz,
         what = numeric()), ncol=2, byrow=T))
   names(ff) <- c("sampleid", "angle", "intensity")
   close(zz)
      
   #zz <- file(datafile, "w") #open connection to datafile
   #write.table(ff, file=datafile, row.names=F, sep=",")
   #close(zz)

   # Return dataframe
   ff
}




##################################################
################### muxd2df ######################
##################################################
muxd2df <- function(uxdfile, range.descriptor) {
   # Function that reads an UXD file which contains several ranges
   # (created in a programmed run, for example)
   # Arguments
   # :: uxdfile (filename with extension)
   # :: range.descriptor (an array with as many elements as
   #    there are ranges in the uxdfile)
   # Returns: dataframe with 3 columns
   
   cchar <- "[;_]" #regexpr matching the comment characters used in Bruker's UXD
   cdata <- "[0-9]" #regexpr matching one character of any digit
   # Create filenames for the output # no longer used, return dataframe instead
   #datafile <- paste(uxdfile,"-",range.descriptor,".data",sep="")
   
   # Read the input multirange file
   ufile <- file(uxdfile, "r")
   f <- readLines(ufile, n=-1) #read _all_ lines from UXD file
   close(ufile)
   
   # This way we identify data rows by looking for numeric characters.
   #wh <- regexpr("[0-9]", f)
   # This way we identify header rows
   # Later we will assume that all other rows are data
   wh <- regexpr(cchar, f)
   
   mh <- wh[1:length(wh)] # this gives you the corresponding index vector
   # the value of each element corresponds to the position of the regexp match.
   # value = 1 means the first character of the row is cchar (row is header)
   # value =-1 means no cchar occur on the row (row is data)
   
   #length(mh[mh == -1]) #total number of datarows in uxdfile
   #mh[mh > 1 | mh < 0] <- 0 #set all header-rows to zero (just to make things easier)
   
   i <- seq(1, length(mh) - 1, 1)
   j <- seq(2, length(mh), 1)
   starts <- which(mh[i] == 1 & mh[j] != 1) + 1 #start indices
   ends   <- which(mh[i] != 1 & mh[j] == 1) #end indices, except the last
   ends   <- c(ends, length(mh)) #fixed the last index of ends   
   
   ff <- data.frame(NULL)
   for (s in 1:length(range.descriptor)) {
      zz <- textConnection(f[starts[s]:ends[s]], "r")
      ff <- rbind(ff, data.frame(range.descriptor[s],
            matrix(scan(zz, what = numeric()), ncol=2, byrow=T)))
      close(zz)
   }
   names(ff) <- c("sampleid", "angle", "intensity")
   
   # Return dataframe
   ff
}




##################################################
################### muxd2mtx #####################
##################################################
muxd2mtx <- function(uxdfile, range) {
   # Function that reads an UXD file which contains several ranges
   # (created in a programmed run, for example)
   # Arguments
   # :: uxdfile (filename with extension)
   # :: range (an integer, the number of ranges in the uxdfile)
   # Returns: matrix with 2 columns
   
   cchar <- "[;_]" #regexpr matching the comment characters used in Bruker's UXD
   cdata <- "[0-9]" #regexpr matching one character of any digit
   # Create filenames for the output # no longer used, return dataframe instead
   #datafile <- paste(uxdfile,"-",range.descriptor,".data",sep="")
   
   # Read the input multirange file
   ufile <- file(uxdfile, "r")
   f <- readLines(ufile, n=-1) #read _all_ lines from UXD file
   close(ufile)
   
   # We identify header rows. We will assume that all other rows are data
   wh <- regexpr(cchar, f)
   # length(wh) equals length(f)
   # wh has either 1 or -1 for each element
   # value = 1 means the first character of that row is cchar (row is header)
   # value =-1 means absence of cchar on that row (row is data)
   
   # Since wh contains some attributes (given by regexpr function), we strip out everything
   # but the index vector and assign it to the new vector mh
   mh <- wh[1:length(wh)] # this gives you the corresponding index vector
   
   #length(mh[mh == -1]) #total number of datarows in uxdfile
   #mh[mh > 1 | mh < 0] <- 0 #set all header-rows to zero (just to make things easier)
   
   # Set counters i and j used in assignment below
   i <- seq(1, length(mh) - 1, 1)
   j <- seq(2, length(mh), 1)
   
   starts <- which(mh[i] == 1 & mh[j] != 1) + 1 #start indices
   ends   <- which(mh[i] != 1 & mh[j] == 1) #end indices, except the last
   ends   <- c(ends, length(mh)) #fixes the last index of ends
   # note that length of starts (or ends) gives the number of ranges
   
   ff <- matrix(NA, length(f), 2)
   for (s in 1:range) {
      zz <- textConnection(f[starts[s]:ends[s]], "r")
      ff <- rbind(ff, matrix(scan(zz, what = numeric()), ncol=2, byrow=T))
      close(zz)
   }
   
   # Clean up matrix: remove extra rows
   ff[apply(ff,1,function(x)any(!is.na(x))), ]

   # Return matrix
   ff
}



##################################################
################### muxd2ls ######################
##################################################
muxd2ls <- function(uxdfile) {
   # Function that reads an UXD file which contains several ranges
   # (created in a programmed run, for example)
   # Arguments
   # :: uxdfile (filename with extension)
   # Returns: List of matrices, as many as there were ranges
   # Requires: ??
   # See extensive comments in muxd2mtx()
   
   cchar <- "[;_]" #comment characters used in Bruker's UXD
   cdata <- "[0-9]"  
   
   ufile <- file(uxdfile, "r")
   f <- readLines(ufile, n=-1) #read _all_ lines from UXD file
   close(ufile)
   
   wh <- regexpr(cchar, f)
   mh <- wh[1:length(wh)]
   
   i <- seq(1, length(mh) - 1, 1)
   j <- seq(2, length(mh), 1)
   starts <- which(mh[i] == 1 & mh[j] != 1) + 1
   ends   <- which(mh[i] != 1 & mh[j] == 1)
   ends   <- c(ends, length(mh))
   
   ff <- list()
   for (s in 1:length(starts)) {
      zz <- textConnection(f[starts[s]:ends[s]], "r")
      ms <- matrix(scan(zz, what = numeric()), ncol = 2, byrow = T)
      close(zz)
      ff[[s]] <- ms
   }
   # Return list of matrices
   return(ff)
}








# -------- ##################### -------- #
# -------- #### REPAIR SHOP #### -------- #
# -------- ##################### -------- #


##################################################
################ pearson.beta ####################
##################################################
pearson.beta <- function(m, w) {
   # Function for calculating integral breadth, \beta, of
   # an X-ray reflection using the half-width of FWHM (w)
   # and m parameter from a Pearson fit.
   ## Args: m - shape parameter from Pearson fit (numeric)
   ##       w - half the FWHM for the reflection (numeric)
   ## **** NOTE: gamma(x) out of range for x > 171
   ## Values: The integral breadth of the reflection (numeric)
   # Ref: Birkholz2006, page 96 (Eq. 3.5)
   beta <- ((pi * 2^(2 * (1 - m)) * gamma(2 * m - 1)) /
      ((2^(1 / m) - 1) * (gamma(m))^2)) * w
}



##################################################
################## strip.ka2 #####################
##################################################
strip.ka <- function(angle, intensity) {
   # Function that strips (removes) the Cu Ka2 contribution
   # from an experimental diffractogram (2Th, I(2Th))
   # using the Rachinger algorithm.
   ## Args: 2theta vector
   ##       intensity vector
   ## Values: intensity vector
   # Ref: Birkholz2006, page 99 (Eq. 3.12), and page 31
   Cu.Kai  <- 0.154059 # nanometre
   Cu.Kaii <- 0.154441 # nanometre
   alpha.int.ratio <- 0.5
   int.stripped <- intensity - alpha.int.ratio*intensity
   return(int.stripped)
}




print.xtable.booktabs <- function(x){
   # Make xtable print booktabs tables
   # ARGS: x (matrix)
   require(xtable)
   print(xtable(x), floating=F, hline.after=NULL,
    add.to.row=list(pos=list(-1,0, nrow(x)),
    command=c("\\toprule ", "\\midrule ", "\\bottomrule ")))
}


split.muxd <- function(uxdfile, range.descriptor) {
   # Fix this some other day!!
   # Function that reads an UXD file which contains several ranges
   # (created in a programmed run, for example) and splits it into
   # its ranges and writes to that number of files
   # Arguments
   # :: uxdfile (filename with extension)
   # :: range.descriptor (an array with as many elements as
   #    there are ranges in the uxdfile)
   # Returns: nothing. Writes files to HDD.
   
   cchar <- "[;_]" #regexpr matching the comment characters used in Bruker's UXD
   cdata <- "[0-9]" #regexpr matching one character of any digit
   # Create filenames for the output # no longer used, return dataframe instead
   datafile <- paste(uxdfile,"-",range.descriptor,".data",sep="")
   
   # Read the input multirange file
   f <- readLines(uxdfile, n=-1)
   
   
   # This way we identify data rows by looking for numeric characters.
   #wh <- regexpr("[0-9]", f)
   # This way we identify header rows
   # Later we will assume that all other rows are data
   wh <- regexpr(cchar, f)
   
   mh <- wh[1:length(wh)] # this gives you the corresponding index vector
   # the value of each element corresponds to the position of the regexp match.
   # value = 1 means the first character of the row is cchar (row is header)
   # value =-1 means no cchar occur on the row (row is data)
   
   #length(mh[mh == -1]) #total number of datarows in uxdfile
   #mh[mh > 1 | mh < 0] <- 0 #set all header-rows to zero (just to make things easier)
   
   i <- seq(1, length(mh) - 1, 1)
   j <- seq(2, length(mh), 1)
   starts <- which(mh[i] == 1 & mh[j] != 1) + 1 #start indices
   ends   <- which(mh[i] != 1 & mh[j] == 1) #end indices, except the last
   ends   <- c(ends, length(mh)) #fixed the last index of ends
   
   #ff <- data.frame(NULL)
   for (s in 1:length(range.descriptor)) {
      matrix(scan(file=textConnection(f[starts[s]:ends[s]]),
            what = numeric()), ncol=2, byrow=T)
   }
   names(ff) <- c("sampleid", "angle", "intensity")
   
   # Return dataframe
   ff
}
