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
#    - EliminateKa2
#    - print.xtable.booktabs
#    - split.muxd
#    - strip.ka2
#    - pearson.beta

## All XRD reading functions need to be re-written to identify the data instead of identifying metadata
## Take inspiration from the functions in CHI.R




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
   #       wavelength      - X-ray wavelength used (default 1.54056 A, Cu Ka)
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
   # For example the PDF files produced by the PDF database at Angstrom's X-ray lab
   # NOTE: sometimes intensity values are specified as less than some value.
   #       In those cases, this function simply strips the less-than character.
   #       (Perhaps not true, see the int.Tex column)
   # ARGS: pdffile (complete path and filename to PDF file)
   # VALUE: dataframe with 9 columns:
   #        thth angles (numeric),
   #        d (numeric),
   #        h index (numeric),
   #        k index (numeric),
   #        l index (numeric),
   #        hkl indices (string),
   #        hkl.TeX indices formatted for LaTeX (string),
   #        intensity (numeric),
   #        int.TeX intensity formatted for LaTeX (string)
   # attr:  This function sets the following attributes:
   #        ApplicationName,
   #        ApplicationVersion,
   #        pdfNumber,
   #        chemicalformula,
   #        empiricalformula,
   #        wavelength
   #
   require(XML)
   doc <- xmlTreeParse(pdffile)
   pdf <- xmlRoot(doc)
   rmchar <- "[^0-9]"
   #
   angles <- data.frame(NULL)
   for (i in 1:length(pdf[["graphs"]][["stick_series"]])) {
   angles <- rbind(angles, data.frame(stringsAsFactors = FALSE,#
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
   attr(angles, "empiricalformula") <- gsub("[ ]", "", xmlValue(pdf[["pdf_data"]][["empirical_formula"]]))
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
   
   rgxp.sampleid <- "[^/]*(?=\\.\\w*)" ## THIS REQUIRES perl=TRUE
   # Regular expression that extracts the filename out of a full path.
   # Matches and extracts everything from the last forward slash (assuming Unix slashes)
   # up until a dot folllowed by an arbitrary number of alphanumeric characters.
   sampleidmtch <- regexpr(rgxp.sampleid, uxdfile, perl=TRUE)
   # Check that there was a match
   if (sampleidmtch < 0) {
      # -1 means no match
      sampleid <- uxdfile
      # If match was unsuccessful we use the argument as passed to this function as sampleid
   }
   sampleid <- substr(uxdfile, sampleidmtch, (sampleidmtch + attr(sampleidmtch, "match.length") - 1))
   
   zz <- textConnection(f, "r")
   ff <- data.frame(sampleid, matrix(scan(zz,
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
################ EliminateKa2 ####################
##################################################
EliminateKa2 <- function(thth, intensity, lever = 2) {
   # Args: 2theta vector, numeric
   #       intensity vector, numeric
   #       parameter pairs, numeric between 1 and 6, default 2
   # "High-quality a2-free patterns can be obtained in most cases
   #  using five (5) or seven (7) pairs of parameters."
   # "When the step-width is less than 0.01 degrees and the
   #  2theta angle is high, a large number of parameter pairs
   #  should be used to get accurate results." {Dong1999a}
   ### 1 2 3 4  5  6 - lever
   ### 3 5 7 9 15 25 - corresponding parameter pairs
   #
   ##### THIS FUNCTION USES THESE OTHER FUNCTIONS #####
   # REQUIRES: common.R :: as.radians() - converts degrees to radians
   # REQUIRES: stats::approx - linear interpolation
   ####################################################
   
   # For test-purposes, use the following data
   # /home/taha/chepec/laboratory/XRD/0103-instrumentbroadening/100917Th2ThLong-counts.UXD
   
   ##### STILL UNDER CONSTRUCTION ####
   ##### STILL UNDER CONSTRUCTION ####
   
   startdatapoint <- 4
   

   #              Ka1a      Ka1b       Ka2a      Ka2b
   CuKa.data <- c(1.534753, 1.540596, 1.541058,  1.544410, 1.544721,
                  3.69,     0.44,     0.60,      0.52,     0.62,
                  1.60,    57.07,     7.64,     25.38,     8.31)
   CuKa <- data.frame(matrix(CuKa.data, ncol=3, byrow=F))
   names(CuKa) <- c("lambda", "w", "E")
   row.names(CuKa) <- c("Satellites", "Ka1a", "Ka1b", "Ka2a", "Ka2b")

   
   # The following lever arm weights are from Dong1999a
   weights <- list()
   # Three-bar weights
   weights[[1]] <- c(0.005134296, 0.491686047, 0.003179657)
   # Five-bar weights   
   weights[[2]] <- c(0.002614410, 0.011928014, 0.480406967, 
                     0.002121807, 0.002928802)
   # Seven-bar weights
   weights[[3]] <- c(0.001580069, 0.003463773, 0.015533472, 0.422601977, 
                     0.053632977, 0.001572467, 0.001615265)
   # Nine-bar weights
   weights[[4]] <- c(0.001138001, 0.00195272,  0.004324464, 
                     0.019246541, 0.394175823, 0.079159001,
                    -0.003591547, 0.002505604, 0.001089392)
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
   
   # The following lever arm lengths are from Dong1999a
   lengths <- list()
   # Three-bar lengths
   lengths[[1]] <- c(0.998506815, 0.997503913, 0.996699460)
   # Five-bar lengths   
   lengths[[2]] <- c(0.998471576, 0.997935524, 0.997503530, 
                     0.997163494, 0.996606519)
   # Seven-bar lengths
   lengths[[3]] <- c(0.998563433, 0.998204025, 0.997825027, 0.997522195, 
                     0.997297615, 0.996844235, 0.996516288)
   # Nine-bar lengths
   lengths[[4]] <- c(0.998609749, 0.998334027, 0.998054914,
                     0.99776062,  0.997527844, 0.997327154, 
                     0.997028978, 0.996734639, 0.99646335)
   # 15-bar lengths
   lengths[[5]] <- c(0.998671599, 0.99850911,  0.998346447, 0.998183442, 0.998019704,
                     0.997854063, 0.997680649, 0.997533314, 0.997377391, 0.997266106,
                     0.997060614, 0.996888005, 0.996741151, 0.996583672, 0.996418168)
   # 25-bar lengths
   lengths[[6]] <- c(0.998706192, 0.998608958, 0.998511721, 0.998414475, 0.998317209,
                     0.998219906, 0.998122538, 0.998025057, 0.997927367, 0.997829244,
                     0.997730044, 0.997626987, 0.997535705, 0.997458223, 0.997346989,
                     0.997277763, 0.997161452, 0.997057942, 0.996982688, 0.99686108,
                     0.996769728, 0.996675255, 0.996578407, 0.996480641, 0.99638324)
   
   ### --- Arguments check
   # Check that "lever" argument is within bounds
   if (lever > length(weights) || lever < 1) {
      # if not, fall back to the default value
      lever <- 2 # corresponds to 5 parameter pairs
   }
   
   ### --- Arguments check
   # Check that vectors are of the same length
   if (!(length(thth) == length(intensity))) {
      # If not the same length, abort with error message
      stop("Arguments thth and intensity have different lengths!")
   }
   
   ### --- Arguments check
   if (any(thth <= 0)) {
      stop("thth vector contains values less-than or equal to zero")
   }
   
   ### --- Arguments check
   if (any(intensity < 0)) {
      stop("intensity vector contains values less than zero")
   }
   
   # THIS IS NECESSARY, but overlooked for the moment
   #int.p.start <- (1 / (1 + (CuKa["Ka2a", "E"] + CuKa["Ka2b", "E"]) /
      #(CuKa["Ka1a", "E"] + CuKa["Ka1b", "E"]))) * intensity[1:startdatapoint-1]
   
   # Redefine everything
   #thth <- thth[startdatapoint:length(thth)]
   #intensity <- intensity[startdatapoint:length(intensity)]
   
   
   # Convert from 2theta to theta scale for use in first step of calculations
   theta <- thth / 2
   sintheta <- sin(as.radians(theta))
   
   # Corresponds to equation 10 in Dong1999a
   # This is based on the assumption that we are supposed to get delta-thth values
   delta.thth.a <- matrix(0, length(theta), length(weights[[lever]]))
   for (j in 1:length(lengths[[lever]])) {
      delta.thth.a[,j] <- 2 * asin(as.radians((lengths[[lever]][j] * sintheta)))
   }
   
   # Add the calculated deltas to the recorded thth values
   # Corresponds to equation 10 in Dong1999a
   thth.a <- matrix(NA, dim(delta.thth.a)[1], dim(delta.thth.a)[2] + 1)
   # Flip the delta.thth.a matrix vertically (just for convenience)
   delta.thth.a <- delta.thth.a[,dim(delta.thth.a)[2]:1]
   for (j in 2:dim(thth.a)[2]) {
      thth.a[,1] <- thth
      thth.a[,j] <- thth + delta.thth.a[,j - 1]
   }
   
   # Intensities with interpolated intensities at the calculated 2theta values
   int.interp <- matrix(NA, dim(thth.a)[1], dim(thth.a)[2])
   for (j in 2:dim(int.interp)[2]) {
      int.interp[,1] <- intensity
      int.interp[,j] <- approx(thth, intensity, xout = thth.a[,j])$y
   }
   # So far, we have just replaced the old thth-scale with a new one,
   # and calculated the intensitites by linearly interpolating from the old intensities.
   
   # Intensities times lever weights, P(j) (this is what you will substract from int.interp)
   int.p <- matrix(NA, length(theta), length(weights[[lever]]))
   for (j in 1:length(lengths[[lever]])) {
      int.p[,j] <- weights[[lever]][j] * int.interp[,j + 1]
   }
   
   # Calculate intensities with Ka2 contribution stripped
   #int.stripped <- c(int.interp[,-1]) - rowSums(int.p)
   int.stripped <- thth - rowSums(int.p)
   #corr.tmp.df <- data.frame(thth = c(thth.a[,-1]), int.corr = int.stripped)
   corr.tmp.df <- data.frame(thth = thth, int.corr = int.stripped)
   corr.df <- corr.tmp.df[order(corr.tmp.df$thth), ]
   #row.names(corr.df) <- seq(1, length(corr.df$thth))
   
   # Make a dataframe of thth.a and int.a and order it by thth
   #df.a <- data.frame(thth = c(thth.a), intensity = c(int.a))
   #df.as <- df.a[order(df.a$thth), ]
   #row.names(df.as) <- seq(1,length(df.as$thth)) # fixes row names order
   # df.as is exactly as the original data, just with the correct number of
   # interpolated thth and intensity values included.
   
   
   return(corr.df)
   
   
   
   
   # Perhaps thth.a and int.p.terms are the new x and y
   
   
   # Intensities in summation term (this is the Ka2 correction)
   #int.p <- rowSums(int.p.terms)
   
   # Collapse the matrix thth.a into a single column
   #thth.ai <- matrix(NA, dim(thth.a)[1] * dim(thth.a)[2], 2)
   #thth.ai[,1] <- sort(c(thth.a))
   #thth.ai[,2] <- approx(thth, intensity, xout = thth.ai[,1])$y
   
   
   # Built-in functions in R to interpolate or extrapolate
   # stats::approx         Linear interpolation
   # Hmisc::approxExtrap   Linear extrapolation
   # go with linear functions for now, see how that works out
   


   
   
   
   
   
   
   
   
   # This is NOT necessary
   # This vector helps convert from lever to actual number of parameter pairs
   #parpairs <- numeric()
   #for (j in 1:length(weights)) {
   #   parpairs[j] <- length(weights[[j]])
   #}
   
   
   ##### STILL UNDER CONSTRUCTION ####
}




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
