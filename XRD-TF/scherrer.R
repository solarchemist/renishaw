source("/home/taha/chepec/chetex/common/R/common/as.radians.R")

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
   ## REQUIRES: as.radians()
   D <- (shapeconstant * wavelength) / (as.radians(integralbreadth) * cos(as.radians(thth)))
   # cos() - angles must be in radians, not degrees!
   return(D) #units of angstrom
}
