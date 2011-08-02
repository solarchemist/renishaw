##################################################
#################### SHE2AVS #####################
##################################################
SHE2AVS <- function(she) {
   # Converts from SHE scale to absolute vacuum (AVS) scale
   avs <- -(4.5 + she)
   return(avs)
}