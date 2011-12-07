##################################################
#################### AVS2SHE #####################
##################################################
AVS2SHE <- function(avs) {
   # Converts from absolute vacuum scale (AVS) to SHE scale
   she <- -(4.5 + avs)
   return(she)
}
