##################################################
################# as.degrees #####################
##################################################
as.degrees <- function(radians) {
   # Converts from radians to degrees
   degrees <- radians * (180 / pi)
   return(degrees)
}