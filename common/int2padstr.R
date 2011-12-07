##################################################
################## int2padstr ####################
##################################################
int2padstr <- function (ii, pchr, w) {
   ## Description:
   ##   Converts an integer or a vector of integers to
   ##   a string padded with characters. 
   ## Usage:
   ##   int2padstr(ii, pchr, w)
   ## Arguments:
   ##     ii: integer or vector of integers
   ##   pchr: a padding character (e.g., "0")
   ##      w: width of the return string (an integer)
   ##         Make sure to set the width longer than
   ##         or equal to the length of the biggest integer.
   ##         For example, if the integers (ii) are
   ##         in the range 1 - 100, set w to at least 3.
   ## Value:
   ##   A character string or a vector of character strings
   gsub(" ", pchr, formatC(ii, format="s", mode="character", width = w))
}
