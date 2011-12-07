# Function that rounds UP to the nearest interval specified by "nearest"
# http://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x

roundup <- function(x, nearest=1000) {
   ceiling(max(x+10^-9)/nearest + 1/nearest)*nearest
}
