#' Import Renishaw Raman spectra files (ASCII)
#'
#' Reads Raman (ASCII) datafile and returns a dataframe with the original data,
#' as well as interpolated wavenumber and counts values
#' (interpolated data is evenly spaced along x-axis, required for later peak fitting
#' with diffractometry package).
#'
#' @param datafilename text string with full path to experimental file
#'
#' @return Dataframe with the following columns (and no extra attributes):
#'    $ sampleid        : char
#'    $ wavenum         : numeric
#'    $ counts          : numeric
#'    $ wnum.interp     : interpolated wavenumber vector (equidistant)
#'    $ cts.interp      : interpolated counts vector
#' @export
Raman2df <- function(datafilename) {

   chifile <- base::readLines(datafilename, n = -1)
   sampleid <- common::ProvideSampleId(datafilename)

   ff <- data.frame(NULL)
   zz <- base::textConnection(chifile, "r")
   ff <-
      base::rbind(
         ff,
         data.frame(
            sampleid,
            base::matrix(
               base::scan(
                  zz,
                  what = numeric(),
                  sep = "\t"),
               ncol = 2,
               byrow = TRUE),
            stringsAsFactors = FALSE))
   close(zz)
   names(ff) <- c("sampleid", "wavenum", "counts")

   # Sort dataframe by increasing wavenumbers
   ff <- ff[order(ff$wavenum), ]
   # ... sort the rownames as well
   row.names(ff) <- seq(1, dim(ff)[1])

   # Create evenly spaced datapoints by interpolating
   ff <-
      base::cbind(
         ff,
         wnum.interp =
            stats::approx(
               x = ff$wavenum,
               y = ff$counts,
               method = "linear",
               n = length(ff$wavenum))$x,
         cts.interp =
            stats::approx(
               x = ff$wavenum,
               y = ff$counts,
               method = "linear",
               n = length(ff$wavenum))$y)

   return(ff)
}
