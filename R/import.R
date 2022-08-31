#' Import Renishaw Raman spectra files (ASCII)
#'
#' Reads Raman (ASCII) datafile and returns a dataframe with the original data,
#' as well as interpolated wavenumber and counts values
#' (interpolated data is evenly spaced along x-axis, required for later peak fitting
#' with diffractometry package).
#'
#' @param datafilename text string with full path to experimental file
#' @param sampleid sampleid (optional, automatically created if not provided)
#'
#' @return Tibble with the following columns (and possibly some attributes):
#'    $ sampleid        : char
#'    $ wavenum         : numeric
#'    $ counts          : numeric
#'    $ wnum.interp     : interpolated wavenumber vector (equidistant)
#'    $ cts.interp      : interpolated counts vector
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
Raman2df <- function(datafilename, sampleid = "") {

   if (sampleid == "") {
      # assume the user did not set sampleid
      this.sampleid <- common::ProvideSampleId(datafilename)
   } else {
      # assume the user explicitly specified sampleid, use it
      this.sampleid <- sampleid
   }

   ff <-
      readr::read_tsv(
         file = datafilename,
         # datafiles expected to contain column labels: "#Wave #Intensity"
         # force read_tsv to discard empty third column that it detects
         # because the header row uses two tabs as separator
         col_select = 1:2,
         col_names = c("wavenum", "counts"),
         # since we set col_names explicitly, we need to skip the first row
         # in the datafile which contains column labels
         skip = 1,
         col_types = "dd") %>%
      tibble::add_column(
         sampleid = this.sampleid,
         .before = "wavenum") %>%
      # sort dataframe by increasing wavenumbers
      dplyr::arrange(.data$wavenum)

   # Create evenly spaced datapoints by interpolating
   ff <- ff %>%
      dplyr::bind_cols(
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
