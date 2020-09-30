#' Reset xtable attributes for printing as LaTeX
#'
#' Warning: this function needs the following variables to exist in the global namespace:
#' - tp
#' - obs
#' (this should be changed!). But for now, watch out.
#'
#' @param table.object  xtable object
#'
#' @return xtable object
#' @export
ResetTabAttributes <- function(table.object) {
   xtable::caption(table.object) <-
      paste0(
        obs$substrateid[tp], "-", common::int2padstr(tp, "0", 3),
        " (area ", obs$sampling.area[tp], ")",
        ", ", obs$objective.mag[tp],
        ", \\SI{", obs$exposure.time[tp], "}{\\second}",
        ", ", obs$accumulations[tp], " acc.",
        ", \\SI{", obs$laser.power[tp], "}{\\percent}")
   xtable::label(table.object) <-
     paste("tab:peak-table-", common::int2padstr(tp, "0", 3), sep = "")
   names(table.object) <-
      c("{Peak}",
        "{Kernel}",
        "{Wavenumber/\\si{\\per\\cm}}",
        "{Height/\\si{\\counts}}",
        "{Area/\\si{\\counts\\per\\cm}}",
        "{FWHM/\\si{\\per\\cm}}",
        "{$m$}",
        "{Accept}")
   xtable::digits(table.object) <-
      c(0, #row.names
        2, #peak
        1, #kernel
        1, #wavenumber
        1, #height
        1, #area
        1, #FWHM
        1, #m
        0) #accept
   xtable::display(table.object) <-
      c("s", #row.names
        "d", #peak
        "d", #kernel
        "f", #wavenumber
        "f", #height
        "f", #area
        "f", #FWHM
        "f", #m
        "s") #accept
   xtable::align(table.object) <-
      c("l", #row.names
        "S[table-format=2.0]", #peak
        "S[table-format=2.0]", #kernel
        "S[table-format=4.1]", #wavenumber
        "S[table-format=4.1]", #height
        "S[table-format=2.2]", #area
        "S[table-format=2.2]", #FWHM
        "S[table-format=3.1]", #m
        "c") #accept
   #
   return(table.object)
}
