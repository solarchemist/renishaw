# This function was copied from R's documentation (see ?is.integer).

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
   abs(x - round(x)) < tol
}
