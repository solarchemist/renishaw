# renishaw

An R package that imports ASCII-formatted Raman spectra files from Renishaw Wire,
saves them as an R tibble along with a user-defined `sampleid`.

The package also contains functions to baseline and peak fit the spectra
using the `diffractometry` package, and a wrapper function to call them
and save the fitting results to disk.


## Notes

+ The `diffractometry` package was archived from CRAN on 2019-02-05.
  Latest available version was 0.1-10.
  I store a copy of the package at https://public.solarchemist.se/diffractometry
+ https://rdrr.io/cran/diffractometry
