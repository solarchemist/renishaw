pdf2df <- function(pdffile) {
   # Function for extracting information from ICDD PDF XML-files
   # For example the PDF files produced by the PDF database at Angstrom's X-ray lab
   # NOTE: sometimes intensity values are specified as less than some value.
   #       In those cases, the lt sign will be preserved in the column int.Tex.
   #       The intensity column, on the other hand, is numeric and so strips off the lt sign.
   # ARGS: pdffile (complete path and filename to PDF file)
   # VALUE: dataframe with 9 columns:
   #        thth angles (numeric),
   #        d (numeric),
   #        h index (numeric),
   #        k index (numeric),
   #        l index (numeric),
   #        hkl indices (string),
   #        hkl.TeX indices formatted for LaTeX (string),
   #        intensity (numeric),
   #        int.TeX intensity formatted for LaTeX (string),
   #        pdfNumber (string)
   # attr:  This function sets the following attributes:
   #        ApplicationName,
   #        ApplicationVersion,
   #        chemicalformula,
   #        empiricalformula,
   #        wavelength
   #
   require(XML)
   doc <- xmlTreeParse(pdffile)
   pdf <- xmlRoot(doc)
   rmchar <- "[^0-9]"
   #
   angles <- data.frame(NULL)
   for (i in 1:length(pdf[["graphs"]][["stick_series"]])) {
   angles <- rbind(angles, data.frame(stringsAsFactors = FALSE,#
      thth      = as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["theta"]])),
      d         = as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["da"]])),
      h         = as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["h"]])),
      k         = as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["k"]])),
      l         = as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["l"]])),
      hkl       = paste(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["h"]]),
                     xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["k"]]),
                     xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["l"]]), sep = ""),
      hkl.TeX   = paste("\\mbox{$", ifelse(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["h"]])) < 0,
                     paste("\\bar{", abs(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["h"]]))),
                        "}", sep = ""),
                     xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["h"]])),
                  "\\,", ifelse(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["k"]])) < 0,
                     paste("\\bar{", abs(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["k"]]))),
                        "}", sep = ""),
                     xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["k"]])),
                  "\\,", ifelse(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["l"]])) < 0,
                     paste("\\bar{", abs(as.numeric(xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["l"]]))),
                        "}", sep = ""),
                     xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["l"]])),
                  "$}", sep = "", collapse = ""),
      intensity = as.numeric(gsub(rmchar, "", xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["intensity"]]))),
      int.TeX   = paste("{", xmlValue(pdf[["graphs"]][["stick_series"]][[i]][["intensity"]]), "}", sep = ""),
      pdfNumber = xmlValue(pdf[["pdf_data"]][["pdf_number"]]),
      formula   = gsub("[ ]", "", xmlValue(pdf[["pdf_data"]][["empirical_formula"]])) 
      ))
   }
   #
   attr(angles, "ApplicationName") <- xmlAttrs(pdf)[[1]]
   attr(angles, "ApplicationVersion") <- xmlAttrs(pdf)[[2]]
   #attr(angles, "pdfNumber") <- xmlValue(pdf[["pdf_data"]][["pdf_number"]])
   attr(angles, "chemicalformula") <- gsub("[ ]", "", xmlValue(pdf[["pdf_data"]][["chemical_formula"]]))
   attr(angles, "empiricalformula") <- gsub("[ ]", "", xmlValue(pdf[["pdf_data"]][["empirical_formula"]]))
   attr(angles, "wavelength") <- as.numeric(xmlValue(pdf[["graphs"]][["wave_length"]]))
   # Caution: Do not subset. Subsetting causes all attributes to be lost.
   return(angles)
}
