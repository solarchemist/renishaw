##################################################
################ tifftags2df #####################
##################################################
tifftags2df <- function(tiffimage) {
   ## Description:
   ##   Extracts all tags from a TIFF image file 
   ##   using the tiffinfo tool and stores
   ##   a selection of the tags in a dataframe.
   ## Usage:
   ##   tifftags2df(fulltiffpath)
   ## Arguments:
   ##    tiffimage: character string, the full filename
   ##               (with path) to one TIFF file.
   ## Value:
   ##   A dataframe with three columns: 
   ##   sampleid, parameter, value.
   #
   substrate.id <- strsplit(basename(tiffimage), "\\.")[[1]][1]
   tifftags <- system(paste("tiffinfo", tiffimage, sep = " "), 
      intern = TRUE, ignore.stderr = TRUE)
   
   # Clean certain special characters from the tiff tags strings
   # If these strings are left untreated, they cause all sorts of weird errors later
   tifftags.clean <- 
      sub("[ ]+$", "",  #trim trailing spaces, if any
      sub("\\r$", "",   #remove trailing \r
      gsub("\\xb9", "", #remove special character
      gsub("\\xb8", "", #remove special character
      gsub("\\xb7", "", #remove special character
      gsub("\\xb6", "", #remove special character
      gsub("\\xb5", "", #remove special character
      gsub("\\xb4", "", #remove special character
      gsub("\\xb3", "", #remove special character
      gsub("\\xb2", "", #remove special character
      gsub("\\xb1", "", #remove special character
      gsub("\\xb0", "", tifftags, 
         useBytes=T), 
         useBytes=T),
         useBytes=T),
         useBytes=T),
         useBytes=T),
         useBytes=T),
         useBytes=T),
         useBytes=T),
         useBytes=T),
         useBytes=T)))
   
   
   
   # These are the tags we are looking for
   tags.df <- data.frame(stringsAsFactors = FALSE, substrate.id, matrix(#
      # Name # REGEXP # splitchar # unit
      c("EHT",                 "AP\\_ACTUALKV",                " = ", "\\kilo\\volt",
        "High current",        "DP\\_HIGH\\_CURRENT",          " = ", "",
        "WD",                  "AP\\_WD",                      " = ", "\\milli\\metre",
        "Magnification",       "AP\\_MAG",                     " = ", "", # could warrant special treatment
        "Brightness",          "AP\\_BRIGHTNESS",              " = ", "\\percent",
        "Contrast",            "AP\\_CONTRAST",                " = ", "\\percent",
        "Signal A",            "DP\\_DETECTOR\\_CHANNEL",      " = ", "",
        "SCM status",          "DP\\_SCM\\_STATUS",            " = ", "",
        "Specimen current",    "AP\\_SCM",                     " = ", "\\femto\\ampere",
       # Beam
       #"Filament curent",     "AP\\_ACTUALCURRENT",           " = ", "\\ampere",
        "Stigmation Y",        "AP\\_STIG\\_Y",                " = ", "\\percent",
        "Stigmation X",        "AP\\_STIG\\_X",                " = ", "\\percent",
        "Aperture align Y",    "AP\\_APERTURE\\_ALIGN\\_Y",    " = ", "\\percent",
        "Aperture align X",    "AP\\_APERTURE\\_ALIGN\\_X",    " = ", "\\percent",
        "Aperture size",       "AP\\_APERTURESIZE",            " = ", "\\micro\\metre",
        "Beam shift Y",        "AP\\_BEAMSHIFT\\_Y",           " = ", "\\percent",
        "Beam shift X",        "AP\\_BEAMSHIFT\\_X",           " = ", "\\percent",
       #"Beam offset Y",       "AP\\_BEAM\\_OFFSET\\_Y",       " = ", "\\nano\\metre",
       #"Beam offset X",       "AP\\_BEAM\\_OFFSET\\_X",       " = ", "\\nano\\metre",
       # Stage
        "Track Z",             "DP\\_TRACK\\_Z",               " = ", "",
        "Stage at Z",          "AP\\_STAGE\\_AT\\_Z",          " = ", "\\milli\\metre",
        "Stage at Y",          "AP\\_STAGE\\_AT\\_Y",          " = ", "\\milli\\metre",
        "Stage at X",          "AP\\_STAGE\\_AT\\_X",          " = ", "\\milli\\metre",
       # Tilt 
        "Stage tilted?",       "DP\\_STAGE\\_TILTED",          " = ", "",
        "Tilt angle",          "AP\\_TILT\\_ANGLE",            " = ", "",
        "Tilt axis",           "AP\\_TILT\\_AXIS",             " = ", "",
       # Image                 
        "Scan speed",          "DP\\_SCANRATE",                " = ", "",
        "Cycle time",          "AP\\_FRAME\\_TIME",            " = ", "\\second",
        "Freeze on",           "DP\\_FREEZE\\_ON",             " = ", "",
        "Dwell time",          "DP\\_DWELL\\_TIME",            " = ", "\\nano\\second",
        "Noise reduction",     "DP\\_NOISE\\_REDUCTION",       " = ", "",
        "Frames to integrate", "AP\\_FRAME\\_INT\\_COUNT",     " = ", "",
        "Frames to average",   "AP\\_FRAME\\_AVERAGE\\_COUNT", " = ", "",
       #"Profile width",       "AP\\_PROFILE\\_W",             " = ", "\\micro\\metre",
        "Pixel size",          "AP\\_PIXEL\\_SIZE",            " = ", "\\nano\\metre",
       # System
        "Gun vacuum",          "AP\\_COLUMN\\_VAC",            " = ", "\\milli\\bar",
        "System vacuum",       "AP\\_SYSTEM\\_VAC",            " = ", "\\milli\\bar",
        "Filament age",        "AP\\_FILAMENT\\_AGE",          " = ", "\\hour",
       # Misc
        "Photo no.",           "AP\\_PHOTO\\_NUMBER",          " = ", "",
        "Date",                "AP\\_DATE",                    " :",  "",
        "Time",                "AP\\_TIME",                    " :",  ""),
      ncol = 4, byrow = T))
   names(tags.df) <- c("sampleid", "name", "regexp", "splitchar", "unit")
   
   
   for (i in 1:dim(tags.df)[1]) {
      current.tag <- which(regexpr(tags.df$regexp[i], tifftags.clean) == 1) + 1
      value.tmp <- strsplit(tifftags.clean[current.tag], split = tags.df$splitchar[i])[[1]][2]
      # Remove leading spaces from tags.df$value
      tags.df$value[i] <- sub("^\\s+", "", value.tmp, useBytes = T)
      if (tags.df$unit[i] != "") {
         tags.df$value[i] <- paste("\\SI{", strsplit(tags.df$value[i], split = " ")[[1]][1], "}{", tags.df$unit[i], "}", sep = "")
      }
   }
      
   
   tags <- data.frame(stringsAsFactors = FALSE, 
                      sampleid = tags.df$sampleid, 
                      parameter = tags.df$name, 
                      value = tags.df$value)

   return(tags)
}