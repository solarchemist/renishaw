source("/home/taha/chepec/chetex/common/R/common.R")

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
   # The first 47 rows of tifftags contains info about
   # the image itself, and is not interesting
   # The last row of tifftags is the empty string
   trimtags <- 
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
      gsub("\\xb0", "", tifftags[48:length(tifftags) - 1], 
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
   oddtags <- trimtags[seq(1, length(trimtags), by = 2)]
   eventags <- trimtags[seq(2, length(trimtags), by = 2)]
   tmp.tagsdf <- 
      data.frame(stringsAsFactors = FALSE,
      sampleid = substrate.id,
      tag = oddtags,
      value = eventags)
   # A vector of displayed tags
   displaytags <- 
    c("AP_ACTUALKV",             # EHT = 5.00 kV
      "DP_HIGH_CURRENT",         # High Current = Off # commonly used for EDS analysis
      "AP_WD",                   # WD = 4.1 mm
      "AP_MAG",                  # Mag = 15.00 K X
      "AP_BRIGHTNESS",           # Brightness = 49.0 %
      "AP_CONTRAST",             # Contrast = 32.4 %
      "DP_DETECTOR_CHANNEL",     # InLens
      "DP_SCM_STATUS",           # SCM Status = Off
      "AP_SCM",                  # Specimen I = 0 fA # zero if SCM Status = Off
      # Beam
      #notrelevant#"AP_ACTUALCURRENT",        # Fil I = 2.370 A
      "AP_STIG_Y",               # Stigmation Y = -0.2 %
      "AP_STIG_X",               # Stigmation X = 1.9 %
      "AP_APERTURE_ALIGN_Y",     # Aperture Align Y = -1.5 % 
      "AP_APERTURE_ALIGN_X",     # Aperture Align X = -2.4 %
      "AP_APERTURESIZE",         # Aperture Size = 30.00 m
      "AP_BEAMSHIFT_Y",          # Beam Shift Y = 0.1%
      "AP_BEAMSHIFT_X",          # Beam Shift X = 0.1%
      #notrelevant#"AP_BEAM_OFFSET_Y",        # Beam Offset Y = 13.89 nm
      #notrelevant#"AP_BEAM_OFFSET_X",        # Beam Offset X = 9.13 nm
      # Stage
      "DP_TRACK_Z",              # Track Z = On
      "AP_STAGE_AT_Z",           # Stage at Z = 25.857 mm
      "AP_STAGE_AT_Y",           # Stage at Y  = 49.0708 mm
      "AP_STAGE_AT_X",           # Stage at X = 50.0753 mm
      # Tilt 
      "DP_STAGE_TILTED",         # Stage Tilt = In Y
      "AP_TILT_ANGLE",           # Tilt Angle = 0.0
      "AP_TILT_AXIS",            # Tilt Axis = 0.0
      # Image                 
      "DP_SCANRATE",             # Scan Speed = 10
      "AP_FRAME_TIME",           # Cycle Time = 48.7 Secs
      "DP_FREEZE_ON",            # Freeze on = End Frame
      "DP_DWELL_TIME",           # Dwell time = 100 ns
      "DP_NOISE_REDUCTION",      # Noise Reduction = Pixel Avg.
      "AP_FRAME_INT_COUNT",      # Frames to integrate = 0
      "AP_FRAME_AVERAGE_COUNT",  # Frames to average = 1
      #notrelevant#"AP_PROFILE_W",            # Profile Width = 71.79 m
      "AP_PIXEL_SIZE",           # Pixel Size = 23.96 nm
      # System
      "AP_COLUMN_VAC",           # Gun Vacuum = 1.74e-009 mbar
      "AP_SYSTEM_VAC",           # System Vacuum = 2.32e-006 mbar
      "AP_FILAMENT_AGE",         # Filament Age = 1378.35 Hours
      # Misc
      "AP_PHOTO_NUMBER",         # Photo No. = 9069
      "AP_DATE",                 # Date :23 Mar 2011 
      "AP_TIME")                 # Time :19:10:46
      #
      extracted.tagsdf <- data.frame()
      for (i in 1:length(displaytags)) {
         rows.extract <- which(displaytags[i] == tmp.tagsdf$tag)
         extracted.tagsdf <- rbind(extracted.tagsdf, 
           data.frame(stringsAsFactors = FALSE,
           tmp.tagsdf[rows.extract, ]))
      }
      # Fix split of time and date fields
      splittags <- strsplit(extracted.tagsdf$value, " = ")
      for (i in 1:length(splittags)) {
         if (length(splittags[[i]]) == 1) {
            splittags[[i]] <- strsplit(splittags[[i]], " :")[[1]]
         }
      }
      tagsdf <- data.frame()
      for (i in 1:length(splittags)) {
         tagsdf <- rbind(tagsdf,
            data.frame(stringsAsFactors = FALSE,
            sampleid = extracted.tagsdf$sampleid[i],
            parameter = splittags[[i]][1],
            value = splittags[[i]][2]))
      }
      # 
      return(tagsdf)
}