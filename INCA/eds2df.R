eds2df <- function(edstxtfile) {
   ## Description:
   ##   Reads EDS textfile from INCA EDS.
   ##   Stores data in data frame.
   ## Usage:
   ##   eds2df(edstxtfile)
   ## Arguments:
   ##    edstxtfile: character string, the full filename
   ##                (with path) to one txt file.
   ## Value:
   ##   A dataframe 
   #
   incatxt <- file(edstxtfile, "r")
   edsfile <- readLines(incatxt, n = -1) #read all lines of input file
   close(incatxt)
   #
   sampleid <- ProvideSampleId(edstxtfile)
   #
   rgxp.comment <- "^\\#"
   #
   numrow.idx <- regexpr(rgxp.comment, edsfile)
   # scrap the match.length attribute
   attr(numrow.idx, "match.length") <- NULL
   #
   i <- seq(1, length(numrow.idx) - 1, 1)
   j <- seq(2, length(numrow.idx), 1)
   # Start index of data range
   start.idx <- which(numrow.idx[i] == 1 & numrow.idx[j] != 1) + 1
   # End index of the data range
   end.idx <- which(numrow.idx[i] != 1 & numrow.idx[j] == 1)
   #
   zz <- textConnection(edsfile[start.idx:end.idx], "r")
   #
   ff <- data.frame()
   ff <- data.frame(stringsAsFactors = FALSE,
            sampleid = sampleid, 
            matrix(scan(zz, what = numeric(), sep = ","), ncol = 2, byrow = T))
   close(zz)
   names(ff) <- c("sampleid", "energy", "counts")
   #
   ### Collect attributes of this experiment
   # XUnit
   position.XUnit <- regexpr("^\\#XUNITS", edsfile)
   XUnit <- as.character(strsplit(edsfile[which(position.XUnit == 1)], ":\\s")[[1]][2])
   ff$XUnit <- XUnit
   # YUnit
   position.YUnit <- regexpr("^\\#YUNITS", edsfile)
   YUnit <- as.character(strsplit(edsfile[which(position.YUnit == 1)], ":\\s")[[1]][2])
   ff$YUnit <- YUnit
   # Date
   position.Date <- regexpr("^\\#DATE", edsfile)
   Date <- strsplit(edsfile[which(position.Date == 1)], ":\\s")[[1]][2]
   ff$Date <- Date
   # Time
   position.Time <- regexpr("^\\#TIME", edsfile)
   Time <- strsplit(edsfile[which(position.Time == 1)], ":\\s")[[1]][2]
   ff$Time <- Time
   # XPerChannel
   position.XPerChannel <- regexpr("^\\#XPERCHAN", edsfile)
   XPerChannel <- as.numeric(strsplit(edsfile[which(position.XPerChannel == 1)], ":\\s")[[1]][2])
   ff$XPerChannel <- XPerChannel
   # Offset
   position.Offset <- regexpr("^\\#OFFSET", edsfile)
   Offset <- as.numeric(strsplit(edsfile[which(position.Offset == 1)], ":\\s")[[1]][2])
   ff$Offset <- Offset
   # ChOffset
   position.ChOffset <- regexpr("^\\#CHOFFSET", edsfile)
   ChOffset <- as.numeric(strsplit(edsfile[which(position.ChOffset == 1)], ":\\s")[[1]][2])
   ff$ChOffset <- ChOffset
   # LiveTime
   position.LiveTime <- regexpr("^\\#LIVETIME", edsfile)
   LiveTime <- as.numeric(strsplit(edsfile[which(position.LiveTime == 1)], ":\\s")[[1]][2])
   ff$LiveTime <- LiveTime
   # DeadTime is calculated from: REALTIME - LIVETIME
   position.RealTime <- regexpr("^\\#REALTIME", edsfile)
   RealTime <- as.numeric(strsplit(edsfile[which(position.RealTime == 1)], ":\\s")[[1]][2])
   DeadTime <- RealTime - LiveTime
   ff$DeadTime <- DeadTime
   # BeamEnergy
   position.BeamEnergy <- regexpr("^\\#BEAMKV", edsfile)
   BeamEnergy <- as.numeric(strsplit(edsfile[which(position.BeamEnergy == 1)], ":\\s")[[1]][2])
   ff$BeamEnergy <- BeamEnergy   
   #
   return(ff)
}
