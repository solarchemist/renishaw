# Function loads R-data file into a variable instead of into the workspace
# Works well when the R-data file contains only ONE variable
# NOT TESTED for when the R-data file contains many variables
LoadRData2Variable <- function(FullPathToRData) {
   return(eval(parse(text = load(FullPathToRData))))
}
