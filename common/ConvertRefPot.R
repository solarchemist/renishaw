source("/home/taha/chepec/chetex/common/R/common/SHE2AVS.R")
source("/home/taha/chepec/chetex/common/R/common/AVS2SHE.R")
source("/home/taha/chepec/chetex/common/R/common/ConvertRefPotEC.R")

##################################################
################# ConvertRefPot ##################
##################################################
ConvertRefPot <- function(argpotential, argrefscale, valuerefscale) {
   # Check that argpotential is valid numeric 
   
   #  IDEA: make a matrix out of these (scale names and flags)
   
   # Valid scales
   scale.names <- list()
   scale.names[["SHE"]] <- c("SHE", "NHE", "she", "nhe")
   scale.names[["AgCl"]] <- c("Ag/AgCl", "AgCl", "ag/agcl", "agcl")
   scale.names[["SCE"]] <- c("SCE", "sce")
   scale.names[["Li"]] <- c("Li/Li+", "Li", "Li+", "li", "li+", "li/li+")
   scale.names[["AVS"]] <- c("AVS", "avs")
   
   # Set flags
   bool.flags <- as.data.frame(matrix(0, nrow = length(scale.names), ncol = 2))
   names(bool.flags) <- c("argref", "valueref")
   row.names(bool.flags) <- names(scale.names)
   
   # argrefscale
   # Check that argrefscale is valid character mode
   # ...
   
   for (j in 1:length(row.names(bool.flags))) {
      if (any(scale.names[[row.names(bool.flags)[j]]] == argrefscale)) {
         bool.flags[row.names(bool.flags)[j], "argref"] <- j
      }
   }
   
   
   # valuerefscale
   # Check that valuerefscale is valid character mode
   # ...
   
   for (k in 1:length(row.names(bool.flags))) {
      if (any(scale.names[[row.names(bool.flags)[k]]] == valuerefscale)) {
         bool.flags[row.names(bool.flags)[k], "valueref"] <- k
      }
   }
   
   # Depending on which flags are set, call the corresponding function
   
   decision.vector <- colSums(bool.flags)
   
   # Check if both scales are the same (no conversion needed). If so, abort gracefully.
   # ...
   
   if (decision.vector["argref"] == 5 || decision.vector["valueref"] == 5) {
      # AVS is requested, deal with it it
      if (decision.vector["argref"] == 5) {
         # Conversion _from_ AVS
         rnpotential <- ConvertRefPotEC(AVS2SHE(argpotential), 
            "SHE", 
            scale.names[[decision.vector["valueref"]]][1])
      } 
      if (decision.vector["valueref"] == 5) {
         # Conversion _to_ AVS
         rnpotential <- SHE2AVS(ConvertRefPotEC(argpotential, 
            scale.names[[decision.vector["argref"]]][1], 
            "SHE"))
      }
   } else {
      rnpotential <- ConvertRefPotEC(argpotential, 
         scale.names[[decision.vector["argref"]]][1], 
         scale.names[[decision.vector["valueref"]]][1])
   }
   return(rnpotential)
}
