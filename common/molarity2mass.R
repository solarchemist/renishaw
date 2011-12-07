##################################################
############### molarity2mass ####################
##################################################
molarity2mass <- function(formulamass, volume, molarity) {
   # Calculates the required mass of
   # the substance to be dissolved.
   # ARGS: formulamass - formula mass of the substance (in gram per mole)
   #       volume      - volume of the final solution (in liters)
   #       molarity    - molarity (in moles per liter)
   # VALUE: mass of substance (in grams)
   #
   mass <- formulamass * volume * molarity 
   # Unit check:
   # [g * mol-1] * [liter] * [mole * liter-1] = [g]
   return(mass)
}
