##################################################################
## Useful functions needed to run the developed methods in this ## 
## collection of code                                           ##
##################################################################

#############################################
# PAR Function Bootstrap & Jackknife Method #
#############################################

PAR2 <- function(a, b, c, d){
  n= a + b + c + d
  ((a + c)/n) - (c/(c + d))
}

#############################################################################
# Transforming by Fisher Z PAR function for Bootstrap and Jackknife Methods #
#############################################################################

Partrans <- function(PAR){
  TPar <- (1/2)*log((1+PAR)/(1-PAR))
}

######################################
# Undoing Fisher's Z transformation  #
######################################

Untrans <- function(z){
  un.trans <- (exp(2*z)-1)/(exp(2*z)+1)
}
