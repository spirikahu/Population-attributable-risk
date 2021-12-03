#############################################################
## Estimating a credible interval for the PAR where a      ##
## Bayesian approach has being adopted.                    ##
## Note that before this function will run the user must   ##
## first source Useful_Functions.R and install the packages##
## MCMC pack and coda.                                     ##
##                                                         ##
## Function inputs: a, b, c and d which represent observed ##
## values in the following 2 x 2 table                     ##
##       |D+ | D-                                          ##
##    -------------                                        ##
##    E+ | a | b                                           ##
##    E- | c | d                                           ##
##                                                         ##
## nMC represents the number of iterations in the Markov   ##
## chain. Given we have an analytic solution for the       ##  
## posterior here there is no need for burnin.             ##
##                                                         ##
## dp is a vector of length four which represents the      ##
## Dirichlet prior on the probability of each of the four  ##
## cross-classifications in the table in this order:       ##
## (x_{11}, x_{12}, x_{21}, x_{22}).                       ##
## by default the prior is selected as the standard        ##
## reference prior Dirichlet(1,1,1,1).                     ##
##                                                         ##
## Function output: PAR, Lower & Upper CI                  ##
#############################################################

# BayesCI(a=22, b=25, c=82, d=251)

BayesCI <- function(a, b, c, d, nMC = 10000, 
                    conf=0.95, dp=c(1,1,1,1)){
  require(MCMCpack)
  MC <- nMC
  x <- c(a,b,c,d)     
  alpha <- dp                     
  posterior.alpha <- x + alpha    
  posterior.iterates <- rdirichlet(MC, posterior.alpha)
  p1 <- posterior.iterates[,1]
  p2 <- posterior.iterates[,2]
  p3 <- posterior.iterates[,3]
  p4 <- posterior.iterates[,4]
  PAR.iterates <- ((p1 + p3)/(p1 + p2 + p3 + p4))-(p3/(p3 + p4))
  PAR.mean <- mean(PAR.iterates)
  PAR.sd <- sd(PAR.iterates)
  PAR.median <- median(PAR.iterates)
  PAR.PI <- quantile(PAR.iterates, c((1-conf)/2, 0.5+conf/2))
  BayesRes <- data.frame(cbind(PAR.mean, PAR.median, PAR.sd, Lower=PAR.PI[1], Upper=PAR.PI[2]), row.names=""); BayesRes
}
