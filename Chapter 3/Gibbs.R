####################################################################################
#  Gibbs Sampler (Gibbs)                                                           #
#                                                                                  #
# This model incorporates infromation about the sensitivity and specificity of a   #
# diagnositc test and calculates the mean/median PAF and PAR with corresponding    #
# CI intervals. Note that this function requires the installation of the           #
# package MCMCpack to run.                                                         #
#                                                                                  #
# Function inputs:                                                                 #
# data - as a vector of the form (a,b,c,d) where we have:                          #  
#                                                                                  #
#       |D+ | D-                                                                   #
#    -------------                                                                 #
#    E+ | a | b                                                                    #  
#    E- | c | d                                                                    #
#                                                                                  #                                                                 #
# prior - a matrix of shape parameters for the priors of each parameter, in the    #
# order p,q,e,Se,Sp where each parameter has an alpha and beta shape parameter.    #
# Note that it is assumed for this program that all priors take the form of a beta #
# distribution. #                                                                  #
# B - the burnin period, by default this is set to 0                               #
#                                                                                  #
# S - the total number of iterates desired for the chain, default = 1000.          #
#                                                                                  #
# thin- the amount of thinning desired, default is no thining.                     #
#                                                                                  #                                                                  #
#                                                                                  #
# OUTPUT: a list where output$Stats provides the posterior mean, median, sd, Lower #
# and Upper bounds of the credible interval for the PAR and HPD Lower and HPD Upper#
# bounds.                                                                          #
#                                                                                  #
# output$Theta provides a dataframe containing output for the Marcov-chain for     #
# each parameter.                                                                  #
#                                                                                  #
####################################################################################

Gibbs <- function(data, prior,S=1000, B=0, thin=1){
  ## Specifying Se and Sp shape parameters ## 
  alpha.e <- prior[4,1]  
  beta.e <- prior[4,2]
  alpha.p <- prior[5,1]
  beta.p <- prior[5,2]
  
  ## Initialisation ##
  theta <- data.frame(matrix(NA, nrow =1, ncol = 11))
  
  ## Re-assigning the data ##
  X11 <- data[1]
  X12 <- data[2]
  X21 <- data[3]
  X22 <- data[4]
  
  ## Inital Conditions ##
  repeat {
    set.seed(123)
    Se <- rbeta(n=1, shape1=alpha.e, shape2=beta.e)
    Sp <- rbeta(n=1, shape1=alpha.p, shape2=beta.p)
    
    # Observed probabilities (described as eta_{11}, ..., eta_{22} in thesis)
    q11 <- data[1]/sum(data) 
    q12 <- data[2]/sum(data)
    q21 <- data[3]/sum(data)
    q22 <- data[4]/sum(data)
    
    # If you use all 4 eqns in your system rather than 3 you also need to test that all the 
    # probabilites sum to 1. This is accounted for in the case where you only consider the 
    # 3 independent eqns but not here. 
    
    # specifying and solving the set of simultaneous equations (3.1) in thesis 
    a <- matrix(c(Se,0,(1-Sp),0,0,Se,0,(1-Sp),(1-Se),0,Sp,0,0,(1-Se),0,Sp),nrow=4,ncol=4, byrow=TRUE) 
    
    b <- matrix(c(q11,q12,q21,q22),nrow=4,ncol=1)
    
    soln <- solve(a,b)
    # true probabilities (described as pi_{11},..,pi_{22} in thesis)
    p11 <- soln[1,]
    p12 <- soln[2,]
    p21 <- soln[3,]
    p22 <- soln[4,]
    sum <- p11 + p12 + p21 + p22
    
    if(length(which(soln >= 0 & soln <= 1)) >= 4 && (sum = 1)){
      e <-  p11 + p12 
      p <- p11/e
      q <- p21/(1-e) 
      PAR <-  (p11 + p21) - (p21/(p21+p22))
      PAF <- PAR/(p11 + p21)
      theta[1,] <- c(p11,p12,p21,p22,p,q,e,Se,Sp,PAR,PAF)
      colnames(theta) <- c("p11","p12","p21","p22", "p", "q", "e", "Se","Sp","PAR", "PAF")
      break
    }
  }
  
  Se <- theta$Se
  Sp <- theta$Sp 
  p11 <- theta$p11
  p12 <- theta$p12
  p21 <- theta$p21
  p22 <- theta$p22
  
  ## Begining of the Gibbs Sampling procedure ##
  for(i in 1:S){  
    print(i)
    # Full conditional posterior distributions for Y and Z 
    # (latent variables for correctly and incorrectly classified respectively)
    Y11 <- sum(rbinom(n=X11, size=1, prob=Se*p11/(Se*p11+(1-Sp)*p21)))
    Z21 <- X11-sum(Y11)
    
    Y12 <- sum(rbinom(n=X12, size=1, prob=Se*p12/(Se*p12+(1-Sp)*p22)))
    Z22 <- X12-sum(Y12)
    
    Y21 <- sum(rbinom(n=X21, size=1, prob=Sp*p21/(Sp*p21 + (1-Se)*p11)))
    Z11 <- X21 - sum(Y21)
    
    Y22 <- sum(rbinom(n=X22, size=1, prob=Sp*p22/((1-Se)*p12 + Sp*p22)))
    Z12 <- X22 - sum(Y22)
    
    # Full conditional posteriors for Se and Sp dependence on Y & Z
    Se <- rbeta(n=1, shape1=Y11 + Y12 + alpha.e, shape2= Z11 + Z12 + beta.e) 
    Sp <- rbeta(n=1, shape1=Y21 + Y22 + alpha.p, shape2 = Z21 + Z22 + beta.p)
    
    require(MCMCpack)
    # Full conditional posterior for vector of true probabilities (described as pi in thesis)
    p.vec <- rdirichlet(n=1, alpha=c(Y11+Z11+1, Y12+Z12+1, Y21+Z21+1, Y22+Z22+1))
    # assigning true probabilities for each cross-classification in 2x2 table
    p11 <- p.vec[,1]; p12 <- p.vec[,2]; p21 <- p.vec[,3]; p22 <- p.vec[,4]
    
    # Calculating parameters of interest
    e <-  p11 + p12 
    p <- p11/e
    q <- p21/(1-e) 
    
    PAR <- e*(p-q)
    PAF <- PAR/(p*e + q*(1-e))
    
    res <- data.frame(p11,p12,p21,p22,Se,Sp, p, q, e, PAR, PAF)
    
    ## Thinning and Burnin statement ##
    if(i >= B & i %% thin == 0){    # %% = modulus arg                   
      theta <- rbind(theta,res, deparse.level=0)
    }
  }
  theta <- data.frame(theta)
  colnames(theta) <- c("p11","p12","p21","p22","Se","Sp", "p", "q", "e",
                       "PAR", "PAF")
  
  # Statistics (mean, CIs, ect.) #
  p11 <- theta$p11
  p12 <- theta$p12
  p21 <- theta$p21
  p22 <- theta$p22
  p <- theta$p
  q <- theta$q
  e <- theta$e
  Se <- theta$Se
  Sp <- theta$Sp
  PAR <- theta$PAR
  PAF <- theta$PAF
  
  Parameter.Names <- c("p", "q", "e","Se", "Sp","p11", "p12", "p21", "p22", "PAR", "PAF")
  mean <- round(c(mean(p), mean(q), mean(e), mean(Se), mean(Sp), mean(p11), mean(p12), mean(p21), mean(p22),
                  mean(PAR), mean(PAF)), digits = 4)
  sd <- round(c(sd(p), sd(q), sd(e), sd(Se), sd(Sp), sd(p11), sd(p12), sd(p21), sd(p22) 
                , sd(PAR), sd(PAF)), digits = 4)
  median <- round(c(median(p), median(q), median(e), median(Se), median(Sp),
                    median(p11), median(p12), median(p21), median(p22), median(PAR), median(PAF)), digits = 4)
  Lower <- round(c(quantile(p,0.025), quantile(q,0.025), quantile(e,0.025), quantile(Se,0.025),
                   quantile(Sp,0.025),quantile(p11,0.025), quantile(p12,0.025), quantile(p21,0.025),
                   quantile(p22,0.025), quantile(PAR,0.025), quantile(PAF,0.025)),digits = 4)
  Upper <-round(c(quantile(p,0.975), quantile(q,0.975), quantile(e,0.975), quantile(Se,0.975),
                  quantile(Sp,0.975),quantile(p11,0.975),quantile(p12,0.975), quantile(p21,0.975), 
                  quantile(p22,0.975), quantile(PAR,0.975), quantile(PAF,0.975)),digits = 4)
  Stats <- as.data.frame(cbind(mean, sd, Lower, median, Upper),
                         row.names=Parameter.Names)
  Results <- list("Stats" = Stats, "Theta" = theta)
  return(Results)
}

## Running the Gibbs Function ##
data <- c(22,25,82,251) # lepto data
prior = matrix(c(1, 1, 1, 1, 2, 2, 25, 3, 30, 1.5), nrow=5, ncol=2, byrow=T) # Setting shape values for priors order: c(p,q,e,Se,Sp)

Results <- Gibbs(data=data, prior=prior, S=10)
