####################################################################################
#  Component-wise Metropolis-Hastings Independence Sampler (CWMH)                  #
#                                                                                  #
# This model incorporates information about the sensitivity and specificity of a   #
# diagnositc test and calculates the mean/median PAF and PAR with corresponding    #
# 95% CI and HPD intervals. Note that this function requires the installation of   #  
# the package coda to run.                                                         #
#                                                                                  #
# Function inputs:                                                                 #
# data - as a vector of the form (a,b,c,d) where we have:                          #  
#                                                                                  #
#       |D+ | D-                                                                   #
#    -------------                                                                 #
#    E+ | a | b                                                                    #  
#    E- | c | d                                                                    #
#                                                                                  #
# prior - a matrix of shape parameters for the priors of each parameter, in the    #
# order p,q,e,Se,Sp where each parameter has an alpha and beta shape parameter.    #
# Note that it is assumed for this program that all priors take the form of a beta #
# distribution.                                                                    #
#                                                                                  #
# B - the burnin period, by default this is set to 0                               #
#                                                                                  #
# S - the total number of iterates desired for the chain, default = 1000.          #
#                                                                                  #
# thin- the amount of thinning desired, default is no thining.                     #
#                                                                                  # 
# OUTPUT: a list where output$Stats provides the posterior mean, median, sd, Lower #
# and Upper bounds of the credible interval for the PAR and HPD Lower and HPD Upper#
# bounds.                                                                          #
#                                                                                  #
# output$Theta provides a dataframe containing output for the Marcov-chain for     #
# each parameter.                                                                  #
#                                                                                  #
# output$AcceptRate provides the acceptance rate for each parameter.               #
####################################################################################

CWMH <- function(data, prior, B=0, S=1000, thin=1){
    ## Specifying Se and Sp shape parameters ## 
    alpha.e <- prior[4,1]  
    beta.e <- prior[4,2]
    alpha.p <- prior[5,1]
    beta.p <- prior[5,2]
    
    ## Initalisation ##
    theta <- matrix(NA, nrow =1, ncol = 5)
    
    ## Calculating Inital Conditions ##
    repeat {
    Se <- rbeta(n=1, shape1=alpha.e, shape2=beta.e)
    Sp <- rbeta(n=1, shape1=alpha.p, shape2=beta.p)
    
    # Observed probabilities (described as eta_{11}, ..., eta_{22} in thesis)
    q11 <- data[1]/sum(data) 
    q12 <- data[2]/sum(data)
    q21 <- data[3]/sum(data)
    q22 <- data[4]/sum(data)
    
    # If using all 4 eqns for observed probabilities in your system rather than 3 make sure all 
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
      e <- p11 + p12 
      p <- p11/e
      q <- p21/(1-e) 
      theta[1,] <- c(p,q,e,Se,Sp)
      break
    }
  }
  
  r <- ncol(theta)                           # Total number of elements in theta
  prev <- theta[1,]
  
  for(t in 1:S){                              # t-loop, iteration counter
    print(t)
    for(i in 1:r){                            # i-loop, parameter indicator 1:5 = (p,q,e,Se,Sp)
      tt <- prev                              # theta(t,i), vector of current state
      tc <- tt                                # Let the candidate vector = the current state which for each [i] will 
      # be replaced with a candidate value 
      tc[i]<- rbeta(n=1,shape1=prior[i,1], shape2=prior[i,2])   # Candidate dist for all 5 parameters is a flat Beta(1,1)
      
      # Candidate probabilites for multinomial distribution
      p1.c <- tc[4]*tc[1]*tc[3] + (1-tc[5])*tc[2]*(1-tc[3])
      p2.c <- tc[4]*(1-tc[1])*tc[3] + (1-tc[5])*(1-tc[2])*(1-tc[3])
      p3.c <- tc[5]*tc[2]*(1-tc[3]) + (1-tc[4])*tc[1]*tc[3]
      p4.c <- tc[5]*(1-tc[2])*(1-tc[3]) + (1-tc[4])*(1-tc[1])*tc[3]
      
      # Probabilites for time step t of multinomial distribution
      p1 <- tt[4]*tt[1]*tt[3] + (1-tt[5])*tt[2]*(1-tt[3])
      p2 <- tt[4]*(1-tt[1])*tt[3] + (1-tt[5])*(1-tt[2])*(1-tt[3])
      p3 <- tt[5]*tt[2]*(1-tt[3]) + (1-tt[4])*tt[1]*tt[3]
      p4 <- tt[5]*(1-tt[2])*(1-tt[3]) + (1-tt[4])*(1-tt[1])*tt[3]
      
      log.lik.c <- dmultinom(data, prob=c(p1.c,p2.c,p3.c,p4.c), log=T) # log-lik candidate contribution                                                                        
      log.lik <- dmultinom(data, prob=c(p1,p2,p3,p4), log=T)           # log-lik previous timestep contribution
      alpha <-  min(1,exp(log.lik.c - log.lik))    # Probability of acceptance (simplified as implementing independence sampler)
      
      # accept or reject candidate value
      u <- runif(1)                           
      prev <- tt
      if(u <= alpha){
        prev <- tc        
      }
    }
    
    # Thinning and Burnin statement 
    if(t > B & t %% thin == 0){    # %% = modulus arg                   
      theta <- rbind(theta,prev, deparse.level=0)
    }  
  }
  theta <- data.frame(theta)
  colnames(theta) <- c("p","q", "e","Se","Sp")
  
  # acceptance and rejection rates
  require(coda)
  rejectR <- rejectionRate(as.mcmc(theta))
  acceptR <- 1- rejectR
  acceptR <- c(acceptR, PAR=NA, PAF=NA)
  
  # Statistics (mean, CIs, ect.) #
  p <- theta$p
  q <- theta$q
  e <- theta$e
  Se <- theta$Se
  Sp <- theta$Sp
  
  # HPD intervals #
  hpd.p <- HPDinterval(as.mcmc(p), prob=0.95)
  hpd.q <- HPDinterval(as.mcmc(q), prob=0.95)
  hpd.e <- HPDinterval(as.mcmc(e), prob=0.95)
  hpd.se <- HPDinterval(as.mcmc(Se), prob=0.95)
  hpd.sp <- HPDinterval(as.mcmc(Sp), prob=0.95)
  
  # Probabilites for all cross-classifications #
  p11 <- Se*p*e + (1-Sp)*q*(1-e)
  p12 <- Se*(1-p)*e + (1-Sp)*(1-q)*(1-e)
  p21 <- Sp*q*(1-e) + (1-Se)*p*e
  p22 <- Sp*(1-q)*(1-e) + (1-Se)*(1-p)*e
  
  # HPD intervals for cross-classification probabilities #
  hpd.p11 <- HPDinterval(as.mcmc(p11), prob=0.95) 
  hpd.p12 <- HPDinterval(as.mcmc(p12), prob=0.95)
  hpd.p21 <- HPDinterval(as.mcmc(p21), prob=0.95)
  hpd.p22 <- HPDinterval(as.mcmc(p22), prob=0.95)
  
  # PAR and PAF estimates #
  PAR <-  (p11 + p21) - (p21/(p21+p22))
  PAF <- PAR /(p11 + p21)
  theta <- cbind(theta,PAR,PAF)
  
  # HPD intervals #
  hpd.PAR <- HPDinterval(as.mcmc(PAR), prob=0.95) 
  hpd.PAF <- HPDinterval(as.mcmc(PAF), prob=0.95)
  
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
  HPDLower <- round(c(hpd.p[,1], hpd.q[,1], hpd.e[,1], hpd.se[,1], hpd.sp[,1], hpd.p11[,1], hpd.p12[,1],
                      hpd.p21[,1], hpd.p22[,1], hpd.PAR[,1], hpd.PAF[,1]), digits = 4)
  HPDUpper <- round(c(hpd.p[,2], hpd.q[,2], hpd.e[,2], hpd.se[,2], hpd.sp[,2], hpd.p11[,2], hpd.p12[,2],
                      hpd.p21[,2], hpd.p22[,2], hpd.PAR[,2], hpd.PAF[,2]), digits = 4)  
  Stats <- as.data.frame(cbind(mean, sd, Lower, median, Upper, HPDLower, HPDUpper, sample),
                         row.names=Parameter.Names)
  Results <- list("Stats" = Stats, "Theta" = theta, "AcceptRate" = acceptR)  
  return(Results)
}  

## Running the CWMH Function ##
data <- c(22,25,82,251) # lepto data
prior = matrix(c(1, 1, 1, 1, 2, 2, 25, 3, 30, 1.5), nrow=5, ncol=2, byrow=T) # Setting shape values for priors order: c(p,q,e,Se,Sp)

Results <- CWMH(data=data, prior=prior, B=4000, S=10000)

