####################################################################################
#  Metropolis Adjusted Langevin Algorithm (MALA)                                   #
#                                                                                  #
# This model incorporates information about the sensitivity and specificity of a   #
# diagnositc test and calculates the mean/median PAF and PAR with corresponding    #
# 95% CI intervals. Note that this function requires the installation of the       #  
# package coda to run.                                                             #
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
# sigma - scale factor tau described in thesis, a small positive constant.         #                                                                   #
#                                                                                  #
# IC - a prespecified vector of inital conditions for each parameter in the order  #
# (p,q,e,Se,Sp). If not specified IC will be calculated.                           #
#                                                                                  #
# OUTPUT: a list where output$Stats provides the posterior mean, median, sd, Lower #
# and Upper bounds of the credible interval for the PAR.                           #
#                                                                                  #
#                                                                                  #
# output$Theta provides a dataframe containing output for the Marcov-chain for     #
# each parameter.                                                                  #
#                                                                                  #
# output$AcceptRate provides the acceptance rate for each parameter.               #
####################################################################################

MALA <- function(data, S=1000, B=0, thin=1, prior, sigma, IC, ...){
  ## Defining elements of the prior ## 
  alpha.p <- prior[1,1]; beta.p <- prior[1,2]
  alpha.q <- prior[2,1]; beta.q <- prior[2,2]
  alpha.e <- prior[3,1]; beta.e <- prior[3,2]
  alpha.Se <- prior[4,1]; beta.Se <- prior[4,2]
  alpha.Sp <- prior[5,1]; beta.Sp <- prior[5,2]
  
  ## Initialization ## 
  theta <- matrix(NA, nrow =1, ncol = 5)
  thetasave <- matrix(NA, nrow =1, ncol = 5)
  
  ## Inital Conditions ##
  # In this function you can pre-specify inital conditions. If this is not done
  # then the same approach that was applied in CWMH/RWMH/Gibbs applied to calculate
  # the inital conditions
  if(missing(IC)){
    repeat {
      Se <- rbeta(n=1, shape1=alpha.Se, shape2=beta.Se)
      Sp <- rbeta(n=1, shape1=alpha.Sp, shape2=beta.Sp)
      
      # Observed probabilities (described as eta_{11}, ..., eta_{22} in thesis)
      q11 <- data[1]/sum(data) 
      q12 <- data[2]/sum(data)
      q21 <- data[3]/sum(data)
      q22 <- data[4]/sum(data)
      # If you use all 4 eqns in your system rather than 3 you also need to check
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
        theta[1,] <- c(p,q,e,Se,Sp)
        break
      }
    }
  }else{
    theta[1,] <- IC
    thetasave[1,] <- IC
  }
  
  ##########################
  ## Vector Norm Function ##
  ##########################
  norm_vec <- function(x) sqrt(sum(x^2))
  
  ########################################################################
  ## Gradient function for MALA specified by eqns (3.27-3.31) in thesis ##
  ########################################################################
  grad <- function(data, theta, prior){ 
    ## Unpacking Inputs ##
    x11 <- data[1]
    x12 <- data[2]
    x21 <- data[3]
    x22 <- data[4]
    
    p <- theta[1]
    q <- theta[2]
    e <- theta[3]
    Se <- theta[4]
    Sp <- theta[5]
    
    alpha.p <- prior[1,1]; beta.p <- prior[1,2]
    alpha.q <- prior[2,1]; beta.q <- prior[2,2]
    alpha.e <- prior[3,1]; beta.e <- prior[3,2]
    alpha.Se <- prior[4,1]; beta.Se <- prior[4,2]
    alpha.Sp <- prior[5,1]; beta.Sp <- prior[5,2]
    
    ## Multinomial Probabilities ##
    q11 <- Se*p*e + (1-Sp)*q*(1-e)
    q12 <- Se*(1-p)*e + (1-Sp)*(1-q)*(1-e)
    q21 <- Sp*(1-e)*q + (1-Se)*p*e
    q22 <- (1-Se)*(1-p)*e + Sp*(1-q)*(1-e) 
    
    # The derivatives of log(f(theta))
    dp <- ((e*Se*x11)/q11) - ((Se*x12*e)/q12) + (((1-Se)*e*x21)/q21) + ((e*(Se-1)*x22)/q22) +
      ((alpha.p-1)/p) - ((beta.p-1)/(1-p))       
    
    dq <- (((1-Sp)*(1-e)*x11)/q11) + (((1-Sp)*(e-1)*x12)/q12) + ((Sp*(1-e)*x21)/q21)+
      ((Sp*(e-1)*x22)/q22) + (alpha.q-1)/q - ((beta.q-1)/(1-q))
    
    de <- (((Se*p + (Sp-1))*x11)/q11) + (((Se*(1-p) + (1-Sp)*(q-1))*x12)/q12) +
      ((((1-Se)*p - Sp*q)*x21)/q21) + ((((1-Se)*(1-p) + Sp*(q-1))*x22)/q22) +
      ((alpha.e - 1)/e) - ((beta.e -1)/(1-e))
    
    dSe <- ((p*e*x11)/q11) + (((1-p)*e*x12)/q12) - ((p*e*x21)/q21) + ((e*(p-1)*x22)/q22) +
      ((alpha.Se-1)/Se) - ((beta.Se-1)/(1-Se))
    
    dSp <- ((q*(e-1)*x11)/q11) + (((1-q)*(e-1)*x12)/q12) + (((1-e)*q*x21)/q21) + 
      (((1-q)*(1-e)*x22)/q22) + ((alpha.Sp-1)/Sp) - ((beta.Sp -1)/(1-Sp))
    grad <- matrix(c(dp, dq, de, dSe, dSp), nrow=5, ncol=1)
    return(grad)
  }
  
  ##############################
  ## Langevin Algorithm Start ## 
  ##############################
  
  for(t in 1:S){   # t-loop, iteration counter
    print(t)
    sigma <- sigma # Initial value of sigma (scalar)
    
    ## Selecting Candidate Value ##  
    theta <- rbind(theta,t(matrix(theta[t,]) + ((sigma^2)/2) * grad(data, theta[t,],prior) + sigma * matrix(mvrnorm(n=1, mu = rep(0,5), Sigma = diag(5)), ncol=1, nrow=5)))
    # The stochastic equation given above can be represented by the Gaussian distribution (3.6) in thesis
    
    if(length(which(theta[t+1,] <= 0 | theta[t+1,] >= 1)) > 0){
      theta[t+1,] <- theta[t,] # if any of the probabilities fall outside [0,1] reject 
    } else{
      #############################
      ## Likelihood Contribution ##
      #############################
      
      # Candidate Values #
      p.c <- theta[t+1,1]; q.c <- theta[t+1,2]
      e.c <- theta[t+1,3]; Se.c <- theta[t+1,4]
      Sp.c <- theta[t+1,5]                      
      
      p1.c <- Se.c*p.c*e.c + (1-Sp.c)*q.c*(1-e.c)
      p2.c <- Se.c*(1-p.c)*e.c + (1-Sp.c)*(1-q.c)*(1-e.c)
      p3.c <- Sp.c*(1-e.c)*q.c + (1-Se.c)*p.c*e.c
      p4.c <- (1-Se.c)*(1-p.c)*e.c + Sp.c*(1-q.c)*(1-e.c)
      
      log.lik.c <- dmultinom(data, prob=c(p1.c,p2.c,p3.c,p4.c), log=T)
      
      # Current Timestep Values #
      p <- theta[t,1]; q <- theta[t,2]
      e <- theta[t,3]; Se <- theta[t,4]
      Sp <- theta[t,5]                          
      
      p1 <- Se*p*e + (1-Sp)*q*(1-e)
      p2 <- Se*(1-p)*e + (1-Sp)*(1-q)*(1-e)
      p3 <- Sp*(1-e)*q + (1-Se)*p*e
      p4 <- (1-Se)*(1-p)*e + Sp*(1-q)*(1-e)
      
      log.lik <- dmultinom(data, prob=c(p1,p2,p3,p4), log=T)
      
      ########################
      ## Prior Contribution ##
      ########################
      
      # Initialisation #
      log.pr.c <- matrix(NA, nrow=1, ncol=5) 
      log.pr <- matrix(NA, nrow=1, ncol=5)
      
      # assumption made that priors for all parameters are beta distributed
      for(i in 1:5){
        log.pr.c[i] <- dbeta(theta[t+1,i], shape1=prior[i,1], shape2=prior[i,2], log=T)
        log.pr[i] <- dbeta(theta[t,i], shape1=prior[i,1], shape2=prior[i,2], log=T)
      }
      
      ###########################
      ## Proposal Contribution ##
      ###########################
      x <- theta[t+1,] - theta[t,] - (sigma^2/2)*grad(data=data, theta=theta[t,], prior=prior)
      y <- theta[t,] - theta[t+1,] - (sigma^2/2)*grad(data=data, theta=theta[t+1,], prior=prior)
      proposal.num <- exp((-(norm_vec(y))^2)/(2*sigma^2))
      proposal.dom <- exp((-(norm_vec(x))^2)/(2*sigma^2))
      ratio.proposal <- proposal.num/proposal.dom
      
      ############################
      ## Acceptance Probability ##
      ############################
      accept.prob <- exp((sum(log.pr.c) + log.lik.c) - (sum(log.pr) + log.lik))* ratio.proposal
      
      ## Updating Process ##  
      alpha <- min(1, accept.prob)
      u <- runif(1)                          
      
      if(u > alpha){   
        # If u <= alpha the current value becomes the candidate. 
        # Therefore, if u > alpha the current value does not become the candidate value.
        theta[t+1,] <- theta[t,]        
      }
    }
    # Thinning and Burnin statement 
    if(t > B && t %% thin == 0){    # %% = modulus arg                   
      thetasave <- data.frame(rbind(thetasave, theta[t+1,], deparse.level=0))
    }
  } 
  colnames(thetasave) <- c("p","q", "e","Se","Sp")
  thetasave <- thetasave[-1,] # removing first row of NA's
  ## acceptance and rejection rates ##
  require(coda)
  rejectR <- rejectionRate(as.mcmc(thetasave))
  acceptR <- 1- rejectR
  acceptR <- c(acceptR, PAR=NA, PAF=NA)  
  
  # Statistics (mean, CIs, ect.) #
  p <- thetasave$p
  q <- thetasave$q
  e <- thetasave$e
  Se <- thetasave$Se
  Sp <- thetasave$Sp
  
  # Probabilites for all cross-classifications #
  p11 <- Se*p*e + (1-Sp)*q*(1-e)
  p12 <- Se*(1-p)*e + (1-Sp)*(1-q)*(1-e)
  p21 <- Sp*q*(1-e) + (1-Se)*p*e
  p22 <- Sp*(1-q)*(1-e) + (1-Se)*(1-p)*e
  
  # PAR and PAF estimates #
  PAR <-  (p11 + p21) - (p21/(p21+p22))
  PAF <- PAR /(p11 + p21)
  thetasave <- cbind(thetasave,PAR,PAF)
  
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
  Results <- list("Stats" = Stats, "Theta" = thetasave, "AcceptRate" = acceptR)
  return(Results)
  
}

## Running the MALA Function ##
data <- c(22,25,82,251) # lepto data
prior = matrix(c(1, 1, 1, 1, 2, 2, 25, 3, 30, 1.5), nrow=5, ncol=2, byrow=T) # Setting shape values for priors order: c(p,q,e,Se,Sp)

Results <- MALA(data=data, prior=prior, sigma=1, S=100, B=10)
