####################################################################################
#  Importance Sampler for non-identifiable models (Imp)                            #
#                                                                                  #
# This model incorporates information about the sensitivity and specificity of a   #
# diagnositc test and calculates the mean/median PAF and PAR with corresponding    #
# 95% CI and HPD intervals. Note that this function requires the installation of   #  
# the packages coda and MCMCpack to run.                                           #
#                                                                                  #
# Function inputs:                                                                 #
# data - as a vector of the form (a,b,c,d) where we have:                          #  
#                                                                                  #
#       |D+ | D-                                                                   #
#    -------------                                                                 #
#    E+ | a | b                                                                    #  
#    E- | c | d                                                                    #
#                                                                                  #
# S - the total number of iterations desired                                       #
#                                                                                  #
# alpha.e & beta.e - Shape parameters for the Se prior to begin initalisation      #
#                                                                                  #
# alpha.p & beta.p - Shape parameters for the Sp prior to begin initalisation      #
#                                                                                  #
####################################################################################                                                                               

Imp <- function(data, S=1000, alpha.e, beta.e, alpha.p, beta.p){
  require(coda)
  require(MCMCpack)
  
  ## Re-labelling the data ##
  X11 <- data[1]
  X12 <- data[2]
  X21 <- data[3]
  X22 <- data[4]
  
  ## Initalisation ##
  FRes <- NULL
  Count <- 0
  Total <- 0
  q <- NULL
  
  ##############################
  ## Importance sampler start ##
  ##############################
  
  while(Total < S){ # In while loop switch Total to Count if want S accepted values
    Total <- Total + 1
    print(Total)
    
    # Draw from priors for Se and Sp 
    Se <- rbeta(n=1, shape1=alpha.e, shape2=beta.e)
    Sp <- rbeta(n=1, shape1=alpha.p, shape2=beta.p)
    # Convience prior on q vector 
    # The q's are not fixed here as they were in the inf case
    q <- rdirichlet(n=1, alpha=c(X11 + 1, X12 + 1, X21 + 1, X22 + 1)) 
    q11 <- q[,1]
    q12 <- q[,2]
    q21 <- q[,3]
    q22 <- q[,4]   
    
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
      Count <- Count + 1
      e <- p11 + p12 
      p <- p11/e
      q <- p21/(1-e)
      PAR <- e*(p-q)
      PAF <- PAR/(p*e + q*(1-e))
      
      Res <- cbind(p11, p12, p21, p22, p, q, e, Se, Sp, PAR, PAF)
      FRes <- rbind(FRes, Res)
    }
  }  
  theta <- data.frame(FRes)
  colnames(theta) <- c("p11", "p12", "p21", "p22", "p", "q",
                         "e", "Se", "Sp", "PAR", "PAF")
  
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
  Results <- list(Stats = Stats, Theta = theta, Accept = 1-rejectionRate(mcmc(theta)))
  return(Results)
}

## Running the Imp function ##
data <- c(22,25,82,251) # lepto data

Results <- Imp(data=data,S=10,alpha.e=25,beta.e=3,alpha.p=30,beta.p=1.5)
