###################################################### 
## In Pirikahu(2016) we compare the performance of  ##
## developed methods for estimating CIs for the PAR ##
## through simulation. The following code provides  ##
## the functions required to carryout simulations   ##
## for each of the discussed methods.               ##
##                                                  ##
## Input variables: p=P(D+|E+), q=P(D+|E-), e=P(E+),##
## n = sample size, N=number of iterations for the  ##
## simulation.                                      ##
##                                                  ##
######################################################


###########################
# Delta Method Simulation # 
###########################

Deltasim <- function(p, q, e, n, N){
  Count <<- Count + 1  # Count must be assigned as a global variable first
  print(Count)
  DelRes <- NULL       # Initalizing 
  a <- rep(0,N)
  b <- rep(0,N)
  c <- rep(0,N)
  d <- rep(0,N)
  
  for (i in 1:N){
    c.table <- ctable(p,q,e,n)
    a[i] <- c.table[1,]
    b[i] <- c.table[3,]
    c[i] <- c.table[2,]
    d[i] <- c.table[4,]
    
    Res <- DeltaPAR(a[i],b[i],c[i],d[i])
    DelRes <- rbind(DelRes, Res)
  }
  UpperCI <- DelRes[,1] + 1.96*DelRes[,2]
  LowerCI <- DelRes[,1] - 1.96*DelRes[,2]
  UpperCI.trans <- DelRes[,3]+1.96*DelRes[,4]
  LowerCI.trans <- DelRes[,3]-1.96*DelRes[,4]
  
  # Un-transforming Fishers CIs
  UpperCI.Fisher <- Untrans(UpperCI.trans)
  LowerCI.Fisher <- Untrans(LowerCI.trans)
  
  Delta <- as.data.frame(cbind(DelRes, a,b,c,d, LowerCI, UpperCI, LowerCI.Fisher, UpperCI.Fisher));Delta
}

###############################
# Jackknife Method Simualtion #
###############################

JackSim <- function(p, q, e, n, N){
  
  JakRes <- NULL    				# Initialising 
  a <- rep(0,N)		
  b <- rep(0,N)
  c <- rep(0,N)
  d <- rep(0,N)
  for(i in 1:N){						# calculating cells of contingency table
    c.table <- ctable(p,q,e,n)
    a[i] <- c.table[1,]
    b[i] <- c.table[3,]
    c[i] <- c.table[2,]
    d[i] <- c.table[4,]
    
    ParO <- PAR2(a=a[i], b=b[i], c=c[i], d=d[i])      # Par original table (MLE estimate of PAR)
    ParO.trans <- Partrans(ParO)                      # Transforming PAR original table (MLE estimate of transformed PAR) 
    Par1 <- PAR2(a=a[i]-1, b=b[i], c=c[i], d=d[i])		# Par for Jackknife sample 1
    Par2 <- PAR2(a=a[i], b=b[i]-1, c=c[i], d=d[i])		# Par for Jackknife sample 2
    Par3 <- PAR2(a=a[i], b=b[i], c=c[i]-1, d=d[i])		# Par for Jackknife sample 3
    Par4 <- PAR2(a=a[i], b=b[i], c=c[i], d=d[i]-1)		# Par for Jackknife sample 4
    
    Par1.trans <- Partrans(Par1)                      # Transformed Jackknife samples Par
    Par2.trans <- Partrans(Par2)
    Par3.trans <- Partrans(Par3)
    Par4.trans <- Partrans(Par4)
    
    meanPar <- (a[i]*Par1 + b[i]*Par2 + c[i]*Par3 + d[i]*Par4)/n
    var <-(a[i]*((Par1-meanPar)^2)+ b[i]*((Par2-meanPar)^2) + c[i]*((Par3-meanPar)^2) + d[i]*((Par4-meanPar)^2))*(n/(n-1))
    
    meanPar.trans <- (a[i]*Par1.trans + b[i]*Par2.trans + c[i]*Par3.trans + d[i]*Par4.trans)/n
    var.trans <-(a[i]*((Par1.trans-meanPar.trans)^2)+ b[i]*((Par2.trans-meanPar.trans)^2) + c[i]*((Par3.trans-meanPar.trans)^2) + d[i]*((Par4.trans-meanPar.trans)^2))*(n/(n-1))
    
    se <- sqrt(var)
    se.trans <- sqrt(var.trans)
    
    JakRes <- rbind(JakRes, c(ParO, se, ParO.trans, se.trans, a[i],b[i],c[i],d[i]))
    colnames(JakRes) <- c("Par", "se", "Par.trans", "se.trans", "a", "b", "c", "d")	
  }
  UpperCI <- JakRes[,1] + 1.96*(JakRes[,2])
  LowerCI <- JakRes[,1] - 1.96*(JakRes[,2])
  
  UpperCI.trans <- JakRes[,3] + 1.96*(JakRes[,4])
  LowerCI.trans <- JakRes[,3] - 1.96*(JakRes[,4]) 
  
  # Un-transforming Fishers CIs
  UpperCI.Fisher <- Untrans(UpperCI.trans)
  LowerCI.Fisher <- Untrans(LowerCI.trans)
  
  Jackknife <- as.data.frame(cbind(JakRes,LowerCI, UpperCI,  LowerCI.trans, UpperCI.trans, LowerCI.Fisher, UpperCI.Fisher))
}

###############################
# Bootstrap Simulation Method # 
###############################

Bootsim <- function(p, q, e, n, N){
  K <- 1000         # Setting number of Bootstrap replicates
  Res <- NULL			  # Initalising Res
  Par <- rep(0,K)		# Initalising PAR
  a <- rep(0,N)
  b <- rep(0,N)
  c <- rep(0,N)
  d <- rep(0,N)
  
  for(i in 1:N){
    c.table <- ctable(p,q,e,n)  # Formulation table 
    a[i] <- c.table[1,]
    b[i] <- c.table[3,]
    c[i] <- c.table[2,]
    d[i] <- c.table[4,]
    
    Par.O <- PAR2(a[i],b[i],c[i],d[i])     # MLE estimate
    Par.O.trans <- Partrans(Par.O)
    prob2 <- c(a[i]/n, b[i]/n, c[i]/n, d[i]/n)
    for(j in 1:K){				# K = Number of Bootstrap Samples
      repeat{
        c.table <- rmultinom(n = 1, size = n, prob = prob2)
        if(c.table[3,] + c.table[4,] > 0){
          break
        }  
      }
      a <- c.table[1,]
      b <- c.table[2,]
      c <- c.table[3,]
      d <- c.table[4,]
      Par[j] <- PAR2(a,b,c,d)
    }
    # Point estimates
    Par.trans <- Partrans(Par)            # MLE estimate for Fisher transformed PAR
    Percentile.mean.PAR <- mean(Par)      # point estimate of the mean from % Bootstrap distribution of PAR
    Percentile.median.PAR <- median(Par)  # point estimate of the median from % Bootstrap distribution of PAR 
    
    PercentileCIs <-  quantile(Par, c(0.5,0.025,0.975))
    PercentLCI <- PercentileCIs[2] 			# Non-symmetric CIs
    PercentUCI <-  PercentileCIs[3] 
    
    var <- sum((Par-mean(Par))^2)/(K-1)
    se <- sqrt(var)
    
    var.trans <- sum(((Par.trans)-mean(Par.trans))^2)/(K-1)
    se.trans <- sqrt(var.trans)
    
    UpperCI <- Par.O + 1.96*se				# Symmetric CIs
    LowerCI <- Par.O - 1.96*se
    
    UpperCI.trans <- Par.O.trans + 1.96*se.trans
    LowerCI.trans <- Par.O.trans - 1.96*se.trans
    # Un-transforming Fishers CIs
    UpperCI.Fisher <- Untrans(UpperCI.trans)
    LowerCI.Fisher <- Untrans(LowerCI.trans)
    
    Res <- as.data.frame(rbind(Res, c(Par.O, se, Par.O.trans, se.trans, a, b, c, d, LowerCI, UpperCI, LowerCI.Fisher, UpperCI.Fisher, PercentLCI, PercentUCI
                                      , Percentile.mean.PAR, Percentile.median.PAR)))
    colnames(Res) <- c("Par", "SE", "Par.trans", "se.trans","a", "b", "c", "d", "LowerCI", "UpperCI", "LowerCI.Fisher", "UpperCI.Fisher","PercentLCI", "PercentUCI"
                       , "Percentile.mean.PAR", "Percentile.median.PAR")
  }
  Res	
}

###===========================================================================###

###############################
# Bayesian Simulation Methods #
###############################

## Bayesian Simulation for Dirchlet(1,1,1,1) prior ##

Dirchsim <- function(p,q,e,n,N){
  Count <<- Count + 1 
  print(Count)
  library(MCMCpack)
  MC <- 10000
  Res <- NULL  
  a <- rep(0,N)
  b <- rep(0,N)
  c <- rep(0,N)
  d <- rep(0,N)
  PAR <- rep(0,N)
  x <- matrix(0, nrow = N, ncol = 4, byrow = TRUE)
  
  for(i in 1:N){
    c.table <- ctable(p,q,e,n)  # Formulation table 
    a[i] <- c.table[1,]
    b[i] <- c.table[3,]
    c[i] <- c.table[2,]
    d[i] <- c.table[4,]
    x[i,] <- c(a[i],b[i],c[i],d[i])
    alpha <- c(1,1,1,1)
    posterior.alpha <- x[i,] + alpha
    
    # MLE of PAR for each table 
    PAR[i] <- PAR2(a[i],b[i],c[i],d[i])
    
    posterior.iterates <- rdirichlet(MC, posterior.alpha)
    p1 <- posterior.iterates[,1]
    p2 <- posterior.iterates[,2]
    p3 <- posterior.iterates[,3]
    p4 <- posterior.iterates[,4]
    PAR.iterates <- ((p1 + p3)/(p1 + p2 + p3 + p4))-(p3/(p3 + p4))
    
    PAR.MLE <- PAR[i]
    PAR.mean <- mean(PAR.iterates)
    PAR.sd <- sd(PAR.iterates)
    PAR.median <- median(PAR.iterates)
    PAR.PI <- quantile(PAR.iterates, c(0.025,0.975))
    
    Res <- as.data.frame(rbind(Res,c(PAR.MLE, PAR.mean, PAR.sd, PAR.median, PAR.PI[1], PAR.PI[2])))
    colnames(Res) <- c("PAR.MLE", "PAR.mean","PAR.sd","PAR.median", "LowerPI", "UpperPI") 
  }
  Res
} 

## Bayesian Simulation for Jeffreys prior ##

Jeffreysim <- function(p,q,e,n,N){
  
  library(MCMCpack)
  MC <- 10000
  Res <- NULL  
  a <- rep(0,N)
  b <- rep(0,N)
  c <- rep(0,N)
  d <- rep(0,N)
  PAR <- rep(0,N)
  PAR.iterates <- rep(0,N)
  x <- matrix(0, nrow = N, ncol = 4, byrow = TRUE)
  
  for(i in 1:N){
    c.table <- ctable(p,q,e,n)  # Formulation table 
    a[i] <- c.table[1,]
    b[i] <- c.table[3,]
    c[i] <- c.table[2,]
    d[i] <- c.table[4,]
    x[i,] <- c(a[i],b[i],c[i],d[i])
    alpha <- c(0.5,0.5,0.5,0.5)
    posterior.alpha <- x[i,] + alpha
    
    # MLE of PAR for each table 
    PAR[i] <- PAR2(a[i],b[i],c[i],d[i])
    posterior.iterates <- rdirichlet(MC, posterior.alpha)
    # This loop detects if c and d are both equal 0 and if so re-draws from the posterior, as 
    # PAR will be indeterminate at this point. 
    repeat {
      defects <- intersect(which(posterior.iterates[,3]==0), which(posterior.iterates[,4]==0))
      if (length(defects) == 0) break;
      for (j in defects) {
        New.iterate <- rdirichlet(1, posterior.alpha)
        posterior.iterates[j,] <- New.iterate
      }
    }
    
    p1 <- posterior.iterates[,1]
    p2 <- posterior.iterates[,2]
    p3 <- posterior.iterates[,3]
    p4 <- posterior.iterates[,4]
    
    PAR.iterates <- ((p1 + p3)/(p1 + p2 + p3 + p4))-(p3/(p3 + p4))
    
    PAR.MLE <- PAR[i]
    PAR.mean <- mean(PAR.iterates)
    PAR.sd <- sd(PAR.iterates)
    PAR.median <- median(PAR.iterates)
    PAR.PI <- quantile(PAR.iterates, c(0.025,0.975))
    
    Res <- as.data.frame(rbind(Res,c(PAR.MLE, PAR.mean, PAR.sd, PAR.median, PAR.PI[1], PAR.PI[2])))
    colnames(Res) <- c("PAR.MLE", "PAR.mean","PAR.sd","PAR.median", "LowerPI", "UpperPI") 
  }
  Res
}  

## Bayesian simulation with improper prior ##

Impropersim <- function(p,q,e,n,N){
  #Count <<- Count + 1 
  #print(Count)
  library(MCMCpack)
  MC <- 10000
  Res <- NULL  
  a <- rep(0,N)
  b <- rep(0,N)
  c <- rep(0,N)
  d <- rep(0,N)
  PAR <- rep(0,N)
  x <- matrix(0, nrow = N, ncol = 4, byrow = TRUE)
  
  for(i in 1:N){
    c.table <- ctable(p,q,e,n)  # Formulation table 
    a[i] <- c.table[1,]
    b[i] <- c.table[3,]
    c[i] <- c.table[2,]
    d[i] <- c.table[4,]
    x[i,] <- c(a[i],b[i],c[i],d[i])
    alpha <- c(0,0,0,0)
    posterior.alpha <- x[i,] + alpha
    
    # MLE of PAR for each table 
    PAR[i] <- PAR2(a[i],b[i],c[i],d[i])
    
    posterior.iterates <- rdirichlet(MC, posterior.alpha)
    p1 <- posterior.iterates[,1]
    p2 <- posterior.iterates[,2]
    p3 <- posterior.iterates[,3]
    p4 <- posterior.iterates[,4]
    
    
    PAR.iterates <- ((p1 + p3)/(p1 + p2 + p3 + p4))-(p3/(p3 + p4))
    
    PAR.MLE <- PAR[i]
    PAR.mean <- mean(PAR.iterates)
    PAR.sd <- sd(PAR.iterates)
    PAR.median <- median(PAR.iterates)
    PAR.PI <- quantile(PAR.iterates, c(0.025,0.975))
    
    Res <- as.data.frame(rbind(Res,c(PAR.MLE, PAR.mean, PAR.sd, PAR.median, PAR.PI[1], PAR.PI[2])))
    colnames(Res) <- c("PAR.MLE", "PAR.mean","PAR.sd","PAR.median", "LowerPI", "UpperPI") 
  }
  Res
}  

###################################
# Flat priors for Bayesian method #
###################################

# flat over [-1,1] #
flat1sim <- function(p,q,e,n,N){
  Count <<- Count + 1 
  print(Count)
  library(MCMCpack)
  MC <- 10000
  Res <- NULL  
  a <- rep(0,N)
  b <- rep(0,N)
  c <- rep(0,N)
  d <- rep(0,N)
  PAR <- rep(0,N)
  PAR.iterates <- rep(0,N)
  x <- matrix(0, nrow = N, ncol = 4, byrow = TRUE)
  
  for(i in 1:N){
    c.table <- ctable(p,q,e,n)  # Formulation table 
    a[i] <- c.table[1,]
    b[i] <- c.table[3,]
    c[i] <- c.table[2,]
    d[i] <- c.table[4,]
    x[i,] <- c(a[i],b[i],c[i],d[i])
    alpha <- c(1,1,0.001,0.001)
    posterior.alpha <- x[i,] + alpha
    
    # MLE of PAR for each table 
    PAR[i] <- PAR2(a[i],b[i],c[i],d[i])
    
    posterior.iterates <- rdirichlet(MC, posterior.alpha)
    # This loop detects if c and d are both equal 0 and if so re-draws from the posterior, as 
    # PAR will be indeterminate at this point. 
    repeat {
      defects <- intersect(which(posterior.iterates[,3]==0), which(posterior.iterates[,4]==0))
      if (length(defects) == 0) break;
      for (j in defects) {
        New.iterate <- rdirichlet(1, posterior.alpha)
        posterior.iterates[j,] <- New.iterate
      }
    }
    
    p1 <- posterior.iterates[,1]
    p2 <- posterior.iterates[,2]
    p3 <- posterior.iterates[,3]
    p4 <- posterior.iterates[,4]
    
    PAR.iterates <- ((p1 + p3)/(p1 + p2 + p3 + p4))-(p3/(p3 + p4))
    
    PAR.MLE <- PAR[i]
    PAR.mean <- mean(PAR.iterates)
    PAR.sd <- sd(PAR.iterates)
    PAR.median <- median(PAR.iterates)
    PAR.PI <- quantile(PAR.iterates, c(0.025,0.975))
    
    Res <- as.data.frame(rbind(Res,c(PAR.MLE, PAR.mean, PAR.sd, PAR.median, PAR.PI[1], PAR.PI[2])))
    colnames(Res) <- c("PAR.MLE", "PAR.mean","PAR.sd","PAR.median", "LowerPI", "UpperPI") 
  }
  Res
}  

# flat over [0,1] # 
flat2sim <- function(p,q,e,n,N){
  Count <<- Count + 1 
  print(Count)
  library(MCMCpack)
  MC <- 10000
  Res <- NULL  
  a <- rep(0,N)
  b <- rep(0,N)
  c <- rep(0,N)
  d <- rep(0,N)
  PAR <- rep(0,N)
  x <- matrix(0, nrow = N, ncol = 4, byrow = TRUE)
  
  for(i in 1:N){
    c.table <- ctable(p,q,e,n)  # Formulation table 
    a[i] <- c.table[1,]
    b[i] <- c.table[3,]
    c[i] <- c.table[2,]
    d[i] <- c.table[4,]
    x[i,] <- c(a[i],b[i],c[i],d[i])
    alpha <- c(1,0.01,0.01,1)
    posterior.alpha <- x[i,] + alpha
    
    # MLE of PAR for each table 
    PAR[i] <- PAR2(a[i],b[i],c[i],d[i])
    
    posterior.iterates <- rdirichlet(MC, posterior.alpha)
    p1 <- posterior.iterates[,1]
    p2 <- posterior.iterates[,2]
    p3 <- posterior.iterates[,3]
    p4 <- posterior.iterates[,4]
    PAR.iterates <- ((p1 + p3)/(p1 + p2 + p3 + p4))-(p3/(p3 + p4))
    
    PAR.MLE <- PAR[i]
    PAR.mean <- mean(PAR.iterates)
    PAR.sd <- sd(PAR.iterates)
    PAR.median <- median(PAR.iterates)
    PAR.PI <- quantile(PAR.iterates, c(0.025,0.975))
    
    Res <- as.data.frame(rbind(Res,c(PAR.MLE, PAR.mean, PAR.sd, PAR.median, PAR.PI[1], PAR.PI[2])))
    colnames(Res) <- c("PAR.MLE", "PAR.mean","PAR.sd","PAR.median", "LowerPI", "UpperPI") 
  }
  Res
}  

# flat over [-1,0] #
flat3sim <- function(p,q,e,n,N){
  #Count <<- Count + 1 
  #print(Count)
  library(MCMCpack)
  MC <- 10000
  Res <- NULL  
  a <- rep(0,N)
  b <- rep(0,N)
  c <- rep(0,N)
  d <- rep(0,N)
  PAR <- rep(0,N)
  x <- matrix(0, nrow = N, ncol = 4, byrow = TRUE)
  
  for(i in 1:N){
    c.table <- table(p,q,e,n)  # Formulation table 
    a[i] <- c.table[1,]
    b[i] <- c.table[3,]
    c[i] <- c.table[2,]
    d[i] <- c.table[4,]
    x[i,] <- c(a[i],b[i],c[i],d[i])
    alpha <- c(1,0.01,0.01,1)
    posterior.alpha <- x[i,] + alpha
    
    # MLE of PAR for each table 
    PAR[i] <- PAR2(a[i],b[i],c[i],d[i])
    
    posterior.iterates <- rdirichlet(MC, posterior.alpha)
    p1 <- posterior.iterates[,1]
    p2 <- posterior.iterates[,2]
    p3 <- posterior.iterates[,3]
    p4 <- posterior.iterates[,4]
    PAR.iterates <- ((p1 + p3)/(p1 + p2 + p3 + p4))-(p3/(p3 + p4))
    
    PAR.MLE <- PAR[i]
    PAR.mean <- mean(PAR.iterates)
    PAR.sd <- sd(PAR.iterates)
    PAR.median <- median(PAR.iterates)
    PAR.PI <- quantile(PAR.iterates, c(0.025,0.975))
    
    Res <- as.data.frame(rbind(Res,c(PAR.MLE, PAR.mean, PAR.sd, PAR.median, PAR.PI[1], PAR.PI[2])))
    colnames(Res) <- c("PAR.MLE", "PAR.mean","PAR.sd","PAR.median", "LowerPI", "UpperPI") 
  }
  Res
}  
