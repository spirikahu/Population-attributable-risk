##############################################################
##  Bayesian Logistic Regression MCMC Sampler               ##
##                                                          ##
## Where the joint distribution of (E,x) is described by a  ##
## Dirichlet distribution and P(D|E) estimated using a      ##
## logistic regression model. This function can be run      ##
## independently but this function is also incorporated in  ##
## simulation code.                                         ##
##                                                          ##
## Input variables:                                         ##
##                                                          ##
## formula - a formula for the glm model in the appropriate ##
## format which the glm function takes. Note that this code ##
## does not currently account for interactions, but maybe   ##
## extended in future.                                      ##
##                                                          ##
## data - a csv file where the data is in the form of a     ##
## 2x2xk table. This will be extended in future to account  ##
## for long format data.                                    ##
##                                                          ##
## prior - specify either "flat" or "cauchy".               ##
## Default = "flat"                                         ##        
##                                                          ##
## tune - scalar tuning parameter                           ##
##                                                          ##
## B - burnin value , default = 0                           ##
##                                                          ##
## thin - amount of thining to apply                        ##
##                                                          ##
## kick - change the inital condition by addition of        ## 
## specified value. Default = 0.                            ##
##                                                          ##
## L - Total number of iterations desired. Default = 1000.  ##
##                                                          ##
##############################################################

logitfun <- function(formula , data, prior="flat", tune, B=0, thin=1, kick = 0, L=1000){
  ## Initialization ##
  beta <- NULL
  p <- NULL
  q <- NULL
  e <- NULL
  accept.rate <- 0
  
  ## Packages ##
  require(MCMCpack)
  
  ## Internal Functions ##
  
  log.lik.data <- function(success,N,p){
    # This if statement deals to the situation where INFs may occur due
    # to probabilities being 0 or 1.
    
    if(sum(is.infinite(log(p))) > 0 | sum(is.infinite(log(1-p))) > 0){
      t1 <- log(p)
      t1[is.infinite(t1)] <- 0 
      
      t2 <- log(1-p)
      t2[is.infinite(t2)] <- 0 
      
      log.lik <- sum(success * t1 + (N-success) * t2) 
      
    } else{
      log.lik <- sum(success * log(p) + (N-success) * log(1-p)) 
    }
    return(log.lik)
  } 
  
  # when there are 0 counts in the table add a fraction of a case
  data.new <- data
  data.new[which(data.new$low == 0),4] <- 0.1
  
  #################################
  ## Assigning Inital Conditions ##
  #################################
  
  ## beta's are assigned from the model on the original data #E
  model <- glm(formula = formula, data=data.new, family=binomial, x=T)
  summary(model)
  initial <- model$coef + kick     # Pulling out model coefficents
  beta <- rbind(beta,initial)
  Cov <- vcov(model)
  
  beta.current <- beta[1,] 
  X <- model$x  # Model matrix which gives every possible combination of betas
  
  p.logit <- X %*% beta.current      # Gives logit(p)
  p.current <- 1/(1+exp(-p.logit))   # Back transforming logit
  
  ## Calculating q P(D+|E-) ##
  X0 <- X
  X0[,2] <- 0
  q.logit <- X0 %*% beta.current
  q.current <- 1/(1+exp(-q.logit)) 
  
  p <- rbind(p,t(p.current))
  q <- rbind(q,t(q.current))
  
  ## Current time steps contributions to alpha ##
  
  # Taking into consideration two different priors here the flat which is equivalent
  # to assigning 1 and the Cauchy prior 
  if(prior=="flat"){
    log.pr.current <- 1 
  } else {
    log.pr.current <- sum(dcauchy(beta.current,location=0, scale=2.5, log=T))
  } 
  
  log.lik.current <- log.lik.data(success=data$low,N=data$N,p=p.current)
  
  e <- rdirichlet(n=1, alpha=data$N + 1)
  
  ## Eigenvalue decomposition ##
  ES<-eigen(Cov,symmetric=T)
  M<-ES$vec
  Sigmasq<-M %*% diag(sqrt(ES$val)) %*% t(M) # sqrt covariance matrix used as sigma in random walk
  
  ########################
  # MCMC algorithm start #
  ########################
  
  for(i in 1:L){
     print(i)
    
    #########################
    #     Gibbs             #
    #########################
    
    e.current <- rdirichlet(n=1, alpha=data$N + 1) # joint distribution (E,x)
    
    #########################
    #   Metropolis-Hastings #
    #########################
    
    ## Block updating approach ##
    beta.cand <- beta.current + tune*Sigmasq %*% rnorm(4)  
    
    p.cand.logit <- X %*% beta.cand
    p.cand <- 1/(1+exp(-p.cand.logit)) 
    
    q.cand.logit <- X0 %*% beta.cand
    q.cand <- 1/(1+exp(-q.cand.logit)) 
    
    ## Candidate contributions to alpha ##
    if(prior=="flat"){
      log.pr.cand <- 1 
    } else {
      log.pr.cand <- sum(dcauchy(beta.cand,location=0, scale=2.5, log=T))
    }
    
    log.lik.cand <- log.lik.data(success=data$low,N=data$N,p=p.cand)
    ## Acceptance Probability ##
    alpha <- min(1,exp(log.pr.cand + log.lik.cand - log.pr.current -log.lik.current))
    
    u <- runif(1)
    
    if(u <= alpha){
      beta.current <- beta.cand
      p.current <- p.cand
      q.current <- q.cand
      log.lik.current <- log.lik.cand
      log.pr.current <- log.pr.cand
      
      accept.rate <- accept.rate + 1
    }
    
    ## Thinning and Burnin Statement ##
    if(i > B && i %% thin == 0){
      beta <- rbind(beta, t(beta.current))
      p <- rbind(p,t(p.current))
      q <- rbind(q,t(q.current))
      e <- rbind(e,e.current)
    }
    
  }
  beta <- data.frame(beta, row.names=NULL)
  p <- data.frame(p, row.names=NULL)
  q <- data.frame(q, row.names = NULL)
  e <- data.frame(e, row.names=NULL)
  
  prop.race <- tapply(data$N, data$race, sum)/sum(data$N)
  smokers <- data[data$smoke =="y",]
  e.fixed <- as.numeric(smokers$N)/tapply(data$N, data$race, sum)    # This is a vector only of the exposure rates E+
  sums <- tapply(data$N, data$race, sum)                      
  ne.fixed <- as.numeric(data$N)/sums[data$race]                     # This is the vector of all exposures E+ and E-
  
  # Calculating PAR/PAF where fixed estimates are used for joing distribution (E,x)
  PAR.fixed <- rowSums((as.matrix((p-q)) %*% diag(e.fixed[data$race])) %*% diag(prop.race[data$race]))
  PAF.fixed <- PAR.fixed/rowSums(as.matrix(p)%*%diag(ne.fixed)%*%diag(prop.race[data$race]))
  
  # Allowing the joint distribution (E,x) to be represented by a Dirchlet distribution 
  PAR.variable <- rowSums(e*(p-q))
  PAF.variable <- PAR.variable/rowSums(e*p)
  
  all <- cbind(beta,PAR.fixed, PAR.variable, PAF.fixed, PAF.variable)
  mean <- colMeans(all)
  median <- apply(all, 2, function(x) median(x))
  sd <- apply(all,2,function(x) sd(x))
  se <- apply(all,2,function(x) sd(x)/sqrt(length(data[,1])))
  Lower <- apply(all,2,function(x) quantile(x,0.025)[[1]])
  Upper <- apply(all,2,function(x) quantile(x,0.975)[[1]])
  Int.Width <- Upper - Lower
  
  stats <- data.frame(cbind(mean, median, sd, se, Lower, Upper, Int.Width))
  
  return(list(beta=beta, stats=stats, accept.rate = accept.rate/L))
}

# Running the logitfun function ##
lbw <- read.csv("lbwdata.csv")
Results <- logitfun(formula = cbind(low,N-low)~ smoke + race, 
                 data = lbw, prior="flat", tune = 1, B = 0,
                 thin = 1, kick = 0 ,L = 10)
