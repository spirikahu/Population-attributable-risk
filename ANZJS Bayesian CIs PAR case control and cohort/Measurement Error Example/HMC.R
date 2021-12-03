#########################################################################################
# Hamiltonian Monte Carlo (HMC)                                                         #
#                                                                                       #
# This program is an adaptation of that which appears in Figure 2 of "MCMC using        #
# Hamiltonian dynamics", Radford M. Neal, 2010 (Handbook of Markov Chain Monte Carlo)   #
# Note that thining and burn in have not been included here and were performed post-hoc #                                                                                                                                                                             #
# when appropriate tuning parameters could be established.                              #                                                                                       #
#                                                                                       #
# The arguments to the HMC function are as follows:                                     #
#                                                                                       #
#   U          A function to evaluate minus the log of the density of the               #
#              distribution to be sampled, plus any constant - ie, the                  #
#              "potential energy".                                                      #
#                                                                                       #
#   grad_U     A function to evaluate the gradient of U.                                #
#                                                                                       #
#   epsilon    The stepsize to use for the leapfrog steps.                              #
#                                                                                       #
#   L          The number of leapfrog steps to do to propose a new state.               #
#                                                                                       #
#   current_q  The current state (position variables only).                             #
#                                                                                       #
#   prior      The beta priors for each parameter                                       #
#                                                                                       #
#   N          Total number of desired iterations                                       #
#                                                                                       #
# Momentum variables are sampled from independent standard normal                       #
# distributions within this function.  The value return is the vector                   #
# of new position variables (equal to current_q if the endpoint of the                  #
# trajectory was rejected).                                                             #
#########################################################################################

###############################################
## Momentum and Gradient functions for the   ##
## leptospirosis example to be used as input ##
## variables to HMC.                         ##
###############################################

## Momentum function ##
U <- function(data, theta, prior){
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
  
  q11 <- Se*p*e + (1-Sp)*q*(1-e)
  q12 <- Se*(1-p)*e + (1-Sp)*(1-q)*(1-e)
  q21 <- Sp*(1-e)*q + (1-Se)*p*e
  q22 <- (1-Se)*(1-p)*e + Sp*(1-q)*(1-e)
  
  # log-likelihood for the posterior distribution we wish to sample #
  loglik = x11*log(q11) +  x12*log(q12) +  x21*log(q21) + x22*log(q22) + 
    (alpha.Se-1)*log(Se) + (beta.Se-1)*log(1-Se) +
    (alpha.Sp-1)*log(Sp) + (beta.Sp-1)*log(1-Sp) +
    (alpha.p-1)*log(p) + (beta.p-1)*log(1-p) + 
    (alpha.q-1)*log(q) + (beta.q-1)*log(1-q) + 
    (alpha.e-1)*log(e) + (beta.e-1)*log(1-e)
  U = -loglik
  return(U)
}

## gradient function ##
grad_U <- function(data, theta, prior){
  ## THIS FUNCTION DIFFERS FROM THAT IN THE MALA ALGORITHM ##
  ## as we are evaluating -log likelihood of the posterior ##
  
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
  dp <- -((e*Se*x11)/q11) + ((Se*x12*e)/q12) - (((1-Se)*e*x21)/q21) + ((e*(1-Se)*x22)/q22) -
    ((alpha.p-1)/p) + ((beta.p-1)/(1-p))       
  
  dq <- -(((1-Sp)*(1-e)*x11)/q11) + (((1-Sp)*(1-e)*x12)/q12) - ((Sp*(1-e)*x21)/q21) +
    ((Sp*(1-e)*x22)/q22) - (alpha.q-1)/q + ((beta.q-1)/(1-q))
  
  de <- -(((Se*p + (Sp-1))*x11)/q11) - (((Se*(1-p) + (1-Sp)*(q-1))*x12)/q12) -
    ((((1-Se)*p - Sp*q)*x21)/q21) - ((((1-Se)*(1-p) + Sp*(q-1))*x22)/q22) -
    ((alpha.e - 1)/e) + ((beta.e -1)/(1-e))
  
  dSe <- -((p*e*x11)/q11) - (((1-p)*e*x12)/q12) + ((p*e*x21)/q21) + ((e*(1-p)*x22)/q22) -
    ((alpha.Se-1)/Se) + ((beta.Se-1)/(1-Se))
  
  dSp <- ((q*(1-e)*x11)/q11) + (((1-q)*(1-e)*x12)/q12) - (((1-e)*q*x21)/q21) - 
    (((1-q)*(1-e)*x22)/q22) - ((alpha.Sp-1)/Sp) + ((beta.Sp -1)/(1-Sp))
  grad <- matrix(c(dp, dq, de, dSe, dSp), nrow=5, ncol=1)
  return(grad)
}

###################
## HMC function  ##
###################
HMC = function (U, grad_U, epsilon, L, N, data, prior){
  ## Defining elements of the prior ## 
  alpha.p <- prior[1,1]; beta.p <- prior[1,2]
  alpha.q <- prior[2,1]; beta.q <- prior[2,2]
  alpha.e <- prior[3,1]; beta.e <- prior[3,2]
  alpha.Se <- prior[4,1]; beta.Se <- prior[4,2]
  alpha.Sp <- prior[5,1]; beta.Sp <- prior[5,2]
  
  ## Initialization ## 
  theta <- matrix(NA, nrow =1, ncol = 5)
  
  repeat {
    Se <- rbeta(n=1, shape1=alpha.Se, shape2=beta.Se)
    Sp <- rbeta(n=1, shape1=alpha.Sp, shape2=beta.Sp)
    
    q11 <- data[1]/sum(data) 
    q12 <- data[2]/sum(data)
    q21 <- data[3]/sum(data)
    q22 <- data[4]/sum(data)
    
    a <- matrix(c(Se,0,(1-Sp),0,0,Se,0,(1-Sp),(1-Se),0,Sp,0,0,(1-Se),0,Sp),nrow=4,ncol=4, byrow=TRUE) 
    
    b <- matrix(c(q11,q12,q21,q22),nrow=4,ncol=1)
    
    soln <- solve(a,b)
    p11 <- soln[1,]
    p12 <- soln[2,]
    p21 <- soln[3,]
    p22 <- soln[4,]
    
    if(length(which(soln >= 0 & soln <= 1)) >= 4){ 
      e <-  p11 + p12 
      p <- p11/e
      q <- p21/(1-e) 
      theta[1,] <- c(p,q,e,Se,Sp)
      break
    }
  }
  current_Q <- theta[1,]
  for(i in 1:N){
    print(i)
    Q = current_Q
    #P = rnorm(length(Q),0,1)  # independent standard normal variates
    set.seed(123)
    P <- matrix(mvrnorm(n=1, mu = rep(0,5), Sigma = diag(5)), ncol=1, nrow=5)
    current_P = P
    
    # Make a half step for momentum at the beginning
    
    P = P - epsilon * grad_U(data,theta=Q,prior) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L){
      # Make a full step for the position
      
      Q = Q + epsilon * P
      
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) P = P - epsilon * grad_U(data,theta=Q,prior)
    }
    
    # Make a half step for momentum at the end.
    
    P = P - epsilon * grad_U(data, theta=Q, prior) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    P = -P
    if(length(which(Q <= 0 | Q >= 1)) > 0){
      theta <- rbind(theta,t(current_Q))  # reject
      current_Q <- current_Q
    } else{
      # Evaluate potential and kinetic energies at start and end of trajectory
      
      current_U = U(data, theta=current_Q, prior)
      current_K = sum(current_P^2) / 2
      proposed_U = U(data, theta=t(Q), prior)
      proposed_K = sum(P^2) / 2
      
      # Accept or reject the state at end of trajectory, returning either
      # the position at the end of the trajectory or the initial position
      
      if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)){
        theta <- rbind(theta, t(Q))  # accept
        current_Q <- Q
      }else{
        theta <- rbind(theta,t(current_Q))  # reject
        current_Q <- current_Q
      }
    }
  }
  theta <- data.frame(theta)
  colnames(theta) <- c("p","q", "e","Se","Sp")
  return(theta)
}

## Running the HMC Function ##
Results <- HMC(U, grad_U, epsilon=0.001, L=50, N=100, data=c(22,25,82,251), 
                        prior=matrix(c(1, 1, 1, 1, 2, 2, 25, 3, 30, 1.5), nrow=5, ncol=2, byrow=T))
