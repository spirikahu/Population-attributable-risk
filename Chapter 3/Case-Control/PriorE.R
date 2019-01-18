##################################################################################
# Case-control study applying prior to P(E+)                                     #
#                                                                                #
# This code makes the assumption that a Beta(1,1) prior has been used on the     #
# P(E+|D+) (u) anda P(E+|D-) (v) and a Beta(1,10) prior used on P(E+) (e).       #
# Note that this sampler may take a long time to run if finding a e between u    #                                                                             
# and v is particularly difficult                                                #
#                                                                                #
# This code will be extended at a later date to allow for any beta priors to be  #
# implemented.                                                                   #
#                                                                                #
# Input variables: data - vector of the order (a,b,c,d) where                    #
#       |D+ | D-                                                                 #
#    -------------                                                               #   
#    E+ | a | b                                                                  #    
#    E- | c | d                                                                  #
#                                                                                #
# N - Total number of desired iterations                                         #
##################################################################################

CCE <- function(data, N=100){
  ## Calculating number of cases and controls ##
  n1 <- data[1] + data[3]
  n2 <- data[2] + data[4]  

  ## Initalisation ##
  u.save <- NULL
  v.save <- NULL
  e.save <- NULL 
  w.save <- NULL
  PAR.save <- NULL
  PAF.save <- NULL
  switch.save <- NA
  
  ## Inital conditions ## 
  e.save <- rbind(e.save, NA)
  u <- rbeta(n=1, shape1=data[1] + 1, shape2=1 + n1 - data[1]); u.save <- rbind(u.save, u)
  v <- rbeta(n=1, shape1=data[2] + 1, shape2=1 + n2 - data[2]); v.save <- rbind(v.save, v)
  switch.count = 1
  num.switch = 0
  
  ##########################
  # Gibbs/MH Sampler Start #
  ##########################
  
  for(i in 1:N){ # i loop for total number of iterations
    print(i)
    
    if(switch.count==1){
      e <- rbeta(n=1, 1, 10)
      while (e < v | u < e){
        e <- rbeta(n=1, 2, 10)
      }
      
      ## Gibbs sampling of u ##
      u <- rbeta(n=1, shape1=data[1] + 1, shape2=1 + n1 - data[1])
      
      while (u < e){ # keep sampling until e > u
        u <- rbeta(n=1, shape1=data[1] + 1, shape2=1 + n1 - data[1])
      }
      
      ## Gibbs sampling of v ##
      v <- rbeta(n=1, shape1=data[2] + 1, shape2=1 + n2 - data[2])
      
      while (v > e){ # keep sampling until v > e
        v <- rbeta(n=1, shape1=data[2] + 1, shape2=1 + n2 - data[2])
      }
    } else{
      ## Gibbs sampling of e ##
      e <- rbeta(n=1, 1, 10)
      while (e > v | u > e){ # keep sampling until a value is between u and v
        e <- rbeta(n=1, 2, 10)
      }
      
      ## Gibbs sampling of u ##
      u <- rbeta(n=1, shape1=data[1] + 1, shape2=1 + n1 - data[1])
      
      while (u > e){ # keep sampling until u > e
        u <- rbeta(n=1, shape1=data[1] + 1, shape2=1 + n1 - data[1])
      }
      
      ## Gibbs sampling of v ## 
      v <- rbeta(n=1, shape1=data[2] + 1, shape2=1 + n2 - data[2])
      
      while (v < e){ # keep sampling until e > v
        v <- rbeta(n=1, shape1=data[2] + 1, shape2=1 + n2 - data[2])
      }
    }
    
    ## Switching u and v ##
    u.switch <- v
    v.switch <- u
    
    ## Log likelihood function - product binomial ##
    log.lik <- function(x11, x12, n1, n2, p11, p12){        # remember p11 = u and p12 = v
      log.lik <- x11*log(p11) + (n1 - x11)* log(1-p11) + x12*log(p12) + (n2-x12)*log(1-p12)
      return(log.lik)
    }
    
    ## Log prior function - if the priors are the same then this will cancel in the acceptance probability ##
    log.prior <- function(u,v,e){
      u.prior <- dbeta(u, shape1 =1, shape2=1, log = TRUE)
      v.prior <- dbeta(v, shape1 =1, shape2=1, log = TRUE)
      e.prior <- dbeta(e, shape1 =1, shape2=10, log = TRUE)
      log.prior <- u.prior + v.prior + e.prior
      return(log.prior)
    }
    
    ## Working out the likelihood contribution ##
    switch.log.lik <- log.lik(x11=22,x12=25,n1=n1,n2=n2,p11=u.switch,p12=v.switch)
    log.lik.current <- log.lik(x11=22,x12=25,n1=n1,n2=n2,p11=u,p12=v)
    
    ## Working out the prior contribution ##
    switch.log.prior <- log.prior(u=u.switch, v=v.switch, e=e) 
    log.prior.current <- log.prior(u=u, v=v, e=e)
    
    ## Acceptance probability - only the ratio of the likelihoods, everything else cancels ##
    
    alpha = min(1, exp(switch.log.lik + switch.log.prior - log.lik.current - log.prior.current))
    r <- runif(n=1)
    
    if(r <= alpha){
      u <- u.switch
      v <- v.switch
      switch.count <- -1*switch.count
      num.switch <- num.switch + 1
    }
    switch.save <- c(switch.save, switch.count)
    u.save <- c(u.save, u)
    v.save <- c(v.save, v)
    e.save <- c(e.save, e)
    
    w <- (e-v)/(u-v) 
    w.save <- c(w.save, w)
    
    t1 <- (1-u)*w*(u*w + v*(1-w))
    t2 <- (1-u)*w + (1-v)*(1-w)
    PAR <- (u*w) - (t1/t2)
    
    PAF <- PAR/w
    PAR.save <- c(PAR.save, PAR)
    PAF.save <- c(PAF.save, PAF)
    
  }
  mean <- round(c(mean(PAR.save), mean(PAF.save)), digits = 4)
  sd <- round(c(sd(PAR.save), sd(PAF.save)), digits = 4)
  median <- round(c(median(PAR.save), median(PAF.save)), digits = 4)
  Lower <- round(c(quantile(PAR.save,0.025), quantile(PAF.save,0.025)),digits = 4)
  Upper <-round(c(quantile(PAR.save,0.975), quantile(PAF.save,0.975)),digits = 4)
  Stats <- data.frame(cbind(mean, sd, Lower, median, Upper), row.names=c("PAR", "PAF"))
  Results <- list(Stats = Stats ,Theta = data.frame(cbind(PAR= PAR.save, PAF=PAF.save)))
  return(Results)
}

## Running the CCE function ##
data <- c(22,25,82,251)

Results <- CCE(data=data, N=1000)
