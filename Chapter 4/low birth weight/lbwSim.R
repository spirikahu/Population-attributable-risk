##############################################################
## Low birth weight example simulation code                 ##
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
## Raceprop - vector of true population proportion for each ## 
## race in the order (white, black, other).                 ##
##                                                          ##
## EGR - vector of true population proportions for each     ## 
## level of exposure given race in the order (E+|white,     ##
## E-|white, E+|black, E-|black, E+|other, E-|other).       ##
##                                                          ##
## N - Sample size interested in. Default the total sum for ##
## the data inserted.                                       ##
##                                                          ##
## tunec - tuning parameter when using the Cauchy prior in  ##
## the MCMC sampler.                                        ##
##                                                          ##
## tunef - tuning parameter when using a flat prior in the  ##
## MCMC sampler.                                            ##
##                                                          ## 
## thin - amount of thinning desired. Default = 1.          ##
##                                                          ##
## kick - change the inital condition by addition of        ## 
## specified value. Default = 0.                            ##
##                                                          ##
## Bc - burnin value for sampler using the Cauchy prior     ##
##                                                          ##
## Bf - burnin value for the sampler using the flat prior   ##
##                                                          ##
## Lf - Total number of iterations desired for the MCMC     ##
## sampler with flat prior.                                 ##
##                                                          ##
## Lc - Total number of iterations desired for the MCMC     ##
## sampler with cauchy prior.                               ##
##                                                          ##
## T - total nunmber of iterations desired for the overall  ##
## simulation.                                              ##
##                                                          ##
##############################################################

lbwSim <- function(formula, data, Raceprop, EGR, N=sum(data$N),
                   tunec=1, tunef=1, Bf=0, Bc=0, 
                     thin=1, kick=0, Lf=1000, Lc=1000, T=100){
  set.seed(2908)
  # MCMC function 
  source("logitSampler.R")
  
  # Initialisation 
  beta <- NULL
  PAR.fixed.f <- NULL
  PAR.fixed.c <- NULL
  
  PAR.variable.f <- NULL
  PAR.variable.c <- NULL
  
  PAF.fixed.f <- NULL
  PAF.fixed.c <- NULL
  
  PAF.variable.f <- NULL
  PAF.variable.c <- NULL
  
  for(i in 1:T){
    print(i)
    # Making a new data set which is identical to the original and then re-write over N 
    # from rmulti, then re-write over low from rbinom
    PropW <- Raceprop[1]
    PropB <- Raceprop[2]
    PropO <- Raceprop[3]
    
    PDWE <- EGR[1]
    PDWNE <- EGR[2]
    PDBE <- EGR[3]
    PDBNE <- EGR[4]
    PDOE <- EGR[5]
    PDONE <- EGR[6]
    
    probs <- c(PropW*PDWE, PropW*PDWNE, PropB*PDBE, PropB*PDBNE, PropO*PDOE, PropO*PDONE)
    
    new.N <-rmultinom(n=1, size=N, prob=probs)
    
    model <- glm(formula = formula, data=data, family=binomial, x=T)
    initial <- model$coef   # Pulling out model coefficents
    beta <- rbind(beta,initial)
    
    X <- model$x  # Model matrix which gives every possible combination of betas
    
    beta.current <- beta[1,] 
    p.logit <- X %*% beta.current  # Gives me logit(p)  
    p.current <- 1/(1+exp(-p.logit))   # Back transforming logit
    
    new.low <- rbinom(n=length(p.current), size=new.N, prob=p.current)
    
    data.new <- data
    data.new$N <- new.N
    data.new$low <- new.low
    print(data.new)
    
    ## Doing the MCMC algorithm ##
    flat <- logitfun(formula = formula, data = data.new, prior="flat", tune = tunef, B = Bf,
                     thin = thin, kick = kick ,L = Lf)
    cauchy <- logitfun(formula = formula, data = data.new, prior="cauchy", tune = tunec, B = Bc,
                       thin = thin, kick = kick ,L = Lc)
    
    ## Extracting stats and CIs for PAR/PAF 
    
    PAR.fixed.f <- rbind(PAR.fixed.f, flat$stats[5,])
    PAR.fixed.c <- rbind(PAR.fixed.c, cauchy$stats[5,])
    
    PAR.variable.f <- rbind(PAR.variable.f, flat$stats[6,])
    PAR.variable.c <- rbind(PAR.variable.c, cauchy$stats[6,])
    
    PAF.fixed.f <- rbind(PAF.fixed.f, flat$stats[7,])
    PAF.fixed.c <- rbind(PAF.fixed.c, cauchy$stats[7,])
    
    PAF.variable.f <- rbind(PAF.variable.f, flat$stats[8,])
    PAF.variable.c <- rbind(PAF.variable.c, cauchy$stats[8,])
  }
  
  Res <- list(PAR.Ffixed=PAR.fixed.f, PAR.Cfixed=PAR.fixed.c, 
              PAR.Fvariable=PAR.variable.f, PAR.Cvariable=PAR.variable.c,
              PAF.Ffixed=PAF.fixed.f, PAF.cfixed=PAF.fixed.c,
              PAF.Fvariable=PAF.variable.f, PAF.Cvariable=PAF.variable.c)
  return(Res)
}

## Running lbwSim function ##
lbw <- read.csv("lbwdata.csv")

Results <- lbwSim(formula = cbind(low,N-low)~ smoke + race, data, Raceprop=c(0.5,0.2,0.3),
                  EGR=c(0.6,0.4,0.4,0.6,0.2,0.8), N=192, tunec=1.2, tunef=1.5, Bf=8000, Bc=6000, 
                  thin=1, kick=0, Lf=18000, Lc=16000, T=1000)

## Calculating True PAR/PAF for comparison ##
model <- glm(formula = cbind(low,N-low)~ smoke + race, data = lbw, family = binomial(link=logit),x=T)
coefs <- model$coef   # Pulling out model coefficents

X <- model$x  # Model matrix which gives every possible combination of betas

p.logit <- X %*% coefs  # Gives me logit(p)  
p <- data.frame(t(1/(1+exp(-p.logit))))   # Back transforming logit
colnames(p) <- c("pW", "qW", "pB", "qB", "pO", "qO")

Raceprop=c(0.5,0.5,0.2,0.2,0.3,0.3)
e=data.frame(t(c(0.6,0.4,0.4,0.6,0.2,0.8)))
colnames(e) <-  c("e11", "e21", "e12", "e22", "e13", "e23")
TruePAR <-  Raceprop[1]*e$e11*(p$pW - p$qW) + Raceprop[3]*e$e12*(p$pB - p$qB) + Raceprop[5]*e$e13*(p$pO - p$qO)
TruePAF <- TruePAR/rowSums(Raceprop*e*p)
