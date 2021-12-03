##################################################################################
# Case-control study applying prior to P(E+)                                     #
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
  e.save <- rbind(e.save, NA)
  
  ## Inital conditions ## 
  u <- rbeta(n=1, shape1=data[1] + 1, shape2=1 + n1 - data[1]); u.save <- rbind(u.save, u)
  v <- rbeta(n=1, shape1=data[2] + 1, shape2=1 + n2 - data[2]); v.save <- rbind(v.save, v)
  
  u.save <- c(u.save, u)
  v.save <- c(v.save, v)
  
  ##########################
  # Gibbs/MH Sampler Start #
  ##########################
  
  for(i in 1:N){ # i loop for total number of iterations
    print(i)
  
    # We want to draw an e for the time step t from a truncated beta(a=1,b=10)  
    # via inverse cdf sampling. To do this we need to first sample d which
    # represents the proposal value for e on the constrained cdf scale.
    # The constraint been that e needs to be between u and v.
    
    # start by determining which is the max and min, u or v? 
    AL <- min(u,v) 
    AU <- max(u,v)
    # calculate the upper/lower bounds for e on the cdf scale
    AL_bounded <- pbeta(AL, shape1 = 1, shape2 = 10)
    AU_bounded <- pbeta(AU, shape1 = 1, shape2 = 10)
    # select a d (that is a proposal e but on the cdf scale)
    d <- runif(1, min=AL_bounded, max=AU_bounded) 
    # now simply take the inverse to get the proposed value of e
    e <- qbeta(d, shape1=1, shape2 = 10)
    e.save <- c(e.save, e)
    
    # sample u an v proposals from beta distributions
    u.proposal <- rbeta(1, shape1=1+data[1], shape2=1+n1-data[1])
    v.proposal <- rbeta(1, shape1=1+data[2], shape2=1+n2-data[2])
    
    # if e is not between the u and v proposals keep sampling u and v until this 
    # condition is met. If the condition is already met then nothing changes 
    while ((u.proposal-e)*(v.proposal-e) >= 0){ 
    u.proposal <- rbeta(1, shape1=1+data[1], shape2=1+n1-data[1])
    v.proposal <- rbeta(1, shape1=1+data[2], shape2=1+n2-data[2])
    }
  
    u <- u.proposal
    v <- v.proposal
    
    u.save <- c(u.save, u)
    v.save <- c(v.save, v)
    w <- (e-v)/(u-v) # P(D+)
    w.save <- c(w.save, w)
    
    # calculation of the PAR 
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

Results <- CCE(data=data, N=2001)
Results2 <- CCE(data=data, N=2001)
Results3 <- CCE(data=data, N=2001)    
    

Results_chain <- Results$Theta
Results_chain$Iteration <- 1:length(Results_chain$PAR)  
# burnin 
burnin_Res1 <- Results_chain %>% filter(Iteration < 1001)
burnin_Res1 %>% summarise(mean_PAR = mean(PAR), sd_PAR = sd(PAR),
                          PAR_CI_LL = quantile(PAR,0.025), 
                          PAR_CI_UL = quantile(PAR,0.975),
                          mean_PAF = mean(PAF), sd_PAF = sd(PAF),
                          PAF_CI_LL = quantile(PAF,0.025), 
                          PAF_CI_UL = quantile(PAF,0.975))

Results2_chain <- Results2$Theta
Results2_chain$Iteration <- 1:length(Results2_chain$PAR)  
burnin_Res2 <- Results2_chain %>% filter(Iteration < 1001)
burnin_Res2 %>% summarise(mean_PAR = mean(PAR), sd_PAR = sd(PAR),
                          PAR_CI_LL = quantile(PAR,0.025), 
                          PAR_CI_UL = quantile(PAR,0.975),
                          mean_PAF = mean(PAF), sd_PAF = sd(PAF),
                          PAF_CI_LL = quantile(PAF,0.025), 
                          PAF_CI_UL = quantile(PAF,0.975))


Results3_chain <- Results3$Theta
Results3_chain$Iteration <- 1:length(Results3_chain$PAR)  
burnin_Res3 <- Results3_chain %>% filter(Iteration < 1001)
burnin_Res3 %>% summarise(mean_PAR = mean(PAR), sd_PAR = sd(PAR),
                          PAR_CI_LL = quantile(PAR,0.025), 
                          PAR_CI_UL = quantile(PAR,0.975),
                          mean_PAF = mean(PAF), sd_PAF = sd(PAF),
                          PAF_CI_LL = quantile(PAF,0.025), 
                          PAF_CI_UL = quantile(PAF,0.975))

# Doing trace plots
Results_chain$PAR_chain2 <- Results2_chain$PAR
Results_chain$PAR_chain3 <- Results3_chain$PAR

Results_chain$PAF_chain2 <- Results2_chain$PAF
Results_chain$PAF_chain3 <- Results3_chain$PAF

# wide to long 
test <- gather(Results_chain, chain, PAR_measurement, c(PAR,PAR_chain2, PAR_chain3, PAF, PAF_chain2, PAF_chain3), factor_key=TRUE)

PAR_res <- test %>% filter(chain == "PAR" | chain=="PAR_chain2" | chain=="PAR_chain3")
PAR_res$chain <- droplevels(PAR_res$chain)

ggplot(aes(x=Iteration, y=PAR_measurement), data=PAR_res) + geom_line() + facet_grid(chain~.)
