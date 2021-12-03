#############################################################
## Estimating a confidence interval for the PAR where      ##
## bootstrap re-sampling is used to estimate the variance. ##
## The percentile bootstrap CI is also provided by this    ## 
## function. Note that before this function will run the   ##
## user must first source Useful_Functions.R               ##
##                                                         ##  
## Function inputs: a, b, c and d which represent observed ##
## values in the following 2 x 2 table                     ##
##       |D+ | D-                                          ##
##    -------------                                        ##
##    E+ | a | b                                           ##
##    E- | c | d                                           ##
##                                                         ##  
## Function output: PAR, Lower & Upper CI (untransformed), ##
## LowerP & UpperP CI (Percentile bootstrap)               ##
#############################################################

# BootCI(a=22, b=25, c=82, d=251)

BootCI <- function(a,b,c,d, conf=0.95){
  n <- a + b + c + d  
  K <- 1000
  Par <- rep(0,K)  	# Initalising PAR
  Par.O <- PAR2(a,b,c,d)    
  prob2 <- c(a/n, b/n, c/n, d/n)
  # Repeat function used to avoid 0's
  for(j in 1:K){  			# K = Number of Bootstrap Samples
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
  Percentile.mean.PAR <- mean(Par)      # point estimate of the mean from % Bootstrap distribution of PAR
  Percentile.median.PAR <- median(Par)  # point estimate of the median from % Bootstrap distribution of PAR 
  
  PercentileCIs <-  quantile(Par, c(0.5,0.025,0.975))
  PercentLCI <- PercentileCIs[2] 			# Non-symmetric CIs
  PercentUCI <-  PercentileCIs[3] 
  
  var <- sum((Par-mean(Par))^2)/(K-1)
  se <- sqrt(var)
  
  multiple <- qnorm(0.5 + conf/2)
  UpperCI <- Par.O + multiple*se				# Symmetric CIs
  LowerCI <- Par.O - multiple*se
  
  Res <- data.frame(cbind(PAR=Par.O, Lower=LowerCI, Upper=UpperCI, LowerP=PercentLCI, UpperP=PercentUCI), row.names="");Res
}