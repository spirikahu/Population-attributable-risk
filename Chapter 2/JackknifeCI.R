#############################################################
## Estimating a confidence interval for the PAR where      ##
## jackknife re-sampling is used to estimate the variance. ##
## A Fisher-z-transformation is also provided by this      ## 
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
## LowerF & UpperF CI (Fisher transformed)                 ##
#############################################################

#JackCI(a=22, b=25, c=82, d=251)

JackCI <- function(a,b,c,d, conf = 0.95){
  n <- a + b + c + d
  ParO <- PAR2(a=a, b=b, c=c, d=d)      # Par original table (MLE estimate of PAR)
  ParO.trans <- Partrans(ParO)          # Transforming PAR original table (MLE estimate of transformed PAR) 
  Par1 <- PAR2(a=a-1, b=b, c=c, d=d)  	# Par for Jackknife sample 1
  Par2 <- PAR2(a=a, b=b-1, c=c, d=d)		# Par for Jackknife sample 2
  Par3 <- PAR2(a=a, b=b, c=c-1, d=d)		# Par for Jackknife sample 3
  Par4 <- PAR2(a=a, b=b, c=c, d=d-1)		# Par for Jackknife sample 4
  
  Par1.trans <- Partrans(Par1)                      # Transformed Jackknife samples Par
  Par2.trans <- Partrans(Par2)
  Par3.trans <- Partrans(Par3)
  Par4.trans <- Partrans(Par4)
  
  meanPar <- (a*Par1 + b*Par2 + c*Par3 + d*Par4)/n
  var <-(a*((Par1-meanPar)^2)+ b*((Par2-meanPar)^2) + c*((Par3-meanPar)^2) + d*((Par4-meanPar)^2))*(n/(n-1))
  
  meanPar.trans <- (a*Par1.trans + b*Par2.trans + c*Par3.trans + d*Par4.trans)/n
  var.trans <-(a*((Par1.trans-meanPar.trans)^2)+ b*((Par2.trans-meanPar.trans)^2) + c*((Par3.trans-meanPar.trans)^2)
               + d*((Par4.trans-meanPar.trans)^2))*(n/(n-1))
  
  se <- sqrt(var)
  se.trans <- sqrt(var.trans)
  
  JakRes <- cbind(ParO, se, ParO.trans, se.trans)
  colnames(JakRes) <- c("Par", "se", "Par.trans", "se.trans")
  
  multiple <- qnorm(0.5 + conf/2)
  UpperCI <- JakRes[,1] + multiple*(JakRes[,2])
  LowerCI <- JakRes[,1] - multiple*(JakRes[,2])
  
  UpperCI.trans <- JakRes[,3] + multiple*(JakRes[,4])
  LowerCI.trans <- JakRes[,3] - multiple*(JakRes[,4]) 
  
  # Un-transforming Fishers CIs
  UpperCI.Fisher <- Untrans(UpperCI.trans)
  LowerCI.Fisher <- Untrans(LowerCI.trans)
  
  Jackknife <- data.frame(cbind(PAR = ParO, Lower=LowerCI,
                                Upper=UpperCI, LowerF=LowerCI.Fisher, UpperF=UpperCI.Fisher), row.names=""); Jackknife
}
