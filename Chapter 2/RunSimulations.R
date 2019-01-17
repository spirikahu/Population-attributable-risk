###########################################################################
## The following code allows the user to carry out the three simulations ##
## described in Pirikahu(2016) or Chapter 2 of the authors PhD thesis.   ##
## Note that these simulations can take a significant amount of time.    ##
##                                                                       ##
## In order for this code to run the user must first have sourced        ##
## Useful_Functions.R and Simulation_Functions.R                         ##
##                                                                       ##
###########################################################################

#####################
# Full Simulation   #
#####################

Simulation <- function(p, q, e, n, N){
  Count = 0    # Needs to be stated as a Global Variable  
  OverallSim <- list(q=q, e=e, p=p, n=n, a=p*e*n, b=((1-p)*e)*n, c=(q*(1-e))*n, d=((1-q)*(1-e))*n, 
                     TruePar=PAR2(a=p*e, b=(1-p)*e, c=q*(1-e), d=(1-q)*(1-e)),
                     Delta=Deltasim(p,q,e,n,N), Jack=JackSim(p,q,e,n,N), Boot=Bootsim(p, q, e, n, N),
                     Bayes=Dirchsim(p,q,e,n,N), Bayes2=Jeffreysim(p,q,e,n,N), Bayes3=Impropersim(p,q,e,n,N))
} 

#########################
## To run simulation 1 ##
#########################

## The parameter space explored for Simulation 1 ##

p <- c(0.01,0.05,0.1,0.2,0.3,0.35,0.4,0.45,0.5)		
q <- c(0.001,0.05,0.1,0.2,0.3,0.35,0.4,0.45,0.5)	
e <- seq(0.1,0.9,0.1)						

xi <- expand.grid(p,q,e)	# In order to find all combinations of e, p and q
colnames(xi) <- c("p", "q", "e")

## running simulation ##
Results <- mapply(Simulation, p=xi$p, q=xi$q, e=xi$e, n=380, N=10)

#########################
## To run simulation 2 ##
#########################

## The parameter space explored for Simulation 2 ##
a <- 5
c <- 5
e <- seq(0.1,0.9,length.out=100)	

# Calculating starting paramters p and q
startpar <- NULL
vals <- NULL
for(i in 1:length(e)){
  n <- 380
  p <- a/(e[i]*n)
  q <- c/((1-e[i])*n)
  vals <- as.data.frame(cbind(p,q,e,n))
  startpar <- as.data.frame(rbind(vals, startpar))
}

## Running simulations ##
Results2 <- mapply(Simulation, p=startpar$p, q=startpar$q, e=startpar$e, n=380, N=10) # n=380
Results2.1 <- mapply(Simulation, p=startpar$p, q=startpar$q, e=startpar$e, n=1000, N=10) # n=1000

#########################
## To run simulation 3 ##
#########################

## The parameter space explored for Simulation 3 ##
a <- seq(1,25,1) 
c <- seq(1,25,1)
acvalues <- expand.grid(a,c) # get all combinations of a and c
colnames(acvalues) <- c("a", "c")

e <- 0.2 # an area of less coverage seen in simulation 2

# Calculating starting paramters p and q
startpar <- NULL
vals <- NULL
for(i in length(acvalues$a):1){
  n <- 380
  p <- acvalues[i,1]/(e*n)
  q <- acvalues[i,2]/((1-e)*n)
  vals <- as.data.frame(cbind(p,q,e,n))
  startpar <- as.data.frame(rbind(vals, startpar))
}

Results3 <- mapply(Simulation, p=startpar$p, q=startpar$q, e=startpar$e, n=380, N=10)

