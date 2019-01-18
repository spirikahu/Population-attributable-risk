##################################################################################
# Case-control study applying prior to P(D+)                                     #
#                                                                                #
# Given that the posterior can be derived analytically for this case we          #
# can sample directly from the posterior distribution and estimate the           #
# PAR and its 95% credible interval. This code is based on the assumption        #
# that a beta distribution is adopted for the priors on u [phi_1 = P(E+|D+)],    #
# v [phi_2 = P(E+|D-)] and w [phi_3 = P(D+)].                                    #
#                                                                                #
##################################################################################

# data #
data <- c(22,25,82,251)
n1 <- data[1] + data[3]
n2 <- data[2] + data[4]  
# Number of replicates #
R <- 10000

# Setting a prior on w=P(D+) #
w <- rbeta(n=R, shape1= 1, shape2=100)  

# Drawing from posteriors for u=P(E+|D+) and v=P(E+|D-) #
u <- rbeta(n=R, shape1=data[1] + 1, shape2=1 + n1 - data[1])
v <- rbeta(n=R, shape1=data[2] + 1, shape2=1 + n2 - data[2])
e <- u*w + v*(1-w)

PAR <- u*w - (((1-u)*w*(u*w+v*(1-w)))/((1-u)*w + (1-v)*(1-w)))
PAF <- PAR/w

means <- c(mean(PAR), mean(PAF))
CIPAR <- quantile(PAR, c(0.5,0.025,0.975))
CIPAF <- quantile(PAF, c(0.5,0.025,0.975))
