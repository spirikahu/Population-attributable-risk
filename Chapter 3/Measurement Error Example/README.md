# Measurement error example

This collection provides code for each of the different MCMC samplers discussed and compared for the cross-sectional study where the uncertainty associated with an imperfect diagnostic test has being considered. 

The code for the Metropolis-Hastings (MH) independence sampler, MH random-walk sampler, Gibbs sampler, Metropolis adjusted Langevin algorithm and Hamiltonian Monte Carlo algorithm are provided in the files CWMH.R, RWMH.R, Gibbs.R, MALA.R and HMC.R respectively. 

The new adapted MH random-walk samplers discussed in this thesis, which take into consideration the shape of the posterior ridge in order to sample more effectively are provided in the files AdMGJTJ.R, AdMHJTDJ.R and AdJTDJprior.R respectively which reprsent the three different forms of Sigma* that have being suggested by equations (3.18-3.20) in the thesis. 

Gustafson (2015) importance sampling approach for non-identifiable models, and by far the superior option when it comes to non-identified models, has been implemented in the file Imp.R. In order to acheive the limiting posterior distribution from this code one can simple set the q vector (described as eta in the thesis) to fixed values based on the data.  
