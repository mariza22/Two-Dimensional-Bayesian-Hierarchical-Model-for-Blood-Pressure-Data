deviance <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\deviance.csv",header=TRUE)[,-1]
deviance <- as.numeric(unlist( c(deviance[1,],deviance[2,],deviance[3,]))) [3001:6001]
theta    <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\theta.csv",header=TRUE)[3001:6001,-1] 
theta_c  <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\theta_c.csv",header=TRUE)[3001:6001,-1]
theta_r  <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\theta_r.csv",header=TRUE)[3001:6001,-1]
theta_s  <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\theta_s.csv",header=TRUE)[3001:6001,-1]
theta_g  <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\theta_g.csv",header=TRUE)[3001:6001,-1]
phi_natl <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\phi_natl.csv",header=TRUE)[3001:6001,-1]
phi_subn <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\phi_subn.csv",header=TRUE)[3001:6001,-1]
phi_comm <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\phi_comm.csv",header=TRUE)[3001:6001,-1]
tau      <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\tau.csv",header=TRUE)[3001:6001,-1]
gamma    <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\gamma.csv",header=TRUE)[3001:6001,-1]
u        <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\u.csv",header=TRUE)[3001:6001,-1]
v        <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\v.csv",header=TRUE)[3001:6001,-1]
sv       <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\sv.csv",header=TRUE)[3001:6001,-1]
w        <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\w.csv",header=TRUE)[3001:6001,-1]
eta_c    <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\eta_c.csv",header=TRUE)[3001:6001,-1]
eta_r    <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\eta_r.csv",header=TRUE)[3001:6001,-1]
eta_s    <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\eta_s.csv",header=TRUE)[3001:6001,-1]
phi_c    <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\phi_c.csv",header=TRUE)[3001:6001,-1]
phi_r    <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\phi_r.csv",header=TRUE)[3001:6001,-1]
phi_s    <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\phi_s.csv",header=TRUE)[3001:6001,-1]
sigma2   <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Results_19-03-22\\sigma2.csv",header=TRUE)[3001:6001,-1]

#########    MCMC.list  ################################################
library(coda)
library("RColorBrewer")
pdf(paste("MCMC.list.pdf", sep=''))
mix          <- list()
title1       <- c(rep("intercept",3),rep("slope",3))
title2       <- rep(c("DBP","SBP","INTER"),2)

#deviance
for ( i in 1:3) mix[[i]] <- mcmc(chains1[[i]]$deviance)
MH.list                  <- mcmc.list(mix)     # add all the chains together
print(gelman.diag(MH.list))      # scale reduction factor diagnostic: it should be close to1.1.05 to 1.1 
gelman.plot(MH.list)             #
#summary(MH.list)
plot(MH.list,col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
#theta
for (j in 1:6){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$theta[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  print(gelman.diag(MH.list))      # scale reduction factor diagnostic: it should be close to1.1.05 to 1.1 
  gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste(title1[j],title2[j]),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
  #traceplot and marginal density plot. Basically, it is the smoothened histogram 
  #of the values in the trace-plot, i.e the distribution of the values of the 
  #parameter in the chain
  #2. Autocorrelation function
  #high correlation is an indication of low mixing.
  autocorr.diag(MH.list)
  autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
  #2.b cross correlation for MCMC output
  crosscorr(MH.list)
  crosscorr.plot(MH.list,col = topo.colors(10))
  #calculates croos correlations between variables in MCMC output, all chainsin MH.list
  #are combined before calculating correlation
  #3 Effective Sample Size
  effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
  lapply(MH.list,effectiveSize)#get the size for each chain individually 
  #the effective sample size indicates how many estimates from the iterations you run are 
  #actually considered good for your mean. For example you might have 1000 iteration and 
  #only 4 of them are depicted the true value.  
  #In STAN they use a fraction between the neff which is the total of the samples close to truth 
  #divided by the total sample size(number of iteration), N. every time we want sth 
  #more than Neff/N > 0.5, around half of them.
  
  #4 Geweke's convergence diagnostic
  geweke.diag(MH.list)
  geweke.plot(MH.list)
}
##theta_c
for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$theta_c)
MH.list    <- mcmc.list(mix)     # add all the chains together
#1 scale reduction factor
print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
gelman.plot(MH.list)             #
#summary(MH.list)
plot(MH.list,col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
##theta_r 
for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$theta_r)
MH.list    <- mcmc.list(mix)     # add all the chains together
#1 scale reduction factor
print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
gelman.plot(MH.list)             #
#summary(MH.list)
plot(MH.list,col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
##theta_s 
for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$theta_r)
MH.list    <- mcmc.list(mix)     # add all the chains together
#1 scale reduction factor
print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
gelman.plot(MH.list)             #
#summary(MH.list)
plot(MH.list,col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
##theta_g
for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$theta_r)
MH.list    <- mcmc.list(mix)     # add all the chains together
#1 scale reduction factor
print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
gelman.plot(MH.list)             #
#summary(MH.list)
plot(MH.list,col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
#phi_natl
for ( j in 1:3){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$phi_natl[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  #1 scale reduction factor
  print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
  gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste(title2[j]),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
#phi_subn
for ( j in 1:3){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$phi_subn[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  #1 scale reduction factor
  print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
  gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste("phi_subn",title2[j]),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
#phi_comm
for ( j in 1:3){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$phi_comm[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  #1 scale reduction factor
  print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
  gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste("phi_comm",title2[j]),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
#tau
for ( j in 1:3){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$tau[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  #1 scale reduction factor
  print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
  gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste("tau",j),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
#gamma
for ( j in 1:15){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$gamma[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  #1 scale reduction factor
  print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
  gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste("gamma",j),col=brewer.pal(n = 10, name = "RdBu"))#Paired"RdBu","Set3" traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
#u
for ( j in 1:250){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$u[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  #1 scale reduction factor
  print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
  #gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste(j),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
#v
for ( j in 1:25){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$v[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  #1 scale reduction factor
  print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
  #gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste(j),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
#sv
for ( j in 1:225){
  for ( i in 1:3) mix[[i]]  <- mcmc(chains1[[i]]$sv[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  #1 scale reduction factor
  print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
  #gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste(j),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
#w
for ( j in 1:25){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$w[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  #1 scale reduction factor
  print(gelman.diag(MH.list))      # diagnostic: it should be close to1.1.05 to 1.1 
  #gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste(j),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
###eta_c
for ( j in 1:3){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$eta_s[,j])#eta_s, eta_c
  MH.list    <- mcmc.list(mix)     # add all the chains together
  print(gelman.diag(MH.list))      #1 scale reduction factor diagnostic: it should be close to1.1.05 to 1.1 
  #gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste(j),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
###phi_s/r/c
for ( j in 1:3){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$phi_r[,j])#phi_s,phi_r
  MH.list    <- mcmc.list(mix)     # add all the chains together
  print(gelman.diag(MH.list))      #1 scale reduction factor diagnostic: it should be close to1.1.05 to 1.1 
  gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste(j),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
###sigma2
for ( j in 1:5){
  for ( i in 1:3) mix[[i]]   <- mcmc(chains1[[i]]$sigma2[,j])
  MH.list    <- mcmc.list(mix)     # add all the chains together
  print(gelman.diag(MH.list))      #1 scale reduction factor diagnostic: it should be close to1.1.05 to 1.1 
  #gelman.plot(MH.list)             #
  #summary(MH.list)
  plot(MH.list,main=paste(j),col=brewer.pal(n = 10, name = "RdBu"))#traceplots, density plots
}
#2. Autocorrelation function
autocorr.diag(MH.list)
autocorr.plot(MH.list)# or autocorr.plot(mix[[1]]) seperately for each chain
#2.b cross correlation for MCMC output
crosscorr(MH.list)
crosscorr.plot(MH.list,col = topo.colors(10))
#3 Effective Sample Size
effectiveSize(MH.list)  #for a mcmc.list object the effective sizes are summed across chains. 
lapply(MH.list,effectiveSize)#get the size for each chain individually 
#4 Geweke's convergence diagnostic
geweke.diag(MH.list)
geweke.plot(MH.list)
######## Combine Chains
sex 			<- "female"
variable 		<- "BP"
##### Combine chains ##################################################
##### Run 1 ###########################################################
burnt             <- seq(1, 2001)
sex 	        		<- "female"
variable 		      <- "bp"
u.combined       	<- chains1[[1]]$u[burnt,]
v.combined       	<- chains1[[1]]$v[burnt,]
sv.combined      	<- chains1[[1]]$sv[burnt,]
w.combined       	<- chains1[[1]]$w[burnt,]
theta.combined   	<- chains1[[1]]$theta[burnt,]
phi_natl.combined <- chains1[[1]]$phi_natl[burnt,]
phi_subn.combined <- chains1[[1]]$phi_subn[burnt,]
phi_comm.combined <- chains1[[1]]$phi_comm[burnt,]
tau.combined 		  <- chains1[[1]]$tau[burnt,]
theta_g.combined 	<- chains1[[1]]$theta_g[burnt]
theta_s.combined 	<- chains1[[1]]$theta_s[burnt]
theta_r.combined 	<- chains1[[1]]$theta_r[burnt]
theta_c.combined 	<- chains1[[1]]$theta_c[burnt]
phi_c.combined   	<- chains1[[1]]$phi_c[burnt,]
phi_r.combined   	<- chains1[[1]]$phi_r[burnt,]
phi_s.combined   	<- chains1[[1]]$phi_s[burnt,]
eta_c.combined   	<- chains1[[1]]$eta_c[burnt,]
eta_r.combined   	<- chains1[[1]]$eta_r[burnt,]
eta_s.combined   	<- chains1[[1]]$eta_s[burnt,]
gamma.combined   	<- chains1[[1]]$gamma[burnt,]
sigma2.combined  	<- chains1[[1]]$sigma2[burnt,]
deviance.combined	<- chains1[[1]]$deviance[burnt]
##### Runs 2 to 10 ####################################################
for(RunNum in 2:3){#34
  print(RunNum)   #seq(1, 4999, length=250)
  u.combined      	 <- rbind(u.combined, chains1[[RunNum]]$u[burnt,])
  v.combined       	 <- rbind(v.combined, chains1[[RunNum]]$v[burnt,])
  sv.combined      	 <- rbind(sv.combined, chains1[[RunNum]]$sv[burnt,])
  w.combined       	 <- rbind(w.combined, chains1[[RunNum]]$w[burnt,])
  theta.combined   	 <- rbind(theta.combined, chains1[[RunNum]]$theta[burnt,])
  phi_natl.combined  <- rbind(phi_natl.combined, chains1[[RunNum]]$phi_natl[burnt,])
  phi_subn.combined  <- rbind(phi_subn.combined, chains1[[RunNum]]$phi_subn[burnt,])
  phi_comm.combined  <- rbind(phi_comm.combined, chains1[[RunNum]]$phi_comm[burnt,])
  tau.combined 		   <- rbind(tau.combined, chains1[[RunNum]]$tau[burnt,])
  theta_g.combined 	 <- c(theta_g.combined, chains1[[RunNum]]$theta_g[burnt])
  theta_s.combined 	 <- c(theta_s.combined, chains1[[RunNum]]$theta_s[burnt])
  theta_r.combined 	 <- c(theta_r.combined, chains1[[RunNum]]$theta_r[burnt])
  theta_c.combined 	 <- c(theta_c.combined, chains1[[RunNum]]$theta_c[burnt])
  phi_c.combined   	 <- rbind(phi_c.combined, chains1[[RunNum]]$phi_c[burnt,])
  phi_r.combined   	 <- rbind(phi_r.combined, chains1[[RunNum]]$phi_r[burnt,])
  phi_s.combined   	 <- rbind(phi_s.combined, chains1[[RunNum]]$phi_s[burnt,])
  eta_c.combined   	 <- rbind(eta_c.combined, chains1[[RunNum]]$eta_c[burnt,])
  eta_r.combined   	 <- rbind(eta_r.combined, chains1[[RunNum]]$eta_r[burnt,])
  eta_s.combined   	 <- rbind(eta_s.combined, chains1[[RunNum]]$eta_s[burnt,])
  gamma.combined   	 <- rbind(gamma.combined, chains1[[RunNum]]$gamma[burnt,])
  sigma2.combined  	 <- rbind(sigma2.combined,chains1[[RunNum]]$sigma2[burnt,])
  deviance.combined	 <- rbind(deviance.combined, chains1[[RunNum]]$deviance[burnt])
}
##### Save chains #####################################################
u  		               <- u.combined;v  		                   <- v.combined
sv 		               <- sv.combined;w  	               	     <- w.combined
theta                <- theta.combined;phi_natl 	           <- phi_natl.combined
phi_subn   	         <- phi_subn.combined;phi_comm   	       <- phi_comm.combined
tau 		             <- tau.combined;theta_g 	               <- theta_g.combined
theta_s      	       <- theta_s.combined;theta_r             <- theta_r.combined
theta_c        	     <- theta_c.combined;phi_c               <- phi_c.combined
phi_r   	           <- phi_r.combined;phi_s   	             <- phi_s.combined
eta_c   	           <- eta_c.combined;eta_r   	             <- eta_r.combined
eta_s   	           <- eta_s.combined;gamma   	             <- gamma.combined
sigma2  	           <- sigma2.combined;deviance	           <- deviance.combined

upto                 <- 2001*3
burnt                <- c(1:upto)
filename <- paste("CV1_",sex,"_",variable,"_Combined", sep='')
rm(u.combined);rm(v.combined);rm(sv.combined);rm(w.combined)
rm(theta.combined);rm(gamma.combined);rm(sigma2.combined)
#save.image(paste(filename,".RData",sep=''));attach(covar);attach(subset)
p=42; multipleRegionsInSregion=19;J=61;N=245;L=9
##### Generate trace plots ############################################
pdf(paste(filename, "bivar_traceplots",".pdf",sep=''),width=10, height=7.5)
par(mfrow=c(3,2), mar=c(0.5,4,0.5,0))
plot(phi_s[burnt,1], type='l', xlab='', ylab='phi_s',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_s[burnt,2], type='l', xlab='', ylab='phi_s',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_s[burnt,3], type='l', xlab='', ylab='phi_s',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(eta_s[burnt,1], type='l', xlab='', ylab='eta_s',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(eta_s[burnt,2], type='l', xlab='', ylab='eta_s',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(eta_s[burnt,3], type='l', xlab='', ylab='eta_s',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_r[burnt,1], type='l', xlab='', ylab='phi_r',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_r[burnt,2], type='l', xlab='', ylab='phi_r',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_r[burnt,3], type='l', xlab='', ylab='phi_r',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(eta_r[burnt,1], type='l', xlab='', ylab='eta_r',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(eta_r[burnt,2], type='l', xlab='', ylab='eta_r',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(eta_r[burnt,3], type='l', xlab='', ylab='eta_r',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_c[burnt,1], type='l', xlab='', ylab='eta_r',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")	
plot(phi_c[burnt,2], type='l', xlab='', ylab='eta_r',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")	
plot(phi_c[burnt,3], type='l', xlab='', ylab='eta_r',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(eta_c[burnt,1], type='l', xlab='', ylab='eta_c',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(eta_c[burnt,2], type='l', xlab='', ylab='eta_c',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(eta_c[burnt,3], type='l', xlab='', ylab='eta_c',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_natl[burnt,1], type='l', xlab='', ylab='phi_natl',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_natl[burnt,2], type='l', xlab='', ylab='phi_natl',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_natl[burnt,3], type='l', xlab='', ylab='phi_natl',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_subn[burnt,1], type='l', xlab='', ylab='phi_subn',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_subn[burnt,2], type='l', xlab='', ylab='phi_subn',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_subn[burnt,3], type='l', xlab='', ylab='phi_subn',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_comm[burnt,1], type='l', xlab='', ylab='phi_comm',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_comm[burnt,2], type='l', xlab='', ylab='phi_comm',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(phi_comm[burnt,3], type='l', xlab='', ylab='phi_comm',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(tau[burnt,1], type='l', xlab='', ylab='tau',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(tau[burnt,2], type='l', xlab='', ylab='tau',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(tau[burnt,3], type='l', xlab='', ylab='tau',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(theta[burnt,1], type='l', xlab='', ylab='a_g',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(theta[burnt,2], type='l', xlab='', ylab='a_g',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(theta[burnt,3], type='l', xlab='', ylab='a_g',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(theta[burnt,4], type='l', xlab='', ylab='b_g',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(theta[burnt,5], type='l', xlab='', ylab='b_g',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(theta[burnt,6], type='l', xlab='', ylab='b_g',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(theta_c[burnt], type='l', xlab='', ylab='theta_c',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(theta_r[burnt], type='l', xlab='', ylab='theta_r',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(theta_s[burnt], type='l', xlab='', ylab='theta_s',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
plot(theta_g[burnt], type='l', xlab='', ylab='theta_g',xaxt='n');abline(v=seq(2001,upto,by=2001),col="grey")
for(i in 1:p) {# p=42
  plot(theta[burnt,3*2+6*L+6*sum(multipleRegionsInSregion)+6*J+3*N+i],type='l', xlab='',ylab=paste('covar', i),xaxt='n')
  abline(v=seq( 2001,upto,by=2001),col="grey")}
for(i in 1:15) { plot(gamma[burnt, i], type='l', xlab='',ylab=paste('gamma', i),xaxt='n')
  abline(v=seq(2001,upto,by=2001),col="grey")}
for(i in 1:15) {plot(log(sigma2[burnt, i]), type='l', xlab='',ylab=paste('log sigma2', i),xaxt='n')
  abline(v=seq(2001,upto,by=2001),col="grey")}
dev.off()

write.csv(deviance, file = paste("deviance.csv"));write.csv(phi_natl,  file = paste("phi_natl.csv"))
write.csv(phi_subn, file = paste("phi_subn.csv"));write.csv(phi_comm,  file = paste("phi_comm.csv"))
write.csv(phi_s,    file = paste("phi_s.csv"))  ;write.csv(phi_r,   file = paste("phi_r.csv"))
write.csv(phi_c,    file = paste("phi_c.csv"))  ;write.csv(eta_s,   file = paste("eta_s.csv")) 
write.csv(eta_r,    file = paste("eta_r.csv"))  ;write.csv(eta_c,   file = paste("eta_c.csv")) 
write.csv(theta,    file = paste("theta.csv"))  ;write.csv(tau,     file = paste("tau.csv"))
write.csv(gamma,    file = paste("gamma.csv"))  ;write.csv(sigma2,  file = paste("sigma2.csv"))
write.csv(theta_c,  file = paste("theta_c.csv"));write.csv(theta_r, file = paste("theta_r.csv"))
write.csv(theta_s,  file = paste("theta_s.csv"));write.csv(theta_g, file = paste("theta_g.csv"))
write.csv(u,        file = paste("u.csv"))      ;write.csv(v,       file = paste("v.csv"))
write.csv(sv,       file = paste("sv.csv"))     ;write.csv(w,       file = paste("w.csv"))

#write.csv(chains1[[1]]$Rage,file = paste("R.age.csv"));write.csv(chains1[[1]]$F.theta.Mm, file = paste("F.theta.Mm.csv"))

