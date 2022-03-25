################################################################################
##### FUNCTIONS ################################################################
################################################################################
##### Gibbs and Gibbs-like update functions ####################################
################################################################################
##### Non-linear trend proposal function #######################################
u.prop.function <- function(i, tPtsPer, V1Inv.diag, Vinv.M, theta.star, theta.old, u.old,P,A,T,eigenNoData,SigmaGenInvNoTheta) {
  #############################################################################
  ### This function is used to propose new values for the non-linear change   #
  ### over time component of the model. It is based on "Gaussian Markov       #
  ### Random Fields: Theory and Applications, Rue H, Held L, Chapman and      #
  ### Hall 2005". A full description of the application to metabolic risks    #
  ### data is provided in Finucane et al and the appendix of Danaei et al.    #
  ### There are four possible scenarios: no data observed in a country, one   #
  ### year of data observed,two years of data observed or more years observed #
  #############################################################################
  val.sw <- paste("val",tPtsPer[i],sep='')                                               # Produces a character string containing the number of time points for a country
  switch(val.sw,                                                                         # Calls appropriate section of function depending on number of time points for a country
         val0 = {                                                                        # No data for that country
           uStar                           <- eigenNoData$vec %*% (1/ sqrt(exp(theta.star)* eigenNoData$val) * rnorm(T^2))
           uStar.old                       <- eigenNoData$vec %*% (1/ sqrt(exp(theta.old) * eigenNoData$val) * rnorm(T^2))
           dens                            <- 1/2 * ((T^2)-3)*theta.star - 1/2 * t(uStar) %*% (exp(theta.star) * SigmaGenInvNoTheta %*% uStar)
           dens.old                        <- 1/2 * ((T^2)-3)*theta.old  - 1/2 * t(u.old[((i-1)*(T^2)+1):(i*(T^2))]) %*% (exp(theta.old) * SigmaGenInvNoTheta %*% u.old[((i-1)*(T^2)+1):(i*(T^2))])
           return(list(u.star=uStar, dens.star=dens, u.star.old=uStar.old, dens.old=dens.old))
         },
         val1 = {                                                                        # One datapoint for that country
           Q                                <- exp(theta.star) * P
           diag(Q)                          <- diag(Q) + V1Inv.diag[((i-1)*(T^2)+1):(i*(T^2))]
           Q.eigen                          <- eigen(Q)
           Q.eigen$values[(T^2-1):(T^2)]    <- Inf
           Qinv                             <- Q.eigen$vectors %*% (1/Q.eigen$values * t(Q.eigen$vectors))
           correct.factor                   <- Qinv %*% matrix(1, T^2, T^2) / sum(Qinv)
           mu                               <- Qinv %*% Vinv.M[((i-1)*(T^2)+1):(i*(T^2))]
           muStar                           <- mu - correct.factor %*% mu
           SigmaStar                        <- Qinv - correct.factor %*% Qinv
           SigmaStar.eigen                  <- eigen(SigmaStar)
           SigmaStar.eigen$val[(T^2-2):(T^2)] <- 0
           #if (any(Re(SigmaStar.eigen$val)<0 & Re(SigmaStar.eigen$val)<=1e-9)){w<-which(Re(SigmaStar.eigen$val)<0 & Re(SigmaStar.eigen$val)<1e-9); SigmaStar.eigen$val[w]<-abs(SigmaStar.eigen$val[w])}
           EigenValGenInv                   <- 1/SigmaStar.eigen$val
           EigenValGenInv[((T^2)-2):(T^2)]  <- 0
           SigmaGenInv                      <- SigmaStar.eigen$vec %*% (EigenValGenInv * t(SigmaStar.eigen$vec))
           uStar                            <- muStar + SigmaStar.eigen$vec %*% (sqrt(SigmaStar.eigen$val) * rnorm(T^2))#T
           #uStar                           <- muStar + SigmaStar.eigen$vec %*% (sqrt(EigenValGenInv) * rnorm(T))
           dens                             <- -1/2 * sum(log(SigmaStar.eigen$val[1:((T^2)-3)])) -1/2 * t(uStar - muStar) %*% (SigmaGenInv %*% (uStar - muStar))
           #dens                            <- -1/2 * sum(log(EigenValGenInv[1:(T-2)])) -1/2 * t(uStar - muStar) %*% (SigmaGenInv %*% (uStar - muStar))
           Q                                <- exp(theta.old) * P
           diag(Q)                          <- diag(Q) + V1Inv.diag[((i-1)*(T^2)+1):(i*(T^2))]#T
           Q.eigen                          <- eigen(Q)
           Q.eigen$values[(T^2-1):(T^2)]    <- Inf
           Qinv                             <- Q.eigen$vectors %*% (1/Q.eigen$values * t(Q.eigen$vectors))
           correct.factor                   <- Qinv %*% matrix(1, (T^2),(T^2)) / sum(Qinv)
           mu                               <- Qinv %*% Vinv.M[((i-1)*(T^2)+1):(i*(T^2))]
           muStar                           <- mu   - correct.factor %*% mu
           SigmaStar                        <- Qinv - correct.factor %*% Qinv
           SigmaStar.eigen                  <- eigen(SigmaStar)
           SigmaStar.eigen$val[(T^2-2):(T^2)] <- 0
           #if (any(Re(SigmaStar.eigen$val)<0 & Re(SigmaStar.eigen$val)<=1e-9)){w<-which(Re(SigmaStar.eigen$val)<0 & Re(SigmaStar.eigen$val)<1e-9); SigmaStar.eigen$val[w]<-abs(SigmaStar.eigen$val[w])}
           EigenValGenInv                   <- 1/SigmaStar.eigen$val
           EigenValGenInv[(T^2-2):(T^2)]    <- 0
           SigmaGenInv                      <- SigmaStar.eigen$vec %*% (EigenValGenInv * t(SigmaStar.eigen$vec))
           uStar.old                        <- muStar + SigmaStar.eigen$vec %*% (sqrt(SigmaStar.eigen$val) * rnorm(T^2))
           #uStar.old                       <- muStaro + SigmaStaro.eigen$vec %*% (sqrt(EigenValGenInvo) * rnorm(T))
           dens.old                         <- -1/2*sum(log(SigmaStar.eigen$val[1:((T^2)-3)])) - 1/2*t(u.old[((i-1)*(T^2)+1):(i*(T^2))] - muStar) %*% (SigmaGenInv %*% (u.old[((i-1)*(T^2)+1):(i*(T^2))]  - muStar))
           #dens.old                        <- -1/2*sum(log(EigenValGenInvo[1:(T-2)])) - 1/2*t(u.old[((i-1)*T+1):(i*T)] - muStaro) %*% (SigmaGenInvo %*% (u.old[((i-1)*T+1):(i*T)] - muStaro))
           return(list(u.star=uStar, dens.star=dens, u.star.old=uStar.old, dens.old=dens.old))
         },
         val2 = {                                                                        # One datapoint for that country
           Q                                <- exp(theta.star) * P
           diag(Q)                          <- diag(Q) + V1Inv.diag[((i-1)*(T^2)+1):(i*(T^2))]
           Q.eigen                          <- eigen(Q)
           Q.eigen$values[T^2]              <- Inf#T
           Qinv                             <- Q.eigen$vectors %*% (1/Q.eigen$values * t(Q.eigen$vectors))
           #correct.factor                  <- Qinv %*% t(A[1:2,]) %*% solve(A[1:2,] %*% Qinv %*% t(A[1:2,])) %*% A[1:2,]
           correct.factor                   <- Qinv %*% t(A[1:2,]) %*% qr.solve(A[1:2,] %*% Qinv %*% t(A[1:2,]),tol=-1e-26) %*% A[1:2,]
           mu                               <- Qinv %*% Vinv.M[((i-1)*(T^2)+1):(i*(T^2))]#T
           muStar                           <- mu - correct.factor %*% mu
           SigmaStar                        <- Qinv - correct.factor %*% Qinv
           SigmaStar.eigen                  <- eigen(SigmaStar)
           SigmaStar.eigen$val[(T^2-2):(T^2)] <- 0
           #if (any(Re(SigmaStar.eigen$val)<0 & Re(SigmaStar.eigen$val)<=1e-9)){w<-which(Re(SigmaStar.eigen$val)<0 & Re(SigmaStar.eigen$val)<1e-9); SigmaStar.eigen$val[w]<-abs(SigmaStar.eigen$val[w])}
           EigenValGenInv                   <- 1/SigmaStar.eigen$val
           EigenValGenInv[((T^2)-2):(T^2)]  <- 0
           SigmaGenInv                      <- SigmaStar.eigen$vec %*% (EigenValGenInv * t(SigmaStar.eigen$vec))
           uStar                            <- muStar + SigmaStar.eigen$vec %*% (sqrt(SigmaStar.eigen$val) * rnorm(T^2))#T
           #uStar                           <- muStar + SigmaStar.eigen$vec %*% (sqrt(EigenValGenInv) * rnorm(T))
           dens                             <- -1/2 * sum(log(SigmaStar.eigen$val[1:((T^2)-3)])) -1/2 * t(uStar - muStar) %*% (SigmaGenInv %*% (uStar - muStar))
           #dens                            <- -1/2 * sum(log(EigenValGenInv[1:(T-2)])) -1/2 * t(uStar - muStar) %*% (SigmaGenInv %*% (uStar - muStar))
           Q                                <- exp(theta.old) * P
           diag(Q)                          <- diag(Q) + V1Inv.diag[((i-1)*(T^2)+1):(i*(T^2))]#T
           Q.eigen                          <- eigen(Q)
           Q.eigen$values[T^2]              <- Inf
           Qinv                             <- Q.eigen$vectors %*% (1/Q.eigen$values * t(Q.eigen$vectors))
           #correct.factor                   <- Qinv %*% t(A[1:2,]) %*% solve(A[1:2,] %*% Qinv %*% t(A[1:2,])) %*% A[1:2,]
           correct.factor                   <- Qinv %*% t(A[1:2,]) %*% qr.solve(A[1:2,] %*% Qinv %*% t(A[1:2,]),tol=-1e-26) %*% A[1:2,]
           mu                               <- Qinv %*% Vinv.M[((i-1)*(T^2)+1):(i*(T^2))]
           muStar                           <- mu   - correct.factor %*% mu
           SigmaStar                        <- Qinv - correct.factor %*% Qinv
           SigmaStar.eigen                  <- eigen(SigmaStar)
           SigmaStar.eigen$val[(T^2-2):(T^2)] <- 0
           #if (any(Re(SigmaStar.eigen$val)<0 & Re(SigmaStar.eigen$val)<=1e-9)){w<-which(Re(SigmaStar.eigen$val)<0 & Re(SigmaStar.eigen$val)<1e-9); SigmaStar.eigen$val[w]<-abs(SigmaStar.eigen$val[w])}
           EigenValGenInv                   <- 1/SigmaStar.eigen$val
           EigenValGenInv[(T^2-2):(T^2)]    <- 0
           SigmaGenInv                      <- SigmaStar.eigen$vec %*% (EigenValGenInv * t(SigmaStar.eigen$vec))
           uStar.old                        <- muStar + SigmaStar.eigen$vec %*% (sqrt(SigmaStar.eigen$val) * rnorm(T^2))
           #uStar.old                       <- muStaro + SigmaStaro.eigen$vec %*% (sqrt(EigenValGenInvo) * rnorm(T))
           dens.old                         <- -1/2*sum(log(SigmaStar.eigen$val[1:((T^2)-3)])) - 1/2*t(u.old[((i-1)*(T^2)+1):(i*(T^2))] - muStar) %*% (SigmaGenInv %*% (u.old[((i-1)*(T^2)+1):(i*(T^2))]  - muStar))
           #dens.old                        <- -1/2*sum(log(EigenValGenInvo[1:(T-2)])) - 1/2*t(u.old[((i-1)*T+1):(i*T)] - muStaro) %*% (SigmaGenInvo %*% (u.old[((i-1)*T+1):(i*T)] - muStaro))
           return(list(u.star=uStar, dens.star=dens, u.star.old=uStar.old, dens.old=dens.old))
         },
         {                                                                           # Two or more datapoints for that country
           Q                               <- exp(theta.star) * P
           diag(Q)                         <- diag(Q) + V1Inv.diag[((i-1)*(T^2)+1):(i*(T^2))]
           Qinv                            <- solve(Q)#qr.solve(Q,tol= 1e-26)
           correct.factor                  <- Qinv %*% t(A) %*% solve(A %*% Qinv %*% t(A)) %*% A
           mu                              <- Qinv %*% Vinv.M[((i-1)*(T^2)+1):(i*(T^2))]
           muStar                          <- mu - correct.factor %*% mu
           SigmaStar                       <- Qinv - correct.factor %*% Qinv
           SigmaStar.eigen                 <- eigen(SigmaStar)
           SigmaStar.eigen$val[(T^2-2):(T^2)]   <- 0 #(T^2-2):(T^2)
           #if (Re(SigmaStar.eigen$val)<0 && Re(SigmaStar.eigen$val)<=1e-9){w<-which(Re(SigmaStar.eigen$val)<0 && Re(SigmaStar.eigen$val)<1e-9); SigmaStar.eigen$val[w]<-abs(SigmaStar.eigen$val[w])}
           EigenValGenInv                  <- 1/(SigmaStar.eigen$val)#Re
           EigenValGenInv[(T^2-2):(T^2)]   <- 0 
           SigmaGenInv                     <- SigmaStar.eigen$vec %*% (EigenValGenInv * t(SigmaStar.eigen$vec))
           uStar                           <- muStar + SigmaStar.eigen$vec %*% (sqrt(SigmaStar.eigen$val) * rnorm(T^2))
           #uStar                          <- muStar + SigmaStar.eigen$vec %*% (sqrt( EigenValGenInv) * rnorm(T))
           dens                            <- -1/2 * sum(log(SigmaStar.eigen$val[1:(T^2-3)])) -1/2 * t(uStar - muStar) %*% (SigmaGenInv %*%  (uStar - muStar))#1:(T^2-3)
           #dens                           <- -1/2 * sum(log( EigenValGenInv[1:(T-2)])) - 1/2 * t(uStar - muStar) %*% (SigmaGenInv %*%  (uStar - muStar))
           Q                               <- exp(theta.old) * P
           diag(Q)                         <- diag(Q) + V1Inv.diag[((i-1)*(T^2)+1):(i*(T^2))]
           Qinv                            <- solve(Q)#qr.solve(Q,tol=1e-26)
           correct.factor                  <- Qinv %*% t(A) %*% solve(A %*% Qinv %*% t(A)) %*% A# qr.solve(aa,tol=1e-26)
           mu                              <- Qinv %*% Vinv.M[((i-1)*(T^2)+1):(i*(T^2))]
           muStar                          <- mu - correct.factor %*% mu
           SigmaStar                       <- Qinv - correct.factor %*% Qinv
           SigmaStar.eigen                 <- eigen(SigmaStar)
           SigmaStar.eigen$val[(T^2-2):(T^2)] <- 0#(T^2-2):(T^2),c(14:17,(T^2-1):(T^2))
           #if (Re(SigmaStar.eigen$val)<0 && Re(SigmaStar.eigen$val)<=1e-9){w<-which(Re(SigmaStar.eigen$val)<0 && Re(SigmaStar.eigen$val)<1e-9); SigmaStar.eigen$val[w]<-abs(SigmaStar.eigen$val[w])}
           EigenValGenInv                  <- 1/(SigmaStar.eigen$val)#Re
           EigenValGenInv[(T^2-2):(T^2)]   <- 0# c(14:17,(T^2-1):(T^2))
           SigmaGenInv                     <- SigmaStar.eigen$vec %*% (EigenValGenInv * t(SigmaStar.eigen$vec))
           uStar.old                       <- muStar + SigmaStar.eigen$vec %*% (sqrt( SigmaStar.eigen$val) * rnorm(T^2))
           #uStar.old                      <- muStar + SigmaStar.eigen$vec %*% (sqrt( EigenValGenInv) * rnorm(T))
           dens.old                        <- -1/2 * sum(log( SigmaStar.eigen$val[1:(T^2-3)]))- 1/2 * t(u.old[((i-1)*(T^2)+1):(i*(T^2))] - muStar) %*% (SigmaGenInv %*% (u.old[((i-1)*(T^2)+1):(i*(T^2))]-muStar))#1:(T^2-3)
           #dens.old                       <- -1/2 * sum(log( EigenValGenInv[1:(T-2)])) - 1/2 * t(u.old[((i-1)*T+1):(i*T)] - muStar) %*% (SigmaGenInv %*% (u.old[((i-1)*T+1):(i*T)] - muStar))
           return(list(u.star=uStar, dens.star=dens, u.star.old=uStar.old, dens.old=dens.old)) } ) }

################################################################################
##### Functions used in Metropolis-Hastings updates                            #
################################################################################

LogLik              <- function(SigmaInv.diag, F.theta.Mm, gamma, R.age,y ){                                        # Function used to calculate log likelihoods
  res3 <- .5*(sum(log(diag(SigmaInv.diag))) - (t(as.vector(y) - F.theta.Mm- R.age %*% gamma))
              %*%as.matrix(diag(SigmaInv.diag)*(as.vector(y) - F.theta.Mm- R.age %*% gamma)) )
  return(res3)
}

ssreLik             <- function(N, phi, ssre) { return((-N*phi - sum(ssre^2)/exp(phi))/2) }
LogPriorPhi         <- function(phi) {  return(phi/2) }

##### Functions used in calculation of the log posteriors for the precision  ###
##### parameters for non-linear trends as described in Danaei et al pp.14-16 ###

CalcThetaPost       <- function(i, uvw,T,P) return(uvw[((i-1)*(T^2)+1):(i*(T^2))] %*% P %*% uvw[((i-1)*(T^2)+1):(i*(T^2))])

LogPostThetaC       <- function(theta_c,u,J,T,P) {
  return((J*(T^2-3)-1) * theta_c / 2 - exp(theta_c)/2 * sum(unlist(lapply(1:J, CalcThetaPost,u,T,P))))
}
LogPostThetaR       <- function(theta_r, v, multipleRegionsInSregion,K,T,P) {
  return((sum(multipleRegionsInSregion)*(T^2-3)-1) * theta_r / 2 - exp(theta_r)/2 * sum(unlist(lapply((1:K)[multipleRegionsInSregion], CalcThetaPost, v,T,P))))
}
LogPostThetaS       <- function(theta_s,sv,L,T,P) {
  return((L*(T^2-3)-1) * theta_s / 2 - exp(theta_s)/2 * sum(unlist(lapply(1:L, CalcThetaPost, sv,T,P))))
}
LogPostThetaG       <- function(theta_g, w,T,P) {
  return((T^2-3-1) * theta_g / 2 - exp(theta_g)/2 * w %*% P %*% w)
}

################################################################################
##### Tuning function to help efficient MCMC updates                           #
##### See Gelman, Roberts, Gilk, Bayesian Statistics 5, pp.599-607             #
################################################################################
adaptJump           <- function(n, pjump, pgoal=NULL, max.mult=2) {#5, me
  pjump[pjump==1]=0.999
  if(is.null(pgoal)) {
    const=rep(log(.44),length(n))                                                   # One dimensional jump
    const[n==2]=log(.35)                                                            # Two dimensional jump
    const[n==3]=log(.32)                                                            # Three dimensional jump
    const[n>3] =log(.25)                                                            # More than three dimensional jump (slightly conservative)
  }
  else const=log(pgoal)
  return(min(max.mult,max(1/max.mult,const/log(pjump))))
}

################################################################################
##### 10 year weighted average function                                        #
################################################################################
TenYearWeightedAvg    <- function(var) {                                                                   # Smooths the country-level covariates using a triangularly-weighted moving average
  avg <- rep(NA, dim(covar)[1])                                                                         # Vector of length equal to the number of rows in the covariate dataset
  for(i in 10:dim(covar)[1]) avg[i] <- weighted.mean(var[(i-9):i], w=(1:10)/sum(1:10))                  # Calculates weighted average
  avg[covar$data_year < (start.year)] <- NA                                                             # Sets to NA for years before the analysis period
  return(avg) }                                                                                         # Returns weighted average and closes function

################################################################################
##### Deviance print function                                                  #
################################################################################
printDeviance         <- function() {
  mean.deviance   <- round(mean(deviance[burnt], na.rm=TRUE),2)                                           # Calculates estimated mean deviance
  p_D             <- round(var(deviance[burnt],  na.rm=TRUE)/2,2)                                         # Calculates effective number of parameters
  dic             <- round(mean(deviance[burnt], na.rm=TRUE) + var(deviance[burnt],na.rm=TRUE)/2, 2)      # Calculates DIC
  cat("mean deviance =", mean.deviance, "\n", "p_d =", p_D, "\n", "p =", p, "\n", "DIC =", dic, "\n") }   # Prints values to screen


################################################################################
##### Trace plot function                                                      #
##### This function prints traceplots for a selection of parameters to a pdf   #
##### document. The list of parameters within the function may be changed      #
##### without affecting the remainder of the code                              #
################################################################################
tracePlots <- function() {
  pdf(paste("_bivariate_traceplots1",".pdf",sep=''), width=10, height=7.5); par(mfrow=c(3,2), mar=c(0.5,4,0.5,0))  # Opens PDF document
  plot(phi_s[take,1],     type='l', xlab='', ylab='phi_s_D')                  # Log variance for normal prior for superregion random intercepts      (log kappa_a^r)
  plot(phi_s[take,2],     type='l', xlab='', ylab='phi_s_S',     xaxt='n')                  # Log variance for normal prior for superregion random intercepts      (log kappa_a^r)
  plot(phi_s[take,3],     type='l', xlab='', ylab='phi_s_INT',   xaxt='n')                  # Log variance for normal prior for superregion random intercepts      (log kappa_a^r)
  plot(phi_r[take,1],     type='l', xlab='', ylab='phi_r_D')                  # Log variance for normal prior for region random intercepts           (log kappa_a^s)
  plot(phi_r[take,2],     type='l', xlab='', ylab='phi_r_S',     xaxt='n')  
  plot(phi_r[take,3],     type='l', xlab='', ylab='phi_r_INT',   xaxt='n')  
  
  plot(phi_c[take,1],     type='l', xlab='', ylab='phi_c_D')                   # Log variance for normal prior for country random intercepts          (log kappa_a^c)
  plot(phi_c[take,2],     type='l', xlab='', ylab='phi_c_S',    xaxt='n')                   # Log variance for normal prior for country random intercepts          (log kappa_a^c)
  plot(phi_c[take,3],     type='l', xlab='', ylab='phi_c_INT',  xaxt='n')                   # Log variance for normal prior for country random intercepts          (log kappa_a^c)
  plot(eta_s[take,1],     type='l', xlab='', ylab='eta_s_D')                   # Log variance for normal prior for superregion random slopes          (log kappa_b^r)
  plot(eta_s[take,2],     type='l', xlab='', ylab='eta_s_S',    xaxt='n')
  plot(eta_s[take,3],     type='l', xlab='', ylab='eta_s_INT',  xaxt='n')
  plot(eta_r[take,1],     type='l', xlab='', ylab='eta_r_D')                   # Log variance for normal prior for region random slopes               (log kappa_b^s)
  plot(eta_r[take,2],     type='l', xlab='', ylab='eta_r_S',    xaxt='n')                   # Log variance for normal prior for region random slopes               (log kappa_b^s)
  plot(eta_r[take,3],     type='l', xlab='', ylab='eta_r_INT',  xaxt='n')                   # Log variance for normal prior for region random slopes               (log kappa_b^s)
  
  plot(eta_c[take,1],     type='l', xlab='', ylab='eta_c_D')                   # Log variance for normal prior for country random slopes              (log kappa_b^c)
  plot(eta_c[take,2],     type='l', xlab='', ylab='eta_c_S',    xaxt='n')                   # Log variance for normal prior for country random slopes              (log kappa_b^c)
  plot(eta_c[take,3],     type='l', xlab='', ylab='eta_c_INT',  xaxt='n')                   # Log variance for normal prior for country random slopes              (log kappa_b^c)
  plot(phi_natl[take,1],  type='l', xlab='', ylab='phi_natl_D')                   # Log variance of random effects for national studies (analogous to collapsed log nu_w and log nu_u)
  plot(phi_natl[take,2],  type='l', xlab='', ylab='phi_natl_S', xaxt='n')                   # Log variance of random effects for national studies (analogous to collapsed log nu_w and log nu_u)
  plot(phi_natl[take,3],  type='l', xlab='', ylab='phi_natl_INT',xaxt='n')                  # Log variance of random effects for national studies (analogous to collapsed log nu_w and log nu_u)
  plot(phi_subn[take,1],  type='l', xlab='', ylab='phi_subn_D')                  # Log variance of random effects for subnational studies               (log nu_s)
  plot(phi_subn[take,2],  type='l', xlab='', ylab='phi_subn_S',  xaxt='n')                  # Log variance of random effects for subnational studies               (log nu_s)
  plot(phi_subn[take,3],  type='l', xlab='', ylab='phi_subn_INT',xaxt='n')                  # Log variance of random effects for subnational studies               (log nu_s)
  plot(phi_comm[take,1],  type='l', xlab='', ylab='phi_comm_D')                  # Log variance of random effects for community studies                 (log nu_c)
  plot(phi_comm[take,2],  type='l', xlab='', ylab='phi_comm_S',  xaxt='n')                  # Log variance of random effects for community studies                 (log nu_c)
  plot(phi_comm[take,3],  type='l', xlab='', ylab='phi_comm_INT',xaxt='n')                  # Log variance of random effects for community studies                 (log nu_c)
  plot(tau[take,1],       type='l', xlab='', ylab='tau_D')                   # Log variance for within-study errors that differ between age groups  (log tau)
  plot(tau[take,2],       type='l', xlab='', ylab='tau_S',      xaxt='n')                   # Log variance for within-study errors that differ between age groups  (log tau)
  plot(tau[take,3],       type='l', xlab='', ylab='tau_INT',    xaxt='n')                   # Log variance for within-study errors that differ between age groups  (log tau)
  plot(theta_c[take],   type='l', xlab='', ylab='theta_c')                   # Log precision parameter for random walk at country level             (log lambda_c)
  plot(theta_r[take],   type='l', xlab='', ylab='theta_r',      xaxt='n')                   # Log precision parameter for random walk at region level              (log lambda_s)
  plot(theta_s[take],   type='l', xlab='', ylab='theta_s')                   # Log precision parameter for random walk at superregion level         (log lambda_r)
  plot(theta_g[take],   type='l', xlab='', ylab='theta_g',      xaxt='n')                   # Log precision parameter for random walk at global level              (log lambda_g)
  plot(deviance[take],  type='l', xlab='', ylab='deviance')                   # Deviance
  plot(theta[take,1], type='l', xlab='', ylab='a_gD')                   # Global random intercept
  plot(theta[take,2], type='l', xlab='', ylab='a_gS',           xaxt='n')                   # Global random intercept
  plot(theta[take,3], type='l', xlab='', ylab='a_gINT',         xaxt='n')                   # Global random intercept
  plot(theta[take,4], type='l', xlab='', ylab='b_gD')                   # Global random slope
  plot(theta[take,5], type='l', xlab='', ylab='b_gS',           xaxt='n')                   # Global random slope
  plot(theta[take,6], type='l', xlab='', ylab='b_gINT',         xaxt='n')                   # Global random slope
  #super-regions_intercept
  for ( i in 7:15) plot(theta[take,i], type='l', xlab=i-6 , ylab= "a_super-region_D")                        # Global random slope
  for ( i in 16:24) plot(theta[take,i], type='l', xlab=i-15 , ylab= "a_super-region_S",   xaxt='n')                      # Global random slope
  for ( i in 25:33) plot(theta[take,i], type='l', xlab=i-24 , ylab= "a_super-region_INT",   xaxt='n')                    # Global random slope
  #super-regions_slope
  for ( i in 34:42) plot(theta[take,i], type='l', xlab=i-33 , ylab= "b_super-region_D")                      # Global random slope
  for ( i in 43:51) plot(theta[take,i], type='l', xlab=i-42 , ylab= "b_super-region_S",   xaxt='n')                      # Global random slope
  for ( i in 52:60) plot(theta[take,i], type='l', xlab=i-51 , ylab= "b_super-region_INT",   xaxt='n')                    # Global random slope
  
  #for(i in 1:p)   plot(theta[,2+2*L+2*sum(multipleRegionsInSregion)+2*J+N+i], type='l', xlab='', ylab=paste('covar', i),xaxt='n')   # Covariate terms  (beta matrix)
  # The theta[,2+2*L+2*sum(multipleRegionsInSregion)+2*J] terms are the linear intercepts and slopes at global, superregion, region and country level respectively
  for(i in 1:12)  plot(gamma[take,i], type='l', xlab='', ylab=paste('gamma', i))                                # Age model parameters (psi and phi terms)
  for(i in c(1,6,11))   plot(log(sigma2[take,i]), type='l', xlab='', ylab=paste('sigma2', i),xaxt='n')                   # Log variances for country-specific random spline coefficients (sigma^2 terms)
  
  for (i in 1:(T^2))     plot(w[take,i],type='l', ylab= paste("w",i))
  for (i in 1:(2*(T^2))) plot(sv[take,i],type='l', ylab=paste("sv",i))
  for (i in 1:(2*(T^2))) plot(v[take,i],type='l', ylab= paste("v",i))
  for (i in 1:(2*(T^2))) plot(u[take,i],type='l', ylab= paste("u",i))
  dev.off()
}
#### END OF FUNCTIONS ##########################################################
##### PRELIMINARY CODE #########################################################
library(MASS)
library(reshape2)
library(Matrix)
library(splines)
library(zoo)
library(Runuran)
library(iterators)
library(parallel)

##### USER DEFINED INPUTS ######################################################
##### Users should set appropriate values for the objects below                #
################################################################################
MCMC <- function(){
  seedVal                                         <- 1                                                                    # Assign each chain a unique seed value
  sex                                             <- "male"                                                               # Set gender: must be either "male" or "female"
  variable                                        <- "bp"                                                                 # Set variable name
  start.year                                      <- 2004                                                # Set start year for analysis
  end.year                                        <- 2014                                                    # Set end year for analysis
  knot1                                           <- 46   #45                                                             # Position of first knot in cubic splines
  middle.age                                      <- 51   #50                                                             # Middle age in analysis
  knot2                                           <- 61   #60                                                             # Position of second knot in cubic splines
  centring.year                                   <- 2008#5 #1997                                                           # Used to centre time
  dbp.centring.value                              <- 78                                                                   # Used to centre DBP
  sbp.centring.value                              <- 130                                                                  # Used to centre SBP
  int.centring.value                              <- 10277#35                                                             # Used to centre SBP
  nLong                                           <- 100000                                                               # Set number of iterations for MCMC loop
  mod.no                                          <- 1                                                                    # Set a unique model identifier: this can be any number or a character string
  epsilon                                         <- .01                                                            # Flat priors on the precision scale
  freq.val                                        <- 200                                                                  # Number of iterations per tuning step in the first 5000 MCMC iterations
  covariate.list                                  <- c("perurb","pc1","pc2","pc3","pc4")                                  #  not taking into account gdp and education
  #set.seed(seedVal)                                                                                                      # Sets random number seed using seedVal object
  sex.val                                         <- switch(sex, male = 1, female = 2)                                    # Sets a numerical variable for sex: 1 for male and 2 for female
  nIts                                            <- nLong/10                                                             # Number of MCMC iterations after thinning
  filename                                        <- paste0("Model",mod.no,"_mean_",sex,"_",variable,"_Seed_",seedVal)    # Sets a base filename used throughout the code
  
  ##### COVARIATES INPUT   #####################################################
  covar                                           <- read.csv("covariates_2017_01_11.csv", header=TRUE)
  covar                                           <- covar[order(covar$region,covar$country,covar$data_year),]                                                  # Sort covariate dataset by region, country and year
  head(covar)
  TenYearWeightedAvg                              <- function(var) {                                                      # Smooths the country-level covariates using a triangularly-weighted moving average
    avg <- rep(NA, dim(covar)[1])                                                                                         # Vector of length equal to the number of rows in the covariate dataset
    for(i in 10:dim(covar)[1]) avg[i] <- weighted.mean(var[(i-9):i], w=(1:10)/sum(1:10))                                  # Calculates weighted average
    avg[covar$data_year < (start.year)] <- NA                                                                             # Sets to NA for years before the analysis period
    return(avg) } 
  Z.weighted                                      <- apply(covar[,noquote(covariate.list)], 2, TenYearWeightedAvg)                                              # Take ten-year weighted average of covariate data
  covar_imp                                       <- data.frame(covar$country,covar$data_year, covar$region, covar$sregion, Z.weighted)                                                    # Create an object with ten-year weighted averages for covariates
  covar                                           <- subset(covar_imp, covar$data_year >= start.year & covar$data_year <= end.year) #start.year                                        # Subset this object by analysis period
  colnames(covar)                                 <- c("country","data_year","region","sregion","perurb","pc1","pc2","pc3","pc4") 
  covar[,noquote(covariate.list)[-1]]             <- scale(covar[,noquote(covariate.list)[-1]], center=TRUE, scale=FALSE)                                       # Centre covariates, except urbanisation for identifiability
  
  ##### DATA INPUT ###############################################################
  
  data                                            <- read.csv("BP_STEPS_dataset_2021_07_21.csv")#Blood Pressure
  data                                            <- data[,c("id_study","sex","age_mean","mean_sbp","mean_dbp","se_sbp","se_dbp","mid_year","survey_type",       # Columns used in analysis
                                                             "urban_rural","Country","Region","Superregion")]
  subset                                          <- data[data$sex==sex.val,]                                                                              # Subset data for gender being modelled
  colnames(subset)                                <- c("uid","sex.final","age","mean_sbp","mean_dbp","se_sbp","se_dbp","data_year","coverage","scope","country","region","sregion")    # Rename columns for analysis
  subset                                          <- subset[(subset$data_year>=start.year) & (subset$data_year<= end.year),]
  
  ### attention in the covariance I will take into account only the countries that are in total in the data as we want to estimate some extra countries beyond 
  ### the dataset
  
  a <- NULL; for (i in 1:length(unique(data$Country))) a <- c(a,c(which(as.character(covar$country)== as.character(unique(data$Country)[i])) ))# countr[i]
  covar         <- covar[a,]
  covar$region  <- as.character(covar$region)
  covar$sregion  <- as.character(covar$sregion)
  ##############################
  
  # Make columns in 'subset' directly available
  y1                                              <- subset$mean_dbp - dbp.centring.value                                                                        # DBP centred for computational stability
  y2                                              <- subset$mean_sbp - sbp.centring.value                                                                        # SBP centred for computational stability
  y3                                              <- subset$mean_dbp*subset$mean_sbp - int.centring.value
  y                                               <- as.vector(cbind(y1,y2,y3)) 
  
  ##### Derive variables related to urbanisation #################################
  
  rural                                           <- as.numeric( subset$scope=="Rural")                                                                          # Indicator variable for rural-only data
  urban                                           <- as.numeric( subset$scope=="Urban")                                                                          # Indicator variable for urban-only data
  ##### Derive variables for the age model #######################################
  t1                                              <- substr(start.year,3,4):substr(end.year,3,4)
  t                                               <- cbind(rep(start.year:end.year,each=length(t1)),rep(start.year:end.year,length(t1)))- centring.year
  t                                               <- cbind(t,1:(length(t1)^2)- ceiling((length(t1)^2)/2))#13
  #time                                           <- subset$data_year - centring.year                                                                            # Centred timepoints in data
  scal                                            <- cbind(subset$data_year,subset$data_year) - centring.year    
  time <- NULL; for (i in 1:dim(scal)[1]) time[i] <- which(t[,1]==scal[i] & t[,2]==scal[i])- ceiling((length(t1)^2)/2)#13
  plot(density(time))                                         #smoother density (non- parametric approach)
  subset$age                                      <- subset$age - middle.age                                                                                     # Age centred by middle age
  plot(density(subset$age))
  age.knot1                                       <- subset$age + middle.age - knot1                                                                             # instead of centralizing the age by the mean we use the knot1
  age.knot1[subset$age< -middle.age + knot1]      <- 0                                                                       # Age using first knot
  plot(density(age.knot1))
  age.knot2                                       <- subset$age + middle.age - knot2
  age.knot2[subset$age< -middle.age + knot2]      <- 0                                                                       # Age using second knot
  plot(density(age.knot2))
  #ageMat                                         <- cbind(subset$age, subset$age^2, subset$age^3, age.knot1^3, age.knot2^3) # Age matrix, mapping datapoints to ages
  ageMat                                          <- bs(subset$age,df = 5)[,1:5]
  ageMat.prime                                    <- t(ageMat)                                                               # Transposed age matrix
  
  ##### Calculate constants ######################################################
  #attach(subset) 
  I                                               <- length(subset$mean_sbp)               # mean_dbp# Number of survey-year-sex-age datapoints
  J                                               <- length(unique(covar$country))                   # Number of countries for which estimates are made (as opposed to number of countries in the dataset)
  K                                               <- length(unique(covar$region))                    # Number of regions for which estimates are made
  L                                               <- length(unique(covar$sregion))                   # Number of superregions for which estimates are made
  N                                               <- length(unique(subset$uid))                      # Number of studies in data
  T                                               <- length(t1)                            # Length of analysis period
  ##### RANDOM WALK CONSTANTS ####################################################
  ##### Penalty matrix ###########################################################
  ##### Non Linear Change at Time ################################################ 
  ##### See Rue and Held, p.110, and Danaei et al Webappendix, p.6
  
  mtx3 <- diag(1,T) ;mtx1 <- mtx2 <- mtx4 <- mtx5 <- mtx6 <- matrix(0,T,T)
  # 2nd order random walk
  mtx1[1,1:3] <- mtx1[T,T:(T-2)]                  <- c(6,-5,1)
  mtx1[2,1:4] <- mtx1[(T-1),T:(T-3)]              <- c(-5,12,-6,1)
  for (i in 3:(T-2)) mtx1[i,(i-2):(i+2)]          <- c(1,-6,12,-6,1)
  #1st block mtx boundary, (1st order random walk)
  mtx2[1,1:2] <-  mtx2[T,T:(T-1)]                 <- c(-5,2)
  for (i in 2:(T-1)) mtx2[i,(i-1):(i+1)]          <- c(2,-7,2)
  # 2nd order random walk
  mtx4[1,1:3] <- mtx4[T,T:(T-2)]                  <- c(12,-7,1)
  mtx4[2,1:4] <- mtx4[(T-1),T:(T-3)]              <- c(-7,20,-8,1)
  for (i in 3:(T-2)) mtx4[i,(i-2):(i+2)]          <- c(1,-8,20,-8,1)
  #1st block mtx boundary, (1st order random walk)
  mtx5[1,1:2] <-  mtx5[T,T:(T-1)]                 <- c(-6,2)
  for (i in 2:(T-1)) mtx5[i,(i-1):(i+1)]          <- c(2,-8,2)
  
  R2D                                             <- matrix(0,T^2,T^2)
  R2D[1:T,1:(3*T)]                                <- cbind(mtx1,mtx2,mtx3)
  R2D[((T-1)*T+1):(T^2),((T-3)*T+1):(T^2)]        <- cbind(mtx3,mtx2,mtx1)
  R2D[(T+1):(2*T),1:(4*T)]                        <- cbind(mtx2,mtx4,mtx5,mtx3)
  R2D[((T-2)*T+1):((T-1)*T),((T-4)*T+1):(T^2)]    <- cbind(mtx3,mtx5,mtx4,mtx2)
  for ( i in 3:(T-2)) R2D[((i-1)*T+1):(i*T),((i-3)*T+1):((i+2)*T)] <- cbind(mtx3,mtx5,mtx4,mtx5,mtx3)
  
  P                                               <- R2D#DD
  ##### Eigenanalysis for countries with no data #################################
  ##### Random walk with no data - see Danaei et al Webappendix p.15
  
  ##Since the random walk is a second order of 2Dim we should restrict the third of the end =0 too 
  ## But it seems as 8 of them should change.They are close to 0.
  eigenNoData                                     <- eigen(P)
  eigenNoData$values[(T^2-2):(T^2)]               <- Inf 
  EigenValGenInv                                  <- eigenNoData$val 
  EigenValGenInv[(T^2-2):(T^2)]                   <- 0  
  SigmaGenInvNoTheta                              <- eigenNoData$vec %*% (EigenValGenInv * t(eigenNoData$vec))
  
  ##### Constraints for random walk ##############################################
  ##### See Danaei et al Webappendix p.15
  ##### Allows identifiabilty of linear random intercepts and slopes
  A                                               <- matrix(1, 3, T^2)                                        # T,49
  #tnew                                           <- rep(t,each=T)                                           # every combination for each year of diastolic to each year of systolic.for the 
  #A[2,]                                          <- tnew - mean(tnew)                                       # when the year is 1990 for diastolic i will repeat it 26 times for 1990-2015 for systolic. 
  # It's always the combination
  A[2,]                                           <- t[,1] - mean(t[,1])                                     # tnew - mean(tnew). I have already dealed with the rep for SBP,DBP and now I ll take 
  A[3,]                                           <- t[,2] - mean(t[,2])                                     # tnew - mean(tnew). I have already dealed with the rep for SBP,DBP and now I ll take 
  #A[2,]                                          <- t - mean(t)                                             # the first column.
  ##### MAPPING MATRICES AND VECTORS #############################################
  #subset
  uid.match                                       <- data.frame(name=unique(subset$uid)[order(unique(subset$uid))],number=1:length(unique(subset$uid)))           # List of all study IDs
  uid.match$coverage                              <- subset$coverage[match(uid.match$name,subset$uid)]                                                            # Adds coverage to list of all study IDs
  #covariates
  country.match                                   <- data.frame(name=unique(covar$country)[order(unique(covar$country))], number=1:length(unique(covar$country)))
  region.match                                    <- data.frame(name=unique(covar$region)[order(unique(covar$region))], number=1:length(unique(covar$region)))
  sregion.match                                   <- data.frame(name=unique(covar$sregion)[order(unique(covar$sregion))], number=1:length(unique(covar$sregion)))
  country.time.match                              <- data.frame(country = rep(as.character(sort(unique(covar$country))),each=T^2),yearD = as.numeric(rep(rep(start.year:end.year,each=T),J)), yearS = as.numeric(rep(start.year:end.year,J)), number = 1:(J*(T^2)))
  region.time.match                               <- data.frame(region = rep(as.character(sort(unique(covar$region))),each=T^2),yearD = as.numeric(rep(rep(start.year:end.year,each=T),K)), yearS = as.numeric(rep(start.year:end.year,K)), number = 1:(K*(T^2)))
  sregion.time.match                              <- data.frame(sregion = rep(as.character(sort(unique(covar$sregion))),each=T^2),yearD = as.numeric(rep(rep(start.year:end.year,each=T),L)), yearS = as.numeric(rep(start.year:end.year,L)), number = 1:(L*(T^2)))
  
  ##### Derive mapping matrices and vectors
  which.uid                                       <- match(subset$uid,uid.match$name)                                                        # Identifies position of each datapoint's study ID in the ID list 
  #p.x sthn deuterh 8esh tou uid.match$name vrisketai to prwto stoixeio tou uid opou emfanizetai 5 fores sthn 1h 8esh
  which.country                                   <- match(subset$country, country.match$name)                                                                 # Identifies position of each datapoint's country in the country list
  which.region                                    <- match(subset$region, region.match$name)                                                                   # Identifies position of each datapoint's region in the region list
  which.sregion                                   <- subset$sregion                                                                                           # Identifies position of each datapoint's superregion in sregion (factor)
  #which.time                                     <- time-min(time)+1
  which.time                                      <- match(paste(subset$data_year,subset$data_year),paste(country.time.match$yearD,country.time.match$yearS))  # Identifies position of each datapoint's year in list of years
  which.countryTime                               <- match(paste(subset$country,subset$data_year,subset$data_year),paste(country.time.match$country,country.time.match$yearD,country.time.match$yearS))  # Identifies position of each datapoint's country-year
  which.countryTime.covar                         <- match(paste(subset$country,subset$data_year,subset$data_year),paste(covar$country,covar$data_year,covar$data_year))                      # Identifies position of each datapoint's country-year in covariates
  which.regionTime                                <- match(paste(subset$region,subset$data_year,subset$data_year),paste(region.time.match$region,region.time.match$yearD,region.time.match$yearS))     # Identifies position of each datapoint's region-year
  which.sregionTime                               <- match(paste(subset$sregion,subset$data_year,subset$data_year),paste(sregion.time.match$sregion,sregion.time.match$yearD,sregion.time.match$yearS))  # Identifies position of each datapoint's superregion-year
  multipleRegionsInSregion                        <- as.logical(!rowSums(table(covar$region,covar$sregion)[,(1:L)[colSums(table(covar$region,covar$sregion)>0)==1]])>0)
  
  # Identifies superregions comprising multiple regions
  country.table                                   <- table(index(which.country),which.country)                                                                # Table identifying countries corresponding to each datapoint
  StudyToCountry                                  <- StudyToCountry.time <- matrix(0, I, J)                                                                   # Matrix mapping datapoints to countries
  StudyToCountry[,as.numeric(colnames(country.table))] <- country.table[,as.numeric(index(colnames(country.table)))]                                     # Matrix mapping datapoints to countries
  StudyToCountry.time                             <- StudyToCountry * time                                                                                    # Matrix mapping datapoint-time combinations to countries
  StudyToRegion                                   <- table(index(subset$region),subset$region)                                                                # Matrix mapping datapoints to regions
  StudyToRegion.time                              <- table(index(subset$region),subset$region) * time                                                         # Matrix mapping datapoint-time combinations to regions
  StudyToSregion                                  <- table(index(subset$sregion),subset$sregion)                                                              # Matrix mapping datapoints to superregions
  StudyToSregion.time                             <- table(index(subset$sregion),subset$sregion) * time                                                       # Matrix mapping datapoint-time combinations to superregions
  StudyToUid                                      <- table(index(which.uid),which.uid)                                                                        # Matrix mapping datapoints to study IDs
  StudyToUid.time                                 <- table(index(which.uid),which.uid) * time                                                                 # Matrix mapping datapoint-time combinations to study IDs
  ##### Create "F" matrix ########################################################
  
  F1                                              <- cbind (
    matrix(1,I,3),                                                                                                     # Column in F matrix corresponding to global intercept DBP, SBP, Interaction
    matrix(time,I,3),                                                                                                  # Column in F matrix corresponding to global slope DBP, SBP, Interaction
    
    matrix(StudyToSregion,I,3*dim(StudyToSregion)[2]),                                                                 # Columns in F matrix corresponding to superregion intercepts INTER
    matrix(StudyToSregion.time,I,3*dim(StudyToSregion.time)[2]),                                                       # Columns in F matrix corresponding to superregion slopes INTER
    
    matrix(StudyToRegion[,multipleRegionsInSregion],I,3*dim(StudyToRegion[,multipleRegionsInSregion])[2]),             # Columns in F matrix corresponding to region intercepts INTER
    matrix(StudyToRegion.time[,multipleRegionsInSregion],I,3*dim(StudyToRegion.time[,multipleRegionsInSregion])[2]),   # Columns in F matrix corresponding to region slopes  DBP
    
    matrix(StudyToCountry,I,3*dim(StudyToCountry)[2]),                                                                 # Columns in F matrix corresponding to country intercepts DBP
    matrix(StudyToCountry.time, I,3*dim(StudyToCountry.time)[2]),                                                      # Columns in F matrix corresponding to country slopes DBP
    
    matrix(StudyToUid,I, 3*dim(StudyToUid)[2]),                                                                        # Columns in F matrix corresponding to study random effects, ei INTER
    
    matrix(as.numeric(subset$coverage=="Subnational"),I,3),                                                            # Column in F matrix corresponding to subnational offset DBP
    matrix(as.numeric(subset$coverage=="Subnational")*time,I,3),                                                       # Column in F matrix corresponding to subnational slope offset DBP
    
    matrix(as.numeric(subset$coverage=="Community"),I,3),
    matrix(as.numeric(subset$coverage=="Community")*time,I,3),                                                         # Column in F matrix corresponding to community slope offset DBP
    
    matrix(covar[,covariate.list[1]][which.countryTime.covar],I,3),                                                    # Column in F matrix corresponding to urbanisation offset DBP
    matrix(covar[,covariate.list[1]][which.countryTime.covar]*time,I,3),                                               # Column in F matrix corresponding to urbanisation slope offset INTER
    
    matrix(covar[,covariate.list[1]][which.countryTime.covar]*rural,I,3),                                              # Column in F matrix corresponding to rural-only offset INTER
    matrix(covar[,covariate.list[1]][which.countryTime.covar]*rural*time,I,3),                                         # Column in F matrix corresponding to rural-only slope offset DBP
    
    matrix((1-covar[,covariate.list[1]][which.countryTime.covar])*urban,I,3),                                          # Column in F matrix corresponding to urban-only offset DBP
    matrix((1-covar[,covariate.list[1]][which.countryTime.covar])*urban*time,I,3),                                     # Column in F matrix corresponding to urban-only slope offset DBP
    
    matrix(covar[,covariate.list[2]][which.countryTime.covar],I,3),                                                    # Column in F matrix corresponding to first food balance covariate offset DBP
    matrix(covar[,covariate.list[3]][which.countryTime.covar],I,3),                                                    # Column in F matrix corresponding to second food balance covariate offset DBP
    matrix(covar[,covariate.list[4]][which.countryTime.covar],I,3),                                                    # Column in F matrix corresponding to third food balance covariate offset DBP
    matrix(covar[,covariate.list[5]][which.countryTime.covar],I,3)                                                     # Column in F matrix corresponding to fourth food balance covariate offset DBP
  )
  F                                               <- rbind(F1,F1,F1)
  all                                             <- 6*(1+L+sum(multipleRegionsInSregion)+J)
  
  #DBP
  F[(dim(F1)[1]+1):(3*dim(F1)[1]),c(1,4,7:(6+L),(7+3*L):(6+4*L),(7+6*L):(6+6*L+sum(multipleRegionsInSregion)),(7+6*L+3*sum(multipleRegionsInSregion)):(6+6*L+4*sum(multipleRegionsInSregion)),
                                    (7+6*L+6*sum(multipleRegionsInSregion)):(6+6*L+6*sum(multipleRegionsInSregion)+J),(7+6*L+6*sum(multipleRegionsInSregion)+3*J):(6+6*L+6*sum(multipleRegionsInSregion)+4*J),
                                    (all+1):(all+N),all+1+3*N+seq(0,3*13,3))]         <- 0 
  #SBP
  F[c(1:dim(F1)[1],(2*dim(F1)[1]+1):(3*dim(F1)[1]) ) ,c(2,5,(7+L):(6+2*L),(7+4*L):(6+5*L),(7+6*L+sum(multipleRegionsInSregion)):(6+6*L+2*sum(multipleRegionsInSregion)),(7+6*L+4*sum(multipleRegionsInSregion)):(6+6*L+5*sum(multipleRegionsInSregion)),
                                                        (7+6*L+6*sum(multipleRegionsInSregion)+J):(6+6*L+6*sum(multipleRegionsInSregion)+2*J),(7+6*L+6*sum(multipleRegionsInSregion)+4*J):(6+6*L+6*sum(multipleRegionsInSregion)+5*J),
                                                        (all+1+N):(all+2*N),all+2+3*N+seq(0,3*13,3))]     <- 0   
  #INTER
  F[1:(2*dim(F1)[1]),c(3,6,(7+2*L):(6+3*L),(7+5*L):(6+6*L),(7+6*L+2*sum(multipleRegionsInSregion)):(6+6*L+3*sum(multipleRegionsInSregion)),(7+6*L+5*sum(multipleRegionsInSregion)):(6+6*L+6*sum(multipleRegionsInSregion)),
                       (7+6*L+6*sum(multipleRegionsInSregion)+2*J):(6+6*L+6*sum(multipleRegionsInSregion)+3*J),(7+6*L+6*sum(multipleRegionsInSregion)+5*J):(6+6*L+6*sum(multipleRegionsInSregion)+6*J),
                       (all+1+2*N):(all+3*N),all+3+3*N + seq(0,3*13,3))] <- 0
  
  F                                               <- Matrix(F, sparse=TRUE)
  p                                               <- dim(F)[2] - (3*2+6*L+6*sum(multipleRegionsInSregion)+6*J+ 3*N)     # Number of covariate parameters in model
  F.prime                                         <- t(F)                                                               # Object derived from the F matrix for use below 
  F.age                                           <- F*subset$age                                                       # Object derived from the F matrix for use below
  F.age.prime                                     <- t(F.age)                                                           # Object derived from the F matrix for use below
  
  ##### INITIAL VALUES ###########################################################
  phi_s <- phi_r <- phi_c                         <- matrix(NA, nIts,3)                                                 # Empty vector for log variance for normal prior for superregion random intercepts (log kappa_a^r)
  eta_s <- eta_r <- eta_c                         <- matrix(NA, nIts,3)                                                 # Empty vector for log variance for normal prior for superregion random slopes (log kappa_b^r)
  phi_natl <- phi_subn <- phi_comm                <- matrix(NA, nIts,3)                                                 # Empty vector for log variance of random effects for national studies SBP   (log nu_"national")
  #### Var(ei)
  tau                                             <- matrix(NA, nIts,3)                                                 # Empty vector for log variance for within-study errors that differ between age groups (log tau) DBP
  deviance                                        <- rep(NA, nIts)                                                      # Empty vector for deviance
  gamma                                           <- matrix(NA, nIts, 3*5+3*5+3*5*J)                                    # Empty matrix for age model parameters
  sigma2                                          <- matrix(NA, nIts, 3*5)                                              # Empty matrix for variances for country-specific random spline coefficients
  theta                                           <- matrix(NA, nIts, dim(F)[2])                                        # Empty matrix for theta matrix (ai,bi,betaX, ei)
  ###non linear change over time
  theta_c                                         <- rep(NA, nIts)                                                      # Empty vector for log precision parameter for random walk at country level     (log lambda_c)
  theta_r                                         <- rep(NA, nIts)                                                      # Empty vector for log precision parameter for random walk at region level      (log lambda_s)
  theta_s                                         <- rep(NA, nIts)                                                      # Empty vector for log precision parameter for random walk at superregion level (log lambda_r)
  theta_g                                         <- rep(NA, nIts)                                                      # Empty vector for log precision parameter for random walk at global level      (log lambda_g)
  # we need to have J*T for each of DBP,SBP,Inter. So, J*T^2
  u                                               <- matrix(NA, nIts, J*(T^2))#J*(T^2)                                  # Empty matrix for component of nonlinear trend at national level
  v                                               <- matrix(NA, nIts, K*(T^2))#T,T*2-3                                  # Empty matrix for component of nonlinear trend at regional level
  sv                                              <- matrix(NA, nIts, L*(T^2))#T,T*2-3                                  # Empty matrix for component of nonlinear trend at superregional level
  w                                               <- matrix(NA, nIts, (T^2))#T,T*2-3                                    # Empty matrix for component of nonlinear trend at global level
  
  ##### The following lines set initial values for model parameters #############
  phi_c.prop.sd <- phi_r.prop.sd <- phi_s.prop.sd <- rep(.6,3)                             # Initial value for proposal SD for log variance for normal prior for country random intercepts (log kappa_a^c)
  eta_c.prop.sd <- eta_r.prop.sd <- eta_s.prop.sd <- rep(.6,3) #0.5                          # Initial value for proposal SD for log variance for normal prior for country random slopes (log kappa_b^c)
  phi_natl.prop.sd <- phi_subn.prop.sd <- phi_comm.prop.sd  <- c(.6,.6,.6)                    # Initial value for proposal SD for log variance of random effects for national studies (log nu_"national")
  tau.prop.sd                                     <- c(.2,.2,.2)#.04,.05,.05                # Initial value for proposal SD for log variance for within-study errors that differ between age groups (log tau) DBP,SBP,INTER
  log.sigma2.prop.sd                              <- rep(.7,5*3)                               # Initial values for proposal SDs for log variances for country-specific random spline coefficients
  theta_c.prop.sd                                 <- 1.83#1.8       #for 5 years .8                                  # Initial value for proposal SD for log precision parameter for random walk at country level (log lambda_c)
  theta_r.prop.sd                                 <- 2.44#2                                   # Initial value for proposal SD for log precision parameter for random walk at region level (log lambda_s)
  theta_s.prop.sd                                 <- 3.24#1.3    #1.5                                     # Initial value for proposal SD for log precision parameter for random walk at superregion level (log lambda_r)
  theta_g.prop.sd                                 <- 7.23#6.5                                         # Initial value for proposal SD for log precision parameter for random walk at global level (log lambda_g)
  phi_s[1,] <- phi_r[1,] <- phi_c[1,]             <- rep(1,3)
  eta_s[1,] <- eta_r[1,] <- eta_c[1,]             <- rep(1,3)
  phi_natl[1,]                                    <- rep(4,3)                                  # Initial value for lohttps://emmy.kent.ac.uk/graphics/8e2ebf57-2b7a-46e0-a822-9fddc9a76f76.pngg variance of random effects for national studies (log nu_"national")
  phi_subn[1,]                                    <- rep(4,3)
  phi_comm[1,]                                    <- rep(4,3)  
  tau[1,1:3]                                      <- rep(1,3)  #2.5                            # Initial value for log variance for within-study errors that differ between age groups (log tau)
  gamma[1,]                                       <- rep(0, 3*5+3*5+3*5*J)                     # Initial values for age model parameters
  sigma2[1,]                                      <- rep(1e-10,3*5)                            # Initial values for variances for country-specific random spline coefficients
  theta[1,]                                       <- rep(0, dim(F)[2])                         # Initial values for theta matrix
  theta_c[1]                                      <- 7#    SBP:7                               # Initial value for log precision parameter for random walk at country level (log lambda_c)
  theta_r[1]                                      <- 8#    SBP:8                               # Initial value for log precision parameter for random walk at region level (log lambda_s)
  theta_s[1]                                      <- 10#,14 SBP:10                             # Initial value for log precision parameter for random walk at superregion level (log lambda_r)
  theta_g[1]                                      <- 12#,15 SBP:12                             # Initial value for log precision parameter for random walk at global level (log lambda_g)
  
  #for each country will have the years for the three DBP,SBP,INTER
  u[1,]                                           <- rep(0, J*(T^2))                           # Initial values for component of nonlinear trend at national level
  v[1,]                                           <- rep(0, K*(T^2))                           # Initial values for component of nonlinear trend at regional level
  sv[1,]                                          <- rep(0, L*(T^2))                           # Initial values for component of nonlinear trend at superregional level
  w[1,]                                           <- rep(0, (T^2))                             # Initial values for component of nonlinear trend at global level
  
  ###### Overdisperse the starting values ########################################
  theta.max                                       <- 20                                                            # Constrain theta values for random walk to 20
  phi_s[1,]                                       <- rnorm(1:3, phi_s[1,], phi_s.prop.sd/2.4*sqrt(6))              # Log variance for normal prior for superregion random intercepts (log kappa_a^r)
  phi_r[1,]                                       <- rnorm(1:3, phi_r[1,], phi_r.prop.sd/2.4*sqrt(6))              # Log variance for normal prior for region random intercepts (log kappa_a^s)
  phi_c[1,]                                       <- rnorm(1:3, phi_c[1,], phi_c.prop.sd/2.4*sqrt(6))              # Log variance for normal prior for country random intercepts (log kappa_a^c)
  eta_s[1,]                                       <- rnorm(1:3, eta_s[1,], eta_s.prop.sd/2.4*sqrt(6))              # Log variance for normal prior for superregion random slopes (log kappa_b^r)
  eta_r[1,]                                       <- rnorm(1:3, eta_r[1,], eta_r.prop.sd/2.4*sqrt(6))              # Log variance for normal prior for region random slopes (log kappa_b^s)
  eta_c[1,]                                       <- rnorm(1:3, eta_c[1,], eta_c.prop.sd/2.4*sqrt(6))              # Log variance for normal prior for country random slopes (log kappa_b^c)
  phi_natl[1,] 	                                  <- rnorm(1:3, phi_natl[1,], phi_natl.prop.sd/2.4*sqrt(3))
  for (ss in 1:3) phi_subn[1,ss] 	                <- ur(pinvd.new(udnorm(mean=phi_subn[1,ss], sd=phi_subn.prop.sd/2.4*sqrt(3), lb=phi_natl[1,ss])),1)
  for (ss in 1:3) phi_comm[1,ss] 	                <- ur(pinvd.new(udnorm(mean=phi_comm[1,ss], sd=phi_comm.prop.sd/2.4*sqrt(3), lb=phi_subn[1,ss])),1)
  tau[1,]                                         <- rnorm(3, tau[1,], tau.prop.sd/2.4)                                # Log variance for within-study errors that differ between age groups (log tau)
  sigma2[1,]                                      <- exp(rnorm(3*5, log(sigma2[1,]), log.sigma2.prop.sd/2.4))           # Variances for country-specific random spline coefficients
  theta_c[1]                                      <- ur(pinvd.new(udnorm(mean=theta_c[1], sd= theta_c.prop.sd/2.4, ub=theta.max)),1)                           # Log precision parameter for random walk at country level (log lambda_c)
  theta_r[1]                                      <- ur(pinvd.new(udnorm(mean=theta_r[1], sd= theta_r.prop.sd/2.4, lb= theta_c[1], ub=theta.max)),1) 
  theta_s[1]                                      <- ur(pinvd.new(udnorm(mean=theta_s[1], sd= theta_s.prop.sd/2.4, lb= theta_r[1], ub=theta.max)),1)
  theta_g[1]                                      <- ur(pinvd.new(udnorm(mean=theta_g[1], sd= theta_g.prop.sd/2.4, lb= theta_s[1], ub=theta.max)),1) 
  
  ##### Diagonal of likelihood variance matrix ###################################
  #I will assume that the rho is the same. 
  # In the covariance matrix we have only take into account the covariance between the dependent variables
  s1                                              <- subset$se_dbp^2 + exp(tau[1,1])                 # Diagonal of variance matrix, including within-study errors that differ between age groups
  s2                                              <- subset$se_sbp^2 + exp(tau[1,2])                 # Diagonal of variance matrix, including within-study errors that differ between age groups
  rho                                             <- cor(y1,y2)
  s12                                             <- rep(rho,I)*subset$se_dbp*subset$se_sbp + exp(tau[1,3]) #covariance: s_(X,Y) = rho_(XY)*s_X*s_Y
  
  ##### alternative way of computing the rho. Rho is variable based on age, so we have 9 different values. ######
  
  # corl <- NULL;subset$rho <- NULL; #subset$age_gr   <- cut(subset$age,breaks = c(unique(subset$age)[1]-1,unique(subset$age)[2:8]) )
  # #for (i in 1:length(levels(subset$age))) corl[i] <- cor(subset[subset$age_gr==levels(subset$age_gr)[i],4],subset[subset$age_gr==levels(subset$age_gr)[i],5])
  # for (i in 1:length(unique(subset$age))) corl[i] <- cor(subset[subset$age==unique(subset$age)[i],4],subset[subset$age==unique(subset$age)[i],5])
  # for (i in 1:length(unique(subset$age)) ) subset$rho[subset$age==unique(subset$age)[i]]            <- abs(corl[i]) ### !
  # s12                                                                                               <- subset$rho*sqrt(s1*s2) + exp(tau[1,3])
  
  Sigma.diag                                      <- matrix(0,3*I, 3*I)                     # Diagonal of variance matrix, including within-study errors that differ between age groups
  #although the means and the variances are different depending on a sample of observations in each study,
  # we will assume that the correlation between DBP and SBP is the same for all the studies.
  Sigma.diag[1:I,1:I]                             <- diag(s1)
  Sigma.diag[(1+I):(2*I),(1+I):(2*I)]             <- diag(s2)
  Sigma.diag[(1+2*I):(3*I),(1+2*I):(3*I)]         <- diag(s12) 
  SigmaInv.diag                                   <- diag(diag(Sigma.diag)^(-1))
  #SigmaInv.diag                                  <- solve(Sigma.diag)                                  # Inverse of diagonal of variance matrix
  ##### Variance of study-specific random effects ################################################################################
  
  V.ssre                                          <- matrix(NA,N,3)          # Empty vector for study-specific random effect variances
  for (ii in 1:3) {#kk <- foreach(ii = 1:3, .combine=rbind) %do% {
    V.ssre[,ii][uid.match$coverage == "National"] <- exp(phi_natl[1,ii])       # Set variances for random effects for national studies
    V.ssre[,ii][uid.match$coverage == "Subnational"]  <- exp(phi_subn[1,ii])       # Set variances for random effects for subnational studies
    V.ssre[,ii][uid.match$coverage == "Community"] <- exp(phi_comm[1,ii])       # Set variances for random effects for community studies
  }
  ##### Diagonal of variance matrices for theta and gamma ########################################################################
  
  V.diag                                          <- c(                                                              # Diagonal of variance for theta matrix
    rep(1/epsilon,3),                                                                                                                             # Flat variance for global intercept DBP
    rep(1/epsilon,3),                                                                                                                             # Flat variance for global slope DBP
    rep(exp(phi_s[1,]), each = L),                                                                                                                # Variances for normal prior for superregion random intercepts DBP
    rep(exp(eta_s[1,]), each = L),                                                                                                                # Variances for normal prior for superregion random slopes DBP
    rep(exp(phi_r[1,]), each = sum(multipleRegionsInSregion)),                                                                                    # Variances for normal prior for region random intercepts DBP
    rep(exp(eta_r[1,]), each = sum(multipleRegionsInSregion)),                                                                                    # Variances for normal prior for region random slopes DBP
    rep(exp(phi_c[1,]), each = J),                                                                                                                # Variances for normal prior for country random intercepts DBBP
    rep(exp(eta_c[1,]), each = J),                                                                                                                # Variances for normal prior for country random slopes DBP
    V.ssre[,1],                                                                                                                                   # Variances for study-specific random effects
    V.ssre[,2],                                                                                                                                   # Variances for study-specific random effects
    V.ssre[,3]
  )
  V.diag                                          <- c(V.diag,rep(1/epsilon, p))                                     # Add flat variances for covariates
  VInv.diag                                       <- 1/V.diag                                                        # Inverse variance
  WInv.diag                                       <- c(rep(epsilon, 3*(5+5)), rep(1/sigma2[1,], each=J) )            # Inverse variance for age model
  ###### Initial values for the theta update #####################################
  
  ageMat3                                         <- matrix(0,3*I,3*dim(ageMat)[2])
  ageMat3[1:I,1:5]                                <- ageMat
  ageMat3[(1+I):(2*I),6:10]                       <- ageMat
  ageMat3[(1+2*I):(3*I),11:15]                    <- ageMat
  Phi.agePlus1                                    <- as.vector(1 + ageMat3 %*% gamma[1,16:30])                       # the first 4335 rows are for theDBP, the next 4335 for SBP etc
  FplusPhiF.age                                   <- F * Phi.agePlus1
  FplusPhiF.age.prime                             <- t(FplusPhiF.age) 
  Q.init                                          <- FplusPhiF.age.prime %*% (diag(SigmaInv.diag) * FplusPhiF.age)
  diag(Q.init)                                    <- diag(Q.init) + VInv.diag
  U.theta.init                                    <- base::chol(Q.init) #Q.init 
  ##################################################################################################
  #Constants of Age for DBP, SBP, INTER
  gamma_const                                     <- matrix(0,3*I,3*dim(ageMat)[2])
  gamma_const[1:I,1:5]                            <- rep(1,I) %*% t(gamma[1,1:5])
  gamma_const[(1+I):(2*I),6:10]                   <- rep(1,I) %*% t(gamma[1,6:10])
  gamma_const[(1+2*I):(3*I),11:15]                <- rep(1,I) %*% t(gamma[1,11:15])
  #StudyTocountry*Agemat
  gamma_countr                                    <- matrix(0,3*I,3*dim(ageMat)[2])
  gamma_countr[1:I,1:5]                           <- matrix(gamma[1, 31:(J*5+30)], J, 5)[which.country,]
  gamma_countr[(1+I):(2*I),6:10]                  <- matrix(gamma[1, (J*5+31):(2*J*5+30)], J, 5)[which.country,]
  gamma_countr[(1+2*I):(3*I),11:15]               <- matrix(gamma[1,(2*J*5+31):(3*J*5+30)], J, 5)[which.country,]
  varphiPlusC.age                                 <- rowSums((gamma_const + gamma_countr) * ageMat3)   # we want to fill the columns fist as its each one of the coefficients for all the countries
  
  ###### Initial values for the update of nonlinear terms ########################
  #I will create again 3 columns for  DBP,SBP,inter
  V1uInv.diag                                     <- Vinv.M.u <- rep(0,J*(T^2))#T      # Country level
  V1vInv.diag                                     <- Vinv.M.v <- rep(0,K*(T^2))#T      # Region level
  V1svInv.diag                                    <- Vinv.M.sv<- rep(0,L*(T^2))#T      # Superregion level
  V1wInv.diag                                     <- Vinv.M.w <- rep(0,(T^2))#T        # Global level
  
  ###### Initial values for the gamma update #####################################
  F.theta.Mm                                      <- y                                 # A vector for the mean [y1,y2]^T
  
  #The computations for the R.age would be easier if we had written y as only one column i.e all the
  #observations one after the other and each three banches would be refering to our group
  
  R                                               <- cbind(                            # Starting matrix mapping ages to gamma parameters
    ageMat,                                                                                                         # Age matrix
    ageMat,
    ageMat,   
    ageMat*F.theta.Mm[1:I],                                                                                         # Age matrix multiplied by mean
    ageMat*F.theta.Mm[(1+I):(2*I)],                                                                                 # Age matrix multiplied by mean
    ageMat*F.theta.Mm[(1+2*I):(3*I)],                                                                               # Age matrix multiplied by mean
    StudyToCountry*ageMat[,1],                                                                                      # Age matrix corresponding to sigma2_1
    StudyToCountry*ageMat[,2],                                                                                      # Age matrix corresponding to sigma2_2
    StudyToCountry*ageMat[,3],                                                                                      # Age matrix corresponding to sigma2_3
    StudyToCountry*ageMat[,4],                                                                                      # Age matrix corresponding to sigma2_4
    StudyToCountry*ageMat[,5],
    
    StudyToCountry*ageMat[,1],                                                                                      # Age matrix corresponding to sigma2_1
    StudyToCountry*ageMat[,2],                                                                                      # Age matrix corresponding to sigma2_2
    StudyToCountry*ageMat[,3],                                                                                      # Age matrix corresponding to sigma2_3
    StudyToCountry*ageMat[,4],                                                                                      # Age matrix corresponding to sigma2_4
    StudyToCountry*ageMat[,5],
    
    StudyToCountry*ageMat[,1],                                                                                      # Age matrix corresponding to sigma2_1
    StudyToCountry*ageMat[,2],                                                                                      # Age matrix corresponding to sigma2_2
    StudyToCountry*ageMat[,3],                                                                                      # Age matrix corresponding to sigma2_3
    StudyToCountry*ageMat[,4],                                                                                      # Age matrix corresponding to sigma2_4
    StudyToCountry*ageMat[,5])                                                                                      # Age matrix corresponding to sigma2_5
  
  R.age                                           <- rbind(R,R,R)
  R.age[1:I,c(6:15,21:30,(J*5+31):(3*J*5+30))]    <- 0 
  R.age[(1+I):(2*I),c(1:5,11:20,26:(J*5+30),(2*J*5+31):(3*J*5+30))]  <- 0 
  R.age[(1+2*I):(3*I),c(1:10,16:25,31:(2*J*5+30))]<- 0 
  R.age                                           <- Matrix(R.age)#, sparse=TRUE
  R.age.prime                                     <- t(R.age) 
  Q.init2                                         <- R.age.prime %*% (diag(SigmaInv.diag)*R.age)                                                                                   # Full conditional precision of gamma
  diag(Q.init2)                                   <- diag(Q.init2) + WInv.diag                                                                                                     # Full conditional precision of gamma
  U.gamma.init                                    <- base::chol(Q.init2)                                                                                                           # cholesky decomposition
  
  ###### Start chains at values consistent with the hyperparameter starting values
  ##### Theta ####################################################################
  
  update.theta                                    <- function(SigmaInv.diag, VInv.diag, u, v, sv, w) {
    Q                                             <- FplusPhiF.age.prime %*% (diag(SigmaInv.diag)* FplusPhiF.age)                                                                  # Full conditional precision of theta
    diag(Q)                                       <- diag(Q) + VInv.diag                                                                                                           # Full conditional precision of theta
    U                                             <- base::chol(Q)
    VinvM                                         <- FplusPhiF.age.prime %*% (diag(SigmaInv.diag) *(y - (u[rep(which.countryTime,3)] +
                                                                                                           v[rep(which.regionTime,3)] +sv[rep(which.sregionTime,3)] + w[rep(which.time,3)]) * Phi.agePlus1 - varphiPlusC.age))
    return(base::backsolve(U, base::forwardsolve(t(U), VinvM) + rnorm(length(VinvM))) )                                                                                # Full conditional precision times full conditional mean
  }
  ###########################
  theta[1,]                                       <- update.theta(SigmaInv.diag, VInv.diag, u[1,], v[1,], sv[1,], w[1,])# First update of theta 
  F.theta.Mm                                      <- as.vector(F%*% theta[1,]+u[1,rep(which.countryTime,3)]+v[1,rep(which.regionTime,3)]+sv[1,rep(which.sregionTime,3)]+w[1,rep(which.time,3)])
  R.age[1:I,16:20]                                <- ageMat*F.theta.Mm[1:I]                                                                                               # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age[(1 + I):(2*I),21:25]                      <- ageMat*F.theta.Mm[(1+I):(2*I)]                                                                                       # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age[(1 + 2*I):(3*I),26:30]                    <- ageMat*F.theta.Mm[(1+2*I):(3*I)]                                                                                     # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age.prime                                     <- t(R.age)                                                                                                             # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  
  ##### Non-linear component at national level ###################################
  #### Users should refer to pp5-7 and pp14-16 of Danaei et al for the line-by-  #
  #### line details of this section of code ######################################
  tPtsPerSregion                                  <- as.numeric(colSums(table(subset$data_year,subset$sregion)>0))                                                                 # Calculates number of years with data present for each superregion
  tPtsPerRegion                                   <- as.numeric(colSums(table(subset$data_year,subset$region)>0))                                                                  # Calculates number of years with data present for each region
  subset$country                                  <- as.character(subset$country)                                                                                                  # Calculates number of years with data present for each country
  country.match$name                              <- as.character(country.match$name)                                                                                              # Calculates number of years with data present for each country
  country.count                                   <- as.numeric(colSums(table(subset$data_year,as.character(subset$country) )>0))                                                               # Calculates number of years with data present for each country
  tPtsPerCountry                                  <- vector("integer", J)                                                                                                          # Calculates number of years with data present for each country
  tPtsPerCountry[as.numeric(colnames(country.table))] <- country.count[index(colnames(country.table))]                                                                                 # Calculates number of years with data present for each country
  #################
  # I will repeat u three times. The u_c, for i country will contain T values and it will be the same for DBP,SBP, INTER
  #since the u value will contain all the info for DBP,SBP,INTER together
  u.Prec.yMinusMean                               <- diag(SigmaInv.diag) * Phi.agePlus1 * (y - (F*Phi.agePlus1) %*% theta[1,] - (v[1,rep(which.regionTime,3)]+sv[1,rep(which.sregionTime,3)]+w[1,][rep(which.time,3)])*Phi.agePlus1-varphiPlusC.age)
  agg1                                            <- aggregate(as.numeric(u.Prec.yMinusMean), list(rep(which.countryTime,3)), sum)
  agg2                                            <- aggregate((diag(SigmaInv.diag) * Phi.agePlus1^2), list(rep(which.countryTime,3)), sum)
  Vinv.M.u[agg1[,1]]                              <- agg1[,2] 
  V1uInv.diag[agg2[,1]]                           <- agg2[,2]                       #for NO-DATA nl=1,2,3 should be the same 
  uStar                                           <- lapply(1:J, u.prop.function, tPtsPerCountry, V1uInv.diag, Vinv.M.u, theta_c[1], theta_c[1], u[1,],P,A,T,eigenNoData,SigmaGenInvNoTheta)
  u.star                                          <- rep(NA, J*(T^2))
  for(j in 1:J) u.star[((j-1)*(T^2)+1):(j*(T^2))] <- uStar[[j]]$u.star
  u[1,]                                           <- Re(u.star)
  u1                                              <- u[1,rep(which.countryTime,3)]
  v1                                              <- v[1,rep(which.regionTime,3)]
  sv1                                             <- sv[1,rep(which.sregionTime,3)]
  w1                                              <- w[1,rep(which.time,3)]
  F.theta.Mm                                      <- as.vector(F%*% theta[1,] + u1 + v1 + sv1 + w1)                                                                                # New estimates of BMI at age 50
  R.age.star                                      <- R.age  
  R.age.star[1:I,16:20]                           <- ageMat*F.theta.Mm[1:I]                                                                                           # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age.star[(1 + I):(2*I),21:25]                 <- ageMat*F.theta.Mm[(1+I):(2*I)]                                                                      # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age.star[(1 + 2*I):(3*I),26:30]               <- ageMat*F.theta.Mm[(1+2*I):(3*I)]                                                                    # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  ##### Non-linear component at regional level ###################################
  ##### See Danaei et al pp5-7 and 14-16 for details #############################
  v.Prec.yMinusMean                               <- diag(SigmaInv.diag) * Phi.agePlus1 * (y - (F*Phi.agePlus1) %*% theta[1,] - (u1 + sv1 + w1) * Phi.agePlus1 - varphiPlusC.age)
  agg3                                            <- aggregate(as.numeric(v.Prec.yMinusMean), list(rep(which.regionTime,3)),sum)
  agg4                                            <- aggregate((diag(SigmaInv.diag) * Phi.agePlus1^2),list(rep(which.regionTime,3)),sum)
  Vinv.M.v[agg3[,1]]                              <- agg3[,2]
  V1vInv.diag[agg4[,1]]                           <- agg4[,2]
  u.star                                          <- rep(NA, K*(T^2))
  u.star[rep(!multipleRegionsInSregion, each=T^2)]<- 0
  uStar                                           <- lapply(1:K, u.prop.function, tPtsPerRegion, V1vInv.diag, Vinv.M.v, theta_r[1],theta_r[1],v[1,],P,A,T,eigenNoData,SigmaGenInvNoTheta)
  for(j in 1:K)  u.star[((j-1)*(T^2)+1):(j*(T^2))]<- uStar[[j]]$u.star
  v[1,]                                           <- Re(u.star)
  v1                                              <- v[1,rep(which.regionTime,3)]
  F.theta.Mm                                      <- as.vector(F%*% theta[1,] + u1 + v1 + sv1 + w1)                                                                                    # New estimates of BMI at age 50
  R.age.star                                      <- R.age 
  R.age.star[1:I,16:20]                           <- ageMat*F.theta.Mm[1:I]                                                                                           # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age.star[(1 + I):(2*I),21:25]                 <- ageMat*F.theta.Mm[(1+I):(2*I)]                                                                      # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age.star[(1 + 2*I):(3*I),26:30]               <- ageMat*F.theta.Mm[(1+2*I):(3*I)]                                                                    # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  ##### Non-linear component at superregional level ##############################
  ##### See Danaei et al pp5-7 and 14-16 for details #############################
  sv.Prec.yMinusMean                              <- diag(SigmaInv.diag) * Phi.agePlus1 * (y - (F*Phi.agePlus1) %*% theta[1,] - (u1 + v1 + w1 ) * Phi.agePlus1 - varphiPlusC.age)
  agg5                                            <- aggregate(as.numeric(sv.Prec.yMinusMean),list(rep(which.sregionTime,3)), sum)
  agg6                                            <- aggregate((diag(SigmaInv.diag) * Phi.agePlus1^2),list(rep(which.sregionTime,3)),sum)
  Vinv.M.sv[agg5[,1]]                             <- agg5[,2]
  V1svInv.diag[agg6[,1]]                          <- agg6[,2]
  u.star                                          <- rep(NA,L*(T^2))
  uStar                                           <- lapply(1:L, u.prop.function, tPtsPerSregion, V1svInv.diag, Vinv.M.sv, theta_s[1], theta_s[1], sv[1,],P,A,T,eigenNoData,SigmaGenInvNoTheta)
  for(j in 1:L) u.star[((j-1)*(T^2)+1):(j*(T^2))] <- uStar[[j]]$u.star
  sv[1,]                                          <- Re(u.star)
  sv1                                             <- sv[1,rep(which.sregionTime,3)]
  F.theta.Mm                                      <- as.vector(F%*% theta[1,] + u1 + v1 + sv1 + w1)                                                                                        # New estimates of BMI at age 50
  R.age.star                                      <- R.age
  R.age.star[1:I,16:20]                           <- ageMat*F.theta.Mm[1:I]                                                                                                                # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age.star[(1 + I):(2*I),21:25]                 <- ageMat*F.theta.Mm[(1+I):(2*I)]                                                                                                        # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age.star[(1 + 2*I):(3*I),26:30]               <- ageMat*F.theta.Mm[(1+2*I):(3*I)]                                                                                                      # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  ##### Non-linear component at global level #####################################
  ##### See Danaei et al pp5-7 and 14-16 for details #############################
  w.Prec.yMinusMean                               <- diag(SigmaInv.diag) * Phi.agePlus1 * (y - (F*Phi.agePlus1) %*% theta[1,] - (u1  + sv1  + v1) *Phi.agePlus1 - varphiPlusC.age)
  agg7                                            <- aggregate(as.numeric(w.Prec.yMinusMean), list(rep(which.time,3)),sum)
  agg8                                            <- aggregate((diag(SigmaInv.diag) * Phi.agePlus1^2), list(rep(which.time,3)),sum)
  Vinv.M.w[agg7[,1]]                              <- agg7[,2]
  V1wInv.diag[agg8[,1]]                           <- agg8[,2]
  u.star                                          <- rep(NA,(T^2))
  uStar                                           <- lapply(1, u.prop.function, 3, V1wInv.diag, Vinv.M.w, theta_g[1], theta_g[1], w[1,],P,A,T,eigenNoData,SigmaGenInvNoTheta)
  u.star                                          <- uStar[[1]]$u.star
  w[1,]                                           <- as.vector(Re(u.star))
  w1                                              <- w[1,rep(which.time,3)]                                                                                                                      #update
  F.theta.Mm                                      <- as.vector(F%*% theta[1,] + u1 + v1 + sv1 + w1)                                                                                        # New estimates of BMI at age 50
  R.age.star                                      <- R.age
  R.age.star[1:I,16:20]                           <- ageMat*F.theta.Mm[1:I]                                                                                                                # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age.star[(1 + I):(2*I),21:25]                 <- ageMat*F.theta.Mm[(1+I):(2*I)]                                                                                                        # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  R.age.star[(1 + 2*I):(3*I),26:30]               <- ageMat*F.theta.Mm[(1+2*I):(3*I)]                                                                                                      # Update gamma mapping matrix based on new mean BMI at age 50 estimates
  ##### Age model ################################################################
  Q                                               <- R.age.prime %*% (diag(SigmaInv.diag) * R.age)                                                                                         # Full conditional precision of gamma
  diag(Q)                                         <- diag(Q) + WInv.diag                                                                                                                   # Full conditional precision of gamma
  U                                               <- base::chol(Q)
  VinvM                                           <- as.matrix(R.age.prime) %*% (as.matrix(diag(SigmaInv.diag) * (y - F.theta.Mm)))                                                         # Full conditional precision x full conditional mean
  gamma[1,]                                       <- base::backsolve(U, base::forwardsolve(t(U),VinvM) + rnorm(length(VinvM)))                        # Returns updated gamma following Gibbs step
  Phi.agePlus1                                    <- as.vector(1 + ageMat3 %*% gamma[1,16:30])                                                                                              # Update derived vector for future calculations
  FplusPhiF.age                                   <- F * Phi.agePlus1                                                                                                                       # Update derived matrix for future calculations
  FplusPhiF.age.prime                             <- t(FplusPhiF.age)                                                                                                                       # Update derived matrix for future calculations
  #Constants of Age for DBP, SBP, INTER
  gamma_const                                     <- matrix(0,3*I,3*dim(ageMat)[2])
  gamma_const[1:I,1:5]                            <- rep(1,I) %*% t(gamma[1,1:5])
  gamma_const[(1+I):(2*I),6:10]                   <- rep(1,I) %*% t(gamma[1,6:10])
  gamma_const[(1+2*I):(3*I),11:15]                <- rep(1,I) %*% t(gamma[1,11:15])
  #StudyTocountry*Agemat
  gamma_countr                                    <- matrix(0,3*I,3*dim(ageMat)[2])
  gamma_countr[1:I,1:5]                           <- matrix(gamma[1,        31:(J*5+30)],   J, 5)[which.country,]
  gamma_countr[(1+I):(2*I),6:10]                  <- matrix(gamma[1,  (J*5+31):(2*J*5+30)], J, 5)[which.country,]
  gamma_countr[(1+2*I):(3*I),11:15]               <- matrix(gamma[1,(2*J*5+31):(3*J*5+30)], J, 5)[which.country,]
  varphiPlusC.age                                 <- rowSums((gamma_const + gamma_countr) * ageMat3)   # we want to fill the columns fist as its each one of the coefficients for all the countries
  ###### Set current value holders ###############################################
  
  phi_c.current                                   <- phi_c[1,]                                             # Current value holder for log variance for normal prior for country random intercepts (log kappa_a^c)
  phi_r.current                                   <- phi_r[1,]                                             # Current value holder for log variance for normal prior for region random intercepts (log kappa_a^s)
  phi_s.current                                   <- phi_s[1,]                                             # Current value holder for log variance for normal prior for superregion random intercepts (log kappa_a^r)
  eta_c.current                                   <- eta_c[1,]                                             # Current value holder for log variance for normal prior for country random slopes (log kappa_b^c)
  eta_r.current                                   <- eta_r[1,]                                             # Current value holder for log variance for normal prior for region random slopes (log kappa_b^s)
  eta_s.current                                   <- eta_s[1,]                                             # Current value holder for log variance for normal prior for superregion random slopes (log kappa_b^r)
  phi_natl.current                                <- phi_natl[1,]                                          # Current value holder for log variance of random effects for national studies (log nu_"national")
  phi_subn.current                                <- phi_subn[1,]                                          # Current value holder for log variance of random effects for subnational studies (log nu_s)
  phi_comm.current                                <- phi_comm[1,]                                          # Current value holder for log variance of random effects for community studies (log nu_c)
  tau.current                                     <- tau[1,]                                               # Current value holder for log variance for within-study errors that differ between age groups (log tau)
  theta_c.current                                 <- theta_c[1]                                           # Current value holder for log precision parameter for random walk at country level (log lambda_c)
  theta_r.current                                 <- theta_r[1]                                           # Current value holder for log precision parameter for random walk at region level (log lambda_s)
  theta_s.current                                 <- theta_s[1]                                           # Current value holder for log precision parameter for random walk at superregion level (log lambda_r)
  theta_g.current                                 <- theta_g[1]                                           # Current value holder for log precision parameter for random walk at global level (log lambda_g)
  gamma.current                                   <- gamma[1,]                                             # Current value holder for age model parameters
  sigma2.current                                  <- sigma2[1,]                                            # Current value holder for variances for country-specific random spline coefficients
  theta.current                                   <- theta[1,]                                             # Current value holder for theta matrix
  u1.current                                      <- u1                                                    # Current value holder for component of nonlinear trend at national level
  v1.current                                      <- v1                                                    # Current value holder for component of nonlinear trend at regional level
  sv1.current                                     <- sv1                                                   # Current value holder for component of nonlinear trend at superregional level
  w1.current                                      <- w1 
  u.current                                       <- u[1,]                                                 # Current value holder for component of nonlinear trend at national level
  v.current                                       <- v[1,]                                                 # Current value holder for component of nonlinear trend at regional level
  sv.current                                      <- sv[1,]                                                # Current value holder for component of nonlinear trend at superregional level
  w.current                                       <- w[1,] 
  deviance[1]                                     <- LogLik(SigmaInv.diag, F.theta.Mm, gamma[1,], R.age,y) * -2 + I * log(2*pi)  # Calculate the initial deviance
  
  ################################################################################
  ##### ***START OF MCMC LOOP*** #################################################
  ################################################################################
  
  # Initialise acceptance counts for Metropolis-Hastings steps
  
  N_comm                                          <- sum(uid.match$coverage=="Community")                      # Number of studies with community coverage
  N_subn                                          <- sum(uid.match$coverage=="Subnational")                    # Number of studies with subnational coverage
  N_natl                                          <- sum(uid.match$coverage=="National")                       # Number of studies with national coverage
  
  # Initialise acceptance counts for Metropolis-Hastings steps
  V.acc                                           <- 0                                                         # Set the acceptance count for the linear random intercepts and slopes to zero
  phi.acc                                         <- rep(0,3) 
  tau.acc                                         <- rep(0,3)                                                  # Set the acceptance count for the log variance for the within-study errors that differ between age groups to zero
  sigma2.acc                                      <- rep(0,15)                                                 # Set the acceptance counts for the variances of the country-specific random spline coefficients to zero
  theta_c.acc                                     <- 0                                                         # Set the acceptance count for the log precision parameter for random walk at country level to zero
  theta_r.acc                                     <- 0                                                         # Set the acceptance count for the log precision parameter for random walk at region level zero
  theta_s.acc                                     <- 0                                                         # Set the acceptance count for the log precision parameter for random walk at superregion level to zero
  theta_g.acc                                     <- 0                                                         # Set the acceptance count for the log precision parameter for random walk at global level to zero
  phi_s.sd.star <- phi_r.sd.star <- phi_c.sd.star <- rep(NA,3)
  eta_s.sd.star <- eta_r.sd.star <- eta_c.sd.star <- rep(NA,3)
  theta_c.star<-theta_r.star<-theta_s.star<-theta_g.star  <- NULL
  start_time1                                     <- Sys.time()
  for(i in 2:nLong ) {                                                                # Loop from 2nd to 55000th MCMC iteration                                                        # Print the number of iterations to screen (every 100 iterations)
    #if(i%%500==0) tracePlots()                                                      # Save traceplots to PDF
    #if(i%%10000==0 | i==nLong) {burnt <- 501:(i/10); printDeviance() }                                       # Print deviance to screen
    start_time2                                   <- Sys.time()
    #if (i%%freq.val==0 & i< 1000){#3000
    #  ww                                        <- c(ww,V.acc)
    #  phi_s.prop.sd                             <- unlist(lapply(1:3, function(tt) phi_s.prop.sd[tt] * adaptJump(6,V.acc/freq.val)))             # Tune the proposal SD for log variance for normal prior for superregion random intercepts (log kappa_a^r)
    #  eta_s.prop.sd                             <- unlist(lapply(1:3, function(tt) eta_s.prop.sd[tt] * adaptJump(6,V.acc/freq.val)))             # Tune the proposal SD for log variance for normal prior for superregion random slopes (log kappa_b^r)
    #  phi_r.prop.sd                             <- unlist(lapply(1:3, function(tt) phi_r.prop.sd[tt] * adaptJump(6,V.acc/freq.val)))             # Tune the proposal SD for log variance for normal prior for region random intercepts (log kappa_a^s)
    #  eta_r.prop.sd                             <- unlist(lapply(1:3, function(tt) eta_r.prop.sd[tt] * adaptJump(6,V.acc/freq.val)))             # Tune the proposal SD for log variance for normal prior for region random slopes (log kappa_b^s)
    #  phi_c.prop.sd                             <- unlist(lapply(1:3, function(tt) phi_c.prop.sd[tt] * adaptJump(6,V.acc/freq.val)))             # Tune the proposal SD for log variance for normal prior for country random intercepts (log kappa_a^c)
    #  eta_c.prop.sd                             <- unlist(lapply(1:3, function(tt) eta_c.prop.sd[tt] * adaptJump(6,V.acc/freq.val)))             # Tune the proposal SD for log variance for normal prior for country random slopes (log kappa_b^c)
    #  theta_c.prop.sd                           <- theta_c.prop.sd * adaptJump(1,theta_c.acc/freq.val)                                          # Tune the proposal SD for log precision parameter for random walk at country level (log lambda_c)
    #  theta_r.prop.sd                           <- theta_r.prop.sd * adaptJump(1,theta_r.acc/freq.val)                                          # Tune the proposal SD for log precision parameter for random walk at region level (log lambda_s)
    #  theta_s.prop.sd                           <- theta_s.prop.sd * adaptJump(1,theta_s.acc/freq.val)                                          # Tune the proposal SD for log precision parameter for random walk at superregion level (log lambda_r)
    #  theta_g.prop.sd                           <- theta_g.prop.sd * adaptJump(1,theta_g.acc/freq.val)                                          # Tune the proposal SD for log precision parameter for random walk at global level (log lambda_g)
    #  tau.prop.sd                               <- unlist(lapply(1:3, function(tt) tau.prop.sd[tt]      * adaptJump(1,tau.acc[tt]/freq.val)))        # Tune the proposal SD for log variance for within-study errors that differ between age groups (log tau)
    #  phi_natl.prop.sd                          <- unlist(lapply(1:3, function(tt) phi_natl.prop.sd[tt] * adaptJump(3,phi.acc[tt]/freq.val)))   # Tune the proposal SD for log variance of random effects for national studies (log nu_"national")
    #  phi_subn.prop.sd                          <- unlist(lapply(1:3, function(tt) phi_subn.prop.sd[tt] * adaptJump(3,phi.acc[tt]/freq.val)))   # Tune the proposal SD for log variance of random effects for subnational studies (log nu_s)
    #  phi_comm.prop.sd                          <- unlist(lapply(1:3, function(tt) phi_comm.prop.sd[tt] * adaptJump(3,phi.acc[tt]/freq.val)))   # Tune the proposal SD for log variance of random effects for community studies (log nu_c)
    #  log.sigma2.prop.sd                        <- unlist(lapply(1:15,function(tt)log.sigma2.prop.sd[tt]* adaptJump(1,sigma2.acc[tt]/freq.val))) # Tune proposal variances for log variances for country-specific random spline coefficients
    #  phi.acc                                   <- rep(0,3)                                                    # Set the acceptance count for the study-specific random effects to zero
    #  tau.acc                                   <- rep(0,3)                                                    # Set the acceptance count for the within-study errors that differ between age groups to zero
    #  theta_c.acc                               <- 0                                                           # Set the acceptance count for the log precision parameter for random walk at country level to zero
    #  theta_r.acc                               <- 0                                                           # Set the acceptance count for the log precision parameter for random walk at region level to zero
    #  theta_s.acc                               <- 0                                                           # Set the acceptance count for the log precision parameter for random walk at superregion level to zero
    #  theta_g.acc                               <- 0                                                           # Set the acceptance count for the log precision parameter for random walk at global level to zero
    #  V.acc                                     <- 0                                                           # Set the acceptance count for the linear random intercepts and slopes to zero
    #  sigma2.acc                                <- rep(0, 15)
    #}
    ssre                                          <- matrix(NA,N,3)
    ssre[,1]                                      <- theta.current[(3+3*L+3*sum(multipleRegionsInSregion)+3*J+1):(3+3*L+3*sum(multipleRegionsInSregion)+3*J+N)]        # we count only the uid vcariables of the F matrix
    ssre[,2]                                      <- theta.current[(3+3*L+3*sum(multipleRegionsInSregion)+3*J+1+N):(3+3*L+3*sum(multipleRegionsInSregion)+3*J+2*N)]    # we count only the uid vcariables of the F matrix
    ssre[,3]                                      <- theta.current[(3+3*L+3*sum(multipleRegionsInSregion)+3*J+1+2*N):(3+3*L+3*sum(multipleRegionsInSregion)+3*J+3*N)]  # we count only the uid vcariables of the F matrix
    phi_natl.star                                 <- phi_subn.star <- phi_comm.star <- NULL
    for (k in 1:3){
      phi_natl.star[k]                            <- rnorm(1, phi_natl.current[k], phi_natl.prop.sd[k])                                     # Propose new value for log variance of national study-specific random effects
      phi_subn.star[k]                            <- rnorm(1, phi_subn.current[k], phi_subn.prop.sd[k])                                     # Propose new value for log variance of subnational study-specific random effects
      phi_comm.star[k]                            <- rnorm(1, phi_comm.current[k], phi_comm.prop.sd[k])                                     # Propose new value for log variance of community study-specific random effects
      
      if(phi_natl.star[k] < phi_subn.star[k] & phi_subn.star[k] < phi_comm.star[k]) {#length(which())!=0                                                 # Constraint: national level has smaller variance than subnational and in turn community
        R1                                        <- exp( 
          ssreLik(N_natl, phi_natl.star[k], ssre[,k][uid.match$coverage=="National"]) +  ##No                                                      # Log likelihood of proposed state for national level
            LogPriorPhi(phi_natl.star[k]) -                                                                                                    # Log prior of proposed state for national level
            ssreLik(N_natl, phi_natl.current[k], ssre[,k][uid.match$coverage=="National"]) -                                                      # Log likelihood of existing state for national level
            LogPriorPhi(phi_natl.current[k]) +                                                                                                    # Log prior of existing state for national level
            ssreLik(N_subn, phi_subn.star[k], ssre[,k][uid.match$coverage=="Subnational"]) + #No,                                                  # Log likelihood of proposed state for subnational level
            LogPriorPhi(phi_subn.star[k]) -                                                                                                    # Log prior of proposed state for subnational level
            ssreLik(N_subn, phi_subn.current[k], ssre[,k][uid.match$coverage=="Subnational"]) -                                                   # Log likelihood of existing state for subnational level
            LogPriorPhi(phi_subn.current[k]) +                                                                                                    # Log prior of existing state for subnational level
            ssreLik(N_comm, phi_comm.star[k], ssre[,k][uid.match$coverage=="Community"]) +                                                     # Log likelihood of proposed state for community level
            LogPriorPhi(phi_comm.star[k]) -                                                                                                    # Log prior of proposed state for community level
            ssreLik(N_comm, phi_comm.current[k], ssre[,k][uid.match$coverage=="Community"]) -                                                     # Log likelihood of existing state for community level
            LogPriorPhi(phi_comm.current[k]) )                                                                                                    # Log prior of existing state for community level
        if(runif(1) < R1[1]) {                                                                                                                              # Metropolis-Hastings update: accept proposed state with probability R
          phi_natl.current[k]                     <- phi_natl.star[k]                                                                     # If accepted, set log variance of random effects for national studies to * value (log nu_natl)
          phi_subn.current[k]                     <- phi_subn.star[k]                                                                     # If accepted, set log variance of random effects for subnational studies to * value (log nu_s)
          phi_comm.current[k]                     <- phi_comm.star[k]                                                                     # If accepted, set log variance of random effects for community studies to * value (log nu_c)
          phi.acc[k]                              <- phi.acc[k] + 1                                                                          # If accepted, increment the acceptance count for the study-specific random effects
          V.ssre[,k][uid.match$coverage == "National"]          <- exp(phi_natl.current[k])                                                                # If accepted, set variance of random effects for national studies to * value
          V.ssre[,k][uid.match$coverage == "Subnational"]       <- exp(phi_subn.current[k])             
          V.ssre[,k][uid.match$coverage == "Community"]         <- exp(phi_comm.current[k])             
          V.diag[(2+2*L+2*sum(multipleRegionsInSregion)+2*J+1+(k-1)*N):(2+2*L+2*sum(multipleRegionsInSregion)+2*J+k*N)]     <- V.ssre[,k]                   # If accepted, update overall variance diagonal vector
          VInv.diag                               <- 1/V.diag                                                                                # If accepted, update overall inverse variance diagonal vector
        }                                                                                                                                                   # Close Metropolis-Hastings update
      }
    }                                                                                                                                                       # Close Metropolis-Hastings update
    #### UPDATE LINEAR RANDOM INTERCEPTS AND SLOPES VARIANCES AND THETA MATRIX #
    #### See pp5 and 13-14 of Danaei et al #####################################
    ############################################################################
    phi_c.sd.star                                 <- unlist(lapply(1:3, function(g) rnorm(1, sqrt(exp(phi_c.current[g])),phi_c.prop.sd[g]) ))
    phi_r.sd.star                                 <- unlist(lapply(1:3, function(g) rnorm(1, sqrt(exp(phi_r.current[g])),phi_r.prop.sd[g]) ))
    phi_s.sd.star                                 <- unlist(lapply(1:3, function(g) rnorm(1, sqrt(exp(phi_s.current[g])),phi_s.prop.sd[g]) ))
    eta_c.sd.star                                 <- unlist(lapply(1:3, function(g) rnorm(1, sqrt(exp(eta_c.current[g])),eta_c.prop.sd[g]) ))
    eta_r.sd.star                                 <- unlist(lapply(1:3, function(g) rnorm(1, sqrt(exp(eta_r.current[g])),eta_r.prop.sd[g]) ))
    eta_s.sd.star                                 <- unlist(lapply(1:3, function(g) rnorm(1, sqrt(exp(eta_s.current[g])),eta_s.prop.sd[g]) ))
    ##### All proposed SD values must be positive ##############################
    if(all(phi_s.sd.star > 0 & eta_s.sd.star > 0 & phi_r.sd.star > 0 & eta_r.sd.star > 0 & phi_c.sd.star > 0 & eta_c.sd.star > 0)) {
      V.diag.star                                 <- c(1/epsilon,1/epsilon,1/epsilon,1/epsilon,1/epsilon,1/epsilon,                                                                                       
                                                       rep(phi_s.sd.star^2, each = L),
                                                       rep(eta_s.sd.star^2, each = L),                                                                         
                                                       rep(phi_r.sd.star^2, each = sum(multipleRegionsInSregion)),
                                                       rep(eta_r.sd.star^2, each = sum(multipleRegionsInSregion)),                                             
                                                       rep(phi_c.sd.star^2, each = J),
                                                       rep(eta_c.sd.star^2, each = J),
                                                       V.ssre[,1],V.ssre[,2],V.ssre[,3],rep(1/epsilon, p) )
      VInv.diag.star                              <- 1/V.diag.star
      ##### Gibbs step to propose new theta based on the new variances #######
      Q1 <- Q.star                                <- FplusPhiF.age.prime %*% (diag(SigmaInv.diag) * FplusPhiF.age)
      VinvM                                       <- FplusPhiF.age.prime %*% as.matrix(diag(SigmaInv.diag) *(as.vector(y)- (u1.current + v1.current + sv1.current + w1.current) 
                                                                                                             * Phi.agePlus1 - varphiPlusC.age)) # Full conditional precision x full conditional mean
      diag(Q.star)                                <- diag(Q.star) + VInv.diag.star
      diag(Q1)                                    <- diag(Q1)     + VInv.diag
      U.star                                      <- base::chol(Q.star)
      U1                                          <- base::chol(Q1)
      M.star                                      <- base::backsolve(U.star,base::forwardsolve(t(U.star),VinvM))
      M                                           <- base::backsolve(U1,base::forwardsolve(t(U1),VinvM))
      theta.star                                  <- M.star + base::backsolve(U.star, rnorm(length(VinvM)))
      F.theta.star.Mm                             <- as.vector(F%*%theta.star + u1.current + v1.current + sv1.current + w1.current) 
      R.age.star                                  <- R.age
      R.age.star[1:I,16:20]                       <- ageMat*F.theta.Mm[1:I]                                                                                           # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star[(1+I):(2*I),21:25]               <- ageMat*F.theta.Mm[(1+I):(2*I)]                                                                      # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star[(1+2*I):(3*I),26:30]             <- ageMat*F.theta.Mm[(1+2*I):(3*I)]                                                                    # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star.prime                            <- t(R.age.star)
      # Metropolis-Hastings acceptance ratio
      (a1                                         <- .5*(sum(log(diag(SigmaInv.diag))) - t(as.vector(y) - F.theta.star.Mm- R.age.star %*% gamma.current)%*%as.matrix(diag(SigmaInv.diag)*(as.vector(y)- F.theta.star.Mm- R.age.star %*% gamma.current)) ) )
      (a2                                         <- .5*(sum(log(VInv.diag.star)) - sum(VInv.diag.star * theta.star^2)) )
      (a3                                         <- - sum(log(diag(U.star))) +.5 * t(theta.star - M.star) %*% Q.star %*% (theta.star - M.star))
      (b1                                         <- .5*(sum(log(diag(SigmaInv.diag))) - t(as.vector(y)- F.theta.Mm- R.age %*% gamma.current)%*%as.matrix(diag(SigmaInv.diag)*(as.vector(y)-F.theta.Mm- R.age %*% gamma.current)) ) )
      (b2                                         <- .5*(sum(log(VInv.diag)) - sum(VInv.diag * theta.current^2)) )
      (b3                                         <- - sum(log(diag(U1))) +.5 * t(theta.current - M) %*% Q1 %*% (theta.current - M))  
      (R2                                         <- exp(a1-b1+a2-b2+a3-b3))
      ##### Accept/reject ###
      if( runif(1) < R2[1] ) {                                                                
        phi_s.current                             <- log(phi_s.sd.star^2)
        eta_s.current                             <- log(eta_s.sd.star^2)
        phi_r.current                             <- log(phi_r.sd.star^2)
        eta_r.current                             <- log(eta_r.sd.star^2)
        phi_c.current                             <- log(phi_c.sd.star^2)
        eta_c.current                             <- log(eta_c.sd.star^2)
        theta.current                             <- theta.star
        VInv.diag                                 <- VInv.diag.star
        V.diag                                    <- V.diag.star
        V.acc                                     <- V.acc + 1
        F.theta.Mm                                <- as.vector(F%*% theta.current + u1.current + v1.current + sv1.current + w1.current)   # If accepted, update
        R.age[1:I,16:20]                          <- ageMat*F.theta.Mm[1:I]                                                                                           # Update gamma mapping matrix based on new mean BMI at age 50 estimates
        R.age[(1 + I):(2*I),21:25]                <- ageMat*F.theta.Mm[(1+I):(2*I)]                                                                      # Update gamma mapping matrix based on new mean BMI at age 50 estimates
        R.age[(1 + 2*I):(3*I),26:30]              <- ageMat*F.theta.Mm[(1+2*I):(3*I)]                                                                    # Update gamma mapping matrix based on new mean BMI at age 50 estimates
        R.age.prime                               <- t(R.age)
      }}
    #### NON-LINEAR TRENDS #####################################################
    #### Users should refer to pp5-7 and pp14-16 of Danaei et al for the line- #
    #### by-line details of this section of code ###############################
    ############################################################################
    # Component of nonlinear trend at national level #
    theta_c.star                                  <- rnorm(1, theta_c.current, theta_c.prop.sd)                              # Propose theta_c*
    if(theta_c.star < theta_r.current & theta_c.star < theta.max) {                         # Only enter update step if constraints are satisfied
      u.Prec.yMinusMean                           <- diag(SigmaInv.diag) * Phi.agePlus1 * (y - (F*Phi.agePlus1) %*% theta.current - (v1.current + sv1.current + w1.current) * Phi.agePlus1 - varphiPlusC.age)
      agg1                                        <- aggregate(as.numeric(u.Prec.yMinusMean), list(rep(which.countryTime,3)), sum)
      agg2                                        <- aggregate((diag(SigmaInv.diag) * Phi.agePlus1^2),list(rep(which.countryTime,3)),sum)
      Vinv.M.u[agg1[,1]]                          <- agg1[,2] 
      V1uInv.diag[agg2[,1]]                       <- agg2[,2]
      u.star1                                     <- rep(NA, J*(T^2)); u.dens.old <- u.dens.star <- rep(NA, J)
      uStar                                       <- lapply(1:J, u.prop.function, tPtsPerCountry, V1uInv.diag, Vinv.M.u, theta_c.star, theta_c.current, u.current,P,A,T,eigenNoData,SigmaGenInvNoTheta)
      for(j in 1:J) u.star1[((j-1)*(T^2)+1):(j*(T^2))]   <- uStar[[j]]$u.star
      u.dens.old                                  <- unlist(lapply(1:J, function(j) uStar[[j]]$dens.old))
      u.dens.star                                 <- unlist(lapply(1:J, function(j) uStar[[j]]$dens.star))
      u.star2                                     <- Re(u.star1); u.dens.old <- Re(u.dens.old) ; u.dens.star <- Re(u.dens.star)
      u.star                                      <- u.star2[rep(which.countryTime,3)]
      F.theta.Mm.star                             <- as.vector(F%*% theta.current + u.star + v1.current + sv1.current + w1.current)
      R.age.star                                  <- R.age
      R.age.star[1:I, 16:20]                      <- ageMat*F.theta.Mm.star[1:I]                                                                                           # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star[(1 + I):(2*I),21:25]             <- ageMat*F.theta.Mm.star[(1+I):(2*I)]                                                                      # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star[(1 + 2*I):(3*I),26:30]           <- ageMat*F.theta.Mm.star[(1+2*I):(3*I)]                                                                    # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R3                                          <- exp( LogLik(SigmaInv.diag, F.theta.Mm.star, gamma.current, R.age.star,y) - 
                                                            LogLik(SigmaInv.diag, F.theta.Mm, gamma.current, R.age,y) +  
                                                            LogPostThetaC(theta_c.star, u.star2,J,T,P) -
                                                            LogPostThetaC(theta_c.current, u.current,J,T,P)+ 
                                                            sum(u.dens.old) - sum(u.dens.star))
      if(runif(1) < R3[1]) {
        theta_c.current                           <- theta_c.star
        u1.current                                <- u.star
        u.current                                 <- u.star2   
        F.theta.Mm                                <- F.theta.Mm.star
        R.age                                     <- R.age.star
        R.age.prime                               <- t(R.age.star)
        theta_c.acc                               <- theta_c.acc + 1
      }}
    # Component of nonlinear trend at regional level #
    theta_r.star                                  <- rnorm(1, theta_r.current, theta_r.prop.sd)
    if(theta_r.star > theta_c.current & theta_r.star < theta_s.current & theta_r.star < theta.max) {
      v.Prec.yMinusMean                           <- diag(SigmaInv.diag) * Phi.agePlus1 * (y - (F*Phi.agePlus1) %*% theta.current - (u1.current + sv1.current + w1.current) * Phi.agePlus1 - varphiPlusC.age)
      agg3                                        <- aggregate(as.numeric(v.Prec.yMinusMean), list(rep(which.regionTime,3)),sum)
      agg4                                        <- aggregate((diag(SigmaInv.diag) * Phi.agePlus1^2), list(rep(which.regionTime,3)),sum)
      Vinv.M.v[agg3[,1]]                          <- agg3[,2]
      V1vInv.diag[agg4[,1]]                       <- agg4[,2]
      u.star1                                     <- rep(NA, K*(T^2)); u.dens.old <- u.dens.star <- rep(NA, K) 
      uStar                                       <- lapply(1:K, u.prop.function, tPtsPerRegion, V1vInv.diag, Vinv.M.v, theta_r.star, theta_r.current, v.current,P,A,T,eigenNoData,SigmaGenInvNoTheta) 
      for(j in 1:K) u.star1[((j-1)*(T^2)+1):(j*(T^2))] <- uStar[[j]]$u.star
      u.dens.old                                  <- unlist(lapply(1:K, function(j) uStar[[j]]$dens.old ))
      u.dens.star                                 <- unlist(lapply(1:K, function(j) uStar[[j]]$dens.star ))
      # Set to zero in those regions that are also super-regions
      u.star1[rep(!multipleRegionsInSregion, each=T^2)]   <- 0
      u.dens.star[!multipleRegionsInSregion]      <- 0
      u.dens.old[!multipleRegionsInSregion]       <- 0
      v.star2                                     <- Re(u.star1) ; u.dens.old  <- Re(u.dens.old) ; u.dens.star <- Re(u.dens.star)
      u.star                                      <- v.star2[rep(which.regionTime, 3)]
      F.theta.Mm.star                             <- as.vector(F%*% theta.current + u1.current + u.star + sv1.current + w1.current)
      R.age.star                                  <- R.age
      R.age.star[1:I, 16:20]                      <- ageMat*F.theta.Mm.star[1:I]                                                                                           # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star[(1 + I):(2*I),21:25]             <- ageMat*F.theta.Mm.star[(1+I):(2*I)]                                                                      # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star[(1 + 2*I):(3*I),26:30]           <- ageMat*F.theta.Mm.star[(1+2*I):(3*I)]                                                                    # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R4                                          <- exp(LogLik(SigmaInv.diag, F.theta.Mm.star, gamma.current, R.age.star,y) -
                                                           LogLik(SigmaInv.diag, F.theta.Mm, gamma.current, R.age,y) +
                                                           LogPostThetaR(theta_r.star, v.star2, multipleRegionsInSregion,K,T,P) - 
                                                           LogPostThetaR(theta_r.current, v.current, multipleRegionsInSregion,K,T,P) +
                                                           sum(u.dens.old) - sum(u.dens.star) )
      if(runif(1) < R4[1]) {
        theta_r.current                           <- theta_r.star
        v1.current                                <- u.star
        v.current                                 <- v.star2
        F.theta.Mm                                <- F.theta.Mm.star
        R.age                                     <- R.age.star
        R.age.prime                               <- t(R.age.star)
        theta_r.acc                               <- theta_r.acc + 1
      }}
    #Component of nonlinear trend at superregional level #
    theta_s.star                                  <- rnorm(1, theta_s.current, theta_s.prop.sd)
    if(theta_s.star > theta_r.current & theta_s.star < theta_g.current & theta_s.star< theta.max) {
      sv.Prec.yMinusMean                          <- diag(SigmaInv.diag )* Phi.agePlus1 * (y - (F*Phi.agePlus1) %*% theta.current- (u1.current + v1.current + w1.current) * Phi.agePlus1 - varphiPlusC.age)
      agg5                                        <- aggregate(as.numeric(sv.Prec.yMinusMean),list(rep(which.sregionTime,3)), sum)
      agg6                                        <- aggregate((diag(SigmaInv.diag) * Phi.agePlus1^2), list(rep(which.sregionTime,3)),sum)
      Vinv.M.sv[agg5[,1]]                         <- agg5[,2]
      V1svInv.diag[agg6[,1]]                      <- agg6[,2]
      u.star1                                     <- rep(NA, L*(T^2)); u.dens.old <- u.dens.star <- rep(NA, L)
      uStar                                       <- lapply(1:L, u.prop.function, tPtsPerSregion, V1svInv.diag, Vinv.M.sv, theta_s.star, theta_s.current, sv.current,P,A,T,eigenNoData,SigmaGenInvNoTheta)
      for(j in 1:L) u.star1[((j-1)*(T^2)+1):(j*(T^2))] <- uStar[[j]]$u.star
      u.dens.old                                  <- unlist(lapply(1:L, function(j) uStar[[j]]$dens.old ))
      u.dens.star                                 <- unlist(lapply(1:L, function(j) uStar[[j]]$dens.star ))
      sv.star2                                    <- Re(u.star1) ; u.dens.old <- Re(u.dens.old) ; u.dens.star <- Re(u.dens.star)
      u.star                                      <- sv.star2[rep(which.sregionTime,3)]
      F.theta.Mm.star                             <- as.vector(F%*% theta.current+ u1.current + v1.current + u.star + w1.current)
      R.age.star                                  <- R.age
      R.age.star[1:I, 16:20]                      <- ageMat*F.theta.Mm.star[1:I]                                                                                           # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star[(1 + I):(2*I),21:25]             <- ageMat*F.theta.Mm.star[(1+I):(2*I)]                                                                      # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star[(1 + 2*I):(3*I),26:30]           <- ageMat*F.theta.Mm.star[(1+2*I):(3*I)]                                                                    # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R5                                          <- exp(LogLik(SigmaInv.diag, F.theta.Mm.star, gamma.current, R.age.star,y) -
                                                           LogLik(SigmaInv.diag, F.theta.Mm, gamma.current, R.age,y) +
                                                           LogPostThetaS(theta_s.star,    sv.star2,L,T,P) -
                                                           LogPostThetaS(theta_s.current, sv.current,L,T,P) +
                                                           sum(u.dens.old) - sum(u.dens.star) )
      if(runif(1) < R5[1]) {
        theta_s.current                           <- theta_s.star
        sv1.current                               <- u.star
        sv.current                                <- sv.star2
        F.theta.Mm                                <- F.theta.Mm.star
        R.age                                     <- R.age.star
        R.age.prime                               <- t(R.age.star)
        theta_s.acc                               <- theta_s.acc + 1
      }}
    # Component of nonlinear trend at global level #
    theta_g.star                                  <- rnorm(1, theta_g.current, theta_g.prop.sd)
    if(theta_g.star > theta_s.current & theta_g.star < theta.max) {
      w.Prec.yMinusMean                           <- diag(SigmaInv.diag) * Phi.agePlus1 * (y - (F*Phi.agePlus1) %*% theta.current - (u1.current + v1.current + sv1.current) * Phi.agePlus1 - varphiPlusC.age)
      agg7                                        <- aggregate(as.numeric(w.Prec.yMinusMean), list(rep(which.time,3)),sum)
      agg8                                        <- aggregate((diag(SigmaInv.diag) * Phi.agePlus1^2), list(rep(which.time,3)),sum)
      Vinv.M.w[agg7[,1]]                          <- agg7[,2]
      V1wInv.diag[agg8[,1]]                       <- agg8[,2]
      # for(j in unique(which.time)) {
      #   Vinv.M.w[j]                             <- sum(w.Prec.yMinusMean[gg][which.time==j])
      #   V1wInv.diag[j]                          <- diag(SigmaInv.diag)[gg][which.time==j] %*% (Phi.agePlus1[gg][which.time==j])^2
      # }
      uStar                                       <- lapply(1, u.prop.function, 3, V1wInv.diag, Vinv.M.w, theta_g.star, theta_g.current, w.current,P,A,T,eigenNoData,SigmaGenInvNoTheta)
      u.star1                                     <- uStar[[1]]$u.star
      u.dens.old                                  <- uStar[[1]]$dens.old
      u.dens.star                                 <- uStar[[1]]$dens.star
      w.star2                                     <- as.vector(Re(u.star1)); u.dens.old <- Re(u.dens.old) ; u.dens.star <- Re(u.dens.star)
      u.star                                      <- w.star2[rep(which.time,3)]
      F.theta.Mm.star                             <- as.vector(F%*% theta.current + u1.current + v1.current + sv1.current + u.star)
      R.age.star                                  <- R.age
      R.age.star[1:I, 16:20]                      <- ageMat*F.theta.Mm.star[1:I]                                                                                           # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star[(1 + I):(2*I),21:25]             <- ageMat*F.theta.Mm.star[(1+I):(2*I)]                                                                      # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R.age.star[(1 + 2*I):(3*I),26:30]           <- ageMat*F.theta.Mm.star[(1+2*I):(3*I)]                                                                    # Update gamma mapping matrix based on new mean BMI at age 50 estimates
      R6                                          <- exp(LogLik(SigmaInv.diag, F.theta.Mm.star, gamma.current, R.age.star,y) - 
                                                           LogLik(SigmaInv.diag, F.theta.Mm, gamma.current, R.age,y) +
                                                           LogPostThetaG(theta_g.star, w.star2,T,P) - 
                                                           LogPostThetaG(theta_g.current, w.current,T,P) +
                                                           sum(u.dens.old) - sum(u.dens.star) )
      if(runif(1) < R6[1]) {
        theta_g.current                           <- theta_g.star
        w1.current                                <- u.star
        w.current                                 <- w.star2
        F.theta.Mm                                <- F.theta.Mm.star
        R.age                                     <- R.age.star
        R.age.prime                               <- t(R.age.star)
        theta_g.acc                               <- theta_g.acc+ 1
      }}
    ##### GAMMA: update using Gibbs ############################################
    ############################################################################
    Q                                             <- R.age.prime %*% (diag(SigmaInv.diag) * R.age)                                                                # Full conditional precision of gamma
    diag(Q)                                       <- diag(Q) + WInv.diag                                                                                              # Full conditional precision of gamma
    U                                             <- base::chol(Q)#Matrix(base:: base::chol(Q, sparse=TRUE))                                                                       # Update base::cholesky decomposition
    VinvM                                         <- as.matrix(R.age.prime) %*% (as.matrix(diag(SigmaInv.diag) * (y - F.theta.Mm)))                                               # Full conditional precision x conditional mean
    gamma.current                                 <- base::backsolve(U, base::forwardsolve(t(U), VinvM) + rnorm(length(VinvM))) # Updated gamma after Gibbs step
    Phi.agePlus1                                  <- as.vector(1 + ageMat3 %*% gamma.current[16:30])                                                                    # Update derived vector for future calculations
    FplusPhiF.age                                 <- F * Phi.agePlus1                                                                                            # Update derived matrix for future calculations
    FplusPhiF.age.prime                           <- t(FplusPhiF.age)                                                                                                 # Update derived matrix for future calculations
    #Constants of Age for DBP, SBP, INTER
    gamma_const[1:I,1:5]                          <- rep(1,I) %*% t(gamma.current[1:5])
    gamma_const[(1+I):(2*I),6:10]                 <- rep(1,I) %*% t(gamma.current[6:10])
    gamma_const[(1+2*I):(3*I),11:15]              <- rep(1,I) %*% t(gamma.current[11:15])
    #StudyTocountry*Agemat
    gamma_countr[1:I,1:5]                         <- matrix(gamma.current[31:(J*5+30)],   J, 5)[which.country,]
    gamma_countr[(1+I):(2*I),6:10]                <- matrix(gamma.current[(J*5+31):(2*J*5+30)], J, 5)[which.country,]
    gamma_countr[(1+2*I):(3*I),11:15]             <- matrix(gamma.current[(2*J*5+31):(3*J*5+30)], J, 5)[which.country,]
    varphiPlusC.age                               <- rowSums((gamma_const + gamma_countr) * ageMat3)                  # we want to fill the columns fist as its each one of the coefficients for all the countries
    
    ##### SIGMA2: update using Metropolis-Hastings #############################
    ############################################################################
    #for DBP
    for(a in 1:15) {                                                    # 5variances for each of our variables                                     # Loop over the five parameters
      log.sigma2.star                             <- rnorm(1, log(sigma2.current[a]), log.sigma2.prop.sd[a])                                 # Propose sigma2* values on log scale
      R7                                          <- exp( -.5*(J*log.sigma2.star + 1/exp(log.sigma2.star) * sum(gamma.current[30+(((a-1)*J+1):(a*J))]^2)) +                                                  # Log likelihood of proposed state
                                                            LogPriorPhi(log.sigma2.star) +                                                                                                            # Log prior of proposed state
                                                            .5*(J*log(sigma2.current[a]) + 1/sigma2.current[a] * sum(gamma.current[30+(((a-1)*J+1):(a*J))]^2)) -                                      # Log likelihood of existing state
                                                            LogPriorPhi(log(sigma2.current[a])) )                                                                                                      # Log prior of existing state # End of calculation of posterior probability
      if(runif(1) < R7[1]) {                                                                                                                             # Metropolis-Hastings update: accept proposed state with probability R
        sigma2.current[a]                         <- exp(log.sigma2.star)                                                                    # If accepted, set sigma2 to sigma2*
        WInv.diag[10+(((a-1)*J+1):(a*J))]         <- rep(1/exp(log.sigma2.star), J)                                                          # If accepted, set variance for sigma2 to variance for sigma2*
        sigma2.acc[a]                             <- sigma2.acc[a]+1                                                                         # If accepted, increment the acceptance count for sigma2
      } }                                                                                                                                             # Close Metropolis-Hastings update                                                                                                                                                 # Close loop over the five parameters
    #### TAU: update using Metropolis-Hastings ################################
    ###########################################################################
    tau.star                                      <- unlist(lapply(1:3, function(l) rnorm(1, tau.current[l], tau.prop.sd[l]) ))              # Propose tau* value on log scale
    s1                                            <- subset$se_dbp^2  + exp(tau.star[1])                                                     # Diagonal of variance matrix, including within-study errors that differ between age groups
    s2                                            <- subset$se_sbp^2  + exp(tau.star[2])                                                     # Diagonal of variance matrix, including within-study errors that differ between age groups
    s12                                           <- rep(rho,I)*subset$se_dbp*subset$se_sbp + exp(tau.star[3])
    Sigma.diag.star                               <- matrix(0,3*I,3*I)
    Sigma.diag.star[1:I,1:I]                      <- diag(s1)
    Sigma.diag.star[(1+I):(2*I),(1+I):(2*I)]      <- diag(s2)
    Sigma.diag.star[(1+2*I):(3*I),(1+2*I):(3*I)]  <- diag(s12)
    SigmaInv.diag.star                            <- diag(diag(Sigma.diag.star)^(-1))                                 # Inverse of diagonal of variance matrix
    for ( l in 1:3){# foreach(l = 1:3, .combine=rbind) %do% {
      R8                                          <- exp(                                                             # Posterior probability R for proposed state
        LogLik(SigmaInv.diag.star, F.theta.Mm, gamma.current, R.age,y) +                                                            # Log likelihood of proposed state
          LogPriorPhi(tau.star[l]) -                                                                                                # Log prior of proposed state
          LogLik(SigmaInv.diag, F.theta.Mm, gamma.current, R.age,y) -                                                                 # Log likelihood of existing state
          LogPriorPhi(tau.current[l]))   # End of calculation of posterior probability
      if(runif(1) < R8[1]) {                                                                                                                              # Metropolis-Hastings update: accept proposed state with probability R
        tau.current[l]                            <- tau.star[l]                                                    # If accepted, set tau to tau star
        tau.acc[l]                                <- tau.acc[l] + 1                                                 # If accepted, increment the acceptance count for tau
        Sigma.diag                                <- Sigma.diag.star                                                # If accepted, set variance to variance*
        SigmaInv.diag                             <- SigmaInv.diag.star                                             # If accepted, set inverse-variance to inverse-variance*
      } }                                                                                                                                                 # Close Metropolis-Hastings update
    ##### Save state every 10th iteration ######################################
    ############################################################################
    if(i%%10==0) {        #10                                                                                                                           # Thinning is carried out here, with parameters only saved every 10 iterations
      j                                           <- i/10                                                                                          # Save iteration in corresponding position/row of vectors/matrices
      phi_natl[j,]                                <- phi_natl.current
      phi_subn[j,]                                <- phi_subn.current
      phi_comm[j,]                                <- phi_comm.current
      phi_s[j,]                                   <- phi_s.current
      phi_r[j,]                                   <- phi_r.current
      phi_c[j,]                                   <- phi_c.current
      eta_s[j,]                                   <- eta_s.current
      eta_r[j,]                                   <- eta_r.current
      eta_c[j,]                                   <- eta_c.current
      theta[j,]                                   <- theta.current
      tau[j,]                                     <- tau.current
      gamma[j,]                                   <- gamma.current
      sigma2[j,]                                  <- sigma2.current
      theta_c[j]                                  <- theta_c.current
      theta_r[j]                                  <- theta_r.current
      theta_s[j]                                  <- theta_s.current
      theta_g[j]                                  <- theta_g.current
      u[j,]                                       <- u.current             
      v[j,]                                       <- v.current        
      sv[j,]                                      <- sv.current   
      w[j,]                                       <- w.current       
      deviance[j]                                 <- LogLik(SigmaInv.diag,F.theta.Mm, gamma.current, R.age,y) * -2 + I * log(2*pi) # Save deviance
    }
    end_time2                                     <- Sys.time()
    print(i)
    print(end_time2-start_time2)
  }
  len                                             <- (nLong/10 - 3000):(nLong/10)
  return(list(deviance=deviance[len],R.age=R.age, F.theta.Mm=F.theta.Mm, phi_natl=phi_natl[len,], 
              phi_subn=phi_subn[len,],phi_comm=phi_comm[len,], phi_s=phi_s[len,], phi_r=phi_r[len,], 
              phi_c=phi_c[len,],eta_s=eta_s[len,],eta_r=eta_r[len,], eta_c=eta_c[len,], theta=theta[len,], 
              tau=tau[len,], gamma=gamma[len,], sigma2=sigma2[len,], theta_c=theta_c[len],theta_r=theta_r[len],
              theta_s=theta_s[len],theta_g=theta_g[len], V.acc=V.acc, phi.acc=phi.acc,
              tau.acc=tau.acc, sigma2.acc=sigma2.acc, theta_c.acc=theta_c.acc,
              theta_r.acc=theta_r.acc, theta_s.acc=theta_s.acc, theta_g.acc=theta_g.acc,
              u=u[len,], v=v[len,], sv=sv[len,], w=w[len,], phi_s.prop.sd=phi_s.prop.sd, 
              eta_s.prop.sd=eta_s.prop.sd, phi_r.prop.sd=phi_r.prop.sd, 
              eta_r.prop.sd=eta_r.prop.sd, phi_c.prop.sd=phi_c.prop.sd, eta_c.prop.sd=eta_c.prop.sd, 
              theta_c.prop.sd=theta_c.prop.sd, theta_r.prop.sd=theta_r.prop.sd,
              theta_s.prop.sd=theta_s.prop.sd, theta_g.prop.sd=theta_g.prop.sd))
  end_time1                                       <- Sys.time() 
  print(end_time1-start_time1)
}
########   Parallel Computing   ###########################################
cl                 <- makeCluster(10) #(detectCores() - 1)
par.setup          <- parLapply( cl, 1:length(cl),function(xx) {c(library(Matrix),library(Runuran),library(zoo),
                                                                  library(MASS),library(reshape2),library(stringr),library(splines))} )
clusterExport( cl, c("MCMC","u.prop.function","LogLik","ssreLik","LogPriorPhi",
                     "CalcThetaPost","LogPostThetaC","LogPostThetaR","LogPostThetaS","LogPostThetaG","adaptJump",
                     "printDeviance","tracePlots"))#"TenYearWeightedAvg",
start_time2        <- Sys.time()
chains             <- parLapply(cl, 1:length(cl), function(xx) MCMC() )
stopCluster(cl)
Sys.time()- start_time2


############# tracePlots() #####
# len                                             <- 9000:10000
# colnames(chains[[10]]$R.age) <-rownames(chains[[10]]$R.age)   <- NULL
# rnd                                             <- rnorm(1,0,20)                                               
# write.csv(c(V.acc=chains[[10]]$V.acc,phi.acc=chains[[10]]$phi.acc,tau.acc=chains[[10]]$tau.acc,sigma2.acc=chains[[10]]$sigma2.acc,
#             theta_c.acc=chains[[10]]$theta_c.acc,theta_r.acc=chains[[10]]$theta_r.acc,theta_s.acc=chains[[10]]$theta_s.acc,
#             theta_g.acc=chains[[10]]$theta_g.acc), file = paste(round(rnd,2),"acceptance1.csv"))
# write.csv(chains[[10]]$deviance[len],  file = paste(round(rnd,2),"deviance1.csv"));write.csv(chains[[10]]$phi_natl[len,], file = paste(round(rnd,2),"phi_natl1.csv"))
# write.csv(chains[[10]]$phi_subn[len,], file = paste(round(rnd,2),"phi_subn1.csv"));write.csv(chains[[10]]$phi_comm[len,], file = paste(round(rnd,2),"phi_comm1.csv"))
# write.csv(as.matrix(chains[[10]]$R.age), file = paste(round(rnd,2),"R.age1.csv")) ;write.csv(chains[[10]]$F.theta.Mm,  file = paste(round(rnd,2),"F.theta.Mm1.csv"))
# write.csv(chains[[10]]$phi_s[len,], file = paste(round(rnd,2),"phi_s1.csv"))      ;write.csv(chains[[10]]$phi_r[len,], file = paste(round(rnd,2),"phi_r1.csv"))
# write.csv(chains[[10]]$phi_c[len,], file = paste(round(rnd,2),"phi_c1.csv"))      ;write.csv(chains[[10]]$eta_s[len,], file = paste(round(rnd,2),"eta_s1.csv")) 
# write.csv(chains[[10]]$eta_r[len,], file = paste(round(rnd,2),"eta_r1.csv"))      ;write.csv(chains[[10]]$eta_c[len,], file = paste(round(rnd,2),"eta_c1.csv")) 
# write.csv(chains[[10]]$theta[len,], file = paste(round(rnd,2),"theta1.csv"))      ;write.csv(chains[[10]]$tau[len,], file = paste(round(rnd,2),"tau1.csv"))
# write.csv(chains[[10]]$gamma[len,], file = paste(round(rnd,2),"gamma1.csv"))      ;write.csv(chains[[10]]$sigma2[len,], file = paste(round(rnd,2),"sigma21.csv"))
# write.csv(chains[[10]]$theta_c[len], file = paste(round(rnd,2),"theta_c1.csv"))   ;write.csv(chains[[10]]$theta_r[len], file = paste(round(rnd,2),"theta_r1.csv"))
# write.csv(chains[[10]]$theta_s[len], file = paste(round(rnd,2),"theta_s1.csv"))   ;write.csv(chains[[10]]$theta_g[len], file = paste(round(rnd,2),"theta_g1.csv"))
# write.csv(chains[[10]]$u[len,], file = paste(round(rnd,2),"u1.csv"))              ;write.csv(chains[[10]]$v[len,], file = paste(round(rnd,2),"v1.csv"))
# write.csv(chains[[10]]$sv[len,], file = paste(round(rnd,2),"sv1.csv"))            ;write.csv(chains[[10]]$w[len,], file = paste(round(rnd,2),"w1.csv"))

