library(doBy)
library(RColorBrewer)
library(ggplot2)
library(maptools)
library(plyr)
library(gpclib);gpclibPermit()##
library(scales)
require(gridExtra)
library(grid)
library(rgdal)
library(dplyr)
library(GiNA)
library(pdp)
library(astrochron)
library(splines)

start.year                                                    <- 2004#1990  #1975, 1951  
end.year                                                      <- 2014#2015  #2016                                                         
centring.year                                                 <- 2008#2002
sex.pop           	                                         	<- substr(sex,1,1)
filename                                                    	<- paste("Model54_mean",sex,sep="_")
dbp.cent.val                                                  <- dbp.centring.value                                                                   
sbp.cent.val                                                  <- sbp.centring.value                                                                  
int.cent.val                                                  <- int.centring.value
cent                                                          <- c(dbp.cent.val,sbp.cent.val,int.cent.val)
burnt                                                         <- seq(1,2001,1)
title                                                         <- c("Diastolic","Systolic","Interaction")

#ch                                                           <- 6

############### Baseline matrix of country-year estimates for 50 y.o.'s #########
theta1                                                        <- as.matrix(theta[burnt,]) #as.matrix(chains[[ch]]$theta[burnt,])#matrix(unlist(read.csv('C:\\Users\\mariz\\Downloads\\c_theta1.csv')),1500,dim(F)[2])
u1                                                            <- as.matrix(u[burnt,]) #as.matrix(chains[[ch]]$u[burnt,])
v1                                                            <- as.matrix(v[burnt,]) #[as.matrix(chains[[ch]]$v[burnt,])
sv1                                                           <- as.matrix(sv[burnt,]) #as.matrix(chains[[ch]]$sv[burnt,])
w1                                                            <- as.matrix(w[burnt,]) #as.matrix(chains[[ch]]$w[burnt,])
baseline                                                      <- matrix(NA, length(burnt), 3*(T^2)*J) #3*T*J
for(k in 1:J){
  country.val 				                                        <- as.character(country.match$name[k])
  country.val.num 	                 	                        <- country.match$number[country.match$name==country.val]
  region.val 				                                          <- as.character(unique(covar$region[covar$country==country.val]))
  region.val.num 		                	                        <- region.match$number[region.match$name==region.val]
  sregion.val 			                                        	<- unique(covar$sregion[covar$country==country.val])
  sregion.val.num 		                	                      <- sregion.match$number[sregion.match$name==sregion.val]
  #sregion
  sregion_pred 	<- sregion.time_pred 	                        <- matrix(0, L, T^2)
  sregion_pred[sregion.val.num,]                              <- 1
  sregion.time_pred[sregion.val.num,]                         <- t[,3]#t
  #region
  region_pred 	<- region.time_pred               	          <- matrix(0, K, T^2)
  region_pred[region.val.num,] 		                            <- 1
  region.time_pred[region.val.num,] 	                        <- t[,3]#t
  region_pred				                                          <- region_pred[multipleRegionsInSregion,]			
  region.time_pred	                        		              <- region.time_pred[multipleRegionsInSregion,]	#16th row won't be into account as it has been counting in regions
  #country
  country_pred 	<- country.time_pred 	                        <- matrix(0, J, T^2)
  country_pred[country.val.num,] 	                         	  <- 1
  country.time_pred[country.val.num,]     	                  <- t[,3]#t
  F_pr                                                        <- rbind(matrix(1,3,T^2), matrix(t[,3],3,T^2,byrow = T), 
                                                                       matrix(sregion_pred,3*L,T^2),rbind(sregion.time_pred,sregion.time_pred,sregion.time_pred), 
                                                                       matrix(region_pred,3*sum(multipleRegionsInSregion),T^2),rbind(region.time_pred,region.time_pred,region.time_pred),
                                                                       matrix(country_pred,3*J,T^2),rbind(country.time_pred,country.time_pred,country.time_pred),
                                                                       matrix(0, 3*N, T^2), 
                                                                       matrix(0, 3*4, T^2),             		
                                                                       t(matrix(covar[,covariate.list[1]][covar$country==country.val],T^2,3)),   	
                                                                       t(matrix(covar[,covariate.list[1]][covar$country==country.val]*t,T^2,3)),   	
                                                                       matrix(0, 3*4, T^2), matrix(0, 3*4, T^2)) 
  F_pred                                                      <- cbind(F_pr,F_pr,F_pr)
  all                                                         <- 6*(1+L+(sum(multipleRegionsInSregion))+J)
  #DBP
  F_pred[c(1,4,7:(6+L),(7+3*L):(6+4*L),(7+6*L):(6+6*L+sum(multipleRegionsInSregion)),(7+6*L+3*sum(multipleRegionsInSregion)):(6+6*L+4*(sum(multipleRegionsInSregion))),
           (7+6*L+6*(sum(multipleRegionsInSregion))):(6+6*L+6*(sum(multipleRegionsInSregion))+J),(7+6*L+6*(sum(multipleRegionsInSregion))+3*J):(6+6*L+6*(sum(multipleRegionsInSregion))+4*J),  
           (all+1):(all+dim(StudyToUid)[2]),all+1+3*dim(StudyToUid)[2]+seq(0,3*13,3)),(dim(F_pr)[2]+1):(3*dim(F_pr)[2])]  <- 0
  #all+1+seq(0,3*13,3)),(dim(F_pr)[2]+1):(3*dim(F_pr)[2])]  <- 0
  #SBP
  F_pred[c(2,5,(7+L):(6+2*L),(7+4*L):(6+5*L),(7+6*L+(sum(multipleRegionsInSregion))):(6+6*L+2*(sum(multipleRegionsInSregion))),
           (7+6*L+4*(sum(multipleRegionsInSregion))):(6+6*L+5*(sum(multipleRegionsInSregion))),(7+6*L+6*(sum(multipleRegionsInSregion))+J):(6+6*L+6*(sum(multipleRegionsInSregion))+2*J),(7+6*L+6*(sum(multipleRegionsInSregion))+4*J):(6+6*L+6*(sum(multipleRegionsInSregion))+5*J),
           (all+1+dim(StudyToUid)[2]):(all+2*dim(StudyToUid)[2]),all+2+3*dim(StudyToUid)[2]+seq(0,3*13,3)),c(1:dim(F_pr)[2],(2*dim(F_pr)[2]+1):(3*dim(F_pr)[2]))] <- 0
  #all+2+seq(0,3*13,3)),c(1:dim(F_pr)[2],(2*dim(F_pr)[2]+1):(3*dim(F_pr)[2]))] <- 0
  #INTER 
  F_pred[c(3,6,(7+2*L):(6+3*L),(7+5*L):(6+6*L),(7+6*L+2*(sum(multipleRegionsInSregion))):(6+6*L+3*(sum(multipleRegionsInSregion))),(7+6*L+5*(sum(multipleRegionsInSregion))):(6+6*L+6*(sum(multipleRegionsInSregion))),
           (7+6*L+6*(sum(multipleRegionsInSregion))+2*J):(6+6*L+6*(sum(multipleRegionsInSregion))+3*J),(7+6*L+6*(sum(multipleRegionsInSregion))+5*J):(6+6*L+6*(sum(multipleRegionsInSregion))+6*J),
           (all+1+2*dim(StudyToUid)[2]):(all+3*dim(StudyToUid)[2]),all+3+3*dim(StudyToUid)[2]+seq(0,3*13,3)),1:(2*dim(F_pr)[2])] <- 0
  #all+3+seq(0,3*13,3)),1:(2*dim(F_pr)[2])] <- 0
  uu                                                           <- rep(which(country.time.match$country == country.val),3) 
  vv                                                           <- rep(which(region.time.match$region == region.val),3) 
  ssv                                                          <- rep(which(sregion.time.match$sregion == sregion.val),3)
  
  baseline[,((k-1)*(3*(T^2))+1):(k*(3*(T^2)))]                 <- theta1 %*% F_pred + u1[,uu] + v1[,vv] + sv1[,ssv] + w1[,rep(1:(T^2),3)]
  # baseline[,((k-1)*3*T+1):(k*(T^2))]                         <- theta1 %*% F_pred + u[burnt, country.time.match$country == country.val] + v[burnt, region.time.match$region == region.val]
  #                                                                                    + sv[burnt, sregion.time.match$sregion == sregion.val] + w[burnt,] 
}
############### Age standardisation ###############################
#the columns inside the baseline are D,S,INT without being all the D,S,INT together etc.
g   <- 1 ;   colD <- colINT <- colS                            <- matrix(0,J*(T^2),1)
for (jj in seq(1,J*(T^2),(T^2))){
  colD[jj:(jj+(T^2)-1)]                                        <- ((g-1)*(T^2)+1):(g*(T^2))         # c(1:T^2,(3*T^2+1):(4*T^2),(6*T+1):(7*T))
  colS[jj:(jj+(T^2)-1)]                                        <- ((g)*(T^2)+1):((g+1)*(T^2))       # c((T^2+1):2*T^2,(4*T+1):(5*T),(7*T+1):(8*T))
  colINT[jj:(jj+(T^2)-1)]                                      <- ((g+1)*(T^2)+1):((g+2)*(T^2))     # c((2*T^2+1):3*T^2,(5*T+1):(6*T),(8*T+1):(9*T))
  g                                                            <- g + 3
}
gamma1                                                         <- as.matrix(gamma[burnt,])#,matrix(unlist(read.csv('C:\\Users\\mariz\\Downloads\\c_gamma1.csv')),1500,3030)
StandardPop 			                                             <- read.csv("C:\\Users\\Mariza\\Desktop\\2D_model_Results_19-03-22\\Data\\StandardPop_WHO.GBD_20plus.bmi.csv") 
age.specific.mean<-Dage.specific.mean <-Sage.specific.mean <-INTage.specific.mean <- list()
for(i in 1:nrow(StandardPop)){ 
  Dage.specific.mean[[i]]                                      <- matrix(NA, length(burnt), J*(T^2)) 
  Sage.specific.mean[[i]]                                      <- matrix(NA, length(burnt), J*(T^2))
  INTage.specific.mean[[i]]                                    <- matrix(NA, length(burnt), J*(T^2))
  age.specific.mean[[i]]                                       <- matrix(NA, length(burnt), 3*J*(T^2)) 
}       
age.vals 		                       	                           <- StandardPop$age_mean - middle.age
ageMat.val                                                     <- bs(age.vals,df = 5)[,1:5]#bs(age.vals,df = 25)[,1:25]
AgeMat31                                                       <- matrix(0,3*dim(ageMat.val)[1],3*dim(ageMat.val)[2])
AgeMat31[1:dim(ageMat.val)[1],1:5]                             <- ageMat.val
AgeMat31[(1+dim(ageMat.val)[1]):(2*dim(ageMat.val)[1]),6:10]   <- ageMat.val
AgeMat31[(1+2*dim(ageMat.val)[1]):(3*dim(ageMat.val)[1]),11:15]<- ageMat.val

for(j in 1:length(age.vals)){
  age.val                                                      <- age.vals[j]
  for(i in 1:length(burnt)){#studies
    c.part                                                     <- rep(NA, 3*J*(T^2))#rep(NA, 3*J*T)
    for(k in 1:J){
      c.part[((k-1)*(T^2)+1):(k*(T^2))]                        <- gamma1[i,  30+k+J*(0:4)] %*% AgeMat31[j,1:5]
      c.part[(J*(T^2)+(k-1)*(T^2)+1):((T^2)*(J+k))]            <- gamma1[i,J*5+30+k+J*(0:4)] %*% AgeMat31[length(age.vals)+j,6:10]
      c.part[(2*J*(T^2)+(k-1)*(T^2)+1):((T^2)*(2*J+k))]        <- gamma1[i,2*J*5+30+k+J*(0:4)] %*% AgeMat31[2*length(age.vals)+j,11:15]
    }
    Dage.specific.mean[[j]][i,]                                <- baseline[i, colD]  + as.vector(AgeMat31[j,1:5] %*% gamma1[i,1:5])  + baseline[i,colD]*as.vector(AgeMat31[j,1:5] %*% gamma1[i,16:20]) + c.part[1:(J*T)]
    Sage.specific.mean[[j]][i,]                                <- baseline[i,colS]   + as.vector(AgeMat31[length(age.vals)+j,6:10] %*% gamma1[i,6:10]) + baseline[i,colS]*as.vector(AgeMat31[length(age.vals)+j,6:10] %*% gamma1[i,21:25]) + c.part[(J*T+1):(2*J*T)]
    INTage.specific.mean[[j]][i,]                              <- baseline[i,colINT] + as.vector(AgeMat31[2*length(age.vals)+j,11:15] %*% gamma1[i,11:15])+ baseline[i,colINT]*as.vector(AgeMat31[2*length(age.vals)+j,11:15] %*% gamma1[i,26:30]) + c.part[(2*J*T+1):(3*J*T)] 
  }}
for (i in 1:nrow(StandardPop) ) age.specific.mean[[i]]         <- cbind(Dage.specific.mean[[i]],Sage.specific.mean[[i]],INTage.specific.mean[[i]])

##### Generate draws #############################################################################
cent.val                                                       <- rep(c(dbp.cent.val,sbp.cent.val,int.cent.val),each=J*(T^2))
for(i in 1:nrow(StandardPop)) {
  ## draws 
  draws.name  	                                               <- paste("draws",i,sep='')
  assign(draws.name, data.frame(rep(country.time.match$country,3),rep(country.time.match$yearD,3),rep(country.time.match$yearS,3),t(age.specific.mean[[i]]) + cent.val))#$year
  draws.temp                                                   <- get(draws.name)
  colnames(draws.temp)                                         <- c("country", "yearD", "yearS", paste("draw.", 1:length(burnt), sep=''))
  draws.temp$age                                               <- StandardPop$age_mean[i]
  assign(draws.name, draws.temp)
  ## mean draws
  mean.draws.name			                                         <- paste("mean.draws",i,sep='')
  assign(mean.draws.name,cbind(draws.temp[,1:3],draws.temp[,(length(burnt)+4)],apply(draws.temp[,4:(length(burnt)+3)],1,mean))) 
  mean.draws.temp                                              <- get(noquote(mean.draws.name))
  colnames(mean.draws.temp)                                    <- c("Country", "YearD", "YearS","Age","Mean")
  assign(mean.draws.name, mean.draws.temp)
  write.csv(mean.draws.temp,file=paste(filename,"_MeanDrawsAge",i,".csv",sep=''),row.names=FALSE)
}
########## Variance components for study-specific random effects ######

# I will find out the mean and the quantile for the 3 of them 
phi_natl1                                                     <- chains[[ch]]$phi_natl[burnt,]
phi_subn1                                                     <- chains[[ch]]$phi_subn[burnt,]
phi_comm1                                                     <- chains[[ch]]$phi_comm[burnt,]

m1                                                            <- apply(sqrt(exp(phi_natl1)),2,mean)
q1                                                            <- apply((sqrt(exp(phi_natl1))),2,function(x) quantile(x,c(.025, .975) ))
m2                                                            <- apply(sqrt(exp(phi_subn1)),2,mean)  
q2                                                            <- apply((sqrt(exp(phi_subn1))),2,function(x) quantile(x,c(.025, .975) ))
m3                                                            <- apply(sqrt(exp(phi_comm1)),2,mean)
q3                                                            <- apply((sqrt(exp(phi_comm1))),2,function(x) quantile(x,c(.025, .975) ))
xlim                                                          <- apply(rbind(q1, q2, q3),2 ,range)
ylim                                                          <- matrix(c(1,3),2,3)

pdf(paste(filename, "_sd_components",".pdf",sep=''))
par(mar=c(4.5,10,1.5,1.5))
for ( i in 1:3){
  plot(xlim[,i], ylim[,i], pch='', xlab='SD',ylab='', yaxt='n',main=title[i])
  axis(2, at=1:3, labels=c('National', 'Sub-national','Community'), las=1)
  lines(q1[,i], c(1,1));points(m1[i], 1); lines(q2[,i], c(2,2));points(m2[i], 2);lines(q3[,i], c(3,3));points(m3[i], 3)
}
dev.off()
foo                                           <- list(3,3,3)
names(foo)                                    <- title
for ( i in 1:3){
  a                                           <- data.frame(c("natl", "subn", "comm"),c(m1[i],m2[i],m3[i]),c(q1[1,i],q2[1,i],q3[1,i]),c(q1[2,i],q2[2,i],q3[2,i]))
  names(a)                                    <- c("coverage","sd","l","u")
  foo[[i]]                                    <- a
}
write.csv(foo,file=paste(filename, "_sd_components.csv", sep=''),row.names=FALSE)

################################### Covar effects ################
covar.names                                   <- c("D:sub-national","S:sub-national","INT:sub-national","D:sub-national x time","S:sub-national x time","INT:sub-national x time",
                                                   "D:community","S:community","INT:community","D:community x time","S:community x time","INT:community x time",
                                                   "D:perurb","S:perurb","INT:perurb","D:perurb x t","S:perurb x t","INT:perurb x t",
                                                   "D:perurb x rural","S:perurb x rural","INT:perurb x rural","D:perurb x rural x t ","S:perurb x rural x t ","INT:perurb x rural x t ",
                                                   "D:(1-perurb) x urban","S:(1-perurb) x urban","INT:(1-perurb) x urban",
                                                   "D:(1-perurb) x urban x t","S:(1-perurb) x urban x t","INT:(1-perurb) x urban x t",
                                                   "D:PC1","S:PC1","INT:PC1","D:PC2","S:PC2","INT:PC2","D:PC3","S:PC3","INT:PC3","D:PC4","S:PC4","INT:PC4") 

pdf(paste(filename, "_covariate_effects",".pdf",sep=''), width=10, height=7.5)
l.covar                                       <- apply(theta1[,(dim(F)[2]-p+1):dim(F)[2]], 2, quantile, .025 )
m.covar                                       <- apply(theta1[,(dim(F)[2]-p+1):dim(F)[2]], 2, mean)
u.covar                                       <- apply(theta1[,(dim(F)[2]-p+1):dim(F)[2]], 2, quantile, .975 )
par(mar=c(3.1,15, 4.1, 2.1))
plot(range(c(l.covar,u.covar)), c(1,p), pch='', yaxt='n', ylab='',main="Covariate effects", xlab='')
abline(v=0, col=grey(.7))
for(i in 1:p){
  lines(c(l.covar[i], u.covar[i]), c((p:1)[i],(p:1)[i]))
  points(m.covar[i], (p:1)[i], pch=19)
}
axis(2, at=p:1, labels=covar.names, las=2);box()
dev.off()
covar.estimates                               <- data.frame(covar.names,m.covar,l.covar,u.covar)
write.csv(covar.estimates, file=paste(filename, "_covariate_effects.csv", sep=''),row.names=FALSE)
###### Urbanisation effects #####
U_D 	                                        <- theta1[, dim(F)[2]-p+13] #D: urbanis
U_S 	                                        <- theta1[, dim(F)[2]-p+14] #S: urbanis
U_INT 	                                      <- theta1[, dim(F)[2]-p+15] #INT: urbanis
UT_D                                          <- theta1[, dim(F)[2]-p+16] #D: urbanis x time
UT_S                                          <- theta1[, dim(F)[2]-p+17] #S: urbanis x time
UT_INT                                        <- theta1[, dim(F)[2]-p+18] #INT: urbanis x time
RU_D 	                                        <- theta1[, dim(F)[2]-p+19] #D: rural x urbanis 
RU_S 	                                        <- theta1[, dim(F)[2]-p+20] #S: rural x urbanis 
RU_INT 	                                      <- theta1[, dim(F)[2]-p+21] #INT: rural x urbanis 
RUT_D 	                                      <- theta1[, dim(F)[2]-p+22] #D: rural x urbanis x time
RUT_S 	                                      <- theta1[, dim(F)[2]-p+23] #S: rural x urbanis x time
RUT_INT 	                                    <- theta1[, dim(F)[2]-p+24] #INT: rural x urbanis x time
NUU_D 	                                      <- theta1[, dim(F)[2]-p+25] #D: (1-urbanis) x urban
NUU_S 	                                      <- theta1[, dim(F)[2]-p+26] #S: (1-urbanis) x urban
NUU_INT 	                                    <- theta1[, dim(F)[2]-p+27] #INT: (1-urbanis) x urban
NUUT_D 	                                      <- theta1[, dim(F)[2]-p+28] #D: (1-urbanis) x urban x time
NUUT_S 	                                      <- theta1[, dim(F)[2]-p+29] #S: (1-urbanis) x urban x time
NUUT_INT 	                                    <- theta1[, dim(F)[2]-p+30] #INT: (1-urbanis) x urban x time
u.vals                                        <- seq(0,1,by=.1) 
###
# t.names                                     <- c(1990,1997,2005,2015) #c(1975,1990,2005,2016)
# t.vals                                      <- c(1990,1997,2005,2015)-centring.year #c(1975,1990,2005,2015)-1997 
#I ll use everywhere the coding years: -13 to 13 
t.names                                       <- c(start.year,start.year+3,end.year-3,end.year) #and I mean 2003.2003,2004.2004,2005.2005,2007.2007 
t.vals                                        <- c(start.year,start.year+3,end.year-3,end.year)-centring.year#2003:2007-centring.year

#####rural
t.val                                         <- t.vals[1]
r80 <- r90 <- r00 <- r08 <- m80 <- m90 <- m00 <- m08 <- u80 <- u90 <- u00 <- u08 <- list(matrix(NA, length(burnt), length(u.vals)),matrix(NA, length(burnt), length(u.vals)),matrix(NA, length(burnt), length(u.vals)))
names(r80)                                    <-  title
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  r80[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val + RU_D*u.val + RUT_D*u.val*t.val
  r80[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val + RU_S*u.val + RUT_S*u.val*t.val
  r80[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val + RU_INT*u.val + RUT_INT*u.val*t.val
}

r80.m<-r80.l<-r80.u<-r90.m <-r90.l<-r90.u<-r00.m<-r00.l<-r00.u<-r08.m<-r08.l<-r08.u <- list(rep(NA,length(u.vals)),rep(NA,length(u.vals)),rep(NA,length(u.vals)))
m80.m<-m80.l<-m80.u<-m90.m <-m90.l<-m90.u<-m00.m<-m00.l<-m00.u<-m08.m<-m08.l<-m08.u <- list(rep(NA,length(u.vals)),rep(NA,length(u.vals)),rep(NA,length(u.vals)))
u80.m<-u80.l<-u80.u<-u90.m <-u90.l<-u90.u<-u00.m<-u00.l<-u00.u<-u08.m<-u08.l<-u08.u <- list(rep(NA,length(u.vals)),rep(NA,length(u.vals)),rep(NA,length(u.vals)))
for ( i in 1:3)  {
  r80.m[[i]]                                  <- colMeans(r80[[i]]) 
  r80.l[[i]]                                  <- apply(r80[[i]], 2, quantile, .025)
  r80.u[[i]]                                  <- apply(r80[[i]], 2, quantile, .975)
}
t.val                                         <- t.vals[2]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  r90[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val + RU_D*u.val + RUT_D*u.val*t.val
  r90[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val + RU_S*u.val + RUT_S*u.val*t.val
  r90[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val + RU_INT*u.val + RUT_INT*u.val*t.val
}
for ( i in 1:3)  {
  r90.m[[i]]                                  <- colMeans(r90[[i]])
  r90.l[[i]]                                  <- apply(r90[[i]], 2, quantile, .025)
  r90.u[[i]]                                  <- apply(r90[[i]], 2, quantile, .975)
}
t.val                                         <- t.vals[3]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  r00[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val + RU_D*u.val + RUT_D*u.val*t.val
  r00[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val + RU_S*u.val + RUT_S*u.val*t.val
  r00[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val + RU_INT*u.val + RUT_INT*u.val*t.val
}
for ( i in 1:3) {
  r00.m[[i]]                                  <- colMeans(r00[[i]])
  r00.l[[i]]                                  <- apply(r00[[i]], 2, quantile, .025)
  r00.u[[i]]                                  <- apply(r00[[i]], 2, quantile, .975)
}
t.val                                         <- t.vals[4]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  r08[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val + RU_D*u.val + RUT_D*u.val*t.val
  r08[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val + RU_S*u.val + RUT_S*u.val*t.val
  r08[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val + RU_INT*u.val + RUT_INT*u.val*t.val
}
for ( i in 1:3) {
  r08.m[[i]]                                  <- colMeans(r08[[i]])
  r08.l[[i]]                                  <- apply(r08[[i]], 2, quantile, .025)
  r08.u[[i]]                                  <- apply(r08[[i]], 2, quantile, .975)
}
####### mixed
t.val                                         <- t.vals[1]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  m80[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val  
  m80[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val  
  m80[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val 
}
for ( i in 1:3) {
  m80.m[[i]]                                  <- colMeans(m80[[i]])
  m80.l[[i]]                                  <- apply(m80[[i]], 2, quantile, .025)
  m80.u[[i]]                                  <- apply(m80[[i]], 2, quantile, .975)
}
t.val                                         <- t.vals[2]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  m90[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val
  m90[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val  
  m90[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val 
}
for ( i in 1:3){
  m90.m[[i]]                                  <- colMeans(m90[[i]])
  m90.l[[i]]                                  <- apply(m90[[i]], 2, quantile, .025)
  m90.u[[i]]                                  <- apply(m90[[i]], 2, quantile, .975)
}
t.val                                         <- t.vals[3]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  m00[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val
  m00[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val
  m00[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val
}
for ( i in 1:3){
  m00.m[[i]]                                  <- colMeans(m00[[i]])
  m00.l[[i]]                                  <- apply(m00[[i]], 2, quantile, .025)
  m00.u[[i]]                                  <- apply(m00[[i]], 2, quantile, .975)
}
t.val                                         <- t.vals[4]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  m08[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val
  m08[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val
  m08[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val
}
for ( i in 1:3){
  m08.m[[i]]                                  <- colMeans(m08[[i]])
  m08.l[[i]]                                  <- apply(m08[[i]], 2, quantile, .025)
  m08.u[[i]]                                  <- apply(m08[[i]], 2, quantile, .975)
}
####### urban
t.val                                         <- t.vals[1]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  u80[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val + NUU_D*(1-u.val) + NUUT_D*(1-u.val)*t.val
  u80[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val + NUU_S*(1-u.val) + NUUT_S*(1-u.val)*t.val
  u80[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val + NUU_INT*(1-u.val) + NUUT_INT*(1-u.val)*t.val
}
for ( i in 1:3){
  u80.m[[i]]                                  <- colMeans(u80[[i]])
  u80.l[[i]]                                  <- apply(u80[[i]], 2, quantile, .025)
  u80.u[[i]]                                  <- apply(u80[[i]], 2, quantile, .975)
}
t.val                                         <- t.vals[2]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  u90[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val + NUU_D*(1-u.val) + NUUT_D*(1-u.val)*t.val
  u90[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val + NUU_S*(1-u.val) + NUUT_S*(1-u.val)*t.val
  u90[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val + NUU_INT*(1-u.val) + NUUT_INT*(1-u.val)*t.val
}
for ( i in 1:3){
  u90.m[[i]]                                  <- colMeans(u90[[i]])
  u90.l[[i]]                                  <- apply(u90[[i]], 2, quantile, .025)
  u90.u[[i]]                                  <- apply(u90[[i]], 2, quantile, .975)
}
t.val                                         <- t.vals[3]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  u00[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val + NUU_D*(1-u.val) + NUUT_D*(1-u.val)*t.val
  u00[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val + NUU_S*(1-u.val) + NUUT_S*(1-u.val)*t.val
  u00[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val + NUU_INT*(1-u.val) + NUUT_INT*(1-u.val)*t.val
}
for ( i in 1:3){
  u00.m[[i]]                                  <- colMeans(u00[[i]])
  u00.l[[i]]                                  <- apply(u00[[i]], 2, quantile, .025)
  u00.u[[i]]                                  <- apply(u00[[i]], 2, quantile, .975)
}
t.val                                         <- t.vals[4]
for(i in 1:length(u.vals)) {
  u.val                                       <- u.vals[i]
  u08[[1]][,i]                                <- U_D*u.val + UT_D*u.val*t.val + NUU_D*(1-u.val) + NUUT_D*(1-u.val)*t.val
  u08[[2]][,i]                                <- U_S*u.val + UT_S*u.val*t.val + NUU_S*(1-u.val) + NUUT_S*(1-u.val)*t.val
  u08[[3]][,i]                                <- U_INT*u.val + UT_INT*u.val*t.val + NUU_INT*(1-u.val) + NUUT_INT*(1-u.val)*t.val
}
for ( i in 1:3){
  u08.m[[i]]                                  <- colMeans(u08[[i]])
  u08.l[[i]]                                  <- apply(u08[[i]], 2, quantile, .025)
  u08.u[[i]]                                  <- apply(u08[[i]], 2, quantile, .975)
}
ylim                                          <- range(c(r80.l[[1]],r80.l[[2]],r80.l[[3]], r90.l[[1]],r90.l[[2]],r90.l[[3]],
                                                         r00.l[[1]],r00.l[[2]],r00.l[[3]], r08.l[[1]],r08.l[[2]],r08.l[[3]],
                                                         m80.l[[1]],m80.l[[2]],m80.l[[3]], m90.l[[1]],m90.l[[2]],m90.l[[3]],
                                                         m00.l[[1]],m08.l[[2]],u80.l[[3]], u90.l[[1]],u90.l[[2]],u90.l[[3]], 
                                                         u00.l[[1]],u00.l[[2]],u00.l[[3]], u08.l[[1]],u08.l[[2]],u08.l[[3]],
                                                         r80.u[[1]],r80.u[[2]],r80.u[[3]], r90.u[[1]],r90.u[[2]],r90.u[[3]],
                                                         r00.u[[1]],r00.u[[2]],r00.u[[3]], r08.u[[1]],r08.u[[2]],r08.u[[3]],
                                                         m80.u[[1]],m80.u[[2]],m80.u[[3]], m90.u[[1]],m90.u[[2]],m90.u[[3]],
                                                         m00.u[[1]],m00.u[[2]],m00.u[[3]], m08.u[[1]],u80.u[[2]],u90.u[[3]],
                                                         u00.u[[1]],u00.u[[2]],u00.u[[3]], u08.u[[1]],u08.u[[2]],u08.u[[3]]))
pdf(paste(filename, "_UrbanisationEffects.pdf",sep=''), width=10, height=7.5)
for (i in 1:3){
  par(mfrow=c(3,4), mar=c(4.6,4.4,2.7,1.3))
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]),c(r80.l[[i]][1], r80.u[[i]], r80.l[[i]][-order(u.vals)+length(u.vals)+1]),
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, r80.m[[i]])
  title(main=t.names[1], ylab="Rural Study", cex.main=1, font.main=1)
  mtext(title[i], side = 3, line = -1.2, outer = TRUE)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]), c(r90.l[[i]][1], r90.u[[i]], r90.l[[i]][-order(u.vals)+length(u.vals)+1]),
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, r90.m[[i]])
  title(main=t.names[2], cex.main=1, font.main=1)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]), c(r00.l[[i]][1], r00.u[[i]], r00.l[[i]][-order(u.vals)+length(u.vals)+1]),#,
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, r00.m[[i]])
  title(main=t.names[3], cex.main=1, font.main=1)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]), c(r08.l[[i]][1], r08.u[[i]], r08.l[[i]][-order(u.vals)+length(u.vals)+1]),#,
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, r08.m[[i]])
  title(main=t.names[4], cex.main=1, font.main=1)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]),c(m80.l[[i]][1], m80.u[[i]], m80.l[[i]][-order(u.vals)+length(u.vals)+1]),#,
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, m80.m[[i]])
  title(ylab="Mixed Study", cex.main=1, font.main=1)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]),c(m90.l[[i]][1], m90.u[[i]], m90.l[[i]][-order(u.vals)+length(u.vals)+1]),#,
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, m90.m[[i]])
  title(cex.main=1, font.main=1)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]),c(m00.l[[i]][1], m00.u[[i]], m00.l[[i]][-order(u.vals)+length(u.vals)+1]),#,
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, m00.m[[i]])
  title(cex.main=1, font.main=1)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]),c(m08.l[[i]][1], m08.u[[i]], m08.l[[i]][-order(u.vals)+length(u.vals)+1]),#,
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, m08.m[[i]])
  title(cex.main=1, font.main=1)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]),c(u80.l[[i]][1], u80.u[[i]], u80.l[[i]][-order(u.vals)+length(u.vals)+1]),#,
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, u80.m[[i]])
  title(ylab="Urban Study", xlab="Country % urban", cex.main=1, font.main=1)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]),c(u90.l[[i]][1], u90.u[[i]], u90.l[[i]][-order(u.vals)+length(u.vals)+1]),#,
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, u90.m[[i]])
  title(xlab="Country % urban", cex.main=1, font.main=1)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]),c(u00.l[[i]][1], u00.u[[i]], u00.l[[i]][-order(u.vals)+length(u.vals)+1]),#,
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, u00.m[[i]])
  title(xlab="Country % urban", cex.main=1, font.main=1)
  
  plot(range(u.vals), ylim, pch='', xaxt='n', xlab='', ylab='')
  axis(1,at=u.vals,seq(0,1,by=.1))
  polygon(c(u.vals[1], u.vals, u.vals[-order(u.vals)+length(u.vals)+1]),c(u08.l[[i]][1], u08.u[[i]], u08.l[[i]][-order(u.vals)+length(u.vals)+1]),#,
          col=grey(.85), border=NA)
  abline(h=0, col="grey")
  lines(u.vals, u08.m[[i]])
  title(xlab="Country % urban", cex.main=1, font.main=1)
}
dev.off()
############# Plot posterior means and intervals for u,v,w, theta1's #####
theta1_c <- chains[[ch]]$theta_c[burnt]
theta1_r <- chains[[ch]]$theta_r[burnt]
theta1_s <- chains[[ch]]$theta_s[burnt]
theta1_g <- chains[[ch]]$theta_g[burnt]
####
u.hat		      <- colMeans(u1)
v.hat  		    <- colMeans(v1)
sv.hat  	    <- colMeans(sv1)
w.hat   	    <- colMeans(w1)
u.l   		    <- apply(u1,  2, quantile, .025)
u.u   		    <- apply(u1,  2, quantile, .975)
v.l   		    <- apply(v1,  2, quantile, .025)
v.u   		    <- apply(v1,  2, quantile, .975)
sv.l   		    <- apply(sv1, 2, quantile, .025)
sv.u   		    <- apply(sv1, 2, quantile, .975)
w.l   		    <- apply(w1,  2, quantile, .025)
w.u   	     	<- apply(w1,  2, quantile, .975)
theta1_c.hat 	<- mean(theta1_c) ; theta1_r.hat <- mean(theta1_r)
theta1_s.hat 	<- mean(theta1_s) ; theta1_g.hat <- mean(theta1_g)
theta1_c.l 	  <- quantile(theta1_c,c(.025,.975))[[1]]; theta1_c.u 	<- quantile(theta1_c,c(.025,.975))[[2]]
theta1_r.l 	  <- quantile(theta1_r,c(.025,.975))[[1]]; theta1_r.u 	<- quantile(theta1_r,c(.025,.975))[[2]]
theta1_s.l  	<- quantile(theta1_s,c(.025,.975))[[1]]; theta1_s.u 	<- quantile(theta1_s,c(.025,.975))[[2]]
theta1_g.l 	  <- quantile(theta1_g,c(.025,.975))[[1]]; theta1_g.u 	<- quantile(theta1_g,c(.025,.975))[[2]]
ylim 	      	<- range(c(u.hat, v.hat,sv.hat, w.hat))
ylim[1]       <- ylim[1]-.1;ylim[2] <- ylim[2]+.1

pdf(paste(filename, "_uvw.hat",".pdf",sep=''), width=10, height=7.5)
par(mar=c(.75,.1,.1,.75), mfrow=c(7,6))
plot(0, 0, pch='', xlim=c(0,T^2), ylim=ylim,xaxt='n', xlab='', yaxt='n', ylab='', main='')
abline(h=seq(-1, 1, by=1) , col='grey')
legend("topright", pch='', legend='w', bty='n')
axis(4)
for(i in 1:(T^2) ) {
  lines(c(i,i), c(w.l[i], w.u[i]))
  points(i, w.hat[i], pch=19)
}
for(l in 1:L) {
  plot(0, 0, pch='', xlim=c((l - 1) * (T^2) + 1,(l*(T^2))), ylim=ylim, yaxt='n', xaxt='n', xlab='', ylab='', main='')
  abline(h = seq(-1, 1, by=1) , col='grey')
  legend("topright", pch='', bty='n',legend=paste("sv:",tolower(region.match$name[l]))) #l
  for(i in ((l - 1) * (T^2) + 1):(l*(T^2))) {
    lines(c(i,i), c(sv.l[i], sv.u[i]))
    points(i, sv.hat[i], pch=19)
  }
}
axis(4)
for(k in 1:K) {
  plot(0, 0, pch='', xlim=c((k - 1) * (T^2) + 1,(k*(T^2))), ylim=ylim, yaxt='n', xaxt='n', xlab='', ylab='', main='')
  abline(h=seq(-1, 1, by=1) , col='grey')
  legend("topright", pch='', bty='n',legend=paste("v:", tolower(region.match$name[k])))
  for(i in ((k - 1) * (T^2) + 1):(k*(T^2))) {
    lines(c(i,i), c(v.l[i], v.u[i]))
    points(i, v.hat[i], pch=19)
  }
}
axis(4)
for(k in 1:J) {
  plot(0, 0, pch='', xlim=c((k - 1) * (T^2) + 1,(k*(T^2))), ylim=ylim, yaxt='n', xaxt='n', xlab='', ylab='', main='')
  abline(h=seq(-1, 1, by=1) , col='grey')
  legend("topright", pch='', bty='n', legend=paste("u:", tolower(country.match$name[k])))
  for(i in ((k - 1) * (T^2) + 1):(k*(T^2))) {
    lines(c(i,i), c(u.l[i], u.u[i]))
    points(i, u.hat[i], pch=19)
  }
}
axis(4)
ylim <- range(theta1_c.l, theta1_r.l,theta1_s.l, theta1_g.l, theta1_c.u, theta1_r.u,theta1_s.u, theta1_g.u)
plot(0, 0, pch='', xlim=c(1,4), ylim=ylim, yaxt='n', xaxt='n', xlab='', ylab='', main='')
legend("topleft", pch='', bty='n', legend="log lambda's")
for ( i in 1:3) {
  lines(c(1,1), c(theta1_c.l[i], theta1_c.u[i]))
  points(1, theta1_c.hat[i], pch=19)
  text(1, theta1_c.hat[i], "c", pos=4)
  lines(c(2,2), c(theta1_r.l[i], theta1_r.u[i]))
  points(2, theta1_r.hat[i], pch=19)
  text(2, theta1_r.hat[i], "r", pos=4)
  lines(c(3,3), c(theta1_s.l[i], theta1_s.u[i]))
  points(3, theta1_s.hat[i], pch=19)
  text(3, theta1_s.hat[i], "s", pos=4)
  lines(c(4,4), c(theta1_g.l[i], theta1_g.u[i]))
  points(4, theta1_g.hat[i], pch=19)
  text(4, theta1_g.hat[i], "g", pos=2)
  axis(4)
}
dev.off()
##### Age standardized country mean 20+ ###############################
cent.val                                             <- matrix(rep(c(dbp.cent.val,sbp.cent.val,int.cent.val),each=J*(T^2)),1,3*J*(T^2) )
age.standardized.mean                                <- age.specific.mean[[1]] * StandardPop$weights[1]
for(j in 2:nrow(StandardPop)) age.standardized.mean  <- age.standardized.mean + age.specific.mean[[j]] * StandardPop$weights[1]
age.standardized.mean 	                             <- age.standardized.mean + cent.val[rep(seq_len(nrow(cent.val)),length(burnt)), ]

age.standard.country.mean                            <- colMeans(age.standardized.mean)
age.standard.country.l                               <- apply(age.standardized.mean, 2, quantile, .025)
age.standard.country.u                               <- apply(age.standardized.mean, 2, quantile, .975)
age.standard.country.se                              <- apply(age.standardized.mean, 2, sd )
country.means                                        <- data.frame(rep(country.time.match$country,3),rep(country.time.match$yearD,3),rep(t[,3],3*J), age.standard.country.mean, age.standard.country.l,age.standard.country.u, age.standard.country.se)#$year
names(country.means)                                 <- c("country", "yearD","cod.year", "mean", "l", "u", "se")
##### Age standardized country mean 20-64 #############################
max.young.age 		                                 	 <- grep("64",StandardPop$age_group)
age.standardized.mean.young                       	 <- age.specific.mean[[1]]*(StandardPop$StandardPop[1]/sum(StandardPop$StandardPop[1:max.young.age]))
for(j in 2:max.young.age) age.standardized.mean.young<- age.standardized.mean.young + age.specific.mean[[j]]*(StandardPop$StandardPop[j]/sum(StandardPop$StandardPop[1:max.young.age]))
age.standardized.mean.young 	                       <- age.standardized.mean.young + cent.val[rep(seq_len(nrow(cent.val)), length(burnt)), ]

age.standard.country.young.mean                      <- colMeans(age.standardized.mean.young )
age.standard.country.young.l                         <- apply(age.standardized.mean.young, 2, quantile,.025 )
age.standard.country.young.u                         <- apply(age.standardized.mean.young, 2, quantile,.975 )
age.standard.country.young.se                        <- apply(age.standardized.mean.young, 2, sd)
country.means.young   	                             <- data.frame(rep(country.time.match$country,3),rep(country.time.match$yearD,3),rep(t[,3],3*J), age.standard.country.young.mean, 
                                                                  age.standard.country.young.l, age.standard.country.young.u, age.standard.country.young.se)
names(country.means.young)                           <- c("country", "yearD", "cod.year", "mean", "l", "u", "se")
##### Age standardized country mean 65+ ###############################
min.old.age 			                                   <- grep("65",StandardPop$age_group)
age.standardized.mean.old                            <- age.specific.mean[[min.old.age]]*(StandardPop$StandardPop[min.old.age]/sum(StandardPop$StandardPop[min.old.age:nrow(StandardPop)]))
for(j in (min.old.age+1):nrow(StandardPop))  age.standardized.mean.old <- age.standardized.mean.old + age.specific.mean[[j]]*(StandardPop$StandardPop[j]/sum(StandardPop$StandardPop[min.old.age:nrow(StandardPop)]))
age.standardized.mean.old 	                         <- age.standardized.mean.old + cent.val[rep(seq_len(nrow(cent.val)), length(burnt)), ]

age.standard.country.old.mean 	                     <- colMeans(age.standardized.mean.old)
age.standard.country.old.l                        	 <- apply(age.standardized.mean.old, 2, quantile,.025 )
age.standard.country.old.u                        	 <- apply(age.standardized.mean.old, 2, quantile,.975 )
age.standard.country.old.se                       	 <- apply(age.standardized.mean.old, 2, sd)
country.means.old 	                      	         <- data.frame(rep(country.time.match$country,3),rep(country.time.match$yearD,3), rep(t[,3],3*J),
                                                                 age.standard.country.old.mean, age.standard.country.old.l, age.standard.country.old.u,age.standard.country.old.se)
names(country.means.old)                             <- c("country", "yearD", "cod.year", "mean", "l", "u", "se")
##### Combine and save country means ##################################
all.means                                            <- cbind(country.means, country.means.young[,4:7], country.means.old[,4:7])
names(all.means)                                     <- c("country", "yearD", "cod.year","mean.all","l.all","u.all","se.all","mean.young","l.young","u.young","se.young","mean.old","l.old","u.old", "se.old")
write.csv(c(all.means)  , paste(filename, "_CountryMeans.csv", sep=''),row.names=FALSE)

##### Age-standardized country-specific linear time slopes ############
slope.vals                                          <- matrix(NA,length(burnt),3)
country.slop                                        <- data.frame(matrix(NA, J, 6))
names(country.slop)                                 <- c("country", "slope.all",  "l.all", "u.all",  "se.all", "p.all")
country.slopes                                      <- list(country.slop,country.slop,country.slop)
names(country.slopes)                               <- title  
for (i in 1:3)      country.slopes[[i]][,1]         <- as.character(country.match$name)
for (i in 1:J){
  country.val                                       <- as.character(country.match$name[i])
  country.specific.age.standardized.meansD          <- age.standardized.mean[,which(country.time.match$country==country.val)]
  country.specific.age.standardized.meansS          <- age.standardized.mean[,(which(country.time.match$country==country.val)+J*(T^2)) ]
  country.specific.age.standardized.meansINT        <- age.standardized.mean[,(which(country.time.match$country==country.val)+2*J*(T^2))]
  
  # slope.vals[,1]                                  <- lm(t(country.specific.age.standardized.meansD)  ~ c(start.year:end.year))$coef[2,] # it is 25 as T^2 repetition for each year*5
  # slope.vals[,2]                                  <- lm(t(country.specific.age.standardized.meansS)  ~ c(start.year:end.year))$coef[2,]
  # slope.vals[,3]                                  <- lm(t(country.specific.age.standardized.meansINT)~ c(start.year:end.year))$coef[2,]
  ####
  # slope.vals[,1]                                  <- lm(t(country.specific.age.standardized.meansD)  ~  c(rep(start.year:end.year,each=T)))$coef[2,]
  # slope.vals[,2]                                  <- lm(t(country.specific.age.standardized.meansS)  ~  c(rep(start.year:end.year,each=T)))$coef[2,]
  # slope.vals[,3]                                  <- lm(t(country.specific.age.standardized.meansINT)~  c(rep(start.year:end.year,each=T)))$coef[2,]
  slope.vals[,1]                                    <- lm(t(country.specific.age.standardized.meansD)  ~  c(1:(T^2)))$coef[2,]# non-central expression of the years
  slope.vals[,2]                                    <- lm(t(country.specific.age.standardized.meansS)  ~  c(1:(T^2)))$coef[2,]
  slope.vals[,3]                                    <- lm(t(country.specific.age.standardized.meansINT) ~  c(1:(T^2)))$coef[2,]
  # slope.vals[,1]                                    <- lm(t(country.specific.age.standardized.meansD)  ~ (cbind(rep(start.year:end.year,each=T),rep(start.year:end.year,T)) ))$coef[2,]
  # slope.vals[,2]                                    <- lm(t(country.specific.age.standardized.meansS)  ~ (cbind(rep(start.year:end.year,each=T),rep(start.year:end.year,T)) ))$coef[2,]
  # slope.vals[,3]                                    <- lm(t(country.specific.age.standardized.meansINT) ~ (cbind(rep(start.year:end.year,each=T),rep(start.year:end.year,T)) ))$coef[2,]
  for (jj in 1:3){  
    country.slopes[[jj]][i,2]                       <- mean(slope.vals[,jj])
    country.slopes[[jj]][i,3]                       <- quantile(slope.vals[,jj], .025)
    country.slopes[[jj]][i,4]                       <- quantile(slope.vals[,jj], .975)
    country.slopes[[jj]][i,5]                       <- sd(slope.vals[,jj])
    country.slopes[[jj]][i,6]                       <- mean(sign(slope.vals[,jj])!=sign(mean(slope.vals[,jj])))
  }}
write.csv(country.slopes, file=paste(filename,"_CountrySlopes.csv", sep=''),row.names=FALSE)
##### Fit plots including age-standardised ###################
#CountryNames                                       <- read.csv("C:\\Users\\mariz\\Dropbox\\PhD_Project\\papers_code_data\\Adult_BMI_model54_test\\country_list2015-09-15.csv")
CountryNames                                        <- subset[,11:13]
CountryNames                                        <- unique(covar$country)[order(unique(covar$country))]#levels(subset[,11])#as.character(CountryNames[order(as.character(CountryNames[,1])),1])
#attach(subset)
coverage.col                                        <- rep(0, length(subset$mean_dbp)) #mean_sbp
coverage.col[subset$coverage=="National"]           <- rainbow(50)[30]
coverage.col[subset$coverage=="Subnational"]        <- rainbow(50)[8]
coverage.col[subset$coverage=="Community"]          <- rainbow(50)[46]
scope.pch                                           <- rep(16, length(subset$mean_dbp)) #
scope.pch[which(subset$scope=="rural")]             <- 15
scope.pch[which(subset$scope=="urban")]             <- 17 
R.age1                                              <- as.matrix(chains[[ch]]$R.age)#,dim(ageMat)[1]*3,30+3*J*5)
countr_means                                        <- list(country.means[1:(J*(T^2)),],country.means[(1+J*(T^2)):(2*J*(T^2)),],country.means[(1+2*J*(T^2)):(3*J*(T^2)),])
se_int                                              <- rho*sqrt(subset$se_dbp*subset$se_sbp)#rep(rho<-cor(y1,y2) ,length(se_dbp))
se                                                  <- list(subset$se_dbp,subset$se_sbp,se_int)
cent                                                <- c(dbp.cent.val,sbp.cent.val,int.cent.val)
minBP                                               <- c(0,0,8000)
maxBP                                               <- c(200,250,12000)
pdf(paste(filename, "_CountryFits.pdf", sep=''), width=7.5, height=10)
for ( d in 1:3){
  for(k in 1:length(unique(covar$country))){
    par(mfrow=c(4,3), mar=c(2.2,2,4,1))
    plot(0,0,pch='',xaxt='n',yaxt='n',bty='n',xlab='',ylab='')
    legend("topleft", pch=c(16,16,16,16,15,17,16), bty='n',c("national", "sub-national", "community", "", "rural", "urban", "mixed"),
           col=c(rainbow(50)[30],rainbow(50)[8],rainbow(50)[46], 'white', 'black', 'black', 'black'), title="Legend")
    country.val                                       <- CountryNames[k];print(country.val)
    sregion.val 	                                    <- unique(covar$sregion[covar$country==country.val])
    region.val                                        <- unique(covar$region[covar$country==country.val])
    l                                                 <- unique(covar$sregion[covar$country==country.val])
    which                                             <- which(country.time.match$country == as.character(country.val))
    plot(c(min(t[,3]),max(t[,3])), c(minBP[d],maxBP[d]), pch='', xlab='', ylab='', xaxt='n', yaxt='n')
    axis(1);axis(2)
    if(sex=="female") foo.sexVal           <- "Female" else if (sex=="male") foo.sexVal <- "Male" 
    mtext(side=3, at= max(countr_means[[d]]$mean) + cent[d], line=0, text=foo.sexVal, cex=.75, adj=1)      
    polygon(c(t[1,3],t[,3], t[(T^2):1,3] ),c(countr_means[[d]]$l[which][1],countr_means[[d]]$u[which][1:(T^2)],countr_means[[d]]$l[which][(T^2):1]), col=grey(.85), border=NA)   
    lines(t[,3], countr_means[[d]]$mean[which][1:(T^2)], lwd=3)  
    uids                                              <- length(unique(subset$uid[as.character(subset$country)==country.val]))
    uidNatl                                           <- length(unique(subset$uid[as.character(subset$country)==country.val & subset$coverage=="National"]))
    title(main=paste(CountryNames[k], "(", tolower(sex), ")", "," ,title[d],sep=''),line=2.6, cex.main=1)
    title(main=sregion.val,font.main=1, cex.main=.9, line=0.6); title(main=region.val, font.main=1, cex.main=.9, line=1.6)
    if(uids == 0) text(max(t[,3]), min(countr_means[[d]]$mean) + cent[d], adj=c(1,0), cex=.8, "No data sources")  
    if(uids == 1) text(max(t[,3]), min(countr_means[[d]]$mean) + cent[d], adj=c(1,0), cex=.8, paste("1 data source, ",uidNatl," national", sep=''))
    if(uids > 1)  text(max(t[,3]), min(countr_means[[d]]$mean) + cent[d], adj=c(1,0), cex=.8, paste(uids," data sources, ",uidNatl," national",sep=''))  
    box()
    plot(0,0,pch='',xaxt='n',yaxt='n',bty='n',xlab='',ylab='')
    counter                                           <- 1
    if(sum(subset$country==as.character(country.val))>0){
      country.val.num                                 <- country.match$number[country.match$name==country.val]
      region.val 	                                    <- as.character(unique(covar$region[covar$country==country.val]))
      region.val.num                                  <- region.match$number[region.match$name==region.val]
      sregion.val       	                            <- unique(covar$sregion[covar$country==country.val])
      sregion.val.num                                 <- sregion.match$number[sregion.match$name==sregion.val]
      #sregion 
      sregion_pred 	<- sregion.time_pred              <- matrix(0, L, T^2)
      sregion_pred[sregion.val.num,]                  <- 1
      sregion.time_pred[sregion.val.num,]             <- t[,3]
      #region 
      region_pred 	<- region.time_pred               <- matrix(0, K, T^2)
      region_pred[region.val.num,]                    <- 1
      region.time_pred[region.val.num,]               <- t[,3]
      region_pred                                     <- region_pred[multipleRegionsInSregion,]			
      region.time_pred                                <- region.time_pred[multipleRegionsInSregion,]	
      #country
      country_pred 	<- country.time_pred              <- matrix(0, J, T^2)
      country_pred[country.val.num,]                  <- 1
      country.time_pred[country.val.num,]             <- t[,3]
      F_pr                                            <- rbind(matrix(1,3,T^2), matrix(t[,3],3,T^2,byrow = T), 
                                                               matrix(sregion_pred,3*L,T^2),rbind(sregion.time_pred,sregion.time_pred,sregion.time_pred), 
                                                               matrix(region_pred,3*(sum(multipleRegionsInSregion)),T^2),rbind(region.time_pred,region.time_pred,region.time_pred),
                                                               matrix(country_pred,3*J,T^2),rbind(country.time_pred,country.time_pred,country.time_pred),
                                                               matrix(0, 3*N,T^2), matrix(0, 3*4,T^2),#
                                                               t(matrix(covar[,covariate.list[1]][covar$country==country.val],T^2,3)),   	
                                                               t(matrix(covar[,covariate.list[1]][covar$country==country.val]*t,T^2,3)),   	
                                                               matrix(0, 3*4,T^2),matrix(0, 3*4,T^2)) #only the first one
      F_pred                                          <- cbind(F_pr,F_pr,F_pr)
      all                                             <- 6*(1+L+(sum(multipleRegionsInSregion))+J)
      #DBP
      F_pred[c(1,4,7:(6+L),(7+3*L):(6+4*L),(7+6*L):(6+6*L+(sum(multipleRegionsInSregion))),(7+6*L+3*(sum(multipleRegionsInSregion))):(6+6*L+4*(sum(multipleRegionsInSregion))),
               (7+6*L+6*(sum(multipleRegionsInSregion))):(6+6*L+6*(sum(multipleRegionsInSregion))+J),(7+6*L+6*(sum(multipleRegionsInSregion))+3*J):(6+6*L+6*(sum(multipleRegionsInSregion))+4*J),  
               (all+1):(all+dim(StudyToUid)[2]),all+1+3*dim(StudyToUid)[2]+seq(0,3*13,3)),(dim(F_pr)[2]+1):(3*dim(F_pr)[2])]  <- 0
      #all+1+seq(0,3*13,3)),(dim(F_pr)[2]+1):(3*dim(F_pr)[2])]  <- 0
      #SBP
      F_pred[c(2,5,(7+L):(6+2*L),(7+4*L):(6+5*L),(7+6*L+(sum(multipleRegionsInSregion))):(6+6*L+2*(sum(multipleRegionsInSregion))),
               (7+6*L+4*(sum(multipleRegionsInSregion))):(6+6*L+5*(sum(multipleRegionsInSregion))),(7+6*L+6*(sum(multipleRegionsInSregion))+J):(6+6*L+6*(sum(multipleRegionsInSregion))+2*J),(7+6*L+6*(sum(multipleRegionsInSregion))+4*J):(6+6*L+6*(sum(multipleRegionsInSregion))+5*J),
               (all+1+dim(StudyToUid)[2]):(all+2*dim(StudyToUid)[2]),all+2+3*dim(StudyToUid)[2]+seq(0,3*13,3)),c(1:dim(F_pr)[2],(2*dim(F_pr)[2]+1):(3*dim(F_pr)[2]))] <- 0
      #all+2+seq(0,3*13,3)),c(1:(dim(F_pr)[2]),(2*dim(F_pr)[2]+1):(3*dim(F_pr)[2]))] <- 0
      #INTER
      F_pred[c(3,6,(7+2*L):(6+3*L),(7+5*L):(6+6*L),(7+6*L+2*(sum(multipleRegionsInSregion))):(6+6*L+3*(sum(multipleRegionsInSregion))),(7+6*L+5*(sum(multipleRegionsInSregion))):(6+6*L+6*(sum(multipleRegionsInSregion))),
               (7+6*L+6*(sum(multipleRegionsInSregion))+2*J):(6+6*L+6*(sum(multipleRegionsInSregion))+3*J),(7+6*L+6*(sum(multipleRegionsInSregion))+5*J):(6+6*L+6*(sum(multipleRegionsInSregion))+6*J),
               (all+1+2*dim(StudyToUid)[2]):(all+3*dim(StudyToUid)[2]),all+3+3*dim(StudyToUid)[2]+seq(0,3*13,3)),1:(2*dim(F_pr)[2])] <- 0
      #all+3+seq(0,3*13,3)),1:(2*dim(F_pr)[2])] <- 0
      uu                                              <- rep(which(country.time.match$country == as.character(country.val)),3)
      vv                                              <- rep(which(region.time.match$region == region.val),3)
      ssv                                             <- rep(which(sregion.time.match$sregion == sregion.val),3)
      ww                                              <- rep(1:T^2,3)
      #The baseline for each country
      baseline.c                                      <- theta1 %*% F_pred + u1[,uu] +v1[,vv] + sv1[,ssv] + w1[,ww]
      baseline.m                                      <- colMeans(baseline.c) # for D,S,INT 
      baselineD.l                                     <- apply(baseline.c[,1:(T^2)], 2, quantile, .025)
      baselineS.l                                     <- apply(baseline.c[,(T^2+1):(2*(T^2))],  2, quantile, .025)
      baselineINT.l                                   <- apply(baseline.c[,(2*(T^2)+1):(3*(T^2))],2, quantile, .025)
      baselineD.u                                     <- apply(baseline.c[,1:(T^2)], 2, quantile, .975)
      baselineS.u                                     <- apply(baseline.c[,(T^2+1):(2*(T^2))],  2, quantile, .975)
      baselineINT.u                                   <- apply(baseline.c[,(2*(T^2)+1):(3*(T^2))],2, quantile, .975)
      age.vals                                        <- unique(subset$age[subset$country==as.character(country.val)])
      for(j in 1:length(age.vals)){
        age.val                                       <- age.vals[order(age.vals)][j] 
        where                                         <- which(subset$country==as.character(country.val) & subset$age==age.val)
        which1 <- list();                      which1 <- list(where,where+I,where+2*I)
        R.age.val                                     <- R.age1[which1[[d]],]
        if(is.matrix(R.age.val)) R.age.val            <- R.age.val[1,] 
        R.age.val                                     <- matrix(R.age.val, length(burnt), length(R.age.val), byrow=T)
        age.specific.m <- age.specific.l              <- age.specific.u <- list(rep(NA, T^2),rep(NA, T^2),rep(NA,T^2)) 
        for(m in 1:(T^2)){
          R.age.val[,16:30]                           <- R.age.val[,1:15] 
          R.age.val[,16:30]                           <- R.age.val[,1:15]*baseline.c[,m]
          Dage.specific                               <- baseline.c[,m] + rowSums(R.age.val[,c(1:5,16:20,31:(30+J*5))]*gamma1[,c(1:5,16:20,31:(30+J*5))])
          Sage.specific                               <- baseline.c[,(T^2)+m]  + rowSums(R.age.val[,c(6:10,21:25,(31+J*5):(30+2*J*5))]*gamma1[,c(6:10,21:25,(31+J*5):(30+2*J*5))])
          INTage.specific                             <- baseline.c[,2*(T^2)+m]+ rowSums(R.age.val[,c(11:15,26:30,(31+2*J*5):(30+3*J*5))]*gamma1[,c(11:15,26:30,(31+2*J*5):(30+3*J*5))])
          age.specific                                <- list(Dage.specific,Sage.specific,INTage.specific)
          for (dd in 1:3){
            age.specific.m[[dd]][m]                   <- mean(age.specific[[dd]])
            age.specific.l[[dd]][m]                   <- quantile(age.specific[[dd]],.025)
            age.specific.u[[dd]][m]                   <- quantile(age.specific[[dd]],.975)
          }
        }  
        plot(c(min(t[,3]),max(t[,3])), c(minBP[d],maxBP[d]), pch='', xlab='', ylab='', xaxt='n', yaxt='n')
        if (counter %% 3 == 1) axis(2)  
        if (counter %% 12 %in% c(0,10,11) | j==length(age.vals) | j==length(age.vals)-1 | j==length(age.vals)-2) axis(1)
        polygon(c(t[1,3],t[,3], t[(T^2):1,3]), c(age.specific.l[[d]][1]+cent[d], age.specific.u[[d]] + cent[d],age.specific.l[[d]][(T^2):1]+cent[d]), col=grey(.85), border=NA)
        lines(t[,3], age.specific.m[[d]] + cent[d], lwd=3)
        which                                         <- which(subset$country==as.character(country.val)&subset$age==age.val) 
        points(time[which], y[((d-1)*I+1):(d*I)][which] + cent[d], col=coverage.col[which], pch=scope.pch[which]) # 
        sem.vals                                      <- se[[d]][which] ## 
        for(i in 1:length(which))
          lines(c(time[which][i], time[which][i]),c(y[((d-1)*I+1):(d*I)][which][i] + cent[d] - 2* sem.vals[i],
                                                    y[((d-1)*I+1):(d*I)][which][i] + cent[d] + 2* sem.vals[i]), col=coverage.col[which][i])
        
        legend("topright", pch='', bty='n', legend=paste("mid-age =", age.val+50))
        if(sex=="female")  foo.sex.val   <- "Female" else if (sex=="male") foo.sex.val <- "Male"  
        if (j > 32 & ((counter+3) %% 12 == 2 | ((counter+3) %% 12 == 1 & j == length(age.vals)))) 
          title(main=paste(CountryNames[k], " (", tolower(sex), "), 4", sep=''),line=2.6, cex.main=1)
        if (j > 21 & j <= 32 & ((counter+3) %% 12 == 2 | ((counter+3) %% 12 == 1 & j == length(age.vals)))) 
          title(main=paste(CountryNames[k], " (", tolower(sex), "), 3", sep=''),line=2.6, cex.main=1)
        if (j<=21 & j > 9 & ((counter+3) %% 12 == 2 | ((counter+3) %% 12 == 1 & j == length(age.vals))))
          title(main=paste(CountryNames[k], " (", tolower(sex), "), 2", sep=''),line=2.6, cex.main=1)   		
        mtext(side=3, at=max(y[((d-1)*I+1):(d*I)]+cent[d]), line=0, text=foo.sex.val, cex=.75, adj=1)
        counter                                       <- counter + 1
        box()
      }}}
}
dev.off()
################# Age trends ##############
age.vals                                              <- seq(18, 80, by = 2) - 50
# ageValMinus45Plus                                   <- age.vals + 50 - 45
# ageValMinus45Plus[age.vals < -50 + 45]              <- 0
# ageValMinus60Plus                                   <- age.vals + 50 - 60
# ageValMinus60Plus[age.vals < -50 + 60]              <- 0
#ageMat.val                                           <- cbind(age.vals, age.vals^2, age.vals^3,ageValMinus45Plus^3,ageValMinus60Plus^3)
ageMat.val                                            <- bs(age.vals,df = 5)[,1:5]
gamma1.hat                                            <- colMeans(gamma1)
y.tilde.vals                                          <- seq(-5,11,by=4)
ylim.valmin                                           <- rep(0,3)
ylim.valmax                                           <- c(160,180,200)
title.val                                             <- "mean BP of 50 y.o.'s"
pred                                                  <- list(matrix(NA, length(y.tilde.vals), length(age.vals)),
                                                              matrix(NA, length(y.tilde.vals), length(age.vals)),
                                                              matrix(NA, length(y.tilde.vals), length(age.vals)) )
dimns                                                 <- list(c(1:5,16:20,31:(30+J*5)),c(6:10,21:25,(31+J*5):(30+2*J*5)),
                                                              c(11:15,26:30,(31+2*J*5):(30+3*J*5)))
for(i in 1:length(y.tilde.vals)) {
  for(j in 1:length(age.vals)) {
    R.age1                                            <- c(ageMat.val[j,], ageMat.val[j,],ageMat.val[j,], ageMat.val[j,],ageMat.val[j,], ageMat.val[j,])
    R.age1[16:20]                                     <- R.age1[16:20]*y.tilde.vals[i]
    R.age1[21:25]                                     <- R.age1[21:25]*y.tilde.vals[i]
    R.age1[26:30]                                     <- R.age1[26:30]*y.tilde.vals[i]
    pred[[1]][i,j]                                    <- y.tilde.vals[i] + R.age1[c(1:5,16:20)] %*% gamma1.hat[c(1:5,16:20)]  + dbp.cent.val 
    pred[[2]][i,j]                                    <- y.tilde.vals[i] + R.age1[c(6:10,21:25)]%*% gamma1.hat[c(6:10,21:25)] + sbp.cent.val 
    pred[[3]][i,j]                                    <- y.tilde.vals[i] + R.age1[c(11:15,26:30)]%*%gamma1.hat[c(11:15,26:30)]+ int.cent.val 
  }}
if(sex=="female")   	main.val                        <- 'Age association with BP in females' else
  if(sex=="male")   	main.val                        <- 'Age association with BP in males'
pdf(paste(filename, "_AgeTrends.pdf", sep=''))
for ( d in 1:3){
  plot(range(age.vals) + 50, c(ylim.valmin[d], ylim.valmax[d]), pch='', xlab='age', ylab= title[d], main=paste(main.val,',',title[d]))
  for(i in 1:length(y.tilde.vals))  lines(age.vals+50, pred[[d]][i,], col=rainbow(length(y.tilde.vals))[i])
  legend("topleft", col=rainbow(length(y.tilde.vals)), legend=y.tilde.vals + cent[d], bty='n', lty=1, title=title.val)
}
dev.off()
pdf(paste(filename, "_CountryAgeTrends.pdf", sep=''))
for (d in 1:3){
  if(sex=="female") main.val                           <- 'Country-specific BP vs. age trends in females' else
    if(sex=="male")   main.val                         <- 'Country-specific BP vs. age trends in males'
    
    plot(range(age.vals) + 50, c(ylim.valmin[d]-50,ylim.valmax[d]+50), pch='',xlab='Age', ylab=title[d], main=main.val)
    for(k in 1:length(country.match$name)) { 
      country.val                                      <- as.character(country.match$name[k])
      country.val.num                                  <- country.match$number[country.match$name==country.val]
      y.hat                                            <- mean(baseline[,country.time.match$country==country.val])
      R.age.v                                          <- cbind(ageMat.val,ageMat.val,ageMat.val, ageMat.val,ageMat.val,ageMat.val,
                                                                matrix(0, length(age.vals),5*J),matrix(0, length(age.vals),5*J),
                                                                matrix(0, length(age.vals),5*J))
      R.age.val                                        <- rbind(R.age.v,R.age.v,R.age.v)
      R.age.val[1:dim(ageMat.val)[1],c(6:15,21:30,(J*5+31):(3*J*5+30))] <- 0
      R.age.val[(1+dim(ageMat.val)[1]):(2*dim(ageMat.val)[1]),c(1:5,11:20,26:(J*5+30),(2*J*5+31):(3*J*5+30))]  <- 0
      R.age.val[(1+2*dim(ageMat.val)[1]):(3*dim(ageMat.val)[1]),c(1:10,16:25,31:(2*J*5+30))]         <- 0
      
      R.age.val[1:dim(ageMat.val)[1],30+country.val.num+J*(0:4)]        <- ageMat.val
      R.age.val[(1+dim(ageMat.val)[1]):(2*dim(ageMat.val)[1]),(30+5*J)+country.val.num+J*(0:4)]      <- ageMat.val
      R.age.val[(1+2*dim(ageMat.val)[1]):(3*dim(ageMat.val)[1]),(30+2*5*J)+ country.val.num+J*(0:4)] <- ageMat.val
      R.age.val[1:dim(ageMat.val)[1],16:20]                             <- R.age.val[1:dim(ageMat.val)[1],16:20]*y.hat
      R.age.val[(1+dim(ageMat.val)[1]):(2*dim(ageMat.val)[1]),21:25]    <- R.age.val[(1+dim(ageMat.val)[1]):(2*dim(ageMat.val)[1]),21:25]*y.hat
      R.age.val[(1+2*dim(ageMat.val)[1]):(3*dim(ageMat.val)[1]),26:30]  <- R.age.val[(1+2*dim(ageMat.val)[1]):(3*dim(ageMat.val)[1]),26:30]*y.hat
      
      if(y.hat < -85) col.val                          <- rainbow(5)[1] 
      if(y.hat > -85 & y.hat < -30) col.val            <- rainbow(5)[2] 
      if(y.hat > -30 & y.hat < -1) col.val             <- rainbow(5)[3] 
      if(y.hat > -1 & y.hat < 216) col.val             <- rainbow(5)[4]
      if(y.hat > 216) col.val                          <- rainbow(5)[5]
      if (d==1) lines(age.vals+50,y.hat + R.age.val[1:length(age.vals),c(1:5,16:20,31:(30+5*J))] %*% gamma1.hat[c(1:5,16:20,31:(30+5*J))] + cent[d],col=col.val) 
      if (d==2) lines(age.vals+50,y.hat + R.age.val[(length(age.vals)+1):(2*length(age.vals)),c(6:10,21:25,(31+5*J):(30+2*5*J))] %*% gamma1.hat[c(6:10,21:25,(31+5*J):(30+2*5*J))]  + cent[d],col=col.val)
      if (d==3) lines(age.vals+50,y.hat + R.age.val[(2*length(age.vals)+1):(3*length(age.vals)),c(11:15,26:30,(31+2*5*J):(30+3*5*J))] %*% gamma1.hat[c(11:15,26:30,(31+2*5*J):(30+3*5*J))]+ cent[d],col=col.val)
    }
    y.hat                                              <- mean(baseline[,country.time.match$country %in% unique(subset$country)])#### subset$country
    R.age.v                                            <- cbind(ageMat.val,ageMat.val,ageMat.val, ageMat.val,ageMat.val,ageMat.val,
                                                                matrix(0, length(age.vals),5*J),matrix(0, length(age.vals),5*J),
                                                                matrix(0, length(age.vals),5*J))
    R.age.val                                          <- rbind(R.age.v,R.age.v,R.age.v)
    R.age.val[1:dim(ageMat.val)[1],c(6:15,21:30,(J*5+31):(3*J*5+30))]   <- 0
    R.age.val[(1+dim(ageMat.val)[1]):(2*dim(ageMat.val)[1]),c(1:5,11:20,26:(J*5+30),(2*J*5+31):(3*J*5+30))]  <- 0
    R.age.val[(1+2*dim(ageMat.val)[1]):(3*dim(ageMat.val)[1]),c(1:10,16:25,31:(2*J*5+30))] <- 0
    
    R.age.val[1:dim(ageMat.val)[1],16:20]                             <- R.age.val[1:dim(ageMat.val)[1],16:20]*y.hat
    R.age.val[(1+dim(ageMat.val)[1]):(2*dim(ageMat.val)[1]),21:25]    <- R.age.val[(1+dim(ageMat.val)[1]):(2*dim(ageMat.val)[1]),21:25]*y.hat
    R.age.val[(1+2*dim(ageMat.val)[1]):(3*dim(ageMat.val)[1]),26:30]  <- R.age.val[(1+2*dim(ageMat.val)[1]):(3*dim(ageMat.val)[1]),26:30]*y.hat
    
    if (d==1) lines(age.vals + 50, y.hat + R.age.val[1:length(age.vals),dimns[[d]]] %*% gamma1.hat[dimns[[d]]] + cent[d],lwd=3)
    if (d==2) lines(age.vals + 50, y.hat + R.age.val[(length(age.vals)+1):(2*length(age.vals)),dimns[[d]]] %*% gamma1.hat[dimns[[d]]] + cent[d],lwd=3)
    if (d==3) lines(age.vals + 50, y.hat + R.age.val[(2*length(age.vals)+1):(3*length(age.vals)),dimns[[d]]] %*% gamma1.hat[dimns[[d]]] + cent[d],lwd=3)  
    box()
}
dev.off()
####### Population weighted for regions and super-regions #######
### Create population dataset ###
data.country                                  <- read.csv("C:\\Users\\Mariza\\Dropbox\\PhD_Project\\papers_code_data\\Adult_BMI_model54_test\\country_list2015-09-15.csv") 
a <- NULL;for ( i in 1:J) a <- c(a,c(which(data.country$Country==country.match$name[i]))) 
data.country                                  <- data.country[a,]

pop.data                                      <- read.csv("C:\\Users\\Mariza\\Dropbox\\PhD_Project\\papers_code_data\\2nd year\\Postprocessing\\popbyagesex_china_revised.csv") 
country.match$name[9]                         <- "Bolivia (Plurinational State of)"  
country.match$name[21]                        <- "Democratic Republic of the Congo"
country.match$name[32]                        <- "Iran (Islamic Republic of)" 
country.match$name[15]                        <- "China, Hong Kong SAR" 
country.match$name[56]                        <- "United Republic of Tanzania"
a <- NULL;for ( i in 1:J) a <- c(a,c(which(pop.data$country==country.match$name[i])))
pop.data                                      <- pop.data[a,]
pop.data                                      <- subset(pop.data,pop.data$sex=='m')#sex.pop = m or f
pop.data                                      <- pop.data[order(pop.data$country,pop.data$year),]
pop.data[pop.data$country=="Sudan",c(8:30)]   <- pop.data[pop.data$country=="Sudan",c(8:30)]+pop.data[pop.data$country=="South Sudan",c(8:30)]
colnames(pop.data)[1]                         <- "iso"
colnames(pop.data)[6]                         <- "data_year"

pop.data                                      <- merge(pop.data,data.country,by="iso",all.y=T)
pop.data                                      <- subset(pop.data,!is.na(pop.data$Region))
pop.data$pop20to24                            <- pop.data$x20
pop.data$pop25to29                            <- pop.data$x25 
pop.data$pop30to34                            <- pop.data$x30
pop.data$pop35to39                            <- pop.data$x35  
pop.data$pop40to44                            <- pop.data$x40
pop.data$pop45to49                            <- pop.data$x45 
pop.data$pop50to54                            <- pop.data$x50
pop.data$pop55to59                            <- pop.data$x55 
pop.data$pop60to64                            <- pop.data$x60
pop.data$pop65to69                            <- pop.data$x65  
pop.data$pop70to74                            <- pop.data$x70
pop.data$pop75to79                            <- pop.data$x75 
pop.data$pop80to84                            <- pop.data$x80
pop.data$pop85plus                            <- pop.data$x85 + pop.data$x90 + pop.data$x95 + pop.data$x100
pop.data                                      <- pop.data[,c("iso","Country","Region","Superregion","data_year","pop20to24","pop25to29",
                                                             "pop30to34","pop35to39","pop40to44","pop45to49","pop50to54","pop55to59","pop60to64","pop65to69",
                                                             "pop70to74","pop75to79","pop80to84","pop85plus")]
age.bands                                     <- c("20to24","25to29","30to34","35to39","40to44","45to49","50to54","55to59","60to64","65to69",
                                                   "70to74","75to79","80to84","85plus")
age.band.names                                <- apply(as.matrix(age.bands),1,function(y) paste("pop",y,sep=''))
other.regions                                 <- read.csv("C:\\Users\\Mariza\\Dropbox\\PhD_Project\\papers_code_data\\2nd year\\Postprocessing\\reglist_complete.csv",header=T)
other.regions$country[25]                     <- "Bolivia (Plurinational State of)"   
other.regions$country[38]                     <- "Democratic Republic of the Congo"
other.regions$country[85]                     <- "Iran (Islamic Republic of)" 
other.regions$country[77]                     <- "China, Hong Kong SAR" 
other.regions$country[186]                    <- "United Republic of Tanzania"
a<-NULL;for ( i in 1:J) a <- c(a,c(which(other.regions$country==country.match$name[i])))
other.regions                                 <- other.regions[a,]
pop.data                                      <- merge(pop.data,other.regions,by="iso",all.x=T)
region.names                                  <- c("Region","Superregion")
pop.data                                      <- subset(pop.data,pop.data$data_year>=start.year&pop.data$data_year<=end.year)

### Population aggregates ###
# Global #
pop.data.global                               <- aggregate(get(paste("pop",age.bands[1],sep='')) ~ data_year ,FUN=sum,data=pop.data)
for(i in 2:length(age.bands)) pop.data.global <- cbind(pop.data.global, aggregate(get(paste("pop",age.bands[i],sep=''))~ data_year,FUN=sum,data=pop.data)[,2])
colnames(pop.data.global)                     <- c("data_year",age.band.names)
# Regions #
for(k in 1:length(region.names)) {
  assign(paste("pop.data.",region.names[k],sep=''),
         aggregate(get(paste("pop",age.bands[1],sep=''))~get(region.names[k]) + data_year,FUN=sum,data=pop.data))
  for(i in 2:length(age.bands)) {
    assign( paste("pop.data.",region.names[k],sep=''),
            cbind(get(paste("pop.data.",region.names[k],sep='')),
                  aggregate(get(paste("pop",age.bands[i],sep=''))~get(region.names[k])+data_year,FUN=sum,data=pop.data)[,3]) )
  } 
}
for(k in 1:length(region.names)) {
  temp.data                                 <- get(paste("pop.data.",region.names[k],sep=''))
  names(temp.data)                          <- c(region.names[k],"data_year",age.band.names)
  assign(paste("pop.data.",region.names[k],sep=''),temp.data)
}
### Merge population data ###
pop.data$merge.global                       <- pop.data$data_year
for(i in 1:length(region.names)) {
  var.temp                                  <- paste("merge.",region.names[i],sep='')
  pop.data[,c(var.temp)]                    <- paste(pop.data[,region.names[i]],pop.data$data_year)
}
pop.data.global$merge.global                <- pop.data.global$data_year
for(i in 1:length(region.names)) {
  data.temp                                 <- get(paste("pop.data.",region.names[i],sep=''))
  var.temp                                  <- paste("merge.",region.names[i],sep='')
  data.temp[,c(var.temp)]                   <- paste(data.temp[,c(region.names[i])],data.temp[,"data_year"]) 
  assign(paste("pop.data.",region.names[i],sep=''),data.temp)
}
pop.data2                                   <- merge(pop.data,pop.data.global,by="merge.global",all.x=T,suffixes=c(".country",".global"))
for(i in 1:length(region.names)) pop.data2  <- merge(pop.data2,get(paste("pop.data.",region.names[i],sep='')),
                                                     by=paste("merge.",region.names[i],sep=''),all.x=T,suffixes=c("",paste(".",region.names[i],sep='')))
### Create final pop data ###
pop.data2                                   <- pop.data2[,c("iso","Country","Region","Superregion","data_year.country",
                                                            paste(age.band.names,".country",sep=''),
                                                            paste(age.band.names,".global",sep=''),
                                                            paste(age.band.names,"",sep=''),
                                                            paste(age.band.names,".Superregion",sep=''))]
colnames(pop.data2)                         <- c("iso","Country","Region","Superregion","data_year",
                                                 paste(age.band.names,".country",sep=''),
                                                 paste(age.band.names,".global",sep=''),
                                                 paste(age.band.names,".Region",sep=''),
                                                 paste(age.band.names,".Superregion",sep=''))
for(i in 1:length(age.band.names)) {
  country.name.temp                         <- paste(age.band.names[i],".country",sep='')
  for(j in 1:length(region.names)) {
    denominator.temp                        <- paste(age.band.names[i],".",region.names[j],sep='')
    prop.temp                               <- pop.data2[,c(country.name.temp)] / ifelse(pop.data2[,denominator.temp]==0,Inf,pop.data2[,denominator.temp])
    pop.data2[,paste(age.band.names[i],".",region.names[j],"prop",sep='')] <- prop.temp
  }
  denominator.global                        <- paste(age.band.names[i],".global",sep='')
  prop.global                               <- pop.data2[,c(country.name.temp)] / ifelse(pop.data2[,denominator.global]==0,Inf,pop.data2[,denominator.global])
  pop.data2[,paste(age.band.names[i],".globalprop",sep='')]<-prop.global
}
pop.weights                                 <- pop.data2[,c("Country","data_year",paste(age.band.names,".globalprop",sep=''),
                                                            paste(age.band.names,".Regionprop",sep=''),paste(age.band.names,".Superregionprop",sep=''))]
pop.weights$merge                           <- paste(pop.weights$Country,pop.weights$data_year)
covar                                       <- merge(covar,other.regions,by="country",all.x=T)
covar$merge                                 <- paste(covar$country,covar$data_year)
covar                                       <- merge(covar,pop.weights,by="merge",all.x=T)
for(i in min(grep("pop",colnames(covar))):ncol(covar)){ covar[is.na(covar[,i]),i] <- 0 }

#### Weighted results ####

for(k in 1:length(region.names)) {
  for(l in 1:length(age.band.names)) {
    assign(paste(region.names[k],".weighted.age.specific.meanD.",age.bands[l],sep=''),
           as.data.frame(t(Dage.specific.mean[[l]]) * covar[rep(seq(nrow(covar)),each=T),paste(age.band.names[l],".",region.names[k],"prop",sep='')]))#
    assign(paste(region.names[k],".weighted.age.specific.meanS.",age.bands[l],sep=''),
           as.data.frame(t(Sage.specific.mean[[l]]) * covar[rep(seq(nrow(covar)),each=T),paste(age.band.names[l],".",region.names[k],"prop",sep='')]))#
    assign(paste(region.names[k],".weighted.age.specific.meanINT.",age.bands[l],sep=''),
           as.data.frame(t(INTage.specific.mean[[l]]) * covar[rep(seq(nrow(covar)),each=T),paste(age.band.names[l],".",region.names[k],"prop",sep='')]))#
  }	
}
for(l in 1:length(age.band.names)) {
  assign(paste("global.weighted.age.specific.meanD.",age.bands[l],sep=''),
         as.data.frame(t(Dage.specific.mean[[l]]) * covar[rep(seq(nrow(covar)),each=T),paste(age.band.names[l],".globalprop",sep='')]))#
  assign(paste("global.weighted.age.specific.meanS.",age.bands[l],sep=''),
         as.data.frame(t(Sage.specific.mean[[l]]) * covar[rep(seq(nrow(covar)),each=T),paste(age.band.names[l],".globalprop",sep='')]))#
  assign(paste("global.weighted.age.specific.meanINT.",age.bands[l],sep=''),
         as.data.frame(t(INTage.specific.mean[[l]]) * covar[rep(seq(nrow(covar)),each=T),paste(age.band.names[l],".globalprop",sep='')]))#
}
region.temp                                   <- c("region","sregion")
for(k in 1:length(region.temp)) {
  for(l in 1:length(age.band.names)) {
    dat.temp1                                 <- get(paste(region.names[k],".weighted.age.specific.meanD.",age.bands[l],sep=''))
    dat.temp1[,"Region"]                      <- covar[rep(seq(nrow(covar)),each=T),region.temp[k]]#
    assign(paste(region.names[k],".weighted.age.specific.meanD.",age.bands[l],sep=''),dat.temp1)
    dat.temp2                                 <- get(paste(region.names[k],".weighted.age.specific.meanS.",age.bands[l],sep=''))
    dat.temp2[,"Region"]                      <- covar[rep(seq(nrow(covar)),each=T),region.temp[k]]#
    assign(paste(region.names[k],".weighted.age.specific.meanS.",age.bands[l],sep=''),dat.temp2)
    dat.temp3                                 <- get(paste(region.names[k],".weighted.age.specific.meanINT.",age.bands[l],sep=''))
    dat.temp3[,"Region"]                      <- covar[rep(seq(nrow(covar)),each=T),region.temp[k]]#
    assign(paste(region.names[k],".weighted.age.specific.meanINT.",age.bands[l],sep=''),dat.temp3)
  }}
for(k in 1:length(region.temp)) {
  for(l in 1:length(age.band.names)) {
    dat.temp1                               <- get(paste(region.names[k],".weighted.age.specific.meanD.",age.bands[l],sep=''))
    dat.temp1[,"Year"]                      <- covar[rep(seq(nrow(covar)),each=T),"data_year.x"]#covar$data_year.x
    assign(paste(region.names[k],".weighted.age.specific.meanD.",age.bands[l],sep=''),dat.temp1)
    dat.temp2                               <- get(paste(region.names[k],".weighted.age.specific.meanS.",age.bands[l],sep=''))
    dat.temp2[,"Year"]                      <- covar[rep(seq(nrow(covar)),each=T),"data_year.x"]#covar$data_year.x
    assign(paste(region.names[k],".weighted.age.specific.meanS.",age.bands[l],sep=''),dat.temp2)
    dat.temp3                               <- get(paste(region.names[k],".weighted.age.specific.meanINT.",age.bands[l],sep=''))
    dat.temp3[,"Year"]                      <- covar[rep(seq(nrow(covar)),each=T),"data_year.x"]#covar$data_year.x
    assign(paste(region.names[k],".weighted.age.specific.meanINT.",age.bands[l],sep=''),dat.temp3)
  }}
for(l in 1:length(age.band.names)) {	
  dat.temp1                                 <- get(paste("global.weighted.age.specific.meanD.",age.bands[l],sep=''))
  dat.temp1[,"Year"]                        <- rep(covar$data_year.x,each=T)#covar$data_year.x
  assign(paste("global.weighted.age.specific.meanD.",age.bands[l],sep=''),dat.temp1)
  dat.temp2                                 <- get(paste("global.weighted.age.specific.meanS.",age.bands[l],sep=''))
  dat.temp2[,"Year"]                        <- rep(covar$data_year.x,each=T)#covar$data_year.x
  assign(paste("global.weighted.age.specific.meanS.",age.bands[l],sep=''),dat.temp2)
  dat.temp3                                 <- get(paste("global.weighted.age.specific.meanINT.",age.bands[l],sep=''))
  dat.temp3[,"Year"]                        <- rep(covar$data_year.x,each=T)#covar$data_year.x
  assign(paste("global.weighted.age.specific.meanINT.",age.bands[l],sep=''),dat.temp3)
}
##me I add the coding in the years for global, region, s-region
for(l in 1:length(age.band.names)) {
  dat.temp1                <- get(paste("global.weighted.age.specific.meanD.",age.bands[l],sep=''))
  dat.temp1[,"cod.year"]   <- as.factor(rep(t[,3],J)); assign(paste("global.weighted.age.specific.meanD.",age.bands[l],sep=''),dat.temp1)
  dat.temp2                <- get(paste("global.weighted.age.specific.meanS.",age.bands[l],sep=''))
  dat.temp2[,"cod.year"]   <- as.factor(rep(t[,3],J)); assign(paste("global.weighted.age.specific.meanS.",age.bands[l],sep=''),dat.temp2)
  dat.temp3                <- get(paste("global.weighted.age.specific.meanINT.",age.bands[l],sep=''))
  dat.temp3[,"cod.year"]   <- as.factor(rep(t[,3],J)); assign(paste("global.weighted.age.specific.meanINT.",age.bands[l],sep=''),dat.temp3)
}
count <- c(K,L)
for(k in 1:length(region.temp)) {
  for(l in 1:length(age.band.names)) {
    dat.temp1              <- get(paste(region.names[k],".weighted.age.specific.meanD.",age.bands[l],sep=''))
    dat.temp1[,"cod.year"] <- as.factor(rep(t[,3],count[k])); assign(paste(region.names[k],".weighted.age.specific.meanD.",age.bands[l],sep=''),dat.temp1)
    dat.temp2              <- get(paste(region.names[k],".weighted.age.specific.meanS.",age.bands[l],sep=''))
    dat.temp2[,"cod.year"] <- as.factor(rep(t[,3],count[k])); assign(paste(region.names[k],".weighted.age.specific.meanS.",age.bands[l],sep=''),dat.temp2)
    dat.temp3              <- get(paste(region.names[k],".weighted.age.specific.meanINT.",age.bands[l],sep=''))
    dat.temp3[,"cod.year"] <- as.factor(rep(t[,3],count[k])); assign(paste(region.names[k],".weighted.age.specific.meanINT.",age.bands[l],sep=''),dat.temp3)
  }}  
##me

for(l in 1:length(age.band.names)) {
  assign(paste("global.weightedD.",age.bands[l],sep=''),
         summaryBy(.~cod.year,data=get(paste("global.weighted.age.specific.meanD.",age.bands[l],sep='')),FUN=sum))#Year
  assign(paste("global.weightedS.",age.bands[l],sep=''),
         summaryBy(.~cod.year,data=get(paste("global.weighted.age.specific.meanS.",age.bands[l],sep='')),FUN=sum))#Year
  assign(paste("global.weightedINT.",age.bands[l],sep=''),
         summaryBy(.~cod.year,data=get(paste("global.weighted.age.specific.meanINT.",age.bands[l],sep='')),FUN=sum))#Year
}
for(k in 1:length(region.temp)) {
  for(l in 1:length(age.band.names)) {
    assign(paste(region.names[k],".weightedD.",age.bands[l],sep=''),
           summaryBy(.~ Region + cod.year, data=get(paste(region.names[k],".weighted.age.specific.meanD.",age.bands[l],sep='')),FUN=sum))#Year
    assign(paste(region.names[k],".weightedS.",age.bands[l],sep=''),
           summaryBy(.~Region + cod.year, data=get(paste(region.names[k],".weighted.age.specific.meanS.",age.bands[l],sep='')),FUN=sum))#Year
    assign(paste(region.names[k],".weightedINT.",age.bands[l],sep=''),
           summaryBy(.~Region + cod.year, data=get(paste(region.names[k],".weighted.age.specific.meanINT.",age.bands[l],sep='')),FUN=sum))#Year
  }}
#### Age-standardised results ####
global.age.standardD         <- (get(paste("global.weightedD.",age.bands[1],sep=''))[,-c(1,203)] * StandardPop$weights[1])#both years out
global.age.standardS         <- (get(paste("global.weightedS.",age.bands[1],sep=''))[,-c(1,203)] * StandardPop$weights[1])#both years out
global.age.standardINT       <- (get(paste("global.weightedINT.",age.bands[1],sep=''))[,-c(1,203)] * StandardPop$weights[1])#both years out
for(i in 2:length(age.band.names)) {
  global.age.standardD       <- global.age.standardD + (get(paste("global.weightedD.",age.bands[i],sep=''))[,-c(1,203)] * StandardPop$weights[i])
  global.age.standardS       <- global.age.standardS + (get(paste("global.weightedS.",age.bands[i],sep=''))[,-c(1,203)] * StandardPop$weights[i])
  global.age.standardINT     <- global.age.standardINT + (get(paste("global.weightedINT.",age.bands[i],sep=''))[,-c(1,203)] * StandardPop$weights[i])
}
for(k in 1:length(region.names)) {
  assign(paste(region.names[k],".age.standardD",sep=''),(get(paste(region.names[k],".weightedD.",age.bands[1],sep=''))[,-c(1:2)] * StandardPop$weights[1]))
  assign(paste(region.names[k],".age.standardS",sep=''),(get(paste(region.names[k],".weightedS.",age.bands[1],sep=''))[,-c(1:2)] * StandardPop$weights[1]))
  assign(paste(region.names[k],".age.standardINT",sep=''),(get(paste(region.names[k],".weightedINT.",age.bands[1],sep=''))[,-c(1:2)] * StandardPop$weights[1]))
}
for(k in 1:length(region.names)) {
  weighted.tempD              <- get(paste(region.names[k],".age.standardD",sep=''))
  weighted.tempS              <- get(paste(region.names[k],".age.standardS",sep=''))
  weighted.tempINT            <- get(paste(region.names[k],".age.standardINT",sep=''))
  
  for(i in 2:length(age.band.names)) {
    weighted.tempD            <- weighted.tempD +
      (get(paste(region.names[k],".weightedD.",age.bands[i],sep=''))[,-c(1:2)] * StandardPop$weights[i])
    weighted.tempS            <- weighted.tempS +
      (get(paste(region.names[k],".weightedS.",age.bands[i],sep=''))[,-c(1:2)] * StandardPop$weights[i])
    weighted.tempINT          <- weighted.tempINT +
      (get(paste(region.names[k],".weightedINT.",age.bands[i],sep=''))[,-c(1:2)] * StandardPop$weights[i])
  }
  assign(paste(region.names[k],".age.standardD",sep=''),weighted.tempD)
  assign(paste(region.names[k],".age.standardS",sep=''),weighted.tempS)
  assign(paste(region.names[k],".age.standardINT",sep=''),weighted.tempINT)
}
global.age.standardD          <- global.age.standardD + cent[1]  
global.age.standardS          <- global.age.standardS + cent[2]
global.age.standardINT        <- global.age.standardINT + cent[3] 
for(k in 1:length(region.names)) {
  weighted.temp1               <- get(paste(region.names[k],".age.standardD",sep='')) + cent[1] 
  weighted.temp2               <- get(paste(region.names[k],".age.standardS",sep='')) + cent[2] 
  weighted.temp3               <- get(paste(region.names[k],".age.standardINT",sep='')) + cent[3] 
  assign(paste(region.names[k],".age.standardD",sep=''),weighted.temp1)
  assign(paste(region.names[k],".age.standardS",sep=''),weighted.temp2)
  assign(paste(region.names[k],".age.standardINT",sep=''),weighted.temp3)
}
### Output summary results ###
# Global #
global.age.standardD.mean      <- rowMeans(global.age.standardD)
global.age.standardS.mean      <- rowMeans(global.age.standardS)
global.age.standardINT.mean    <- rowMeans(global.age.standardINT)
global.age.standardD.l         <- apply(global.age.standardD, 1, quantile, .025)
global.age.standardS.l         <- apply(global.age.standardS, 1, quantile, .025)
global.age.standardINT.l       <- apply(global.age.standardINT, 1, quantile, .025)
global.age.standardD.u         <- apply(global.age.standardD, 1, quantile, .975)
global.age.standardS.u         <- apply(global.age.standardS, 1, quantile, .975)
global.age.standardINT.u       <- apply(global.age.standardINT, 1, quantile, .975)
global.age.standardD.se        <- apply(global.age.standardD, 1, sd)
global.age.standardS.se        <- apply(global.age.standardS, 1, sd)
global.age.standardINT.se      <- apply(global.age.standardINT, 1, sd)
global.means                   <- data.frame(rep(start.year:end.year,each=T), global.weightedD.20to24$cod.year,
                                             global.age.standardD.mean,global.age.standardD.l,global.age.standardD.u,
                                             global.age.standardD.se,
                                             global.age.standardS.mean,global.age.standardS.l,global.age.standardS.u,
                                             global.age.standardS.se,
                                             global.age.standardINT.mean,global.age.standardINT.l,global.age.standardINT.u,
                                             global.age.standardINT.se)
names(global.means)            <- c("year","cod.year", "meanD", "l.D", "u.D", "se.D", "meanS", "l.S", "u.S", "se.S", "meanINT", "l.INT", "u.INT", "se.INT")
write.csv(global.means, paste(filename, "_GlobalMeans.csv", sep=''),row.names=FALSE)
# Regional #
for(k in 1:length(region.names)) {
  regional.temp                <- get(paste(region.names[k],".age.standardD",sep=''))
  regional.temp                <- get(paste(region.names[k],".age.standardS",sep=''))
  regional.temp                <- get(paste(region.names[k],".age.standardINT",sep=''))
  assign(paste(region.names[k],".age.standardD.mean",sep=''), rowMeans(regional.temp))
  assign(paste(region.names[k],".age.standardS.mean",sep=''), rowMeans(regional.temp))
  assign(paste(region.names[k],".age.standardINT.mean",sep=''), rowMeans(regional.temp))
  assign(paste(region.names[k],".age.standardD.l",sep=''),    apply(regional.temp, 1, quantile, .025))
  assign(paste(region.names[k],".age.standardS.l",sep=''),    apply(regional.temp, 1, quantile, .025))
  assign(paste(region.names[k],".age.standardINT.l",sep=''),  apply(regional.temp, 1, quantile, .025))
  assign(paste(region.names[k],".age.standardD.u",sep=''),    apply(regional.temp, 1, quantile, .975))
  assign(paste(region.names[k],".age.standardS.u",sep=''),    apply(regional.temp, 1, quantile, .975))
  assign(paste(region.names[k],".age.standardINT.u",sep=''),  apply(regional.temp, 1, quantile, .975))
  assign(paste(region.names[k],".age.standardD.se",sep=''),   apply(regional.temp, 1, sd))
  assign(paste(region.names[k],".age.standardS.se",sep=''),   apply(regional.temp, 1, sd))
  assign(paste(region.names[k],".age.standardINT.se",sep=''), apply(regional.temp, 1, sd))
  
  assign(paste(region.names[k],".means",sep=''),
         data.frame(get(paste(region.names[k],".weightedD.20to24",sep=''))[,"Region"],get(paste(region.names[k],".weightedD.20to24",sep=''))[,"cod.year"],#Year
                    get(paste(region.names[k],".age.standardD.mean",sep='')),get(paste(region.names[k],".age.standardS.mean",sep='')),get(paste(region.names[k],".age.standardINT.mean",sep='')),
                    get(paste(region.names[k],".age.standardD.l",sep='')),get(paste(region.names[k],".age.standardS.l",sep='')),get(paste(region.names[k],".age.standardINT.l",sep='')),
                    get(paste(region.names[k],".age.standardD.u",sep='')),get(paste(region.names[k],".age.standardS.u",sep='')),get(paste(region.names[k],".age.standardINT.u",sep='')),
                    get(paste(region.names[k],".age.standardD.se",sep='')),get(paste(region.names[k],".age.standardS.se",sep='')),get(paste(region.names[k],".age.standardINT.se",sep='')) ))
}
for(k in 1:length(region.names)) {
  regional.temp                 <- get(paste(region.names[k],".means",sep=''))
  names(regional.temp)          <- c("region", "year", "meanD", "l.D", "u.D", "se.D", "meanS", "l.S", "u.S", "se.S"
                                     ,"meanINT", "l.INT", "u.INT", "se.INT")
  regional.temp                 <- subset(regional.temp,regional.temp$region!="NO REGION")
  write.csv(regional.temp,paste(filename,"_", region.names[k],"_Means.csv", sep=''),row.names=FALSE)
  assign(paste(region.names[k],".means",sep=''),regional.temp)
}
#### Slopes ####
# Global #
slope.valsD <- slope.valsS      <- slope.valsDINT <- rep(NA, length(burnt))
global.slopes                   <- as.data.frame(matrix(NA,1,16))
names(global.slopes)            <- c("region", "meanD.all", "meanS.all", "meanINT.all","l_D.all","l_S.all","l_INT.all","u_D.all",
                                     "u_S.all","u_INT.all", "se_D.all","se_S.all","se_INT.all", "p_D.all", "p_S.all", "p_INT.all")
global.slopes[,1]               <- "Global"
slope.valsD 	                  <- lm(as.matrix(global.age.standardD) ~ t[,3])$coef[2,]#start.year:end.year
slope.valsS                     <- lm(as.matrix(global.age.standardS) ~ t[,3])$coef[2,]#start.year:end.year
slope.valsINT	                  <- lm(as.matrix(global.age.standardINT) ~ t[,3])$coef[2,]#start.year:end.year

global.slopes[1,2]              <- mean(slope.valsD)
global.slopes[1,3]              <- mean(slope.valsS)
global.slopes[1,4]              <- mean(slope.valsINT)

global.slopes[1,5] 	            <- quantile(slope.valsD, .025)
global.slopes[1,6] 	            <- quantile(slope.valsS, .025)
global.slopes[1,7] 	            <- quantile(slope.valsINT, .025)

global.slopes[1,8] 	            <- quantile(slope.valsD, .975)
global.slopes[1,9] 	            <- quantile(slope.valsS, .975)
global.slopes[1,10] 	          <- quantile(slope.valsINT, .975)

global.slopes[1,11] 	          <- sd(slope.valsD)
global.slopes[1,12] 	          <- sd(slope.valsS)
global.slopes[1,13] 	          <- sd(slope.valsINT)

global.slopes[1,14] 	          <- mean(sign(slope.valsD)!=sign(mean(slope.valsD)))
global.slopes[1,15] 	          <- mean(sign(slope.valsS)!=sign(mean(slope.valsS)))
global.slopes[1,16] 	          <- mean(sign(slope.valsINT)!=sign(mean(slope.valsINT)))

write.csv(global.slopes, file=paste(filename,"_GlobalSlopes.csv", sep=''),row.names=FALSE)
### Regions ###
for(k in 1:length(region.names)) {
  slope.valsD <- slope.valsS      <- slope.valsINT <- rep(NA, length(burnt))
  region.list                     <- unique(get(paste(region.names[k],".weightedD.20to24",sep=''))[,"Region"])
  region.list                     <- as.character(subset(region.list,region.list!="NO REGION"))
  number.regions                  <- length(region.list)
  assign(paste(region.names[k],".slopes",sep=''),as.data.frame(matrix(NA,number.regions,6)))
  results.mat                     <- get(paste(region.names[k],".slopes",sep=''))
  results.mat[,1]                 <- region.list
  for(i in 1:number.regions) {
    dat.subD                      <- subset(get(paste(region.names[k],".age.standardD",sep='')),get(paste(region.names[k],".weightedD.20to24",sep=''))[,"Region"]==region.list[i])
    dat.subS                      <- subset(get(paste(region.names[k],".age.standardS",sep='')),get(paste(region.names[k],".weightedS.20to24",sep=''))[,"Region"]==region.list[i])
    dat.subINT                    <- subset(get(paste(region.names[k],".age.standardINT",sep='')),get(paste(region.names[k],".weightedINT.20to24",sep=''))[,"Region"]==region.list[i])
    
    slope.valsD                   <- lm(as.matrix(dat.subD) ~ t[,3])$coef[2,]#start.year:end.year
    slope.valsS                   <- lm(as.matrix(dat.subS) ~ t[,3])$coef[2,]#start.year:end.year
    slope.valsINT                 <- lm(as.matrix(dat.subINT) ~ t[,3])$coef[2,]#start.year:end.year
    
    results.mat[i,2] 	            <- mean(slope.valsD)
    results.mat[i,3] 	            <- mean(slope.valsS)
    results.mat[i,4] 	            <- mean(slope.valsINT)
    
    results.mat[i,5]          	  <- quantile(slope.valsD, .025)
    results.mat[i,6]          	  <- quantile(slope.valsS, .025)
    results.mat[i,7]          	  <- quantile(slope.valsINT, .025)
    
    results.mat[i,8]          	  <- quantile(slope.valsD, .975)
    results.mat[i,9]          	  <- quantile(slope.valsS, .975)
    results.mat[i,10]          	  <- quantile(slope.valsINT, .975)
    
    results.mat[i,11] 	          <- sd(slope.valsD)
    results.mat[i,12] 	          <- sd(slope.valsS)
    results.mat[i,13] 	          <- sd(slope.valsINT)
    
    results.mat[i,14]	            <- mean(sign(slope.valsD)!=sign(mean(slope.valsD)))
    results.mat[i,15]	            <- mean(sign(slope.valsS)!=sign(mean(slope.valsS)))
    results.mat[i,16]	            <- mean(sign(slope.valsINT)!=sign(mean(slope.valsINT)))
  }
  names(results.mat)              <- c("region","meanD.all","meanS.all","meanINT.all","l_D.all","l_S.all","l_INT.all","u_D.all",
                                       "u_S.all","u_INT.all","se_D.all","se_S.all","se_INT.all", "p_D.all", "p_S.all", "p_INT.all")
  write.csv(results.mat,paste(filename,"_", region.names[k],"_Slopes.csv", sep=''),row.names=FALSE)
  assign(paste(region.names[k],".slopes",sep=''),results.mat)
}
######## Slope changes ###########
country.slopeD1 <- country.slopeS1 <- country.slopeINT1 <- as.data.frame(matrix(NA, J,length(burnt)+1))
country.slopeD2 <- country.slopeS2 <- country.slopeINT2 <- as.data.frame(matrix(NA, J,length(burnt)))
country.slopeD1[,1] <- country.slopeS1[,1]   <- country.slopeINT1[,1]<-as.character(country.match$name)
country.slopeD2[,1] <- country.slopeS2[,1]   <- country.slopeINT2[,1]<-as.character(country.match$name)
start.slope1 <- start.year;     end.slope1   <- centring.year  # 1) start.year, 2) 2003
start.slope2 <- centring.year+1;end.slope2   <- end.year #1) 2004 2)end.year
for(i in 1:J){
  country.val                                <- as.character(country.match$name[i])
  country.resultsD1                          <- as.data.frame(t(age.standardized.mean[,1:(J*(T^2))][,country.time.match$country==country.val &
                                                                                                      country.time.match$yearD>=start.slope1 & country.time.match$yearD<=end.slope1]+cent[1]))#
  country.resultsS1                          <- as.data.frame(t(age.standardized.mean[,(1+J*(T^2)):(2*J*(T^2))][,country.time.match$country==country.val &
                                                                                                                  country.time.match$yearD>=start.slope1 & country.time.match$yearD<=end.slope1]+cent[2]))#
  country.resultsINT1                        <- as.data.frame(t(age.standardized.mean[,(1+2*J*(T^2)):(3*J*(T^2))][,country.time.match$country==country.val &
                                                                                                                    country.time.match$yearD>=start.slope1 & country.time.match$yearD<=end.slope1]+cent[3]))#
  
  country.resultsD2                          <- as.data.frame(t(age.standardized.mean[,1:(J*(T^2))][,country.time.match$country==country.val &
                                                                                                      country.time.match$yearD>=start.slope2 & country.time.match$yearD<=end.slope2]+cent[1]))#
  country.resultsS2                          <- as.data.frame(t(age.standardized.mean[,(1+J*(T^2)):(2*J*(T^2))][,country.time.match$country==country.val &
                                                                                                                  country.time.match$yearD>=start.slope2 & country.time.match$yearD<=end.slope2]+cent[2]))#
  country.resultsINT2                        <- as.data.frame(t(age.standardized.mean[,(1+2*J*(T^2)):(3*J*(T^2))][,country.time.match$country==country.val &
                                                                                                                    country.time.match$yearD>=start.slope2 & country.time.match$yearD<=end.slope2]+cent[3]))#
  
  year1                                      <- t[1:((T-round(T/2))*T),3] #start.slope1:end.slope1. I will use the coding of the years instead
  year2	                                     <- t[(1+(T-round(T/2))*T):(T^2),3]
  country.slopeD1[i,2:(length(burnt)+1)]     <- as.numeric(lm(t(t(country.resultsD1))~year1)$coef[2,])
  country.slopeS1[i,2:(length(burnt)+1)]     <- as.numeric(lm(t(t(country.resultsS1))~year1)$coef[2,])
  country.slopeINT1[i,2:(length(burnt)+1)]   <- as.numeric(lm(t(t(country.resultsINT1))~year1)$coef[2,])
  
  country.slopeD2[i,2:(length(burnt)+1)]     <- as.numeric(lm(t(t(country.resultsD2))~year2)$coef[2,])
  country.slopeS2[i,2:(length(burnt)+1)]     <- as.numeric(lm(t(t(country.resultsS2))~year2)$coef[2,])
  country.slopeINT2[i,2:(length(burnt)+1)]   <- as.numeric(lm(t(t(country.resultsINT2))~year2)$coef[2,])
}
country.slope.sign1 <- country.slope.sign2   <- list()
country.slope1                               <- list(country.slopeD1,country.slopeS1,country.slopeINT1)
country.slope2                               <- list(country.slopeD2,country.slopeS2,country.slopeINT2)
country.slope.equals 		                     <- list()
country.slope.equals.mean                    <- list(matrix(NA, J, 2),matrix(NA, J, 2),matrix(NA, J, 2))	
names(country.slope.equals.mean)             <- title
for (d in 1:3){
  country.slope.sign1[[d]]  	               <- sign(country.slope1[[d]][,-1])
  country.slope.sign2[[d]]             	     <- sign(country.slope2[[d]][,-1])
  country.slope.equals[[d]]                  <- country.slope.sign1[[d]]!=country.slope.sign2[[d]]
  country.slope.equals.mean[[d]][,1]	       <- as.character(country.match$name)
  country.slope.equals.mean[[d]][,2]	       <- rowMeans(country.slope.equals[[d]])
  country.slope.equals.mean[[d]] 	           <- as.data.frame(country.slope.equals.mean[[d]])
  names(country.slope.equals.mean[[d]])      <- c("Country","Probability")
}
write.csv(country.slope.equals.mean,file=paste(filename,"_Country_slope-change.csv",sep=''),row.names=FALSE)
countries                                    <- read.csv("C:\\Users\\mariz\\Dropbox\\PhD_Project\\papers_code_data\\Adult_BMI_model54_test\\country_list2015-09-15.csv")
country.slope.mean1 <- country.slope.mean2   <- list(as.vector(1:J),as.vector(1:J),as.vector(1:J))
global.slope.mean1  <- global.slope.mean2    <- list(1,1,1) 
global.slope1       <- global.slope2         <- list(as.data.frame(matrix(NA,1,length(burnt))),as.data.frame(matrix(NA,1,length(burnt))),as.data.frame(matrix(NA,1,length(burnt))))
country.slope.mean  <- global.slope.mean     <- list(matrix(NA,1,3),matrix(NA,1,3),matrix(NA,1,3))
reg.temp                                     <- list(NULL,NULL,NULL)
global.age.standard                          <- list(global.age.standardD,global.age.standardS,global.age.standardINT)
for ( d in 1:3){
  country.slope.mean1[[d]]                   <- rowMeans(country.slope1[[d]][,-c(1)])
  country.slope.mean2[[d]]                   <- rowMeans(country.slope2[[d]][,-c(1)])
  country.slope.mean[[d]]                    <- as.data.frame(cbind(country.slope1[[d]][,1], country.slope.mean1[[d]], country.slope.mean2[[d]]))
  names(country.slope.mean[[d]])             <- c("Country","Early","Late")
  country.slope.mean[[d]]                    <- merge(country.slope.mean[[d]],countries,by="Country")
  global.slope1[[d]][,1]<-global.slope2[[d]][,1] <- "World"
  global.results1                            <- as.data.frame(t(global.age.standard[[d]])[,1:((T-round(T/2))*T)])# 1:15 for T^2=25
  global.results2                            <- as.data.frame(t(global.age.standard[[d]])[,(1+(T-round(T/2))*T):(T^2)]) #16:25
  global.slope1[[d]][,2:((length(burnt)+1))] <- as.numeric(lm((t(global.results1))~year1)$coef[2,])
  global.slope2[[d]][,2:((length(burnt)+1))] <- as.numeric(lm((t(global.results2))~year2)$coef[2,])
  global.slope.mean1[[d]]                    <- rowMeans(global.slope1[[d]][,-c(1)])
  global.slope.mean2[[d]]                    <- rowMeans(global.slope2[[d]][,-c(1)])
  global.slope.mean[[d]]                     <- as.data.frame(cbind(global.slope1[[d]][,1],global.slope.mean1[[d]],global.slope.mean2[[d]]))
  names(global.slope.mean[[d]])              <- c("Country","Early","Late")
  reg.temp[[d]]                              <- unique(country.slope.mean[[d]][,5:6])
  reg.temp[[d]]                              <- reg.temp[[d]][order(reg.temp[[d]][,2]),]
  reg.order                                  <- c(19,18,20,21,3,4,17,5,6,7,15,16,11,13,14,12,9,10,8,1,2)
  reg.temp[[d]]                              <- reg.temp[[d]][c(reg.order),]
  reg.temp[[d]]                              <- cbind(reg.temp[[d]],c(15,16,17,18,15,16,15,15,16,15,15,16,15,16,17,18,15,16,17,15,16))
  names(reg.temp[[d]])[3]                    <- "pch"
}
palette(brewer.pal(12,"Paired")) ;dev.off()
palette(c(palette()[c(2,6,8,12,1,4,10,7,5)]));dev.off()
pdf(paste(filename, "_SlopeChangePlot.pdf", sep=''), width=9, height=9)
par(mar = c(5,4.5,4,2) + 0.1)
for ( d in 1:3){
  plot(as.numeric(as.matrix(country.slope.mean[[d]]$Early))*10,as.numeric(as.matrix(country.slope.mean[[d]]$Late))*10,
       xlab=expression("Change in age-standardised BP before 2000 (kg/"*m^2*" per decade)"),ylab= title[d],yaxt="n",
       col=as.numeric(as.factor(country.slope.mean[[d]]$Superregion)),pch=20,xlim=c(-10,13),ylim=c(-10,13.5))
  abline(0,1,lty=2,col="lightgrey")
  abline(0,0,lty=2,col="lightgrey")
  abline(v=0,lty=2,col="lightgrey")
  text(-10,12,"F",cex=.8)
  text(10,13.5,"A",cex=.8)
  text(12.5,4,"B",cex=.8)
  text(12.5,-4,"C",cex=.8)
  text(-1,-10,"D",cex=.8)
  text(-10,-1,"E",cex=.8)
  points(as.numeric(as.matrix(global.slope.mean[[d]][,"Early"]))*10, as.numeric(as.matrix(global.slope.mean[[d]][,"Late"]))*10, pch=20,col="black",cex=2)
  axis(2,labels=TRUE)
  text(-.95,2.5,"Women")
}
legend("bottomright",c(unique(as.character(reg.temp$Superregion)),"World"),col=c(as.numeric(unique(reg.temp$Superregion)),"black"),pch=20,cex=.9,pt.cex=c(rep(.9,9),1.5),box.lty=0)
dev.off()
write.csv(country.slope.mean,file=paste(filename, "_SlopeChanges.csv", sep=''),row.names=FALSE)
######## Ranking means ############
CountryNames                            <- data.frame(unique(covar$country));names(CountryNames)<-c("country")
region.names                            <- unique(covar$region)
covar2                                  <- subset(covar, covar$data_year.x==start.year)
covar2                                  <- covar2[,c("country","region","sregion")]
colnames(covar2)                        <- c("covar.country","covar.region","covar.sregion")
covar2$covar.country                    <- as.character(covar2$covar.country)
#for (i in 1:3) {countr_means[[i]]$yearS <- rep(t[,3],21); names(countr_means[[i]])[3]<-"cod.year"}
foo        <- foo.y                     <- countr_means  # I will need to add the coding year(-12 to 12)
names(foo) <- names(foo.y)              <- title
reg.temp                                <- list(NULL,NULL,NULL)
for (d in 1:3){
  for(year.val in c(min(t[,3]),max(t[,3])) ) { #start.year,end.year
    foo.y[[d]]                            <- foo[[d]][foo[[d]][,"cod.year"]==year.val,]#year
    foo.y[[d]]                            <- foo.y[[d]][,c(1,4,5,6)]#c(1,3,4,5) as we had only 1 column about year
    names(foo.y[[d]])                     <- c("f.country", "f.mean", "f.l", "f.u")
    foo.y[[d]][,c(3,4)]                   <- foo.y[[d]][,c(3,4)]
    foo.y[[d]]                            <- merge(covar2, foo.y[[d]], by.x="covar.country",by.y="f.country")
    foo.y[[d]]                            <- merge(foo.y[[d]], CountryNames, by.x="covar.country", by.y="country")
    foo.y[[d]]                            <- rbind(foo.y[[d]], rep(NA, dim(foo.y[[d]])[2]))
    for(i in 1:3)       foo.y[[d]][,i]    <- as.character(foo.y[[d]][,i])
    foo.y[[d]][nrow(foo.y[[d]]),1:3]      <-"World"
    foo.y[[d]]$f.mean[dim(foo.y[[d]])[1]] <- global.means[,((d-1)*4+3)][global.means$cod.year==year.val]#$mean
    foo.y[[d]]$f.l[dim(foo.y[[d]])[1]]    <- global.means[,((d-1)*4+4)][global.means$cod.year==year.val] #$l
    foo.y[[d]]$f.u[dim(foo.y[[d]])[1]]    <- global.means[,((d-1)*4+5)][global.means$cod.year==year.val]#$u
    Jnew                                  <- dim(foo.y[[d]])[1]
    order                                 <- order(foo.y[[d]]$f.mean)
    reg.temp[[d]]                         <- unique(foo.y[[d]][,2:3])
    reg.temp[[d]]                         <- reg.temp[[d]][order(reg.temp[[d]][,2]),]
    length.reg                            <- as.numeric(table(reg.temp[[d]][,2]))
    pch.temp                              <- NA
    for(i in 1:length(length.reg)) pch.temp <- c(pch.temp,c(15:19)[1:length.reg[i]])
    
    reg.temp[[d]]                         <- cbind(reg.temp[[d]],pch.temp[-1])
    pch.use                               <- reg.temp[[d]][match(foo.y[[d]]$covar.region,reg.temp[[d]]$covar.region),3]
    pch.use[J+2]                          <- 20
    palette(c(palette(),"black"));dev.off()
    pdf(paste(d,filename, year.val, "Means.pdf", sep=''), height=20, width=6)
    par(mar=c(3.6,9,1.1,1.1))
    plot(range(c(foo.y[[d]]$f.l,foo.y[[d]]$f.u)), c(1+6.25,Jnew-6.25), yaxt='n', xlab='', pch='',ylab='', cex.axis=.5, cex.lab=.5, xaxt='n')
    for(odd in seq(1,Jnew,by=2)) polygon(c(min(foo.y[[d]]$f.l),min(foo.y[[d]]$f.l),max(foo.y[[d]]$f.u),max(foo.y[[d]]$f.u),min(foo.y[[d]]$f.l)), c(odd-.5,odd+.5,odd+.5,odd-.5,odd-.5),col=grey(.95), border=NA)
    foo3                                  <- axis(1, cex.axis=.5, cex.lab=.5)
    abline(v=foo3, col=grey(.8))
    axis(2, at=1:Jnew, labels=covar2$covar.country[order], las=2, cex.axis=.5)
    for(i in 1:Jnew) {
      lines(c(foo.y[[d]]$f.l[order][i],foo.y[[d]]$f.u[order][i]), c(i,i), col=as.numeric(as.factor(foo.y[[d]]$covar.sregion))[order][i])
      points(foo.y[[d]]$f.mean[order][i], i, col=as.numeric(as.factor(foo.y[[d]]$covar.sregion))[order][i],
             pch=pch.use[order][i], cex=c(rep(.5,201),1)[order][i])
    }
    if(sex=="female") sex.space <- " female " else if(sex=="male") sex.space <- " male "
    title(xlab=bquote(.(year.val)*.(sex.space)*.(toupper(title[d]))*" means and uncertainty intervals (kg/"*m^2*")"), line=2, cex.lab=.5)
    box()
    dev.off()
  }} 
######## Ranking slopes ##########
foo                                 <- country.slopes
foo.y                               <- list()
for ( d in 1:3){
  foo.y[[d]]                        <- foo[[d]][,c(1,2,3,4)]
  names(foo.y[[d]])                 <- c("f.country", "f.mean", "f.l", "f.u")
  foo.y[[d]][,c(2,3,4)]             <- foo.y[[d]][,c(2,3,4)]*10
  foo.y[[d]]                        <- merge(covar2, foo.y[[d]], by.x="covar.country", by.y="f.country") 
  foo.y[[d]]                        <- merge(foo.y[[d]], CountryNames, by.x="covar.country", by.y="country")
  foo.y[[d]]                        <- rbind(foo.y[[d]], rep(NA, dim(foo.y[[d]])[2]))
  for(i in 1:3) foo.y[[d]][,i]      <- as.character(foo.y[[d]][,i])
  foo.y[[d]][nrow(foo.y[[d]]),c(1:3)]<- "World"
  foo.y[[d]]$f.mean[dim(foo.y[[d]])[1]] <- global.slopes[,1+d]*10  #"mean.all"
  foo.y[[d]]$f.l[dim(foo.y[[d]])[1]]<- global.slopes[,4+d]*10  #"l.all"
  foo.y[[d]]$f.u[dim(foo.y[[d]])[1]]<- global.slopes[,7+d]*10  #"u.all"
  reg.temp[[d]]                     <- unique(foo.y[[d]][,2:3])
  reg.temp[[d]]                     <- reg.temp[[d]][order(reg.temp[[d]][,2]),]
  length.reg                        <- as.numeric(table(reg.temp[[d]][,2]))
  pch.temp                          <- NA
  for(i in 1:length(length.reg))  pch.temp <- c(pch.temp,c(15:19)[1:length.reg[i]])
  reg.temp[[d]]                     <- cbind(reg.temp[[d]],pch.temp[-1])
  pch.use                           <- reg.temp[[d]][match(foo.y[[d]]$covar.region,reg.temp[[d]]$covar.region),3]
  pch.use[J+2]                      <- 20
  Jnew                              <- dim(foo.y[[d]])[1]
  order                             <- order(foo.y[[d]]$f.mean)
  pdf(paste(d,filename, "Slopes.pdf", sep=''), height=20, width=6)
  par(mar=c(3.6,9,1.1,1.1))
  plot(range(c(foo.y[[d]]$f.l,foo.y[[d]]$f.u)), c(1+6.25,Jnew-6.25), yaxt='n', xlab='', pch='',ylab='', cex.axis=.5, cex.lab=.5, xaxt='n')
  for(odd in seq(1,Jnew,by=2)) polygon(c(min(foo.y[[d]]$f.l),min(foo.y[[d]]$f.l),max(foo.y[[d]]$f.u),max(foo.y[[d]]$f.u),min(foo.y[[d]]$f.l)),c(odd-.5,odd+.5,odd+.5,odd-.5,odd-.5), col=grey(.95), border=NA)
  foo3                              <- axis(1, cex.axis=.5, cex.lab=.5)
  abline(v=foo3, col=grey(.8))
  axis(2, at=1:Jnew, labels=covar2$covar.country[order], las=2, cex.axis=.5)
  for(i in 1:Jnew){
    lines(c(foo.y[[d]]$f.l[order][i], foo.y[[d]]$f.u[order][i]), c(i,i), col=as.numeric(as.factor(foo.y$covar.sregion))[order][i])
    points(foo.y[[d]]$f.mean[order][i], i, col=as.numeric(as.factor(foo.y[[d]]$covar.sregion))[order][i], pch=pch.use[order][i], cex=c(rep(.5,201),1)[order][i])
  }
  title(xlab=bquote("Change in age-standardised "*.(tolower(sex))*" "*.(toupper(title[d]))*" (kg/"*m^2*"/decade)"), line=2, cex.lab=.5)
  box()
  dev.off()
}
####### Country means map #####

#rm(aa,drawsD1,drawsD10,drawsD11,drawsD12,drawsD13,drawsD14) 
#rm(aa,drawsD2,drawsD3,drawsD4,drawsD5,drawsD6,drawsD7,draws.nameD9)
map.year                                   <- 0#2015#2016 i ll use the codding now

country.means.mapD                         <- subset(countr_means[[1]],countr_means[[1]]$cod.year==map.year)[,c(1,3)]
country.means.mapS                         <- subset(countr_means[[2]],countr_means[[2]]$cod.year==map.year)[,3]
country.means.mapINT                       <- subset(countr_means[[3]],countr_means[[3]]$cod.year==map.year)[,3]
country.means.map                          <- data.frame("id" = NA,"names" = country.means.mapD[,1],"mean_D" = country.means.mapD[,2],
                                                         "mean_S"=country.means.mapS,"mean_INT"=country.means.mapINT)

levels(country.means.map$names)[122]       <- "Burma"
levels(country.means.map$names)[166]       <- "Korea, Republic of"
levels(country.means.map$names)[132]       <- "Korea, Democratic People's Republic of"
levels(country.means.map$names)[102]       <- "Libyan Arab Jamahiriya"
levels(country.means.map$names)[148]       <- "Russia"
levels(country.means.map$names)[134]       <- "Palestine"
levels(country.means.map$names)[192]       <- "United States"
levels(country.means.map$names)[54]        <- "Democratic Republic of the Congo"
levels(country.means.map$names)[177]       <- "United Republic of Tanzania"
levels(country.means.map$names)[84]        <- "Iran (Islamic Republic of)"
levels(country.means.map$names)[97]        <- "Lao People's Democratic Republic"
levels(country.means.map$names)[76]        <- "Guinea-Bissau"
levels(country.means.map$names)[105]       <- "North Macedonia"
levels(country.means.map$names)[117]       <- "Republic of Moldova"
levels(country.means.map$names)[116]       <- "Micronesia, Federated States of"
#levels(country.means.map$names)[31]        <- "Cape Verde"

country.means.map[,1]                      <- c(30,3,1,6,132,5,0,7,4,8,56,2,12,9,13,10,101,128,14,18,11,31,16,15,
                                                180,20,22,21,208,28,41,24,34,23,39,35,32,29,129,37,36,26,42,38,
                                                86,79,40,43,57,44,45,46,47,27,48,49,54,51,53,52,55,60,59,64,63,66,
                                                65,67,71,68,73,70,69,74,75,166,76,77,78,80,81,82,223,83,87,50,84,85,89,88,90,97,91,
                                                94,96,92,98,100,99,179,103,106,102,135,111,107,122,120,118,112,116,230,115,114,119,
                                                62,170,109,139,113,121,17,216,155,154,152,158,157,125,151,124,93,153,138,117,161,
                                                229,163,165,159,160,171,162,164,172,167,169,173,174,176,187,211,218,198,175,181,236,
                                                177,183,184,104,182,19,185,178,95,186,25,188,156,219,189,191,190,244,194,203,
                                                193,227,197,195,196,192,199,200,202,201,204,206,226,205,207,209,210,150,212,214,220,221,222)

dsn                                       <- "C:\\Users\\mariz\\Dropbox\\PhD_Project\\papers_code_data\\2nd year\\Postprocessing"
country.means.map$id                      <- as.character(country.means.map$id)
world.map                                 <- readOGR(dsn=dsn,layer="World_map")
world.ggmap                               <- fortify(model = world.map,Regions="NAME")
#world.ggmap$id[world.ggmap$nam=="Svalbard"]<- "154"
world.ggmapsr                             <- merge(world.ggmap,country.means.map,by="id",all=TRUE)
world.ggmapsr                             <- arrange(world.ggmapsr,order)

#world.ggmapsr$mean_D[world.ggmapsr$id=="179"& world.ggmapsr$hole==TRUE] <- country.means.map$mean_D[country.means.map$names=="Lesotho"]
colors                                    <- brewer.pal(9,"YlOrRd")#
colors2                                   <- c('yellow','green','red', 'purple','cyan', 'pink', 'orange', 'blue','brown')

world.plot                                <- ggplot(data=world.ggmapsr, aes(x=long, y=lat, group=group),size=1e-100) + geom_path(size=1e-100) +
  scale_y_continuous(breaks=(-2:2) * 30) + scale_x_continuous(breaks=(-4:4) * 45) + 
  theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# global.plotD                              <- world.plot + geom_polygon(aes(fill = mean_D), size=1e-100,colour="black") + annotate("text",x=-183,y=90,label="Women",size=4.25) + 
#                                               scale_fill_gradientn(colours=colors,name=expression(atop("Age-standardised ","Diastolic BP(kg/"*m^2*")")),na.value="darkgrey",lim=c(min(country.means.map[,3]),max(country.means.map[,3])),
#                                               breaks=c(50,65,75,85,95),labels=c(50,65,75,85,95)) + theme(legend.title=element_text(hjust=0))
# global.plotS                              <- world.plot + geom_polygon(aes(fill = mean_S), size=1e-100,colour="black") + annotate("text",x=-183,y=90,label="Women",size=4.25) + 
#                                               scale_fill_gradientn(colours=colors,name=expression(atop("Age-standardised ","Systolic BP(kg/"*m^2*")")),na.value="darkgrey",lim=c(min(country.means.map[,4]),max(country.means.map[,4])),
#                                               breaks=c(80,100,120,140,155),labels=c(80,100,120,140,155)) + theme(legend.title=element_text(hjust=0))
# global.plotINT                            <- world.plot + geom_polygon(aes(fill = mean_S), size=1e-100,colour="black") + annotate("text",x=-183,y=90,label="Women",size=4.25) + 
#                                               scale_fill_gradientn(colours=colors2,name=expression(atop("Age-standardised ","Interaction BP(kg/"*m^2*")")),na.value="darkgrey",lim=c(min(country.means.map[,5]),max(country.means.map[,5])),
#                                               breaks=c(-100,-50,0,50,100,150,200),labels=c(-100,-50,0,50,100,150,200)) + theme(legend.title=element_text(hjust=0))

caribbean                                 <- c("Antigua and Barbuda","Bahamas","Barbados","Cuba","Dominica","Dominican Republic","Grenada","Haiti","Jamaica","Saint Kitts and Nevis","Saint Lucia",
                                               "Puerto Rico","Trinidad and Tobago","Saint Vincent and the Grenadines")
world.map2                                <- world.map[world.map@data$NAME%in%caribbean,]
world.ggmap                               <- fortify(model = world.map2,Regions="NAME")
world.ggmapsr                             <- merge(world.ggmap,country.means.map,by="id",all=TRUE)
world.ggmapsr                             <- arrange(world.ggmapsr,order)
world.plot1                               <- ggplot(data=world.ggmapsr, aes(x=long, y=lat, group=group)) + geom_path(size=1e-100) + scale_y_continuous(breaks=(-2:2) * 30) + scale_x_continuous(breaks=(-4:4) * 45)  

# caribbean.plotD                         <- world.plot1 + geom_polygon(aes(fill = mean_D), size=1e-100,colour="black") + scale_fill_gradientn(colours=colors,na.value="darkgrey",
#                                               lim=c(min(country.means.map[,3]),max(country.means.map[,3]))) + ggtitle(paste("Caribbean")) + theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
#                                               axis.ticks = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.background = element_rect(fill = "white"),
#                                               panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(plot.title=element_text(hjust=0.5))
# caribbean.plotS                         <- world.plot1 + geom_polygon(aes(fill = mean_S), size=1e-100,colour="black") + scale_fill_gradientn(colours=colors,na.value="darkgrey",
#                                               lim=c(min(country.means.map[,4]),max(country.means.map[,4]))) + ggtitle(paste("Caribbean")) + theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
#                                               axis.ticks = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.background = element_rect(fill = "white"),
#                                               panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(plot.title=element_text(hjust=0.5))
# caribbean.plotINT                       <- world.plot1 + geom_polygon(aes(fill = mean_S), size=1e-100,colour="black") + scale_fill_gradientn(colours=colors,na.value="darkgrey",
#                                               lim=c(min(country.means.map[,5]),max(country.means.map[,5]))) + ggtitle(paste("Caribbean")) + theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
#                                               axis.ticks = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.background = element_rect(fill = "white"),
#                                               panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(plot.title=element_text(hjust=0.5))
country.list                              <- read.csv("C:\\Users\\mariz\\Dropbox\\PhD_Project\\papers_code_data\\Adult_BMI_model54_test\\country_list2015-09-15.csv")[,1:2]
country.list[,2]                          <- as.character(country.list[,2])
country.list[63,2]                        <- "Micronesia, Federated States of"
inset.countries                           <- c("American Samoa","Bahrain","Bermuda","Brunei Darussalam","Cabo Verde","Comoros","Cook Islands","Fiji","French Polynesia","Kiribati","Maldives",
                                               "Marshall Islands","Mauritius","Micronesia, Federated States of","Montenegro","Nauru","Niue","Palau","Samoa","Sao Tome and Principe","Seychelles",
                                               "Solomon Islands","Tokelau","Tonga","Tuvalu","Vanuatu")
inset.iso                                 <- NA
for(i in 1:length(inset.countries)) inset.iso[i] <- as.character(country.list[country.list[,2]==inset.countries[i],1])
inset.title                               <- c("American Samoa","Bahrain","Bermuda","Brunei Darussalam","Cabo Verde","Comoros","Cook Islands","Fiji","French Polynesia","Kiribati","Maldives",
                                               "Marshall Islands","Mauritius","Micronesia F.S.","Montenegro","Nauru","Niue","Palau","Samoa","Sao Tome and Principe","Seychelles","Solomon Islands",                
                                               "Tokelau","Tonga","Tuvalu","Vanuatu")   
# lim.min                                   <- as.numeric(apply(country.means.map,2,min)[3:5])
# lim.max                                   <- as.numeric(apply(country.means.map,2,max)[3:5])
lim.min<- lim.max <- NULL
for (i in 3:5) {
  lim.min[i-2]<-summary(country.means.map[,i])[1]
  lim.max[i-2]<-summary(country.means.map[,i])[6]
}
pdf(paste(filename,"_BP_map_",map.year,".pdf",sep=''),width=15,height=10)
for (m in 3:5)   {
  for(i in 1:length(inset.countries)) {
    country.name.temp                       <- inset.countries[i]
    country.iso.temp                        <- inset.iso[i]
    country.title.temp                      <- inset.title[i]
    dat.temp                                <- subset(country.means.map,country.means.map[,2]==country.name.temp)
    square                                  <- as.data.frame(cbind(as.vector(c(0,0)),as.vector(c(1,1))))
    square                                  <- cbind(square,dat.temp[,m])
    names(square)                           <- c("x","y","col")
    assign(paste(country.iso.temp,".plot",sep=''),ggplot(aes(x=x,y=y),data=square) + geom_rect(aes(fill = col),xmin=-1,xmax=1,ymin=-1,ymax=1.5) + 
             scale_fill_gradientn(colours=colors,lim=c(lim.min[m-2],lim.max[m-2])) + ggtitle(country.title.temp) + theme(
               text=element_text(size=7),plot.margin = unit(c(0.2,1,0.2,1), "lines"),axis.line = element_blank(),axis.text.x = element_blank(), 
               axis.text.y = element_blank(),axis.ticks = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border=
                 element_rect(colour = "black", fill = NA, size=0.2),legend.position = "none")+theme(plot.title=element_text(hjust=0.5)) )
  }
  if (m==3) {
    caribbean.plot     <- world.plot1 + geom_polygon(aes(fill = mean_D), size=1e-100,colour="black") +
      scale_fill_gradientn(colours=colors,na.value="darkgrey",lim=c(min(country.means.map[,3]),max(country.means.map[,3])))+
      ggtitle(paste("Caribbean")) + theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
                                          axis.ticks = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
                                          panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) + theme(plot.title=element_text(hjust=0.5))
    
    global.plot       <-  world.plot + geom_polygon(aes(fill = mean_D), size=1e-100,colour="black") +
      annotate("text",x=-183,y=90,label="Women",size=4.25) + scale_fill_gradientn(colours=colors,
                                                                                  name=expression(atop("Age-standardised ","Diastolic BP(kg/"*m^2*")")),na.value="darkgrey",
                                                                                  lim=c(min(country.means.map[,3]),max(country.means.map[,3])),breaks=c(50,65,75,85,95),
                                                                                  labels=c(50,65,75,85,95)) + theme(legend.title=element_text(hjust=0))
    
  }else if (m==4){
    caribbean.plot    <- world.plot1 + geom_polygon(aes(fill = mean_S), size=1e-100,colour="black") + 
      scale_fill_gradientn(colours=colors,na.value="darkgrey",lim=c(min(country.means.map[,4]),max(country.means.map[,4]))) +
      ggtitle(paste("Caribbean")) + theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
                                          axis.ticks = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
                                          panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) + theme(plot.title=element_text(hjust=0.5))
    
    global.plot       <- world.plot + geom_polygon(aes(fill = mean_S), size=1e-100,colour="black") +
      annotate("text",x=-183,y=90,label="Women",size=4.25) +  scale_fill_gradientn(colours=colors,
                                                                                   name=expression(atop("Age-standardised ","Systolic BP(kg/"*m^2*")")),na.value="darkgrey",
                                                                                   lim=c(min(country.means.map[,4]),max(country.means.map[,4])),breaks=c(80,100,120,140,155),
                                                                                   labels=c(80,100,120,140,155)) + theme(legend.title=element_text(hjust=0))
  }else if (m==5){
    caribbean.plot    <- world.plot1 + geom_polygon(aes(fill = mean_S), size=1e-100,colour="black") + 
      scale_fill_gradientn(colours=colors,na.value="darkgrey",lim=c(min(country.means.map[,5]),max(country.means.map[,5]))) + 
      ggtitle(paste("Caribbean")) + theme(axis.line = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),
                                          axis.ticks = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",
                                          panel.background = element_rect(fill = "white"),panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) + theme(plot.title=element_text(hjust=0.5))
    
    global.plot       <- world.plot + geom_polygon(aes(fill = mean_S), size=1e-100,colour="black") +
      annotate("text",x=-183,y=90,label="Women",size=4.25) + scale_fill_gradientn(colours=colors,
                                                                                  name=expression(atop("Age-standardised ","Interaction BP(kg/"*m^2*")")),na.value="darkgrey",
                                                                                  lim=c(min(country.means.map[,5]),max(country.means.map[,5])),breaks=c(-100,-50,0,50,100,150,200),
                                                                                  labels=c(-100,-50,0,50,100,150,200)) + theme(legend.title=element_text(hjust=0))
    
  }
  grid.arrange(arrangeGrob(global.plot,arrangeGrob(caribbean.plot,arrangeGrob(
    ASM.plot,BHR.plot,BMU.plot,BRN.plot,CPV.plot,COM.plot,COK.plot,FJI.plot,PYF.plot,
    KIR.plot,MDV.plot,MHL.plot,MUS.plot,FSM.plot,MNE.plot,NRU.plot,NIU.plot,PLW.plot,
    WSM.plot,STP.plot,SYC.plot,SLB.plot,TKL.plot,TON.plot,TUV.plot,VUT.plot,
    nrow=4,ncol=7,widths=unit.c(rep(unit(1/7, "npc"),7)),
    heights=unit.c(unit(.25, "npc"),unit(.25, "npc"),unit(.25, "npc"),unit(.25, "npc")) ),
    nrow=1,ncol=2,widths=unit.c(unit(.2, "npc"),unit(.8, "npc")) 	),
    nrow=2,ncol=1,heights=unit.c(unit(.8, "npc"),unit(.2, "npc")) ) )
}
dev.off() 
#############
