
# Note: This code checks whether a configuration has already been fitted
# and replaced then

###################
# SOURCE DATA ##### 
# BONEMARROW (bmt)#
###################

source("cure_utils.R")

#install.packages("smcure")
library("smcure")
data(bmt)
colnames(bmt) <- c("time", "delta","Auto")

#Re-scale data
bmt$time <- bmt$time / max(bmt$time)


################
# INLA EXTENSION
################
#install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
#update.packages("INLA")
#install.packages("INLABMA")
library(INLABMA)
library(INLA)
#############


#################
##  ALGORITHM  ##
#################
niter <- 200
#niter <- 40
##########
# STEP 0.#
##########

#z: matrix to save the different z vector configurations 
#explored.
#z dimension: number of individuals * number of iterations
z <- matrix(NA, nrow = nrow(bmt), ncol = niter + 1) 
#z[,1] <- sample(0:1, nrow(bmt), rep = TRUE)
#
##Not censoring observartions always  z=0 (uncured group)
z[, 1][which(bmt$delta == 1)] <- 0 

z[, 1][which(bmt$delta == 0)] <- sample(0:1, sum(bmt$delta == 0), rep = TRUE)

# List of different configurations
z.list <- list()
logistic.list <- list()
survival.list <-list()
n.z <- 0 #Number of different values of z altogether

z.idx <- rep(NA, niter) #Index to indicate the z, logictic, etc. obtained in a given iteration

##############################
## USEFUL OBJECTS FOR STEP 1.#
##############################

# LOGISTIC REGRESSION MODEL

#d.logistic: array to save each iteration database.
#Note that d.logistic databases only differ in z values 
#in each iteration step.
d.logistic <- bmt[, c("time", "delta", "Auto")]
d.logistic$z <- z[, 1] #Starting point

#logistic.inla: list to save each iteration logistic model 
#posterior distribution.
logistic.inla <- list()

#mliklogistic: matrix to save the logistic model marginal 
#log-likelihood.
mliklogistic <- matrix(NA, nrow = 2, ncol = niter)

#beta1: matrix to save the modes of the conditional 
#posterior marginals of the incidence model.
beta1 <- matrix(NA, nrow = 2, ncol = niter)

# SURVIVAL REGRESSION MODEL

#d.survival: list to save each iteration survival database.
d.survival <- list()

#survival.inla: list to save each iteration survival model 
#posterior distribution.
survival.inla <- list()

#mliklogistic: matrix to save the survival model marginal 
#log-likelihood
mliksurvival <- matrix(NA, nrow = 2, ncol = niter)

#beta2 and alpha: matrix to save each iteration survival model 
#posterior conditional models.
beta2 <- matrix(NA, nrow = 2, ncol = niter)
alpha <- matrix(NA, nrow = 1, ncol = niter)



##############################
## USEFUL OBJECTS FOR STEP 2.#
##############################

x1 <- as.matrix(cbind(rep(1, nrow(bmt)), bmt$Auto))
x2 <- x1

##############################
## USEFUL OBJECTS FOR STEP 3.#
##############################

eta <- matrix(NA, nrow = nrow(bmt), ncol = niter) 
su <- matrix(NA, nrow = nrow(bmt), ncol = niter)
pz <- matrix(NA, nrow = nrow(bmt), ncol = niter)

#Setup using index
n.z <- 0

# Just create a name for the configuration z (a better algorithm can be used)
z.id <- function(z) {
  return(paste(z, collapse = ""))
}

Sys.time()->start;
i <- 1 
while(i <= niter) {

  print(paste0("**ITERATION: ", i ))

  # CHeck wether has already been sampled
  zz <- z.id(z[, i])
  aux.idx <- which(zz == names(z.list))

  if(length(aux.idx) == 0) { #Fit models 
  
  #########
  #STEP 1.#
  #########
  
  
  #FIT LOGISTIC REGRESSION MODEL WITH INLA 
  d.logistic$z <- z[, i]
  
  logistic.inla <- inla(z ~ 1 + Auto, family = "binomial", 
    data = as.data.frame(d.logistic), Ntrials = 1,
    control.predictor = list(link = 1),
    control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
      mean = 0, prec = 0.001)) 
  
  #FIT SURVIVAL MODEL WITH INLA 
  
  d.survival <- subset(d.logistic, z == 0)
  
  survival.inla <- inla(inla.surv(time, delta) ~ 1 + Auto, 
    data = d.survival, family = "weibullsurv", 
    control.predictor = list(link = 1),
    control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
      mean = 0, prec = 0.001),
    control.mode =  list(theta = 0.1, restart = TRUE),
    #control.mode =  list(theta = ifelse(i == 1, 0, survival.inla[[i-1]]$mode$theta), restart = TRUE),
    control.family = list(prior = "loggamma", param = c(0.1, 0.1)))
  

  print(table(d.survival$delta))

    #Add configuration and results
    n.z <- n.z + 1 
    aux.idx <- n.z

    z.list[[n.z]] <- z[, i]
    z.idx[i] <- aux.idx
    names(z.list)[n.z] <- z.id(z[, i])
    logistic.list[[n.z]] <- logistic.inla
    survival.list[[n.z]] <- survival.inla

  } else { #z has already appeared
    z.idx[i] <- aux.idx
  }

  ###########  
  #STEP 2.###
  ###########
  
  beta1[, i] <- logistic.list[[aux.idx]]$summary.fixed[, "mode"]
  beta2[, i] <- survival.list[[aux.idx]]$summary.fixed[, "mode"]
  alpha[, i] <- survival.list[[aux.idx]]$summary.hyperpar[, "mode"]
  
  #marginal log-likelihoods to check convergence
  mliklogistic[, i] <- logistic.list[[aux.idx]]$mlik
  mliksurvival[, i] <- survival.list[[aux.idx]]$mlik
  
  ###########  
  #STEP 3.###
  ###########
  
  # Compute cure proportion
  eta[, i] <- cureprop(x1, beta1[, i])
  
  # Compute survival
  su[, i] <- survunc(bmt$time, alpha[, i], x2, beta2[, i])
  
  # Sample z
  
  pz[, i] <- piz(bmt$delta, eta[, i], su[, i])  
  
  z[, i + 1] <- upz(pz[, i])

  i <- i + 1
} 
print(Sys.time()-start);
#Time difference of 1.240756 mins

###Burning and thinning
nthining <- seq(21, niter, by = 2)
weights <- table(z.idx[nthining]) / length(nthining)

z.inla <- z[, as.integer(names(weights))]
survival.inla <- survival.list[as.integer(names(weights))]
logistic.inla <- logistic.list[as.integer(names(weights))]

#BMA averaging
length(nthining)
bma.logistic <- INLABMA:::fitmargBMA2(logistic.inla, weights,
  "marginals.fixed")
bma.survival <- INLABMA:::fitmargBMA2(survival.inla, weights,
  "marginals.fixed")
bma.survival.alpha <- INLABMA:::fitmargBMA2(survival.inla,
  weights, "marginals.hyperpar")

###SUMMARY POSTERIORS
inla.zmarginal(bma.logistic[[1]])#ball1
inla.zmarginal(bma.logistic[[2]])#baut1

inla.zmarginal(bma.survival[[1]])#ball2
inla.zmarginal(bma.survival[[2]])#baut2





###########
# JAGS CODE
############
library(rjags)
library(random)
library(R2WinBUGS)
library(doParallel)
#############


###################
#    MIXTURE      #
#CURE RATE MODEL-II#
#WEIBULL COVARIATES#
###################

#############

# Parameters for simulation
ni<-200000
#ni<-20000
na<-0.2*ni; nt<-max(1, floor(2*ni/1000)) #1500 sample size


bonemarrow_mixtureII <- function(){

  for(i in 1:n.z) {
    # Non-observed events
    z[idx.z[i]] ~ dbern(eta[idx.z[i]])
    #logit(eta[idx.z[i]]) <- betaAll1 + betaAut1 * Auto[idx.z[i]]
  }

  for(i in 1:n.zfixed) {
    # Observed events: non-cured
    z[idx.zfixed[i]] <- 0
    #logit(eta[idx.zfixed[i]]) <- betaAll1 + betaAut1 * Auto[idx.zfixed[i]]

  }

  for(i in 1: N)
  {
    # LOGISTIC REGRESSION MODEL 
    
    #z[i] ~ dbern(eta[i])
    
    logit(eta[i]) <- betaAll1+betaAut1*Auto[i]
    
    # WEIBULL AFT MODEL
    
    #scale parameter
    lambda[i] <- exp(betaAll2 + betaAut2 * Auto[i])#beta0
    
    #survival function
    s[i] <- exp(-lambda[i]*pow(time[i],alpha))
    
    #hazard function
    h[i] <- lambda[i]*alpha*pow(time[i],alpha-1)
    
    #density function
    f[i] <- lambda[i]*alpha*pow(time[i],alpha-1)*exp(-lambda[i]*pow(time[i],alpha))
    
    
    # DEFINITION OF THE LOG-LIKELIHOOD USING ZEROS TRICK
    #delta: not censoring indicator
    
    L[i] <- pow((1-eta[i])*f[i],delta[i])*pow((eta[i])+(1-eta[i])*s[i],1-delta[i])
    
    logL[i] <- log(L[i])
    
    phi[i] <- C - logL[i]
    
    zeros[i] ~ dpois(phi[i])
    
  }
  
  # MARGINAL PRIOR DISTRIBUTIONS
  
  betaAll1 ~ dnorm(0,0.001)
  betaAut1 ~ dnorm(0,0.001)
  betaAll2 ~ dnorm(0,0.001)
  betaAut2 ~ dnorm(0,0.001)
  sigma ~ dunif(0,10)
  alpha <- 1 / sigma
  
}

filename.JAGS <- file.path( paste("bonemarrow_mixtureII.txt"))
write.model(bonemarrow_mixtureII, filename.JAGS)




data.JAGS <- list(N = nrow(bmt), time = bmt$time, delta = bmt$delta,
  C = 500000, zeros = rep(0, nrow(bmt)), Auto = bmt$Auto, 
  z = rep(NA, length = nrow(bmt))) #age=e1684$age, gender=e1684$gender, ps=e1684$ps)

#Add fixed and nonfixed z's
data.JAGS$idx.z <- which(data.JAGS$delta == 0)
data.JAGS$idx.zfixed <- which(data.JAGS$delta == 1)
data.JAGS$n.z <- length(data.JAGS$idx.z)
data.JAGS$n.zfixed <- length(data.JAGS$idx.zfixed)

registerDoParallel(3)
seeds <- c(34,73,1500)
jags.inits <- function(i){
  return(list(.RNG.name="lecuyer::RngStream", .RNG.seed=seeds[i]))
}
mcmc.combine <- function(...){
  return(as.mcmc.list(sapply(list(...), mcmc)))
}

parameters.JAGS <- c("z","eta","betaAll1","betaAut1", "alpha","betaAll2","betaAut2")#

Sys.time()->start;
jags.parsamples <- foreach( 
  i=1:getDoParWorkers(), .inorder=F, .packages=c('rjags','random'), .combine="mcmc.combine", .multicombine=T) %dopar% {
    load.module("lecuyer")
    model.jags <- jags.model(data=data.JAGS, file=filename.JAGS, inits=jags.inits(i), n.adapt=na)
    result <- coda.samples(model.jags, variable.names=parameters.JAGS, n.iter=ni, thin=nt)
    return(result)
  } 
print(Sys.time()-start);
summary(jags.parsamples,dig=3)

#-----------MCMC Convergence----------
# plot(jags.parsamples)
gelman.diag(jags.parsamples, multivariate = FALSE)
#-------------------------------------
#---------------------------------------------------------------#
#                           RESULTS                             #
#---------------------------------------------------------------#


## Prepare MCMC outcomes - Multiplicative frailties
result1 <- do.call(rbind, jags.parsamples)

sims.list <- vector("list", length(parameters.JAGS))

names(sims.list) <- parameters.JAGS
# names(sims.list2) <- parameters2

## Put the posterior sample of each parameter in a list
for(p1 in seq_along(parameters.JAGS)){
  iik <- grep(paste("^", parameters.JAGS[p1], sep=""), colnames(result1))
  sims.list[[p1]] <- result1[,iik]
}

## PARAMETERS AND DERIVED QUANTITIES SAMPLES #1500 sample size
parameters.JAGS <- c("z", "eta", "betaAll1", "betaAut1", "alpha",
  "betaAll2", "betaAut2")#
z<- sims.list[[1]]  
eta <- sims.list[[2]] 
betaAll1 <- sims.list[[3]] 
betaAut1<- sims.list[[4]]
alpha<-sims.list[[5]]
betaAll2<-sims.list[[6]]
betaAut2<-sims.list[[7]]


#
#  SUMMARY OF RESULTS
#



###SUMMARY POSTERIORS
inla.zmarginal(bma.logistic[[1]])#ball1
inla.zmarginal(bma.logistic[[2]])#baut1

inla.zmarginal(bma.survival[[1]])#ball2
inla.zmarginal(bma.survival[[2]])#baut2

#JAGS
summary(betaAll1)
summary(betaAut1)
summary(betaAll2)
summary(betaAut2)


#PLot stuff
par(mfrow = c(2, 2))

plot(bma.logistic[[1]], type = "l", main = "ball1")
lines(density(betaAll1), col = "red")
plot(bma.logistic[[2]], type = "l", main = "baut1")
lines(density(betaAut1), col = "red")
plot(bma.survival[[1]], type = "l", main = "ball2")
lines(density(betaAll2), col = "red")
plot(bma.survival[[2]], type = "l", main = "baut2")
lines(density(betaAut2), col = "red")


# INLA DERIVED QUANTITIES---> Cure proportion and Survival curves por susceptible individuals

#Check convergence and configuration of z
mlik.s <- unlist(lapply(survival.inla, function(X){X$mlik[1, 1]}))
mlik.l <- unlist(lapply(logistic.inla, function(X){X$mlik[1, 1]}))
which.max(mlik.s+mlik.l)#41
hist(mlik.s+mlik.l)#models probbaility
plot(mlik.s+mlik.l)#check convergence


#Model adjustment with z most likely configuration

d.logistic.mlm<-cbind.data.frame(time=bmt$time, delta=bmt$delta, Auto=bmt$Auto, z=z.inla[, as.integer(which.max(mlik.s+mlik.l))])

logistic.mlm <- inla(z ~ 1 + Auto, family = "binomial", 
                      data = as.data.frame(d.logistic.mlm), Ntrials = 1,
                      control.compute=list(config = TRUE),
                      control.predictor = list(link = 1),
                      control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
                                           mean = 0, prec = 0.001)) 

logistic.mlm.sample<-inla.posterior.sample(3000,logistic.mlm)
str(logistic.mlm.sample[[1]])
length(logistic.mlm.sample)
coefs <- do.call(cbind, lapply(logistic.mlm.sample, function(X){X$latent[92:93, 1]}))

#Cured proportion
#i1=allogeneic
#i2=autologous
xx1 <- as.vector(t(coefs) %*% matrix(c(1, 0), ncol = 1))
xx2<-as.vector(t(coefs) %*% matrix(c(1, 1), ncol = 1))


pi1in<-1/(1+exp(-xx1)) # allogenic cure proportion
pi2in<-1/(1+exp(-xx2)) # autologous cure proportion




d.survival.mlm<-subset(d.logistic.mlm[,],d.logistic.mlm[,4]==0)

survival.mlm <- inla(inla.surv(time, delta) ~ 1 + Auto, 
                      data = d.survival.mlm, family = "weibullsurv", 
                      control.compute=list(config=TRUE),
                      control.predictor = list(link = 1),
                      control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
                                           mean = 0, prec = 0.001),
                      control.mode =  list(theta = 0.1, restart = TRUE),
                      control.family = list(prior = "loggamma", param = c(0.1, 0.1)))

survival.mlm.sample<-inla.posterior.sample(3000,survival.mlm)
str(survival.mlm.sample[[1]])
length(survival.mlm.sample)
coefs <- do.call(cbind, lapply(survival.mlm.sample, function(X){X$latent[70:71, 1]}))
alphain <-lapply(survival.mlm.sample, function(X){X$hyperpar[1]})
xx1 <- as.vector(t(coefs) %*% matrix(c(1, 0), ncol = 1))
xx2<-as.vector(t(coefs) %*% matrix(c(1, 1), ncol = 1))

lpi1in<-exp(xx1)
lpi2in<-exp(xx2)



t<-seq(0,max(bmt$time),by=0.01)
spi1in=matrix(NA,nrow=length(xx1),ncol=length(t))
spi2in=matrix(NA,nrow=length(xx1),ncol=length(t))

for(i in 1:length(xx1)){
  for(j in 1:length(t)){
    
    spi1in[i,j]=exp(-lpi1in[i]*(t[j]^alphain[i][[1]]))
    spi2in[i,j]=exp(-lpi2in[i]*(t[j]^alphain[i][[1]]))
    
  }
}

# Mean of the posterior distribution

mspi1in <- apply(spi1in,2,mean)# allogeneic
mspi2in <- apply(spi2in,2,mean)# autologous


# MCMC DERIVED QUANTITIES---> Cure proportion and Survival curves por susceptible individuals

####POSTERIOR SUMMARY OF DERIVED QUANTITIES MCMC
#Cured proportion
#i1=allogeneic
#i2=autologous

pi1=exp(betaAll1)/(1+exp(betaAll1))
pi2=exp(betaAll1+betaAut1)/(1+exp(betaAll1+betaAut1))

#Survival uncured proportion

lpi1=vector()
lpi2=vector()

for(i in 1:length(betaAll1)){
  
  lpi1[i]<-exp(betaAll2[i])
  lpi2[i]=exp(betaAll2[i]+betaAut2[i])
  
}

t<-seq(0,max(bmt$time),by=0.01)
spi1=matrix(NA,nrow=length(betaAll1),ncol=length(t))
spi2=matrix(NA,nrow=length(betaAll1),ncol=length(t))


for(i in 1:length(betaAll1)){
  for(j in 1:length(t)){
    
    spi1[i,j]=exp(-lpi1[i]*(t[j]^alpha[i]))
    spi2[i,j]=exp(-lpi2[i]*(t[j]^alpha[i]))
    
    
  }
}

mspi1<-apply(spi1,2,mean)
mspi2<-apply(spi2,2,mean)



par(mar=c(8,4,4,4))
plot(t,mspi1in,type="l",lwd=10,ylab="",xlab="",cex.axis=3,cex.lab=4)
lines(t,mspi1,col="red",lty=2,lwd=10)
mtext(text="Time", side=1, line=3.7, las=1,cex=3,at=830)
grid(nx = NULL, ny = NULL, col = "gray50", lty = "solid",
     lwd = par("lwd"), equilogs = TRUE)
box(col = 'gray50')




par(mar=c(8,4,4,4))
plot(t,mspi2in,type="l",lwd=10,ylab="",xlab="",cex.axis=3,cex.lab=4)
lines(t,mspi2,col="red",lty=2,lwd=10)
mtext(text="Time (days)", side=1, line=3.7, las=1,cex=3,at=830)
grid(nx = NULL, ny = NULL, col = "gray50", lty = "solid",
     lwd = par("lwd"), equilogs = TRUE)
box(col = 'gray50')








#Estimates of z
dev.new()

plot((z.inla %*% weights)[bmt$delta == 0],
  apply(z, 2, mean)[bmt$delta == 0],
  ylim = c(0,1), xlim = c(0,1), xlab = "INLA", ylab = "MCMC")
abline(0, 1)

# Fit model with two extreme cases:
# (1) all patients in survival model
# (2) only patients with observed event in survival model

m1 <- inla(inla.surv(time, delta) ~ 1 + Auto,
  data = bmt, family = "weibullsurv",
  control.predictor = list(link = 1),
  control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
    mean = 0, prec = 0.001),
  control.mode =  list(theta = 0.1, restart = TRUE),
  control.family = list(prior = "loggamma", param = c(0.1, 0.1)))

m2 <- inla(inla.surv(time, delta) ~ 1 + Auto,
  data = bmt[bmt$delta == 1, ], family = "weibullsurv",
  control.predictor = list(link = 1),
  control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
    mean = 0, prec = 0.001),
  control.mode =  list(theta = 0.1, restart = TRUE),
  control.family = list(prior = "loggamma", param = c(0.1, 0.1)))

#PLot stuff
pdf(file = "bonemarrow-margs.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))

plot(bma.survival[[1]], type = "l", main = "ball2")
lines(density(betaAll2), col = "red")
lines(m1$marginals.fixed[[1]], col = "blue")
lines(m2$marginals.fixed[[1]], col = "green")
plot(bma.survival[[2]], type = "l", main = "baut2")
lines(density(betaAut2), col = "red")
lines(m1$marginals.fixed[[2]], col = "blue")
lines(m2$marginals.fixed[[2]], col = "green")
dev.off()



#save(file = "code_bonemarrow-speedup.RData", list = ls())

