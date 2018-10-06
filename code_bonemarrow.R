
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


##############################
## USEFUL OBJECTS FOR STEP 1.#
##############################

# LOGISTIC REGRESSION MODEL

#d.logistic: array to save each iteration database.
#Note that d.logistic databases only differ in z values 
#in each iteration step.
d.logistic <- array(NA, dim = c(nrow(bmt), ncol(bmt) + 1, niter + 1),
  dimnames = list(NULL, c("time", "delta", "Auto", "z"), NULL))

d.logistic[ , , 1] <- matrix(c(bmt$time, bmt$delta, bmt$Auto,
  z[, 1]), byrow = FALSE)

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


i <- 1 
while(i <= niter) {

  print(paste0("**ITERATION: ", i ))
  
  #########
  #STEP 1.#
  #########
  
  
  #FIT LOGISTIC REGRESSION MODEL WITH INLA 
  
  logistic.inla[[i]] <- inla(z ~ 1 + Auto, family = "binomial", 
    data = as.data.frame(d.logistic[ , , i]), Ntrials = 1,
    control.predictor = list(link = 1),
    control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
      mean = 0, prec = 0.001)) 
  
  #FIT SURVIVAL MODEL WITH INLA 
  
  d.survival[[i]] <- array(NA, dim = c(length(which(d.logistic[, 4, i] == 0)),
    ncol(bmt) + 1, niter),
    dimnames = list(NULL, c("time", "delta", "Auto", "z"), NULL))
  
  
  d.survival[[i]] <- as.data.frame(subset(d.logistic[ , , i],
    d.logistic[, 4, i] == 0))
  
  survival.inla[[i]] <- inla(inla.surv(time, delta) ~ 1 + Auto, 
    data = d.survival[[i]], family = "weibullsurv", 
    control.predictor = list(link = 1),
    control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
      mean = 0, prec = 0.001),
    control.mode =  list(theta = 0.1, restart = TRUE),
    #control.mode =  list(theta = ifelse(i == 1, 0, survival.inla[[i-1]]$mode$theta), restart = TRUE),
    control.family = list(prior = "loggamma", param = c(0.1, 0.1)))
  

  print(table(d.survival[[i]]$delta))

  ###########  
  #STEP 2.###
  ###########
  
  beta1[, i] <- logistic.inla[[i]]$summary.fixed[, "mode"]
  beta2[, i] <- survival.inla[[i]]$summary.fixed[, "mode"]
  alpha[, i] <- survival.inla[[i]]$summary.hyperpar[, "mode"]
  
  #marginal log-likelihoods to check convergence
  mliklogistic[, i] <- logistic.inla[[i]]$mlik
  mliksurvival[, i] <- survival.inla[[i]]$mlik
  
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
  d.logistic[ , , i + 1] <- matrix(c(bmt$time, bmt$delta, bmt$Auto,
    z[, i + 1]), byrow = TRUE)

  i <- i + 1
} 

###Burning and thinning
nthining <- seq(21, niter, by = 2)
survival.inla <- survival.inla[c(nthining)]
logistic.inla <- logistic.inla[c(nthining)]
z.inla <-z[, c(nthining)]

#BMA averaging
length(nthining)
bma.logistic <- INLABMA:::fitmargBMA2(logistic.inla, rep(1 / length(nthining),
  length(nthining)), "marginals.fixed")
bma.survival <- INLABMA:::fitmargBMA2(survival.inla, rep(1 / length(nthining),
  length(nthining)), "marginals.fixed")
bma.survival.alpha <- INLABMA:::fitmargBMA2(survival.inla,
  rep(1 / length(nthining), length(nthining)), "marginals.hyperpar")

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
plot(jags.parsamples)
gelman.diag(jags.parsamples,multivariate=F)
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

#Estimates of z
dev.new()

plot(apply(z.inla, 1, mean)[bmt$delta == 0], apply(z, 2, mean)[bmt$delta == 0],
  ylim = c(0,1), xlim = c(0,1), xlab = "INLA", ylab = "MCMC")
abline(0, 1)


