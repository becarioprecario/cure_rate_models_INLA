
###################
# SOURCE DATA ##### 
# MELANOMA (e1684)#
###################
#setwd("C:/Users/Elena/Dropbox/TESIS/models_ch5")
#setwd("C:/Users/MICOLOGIA_JOAQUIN/Dropbox/Dropbox/cure_rate_models/paper_cure_models/code")
library("smcure")
data(e1684)
e1684 <- e1684[-c(37), ] # SEX=NA; AGE=NA
colnames(e1684) <- c("IFN", "time", "delta", "A", "W")

#Re-scale time
e1684$time <- e1684$time / max(e1684$time)



################
# INLA EXTENSION
################
#install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
#update.packages("INLA")
#install.packages("INLABMA")
library(INLABMA)
library(INLA)
#############

source("cure_utils.R")

#################
##  ALGORITHM  ##
#################
niter <- 500
#niter <- 40

##########
# STEP 0.#
##########

#z: matrix to save the different z vector configurations 
#explored.
#z dimension: number of individuals * number of iterations
z <- matrix(NA, nrow = nrow(e1684), ncol = niter + 1) 
#z[, 1] <- sample(0:1, nrow(e1684), rep = TRUE)

#Not censoring observartions always  z=0 (uncured group)
z[, 1][which(e1684$delta == 1)] <- 0 

z[, 1][which(e1684$delta == 0)] <- sample(0:1, sum(e1684$delta == 0),
  rep = TRUE)

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
d.logistic <- e1684[, c("IFN", "time", "delta", "A", "W")]
d.logistic$z <- z[, 1] #Starting point

#logistic.inla: list to save each iteration logistic model 
#posterior distribution.
logistic.inla <- list()

#mliklogistic: matrix to save the logistic model marginal 
#log-likelihood.
mliklogistic <- matrix(NA, nrow = 2, ncol = niter)

#beta1: matrix to save the modes of the conditional 
#posterior marginals of the incidence model.
beta1 <- matrix(NA, nrow = 4, ncol = niter)

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
beta2 <- matrix(NA, nrow = 4, ncol = niter)
alpha <- matrix(NA, nrow = 1, ncol = niter)



##############################
## USEFUL OBJECTS FOR STEP 2.#
##############################

x1 <- as.matrix(cbind(rep(1, nrow(e1684)), e1684$W, e1684$IFN, e1684$A))
x2 <- x1

##############################
## USEFUL OBJECTS FOR STEP 3.#
##############################

eta <- matrix(NA, nrow = nrow(e1684), ncol = niter) 
su <- matrix(NA, nrow = nrow(e1684), ncol = niter)
pz <- matrix(NA, nrow = nrow(e1684), ncol = niter)

#Setup using index
n.z <- 0

# Just create a name for the configuration z (a better algorithm can be used)
z.id <- function(z) {
  return(paste(z, collapse = ""))
}

i <- 1
Sys.time() -> start;
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

  
 logistic.inla = inla(z ~ 1 + W + IFN + A, family = "binomial", 
    data = as.data.frame(d.logistic), Ntrials = 1,
    control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
      mean = 0, prec = 0.001)) 
  
  #FIT SURVIVAL MODEL WITH INLA 
  d.survival <- subset(d.logistic, z == 0)

  survival.inla <- inla(inla.surv(time, delta) ~ 1 + W + IFN + A, 
    data = d.survival, family = "weibullsurv", 
    control.predictor = list(link = 1), control.fixed = list(mean.intercept = 0,
      prec.intercept = 0.001, mean = 0, prec = 0.001),
    control.mode = list(theta = 0.15, restart = TRUE),
    #control.mode =  list(theta = ifelse(i == 1, 0, survival.inla[[i-1]]$mode$theta), restart = TRUE),

    control.family = list(prior = "loggamma", param = c(0.01, 0.01)),
    verbose = FALSE)

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
  
  beta1[,i] <- logistic.list[[aux.idx]]$summary.fixed[, "mode"]
  beta2[,i] <- survival.list[[aux.idx]]$summary.fixed[, "mode"]
  alpha[,i] <- survival.list[[aux.idx]]$summary.hyperpar[, "mode"]
  
  #marginal log-likelihoods to check convergence
  mliklogistic[,i] <- logistic.list[[aux.idx]]$mlik
  mliksurvival[,i] <- survival.list[[aux.idx]]$mlik
  
  ###########  
  #STEP 3.###
  ###########
  
  # Compute cure proportion
  eta[, i] <- cureprop(x1, beta1[, i])
  
  # Compute survival
  su[, i] <- survunc(e1684$time, alpha[, i], x2, beta2[, i])
  
  # Sample z
  
  pz[, i] <- piz(e1684$delta, eta[, i], su[, i])  
  
  z[, i + 1] <- upz(pz[, i])

  i <- i + 1
} 
print(Sys.time() - start);


###Burning and thinning
nthining <- seq(51, niter, by = 5)
weights <- table(z.idx[nthining]) / length(nthining)

z.inla <- z[, as.integer(names(weights))]
survival.inla <- survival.list[ as.integer(names(weights))]
logistic.inla <- logistic.list[ as.integer(names(weights))]

#BMA averaging
length(nthining)

#BMA averaging
bma.logistic <- INLABMA:::fitmargBMA2(logistic.inla, weights,
  "marginals.fixed")
bma.survival <- INLABMA:::fitmargBMA2(survival.inla, weights,
  "marginals.fixed")
bma.survival.alpha <- INLABMA:::fitmargBMA2(survival.inla,
  weights, "marginals.hyperpar")

#save(file = "melanoma_mgs_outcomes.RData", list = ls())

###SUMMARY POSTERIORS
inla.zmarginal(bma.logistic[[1]])#beta01
inla.zmarginal(bma.logistic[[2]])#betaW1
inla.zmarginal(bma.logistic[[3]])#betaIFN1
inla.zmarginal(bma.logistic[[4]])#betaA1
lapply(bma.logistic, function(X){
  1 - inla.pmarginal(0, X)
})


inla.zmarginal(bma.survival[[2]])#betaW2
inla.zmarginal(bma.survival[[3]])#betaIFN2
inla.zmarginal(bma.survival[[4]])#betaA2
lambda <- inla.tmarginal(exp,bma.survival[[1]])
inla.zmarginal(lambda)#lambda
inla.zmarginal(bma.survival.alpha[[1]])#alpha


############
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
ni <- 100000
#ni <- 10000
na <- 0.2 * ni
nt <- max(1, floor(2 * ni / 1000)) #1500 sample size

melanoma_mixtureII <- function(){
  
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
    
    logit(eta[i]) <- beta01 + betaW1 * W[i] + betaIFN1 * IFN[i] + betaA1 * A[i]
    
    # COX PROPORTIONAL HAZARDS MODEL
    #lambba[i]
    lambda[i] <- exp(beta02 + betaW2 * W[i] + betaIFN2 * IFN[i] + betaA2 * A[i])
    
    #survival function
    s[i] <- exp(-lambda[i] * pow(time[i], alpha))
    
    #hazard function
    h[i] <- lambda[i] * alpha * pow(time[i], alpha - 1)
    
    #density function
    f[i]<-lambda[i] * alpha * pow(time[i], alpha - 1) * 
      exp(-lambda[i] * pow(time[i], alpha)) 
    
    
    # DEFINITION OF THE LOG-LIKELIHOOD USING ZEROS TRICK
    #delta: not censoring indicator
    
    L[i] <- pow((1 - eta[i]) * f[i], delta[i]) *
      pow((eta[i]) + (1 - eta[i]) * s[i],1-delta[i])
    
    logL[i] <- log(L[i])
    
    phi[i] <- C - logL[i]
    
    zeros[i] ~ dpois(phi[i])
    
  }
  
  # MARGINAL PRIOR DISTRIBUTIONS
  
  beta01 ~ dnorm(0, 0.001)
  betaW1 ~ dnorm(0, 0.001)
  betaIFN1 ~ dnorm(0, 0.001)
  betaA1 ~ dnorm(0, 0.001)
  beta02 ~ dnorm(0, 0.001)
  betaW2 ~ dnorm(0, 0.001)
  betaIFN2 ~ dnorm(0, 0.001)
  betaA2 ~ dnorm(0, 0.001)
  alpha ~ dgamma(0.01, 0.01)
 
}


filename.JAGS <- file.path("melanoma_mixtureII.txt")
write.model(melanoma_mixtureII, filename.JAGS)




data.JAGS <- list(N = 284, time = e1684$time, delta = e1684$delta,
  C = 500000, zeros = rep(0, 284), W = e1684$W, IFN = e1684$IFN,
  A = e1684$A, z = rep(NA, length = 284)) 

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

parameters.JAGS <- c("z","eta","beta01","betaW1","betaIFN1","betaA1", 
 "alpha", "beta02", "betaW2", "betaIFN2", "betaA2")

Sys.time() -> start;
jags.parsamples <- foreach( 
  i=1:getDoParWorkers(), .inorder=F, .packages=c('rjags','random'),
    .combine="mcmc.combine", .multicombine=T) %dopar% {
      load.module("lecuyer")
      model.jags <- jags.model(data = data.JAGS, file = filename.JAGS,
        inits = jags.inits(i), n.adapt = na)
    result <- coda.samples(model.jags, variable.names = parameters.JAGS,
       n.iter = ni, thin = nt)
    return(result)
  } 
print(Sys.time() - start);
summary(jags.parsamples, dig = 3)

#-----------MCMC Convergence----------
#plot(jags.parsamples)
gelman.diag(jags.parsamples, multivariate = FALSE)
#-------------------------------------
#---------------------------------------------------------------#
#                           RESULTS                             #
#---------------------------------------------------------------#


## Prepare MCMC outcomes - Multiplicative frailties
result1 <- do.call(rbind, jags.parsamples)

sims.list <- vector("list", length(parameters.JAGS))

names(sims.list) <- parameters.JAGS

## Put the posterior sample of each parameter in a list
for(p1 in seq_along(parameters.JAGS)){
  iik <- grep(paste("^", parameters.JAGS[p1], sep = ""), colnames(result1))
  sims.list[[p1]] <- result1[, iik]
}

## PARAMETERS AND DERIVED QUANTITIES SAMPLES #1500 sample size
parameters.JAGS <- c("z", "eta", "beta01", "betaW1", "betaIFN1", "betaA1",
  "alpha", "beta02", "betaW2", "betaIFN2", "betaA2")
#z<- sims.list[[1]]  
eta <- sims.list[[2]] 
beta01 <- sims.list[[3]] 
betaW1 <- sims.list[[4]]
betaIFN1 <- sims.list[[5]]
betaA1 <- sims.list[[6]]
alpha <- sims.list[[7]]
beta02 <- sims.list[[8]]
betaW2 <- sims.list[[9]]
betaIFN2 <- sims.list[[10]]
betaA2 <- sims.list[[11]]

#save(file = "melanoma_mcmc_outcomes.RData", list = ls())


#
#  SUMMARY OF RESULTS
#



###SUMMARY POSTERIORS
inla.zmarginal(bma.logistic[[1]])
inla.zmarginal(bma.logistic[[2]])
inla.zmarginal(bma.logistic[[3]])
inla.zmarginal(bma.logistic[[4]])

inla.zmarginal(bma.survival[[1]])
inla.zmarginal(bma.survival[[2]])
inla.zmarginal(bma.survival[[3]])
inla.zmarginal(bma.survival[[4]])

#JAGS
summary(beta01)
summary(betaW1)
summary(betaIFN1)
summary(betaA1)
#summary(alpha)
summary(beta02)
summary(betaW2)
summary(betaIFN2)
summary(betaA2)



#PLot stuff
par(mfrow = c(2, 4))

plot(bma.logistic[[1]], type = "l", main = "beta01")
lines(density(beta01), col = "red")
plot(bma.logistic[[2]], type = "l", main = "betaW1")
lines(density(betaW1), col = "red")
plot(bma.logistic[[3]], type = "l", main = "betaIFN1")
lines(density(betaIFN1), col = "red")
plot(bma.logistic[[4]], type = "l", main = "betaA1")
lines(density(betaA1), col = "red")



plot(bma.survival[[1]], type = "l", main = "beta02")
lines(density(beta02), col = "red")
plot(bma.survival[[2]], type = "l", main = "betaW2")
lines(density(betaW2), col = "red")
plot(bma.survival[[3]], type = "l", main = "betaIFN2")
lines(density(betaIFN2), col = "red")
plot(bma.survival[[4]], type = "l", main = "betaA2")
lines(density(betaA2), col = "red")



# INLA DERIVED QUANTITIES---> Cure proportion and Survival curves por susceptible individuals

#Check convergence and configuration of z
mlik.s <- unlist(lapply(survival.inla, function(X) {
  X$mlik[1, 1]
}))
mlik.l <- unlist(lapply(logistic.inla, function(X) {
  X$mlik[1, 1]
}))
which.max(mlik.s + mlik.l)#
hist(mlik.s + mlik.l)#models probbaility
plot(mlik.s + mlik.l)#check convergence


#Model adjustment with z most likely configuration

d.logistic.mlm <- cbind.data.frame(time = e1684$time, delta = e1684$delta,
  W = e1684$W, IFN = e1684$IFN, A = e1684$A,
  z = z.inla[, as.integer(which.max(mlik.s + mlik.l))])


logistic.mlm = inla(z ~ 1 + W + IFN + A, family = "binomial", 
  data = as.data.frame(d.logistic.mlm), Ntrials = 1,
  control.compute = list(config = TRUE),
  control.fixed = list(mean.intercept = 0, prec.intercept = 0.001,
    mean = 0, prec = 0.001)) 


logistic.mlm.sample <- inla.posterior.sample(3000, logistic.mlm)
str(logistic.mlm.sample[[1]])
length(logistic.mlm.sample)
coefs <- do.call(cbind, lapply(logistic.mlm.sample, function(X){
  X$latent[285:288, 1]
}))

#Cured proportion
#i1=male,no treated,standarized age
#i2=male, treated, standarized age
#i3=female,no treated, standarized age
#i4,female, treated, standarized age

xx1 <- as.vector(t(coefs) %*% matrix(c(1, 0, 0, mean(e1684$A)), ncol = 1))
xx2 <- as.vector(t(coefs) %*% matrix(c(1, 0, 1, mean(e1684$A)), ncol = 1))
xx3 <- as.vector(t(coefs) %*% matrix(c(1, 1, 0, mean(e1684$A)), ncol = 1))
xx4 <- as.vector(t(coefs) %*% matrix(c(1, 1, 1, mean(e1684$A)), ncol = 1))

pi1in <- 1 / (1 + exp(-xx1))
pi2in <- 1 / (1 + exp(-xx2))
pi3in <- 1 / (1 + exp(-xx3))
pi4in <- 1 / (1 + exp(-xx4))




d.survival.mlm <- subset(d.logistic.mlm[, ], d.logistic.mlm[, 6] == 0)



survival.mlm <- inla(inla.surv(time, delta) ~ 1 + W + IFN + A, 
  data = d.survival.mlm, family = "weibullsurv", 
  control.compute = list(config = TRUE),
  control.predictor = list(link = 1), control.fixed = list(mean.intercept = 0,
    prec.intercept = 0.001, mean = 0, prec = 0.001),
  control.mode = list(theta = 0.15, restart = TRUE),
  control.family = list(prior = "loggamma", param = c(0.01, 0.01)),
  verbose = FALSE)

survival.mlm.sample <- inla.posterior.sample(3000, survival.mlm)
str(survival.mlm.sample[[1]])
length(survival.mlm.sample)
coefs <- do.call(cbind, lapply(survival.mlm.sample, function(X) {
  X$latent[197:200, 1]
}))
alphain <-lapply(survival.mlm.sample, function(X){
  X$hyperpar[1]
})
xx1 <- as.vector(t(coefs) %*% matrix(c(1, 0, 0, mean(e1684$A)), ncol = 1))
xx2 <- as.vector(t(coefs) %*% matrix(c(1, 0, 1, mean(e1684$A)), ncol = 1))
xx3 <- as.vector(t(coefs) %*% matrix(c(1, 1, 0, mean(e1684$A)), ncol = 1))
xx4 <- as.vector(t(coefs) %*% matrix(c(1, 1, 1, mean(e1684$A)), ncol = 1))

lpi1in <- exp(xx1)
lpi2in <- exp(xx2)
lpi3in <- exp(xx3)
lpi4in <- exp(xx4)


t <-  seq(0, max(e1684$time), by = 0.01)
spi1in = matrix(NA, nrow = length(xx1), ncol = length(t))
spi2in = matrix(NA, nrow = length(xx1), ncol = length(t))
spi3in = matrix(NA, nrow = length(xx1), ncol = length(t))
spi4in = matrix(NA, nrow = length(xx1), ncol = length(t))

for(i in 1:length(xx1)) {
  for(j in 1:length(t)) {
    
    spi1in[i,j] = exp(-lpi1in[i] * (t[j]^alphain[i][[1]]))
    spi2in[i,j] = exp(-lpi2in[i] * (t[j]^alphain[i][[1]]))
    spi3in[i,j] = exp(-lpi3in[i] * (t[j]^alphain[i][[1]]))
    spi4in[i,j] = exp(-lpi4in[i] * (t[j]^alphain[i][[1]]))
    
  }
}

mspi1in <- apply(spi1in,2,mean)
mspi2in <- apply(spi2in,2,mean)
mspi3in <- apply(spi3in,2,mean)
mspi4in <- apply(spi4in,2,mean)



# MCMC DERIVED QUANTITIES---> Cure proportion and Survival curves por susceptible individuals

####POSTERIOR SUMMARY OF DERIVED QUANTITIES MCMC
#Cured proportion
#i1=male,no treated,standarized age
#i2=male, treated, standarized age
#i3=female,no treated, standarized age
#i4,female, treated, standarized age

summary(beta01)
summary(betaW1)
summary(betaIFN1)
summary(betaA1)
#summary(alpha)
summary(beta02)
summary(betaW2)
summary(betaIFN2)
summary(betaA2)

pi1 = exp(beta01 + betaA1 * mean(e1684$A)) /
  (1 + exp(beta01 + betaA1 * mean(e1684$A)))
pi2 = exp(beta01 + betaIFN1 + betaA1 * mean(e1684$A)) / 
  (1 + exp(beta01 + betaIFN1 + betaA1 * mean(e1684$A)))
pi3 = exp(beta01 + betaW1 + betaA1 * mean(e1684$A)) /
  (1 + exp(beta01 + betaW1 + betaA1 * mean(e1684$A)))
pi4 = exp(beta01 + betaW1 + betaIFN1 + betaA1 * mean(e1684$A)) / 
  (1 + exp(beta01 + betaW1 + betaIFN1 + betaA1 * mean(e1684$A)))



#Survival uncured proportion

lpi1 = vector()
lpi2 = vector()
lpi3 = vector()
lpi4 = vector()

for(i in 1:length(beta02)) {
  
  lpi1[i] <- exp(beta02[i] + betaA2[i] * mean(e1684$A))
  lpi2[i] <- exp(beta02[i] + betaIFN2[i] + betaA2[i] * mean(e1684$A))
  lpi3[i] <- exp(beta02[i] + betaW2[i] + betaA2[i] * mean(e1684$A))
  lpi4[i] <- exp(beta02[i] + betaW2[i] + betaIFN2[i] + betaA2[i] * mean(e1684$A))
}

t <- seq(0, max(e1684$time), by = 0.01)
spi1 = matrix(NA, nrow = length(beta02), ncol = length(t))
spi2 = matrix(NA, nrow = length(beta02), ncol = length(t))
spi3 = matrix(NA, nrow = length(beta02), ncol = length(t))
spi4 = matrix(NA, nrow = length(beta02), ncol = length(t))

for(i in 1:length(beta02)){
  for(j in 1:length(t)){
    
    spi1[i,j] = exp(-lpi1[i] * (t[j]^alpha[i]))
    spi2[i,j] = exp(-lpi2[i] * (t[j]^alpha[i]))
    spi3[i,j] = exp(-lpi3[i] * (t[j]^alpha[i]))
    spi4[i,j] = exp(-lpi4[i] * (t[j]^alpha[i]))
    
  }
}

mspi1 <- apply(spi1, 2, mean)
mspi2 <- apply(spi2, 2, mean)
mspi3 <- apply(spi3, 2, mean)
mspi4 <- apply(spi4, 2, mean)


par(mfrow=c(2,2))

postscript(file = "./survivalmn.eps", horizontal = FALSE, width = 400,
  height = 400)
par(mar = c(8, 4, 4, 4))
plot(t, mspi1in, type = "l", lwd = 10,  ylab = "", xlab = "", cex.axis = 3,
  cex.lab = 4)
lines(t, mspi1, col = "red", lty = 2, lwd = 10)
mtext(text = "Time", side = 1, line = 3.7, las = 1, cex = 3, at = 0.5)
grid(nx = NULL, ny = NULL, col = "gray50", lty = "solid",
  lwd = par("lwd"), equilogs = TRUE)
box(col = 'gray50')
dev.off()



postscript(file = "./survivalmt.eps", horizontal = FALSE, width = 400,
  height = 400)
par(mar = c(8, 4, 4, 4))
plot(t, mspi2in, type = "l", lwd = 10, ylab = "", xlab = "", cex.axis = 3,
  cex.lab = 4)
lines(t, mspi1, col = "red", lty = 2, lwd = 10)
mtext(text = "Time", side = 1, line = 3.7, las = 1, cex = 3, at = 0.5)
grid(nx = NULL, ny = NULL, col = "gray50", lty = "solid",
  lwd = par("lwd"), equilogs = TRUE)
box(col = 'gray50')
dev.off()


postscript(file = "./survivalfn.eps", horizontal = FALSE, width = 400,
  height = 400)
par(mar = c(8, 4, 4, 4))
plot(t, mspi3in, type = "l", lwd = 10, ylab = "", xlab = "", cex.axis = 3,
  cex.lab = 4)
lines(t, mspi3, col = "red", lty = 2, lwd = 10)
mtext(text = "Time", side = 1, line = 3.7, las = 1, cex = 3, at = 0.5)
grid(nx = NULL, ny = NULL, col = "gray50", lty = "solid",
  lwd = par("lwd"), equilogs = TRUE)
box(col = 'gray50')
dev.off()


postscript(file = "./survivalft.eps", horizontal = FALSE, width = 400,
  height = 400)
par(mar = c(8, 4, 4, 4))
plot(t, mspi4in, type = "l", lwd = 10, ylab = "", xlab = "", cex.axis = 3,
  cex.lab = 4)
lines(t, mspi4, col = "red", lty = 2,lwd = 10)
mtext(text = "Time", side = 1, line = 3.7, las = 1, cex = 3, at = 0.5)
grid(nx = NULL, ny = NULL, col = "gray50", lty = "solid",
  lwd = par("lwd"), equilogs = TRUE)
box(col = 'gray50')
dev.off()


