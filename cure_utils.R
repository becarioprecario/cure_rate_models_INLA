# Common functions to implement cure models with INLA

####################
# USEFUL FUNCTIONS #
####################

# CURE PROPORTION (eta)  

#x1: Design matrix of covariates
#beta1: Column-vector of coefficients
cureprop <- function(x1, beta1){
  mue <- (x1 %*% beta1)[, 1]
  p <- exp(mue) / (1 + exp(mue))
  return(p)
}

#SURVIVAL UNCURED (Su(t))
#Note that in INLA Weibull likelihood  parameterisation
#contains the "true" lambda scale parameter inside the
#linear predictor as the intercept term.

#time: Time
#alpha: Alpha parameter of weibull
#x2: Design matrix of covariates
#beta2: Column-vector of coefficients
survunc <- function(time, alpha, x2, beta2){
  mus <- (x2 %*% beta2)[, 1]
  lambda <- exp(mus)
  s <- exp(-lambda * time^alpha)
  return(s)
}

# FULL CONDITIONAL DISTRIBUTION OF THE LATENT VARIABLE Z
#d: not censoring indicator variable

#delta: Censorship indicator (1=observed, 0=censored)
#eta: Lineaer predictor(?)
#su: Survival function (?)
piz<-function(delta, eta, su){
  A <- (1 - eta) * su
  B <- eta + ( (1 - eta) * su)
  aux <- (1 - (A / B)) + delta * ((A / B) - 1)
  return(aux)
}

# SAMPLE Z
#z=1 cured 
#z=0 uncured

#w: weights for sampling, i.e., w = Pr(z = 1)
upz <- function(w){
  aux <- rep(NA, length(w))
  for(i in 1:length(w)){
    aux[i] <- sample(0:1, 1, prob = c(1 - w[i], w[i]))
  }
  return(aux)
}


