library(mrgsolve)
library(magrittr)


# Load somemodels. To get the fixed effects right, their respective covariates
# *must* be param'ed in

# for(filename in list.files(pattern="\\.cpp$")) {
#   noext <- sub("\\.\\w*$", "", filename)
#   assign(paste("mod.", noext, sep=""), mread(noext))
# }

mod.vanc_roberts2011 <- mread("vanc_roberts2011")
mod.vanc_thomson2009 <- mread("vanc_thomson2009")


# Note that dosing is not mentioned below. It is the caller's responsibility.


# Mutate a model to suit an individual's etas

ipred <- function(mod, eta) {
  # "individual prediction"
  # mutate mod with the provided etas and zero omega matrix
  
  # zero matrices
  zomat <- (mod %>% omat %>% as.matrix) * 0
  mod %<>% omat(zomat)
  
  # set etas
  etalist <- as.list(eta)
  names(etalist) <- paste('FORCEETA', seq(length(eta)), sep='')
  mod %<>% param(etalist)
  
  mod
}


# Objective function

bayesian.ofv <- function(eta, mod, t, y) {
  # simulate ipred, no RUV, at just the times we ask
  y2 <- mod %>% ipred(eta) %>% smat(dmat(0)) %>% obsonly %>% mrgsim(end=-1, add=t)
  y2 <- y2$DV

  # extract a few bits from the model
  sigma.num <- mod %>% smat %>% as.matrix %>% as.numeric
  omega.mat <- mod %>% omat %>% as.matrix
  eta.mat.horiz <- eta %>% matrix(nrow=1)

  # calculate the RUV part of the objective function
  actual.variances <- y2^2 * sigma.num # sig2j
  sqwres <- log(actual.variances) + (y2-y)^2/actual.variances
  
  # and the ETA part of the objective function
  # what's happening here?
  n0n <- diag(eta.mat.horiz %*% solve(omega.mat) %*% t(eta.mat.horiz))
  
  sum(sqwres) + n0n
}


remove.mod.uncertainty <- function(mod) {
  zeromat <- mod %>% omat %>% as.matrix * 0
  mod %<>% omat(zeromat)
  
  zeromat <- mod %>% smat %>% as.matrix * 0
  mod %<>% smat(zeromat)
}


# Change the model to the population of people indistinguishable from our individual

hack.mod.for.fit <- function(mod, fit)
{
  # force our estimated etas on the model
  etalist <- as.list(fit$par)
  names(etalist) <- paste('FORCEETA', seq(length(etalist)), sep='')

  mod %<>% param(etalist)
  
  # misuse the omega matrix 
  matrix.describing.uncertainty.in.our.estimates.etas <- solve(fit$hessian)

  mod %<>% omat(matrix.describing.uncertainty.in.our.estimates.etas)
}


# Fit function, which returns my "suspected" etas and their Hessian

tdm <- function(mod, t, y)
{
  n.eta <- mod %>% omat %>% as.matrix %>% nrow
  init.eta <- rep(0, n.eta)
  
  optim(init.eta, bayesian.ofv, hessian=T, mod=mod, t=t, y=y)
}
