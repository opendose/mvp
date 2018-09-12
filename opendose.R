library(mvtnorm) # multivariate normal distribution

# Extract a line of R source starting with "OpenDose_" from an mrgsolve model
opendose_modextractR <- function(mod, name) {
  matchre <- sprintf("OpenDose_%s=", name)
  cutre <- sprintf(".*%s", matchre)
  match <- grep(matchre, mod@code, value=T)
  stopifnot(length(match) == 1)
  cut <- sub(cutre, "", match)
  eval(parse(text=cut))
}

# Rewrite an mrgsolve model to take its ETAs from PARAMs called ETA1...n,
# instead of the random variables ETA(1...n)
opendose_modrewriteETA <- function(mod) {
  etacount <- nrow(omat(mod))
  code <- mod@code
  code <- gsub("\\bETA\\(([0-9]+)\\)", "ETA\\1", code)
  etanames <- sprintf("ETA%d", seq(etacount))
  etadecl <- paste(sprintf("%s=0", etanames), collapse=", ")
  code <- sub("^(\\$PARAM.*)", sprintf("\\1, %s", etadecl), code)
  code <- paste(code, collapse="\n")
  mcode("mynewmodel", code)
}

opendose_residual_zscores <- function(sim) {
  if("FLOGVAR" %in% colnames(sim)) {
    residuals.zscores <- log(sim$Y/sim$F)/sqrt(sim$FLOGVAR)
  } else {
    residuals.zscores <- (sim$Y-sim$F)/sqrt(sim$FVAR)
  }
}

# Perform a Monte Carlo simulation to calculate the a posteriori
# probability distribution of the ETAs of an mrgsolve model
# (Depends on an OpenDose_ScaledResidualPDF function in the model code)
opendose_gamble <- function(mod, tdmT, tdmY, methods="ebe") {
  stopifnot(all(methods %in% c("ebe", "montecarlo_posterior", "montecarlo_prior")))
  stopifnot(length(methods) > 0)

  omega <- mod %>% omat %>% as.matrix
  etacount <- nrow(omega)
  etanames <- paste("ETA", seq(etacount), sep="")
  
  output <- list()

  if("ebe" %in% methods) {
    # minimise (eta pdf) * (residual pdf)
    set.seed(1)

    ebe.solution <- optim(rep(0, nrow(omega)), function(etavec) {
      eta.logpdf <- dmvnorm(etavec, sigma=omega, log=T)
      
      etalist <- as.list(etavec)
      names(etalist) <- etanames
      sim <- mod %>% param(etalist) %>% mrgsim(end=-1, add=tdmT) %>% filter(time %in% tdmT) %>% as.data.frame
      sim$Y <- tdmY
      residuals.zscores <- opendose_residual_zscores(sim)
      residuals.logpdf <- dnorm(residuals.zscores, log=T)
      
      -(eta.logpdf + sum(residuals.logpdf))
    })
    
    ebe <- ebe.solution$par %>% as.list
    names(ebe) <- etanames
    output$ebe <- ebe
  }
  
  if(any(c("montecarlo_posterior", "montecarlo_prior") %in% methods)) {
    set.seed(1)
    df <- rmvnorm(1000, sigma=omega) %>% as.data.frame
    colnames(df) <- etanames

    if("montecarlo_prior" %in% methods) output$montecarlo_prior <- df

    if("montecarlo_posterior" %in% methods) {
      likelihood <- apply(df, 1, function(etalist) { # iterate over rows and calculate a likelihood for each
        names(etalist) <- etanames
        sim <- mod %>% param(etalist) %>% mrgsim(end=-1, add=tdmT) %>% filter(time %in% tdmT) %>% as.data.frame
        sim$Y <- tdmY
        residuals.zscores <- opendose_residual_zscores(sim)
        residuals.pdf <- dnorm(residuals.zscores) / dnorm(0) # this might overestimate how reasonable the TDM figure is
        prod(residuals.pdf) # multiply the probabilities together!
      })
      
      keepvector <- rbinom(n=length(likelihood), size=1, prob=likelihood)
      df <- df[keepvector==1, ]
      output$montecarlo_posterior <- df
    }
  }
  
  output
}
