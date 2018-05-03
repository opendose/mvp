library(mrgsolve)
library(magrittr)



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



code <- '
$PARAM TBW=74.8, CRCL=1, FORCEETA1=0, FORCEETA2=0

$CMT CENT

$PKMODEL ncmt=1

$OMEGA 0.139876 0.151321
$SIGMA 0.039601

$MAIN
double V = TBW * 1.53 * (1 + FORCEETA1 + ETA(1));
double CL = CRCL * 4.58 * (1 + FORCEETA2 + ETA(2));

$TABLE
double DV = (CENT/V)*(1+EPS(1)); // observed, only fuzzes when we set sigma
double ET1 = ETA(1);
double ET2 = ETA(2);

$CAPTURE DV ET1 ET2 CL
'

mod <- mcode("roberts", code)



# LET'S TALK TDM: need to know inputs and outputs

# INPUTS
doses <- c(
  ev(time=0, amt=1500)
)

# OUTPUTS (what you measured and when you measured it)

# vector of times that TDM measurements were taken
# the values of those measurements

t <- c(12)
y <- c(20)



mod %<>% ev(doses)
# here I should also resolve the fixed effects of the model




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




fit <- optim(c(0.1,0.1), bayesian.ofv, hessian=T, mod=mod, t=t, y=y)



# okay, now to graph the bastard...

