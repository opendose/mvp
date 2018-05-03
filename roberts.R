library(magrittr)
library(mrgsolve)



# LET'S TALK TDM: need to know inputs and outputs

# INPUTS
doses <- c(
  ev(time=0, amt=1000),
  ev(time=12, amt=1000)
)

# OUTPUTS (what you measured and when you measured it)

# vector of times that TDM measurements were taken
# the values of those measurements

x <- c(12, 24, 36)
y <- c(15, 20, 22)



# Create some models...

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

$CAPTURE DV ET1 ET2
'

mod <- mcode("roberts", code)
mod.fe <- mod %>% ev(doses) # add fixed effects in here later (they are currently not calculated from the model)

omega <- omat(mod) %>% as.matrix
sigma <- smat(mod) %>% as.matrix


objective <- function(eta, mod, x, y) {
  # remove omega and sigma matrices
  mod %<>% omat(cmat(0,0,0)) %>% smat(cmat(0))
  
  # add our own etas instead
  mod %<>% param(FORCEETA1=eta[1], FORCEETA2=eta[2])

  # get observation records
  y.sim <- mod %>% mrgsim(end=-1, add=x)
  y.sim <- y.sim$DV
  
  # now calculate an objective function (eek!)
  eta_m <- matrix(eta, nrow=1)

  # RUV component of objective function
  sig2j <- y.sim^2 * as.numeric(sigma)
  sqwres <- log(sig2j) + (1/sig2j) * (y.sim)^2

  # ETA component of objective function
  n0n <- diag(eta_m %*% solve(omega) %*% t(eta_m))
  
  # Sum them
  sum(sqwres) + n0n
}

fit <- newuoa(c(1,1), objective, mod=mod.fe, x=x, y=y)
