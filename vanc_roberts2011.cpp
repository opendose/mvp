$PARAM TBW=74.8, CRCL=90.7, FORCEETA1=0, FORCEETA2=0

$CMT CENT

$PKMODEL ncmt=1

$OMEGA 0.139876 0.151321
$SIGMA 0.039601

$MAIN
double V = TBW * 1.53 * exp(FORCEETA1 + ETA(1));
double CL = CRCL * 0.0458 * exp(FORCEETA2 + ETA(2));

$TABLE
double DV = (CENT/V)*(1+EPS(1)); // observed, only fuzzes when we set sigma
double ET1 = ETA(1);
double ET2 = ETA(2);

$CAPTURE DV ET1 ET2 CL
