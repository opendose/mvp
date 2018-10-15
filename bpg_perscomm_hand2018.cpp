$PROB
Unpublished model

$PARAM FFM=40, OVWT=0

$CMT ABS1 ABS2 CENT

$MAIN
double POPKEL    = 1.32 * 24 * pow(FFM/70, -0.25); // times in days!!
double POPV      = 72.2 * FFM/56.1; // 56.1 = FFM of average (70kg) person
double POPT12PRI = 0.455; // times in days
double POPT12SEC = 8.88 * (OVWT ? 1.865 : 1);

double INDKEL    = POPKEL;
double INDV      = POPV      * exp(ETA(2));
double INDT12PRI = POPT12PRI * exp(ETA(1) * -1.24);
double INDT12SEC = POPT12SEC * exp(ETA(1));

$ODE
dxdt_ABS1 = -(log(2)/INDT12PRI)*ABS1;
dxdt_ABS2 = (log(2)/INDT12PRI)*ABS1 - (log(2)/INDT12SEC)*ABS2;
dxdt_CENT = (log(2)/INDT12SEC)*ABS2 - INDKEL*CENT;

$TABLE
double F = CENT/INDV;
double FLOGVAR = 0.1204;
double Y = F * exp(EPS(1));

$CAPTURE POPKEL, POPV, POPT12PRI, POPT12SEC, INDKEL, INDV, INDT12PRI, INDT12SEC, F, FLOGVAR, Y, 

$OMEGA @block
0.4842
-0.132 0.0648

$SIGMA 0.1204
