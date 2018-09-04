$PARAM FFM=40, BMI=20, FORCEETA1=0, FORCEETA2=0

$SET delta=1/24/6

$CMT ABS1 ABS2 CENT

$MAIN
double POPKEL = 1.32 * 24 * pow(FFM/70, -0.25); // times in days!!
double POPV = 72.2 * FFM/56.1; // 56.1 = FFM of average (70kg) person
double POPT12ABS = 0.455; // times in days
double POPT12DIS = 8.88 * (BMI>25 ? 1.865 : 1);

double INDKEL = POPKEL;
double INDV = POPV * exp(ETA(2) + FORCEETA2);
double INDT12ABS = POPT12ABS * exp((ETA(1) + FORCEETA1) * -1.24);
double INDT12DIS = POPT12DIS * exp(ETA(1) + FORCEETA1);

double KEL = INDKEL;
double KA = log(2)/INDT12ABS;
double KD = log(2)/INDT12DIS;

$ODE
dxdt_ABS1 = -KA*ABS1;
dxdt_ABS2 = KA*ABS1 - KD*ABS2;
dxdt_CENT = KD*ABS2 - KEL*CENT;

$OMEGA >> block=TRUE
0.4842
-0.132 0.0648

$SIGMA 0.1204

$TABLE
double DV = (CENT/INDV)*exp(EPS(1)); // observed, only fuzzes when we set sigma
double ET1 = ETA(1);
double ET2 = ETA(2);

$CAPTURE DV ET1 ET2 POPKEL INDKEL
