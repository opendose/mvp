$MODEL
COMP=(ABS1,DEFDOSE) ; 1st absorption compartment
COMP=(ABS2) ; 2nd absorption compartment
COMP=(PLASMA) ; PLASMA COMPARTMENT

$PK
KEL=THETA(1)*((FFM/70)**(-0.25))
V=THETA(2)*(FFM/56.1)*EXP(ETA(2))
ETARAT=THETA(6)
ABS=THETA(3)*EXP(ETA(1)*ETARAT*(-1))    
KA=LOG(2)/(ABS*24)
OBES=0                           
IF(BMI.GT.25) OBES=1                           
DIS=THETA(4)*EXP(ETA(1)) *(1+((OBES)*THETA(5)))
KD=LOG(2)/(DIS*24)
S2=V
F1=1
ALAG1=0

$DES
DADT(1)= -KD*A(1)
DADT(2)=KD*A(1)-KA*A(2)
DADT(3)= KA*A(2) - KEL*A(3)

$ERROR
IPRED=LOG(F+0.001)
IRES=DV-IPRED
IWRES=IRES/1
Y=LOG(F)+ERR(1)

$THETA
(0, 1.32) FIX ;KEL 
(0, 72.2) ; V2
(0, 0.455) ; KA
(0, 8.88) ; KD
(0, 0.865) ; KD-OBES
(0, 1.24) ; ETARAT

$OMEGA BLOCK(2)
 0.4842  ; KD
 -0.132   ; r(kD, V)
0.0648  ; V

$SIGMA 0.1204  ; Properr
