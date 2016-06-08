************************************************************************************************** ;
*	ASSESSING NATURAL DIRECT AND INDIRECT EFFECTS FOR A CONTINUOUS EXPOSURE AND A                  ;
*	DICHOTOMOUS OUTCOME                                                                            ;
*   Version 1.1, Uploaded 06/08/2016        										                   ;
*																					               ;
*	Code written by Wei Wang (2015)													               ;
* 																					               ;
*	Reference:																		               ;
* 	Wang W, and Zhang B (2016). Assessing Natural Direct and Indirect Effects for a Continuous     ; 
*   Exposure and a Dichotomous Outcome. Journal of Statistical Theory and Practice. Submitted.     ;
* 																					               ;
***************************************************************************************************;
*	PROGRAM DESCRIPTION:															               ;
*   -----------------------------------------------------------------------------------------------;
*	THIS PROGRAM PRODUCES THE ESTIMATION OF THE NATURAL INDIRECT EFFECT, DIRECT                    ; 
*   EFFECT AND TOTAL EFFECT FOR A CONTINUOUS EXPOSURE AND A DICHOTOMOUS OUTCOME IN                 ; 
*   TWO-STAGE MEDIATION MODELS. IT APPLIED AND WORKED FOR ONE CONTINUOUS EXPOSURE, ONE             ; 
*   DICHOTOMOUS OUTCOME AND ONE BINARY OR CONTINUOUS MEDIATOR AND ANY NUMBER OF                    ; 
*   BASELINE COVARIATES WHICH MAY AFFECT THE MEDIATOR AND THE DICHOTOMOUS OUTCOME. THE NATURAL     ; 
*	DIRECT AND INDIRECT EFFECTS ON BOTH THE RISK DIFFERENCE AND ODDS RATIO SCALES WERE DEFINED     ;
*	AND ESTIMATED UTILIZING THE EMPIRICAL JOINT DISTRIBUTION OF THE EXPOSURE AND BASELINE   	   ;
*	COVARIATES FROM THE WHOEL SAMPLE ANALYSIS POPULATION.       								   ;
*	                                          													   ;
*	DATASET REQUIREMENTS:															               ;
*   -----------------------------------------------------------------------------------------------;
*	THE DATASET MUST HAVE ONE LINE PER SUBJECT WHERE EACH SUBJECT MUST CONTAIN   	               ;
*	ONE CONTINUOUS EXPOSURE VARIABLE, ONE MEDIATOR, EITHER BINARY OR                               ;
*   CONTINUOUS, ONE DICHOTOMOUS OUTCOME AND ANY NUMBER OF  BASELINE COVARIATES. OF NOTE,           ;
*   BASELINE COVARIATE CAN BE DISCRETE OR CONTINUOUS.                                              ;
*                                                                                                  ;
*	THIS MACRO CAN ONLY WORK FOR DATASET WITHOUT MISSING VALUES. THE MISSING VALUES SHOULD BE      ;  
*   IMPUTED OR ANY SUBJECT WITH MISSING DATA SHOULD BE EXCLUDED. THE SAMPLE                        ;
*	DATA SET WITH ONE COVARIATE AND ONE CONTINUOUS MEDIATOR IS LISTED AS FOLLOWING,                ;
*			                                                                                       ;
*			SUBJECT  EXPOSURE  MEDIATOR   OUTCOME  COVARIATE              		                   ;
*				1		21.2		21		  1		   65   					                   ;
*				2		25.6		28		  0		   40   					                   ;
*				3		30.4		25		  1		   37  	    		            	           ;
*				4		35.7		23		  1		   77  	               				           ;
*																					               ;
*																					               ;
*	MODEL SPECIFICATION																               ;
*   -----------------------------------------------------------------------------------------------;
*	X: CONTINUOUS EXPOSURE, M: NORMALIY DISTRIBUTED MEDIATOR OR BINARY MEDIATOR, Y: BINARY OUTCOME ;
*   W: BASELINE COVARIATES THAT AFFECT THE MEDIATOR AND THE DICHOTOMOUS OUTCOME.                   ;
*																					               ;
*	OUTCOME LOGIT(P(Y = 1)) = BETA0 + BETA1 * X + BETA2 * M + BEAT3 * W                            ;
*	FOR CONTINUOUS MEDIATOR M, THE FOLLOWING MEDIATOR MODEL IS SPECIFIED   			               ;
*	M = ALPHA0 + ALPHA1 * X + ALPHA2 * W +  EPSILON	   		                            	       ;
*	WHERE EPSILON~N(0, SIGMASQUARE)	   				                                               ;
*																					               ;
*	FOR CONTINUOUS MEDIATOR	M, THE FOLLOWING MEDIATION MODEL IS SPECIFIED 			               ;
*	LOGIT(P(M = 1)) = ALPHA0 + ALPHA1 * X + ALPHA2 * W	   		                        	       ;
*	                        										    			               ;
*	MACRO VARIABLES:																               ;
*   -----------------------------------------------------------------------------------------------;
*																					               ;
*	DATASET:        MAIN INPUT DATASET												               ;
*																					               ;
*	ANALYSIS: 		NUMERIC VALUE 1-2 INDICATING MEDIATOR TYPE IN THE DATASET                      ;
*												                    							   ;
*					1 = BINARY MEDIATOR                                                            ;
*																					               ;
*					2 = CONTINOUS MEDIATOR           	                                           ;
*																					               ;
*	Y:		SHOULD BE BINARY VARIABLE                                           	               ;
*																					               ;
*	X:		CONTINUOUS EXPOSURE VARIABLE AFTER SCALING IN ORDER TO PROVIDE APPROPRIATE             ; 
*           EXPOSURE EFFET AND MEIDATION EFFECT ESTIMATES CORRESPONDING 'ONE' UNIT INCREASE        ; 
*           OF THE EXPOSURE VARIABLE                                                               ; 
*	                                                                                               ;
*	M:		EITHER BINARY OR NORMALLY DISTRIBUTED CONTINUOUS MEDIATOR IS ALLOWED (SHOULD BE        ;
*           CONSISTENT WITH MACRO VARIABLE ANALYSIS)                                               ;
*																					               ;
*   W:      THE MACRO ALLOWS ANY NUMBER OF BASELINE COVARIATES (CAN BE MISSING), AND ALL SPECIFIED ; 
*           BASELINE COVARIATES WILL AFFECT BOTH MEDIATOR AND OUTCOMES. IN ADDITION, CATEGORICAL   ; 
*           COVARIATES WITH MORE THAN TWO DISTINCT VALUES SHOULD BE TRANSFORMED TO DUMMY VARIABLES ;
*			FIRST BEFORE RUNNING THIS MACRO.            							               ;
*																					               ;
*	OUTRD:  NAME OF THE OUTPUT DATASET CONTAINING THE RISK DIFFERENCE SCALE CAUSAL EFFECTS         ;
*			ESTIMATE. THIS DATASET CONTAINS THE FOLLOWING CAUSAL QUATITIES AND THEIR STANDARD      ;
*			ERRORS FROM DELTA METHOD (FORMULA 17 IN REFERENCE GIVES DEFINITION                     ;
*           AND ESTIMATION METHODS)                                                                ;
*			TE = E(Y(X+1, M(X+1))) - E(Y(X, M(X)))     STANDARD ERROR VARIABLE NAME: GTESD         ;
*			IE0 = E(Y(X, M(X+1))) - E(Y(X, M(X)))	   STANDARD ERROR VARIABLE NAME: GIE0SD        ;
*			DE1 = E(Y(X+1, M(X+1))) - E(Y(X, M(X+1)))  STANDARD ERROR VARIABLE NAME: GDE1SD        ;
*			IE1 = E(Y(X+1, M(X+1))) - E(Y(X+1, M(X)))  STANDARD ERROR VARIABLE NAME: GIE1SD        ;
*			DE0 = E(Y(X+1, M(X))) - E(Y(X, M(X)))      STANDARD ERROR VARIABLE NAME: GDE0SD        ;
*			R1 = IE1/TE								   STANDARD ERROR VARIABLE NAME: GR1SD         ;
*			R0 = IE0/TE								   STANDARD ERROR VARIABLE NAME: GR0SD         ;
*																					               ;
*	OUTOR:  NAME OF THE OUTPUT DATASET CONTAINING THE ODDS RATIO SCALE CAUSAL EFFECTS ESTIMATE.    ;
*			THIS DATASET CONTAINS THE FOLLOWING CAUSAL QUATITIES AND THEIR STANDARD ERRORS FROM    ;
*			DELTA METHOD. OF NOTE, THE FINAL OUTPUT DATA SET PRESENTED LOG ODDS RATIO SCALE        ;
*			CAUSAL QUANTITIES AS FOLLOWING,                                                        ;
*																					               ;
*			LOGTEOR = LOG{E(Y(X+1, M(X+1)))/(1 - E(Y(X+1, M(X+1))))}                               ; 
*                   - LOG{E(Y(X, M(X)))/(1 - E(Y(X, M(X))))}       SE VARIABLE NAME: GTEORSD       ;
*			LOGIE0OR = LOG{E(Y(X, M(X+1)))/(1 - E(Y(X, M(X+1))))}                                  ; 
*                   - LOG{E(Y(X, M(X)))/(1 - E(Y(X, M(X))))}   	   SE VARIABLE NAME: GIE0ORSD      ;
*			LOGDE1OR = LOG{E(Y(X+1, M(X+1)))/(1 - E(Y(X+1, M(X+1))))}                              ; 
*                   - LOG{E(Y(X, M(X+1)))/(1 - E(Y(X, M(X+1))))}   SE VARIABLE NAME: GDE1ORSD      ;
*			LOGIE1OR = LOG{E(Y(X+1, M(X+1)))/(1 - E(Y(X+1, M(X+1))))}                              ; 
*                   - LOG{E(Y(X+1, M(X)))/(1 - E(Y(X+1, M(X))))}   SE VARIABLE NAME: GIE1ORSD      ;
*			LOGDE0OR = LOG{E(Y(X+1, M(X)))/(1 - E(Y(X+1, M(X))))}                                  ; 
*                   - LOG{E(Y(X, M(X)))/(1 - E(Y(X, M(X))))}   	   SE VARIABLE NAME: GDE0ORSD      ;
*																					               ;
*			R1OR = LOGIE1OR/LOGTEOR								   SE VARIABLE NAME: GR1ORSD       ;
*			R0OR = LOGIE0OR/LOGTEOR								   SE VARIABLE NAME: GR0ORSD       ;
*																					               ;
*	EXAMPLE CODE:																	               ;
*   %include 'CONMED.sas' 														                   ;
*																					               ;
*	BINARY MEDIATOR:    															               ;
*   %CONMED(dataset=DSN, analysis = 1, y=RESPONSE, x=PRED, m = MED, w = COV1 COV2, out = OUT1)     ;     
*	    						                 									               ;
*	    											                				               ;
*	CONTINUOUS MEDIATOR:    															           ;
*   %CONMED(dataset=DSN, analysis = 2, y=RESPONSE, x=PRED, m = MED, w = COV1 COV2, out = OUT1)     ; 
*																					               ;
***************************************************************************************************;

options ls=90 ps=64;
dm log 'clear'; dm list 'clear' continue;

***************************************************************************************************;
*************************************Final Macro***************************************************;
***************************************************************************************************;

/* Continuous Mediator Linear Regression Model with PROC IML*/ 

%macro tobit0(dsn=, y=,x=, theta= 0, out1 =, out2=);
proc iml;
reset noname;
use &dsn;
read all var {&x} into x;
read all var {&y} into y;
yx=y || x;
* a = lower truncation point;
* define var indicating censored ohs;
limit = y;
start LL(theta) global (yx, limit);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
sum1 = 0;
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do m = 1 to n0;
   sum1 = sum1 -(1/2) *(log(2*pi) + log(sigma2))
- (y0[m] - (x0[m,]*beta))**2/(2*sigma2);
end;

f = sum1;
return(f);
finish LL ;

* define gradient;
start GRAD(theta) global (yx, limit);
* split data into censored and uncensored ohs;
yx0 = yx[loc((limit ^= .)),];
n0 = nrow(yx0);
k = ncol(yx0);
y0 = yx0[,1];
x0 = j(n0,1,1)|| yx0[,2:k];

*print yx0, n0, k, y0, x0, yx1, n1, y1, x1, yx2, n2, y2, x2, theta;
* separate param vector into beta and sigma;
k = ncol(theta);
k1 = k-1;
beta = theta [1:k1];
sigma = theta [k];
sigma2 = sigma* sigma;
g = j(1,k,0);
pi = arcos(-1);
*print k, k1, beta, sigma, sigma2, sum1;

do m = 1 to n0;
   mu = x0[m,]*beta;

   g[1:k1] = g[1:k1] + (y0[m] - mu)/(max(sigma2, exp(-150))) * t(x0[m,]); 
   g[k] = g[k] - 1/sigma + (y0[m] - mu) ** 2/sigma**3;
end;
return(g);
finish GRAD ;

optn={1 2};
*tc={3000 15000 . 1E-12};
tc={3000 15000 . 1E-6};
*tc={3000 15000};
* define starting va1ues for theta: either user-provided or O1S;
theta0 = &theta;
if theta0 = 0 then
do;
k = ncol(yx) + 1;
n = nrow(y);
theta0 = j(1,k,0);
xx = j(n,1,1) || x;
beta = inv(t(xx)*xx)*t(xx) *y;
e = y - xx*beta;
s = sqrt(ssq(e)/(n- k));
theta0[1:(k-1)] = t(beta);
theta0[k] = s;
print 'OLS starting values for theta = '
theta0;
end;
else print 'User-supplied starting values for
theta = ' theta0;
       con1 = {  .    . , 
                .    . };
*				theta = {1 1 1 1};
t1 = j(1,ncol(theta0),.);
t2 = t1;
t2[ncol(theta0)] = 0;
t3 = t2//t1;
con = t3;
*call nlpqn(rc, theta, 'LL',theta0,optn, con, tc);
call nlpqn(rc, theta, 'LL',theta0,optn, con, tc) grd='GRAD';
call nlpfdd(f, g, h,'LL', theta);
var = inv(-h);
sd = t(j(1,ncol(theta0),.));
temp = vecdiag(var);
do i = 1 to nrow(temp);
if temp[i] >= 0 then sd[i] = sqrt(temp[i]);
else sd[i] = .;
end; 
*sd = sqrt(vecdiag(var));
print 'Hessian = ' h ' Covariance matrix = ' var
' Standard errors = ' sd;
*print 'Hessian = ' h ' Covariance matrix = ' var
;
print rc theta;

postprobs=theta;

create &out1 from postprobs;
append from postprobs;

postprobs1=var;

create &out2 from postprobs1;
append from postprobs1;

quit;
run;
%mend;

/* Continous exposure and continuous mediator with at least one covariate */

%macro conmedcon(dataset=' ', y=' ', x=' ', m = '', w = '', out =' ');

proc logistic data = &dataset;
model &y (event = '1') = &x &m &w/COVB;
ods output covb = aa1
ParameterEstimates = aa2;
run; 

proc transpose data = aa2 out = aa3;
id variable;
var estimate;
run;

data aa4;
merge ff11 aa3;
drop _name_;
run;

data ff13;
set ff12;
id = _n_;
run;

proc means data = ff12;
var col1;
output out = aa11 n(col1) = np;
run;

data aa12;
set aa11;
nb = 2 * np;
run;

data _null_;
set aa11;
call symput('NP', np);
call symput('NB', nb);
run;

data aa21;
set aa1;
id = _n_ + &NP;
run;

data aa5;
merge ff13 aa21;
by id;
drop id parameter;
run;

data aa6;
set aa5;
array x[*] col1 -- &w;
do i = 1 to dim(x);
if x[i] = . then x[i] = 0;
end;
drop i;
run;

proc iml;
reset noname;
use &dataset;
read all var {&x &w} into x;
use aa6;
read all into var;
use abswei;
read all var {abs weight} into w;
use aa4;
read all into theta;

n = nrow(x);
k = ncol(theta)/2;
alpha = theta [1:k-1];
sigma = theta[k];
beta = theta[k+1:2*k];
nw = ncol(x) - 1;
mx0 = j(n,1,1) || x;
mx1 = j(n,1,1) || (x[,1] + 1) || x[,2:(nw+1)];

yx0 = j(n,1,1) || x[,1] || j(n,1,0) || x[,2:(nw+1)];
yx1 = j(n,1,1) || (x[,1] + 1) || j(n,1,0) || x[,2:(nw+1)];
/*x2 = yx0[1:3,];*/
/*x3 = yx1[1:3,];*/
/**/
/*print theta k alpha beta x2 x3;*/
sum11 = 0;
sum10 = 0;
sum01 = 0;
sum00 = 0;

TE = 0;
IE1 = 0;
DE0 = 0;
IE0 = 0;
DE1 = 0;

TEOR = 0;
IE1OR = 0;
DE0OR = 0;
IE0OR = 0;
DE1OR = 0;

gte = j(1,k*2,0);
gie1 = j(1,k*2,0);
gde0 = j(1,k*2,0);
gie0 = j(1,k*2,0);
gde1 = j(1,k*2,0);

gteor = j(1,k*2,0);
gie1or = j(1,k*2,0);
gde0or = j(1,k*2,0);
gie0or = j(1,k*2,0);
gde1or = j(1,k*2,0);

do i = 1 to n;
   mum1 = mx1[i,]*alpha; 
   mum0 = mx0[i,]*alpha; 
   muy11 = yx1[i,]*beta + beta[3] * (sqrt(2) * sigma * w[,1] + mum1); 
   muy10 = yx1[i,]*beta + beta[3] * (sqrt(2) * sigma * w[,1] + mum0); 
   muy01 = yx0[i,]*beta + beta[3] * (sqrt(2) * sigma * w[,1] + mum1); 
   muy00 = yx0[i,]*beta + beta[3] * (sqrt(2) * sigma * w[,1] + mum0); 

   E11 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy11))))* w[,2]; 
   E10 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy10))))* w[,2]; 
   E01 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy01))))* w[,2]; 
   E00 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy00))))* w[,2]; 

   sum11 = sum11 + E11;
   sum10 = sum10 + E10;
   sum01 = sum01 + E01;
   sum00 = sum00 + E00;

   TE = TE + (E11 - E00);
   IE1 = IE1 + (E11 - E10);
   DE0 = DE0 + (E10 - E00);
   IE0 = IE0 + (E01 - E00);
   DE1 = DE1 + (E11 - E01);

   TEOR = TEOR +  LOG((E11/(1-E11))/(E00/(1-E00)));
   IE1OR = IE1OR + LOG((E11/(1-E11))/(E10/(1-E10)));
   DE0OR = DE0OR + LOG((E10/(1-E10))/(E00/(1-E00)));
   IE0OR = IE0OR + LOG((E01/(1-E01))/(E00/(1-E00)));
   DE1OR = DE1OR + LOG((E11/(1-E11))/(E01/(1-E01)));

/* Derivative for alpha coefficient */

   E11_1 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy11))) # (sqrt(2) * w[,1]/sigma))* w[,2]; 
   E10_1 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy10))) # (sqrt(2) * w[,1]/sigma))* w[,2]; 
   E01_1 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy01))) # (sqrt(2) * w[,1]/sigma))* w[,2]; 
   E00_1 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy00))) # (sqrt(2) * w[,1]/sigma))* w[,2]; 

/* Derivative for sigma */

   E11_2 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy11))) # (2 * w[,1] # w[,1]/sigma - 1/sigma))* w[,2]; 
   E10_2 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy10))) # (2 * w[,1] # w[,1]/sigma - 1/sigma))* w[,2]; 
   E01_2 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy01))) # (2 * w[,1] # w[,1]/sigma - 1/sigma))* w[,2]; 
   E00_2 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy00))) # (2 * w[,1] # w[,1]/sigma - 1/sigma))* w[,2]; 

/* Derivative for other beta except for beta2 */

   E11_3 = t(1/sqrt(constant('pi')) * (exp((-1) * muy11)/(1 + exp((-1) * muy11))/(1 + exp((-1) * muy11)) ))* w[,2]; 
   E10_3 = t(1/sqrt(constant('pi')) * (exp((-1) * muy10)/(1 + exp((-1) * muy10))/(1 + exp((-1) * muy10)) ))* w[,2]; 
   E01_3 = t(1/sqrt(constant('pi')) * (exp((-1) * muy01)/(1 + exp((-1) * muy01))/(1 + exp((-1) * muy01)) ))* w[,2]; 
   E00_3 = t(1/sqrt(constant('pi')) * (exp((-1) * muy00)/(1 + exp((-1) * muy00))/(1 + exp((-1) * muy00)) ))* w[,2]; 

/* Derivative for beta2 */

   E11_4 = t(1/sqrt(constant('pi')) * (exp((-1) * muy11)/(1 + exp((-1) * muy11))/(1 + exp((-1) * muy11)) # (sqrt(2) * w[,1] * sigma + mum1) ))* w[,2]; 
   E10_4 = t(1/sqrt(constant('pi')) * (exp((-1) * muy10)/(1 + exp((-1) * muy10))/(1 + exp((-1) * muy10)) # (sqrt(2) * w[,1] * sigma + mum0) ))* w[,2]; 
   E01_4 = t(1/sqrt(constant('pi')) * (exp((-1) * muy01)/(1 + exp((-1) * muy01))/(1 + exp((-1) * muy01)) # (sqrt(2) * w[,1] * sigma + mum1) ))* w[,2]; 
   E00_4 = t(1/sqrt(constant('pi')) * (exp((-1) * muy00)/(1 + exp((-1) * muy00))/(1 + exp((-1) * muy00)) # (sqrt(2) * w[,1] * sigma + mum0) ))* w[,2]; 

   gte[1:(k-1)] = gte[1:(k-1)] + E11_1 * t(mx1[i,]) - E00_1* t(mx0[i,]);
   gte[k] = gte[k] + E11_2 - E00_2;
   gte[(k+1):(k+2)] = gte[(k+1):(k+2)] + E11_3 * t(yx1[i,1:2]) - E00_3 * t(yx0[i,1:2]);
   gte[k+3] = gte[k+3] + E11_4 - E00_4;
   gte[(k+4):2*k] = gte[(k+4):2*k] + E11_3 * t(yx1[i,4:k]) - E00_3 * t(yx0[i,4:k]);

   gie1[1:(k-1)] = gie1[1:(k-1)] + E11_1 * t(mx1[i,]) - E10_1* t(mx0[i,]);
   gie1[k] = gie1[k] + E11_2 - E10_2;
   gie1[(k+1):(k+2)] = gie1[(k+1):(k+2)] + E11_3 * t(yx1[i,1:2]) - E10_3 * t(yx1[i,1:2]);
   gie1[k+3] = gie1[k+3] + E11_4 - E10_4;
   gie1[(k+4):2*k] = gie1[(k+4):2*k] + E11_3 * t(yx1[i,4:k]) - E10_3 * t(yx1[i,4:k]);

   gde0[1:(k-1)] = gde0[1:(k-1)] + E10_1 * t(mx0[i,]) - E00_1* t(mx0[i,]);
   gde0[k] = gde0[k] + E10_2 - E00_2;
   gde0[(k+1):(k+2)] = gde0[(k+1):(k+2)] + E10_3 * t(yx1[i,1:2]) - E00_3 * t(yx0[i,1:2]);
   gde0[k+3] = gde0[k+3] + E10_4 - E00_4;
   gde0[(k+4):2*k] = gde0[(k+4):2*k] + E10_3 * t(yx1[i,4:k]) - E00_3 * t(yx0[i,4:k]);

   gie0[1:(k-1)] = gie0[1:(k-1)] + E01_1 * t(mx1[i,]) - E00_1* t(mx0[i,]);
   gie0[k] = gie0[k] + E01_2 - E00_2;
   gie0[(k+1):(k+2)] = gie0[(k+1):(k+2)] + E01_3 * t(yx0[i,1:2]) - E00_3 * t(yx0[i,1:2]);
   gie0[k+3] = gie0[k+3] + E01_4 - E00_4;
   gie0[(k+4):2*k] = gie0[(k+4):2*k] + E01_3 * t(yx0[i,4:k]) - E00_3 * t(yx0[i,4:k]);

   gde1[1:(k-1)] = gde1[1:(k-1)] + E11_1 * t(mx1[i,]) - E01_1* t(mx1[i,]);
   gde1[k] = gde1[k] + E11_2 - E01_2;
   gde1[(k+1):(k+2)] = gde1[(k+1):(k+2)] + E11_3 * t(yx1[i,1:2]) - E01_3 * t(yx0[i,1:2]);
   gde1[k+3] = gde1[k+3] + E11_4 - E01_4;
   gde1[(k+4):2*k] = gde1[(k+4):2*k] + E11_3 * t(yx1[i,4:k]) - E01_3 * t(yx0[i,4:k]);

   gteor[1:(k-1)] = gteor[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E00 - E11 * E00_1* t(mx0[i,]))*(E00 - E00 * E11)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E11 - E00 * E11_1 * t(mx1[i,]))*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[k] = gteor[k] + ((E11_2 - E11_2 * E00 - E11 * E00_2)*(E00 - E00 * E11)-(E00_2 - E00_2 * E11 - E00 * E11_2)*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[(k+1):(k+2)] = gteor[(k+1):(k+2)] + ((E11_3 * t(yx1[i,1:2]) - E11_3 * t(yx1[i,1:2]) * E00 - E11 * E00_3 * t(yx0[i,1:2]))*(E00 - E00 * E11)-(E00_3 * t(yx0[i,1:2]) - E00_3 * t(yx0[i,1:2]) * E11 - E00 * E11_3 * t(yx1[i,1:2]))*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[k+3] = gteor[k+3] + ((E11_4 - E11_4 * E00 - E11 * E00_4)*(E00 - E00 * E11)-(E00_4 - E00_4 * E11 - E00 * E11_4)*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[(k+4):2*k] = gteor[(k+4):2*k] + ((E11_3 * t(yx1[i,4:k]) - E11_3 * t(yx1[i,4:k]) * E00 - E11 * E00_3 * t(yx0[i,4:k]))*(E00 - E00 * E11)-(E00_3 * t(yx0[i,4:k]) - E00_3 * t(yx0[i,4:k]) * E11 - E00 * E11_3 * t(yx1[i,4:k]))*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;

   gie1or[1:(k-1)] = gie1or[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E10 - E11 * E10_1* t(mx0[i,]))*(E10 - E10 * E11)-(E10_1* t(mx0[i,]) - E10_1* t(mx0[i,]) * E11 - E10 * E11_1 * t(mx1[i,]))*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[k] = gie1or[k] + ((E11_2 - E11_2 * E10 - E11 * E10_2)*(E10 - E10 * E11)-(E10_2 - E10_2 * E11 - E10 * E11_2)*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[(k+1):(k+2)] = gie1or[(k+1):(k+2)] + ((E11_3 * t(yx1[i,1:2]) - E11_3 * t(yx1[i,1:2]) * E10 - E11 * E10_3 * t(yx1[i,1:2]))*(E10 - E10 * E11)-(E10_3 * t(yx1[i,1:2]) - E10_3 * t(yx1[i,1:2]) * E11 - E10 * E11_3 * t(yx1[i,1:2]))*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[k+3] = gie1or[k+3] + ((E11_4 - E11_4 * E10 - E11 * E10_4)*(E10 - E10 * E11)-(E10_4 - E10_4 * E11 - E10 * E11_4)*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[(k+4):2*k] = gie1or[(k+4):2*k] + ((E11_3 * t(yx1[i,4:k]) - E11_3 * t(yx1[i,4:k]) * E10 - E11 * E10_3 * t(yx1[i,4:k]))*(E10 - E10 * E11)-(E10_3 * t(yx1[i,4:k]) - E10_3 * t(yx1[i,4:k]) * E11 - E10 * E11_3 * t(yx1[i,4:k]))*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;

   gde0or[1:(k-1)] = gde0or[1:(k-1)] + ((E10_1 * t(mx0[i,]) - E10_1 * t(mx0[i,]) * E00 - E10 * E00_1* t(mx0[i,]))*(E00 - E00 * E10)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E10 - E00 * E10_1 * t(mx0[i,]))*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[k] = gde0or[k] + ((E10_2 - E10_2 * E00 - E10 * E00_2)*(E00 - E00 * E10)-(E00_2 - E00_2 * E10 - E00 * E10_2)*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[(k+1):(k+2)] = gde0or[(k+1):(k+2)] + ((E10_3 * t(yx1[i,1:2]) - E10_3 * t(yx1[i,1:2]) * E00 - E10 * E00_3 * t(yx0[i,1:2]))*(E00 - E00 * E10)-(E00_3 * t(yx0[i,1:2]) - E00_3 * t(yx0[i,1:2]) * E10 - E00 * E10_3 * t(yx1[i,1:2]))*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[k+3] = gde0or[k+3] + ((E10_4 - E10_4 * E00 - E10 * E00_4)*(E00 - E00 * E10)-(E00_4 - E00_4 * E10 - E00 * E10_4)*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[(k+4):2*k] = gde0or[(k+4):2*k] + ((E10_3 * t(yx1[i,4:k]) - E10_3 * t(yx1[i,4:k]) * E00 - E10 * E00_3 * t(yx0[i,4:k]))*(E00 - E00 * E10)-(E00_3 * t(yx0[i,4:k]) - E00_3 * t(yx0[i,4:k]) * E10 - E00 * E10_3 * t(yx1[i,4:k]))*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;

   gie0or[1:(k-1)] = gie0or[1:(k-1)] + ((E01_1 * t(mx1[i,]) - E01_1 * t(mx1[i,]) * E00 - E01 * E00_1* t(mx0[i,]))*(E00 - E00 * E01)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E01 - E00 * E01_1 * t(mx1[i,]))*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[k] = gie0or[k] + ((E01_2 - E01_2 * E00 - E01 * E00_2)*(E00 - E00 * E01)-(E00_2 - E00_2 * E01 - E00 * E01_2)*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[(k+1):(k+2)] = gie0or[(k+1):(k+2)] + ((E01_3 * t(yx0[i,1:2]) - E01_3 * t(yx0[i,1:2]) * E00 - E01 * E00_3 * t(yx0[i,1:2]))*(E00 - E00 * E01)-(E00_3 * t(yx0[i,1:2]) - E00_3 * t(yx0[i,1:2]) * E01 - E00 * E01_3 * t(yx0[i,1:2]))*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[k+3] = gie0or[k+3] + ((E01_4 - E01_4 * E00 - E01 * E00_4)*(E00 - E00 * E01)-(E00_4 - E00_4 * E01 - E00 * E01_4)*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[(k+4):2*k] = gie0or[(k+4):2*k] + ((E01_3 * t(yx0[i,4:k]) - E01_3 * t(yx0[i,4:k]) * E00 - E01 * E00_3 * t(yx0[i,4:k]))*(E00 - E00 * E01)-(E00_3 * t(yx0[i,4:k]) - E00_3 * t(yx0[i,4:k]) * E01 - E00 * E01_3 * t(yx0[i,4:k]))*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;

   gde1or[1:(k-1)] = gde1or[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E01 - E11 * E01_1* t(mx1[i,]))*(E01 - E01 * E11)-(E01_1* t(mx1[i,]) - E01_1* t(mx1[i,]) * E11 - E01 * E11_1 * t(mx1[i,]))*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[k] = gde1or[k] + ((E11_2 - E11_2 * E01 - E11 * E01_2)*(E01 - E01 * E11)-(E01_2 - E01_2 * E11 - E01 * E11_2)*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[(k+1):(k+2)] = gde1or[(k+1):(k+2)] + ((E11_3 * t(yx1[i,1:2]) - E11_3 * t(yx1[i,1:2]) * E01 - E11 * E01_3 * t(yx0[i,1:2]))*(E01 - E01 * E11)-(E01_3 * t(yx0[i,1:2]) - E01_3 * t(yx0[i,1:2]) * E11 - E01 * E11_3 * t(yx1[i,1:2]))*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[k+3] = gde1or[k+3] + ((E11_4 - E11_4 * E01 - E11 * E01_4)*(E01 - E01 * E11)-(E01_4 - E01_4 * E11 - E01 * E11_4)*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[(k+4):2*k] = gde1or[(k+4):2*k] + ((E11_3 * t(yx1[i,4:k]) - E11_3 * t(yx1[i,4:k]) * E01 - E11 * E01_3 * t(yx0[i,4:k]))*(E01 - E01 * E11)-(E01_3 * t(yx0[i,4:k]) - E01_3 * t(yx0[i,4:k]) * E11 - E01 * E11_3 * t(yx1[i,4:k]))*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;

end;

gr1 = (gie1 * te - gte * ie1)/te ** 2;
gr1var = gr1 * var  * t(gr1);
if gr1var >= 0 then gr1sd = sqrt(gr1var);

gr0 = (gie0 * te - gte * ie0)/te ** 2;
gr0var = gr0 * var * t(gr0);
if gr0var >= 0 then gr0sd = sqrt(gr0var);

gr1or = (gie1or * teor - gteor * ie1or)/teor ** 2;
gr1orvar = gr1or * var  * t(gr1or);
if gr1orvar >= 0 then gr1orsd = sqrt(gr1orvar);

gr0or = (gie0or * teor - gteor * ie0or)/teor ** 2;
gr0orvar = gr0or * var  * t(gr0or);
if gr0orvar >= 0 then gr0orsd = sqrt(gr0orvar);

sum11 = sum11/n;
sum10 = sum10/n;
sum01 = sum01/n;
sum00 = sum00/n;

TE = TE/n;
IE1 = IE1/n;
DE0 = DE0/n;
IE0 = IE0/n;
DE1 = DE1/n;

R1 = IE1/TE;
R0 = IE0/TE;

TEOR = TEOR * (1/n);
IE1OR = IE1OR * (1/n);
DE0OR = DE0OR * (1/n);
IE0OR = IE0OR * (1/n);
DE1OR = DE1OR * (1/n);

R1OR = IE1OR/TEOR;
R0OR = IE0OR/TEOR;

gte = gte/n;
gtevar = gte * var * t(gte);
if gtevar >= 0 then gtesd = sqrt(gtevar);
gie1 = gie1/n;
gie1var = gie1 * var * t(gie1);
if gie1var >= 0 then gie1sd = sqrt(gie1var);
gde0 = gde0/n;
gde0var = gde0 * var * t(gde0);
if gde0var >= 0 then gde0sd = sqrt(gde0var);
gie0 = gie0/n;
gie0var = gie0 * var * t(gie0);
if gie0var >= 0 then gie0sd = sqrt(gie0var);
gde1 = gde1/n;
gde1var = gde1 * var * t(gde1);
if gde1var >= 0 then gde1sd = sqrt(gde1var);

gteor = gteor/n;
gteorvar = gteor * var * t(gteor);
if gteorvar >= 0 then gteorsd = sqrt(gteorvar);
gie1or = gie1or/n;
gie1orvar = gie1or * var * t(gie1or);
if gie1orvar >= 0 then gie1orsd = sqrt(gie1orvar);
gde0or = gde0or/n;
gde0orvar = gde0or * var * t(gde0or);
if gde0orvar >= 0 then gde0orsd = sqrt(gde0orvar);
gie0or = gie0or/n;
gie0orvar = gie0or * var * t(gie0or);
if gie0orvar >= 0 then gie0orsd = sqrt(gie0orvar);
gde1or = gde1or/n;
gde1orvar = gde1or * var * t(gde1or);
if gde1orvar >= 0 then gde1orsd = sqrt(gde1orvar);

print sum11 sum10 sum01 sum00 TE IE1 DE0 IE0 DE1 TEOR IE1OR DE0OR IE0OR DE1OR gtesd gie1sd gde0sd gie0sd gde1sd
gteorsd gie1orsd gde0orsd gie0orsd gde1orsd;

postprobs=sum11||sum10||sum01||sum00||TE||IE1||DE0||IE0||DE1||TEOR||IE1OR||DE0OR||IE0OR||DE1OR||gtesd||gie1sd||gde0sd||gie0sd||gde1sd
||gteorsd||gie1orsd||gde0orsd||gie0orsd||gde1orsd||r1||r0||r1or||r0or||gr1sd||gr0sd||gr1orsd||gr0orsd;
 cname = {"sum11" "sum10" "sum01" "sum00" "TE" "IE1" "DE0" "IE0" "DE1" "LOGTEOR" "LOGIE1OR" "LOGDE0OR" "LOGIE0OR" "LOGDE1OR" "gtesd" "gie1sd" "gde0sd" "gie0sd" 
"gde1sd" "gteorsd" "gie1orsd" "gde0orsd" "gie0orsd" "gde1orsd" "r1" "r0" "r1or" "r0or" "gr1sd" "gr0sd" "gr1orsd" "gr0orsd"};

create &out from postprobs  [ colname=cname ];
append from postprobs;

quit;

%mend;

/* Continous exposure and binary mediator with at least one covariate */

%macro conmedbi(dataset=' ', y=' ', x=' ', m = '', w = '', out =' ');

proc logistic data = &dataset;
model &m (event = '1') = &x &w/COVB;
ods output covb = bb1
ParameterEstimates = bb2;
run; 

proc transpose data = bb2 out = bb3;
var estimate;
run;

proc logistic data = &dataset;
model &y (event = '1') = &x &m &w/COVB;
ods output covb = aa1
ParameterEstimates = aa2;
run; 

proc transpose data = aa2 out = aa3;
id variable;
var estimate;
run;

data aa4;
merge bb3 aa3;
drop _name_;
run;

proc transpose data = bb1 out = bb6;
run;

data bb4;
set bb6;
id = _n_;
run;

proc means data = bb4;
var id;
output out = bb5 n(id) = np;
run;

data _null_;
set bb5;
call symput('NP', np);
run;

data aa21;
set aa1;
id = _n_ + &NP;
run;

data aa5;
merge bb4 aa21;
by id;
drop id parameter _name_;
run;

data aa6;
set aa5;
array x[*] col1 -- &w;
do i = 1 to dim(x);
if x[i] = . then x[i] = 0;
end;
drop i;
run;

proc iml;
reset noname;
use &dataset;
read all var {&x &w} into x;
use aa6;
read all into var;
use abswei;
read all var {abs weight} into w;
use aa4;
read all into theta;

n = nrow(x);
k = (ncol(theta) + 1)/2;
alpha = theta [1:k-1];
beta = theta[k:2*k - 1];
nw = ncol(x) - 1;
mx0 = j(n,1,1) || x;
mx1 = j(n,1,1) || (x[,1] + 1) || x[,2:(nw+1)];

yx00 = j(n,1,1) || x[,1] || j(n,1,0) || x[,2:(nw+1)];
yx10 = j(n,1,1) || (x[,1] + 1) || j(n,1,0) || x[,2:(nw+1)];
yx01 = j(n,1,1) || x[,1] || j(n,1,1) || x[,2:(nw+1)];
yx11 = j(n,1,1) || (x[,1] + 1) || j(n,1,1) || x[,2:(nw+1)];
/*x2 = yx0[1:3,];*/
/*x3 = yx1[1:3,];*/
/**/
/*print theta k alpha beta x2 x3;*/
sum11 = 0;
sum10 = 0;
sum01 = 0;
sum00 = 0;

TE = 0;
IE1 = 0;
DE0 = 0;
IE0 = 0;
DE1 = 0;

TEOR = 0;
IE1OR = 0;
DE0OR = 0;
IE0OR = 0;
DE1OR = 0;

gte = j(1,k*2-1,0);
gie1 = j(1,k*2-1,0);
gde0 = j(1,k*2-1,0);
gie0 = j(1,k*2-1,0);
gde1 = j(1,k*2-1,0);

gteor = j(1,k*2-1,0);
gie1or = j(1,k*2-1,0);
gde0or = j(1,k*2-1,0);
gie0or = j(1,k*2-1,0);
gde1or = j(1,k*2-1,0);

do i = 1 to n;
   mum1 = mx1[i,]*alpha; 
   mum0 = mx0[i,]*alpha; 
   muy11 = yx11[i,]*beta; 
   muy10 = yx10[i,]*beta; 
   muy01 = yx01[i,]*beta; 
   muy00 = yx00[i,]*beta; 

   E11 = 1/(1 + exp((-1)*muy11))/(1 + exp((-1)*mum1)) + 1/(1 + exp((-1)*muy10))*(1-1/(1 + exp((-1)*mum1))); 
   E10 = 1/(1 + exp((-1)*muy11))/(1 + exp((-1)*mum0)) + 1/(1 + exp((-1)*muy10))*(1-1/(1 + exp((-1)*mum0))); 
   E01 = 1/(1 + exp((-1)*muy01))/(1 + exp((-1)*mum1)) + 1/(1 + exp((-1)*muy00))*(1-1/(1 + exp((-1)*mum1))); 
   E00 = 1/(1 + exp((-1)*muy01))/(1 + exp((-1)*mum0)) + 1/(1 + exp((-1)*muy00))*(1-1/(1 + exp((-1)*mum0))); 

   sum11 = sum11 + E11;
   sum10 = sum10 + E10;
   sum01 = sum01 + E01;
   sum00 = sum00 + E00;

   TE = TE + (E11 - E00);
   IE1 = IE1 + (E11 - E10);
   DE0 = DE0 + (E10 - E00);
   IE0 = IE0 + (E01 - E00);
   DE1 = DE1 + (E11 - E01);

   TEOR = TEOR +  LOG((E11/(1-E11))/(E00/(1-E00)));
   IE1OR = IE1OR + LOG((E11/(1-E11))/(E10/(1-E10)));
   DE0OR = DE0OR + LOG((E10/(1-E10))/(E00/(1-E00)));
   IE0OR = IE0OR + LOG((E01/(1-E01))/(E00/(1-E00)));
   DE1OR = DE1OR + LOG((E11/(1-E11))/(E01/(1-E01)));

/* Derivative for alpha coefficient */

   E11_1 = (1/(1 + exp((-1)*muy11)) - 1/(1 + exp((-1)*muy10)))*exp((-1) * mum1)/(1 + exp((-1) * mum1))**2; 
   E10_1 = (1/(1 + exp((-1)*muy11)) - 1/(1 + exp((-1)*muy10)))*exp((-1) * mum0)/(1 + exp((-1) * mum0))**2; 
   E01_1 = (1/(1 + exp((-1)*muy01)) - 1/(1 + exp((-1)*muy00)))*exp((-1) * mum1)/(1 + exp((-1) * mum1))**2; 
   E00_1 = (1/(1 + exp((-1)*muy01)) - 1/(1 + exp((-1)*muy00)))*exp((-1) * mum0)/(1 + exp((-1) * mum0))**2; 

/* Derivative for other beta except for beta2 */

   E11_3 = 1/(1 + exp((-1) * mum1)) * exp((-1)*muy11)/(1 + exp((-1)*muy11))**2 + (1 - 1/(1 + exp((-1) * mum1))) * exp((-1)*muy10)/(1 + exp((-1)*muy10))**2;
   E10_3 = 1/(1 + exp((-1) * mum0)) * exp((-1)*muy11)/(1 + exp((-1)*muy11))**2 + (1 - 1/(1 + exp((-1) * mum0))) * exp((-1)*muy10)/(1 + exp((-1)*muy10))**2; 
   E01_3 = 1/(1 + exp((-1) * mum1)) * exp((-1)*muy01)/(1 + exp((-1)*muy01))**2 + (1 - 1/(1 + exp((-1) * mum1))) * exp((-1)*muy00)/(1 + exp((-1)*muy00))**2; 
   E00_3 = 1/(1 + exp((-1) * mum0)) * exp((-1)*muy01)/(1 + exp((-1)*muy01))**2 + (1 - 1/(1 + exp((-1) * mum0))) * exp((-1)*muy00)/(1 + exp((-1)*muy00))**2; 

/* Derivative for beta2 */

   E11_4 = 1/(1 + exp((-1) * mum1)) * exp((-1)*muy11)/(1 + exp((-1)*muy11))**2;
   E10_4 = 1/(1 + exp((-1) * mum0)) * exp((-1)*muy11)/(1 + exp((-1)*muy11))**2; 
   E01_4 = 1/(1 + exp((-1) * mum1)) * exp((-1)*muy01)/(1 + exp((-1)*muy01))**2; 
   E00_4 = 1/(1 + exp((-1) * mum0)) * exp((-1)*muy01)/(1 + exp((-1)*muy01))**2; 

   gte[1:(k-1)] = gte[1:(k-1)] + E11_1 * t(mx1[i,]) - E00_1* t(mx0[i,]);
   gte[k:(k+1)] = gte[k:(k+1)] + E11_3 * t(yx11[i,1:2]) - E00_3 * t(yx00[i,1:2]);
   gte[k+2] = gte[k+2] + E11_4 - E00_4;
   gte[(k+3):2*k-1] = gte[(k+3):2*k-1] + E11_3 * t(yx11[i,4:k]) - E00_3 * t(yx00[i,4:k]);

   gie1[1:(k-1)] = gie1[1:(k-1)] + E11_1 * t(mx1[i,]) - E10_1* t(mx0[i,]);
   gie1[k:(k+1)] = gie1[k:(k+1)] + E11_3 * t(yx11[i,1:2]) - E10_3 * t(yx11[i,1:2]);
   gie1[k+2] = gie1[k+2] + E11_4 - E10_4;
   gie1[(k+3):2*k-1] = gie1[(k+3):2*k-1] + E11_3 * t(yx11[i,4:k]) - E10_3 * t(yx11[i,4:k]);

   gde0[1:(k-1)] = gde0[1:(k-1)] + E10_1 * t(mx0[i,]) - E00_1* t(mx0[i,]);
   gde0[k:(k+1)] = gde0[k:(k+1)] + E10_3 * t(yx11[i,1:2]) - E00_3 * t(yx00[i,1:2]);
   gde0[k+2] = gde0[k+2] + E10_4 - E00_4;
   gde0[(k+3):2*k-1] = gde0[(k+3):2*k-1] + E10_3 * t(yx11[i,4:k]) - E00_3 * t(yx00[i,4:k]);

   gie0[1:(k-1)] = gie0[1:(k-1)] + E01_1 * t(mx1[i,]) - E00_1* t(mx0[i,]);
   gie0[k:(k+1)] = gie0[k:(k+1)] + E01_3 * t(yx00[i,1:2]) - E00_3 * t(yx00[i,1:2]);
   gie0[k+2] = gie0[k+2] + E01_4 - E00_4;
   gie0[(k+3):2*k-1] = gie0[(k+3):2*k-1] + E01_3 * t(yx00[i,4:k]) - E00_3 * t(yx00[i,4:k]);

   gde1[1:(k-1)] = gde1[1:(k-1)] + E11_1 * t(mx1[i,]) - E01_1* t(mx1[i,]);
   gde1[k:(k+1)] = gde1[k:(k+1)] + E11_3 * t(yx11[i,1:2]) - E01_3 * t(yx00[i,1:2]);
   gde1[k+2] = gde1[k+2] + E11_4 - E01_4;
   gde1[(k+3):2*k-1] = gde1[(k+3):2*k-1] + E11_3 * t(yx11[i,4:k]) - E01_3 * t(yx00[i,4:k]);

   gteor[1:(k-1)] = gteor[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E00 - E11 * E00_1* t(mx0[i,]))*(E00 - E00 * E11)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E11 - E00 * E11_1 * t(mx1[i,]))*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[k:(k+1)] = gteor[k:(k+1)] + ((E11_3 * t(yx11[i,1:2]) - E11_3 * t(yx11[i,1:2]) * E00 - E11 * E00_3 * t(yx00[i,1:2]))*(E00 - E00 * E11)-(E00_3 * t(yx00[i,1:2]) - E00_3 * t(yx00[i,1:2]) * E11 - E00 * E11_3 * t(yx11[i,1:2]))*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[k+2] = gteor[k+2] + ((E11_4 - E11_4 * E00 - E11 * E00_4)*(E00 - E00 * E11)-(E00_4 - E00_4 * E11 - E00 * E11_4)*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[(k+3):2*k-1] = gteor[(k+3):2*k-1] + ((E11_3 * t(yx11[i,4:k]) - E11_3 * t(yx11[i,4:k]) * E00 - E11 * E00_3 * t(yx00[i,4:k]))*(E00 - E00 * E11)-(E00_3 * t(yx00[i,4:k]) - E00_3 * t(yx00[i,4:k]) * E11 - E00 * E11_3 * t(yx11[i,4:k]))*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;

   gie1or[1:(k-1)] = gie1or[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E10 - E11 * E10_1* t(mx0[i,]))*(E10 - E10 * E11)-(E10_1* t(mx0[i,]) - E10_1* t(mx0[i,]) * E11 - E10 * E11_1 * t(mx1[i,]))*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[k:(k+1)] = gie1or[k:(k+1)] + ((E11_3 * t(yx11[i,1:2]) - E11_3 * t(yx11[i,1:2]) * E10 - E11 * E10_3 * t(yx11[i,1:2]))*(E10 - E10 * E11)-(E10_3 * t(yx11[i,1:2]) - E10_3 * t(yx11[i,1:2]) * E11 - E10 * E11_3 * t(yx11[i,1:2]))*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[k+2] = gie1or[k+2] + ((E11_4 - E11_4 * E10 - E11 * E10_4)*(E10 - E10 * E11)-(E10_4 - E10_4 * E11 - E10 * E11_4)*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[(k+3):2*k-1] = gie1or[(k+3):2*k-1] + ((E11_3 * t(yx11[i,4:k]) - E11_3 * t(yx11[i,4:k]) * E10 - E11 * E10_3 * t(yx11[i,4:k]))*(E10 - E10 * E11)-(E10_3 * t(yx11[i,4:k]) - E10_3 * t(yx11[i,4:k]) * E11 - E10 * E11_3 * t(yx11[i,4:k]))*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;

   gde0or[1:(k-1)] = gde0or[1:(k-1)] + ((E10_1 * t(mx0[i,]) - E10_1 * t(mx0[i,]) * E00 - E10 * E00_1* t(mx0[i,]))*(E00 - E00 * E10)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E10 - E00 * E10_1 * t(mx0[i,]))*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[k:(k+1)] = gde0or[k:(k+1)] + ((E10_3 * t(yx11[i,1:2]) - E10_3 * t(yx11[i,1:2]) * E00 - E10 * E00_3 * t(yx00[i,1:2]))*(E00 - E00 * E10)-(E00_3 * t(yx00[i,1:2]) - E00_3 * t(yx00[i,1:2]) * E10 - E00 * E10_3 * t(yx11[i,1:2]))*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[k+2] = gde0or[k+2] + ((E10_4 - E10_4 * E00 - E10 * E00_4)*(E00 - E00 * E10)-(E00_4 - E00_4 * E10 - E00 * E10_4)*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[(k+3):2*k-1] = gde0or[(k+3):2*k-1] + ((E10_3 * t(yx11[i,4:k]) - E10_3 * t(yx11[i,4:k]) * E00 - E10 * E00_3 * t(yx00[i,4:k]))*(E00 - E00 * E10)-(E00_3 * t(yx00[i,4:k]) - E00_3 * t(yx00[i,4:k]) * E10 - E00 * E10_3 * t(yx11[i,4:k]))*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;

   gie0or[1:(k-1)] = gie0or[1:(k-1)] + ((E01_1 * t(mx1[i,]) - E01_1 * t(mx1[i,]) * E00 - E01 * E00_1* t(mx0[i,]))*(E00 - E00 * E01)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E01 - E00 * E01_1 * t(mx1[i,]))*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[k:(k+1)] = gie0or[k:(k+1)] + ((E01_3 * t(yx00[i,1:2]) - E01_3 * t(yx00[i,1:2]) * E00 - E01 * E00_3 * t(yx00[i,1:2]))*(E00 - E00 * E01)-(E00_3 * t(yx00[i,1:2]) - E00_3 * t(yx00[i,1:2]) * E01 - E00 * E01_3 * t(yx00[i,1:2]))*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[k+2] = gie0or[k+2] + ((E01_4 - E01_4 * E00 - E01 * E00_4)*(E00 - E00 * E01)-(E00_4 - E00_4 * E01 - E00 * E01_4)*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[(k+3):2*k-1] = gie0or[(k+3):2*k-1] + ((E01_3 * t(yx00[i,4:k]) - E01_3 * t(yx00[i,4:k]) * E00 - E01 * E00_3 * t(yx00[i,4:k]))*(E00 - E00 * E01)-(E00_3 * t(yx00[i,4:k]) - E00_3 * t(yx00[i,4:k]) * E01 - E00 * E01_3 * t(yx00[i,4:k]))*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;

   gde1or[1:(k-1)] = gde1or[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E01 - E11 * E01_1* t(mx1[i,]))*(E01 - E01 * E11)-(E01_1* t(mx1[i,]) - E01_1* t(mx1[i,]) * E11 - E01 * E11_1 * t(mx1[i,]))*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[k:(k+1)] = gde1or[k:(k+1)] + ((E11_3 * t(yx11[i,1:2]) - E11_3 * t(yx11[i,1:2]) * E01 - E11 * E01_3 * t(yx00[i,1:2]))*(E01 - E01 * E11)-(E01_3 * t(yx00[i,1:2]) - E01_3 * t(yx00[i,1:2]) * E11 - E01 * E11_3 * t(yx11[i,1:2]))*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[k+2] = gde1or[k+2] + ((E11_4 - E11_4 * E01 - E11 * E01_4)*(E01 - E01 * E11)-(E01_4 - E01_4 * E11 - E01 * E11_4)*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[(k+3):2*k-1] = gde1or[(k+3):2*k-1] + ((E11_3 * t(yx11[i,4:k]) - E11_3 * t(yx11[i,4:k]) * E01 - E11 * E01_3 * t(yx00[i,4:k]))*(E01 - E01 * E11)-(E01_3 * t(yx00[i,4:k]) - E01_3 * t(yx00[i,4:k]) * E11 - E01 * E11_3 * t(yx11[i,4:k]))*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;

end;

gr1 = (gie1 * te - gte * ie1)/te ** 2;
gr1var = gr1 * var  * t(gr1);
if gr1var >= 0 then gr1sd = sqrt(gr1var);

gr0 = (gie0 * te - gte * ie0)/te ** 2;
gr0var = gr0 * var * t(gr0);
if gr0var >= 0 then gr0sd = sqrt(gr0var);

gr1or = (gie1or * teor - gteor * ie1or)/teor ** 2;
gr1orvar = gr1or * var  * t(gr1or);
if gr1orvar >= 0 then gr1orsd = sqrt(gr1orvar);

gr0or = (gie0or * teor - gteor * ie0or)/teor ** 2;
gr0orvar = gr0or * var  * t(gr0or);
if gr0orvar >= 0 then gr0orsd = sqrt(gr0orvar);

sum11 = sum11/n;
sum10 = sum10/n;
sum01 = sum01/n;
sum00 = sum00/n;

TE = TE/n;
IE1 = IE1/n;
DE0 = DE0/n;
IE0 = IE0/n;
DE1 = DE1/n;

R1 = IE1/TE;
R0 = IE0/TE;

TEOR = TEOR * (1/n);
IE1OR = IE1OR * (1/n);
DE0OR = DE0OR * (1/n);
IE0OR = IE0OR * (1/n);
DE1OR = DE1OR * (1/n);

R1OR = IE1OR/TEOR;
R0OR = IE0OR/TEOR;

gte = gte/n;
gtevar = gte * var * t(gte);
if gtevar >= 0 then gtesd = sqrt(gtevar);
gie1 = gie1/n;
gie1var = gie1 * var * t(gie1);
if gie1var >= 0 then gie1sd = sqrt(gie1var);
gde0 = gde0/n;
gde0var = gde0 * var * t(gde0);
if gde0var >= 0 then gde0sd = sqrt(gde0var);
gie0 = gie0/n;
gie0var = gie0 * var * t(gie0);
if gie0var >= 0 then gie0sd = sqrt(gie0var);
gde1 = gde1/n;
gde1var = gde1 * var * t(gde1);
if gde1var >= 0 then gde1sd = sqrt(gde1var);

gteor = gteor/n;
gteorvar = gteor * var * t(gteor);
if gteorvar >= 0 then gteorsd = sqrt(gteorvar);
gie1or = gie1or/n;
gie1orvar = gie1or * var * t(gie1or);
if gie1orvar >= 0 then gie1orsd = sqrt(gie1orvar);
gde0or = gde0or/n;
gde0orvar = gde0or * var * t(gde0or);
if gde0orvar >= 0 then gde0orsd = sqrt(gde0orvar);
gie0or = gie0or/n;
gie0orvar = gie0or * var * t(gie0or);
if gie0orvar >= 0 then gie0orsd = sqrt(gie0orvar);
gde1or = gde1or/n;
gde1orvar = gde1or * var * t(gde1or);
if gde1orvar >= 0 then gde1orsd = sqrt(gde1orvar);

print sum11 sum10 sum01 sum00 TE IE1 DE0 IE0 DE1 TEOR IE1OR DE0OR IE0OR DE1OR gtesd gie1sd gde0sd gie0sd gde1sd
gteorsd gie1orsd gde0orsd gie0orsd gde1orsd;

postprobs=sum11||sum10||sum01||sum00||TE||IE1||DE0||IE0||DE1||TEOR||IE1OR||DE0OR||IE0OR||DE1OR||gtesd||gie1sd||gde0sd||gie0sd||gde1sd
||gteorsd||gie1orsd||gde0orsd||gie0orsd||gde1orsd||r1||r0||r1or||r0or||gr1sd||gr0sd||gr1orsd||gr0orsd;
 cname = {"sum11" "sum10" "sum01" "sum00" "TE" "IE1" "DE0" "IE0" "DE1" "LOGTEOR" "LOGIE1OR" "LOGDE0OR" "LOGIE0OR" "LOGDE1OR" "gtesd" "gie1sd" "gde0sd" "gie0sd" 
"gde1sd" "gteorsd" "gie1orsd" "gde0orsd" "gie0orsd" "gde1orsd" "r1" "r0" "r1or" "r0or" "gr1sd" "gr0sd" "gr1orsd" "gr0orsd"};

create &out from postprobs  [ colname=cname ];
append from postprobs;

quit;

%mend;

/* Continous exposure and continuous mediator without any covariate */

%macro conmedcon1(dataset=' ', y=' ', x=' ', m = '', out =' ');

proc logistic data = &dataset;
model &y (event = '1') = &x &m/COVB;
ods output covb = aa1
ParameterEstimates = aa2;
run; 

proc transpose data = aa2 out = aa3;
id variable;
var estimate;
run;

data aa4;
merge ff11 aa3;
drop _name_;
run;

data ff13;
set ff12;
id = _n_;
run;

proc means data = ff12;
var col1;
output out = aa11 n(col1) = np;
run;

data aa12;
set aa11;
nb = 2 * np;
run;

data _null_;
set aa11;
call symput('NP', np);
call symput('NB', nb);
run;

data aa21;
set aa1;
id = _n_ + &NP;
run;

data aa5;
merge ff13 aa21;
by id;
drop id parameter;
run;

data aa6;
set aa5;
array x[*] col1 -- &m;
do i = 1 to dim(x);
if x[i] = . then x[i] = 0;
end;
drop i;
run;

proc iml;
reset noname;
use &dataset;
read all var {&x} into x;
use aa6;
read all into var;
use abswei;
read all var {abs weight} into w;
use aa4;
read all into theta;

n = nrow(x);
k = ncol(theta)/2;
alpha = theta [1:k-1];
sigma = theta[k];
beta = theta[k+1:2*k];
mx0 = j(n,1,1) || x;
mx1 = j(n,1,1) || (x[,1] + 1) ;

yx0 = j(n,1,1) || x[,1] || j(n,1,0);
yx1 = j(n,1,1) || (x[,1] + 1) || j(n,1,0);
/*x2 = yx0[1:3,];*/
/*x3 = yx1[1:3,];*/
/**/
/*print theta k alpha beta x2 x3;*/
sum11 = 0;
sum10 = 0;
sum01 = 0;
sum00 = 0;

TE = 0;
IE1 = 0;
DE0 = 0;
IE0 = 0;
DE1 = 0;

TEOR = 0;
IE1OR = 0;
DE0OR = 0;
IE0OR = 0;
DE1OR = 0;

gte = j(1,k*2,0);
gie1 = j(1,k*2,0);
gde0 = j(1,k*2,0);
gie0 = j(1,k*2,0);
gde1 = j(1,k*2,0);

gteor = j(1,k*2,0);
gie1or = j(1,k*2,0);
gde0or = j(1,k*2,0);
gie0or = j(1,k*2,0);
gde1or = j(1,k*2,0);

do i = 1 to n;
   mum1 = mx1[i,]*alpha; 
   mum0 = mx0[i,]*alpha; 
   muy11 = yx1[i,]*beta + beta[3] * (sqrt(2) * sigma * w[,1] + mum1); 
   muy10 = yx1[i,]*beta + beta[3] * (sqrt(2) * sigma * w[,1] + mum0); 
   muy01 = yx0[i,]*beta + beta[3] * (sqrt(2) * sigma * w[,1] + mum1); 
   muy00 = yx0[i,]*beta + beta[3] * (sqrt(2) * sigma * w[,1] + mum0); 

   E11 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy11))))* w[,2]; 
   E10 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy10))))* w[,2]; 
   E01 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy01))))* w[,2]; 
   E00 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy00))))* w[,2]; 

   sum11 = sum11 + E11;
   sum10 = sum10 + E10;
   sum01 = sum01 + E01;
   sum00 = sum00 + E00;

   TE = TE + (E11 - E00);
   IE1 = IE1 + (E11 - E10);
   DE0 = DE0 + (E10 - E00);
   IE0 = IE0 + (E01 - E00);
   DE1 = DE1 + (E11 - E01);

   TEOR = TEOR +  LOG((E11/(1-E11))/(E00/(1-E00)));
   IE1OR = IE1OR + LOG((E11/(1-E11))/(E10/(1-E10)));
   DE0OR = DE0OR + LOG((E10/(1-E10))/(E00/(1-E00)));
   IE0OR = IE0OR + LOG((E01/(1-E01))/(E00/(1-E00)));
   DE1OR = DE1OR + LOG((E11/(1-E11))/(E01/(1-E01)));

/* Derivative for alpha coefficient */

   E11_1 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy11))) # (sqrt(2) * w[,1]/sigma))* w[,2]; 
   E10_1 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy10))) # (sqrt(2) * w[,1]/sigma))* w[,2]; 
   E01_1 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy01))) # (sqrt(2) * w[,1]/sigma))* w[,2]; 
   E00_1 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy00))) # (sqrt(2) * w[,1]/sigma))* w[,2]; 

/* Derivative for sigma */

   E11_2 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy11))) # (2 * w[,1] # w[,1]/sigma - 1/sigma))* w[,2]; 
   E10_2 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy10))) # (2 * w[,1] # w[,1]/sigma - 1/sigma))* w[,2]; 
   E01_2 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy01))) # (2 * w[,1] # w[,1]/sigma - 1/sigma))* w[,2]; 
   E00_2 = t(1/sqrt(constant('pi')) * (1/(1 + exp((-1) * muy00))) # (2 * w[,1] # w[,1]/sigma - 1/sigma))* w[,2]; 

/* Derivative for other beta except for beta2 */

   E11_3 = t(1/sqrt(constant('pi')) * (exp((-1) * muy11)/(1 + exp((-1) * muy11))/(1 + exp((-1) * muy11)) ))* w[,2]; 
   E10_3 = t(1/sqrt(constant('pi')) * (exp((-1) * muy10)/(1 + exp((-1) * muy10))/(1 + exp((-1) * muy10)) ))* w[,2]; 
   E01_3 = t(1/sqrt(constant('pi')) * (exp((-1) * muy01)/(1 + exp((-1) * muy01))/(1 + exp((-1) * muy01)) ))* w[,2]; 
   E00_3 = t(1/sqrt(constant('pi')) * (exp((-1) * muy00)/(1 + exp((-1) * muy00))/(1 + exp((-1) * muy00)) ))* w[,2]; 

/* Derivative for beta2 */

   E11_4 = t(1/sqrt(constant('pi')) * (exp((-1) * muy11)/(1 + exp((-1) * muy11))/(1 + exp((-1) * muy11)) # (sqrt(2) * w[,1] * sigma + mum1) ))* w[,2]; 
   E10_4 = t(1/sqrt(constant('pi')) * (exp((-1) * muy10)/(1 + exp((-1) * muy10))/(1 + exp((-1) * muy10)) # (sqrt(2) * w[,1] * sigma + mum0) ))* w[,2]; 
   E01_4 = t(1/sqrt(constant('pi')) * (exp((-1) * muy01)/(1 + exp((-1) * muy01))/(1 + exp((-1) * muy01)) # (sqrt(2) * w[,1] * sigma + mum1) ))* w[,2]; 
   E00_4 = t(1/sqrt(constant('pi')) * (exp((-1) * muy00)/(1 + exp((-1) * muy00))/(1 + exp((-1) * muy00)) # (sqrt(2) * w[,1] * sigma + mum0) ))* w[,2]; 

   gte[1:(k-1)] = gte[1:(k-1)] + E11_1 * t(mx1[i,]) - E00_1* t(mx0[i,]);
   gte[k] = gte[k] + E11_2 - E00_2;
   gte[(k+1):(k+2)] = gte[(k+1):(k+2)] + E11_3 * t(yx1[i,1:2]) - E00_3 * t(yx0[i,1:2]);
   gte[k+3] = gte[k+3] + E11_4 - E00_4;

   gie1[1:(k-1)] = gie1[1:(k-1)] + E11_1 * t(mx1[i,]) - E10_1* t(mx0[i,]);
   gie1[k] = gie1[k] + E11_2 - E10_2;
   gie1[(k+1):(k+2)] = gie1[(k+1):(k+2)] + E11_3 * t(yx1[i,1:2]) - E10_3 * t(yx1[i,1:2]);
   gie1[k+3] = gie1[k+3] + E11_4 - E10_4;

   gde0[1:(k-1)] = gde0[1:(k-1)] + E10_1 * t(mx0[i,]) - E00_1* t(mx0[i,]);
   gde0[k] = gde0[k] + E10_2 - E00_2;
   gde0[(k+1):(k+2)] = gde0[(k+1):(k+2)] + E10_3 * t(yx1[i,1:2]) - E00_3 * t(yx0[i,1:2]);
   gde0[k+3] = gde0[k+3] + E10_4 - E00_4;

   gie0[1:(k-1)] = gie0[1:(k-1)] + E01_1 * t(mx1[i,]) - E00_1* t(mx0[i,]);
   gie0[k] = gie0[k] + E01_2 - E00_2;
   gie0[(k+1):(k+2)] = gie0[(k+1):(k+2)] + E01_3 * t(yx0[i,1:2]) - E00_3 * t(yx0[i,1:2]);
   gie0[k+3] = gie0[k+3] + E01_4 - E00_4;

   gde1[1:(k-1)] = gde1[1:(k-1)] + E11_1 * t(mx1[i,]) - E01_1* t(mx1[i,]);
   gde1[k] = gde1[k] + E11_2 - E01_2;
   gde1[(k+1):(k+2)] = gde1[(k+1):(k+2)] + E11_3 * t(yx1[i,1:2]) - E01_3 * t(yx0[i,1:2]);
   gde1[k+3] = gde1[k+3] + E11_4 - E01_4;

   gteor[1:(k-1)] = gteor[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E00 - E11 * E00_1* t(mx0[i,]))*(E00 - E00 * E11)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E11 - E00 * E11_1 * t(mx1[i,]))*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[k] = gteor[k] + ((E11_2 - E11_2 * E00 - E11 * E00_2)*(E00 - E00 * E11)-(E00_2 - E00_2 * E11 - E00 * E11_2)*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[(k+1):(k+2)] = gteor[(k+1):(k+2)] + ((E11_3 * t(yx1[i,1:2]) - E11_3 * t(yx1[i,1:2]) * E00 - E11 * E00_3 * t(yx0[i,1:2]))*(E00 - E00 * E11)-(E00_3 * t(yx0[i,1:2]) - E00_3 * t(yx0[i,1:2]) * E11 - E00 * E11_3 * t(yx1[i,1:2]))*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[k+3] = gteor[k+3] + ((E11_4 - E11_4 * E00 - E11 * E00_4)*(E00 - E00 * E11)-(E00_4 - E00_4 * E11 - E00 * E11_4)*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;

   gie1or[1:(k-1)] = gie1or[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E10 - E11 * E10_1* t(mx0[i,]))*(E10 - E10 * E11)-(E10_1* t(mx0[i,]) - E10_1* t(mx0[i,]) * E11 - E10 * E11_1 * t(mx1[i,]))*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[k] = gie1or[k] + ((E11_2 - E11_2 * E10 - E11 * E10_2)*(E10 - E10 * E11)-(E10_2 - E10_2 * E11 - E10 * E11_2)*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[(k+1):(k+2)] = gie1or[(k+1):(k+2)] + ((E11_3 * t(yx1[i,1:2]) - E11_3 * t(yx1[i,1:2]) * E10 - E11 * E10_3 * t(yx1[i,1:2]))*(E10 - E10 * E11)-(E10_3 * t(yx1[i,1:2]) - E10_3 * t(yx1[i,1:2]) * E11 - E10 * E11_3 * t(yx1[i,1:2]))*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[k+3] = gie1or[k+3] + ((E11_4 - E11_4 * E10 - E11 * E10_4)*(E10 - E10 * E11)-(E10_4 - E10_4 * E11 - E10 * E11_4)*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;

   gde0or[1:(k-1)] = gde0or[1:(k-1)] + ((E10_1 * t(mx0[i,]) - E10_1 * t(mx0[i,]) * E00 - E10 * E00_1* t(mx0[i,]))*(E00 - E00 * E10)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E10 - E00 * E10_1 * t(mx0[i,]))*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[k] = gde0or[k] + ((E10_2 - E10_2 * E00 - E10 * E00_2)*(E00 - E00 * E10)-(E00_2 - E00_2 * E10 - E00 * E10_2)*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[(k+1):(k+2)] = gde0or[(k+1):(k+2)] + ((E10_3 * t(yx1[i,1:2]) - E10_3 * t(yx1[i,1:2]) * E00 - E10 * E00_3 * t(yx0[i,1:2]))*(E00 - E00 * E10)-(E00_3 * t(yx0[i,1:2]) - E00_3 * t(yx0[i,1:2]) * E10 - E00 * E10_3 * t(yx1[i,1:2]))*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[k+3] = gde0or[k+3] + ((E10_4 - E10_4 * E00 - E10 * E00_4)*(E00 - E00 * E10)-(E00_4 - E00_4 * E10 - E00 * E10_4)*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;

   gie0or[1:(k-1)] = gie0or[1:(k-1)] + ((E01_1 * t(mx1[i,]) - E01_1 * t(mx1[i,]) * E00 - E01 * E00_1* t(mx0[i,]))*(E00 - E00 * E01)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E01 - E00 * E01_1 * t(mx1[i,]))*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[k] = gie0or[k] + ((E01_2 - E01_2 * E00 - E01 * E00_2)*(E00 - E00 * E01)-(E00_2 - E00_2 * E01 - E00 * E01_2)*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[(k+1):(k+2)] = gie0or[(k+1):(k+2)] + ((E01_3 * t(yx0[i,1:2]) - E01_3 * t(yx0[i,1:2]) * E00 - E01 * E00_3 * t(yx0[i,1:2]))*(E00 - E00 * E01)-(E00_3 * t(yx0[i,1:2]) - E00_3 * t(yx0[i,1:2]) * E01 - E00 * E01_3 * t(yx0[i,1:2]))*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[k+3] = gie0or[k+3] + ((E01_4 - E01_4 * E00 - E01 * E00_4)*(E00 - E00 * E01)-(E00_4 - E00_4 * E01 - E00 * E01_4)*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;

   gde1or[1:(k-1)] = gde1or[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E01 - E11 * E01_1* t(mx1[i,]))*(E01 - E01 * E11)-(E01_1* t(mx1[i,]) - E01_1* t(mx1[i,]) * E11 - E01 * E11_1 * t(mx1[i,]))*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[k] = gde1or[k] + ((E11_2 - E11_2 * E01 - E11 * E01_2)*(E01 - E01 * E11)-(E01_2 - E01_2 * E11 - E01 * E11_2)*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[(k+1):(k+2)] = gde1or[(k+1):(k+2)] + ((E11_3 * t(yx1[i,1:2]) - E11_3 * t(yx1[i,1:2]) * E01 - E11 * E01_3 * t(yx0[i,1:2]))*(E01 - E01 * E11)-(E01_3 * t(yx0[i,1:2]) - E01_3 * t(yx0[i,1:2]) * E11 - E01 * E11_3 * t(yx1[i,1:2]))*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[k+3] = gde1or[k+3] + ((E11_4 - E11_4 * E01 - E11 * E01_4)*(E01 - E01 * E11)-(E01_4 - E01_4 * E11 - E01 * E11_4)*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;

end;

gr1 = (gie1 * te - gte * ie1)/te ** 2;
gr1var = gr1 * var  * t(gr1);
if gr1var >= 0 then gr1sd = sqrt(gr1var);

gr0 = (gie0 * te - gte * ie0)/te ** 2;
gr0var = gr0 * var * t(gr0);
if gr0var >= 0 then gr0sd = sqrt(gr0var);

gr1or = (gie1or * teor - gteor * ie1or)/teor ** 2;
gr1orvar = gr1or * var  * t(gr1or);
if gr1orvar >= 0 then gr1orsd = sqrt(gr1orvar);

gr0or = (gie0or * teor - gteor * ie0or)/teor ** 2;
gr0orvar = gr0or * var  * t(gr0or);
if gr0orvar >= 0 then gr0orsd = sqrt(gr0orvar);

sum11 = sum11/n;
sum10 = sum10/n;
sum01 = sum01/n;
sum00 = sum00/n;

TE = TE/n;
IE1 = IE1/n;
DE0 = DE0/n;
IE0 = IE0/n;
DE1 = DE1/n;

R1 = IE1/TE;
R0 = IE0/TE;

TEOR = TEOR * (1/n);
IE1OR = IE1OR * (1/n);
DE0OR = DE0OR * (1/n);
IE0OR = IE0OR * (1/n);
DE1OR = DE1OR * (1/n);

R1OR = IE1OR/TEOR;
R0OR = IE0OR/TEOR;

gte = gte/n;
gtevar = gte * var * t(gte);
if gtevar >= 0 then gtesd = sqrt(gtevar);
gie1 = gie1/n;
gie1var = gie1 * var * t(gie1);
if gie1var >= 0 then gie1sd = sqrt(gie1var);
gde0 = gde0/n;
gde0var = gde0 * var * t(gde0);
if gde0var >= 0 then gde0sd = sqrt(gde0var);
gie0 = gie0/n;
gie0var = gie0 * var * t(gie0);
if gie0var >= 0 then gie0sd = sqrt(gie0var);
gde1 = gde1/n;
gde1var = gde1 * var * t(gde1);
if gde1var >= 0 then gde1sd = sqrt(gde1var);

gteor = gteor/n;
gteorvar = gteor * var * t(gteor);
if gteorvar >= 0 then gteorsd = sqrt(gteorvar);
gie1or = gie1or/n;
gie1orvar = gie1or * var * t(gie1or);
if gie1orvar >= 0 then gie1orsd = sqrt(gie1orvar);
gde0or = gde0or/n;
gde0orvar = gde0or * var * t(gde0or);
if gde0orvar >= 0 then gde0orsd = sqrt(gde0orvar);
gie0or = gie0or/n;
gie0orvar = gie0or * var * t(gie0or);
if gie0orvar >= 0 then gie0orsd = sqrt(gie0orvar);
gde1or = gde1or/n;
gde1orvar = gde1or * var * t(gde1or);
if gde1orvar >= 0 then gde1orsd = sqrt(gde1orvar);

print sum11 sum10 sum01 sum00 TE IE1 DE0 IE0 DE1 TEOR IE1OR DE0OR IE0OR DE1OR gtesd gie1sd gde0sd gie0sd gde1sd
gteorsd gie1orsd gde0orsd gie0orsd gde1orsd;

postprobs=sum11||sum10||sum01||sum00||TE||IE1||DE0||IE0||DE1||TEOR||IE1OR||DE0OR||IE0OR||DE1OR||gtesd||gie1sd||gde0sd||gie0sd||gde1sd
||gteorsd||gie1orsd||gde0orsd||gie0orsd||gde1orsd||r1||r0||r1or||r0or||gr1sd||gr0sd||gr1orsd||gr0orsd;
 cname = {"sum11" "sum10" "sum01" "sum00" "TE" "IE1" "DE0" "IE0" "DE1" "LOGTEOR" "LOGIE1OR" "LOGDE0OR" "LOGIE0OR" "LOGDE1OR" "gtesd" "gie1sd" "gde0sd" "gie0sd" 
"gde1sd" "gteorsd" "gie1orsd" "gde0orsd" "gie0orsd" "gde1orsd" "r1" "r0" "r1or" "r0or" "gr1sd" "gr0sd" "gr1orsd" "gr0orsd"};

create &out from postprobs  [ colname=cname ];
append from postprobs;

quit;

%mend;

/* Continous exposure and binary mediator without any covariate */

%macro conmedbi1(dataset=' ', y=' ', x=' ', m = '', out =' ');

proc logistic data = &dataset;
model &m (event = '1') = &x/COVB;
ods output covb = bb1
ParameterEstimates = bb2;
run; 

proc transpose data = bb2 out = bb3;
var estimate;
run;

proc logistic data = &dataset;
model &y (event = '1') = &x &m/COVB;
ods output covb = aa1
ParameterEstimates = aa2;
run; 

proc transpose data = aa2 out = aa3;
id variable;
var estimate;
run;

data aa4;
merge bb3 aa3;
drop _name_;
run;

proc transpose data = bb1 out = bb6;
run;

data bb4;
set bb6;
id = _n_;
run;

proc means data = bb4;
var id;
output out = bb5 n(id) = np;
run;

data _null_;
set bb5;
call symput('NP', np);
run;

data aa21;
set aa1;
id = _n_ + &NP;
run;

data aa5;
merge bb4 aa21;
by id;
drop id parameter _name_;
run;

data aa6;
set aa5;
array x[*] col1 -- &m;
do i = 1 to dim(x);
if x[i] = . then x[i] = 0;
end;
drop i;
run;

proc iml;
reset noname;
use &dataset;
read all var {&x} into x;
use aa6;
read all into var;
use abswei;
read all var {abs weight} into w;
use aa4;
read all into theta;

n = nrow(x);
k = (ncol(theta) + 1)/2;
alpha = theta [1:k-1];
beta = theta[k:2*k - 1];
mx0 = j(n,1,1) || x;
mx1 = j(n,1,1) || (x[,1] + 1);

yx00 = j(n,1,1) || x[,1] || j(n,1,0);
yx10 = j(n,1,1) || (x[,1] + 1) || j(n,1,0);
yx01 = j(n,1,1) || x[,1] || j(n,1,1);
yx11 = j(n,1,1) || (x[,1] + 1) || j(n,1,1);
/*x2 = yx0[1:3,];*/
/*x3 = yx1[1:3,];*/
/**/
/*print theta k alpha beta x2 x3;*/
sum11 = 0;
sum10 = 0;
sum01 = 0;
sum00 = 0;

TE = 0;
IE1 = 0;
DE0 = 0;
IE0 = 0;
DE1 = 0;

TEOR = 0;
IE1OR = 0;
DE0OR = 0;
IE0OR = 0;
DE1OR = 0;

gte = j(1,k*2-1,0);
gie1 = j(1,k*2-1,0);
gde0 = j(1,k*2-1,0);
gie0 = j(1,k*2-1,0);
gde1 = j(1,k*2-1,0);

gteor = j(1,k*2-1,0);
gie1or = j(1,k*2-1,0);
gde0or = j(1,k*2-1,0);
gie0or = j(1,k*2-1,0);
gde1or = j(1,k*2-1,0);

do i = 1 to n;
   mum1 = mx1[i,]*alpha; 
   mum0 = mx0[i,]*alpha; 
   muy11 = yx11[i,]*beta; 
   muy10 = yx10[i,]*beta; 
   muy01 = yx01[i,]*beta; 
   muy00 = yx00[i,]*beta; 

   E11 = 1/(1 + exp((-1)*muy11))/(1 + exp((-1)*mum1)) + 1/(1 + exp((-1)*muy10))*(1-1/(1 + exp((-1)*mum1))); 
   E10 = 1/(1 + exp((-1)*muy11))/(1 + exp((-1)*mum0)) + 1/(1 + exp((-1)*muy10))*(1-1/(1 + exp((-1)*mum0))); 
   E01 = 1/(1 + exp((-1)*muy01))/(1 + exp((-1)*mum1)) + 1/(1 + exp((-1)*muy00))*(1-1/(1 + exp((-1)*mum1))); 
   E00 = 1/(1 + exp((-1)*muy01))/(1 + exp((-1)*mum0)) + 1/(1 + exp((-1)*muy00))*(1-1/(1 + exp((-1)*mum0))); 

   sum11 = sum11 + E11;
   sum10 = sum10 + E10;
   sum01 = sum01 + E01;
   sum00 = sum00 + E00;

   TE = TE + (E11 - E00);
   IE1 = IE1 + (E11 - E10);
   DE0 = DE0 + (E10 - E00);
   IE0 = IE0 + (E01 - E00);
   DE1 = DE1 + (E11 - E01);

   TEOR = TEOR +  LOG((E11/(1-E11))/(E00/(1-E00)));
   IE1OR = IE1OR + LOG((E11/(1-E11))/(E10/(1-E10)));
   DE0OR = DE0OR + LOG((E10/(1-E10))/(E00/(1-E00)));
   IE0OR = IE0OR + LOG((E01/(1-E01))/(E00/(1-E00)));
   DE1OR = DE1OR + LOG((E11/(1-E11))/(E01/(1-E01)));

/* Derivative for alpha coefficient */

   E11_1 = (1/(1 + exp((-1)*muy11)) - 1/(1 + exp((-1)*muy10)))*exp((-1) * mum1)/(1 + exp((-1) * mum1))**2; 
   E10_1 = (1/(1 + exp((-1)*muy11)) - 1/(1 + exp((-1)*muy10)))*exp((-1) * mum0)/(1 + exp((-1) * mum0))**2; 
   E01_1 = (1/(1 + exp((-1)*muy01)) - 1/(1 + exp((-1)*muy00)))*exp((-1) * mum1)/(1 + exp((-1) * mum1))**2; 
   E00_1 = (1/(1 + exp((-1)*muy01)) - 1/(1 + exp((-1)*muy00)))*exp((-1) * mum0)/(1 + exp((-1) * mum0))**2; 

/* Derivative for other beta except for beta2 */

   E11_3 = 1/(1 + exp((-1) * mum1)) * exp((-1)*muy11)/(1 + exp((-1)*muy11))**2 + (1 - 1/(1 + exp((-1) * mum1))) * exp((-1)*muy10)/(1 + exp((-1)*muy10))**2;
   E10_3 = 1/(1 + exp((-1) * mum0)) * exp((-1)*muy11)/(1 + exp((-1)*muy11))**2 + (1 - 1/(1 + exp((-1) * mum0))) * exp((-1)*muy10)/(1 + exp((-1)*muy10))**2; 
   E01_3 = 1/(1 + exp((-1) * mum1)) * exp((-1)*muy01)/(1 + exp((-1)*muy01))**2 + (1 - 1/(1 + exp((-1) * mum1))) * exp((-1)*muy00)/(1 + exp((-1)*muy00))**2; 
   E00_3 = 1/(1 + exp((-1) * mum0)) * exp((-1)*muy01)/(1 + exp((-1)*muy01))**2 + (1 - 1/(1 + exp((-1) * mum0))) * exp((-1)*muy00)/(1 + exp((-1)*muy00))**2; 

/* Derivative for beta2 */

   E11_4 = 1/(1 + exp((-1) * mum1)) * exp((-1)*muy11)/(1 + exp((-1)*muy11))**2;
   E10_4 = 1/(1 + exp((-1) * mum0)) * exp((-1)*muy11)/(1 + exp((-1)*muy11))**2; 
   E01_4 = 1/(1 + exp((-1) * mum1)) * exp((-1)*muy01)/(1 + exp((-1)*muy01))**2; 
   E00_4 = 1/(1 + exp((-1) * mum0)) * exp((-1)*muy01)/(1 + exp((-1)*muy01))**2; 

   gte[1:(k-1)] = gte[1:(k-1)] + E11_1 * t(mx1[i,]) - E00_1* t(mx0[i,]);
   gte[k:(k+1)] = gte[k:(k+1)] + E11_3 * t(yx11[i,1:2]) - E00_3 * t(yx00[i,1:2]);
   gte[k+2] = gte[k+2] + E11_4 - E00_4;

   gie1[1:(k-1)] = gie1[1:(k-1)] + E11_1 * t(mx1[i,]) - E10_1* t(mx0[i,]);
   gie1[k:(k+1)] = gie1[k:(k+1)] + E11_3 * t(yx11[i,1:2]) - E10_3 * t(yx11[i,1:2]);
   gie1[k+2] = gie1[k+2] + E11_4 - E10_4;

   gde0[1:(k-1)] = gde0[1:(k-1)] + E10_1 * t(mx0[i,]) - E00_1* t(mx0[i,]);
   gde0[k:(k+1)] = gde0[k:(k+1)] + E10_3 * t(yx11[i,1:2]) - E00_3 * t(yx00[i,1:2]);
   gde0[k+2] = gde0[k+2] + E10_4 - E00_4;

   gie0[1:(k-1)] = gie0[1:(k-1)] + E01_1 * t(mx1[i,]) - E00_1* t(mx0[i,]);
   gie0[k:(k+1)] = gie0[k:(k+1)] + E01_3 * t(yx00[i,1:2]) - E00_3 * t(yx00[i,1:2]);
   gie0[k+2] = gie0[k+2] + E01_4 - E00_4;

   gde1[1:(k-1)] = gde1[1:(k-1)] + E11_1 * t(mx1[i,]) - E01_1* t(mx1[i,]);
   gde1[k:(k+1)] = gde1[k:(k+1)] + E11_3 * t(yx11[i,1:2]) - E01_3 * t(yx00[i,1:2]);
   gde1[k+2] = gde1[k+2] + E11_4 - E01_4;

   gteor[1:(k-1)] = gteor[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E00 - E11 * E00_1* t(mx0[i,]))*(E00 - E00 * E11)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E11 - E00 * E11_1 * t(mx1[i,]))*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[k:(k+1)] = gteor[k:(k+1)] + ((E11_3 * t(yx11[i,1:2]) - E11_3 * t(yx11[i,1:2]) * E00 - E11 * E00_3 * t(yx00[i,1:2]))*(E00 - E00 * E11)-(E00_3 * t(yx00[i,1:2]) - E00_3 * t(yx00[i,1:2]) * E11 - E00 * E11_3 * t(yx11[i,1:2]))*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;
   gteor[k+2] = gteor[k+2] + ((E11_4 - E11_4 * E00 - E11 * E00_4)*(E00 - E00 * E11)-(E00_4 - E00_4 * E11 - E00 * E11_4)*(E11 - E00 * E11))/((E11/(1-E11))/(E00/(1-E00)))/(E00*(1-E11))**2;

   gie1or[1:(k-1)] = gie1or[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E10 - E11 * E10_1* t(mx0[i,]))*(E10 - E10 * E11)-(E10_1* t(mx0[i,]) - E10_1* t(mx0[i,]) * E11 - E10 * E11_1 * t(mx1[i,]))*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[k:(k+1)] = gie1or[k:(k+1)] + ((E11_3 * t(yx11[i,1:2]) - E11_3 * t(yx11[i,1:2]) * E10 - E11 * E10_3 * t(yx11[i,1:2]))*(E10 - E10 * E11)-(E10_3 * t(yx11[i,1:2]) - E10_3 * t(yx11[i,1:2]) * E11 - E10 * E11_3 * t(yx11[i,1:2]))*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;
   gie1or[k+2] = gie1or[k+2] + ((E11_4 - E11_4 * E10 - E11 * E10_4)*(E10 - E10 * E11)-(E10_4 - E10_4 * E11 - E10 * E11_4)*(E11 - E10 * E11))/((E11/(1-E11))/(E10/(1-E10)))/(E10*(1-E11))**2;

   gde0or[1:(k-1)] = gde0or[1:(k-1)] + ((E10_1 * t(mx0[i,]) - E10_1 * t(mx0[i,]) * E00 - E10 * E00_1* t(mx0[i,]))*(E00 - E00 * E10)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E10 - E00 * E10_1 * t(mx0[i,]))*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[k:(k+1)] = gde0or[k:(k+1)] + ((E10_3 * t(yx11[i,1:2]) - E10_3 * t(yx11[i,1:2]) * E00 - E10 * E00_3 * t(yx00[i,1:2]))*(E00 - E00 * E10)-(E00_3 * t(yx00[i,1:2]) - E00_3 * t(yx00[i,1:2]) * E10 - E00 * E10_3 * t(yx11[i,1:2]))*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;
   gde0or[k+2] = gde0or[k+2] + ((E10_4 - E10_4 * E00 - E10 * E00_4)*(E00 - E00 * E10)-(E00_4 - E00_4 * E10 - E00 * E10_4)*(E10 - E00 * E10))/((E10/(1-E10))/(E00/(1-E00)))/(E00*(1-E10))**2;

   gie0or[1:(k-1)] = gie0or[1:(k-1)] + ((E01_1 * t(mx1[i,]) - E01_1 * t(mx1[i,]) * E00 - E01 * E00_1* t(mx0[i,]))*(E00 - E00 * E01)-(E00_1* t(mx0[i,]) - E00_1* t(mx0[i,]) * E01 - E00 * E01_1 * t(mx1[i,]))*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[k:(k+1)] = gie0or[k:(k+1)] + ((E01_3 * t(yx00[i,1:2]) - E01_3 * t(yx00[i,1:2]) * E00 - E01 * E00_3 * t(yx00[i,1:2]))*(E00 - E00 * E01)-(E00_3 * t(yx00[i,1:2]) - E00_3 * t(yx00[i,1:2]) * E01 - E00 * E01_3 * t(yx00[i,1:2]))*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;
   gie0or[k+2] = gie0or[k+2] + ((E01_4 - E01_4 * E00 - E01 * E00_4)*(E00 - E00 * E01)-(E00_4 - E00_4 * E01 - E00 * E01_4)*(E01 - E00 * E01))/((E01/(1-E01))/(E00/(1-E00)))/(E00*(1-E01))**2;

   gde1or[1:(k-1)] = gde1or[1:(k-1)] + ((E11_1 * t(mx1[i,]) - E11_1 * t(mx1[i,]) * E01 - E11 * E01_1* t(mx1[i,]))*(E01 - E01 * E11)-(E01_1* t(mx1[i,]) - E01_1* t(mx1[i,]) * E11 - E01 * E11_1 * t(mx1[i,]))*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[k:(k+1)] = gde1or[k:(k+1)] + ((E11_3 * t(yx11[i,1:2]) - E11_3 * t(yx11[i,1:2]) * E01 - E11 * E01_3 * t(yx00[i,1:2]))*(E01 - E01 * E11)-(E01_3 * t(yx00[i,1:2]) - E01_3 * t(yx00[i,1:2]) * E11 - E01 * E11_3 * t(yx11[i,1:2]))*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;
   gde1or[k+2] = gde1or[k+2] + ((E11_4 - E11_4 * E01 - E11 * E01_4)*(E01 - E01 * E11)-(E01_4 - E01_4 * E11 - E01 * E11_4)*(E11 - E01 * E11))/((E11/(1-E11))/(E01/(1-E01)))/(E01*(1-E11))**2;

end;

gr1 = (gie1 * te - gte * ie1)/te ** 2;
gr1var = gr1 * var  * t(gr1);
if gr1var >= 0 then gr1sd = sqrt(gr1var);

gr0 = (gie0 * te - gte * ie0)/te ** 2;
gr0var = gr0 * var * t(gr0);
if gr0var >= 0 then gr0sd = sqrt(gr0var);

gr1or = (gie1or * teor - gteor * ie1or)/teor ** 2;
gr1orvar = gr1or * var  * t(gr1or);
if gr1orvar >= 0 then gr1orsd = sqrt(gr1orvar);

gr0or = (gie0or * teor - gteor * ie0or)/teor ** 2;
gr0orvar = gr0or * var  * t(gr0or);
if gr0orvar >= 0 then gr0orsd = sqrt(gr0orvar);

sum11 = sum11/n;
sum10 = sum10/n;
sum01 = sum01/n;
sum00 = sum00/n;

TE = TE/n;
IE1 = IE1/n;
DE0 = DE0/n;
IE0 = IE0/n;
DE1 = DE1/n;

R1 = IE1/TE;
R0 = IE0/TE;

TEOR = TEOR * (1/n);
IE1OR = IE1OR * (1/n);
DE0OR = DE0OR * (1/n);
IE0OR = IE0OR * (1/n);
DE1OR = DE1OR * (1/n);

R1OR = IE1OR/TEOR;
R0OR = IE0OR/TEOR;

gte = gte/n;
gtevar = gte * var * t(gte);
if gtevar >= 0 then gtesd = sqrt(gtevar);
gie1 = gie1/n;
gie1var = gie1 * var * t(gie1);
if gie1var >= 0 then gie1sd = sqrt(gie1var);
gde0 = gde0/n;
gde0var = gde0 * var * t(gde0);
if gde0var >= 0 then gde0sd = sqrt(gde0var);
gie0 = gie0/n;
gie0var = gie0 * var * t(gie0);
if gie0var >= 0 then gie0sd = sqrt(gie0var);
gde1 = gde1/n;
gde1var = gde1 * var * t(gde1);
if gde1var >= 0 then gde1sd = sqrt(gde1var);

gteor = gteor/n;
gteorvar = gteor * var * t(gteor);
if gteorvar >= 0 then gteorsd = sqrt(gteorvar);
gie1or = gie1or/n;
gie1orvar = gie1or * var * t(gie1or);
if gie1orvar >= 0 then gie1orsd = sqrt(gie1orvar);
gde0or = gde0or/n;
gde0orvar = gde0or * var * t(gde0or);
if gde0orvar >= 0 then gde0orsd = sqrt(gde0orvar);
gie0or = gie0or/n;
gie0orvar = gie0or * var * t(gie0or);
if gie0orvar >= 0 then gie0orsd = sqrt(gie0orvar);
gde1or = gde1or/n;
gde1orvar = gde1or * var * t(gde1or);
if gde1orvar >= 0 then gde1orsd = sqrt(gde1orvar);

print sum11 sum10 sum01 sum00 TE IE1 DE0 IE0 DE1 TEOR IE1OR DE0OR IE0OR DE1OR gtesd gie1sd gde0sd gie0sd gde1sd
gteorsd gie1orsd gde0orsd gie0orsd gde1orsd;

postprobs=sum11||sum10||sum01||sum00||TE||IE1||DE0||IE0||DE1||TEOR||IE1OR||DE0OR||IE0OR||DE1OR||gtesd||gie1sd||gde0sd||gie0sd||gde1sd
||gteorsd||gie1orsd||gde0orsd||gie0orsd||gde1orsd||r1||r0||r1or||r0or||gr1sd||gr0sd||gr1orsd||gr0orsd;
 cname = {"sum11" "sum10" "sum01" "sum00" "TE" "IE1" "DE0" "IE0" "DE1" "LOGTEOR" "LOGIE1OR" "LOGDE0OR" "LOGIE0OR" "LOGDE1OR" "gtesd" "gie1sd" "gde0sd" "gie0sd" 
"gde1sd" "gteorsd" "gie1orsd" "gde0orsd" "gie0orsd" "gde1orsd" "r1" "r0" "r1or" "r0or" "gr1sd" "gr0sd" "gr1orsd" "gr0orsd"};

create &out from postprobs  [ colname=cname ];
append from postprobs;

quit;

%mend;

%macro conmed(analysis = ., dataset=' ', y=' ', x=' ', m = '', w = '', outrd =' ', outor = '');

%if &analysis ne 1 and &analysis ne 2 %then %do;
  proc iml;
    print,"WARNING: NEED TO SPECIFY CORRECT NUMBER FOR VARIABLE ANALYSIS, 1 OR 2","PROGRAM WILL TERMINATE",;
  quit;
%end;

%else %if &dataset= | &y= | &x= | &m = %then %do;
    proc iml;
      print,"WARNING: NEED TO SPECIFY DATASET, PREDICTOR, MEDIATOR AND OUTCOME","PROGRAM WILL TERMINATE",;
	quit;
%end;

%else %if &w= and &analysis=1 %then %do;

%tobit0(dsn=&dataset, y=&m, x=&x, out1 = ff11, out2 = ff12);

%conmedcon1(dataset=&dataset, y=&y, x=&x, m =&m, out = out);

data &outrd;
retain te gtesd ie0 gie0sd de1 gde1sd ie1 gie1sd de0 gde0sd r1 gr1sd r0 gr0sd;
set out;
keep te gtesd ie0 gie0sd de1 gde1sd ie1 gie1sd de0 gde0sd r1 gr1sd r0 gr0sd;
run;

data &outor;
retain logteor gteorsd logie0or gie0orsd logde1or gde1orsd logie1or gie1orsd logde0or gde0orsd r1or gr1orsd r0or gr0orsd;
set out;
keep logteor gteorsd logie0or gie0orsd logde1or gde1orsd logie1or gie1orsd logde0or gde0orsd r1or gr1orsd r0or gr0orsd;
run;

proc datasets lib = work;
delete ff11 ff12 ff13 aa1 aa2 aa11 aa12 aa21 aa3 aa4 aa5 aa6 out;
run;

%end;

%else %if &w= and &analysis=2 %then %do;

%conmedbi1(dataset=&dataset, y=&y, x=&x, m =&m, out = out);

data &outrd;
retain te gtesd ie0 gie0sd de1 gde1sd ie1 gie1sd de0 gde0sd r1 gr1sd r0 gr0sd;
set out;
keep te gtesd ie0 gie0sd de1 gde1sd ie1 gie1sd de0 gde0sd r1 gr1sd r0 gr0sd;
run;

data &outor;
retain logteor gteorsd logie0or gie0orsd logde1or gde1orsd logie1or gie1orsd logde0or gde0orsd r1or gr1orsd r0or gr0orsd;
set out;
keep logteor gteorsd logie0or gie0orsd logde1or gde1orsd logie1or gie1orsd logde0or gde0orsd r1or gr1orsd r0or gr0orsd;
run;

proc datasets lib = work;
delete aa1 aa2 aa3 aa21 aa4 aa5 aa6 bb1 bb2 bb3 bb4 bb5 bb6 out;
run;

%end;

%else %if &analysis=1 %then %do;

%tobit0(dsn=&dataset, y=&m, x=&x &w, out1 = ff11, out2 = ff12);

%conmedcon(dataset=&dataset, y=&y, x=&x, m =&m, w =&w, out = out);

data &outrd;
retain te gtesd ie0 gie0sd de1 gde1sd ie1 gie1sd de0 gde0sd r1 gr1sd r0 gr0sd;
set out;
keep te gtesd ie0 gie0sd de1 gde1sd ie1 gie1sd de0 gde0sd r1 gr1sd r0 gr0sd;
run;

data &outor;
retain logteor gteorsd logie0or gie0orsd logde1or gde1orsd logie1or gie1orsd logde0or gde0orsd r1or gr1orsd r0or gr0orsd;
set out;
keep logteor gteorsd logie0or gie0orsd logde1or gde1orsd logie1or gie1orsd logde0or gde0orsd r1or gr1orsd r0or gr0orsd;
run;

proc datasets lib = work;
delete ff11 ff12 ff13 aa1 aa2 aa11 aa12 aa21 aa3 aa4 aa5 aa6 out;
run;

%end;

%else %if &analysis=2 %then %do;

%conmedbi(dataset=&dataset, y=&y, x=&x, m =&m, w =&w, out = out);

data &outrd;
retain te gtesd ie0 gie0sd de1 gde1sd ie1 gie1sd de0 gde0sd r1 gr1sd r0 gr0sd;
set out;
keep te gtesd ie0 gie0sd de1 gde1sd ie1 gie1sd de0 gde0sd r1 gr1sd r0 gr0sd;
run;

data &outor;
retain logteor gteorsd logie0or gie0orsd logde1or gde1orsd logie1or gie1orsd logde0or gde0orsd r1or gr1orsd r0or gr0orsd;
set out;
keep logteor gteorsd logie0or gie0orsd logde1or gde1orsd logie1or gie1orsd logde0or gde0orsd r1or gr1orsd r0or gr0orsd;
run;

proc datasets lib = work;
delete aa1 aa2 aa3 aa21 aa4 aa5 aa6 bb1 bb2 bb3 bb4 bb5 bb6 out;
run;

%end;

%mend;

* Data set ABSWEI contains the nodes and weights of the Gauss-Hermite quadratrue used for the      ;
* integration of normally distributed continuous mediator                                          ;

data abswei;
input k abs weight;
datalines;
1	-8.098761139250850052013	2.5910437138470814735E-29
2	-7.411582531485468809439	8.544056963775510774E-25
3	-6.84023730524935541785	2.5675933654116696605E-21
4	-6.328255351220081955657	1.9891810121165024856E-18
5	-5.85409505603040010804	6.0083587894908166903E-16
6	-5.4066542479701276084	8.80570764521613225662E-14
7	-4.979260978545255871627	7.1565280526903187084E-12
8	-4.567502072844394855169	3.52562079136541190292E-10
9	-4.16825706683250020154	1.12123608322758101745E-8
10	-3.779206753435223493119	2.41114416367052344179E-7
11	-3.398558265859628346294	3.6315761506930235118E-6
12	-3.024879883901284437677	3.93693398109249277043E-5
13	-2.656995998442895794981	3.13853594541331475647E-4
14	-2.293917141875083421885	0.001871496829597952779484
15	-1.934791472282295793298	0.00846088800825813243994
16	-1.578869894931613886258	0.0293125655361723698457
17	-1.225480109046289030949	0.0784746058654043913089
18	-0.8740066123570880774379	0.1633787327132714570815
19	-0.5238747138322771926149	0.265728251877377076143
20	-0.1745372145975823834895	0.338643277425589218202
21	0.174537214597582383489	0.3386432774255892182024
22	0.523874713832277192615	0.265728251877377076143
23	0.874006612357088077438	0.163378732713271457082
24	1.225480109046289030949	0.0784746058654043913089
25	1.578869894931613886258	0.0293125655361723698457
26	1.934791472282295793298	0.00846088800825813243994
27	2.293917141875083421885	0.00187149682959795277948
28	2.656995998442895794981	3.13853594541331475647E-4
29	3.024879883901284437677	3.93693398109249277043E-5
30	3.398558265859628346294	3.6315761506930235118E-6
31	3.77920675343522349312	2.41114416367052344179E-7
32	4.168257066832500201536	1.12123608322758101745E-8
33	4.56750207284439485517	3.5256207913654119029E-10
34	4.97926097854525587163	7.1565280526903187084E-12
35	5.4066542479701276084	8.8057076452161322566E-14
36	5.85409505603040010804	6.0083587894908166903E-16
37	6.328255351220081955657	1.9891810121165024856E-18
38	6.84023730524935541785	2.5675933654116696605E-21
39	7.41158253148546880944	8.54405696377551077388E-25
40	8.098761139250850052013	2.5910437138470814735E-29
;
run;

proc import datafile = 'C:\Users\Wei Wang\Box Sync\Con Exp Med\2-programs\SAS Macro Upload\example.xls' out = ttemp1 replace;
run; 

%conmed(analysis = 1, dataset=ttemp1, y=outcome, x=exposure, m =mediator1, w =cov1 cov2 cov3, outrd =outrd1, outor = outor1);
%conmed(analysis = 1, dataset=ttemp1, y=outcome, x=exposure, m =mediator1, w =, outrd =outrd2, outor = outor2);
%conmed(analysis = 2, dataset=ttemp1, y=outcome, x=exposure, m =mediator2, w =cov1 cov2 cov3, outrd =outrd3, outor = outor3);
%conmed(analysis = 2, dataset=ttemp1, y=outcome, x=exposure, m =mediator2, w =, outrd =outrd4, outor = outor4);
