PROC IMPORT OUT= WORK.PRIAS_ID
            DATAFILE= "C:\Users\838035\Dropbox\PhD\src\jmphd\prias\prias
_id.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

proc phreg data=prias_id;
model year_discontinued*event_indicator(0) = Age Age*Age/ type3;
run;

/* Cox snell */
proc phreg data=prias_id;
model year_discontinued*event_indicator(0) = Age Age*Age/ type3;
output out = phreg_coxsnell LOGSURV = h;
run;

data phreg_coxsnell;
  set phreg_coxsnell;
  h = -h;  
run;

proc lifetest data=phreg_coxsnell nelson outsurv=survres;
time h*event_indicator(0);
run;

data survres;
set survres;
ls=(-log(survival));
run;

proc sgplot data=survres NOAUTOLEGEND;
step x=h y=ls;
xaxis label="Cox Snell residuals";
yaxis label="Estimated cumulative hazard rate";
lineparm x=0 y=0 slope=1 /lineattrs=(color=black);
run;

/*Longitudinal analysis of gleason */
PROC IMPORT OUT= WORK.gleason_long 
            DATAFILE= "D:\Dropbox\PhD\src\prias\gleason_long.csv"
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

proc glimmix data=gleason_long method=quad(qpoints=30);
class P_ID;
nloptions maxiter=1000;
model digleason(event="High") = Age segmentVisitTimeYears visitTimeYears/ dist=BINARY link=logit solution corrb;
random intercept/ subject=P_ID type=un;
title "Random Intercept and no splines";
run;

proc glimmix data=gleason_long method=quad(qpoints=5);
class P_ID;
nloptions maxiter=1000;
model digleason(event="Low") = Age fixedSplines_1 fixedSplines_2 fixedSplines_3/ dist=BINARY link=logit solution corrb;
random intercept/ subject=P_ID type=un;
title "Random Intercept and no splines";
contrast 'Test main effect splines' fixedSplines_1 -1 1, fixedSplines_2 -1 1, fixedSplines_3 -1 1;
run;

/*Longitudinal analysis of gleason */
PROC IMPORT OUT= WORK.dre_long 
            DATAFILE= "D:\Dropbox\PhD\src\prias\dre_long.csv"
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

proc glimmix data=dre_long method=quad(qpoints=5);
class P_ID;
nloptions maxiter=1000;
model didre(event=">T1") = Age segmentVisitTimeYears visitTimeYears/ dist=BINARY link=logit solution corrb;
random intercept/ subject=P_ID type=un;
title "Random Intercept and no splines";
run;

/*Longitudinal analysis of PSA*/
proc mixed data=prias_long method=REML;
class P_ID;
model logpsa1 = polyage_1 polyage_2 nsFixed_1-nsFixed_4 
nsFixed_1*polyage_1 nsFixed_2*polyage_1 nsFixed_3*polyage_1 nsFixed_4*polyage_1/solution;
random intercept nsRandom_1 nsRandom_2/type=un subject=P_ID vcorr=1 v=1 g gcorr;
/*random intercept/type=un subject=P_ID vcorr=1 v=1 g gcorr;*/
contrast 'Test main effect splines' nsFixed_1 -1 1, nsFixed_2 -1 1, nsFixed_3 -1 1, nsFixed_4 -1 1;
contrast 'Test interaction effect splines' nsFixed_1*polyage_1 -1 1, nsFixed_2*polyage_1 -1 1, nsFixed_3*polyage_1 -1 1, nsFixed_4*polyage_1 -1 1;
run;

PROC IMPORT OUT= WORK.PRIAS_LONG 
            DATAFILE= "C:\Users\838035\Dropbox\PhD\src\jmphd\prias\prias
_long.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;


data prias_long;
set prias_long;
if didre = "NA" then didre="";
if digleason = "NA" then digleason="";
intercept100 = 0.1;
run;

proc glimmix data=prias_long method=quad(qpoints=5);
class P_ID;
nloptions maxiter=200;
model digleason(event="Low") = intercept100 Age visitTimeYears/ dist=BINARY link=logit solution corrb noint;
random intercept100 visitTimeYears/ subject=P_ID type=un;
title "Gleason: Random Slope and no splines";
run;

proc glimmix data=prias_long method=quad(qpoints=5);
class P_ID;
nloptions maxiter=100;
model digleason(event="Low") = Age visitTimeYears/ dist=BINARY link=logit solution corrb;
random intercept / subject=P_ID type=un;
title "Random Intercept and no splines";
run;

proc glimmix data=prias_long method=quad(qpoints=5);
class P_ID;
nloptions maxiter=100;
model digleason(event="Low") = Age visitTimeYears/ dist=BINARY link=logit solution corrb;
random intercept visitTimeYears/ subject=P_ID type=un;
title "Random Intercept and no splines";
run;
