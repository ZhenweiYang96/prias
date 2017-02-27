PROC IMPORT OUT= WORK.prias 
            DATAFILE= "C:\Users\838035\Documents\ErasmusMC_datasets\2016-10 PR
IAS DB - DRE, PSA and biopsy results..sav" 
            DBMS=SPSS REPLACE;

RUN;

PROC IMPORT OUT= WORK.prias 
            DATAFILE= "C:\Users\Anirudh\Documents\Office Docs\2016-10 PR
IAS DB - DRE, PSA and biopsy results..sav" 
            DBMS=SPSS REPLACE;

RUN;

data prias;
set prias;
Gleason_sum_2 = Gleason1_2 + Gleason2_2;
if Age ~= . then output;
run;

data prias;
set prias;
if Gleason_sum >= 8 then Gl_biop1_grp = 'High Risk (Gl >= 8)';
else if Gleason_sum = 7 then Gl_biop1_grp = 'Med Risk (Gl = 7)';
else if Gleason_sum = 0 then Gl_biop1_grp = 'Missing';
else Gl_biop1_grp = 'Low Risk (Gl <=6)';

if Gleason_sum_2 >= 8 then Gl_biop2_grp = 'High Risk (Gl >= 8)';
else if Gleason_sum_2 = 7 then Gl_biop2_grp = 'Med Risk (Gl = 7)';
else if Gleason_sum_2 = 0 then Gl_biop2_grp = 'Missing';
else Gl_biop2_grp = 'Low Risk (Gl <=6)';

run;


proc tabulate data = prias;
class Gl_biop1_grp Gl_biop2_grp;
table Gl_biop1_grp all='Total', Gl_biop2_grp all='Total';
run;


data prias_both_biopsy;
set prias;
if Gleason_sum_2>0 and Gleason_sum>0 then output;
run;
