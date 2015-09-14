*--------------------------------------------------------------;
* PROGRAM TO CALCULATE RESIDUAL FEED INTAKE                    ;
* PROJECT WITH ERIN CONNOR, BFGL                               ;
*--------------------------------------------------------------;
data all_data;
  infile 'data.prn'  missover firstobs=2;
  input @5 cow_id  4.  @16 trial 1.  @19 reg_num 9.  @33 dim 3.  @41 length 3.  @46 date mmddyy8.  @58 dmi 7.4  @72 lact 1.
        @78 countofdim 3.  @87 minofdim 2.  @94 maxofdim 3.  @104 dailymilk 7.4  @119engy_intake 8.4  @132 weight 6.2  
	@142 pedometer 4.1  @151 time $2.  @155 milk 7.4  @164 fat 6.2  @174 prot 4.2;
  drop pedometer dailymilk length countofdim minofdim maxofdim;
  if dim le 90;
run;

proc sort data=all_data; by cow_id lact; run;    
proc means data=all_data; run;
proc print data=all_data (obs=500); run;


*-------------------------------------------------------------------------;
*  CALCULATE CALVING DATE                                                  ;
*-------------------------------------------------------------------------;
proc sort data=all_data; by cow_id lact date; run;
proc sort data=all_data out=one nodupkeys; by cow_id lact; run;
data one;
  set one;
  calvdt = date - dim;
  keep cow_id lact calvdt;
  run;
  
data all_data;
  merge one(in=in1) all_data(in=in2);
  by cow_id lact;
  if in1 and in2;
  run;


*-------------------------------------------------------------------------;
*  RUN PROC GLM TO OBTAIN BODY WEIGHT REGRESSIONS                         ;
*-------------------------------------------------------------------------;

proc sort data=all_data out=a nodupkeys; by cow_id date; run;

data wgt;
  set a;
  if weight ne .;
  run;
   
proc sort data=wgt; by cow_id lact; run;

ods trace on;
ods listing close;
ods output ParameterEstimates=est (drop=stderr tvalue probt);

proc glm data=wgt;
  model weight = dim dim*dim;
  by cow_id lact;
  run;

ods _all_ close;
ods listing;
ods trace off;

proc transpose data=est out=wgt_est;
  by cow_id lact;
  idlabel parameter;
  var estimate;
  run;

data wgt_est; 
  set wgt_est; 
  Intercept=COL1; EST_dim=COL2; EST_dim2=COL3; 
  drop COL1 COL2 COL3 _NAME_;
run;

proc datasets; delete a wgt est; run; quit;


*-------------------------------------------------------------------------;
*  CALCULATE DAILY MILK                                           ;
*  MAKE ALL_DATA FILE WITH ONLY ONE OBSERVATION PER DATE                  ;                                                
*-------------------------------------------------------------------------;

data a;
  set all_data;
  keep cow_id date milk;
  run;
  
proc sort data=a; by cow_id date; run;
proc means data=a noprint sum; var milk; by cow_id date; output out=milkavg sum=Dailymilk; run;
data milkavg; set milkavg; keep cow_id date dailymilk; run;

proc sort data=all_data out=a1 nodupkeys; by cow_id date; run;

data a1;
  merge a1 (in=in1) milkavg(in=in2);
  by cow_id date;
  drop time milk;
  run;

proc datasets; delete a milkavg; run; quit;


*-------------------------------------------------------------------------;
*  RUN PROC GLM TO OBTAIN MILK COMPONENT REGRESSIONS                      ;
*-------------------------------------------------------------------------;

proc sort data=a1 out=a(keep=cow_id date lact dim dailymilk) nodupkeys; by cow_id date; run;

data milk;
  set all_data;
  if time='AM' then am_pm=1;
  else am_pm=2;
  drop milk time;
  if fat eq . and prot eq . then delete;
run;

proc sort data=milk; by cow_id date; run;
proc sort data=a; by cow_id date; run;

data comp;
  merge a milk(in=in1); 
  by cow_id date;
  if in1;
  
  milk_fat = (fat/100)*dailymilk;
  milk_prot = (prot/100)*dailymilk;
  keep cow_id lact am_pm dim milk_fat milk_prot; 
run;

proc sort data=comp; by cow_id lact; run;
/*proc sort data=a; by cow_id lact; run;*/

ods trace on;
ods listing close;
ods output ParameterEstimates=est (drop=stderr tvalue probt Biased);

proc glm data=comp;
  class am_pm;
  model milk_fat milk_prot = am_pm dim dim*dim / solution;
  by cow_id lact;
  run;
ods _all_ close;
ods listing;
ods trace off;

data milkfatest milkprotest;
  set est; 
  drop dependent; 
  if dependent='milk_fat' then output milkfatest;
  else output milkprotest;
  run;

*MILK FAT ESTIMATES;  
proc transpose data=milkfatest out=fat_est;
  by cow_id lact;
  idlabel parameter;
  var estimate;
  run;

data fat_est; 
  set fat_est; 
  Intercept=COL1; EST_ampm=COL2; EST_dim=COL4; EST_dim2=COL5; 
  drop COL1 COL2 COL3 COL4 COL5 _NAME_;
run;

data fat_est;
  merge a fat_est;
  by cow_id lact;
run;

data fat_est;
  set fat_est;
  EST_ampm = EST_ampm/2;
  pred_milk_fat=intercept + EST_ampm + (EST_dim*dim) + (EST_dim2*dim*dim);
  drop intercept est_ampm est_dim est_dim2;
  run;


*MILK PROTEIN ESTIMATES; 
proc transpose data=milkprotest out=prot_est;
  by cow_id lact;
  idlabel parameter;
  var estimate;
  run;

data prot_est; 
  set prot_est; 
  Intercept=COL1; EST_ampm=COL2; EST_dim=COL4; EST_dim2=COL5; 
  drop COL1 COL2 COL3 COL4 COL5 _NAME_;
run;
 
data prot_est;
  merge a prot_est;
  by cow_id lact;
run;

data prot_est;
  set prot_est;
  EST_ampm = EST_ampm/2;
  pred_milk_prot=intercept + EST_ampm + (EST_dim*dim) + (EST_dim2*dim*dim);
  drop intercept est_ampm est_dim est_dim2;
  run;   

data comp_est;
  merge fat_est prot_est;
  by cow_id lact date; 
  run;

proc datasets; delete a all_data milk est milkfatest milkprotest comp fat_est prot_est; run; quit;  



*-------------------------------------------------------------------------;
*  CALCULATE ECM                                                          ;                                                      
*-------------------------------------------------------------------------;

proc sort data=a1; by cow_id lact date; run;
proc sort data=comp_est; by cow_id lact date; run;

data all_data;
  merge a1(in=in1) comp_est (in=in2);
  by cow_id lact date;
  if in1 and in2;
  run;

     
data all_data;
  set all_data;
  ecm= (0.327*dailymilk) + (12.95*pred_milk_fat) + (7.2*pred_milk_prot); 
  run;
  
proc datasets; delete comp_est a1; run; quit;
  
data _null_;
  set all_data; 
  ecm = round(ecm,0.0001);
  file 'allvars';
  put @6 cow_id    @21 trial  @31 reg_num  @46 dim  @52 date   @69 engy_intake  @86 dmi 
      @96 weight   @132 lact @142 ecm;
run;
 
   
*-------------------------------------------------------------------------;
*  CALCULATE PREDICTED WEIGHT                                             ;                                                  
*-------------------------------------------------------------------------;
proc means data=all_data noprint mean; var weight; by cow_id lact; output out=avgwgt mean=Avg_actualwt; run; 
data avgwgt; set avgwgt; keep cow_id lact Avg_actualwt; run;

data a;
  merge all_data avgwgt;
  by cow_id lact;
  run;
  
proc sort data=a; by cow_id lact; run;
proc sort data=wgt_est; by cow_id lact; run;
data a noid;
  merge a(in=in1) wgt_est(in=in2);
  by cow_id lact;
  if in1 and in2 then output a;
  if in1 and not in2 then output noid;
  run;

proc print data=noid; run;

data all_data;
  set a;
  pred_wt = intercept + (EST_dim*dim) + (EST_dim2*dim*dim);
  drop intercept EST_dim EST_dim2;
  run;
proc datasets; delete wgt_est noid avgwgt a cows; run; quit;

proc sort data=all_data; by dim; run;
proc means data=all_data; by dim; run;
proc means data=all_data noprint mean; var dailymilk; by dim; output out=mean_milk mean=AvgMilk; run;
proc means data=all_data noprint mean; var ecm; by dim; output out=mean_ecm mean=AvgECM; run;
proc print data=mean_milk; run;
proc print data=mean_ecm; run;

proc sort data=all_data; by cow_id dim; run;

data _null_;
  set all_data; 
  if dailymilk ne .;
  ecm = round(ecm,0.0001);
  dailymilk=round(dailymilk,0.0001);
  pred_milk_fat=round(pred_milk_fat,0.0001);
  pred_milk_prot=round(pred_milk_prot,0.0001);
  pred_wt=round(pred_wt,0.0001);
  
  file 'pred_vars';
  put @1 cow_id    @8 reg_num  @20 lact  @25 dim @30 dailymilk @40 pred_milk_fat @50 pred_milk_prot @60 ecm @70 pred_wt;
run;


*-------------------------------------------------------------------------;
*  CALCULATE GAIN WITHIN TIME FRAME                                       ;                                                  
*-------------------------------------------------------------------------;
proc sort data=all_data; by cow_id lact dim; run;
proc sort data=all_data out=first nodupkeys; by cow_id lact; run;
data first; set first; first_weight=pred_wt; first_day=dim; keep cow_id lact first_weight first_day; run;

proc sort data=all_data; by cow_id lact descending dim ; run;
proc sort data=all_data out=last nodupkeys; by cow_id lact; run;
data last; set last; last_weight=pred_wt; last_day=dim; keep cow_id lact last_weight last_day; run;

data weight;
  merge first last;
  by cow_id lact;
  gain=(last_weight - first_weight)/(last_day - first_day);
  keep cow_id lact gain;
  run;

proc print data=weight(obs=10); run;

proc sort data=weight nodupkeys; by cow_id lact; run;

proc datasets; delete first last; run; quit;

data all_data;
  set all_data;
  if date ne .;
  run;
  
*-------------------------------------------------------------------------;
*  CALCULATE AVERAGES                                                     ;                                                  
*-------------------------------------------------------------------------;
  
proc sort data=all_data; by cow_id lact; run;
proc means noprint data=all_data;
  var pred_wt engy_intake dmi ecm;
  by cow_id lact; 
  output out=avgs mean = avg_predwt avg_engy_intake avg_dmi avg_ecm n = n;
  run;


*-------------------------------------------------------------------------;
*  MERGE AVERAGES AND WEIGHT TOGETHER                                     ;                                                  
*-------------------------------------------------------------------------;

data avgs;
  merge avgs weight;
  by cow_id lact;
  run;
  
proc sort data=all_data out=actavg(keep=cow_id lact Avg_actualwt calvdt reg_num) nodupkeys; by cow_id lact; run;

data avgs;
  merge avgs actavg;
  by cow_id lact;
  drop _TYPE_ _FREQ_;
  run;

proc corr data=avgs; var avg_predwt Avg_actualwt; run;

data avgs;
  set avgs;
  avg_predwt75 = avg_predwt**0.75;
  month_calvdt=month(calvdt);
 
   calvyr = year(calvdt);
   if 1 <= month(calvdt) <= 2 then calving_season = 1;
   else if 3 <= month(calvdt) <= 5 then calving_season = 2;
   else if 6 <= month(calvdt) <= 8 then calving_season = 3;
   else if 9 <= month(calvdt) <= 11 then calving_season = 4;
   else calving_season = 1;
   ys=strip(calvyr)!!strip(calving_season);
   drop month_calvdt calvyr calving_season;
  run;

proc freq data=avgs; tables ys; run;
proc datasets; delete weight actavg; run; quit;

*-------------------------------------------------------------------------;
*  ADD LACTATION TO MODEL TO SEE IF SIGNIFICANT                           ;
*  PROC GLM FOR AVERAGE DMI, ALL DATA                                     ;                                                  
*   SEE IF WE CAN USE 1, 2, 3+ cateogires                                 ;
*-------------------------------------------------------------------------;

data lact;
  set avgs;

  if lact=1 then parity=1;
  else if lact=2 then parity=2;
  else parity=3;
  if cow_id=3466 and lact=2 then ys='20124';
run;
proc freq data=lact; tables parity; run;

proc sort data=lact; by lact cow_id; run;


*-------------------------------------------------------------------------;
*  ADD YEAR-SEASON of CALVING TO MODEL TO SEE IF SIGNIFICANT                ;                                                  
*-------------------------------------------------------------------------;

proc glm data=lact ;
  class parity;
  model  avg_engy_intake = parity avg_predwt75 gain avg_ecm / ss3 solution;
  output out=pred_3lacts_ys p = pred r=resid;
  
  estimate  'Diff12' parity 1 -1;
  estimate  'Diff13' parity 1  0 -1;
  estimate  'Diff23' parity 0  1 -1; 

  run;
  
*ys;

proc glm data=lact ;
  class parity ;
  model  avg_dmi = parity avg_predwt75 gain avg_ecm / ss3 solution;
  output out=pred_dmi_3lacts_ys p = pred_dmi r=resid_dmi;

  estimate  'Diff12' parity 1 -1;
  estimate  'Diff13' parity 1  0 -1;
  estimate  'Diff23' parity 0  1 -1; 
  run;




*-------------------------------------------------------------------------;
*  MERGE RESULTS                                                          ;                                                  
*-------------------------------------------------------------------------;



proc sort data=pred_3lacts_ys; by cow_id lact; run;
proc sort data=pred_dmi_3lacts_ys; by cow_id lact; run;
data final;
  merge pred_dmi_3lacts_ys pred_3lacts_ys;
  by cow_id lact;
  pred=round(pred, .001);
  resid=round(resid, .001);
  pred_dmi=round(pred_dmi, .001);
  resid_dmi=round(resid_dmi, .001);
  format calvdt mmddyy10.;
  drop avg_actualwt avg_predwt;
  run;

proc print data=final; title '3 PARITIES INCLUDED IN MODEL:RFI using 3 lactation groups, YEAR-SEASON OF CALVING'; run; title;

data _null_;
  set final;    
  file 'rfi.dat';
  put @1 cow_id @10 reg_num @20 lact @25 pred @35 resid @45 pred_dmi @55 resid_dmi;
run;  
  




quit;
