# PRIAS

PRIAS (https://www.prias-project.org/), is an abbreviation for Prostate Cancer Research International Active Surveillance. Curently it is the world's largest active surveillance program. As part of a 4 year PhD program in Biostatistics I am working towards creating a personalized schedule of biopsies for prostate cancer patients. To this end, I use a dataset which contains anonymized information of patients enrolled in PRIAS. In this repository, I am sharing the code for the statistical model that I use for the aforementioend goal. Furthermore, I am also sharing the code for the personalized schedules that I have created, and code for the simualtion study I performed to evaluate the performance of various schedules.

## Code for the Model
I use a joint model for time to event and longitudinal data to model the evolution of prostate specific antigen levels over time, and to simulatenouly model their association with the hazard of Gleason score reclassification. The R package I use for this purpose is called JMbayes (https://cran.r-project.org/web/packages/JMbayes/JMbayes.pdf). The API I use however are currently not hosted on CRAN, and can be found here:
https://github.com/drizopoulos/JMbayes/blob/master/man/mvJointModelBayes.Rd

Now since public access to PRIAS data is forbidden, you will most probably be not able to run my code. Having said that, the code for the model that I fit can be found here:
For the relative risk part, I first create a model using the well known survreg API (https://stat.ethz.ch/R-manual/R-devel/library/survival/html/survreg.html)
https://github.com/anirudhtomer/prias/blob/master/src/R/Gleason%20as%20event/coxAnalysis_intervalcensored.R

Using the definitions from the survreg model object, I then create a joint model.
https://github.com/anirudhtomer/prias/blob/master/src/R/Gleason%20as%20event/psaLongAnalysis_intervalcensored.R
As you can see I fit multiple models, however the one I use finally is called "joint_psa_spline_pt1pt54_pt1_tdboth".

## Code for the Personalized Schedules, and Simulation Study
Using the fitted model, I simulate evolution of PSA for 1000 patients x 500 datasets. I also generate the time of Gleason reclassification for these simulated patients. The simulation code can be found here (start reading from the bottom to the top):
https://github.com/anirudhtomer/prias/blob/master/src/R/Sim%20Study%20Mixture/simCommon.R

Now for these patients I also generate personalized schedules, yearly/annual schedule and PRIAS schedule of biopsies, which can be found here:
https://github.com/anirudhtomer/prias/blob/master/src/R/Sim%20Study%20Mixture/nbAndOffset.R

Then I run simulations for 500 datasets to check which schedule conducts how many biopsies to detect Gleason reclassification and how much do they overshoot the real Gleason reclassification time of patients. This can be found here:
https://github.com/anirudhtomer/prias/blob/master/src/R/Sim%20Study%20Mixture/SimulateJM.R

In addition I also create personalized schedules for 3 patients from PRIAS, to see how the personalized schedules work with real world data. Ufortunately we do not know the real Gleason reclassification time for any of the PRIAS patients and hence the simualtion study. Nevertheless the code is here:
https://github.com/anirudhtomer/prias/blob/master/src/R/Gleason%20as%20event/biopsyPRIAS.R

To calculate the expected Gleason reclassification time and variance of posterior predictive distribution of Gleason reclassification time, one can use the following code:
https://github.com/anirudhtomer/JMbayes/blob/master/dev/expectedCondFailureTime.R
https://github.com/anirudhtomer/JMbayes/blob/master/dev/varCondFailureTime.R

A more optimized version of these however is present in the file (I use this because it is faster; and only a tiny tiny bit less accurate):
https://github.com/anirudhtomer/prias/blob/master/src/R/Sim%20Study%20Mixture/simCommon.R

Lastly, I generate graphs to see my results, which can be found here:
https://github.com/anirudhtomer/prias/blob/master/src/R/Sim%20Study%20Mixture/produceResults.R





