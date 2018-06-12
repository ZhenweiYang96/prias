#Load the dataset containing a training datset and test dataset
load("../dataset.Rdata")

#Install the required packages
install.packages("JMbayes")
library(JMbayes)
library(splines)
library(survival)

#We need to fit a cox model, to extract the model object in the joint model later
cox_Model_training = coxph(Surv(progression_time, progressed) ~ I(Age - 70) +  I((Age - 70)^2),
                           data = trainingDs.id, x=T, model=T)

#mvglmer fits multivariate generalized linear mixed models (Bayesian) using JAGS as well STAN. 
#Default is JAGS
mvglmer_psa_training=mvglmer(list(logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + 
                                    ns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)) + 
                                    (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7))|P_ID)), 
                             data = trainingDs, families = list(gaussian))

#We specify the formula for calculation of velocity of PSA
forms_psa_training <- list("logpsa1" = "value",
                           "logpsa1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)),
                                            random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7)),
                                            indFixed = 4:7, indRandom=2:3, name = "slope"))

#Fitting a joint model
mvJoint_psa_tdboth_training <- mvJointModelBayes(mvglmer_psa_training, cox_Model_training, timeVar = "visitTimeYears",
                                                 Formulas = forms_psa_training)


#Plot the model object to analyze traceplots
plot(mvJoint_psa_tdboth_training)

#Check the summary for parameter estimates
summary(mvJoint_psa_tdboth_training)

