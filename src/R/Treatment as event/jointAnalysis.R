#DRE PSA long, and gleason as survival
dre_psa_data_set = prias_long_survival_eligible[!is.na(prias_long_survival_eligible$dre) & !is.na(prias_long_survival_eligible$psa),]
dre_psa_data_set$P_ID = droplevels(dre_psa_data_set$P_ID)
dre_psa_data_set$didre_num = ifelse(dre_psa_data_set$didre=="T1", 0, 1)
dre_psa_data_set$logpsa1 = log(dre_psa_data_set$psa + 1)
save.image("dre_psa_gl.Rdata")

mvglm_psa_dre_segmented_jags = mvglmer(list(didre_num~Age + visitTimeYears +
                                              I((visitTimeYears-1)*(visitTimeYears>1)) + 
                                              (1|P_ID),

                                              logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + visitTimeYears + 
    											(visitTimeYears|P_ID)),
                                       data = dre_psa_data_set, families = list(binomial, gaussian), engine="JAGS")
save.image(file = "dre_psa_gl.Rdata")
mvJoint_psa_dre_segmented_jags = mvJointModelBayes(mvglm_psa_dre_segmented_jags, coxModel, timeVar = "visitTimeYears")
save.image(file = "dre_psa_gl.Rdata")

mvglm_psa_spline_dre_segmented_jags = mvglmer(list(didre_num~Age + visitTimeYears +
                                              I((visitTimeYears-1)*(visitTimeYears>1)) + 
                                              (1|P_ID),
                                            
                                              logpsa1 ~  I(Age - 70) +  I((Age - 70)^2) + 
                                                ns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0,15)) + 
                                                (ns(visitTimeYears, knots=c(1), Boundary.knots=c(0,15))|P_ID)),
                                       data = dre_psa_data_set, families = list(binomial, gaussian), engine="JAGS")
save.image("dre_psa_gl.Rdata")

mvJoint_psa_spline_dre_segmented_jags = mvJointModelBayes(mvglm_psa_spline_dre_segmented_jags, coxModel, timeVar = "visitTimeYears", 
                                Formulas = list("didre_num"="value",
                                                "logpsa1" = "value",
                                                "logpsa1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0,15)),
                                                                 random=~0 + dns(visitTimeYears, knots=c(1), Boundary.knots=c(0,15)),
                                                                 indFixed = 4:7, indRandom=2:3, name = "slope")))



##################
gleason_dre_dataset = prias_long[!is.na(prias_long$didre) & !is.na(prias_long$digleason),]
gleason_dre_dataset$P_ID = droplevels(gleason_dre_dataset$P_ID)
gleason_dre_dataset$digleason_num = ifelse(gleason_dre_dataset$digleason=="Low", 0, 1)
gleason_dre_dataset$didre_num = ifelse(gleason_dre_dataset$didre=="T1", 0, 1)

mvglm_dre_gleason_segmented = mvglmer(list(
  
  didre_num ~ Age + visitTimeYears + I((visitTimeYears-0.15)*(visitTimeYears>0.15)) + 
    (1|P_ID),
  
  digleason_num ~ Age + visitTimeYears + I((visitTimeYears-0.05)*(visitTimeYears>0.05)) + 
    (1|P_ID)
  
),data = gleason_dre_dataset, 
families = list(binomial, binomial))

jointModel_prias_dre_gleason = mvJointModelBayes(mvglm_dre_gleason_segmented, coxModel, timeVar = "visitTimeYears")

save.image("gleason_dre.Rdata")

model_all = mvglmer(list(
  
  didre ~ Age + ns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5)) +
    (ns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5))|P_ID),
  
  digleason ~ Age + ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05)) +
    (ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05))|P_ID),
  
  logpsa1 ~ Age  +
    ns(visitTimeYears, knots=c(0.1, 0.35, 0.6, 0.85), Boundary.knots = c(0,1.5)) * Age+
    (ns(visitTimeYears, knots=c(0.1), Boundary.knots = c(0,1.5))|P_ID)
  
),data = prias_long[prias_long$P_ID!=957,], 
families = list(binomial, binomial, gaussian))



jointModel_prias = mvJointModelBayes(model_all, coxModel, timeVar = "visitTimeYears")

forms_prias <- list("logpsa1" = "value",
                             "logpsa1" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.35, 0.6, 0.85), Boundary.knots = c(0,1.5)),
                                                      random = ~ 0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots = c(0,1.5)), indFixed = 3:7, indRandom = 2:3,
                                                      name = "slope"),

                             "digleason" = "value",
                             "digleason" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05)), 
                                               random = ~ 0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05)), indFixed = 3:4, indRandom = 2:3, 
                                               name = "slope"),
                             
                             "didre" = "value",
                             "didre" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5)), 
                                                random = ~ 0 + dns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5)), 
                                             indFixed = 3:4, indRandom = 2:3, 
                                                name = "slope"))

jointModel_slope_prias <- update(jointModel_prias, Formulas = forms_prias,
                                         priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))
