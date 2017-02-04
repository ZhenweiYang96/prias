#subject 957 doesn't have any gleason. removing it

save.image()
model_binary = mvglmer(list(
  
  didre ~ Age + ns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5)) +
    (1|P_ID),
  
  digleason ~ Age + ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05)) +
    (1|P_ID)
  
),data = prias_long[prias_long$P_ID!=957,], 
families = list(binomial, binomial))

save.image()

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

jointModel_prias_gleason = mvJointModelBayes(mvglm_gleason, coxModel, timeVar = "visitTimeYears")

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
