library(lme4)
library(splines)
library(ggplot2)

quantiles_time_years = round(quantile(prias_long[!is.na(prias_long$dre),]$visitTimeYears, probs = seq(0,1, by=0.01)), 3)
quantiles_time_years = unique(quantiles_time_years)

odds = vector("numeric", length(quantiles_time_years))
quantiles_time_years = c(-1, quantiles_time_years)
for(i in 2:length(quantiles_time_years)){
  prias_df_temp = prias_long[prias_long$visitTimeYears<=quantiles_time_years[i] & 
                               prias_long$visitTimeYears > quantiles_time_years[i-1],]
  odds[i-1] = table(prias_df_temp$didre)["T1"]/table(prias_df_temp$didre)[">T1"]
}

qplot(y=log(odds), x=quantiles_time_years[-1], geom=c("point","line", "smooth")) + 
  ticksX(1.0, 0.05) + xlim(0,0.9) + xlab("Time (years)")

ggplot(data = prias_long, aes(y=as.numeric(didre),x=visitTimeYears)) + 
  geom_point() + stat_smooth() + ticksX(1.5, 0.05)

prias_long = cbind(prias_long, nsFixed_dre=ns(prias_long$visitTimeYears, knots=0.3),
                   nsRandom_dre=ns(prias_long$visitTimeYears, knots=0.3))
prias_long = cbind(prias_long, nsFixed_auto_dre=ns(prias_long$visitTimeYears, df=2),
                   nsRandom_auto_dre=ns(prias_long$visitTimeYears, df=2))
prias_long = cbind(prias_long, nsFixed_auto_dre=ns(prias_long$visitTimeYears, knots=0.25),
                   nsRandom_auto_dre=ns(prias_long$visitTimeYears, knots=0.25))
write.csv(prias_long, "prias_long.csv", row.names = F)

#Model fitting
prias_long$firstvisit = ifelse(prias_long$visit_number==1, 1, 0)
prias_long$didre_num = ifelse(prias_long$didre=="Low", 0, 1)

mvglm_dre = mvglmer(list(didre ~ Age + ns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5)) + 
                           (1|P_ID)),
                    data = prias_long[prias_long$P_ID!=957,], families = list(binomial), engine="STAN")
jointModel_prias_dre = mvJointModelBayes(mvglm_dre, coxModel, timeVar = "visitTimeYears")

forms_prias_dre <- list("didre" = "value",
                    "didre" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5)), 
                                   random = ~0,
                                   indFixed = 3:4, 
                                   name = "slope"))

jointModel_slope_prias_dre <- update(jointModel_prias_dre, Formulas = forms_prias_dre,
                                 priors = list(shrink_gammas = TRUE, shrink_alphas = TRUE))



glmer(didre ~ I(Age/10) + ns(visitTimeYears, knots=c(0.25), Boundary.knots=c(0,1.5)) + (1|P_ID),
        data=prias_long[!is.na(prias_long$didre),], family="binomial")

model_dre_1 = glmmPQL(didre_num ~  Age + visitTimeDays, random = ~ 1 | P_ID,
        family = "binomial", data = prias_long[!is.na(prias_long$didre),])



anova(model_dre_1, type="marginal")

qplot(y=psaModel$residuals[,2], x=psaModel$fitted[,2])
model.matrix(formula(psaModel), getData(psaModel))
View(round((cor(model.matrix(formula(psaModel), getData(psaModel))[,-1])),3))
kappa(cor(model.matrix(formula(psaModel), getData(psaModel))[,-1]))
