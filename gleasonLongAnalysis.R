library(lme4)
library(splines)
library(ggplot2)

quantiles_time_years = round(quantile(prias_long[!is.na(prias_long$gleason),]$visitTimeYears, probs = seq(0,1, by=0.005)), 3)
quantiles_time_years = unique(quantiles_time_years)

odds = vector("numeric", length(quantiles_time_years))
quantiles_time_years = c(-1, quantiles_time_years)
for(i in 2:length(quantiles_time_years)){
  prias_df_temp = prias_long[prias_long$visitTimeYears<=quantiles_time_years[i] & 
                               prias_long$visitTimeYears > quantiles_time_years[i-1],]
  odds[i-1] = table(prias_df_temp$digleason)["Low"]/table(prias_df_temp$digleason)["High"]
}

qplot(y=log(odds), x=quantiles_time_years[-1], geom=c("point","line", "smooth")) + 
  ticksX(1.5, 0.05) + xlab("Time (years)")

ggplot(data = prias_long, aes(y=as.numeric(digleason),x=visitTimeYears)) + geom_point() + stat_smooth()

#Model fitting
prias_long$firstvisit = ifelse(prias_long$visit_number==1, 1, 0)
prias_long$digleason_num = ifelse(prias_long$digleason=="Low", 0, 1)

glmer_gleason_int = glmer(digleason ~I(Age/10) + ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05)) + 
        (1|P_ID), data =  prias_long[prias_long$P_ID!=957 & !is.na(prias_long$digleason),], family = binomial)


glmer_gleason_int_glmmpql = glmmPQL(digleason_num ~  I(Age/10) + firstvisit + visitTimeYears, random = ~ 1 | P_ID,
        family = "binomial", data = prias_long[prias_long$P_ID!=957 & !is.na(prias_long$digleason),])

dLongBin <- function (y, eta.y, scale, log = FALSE, data) {
  dbinom(x = y, size = 1, prob = plogis(eta.y), log = log)
}

jointFit_gleason <- jointModelBayes(glmer_gleason_int_glmmpql, coxModel, timeVar = "visitTimeYears", 
                                 densLong = dLongBin)

mvglm_gleason = mvglmer(list(digleason ~Age + ns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05)) + 
                       (1|P_ID)),
        data = prias_long[prias_long$P_ID!=957,], families = list(binomial), engine="STAN")

save.image()
mvglm_gleason_slope = mvglmer(list(digleason ~Age +visitTimeYears + 
                               (1|P_ID)),
                        data = prias_long[prias_long$P_ID!=957,], families = list(binomial), engine="STAN")


joinmodel_gleason = mvJointModelBayes(mvglm_gleason, coxModel, timeVar = "visitTimeYears", 
                  Formulas = list("digleason" = "value",
                                  "digleason" = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0,1.05)),
                                                     random=~0,
                                                     indFixed = 3:4, name = "slope")))

anova(model_dre_1, type="marginal")

qplot(y=psaModel$residuals[,2], x=psaModel$fitted[,2])
model.matrix(formula(psaModel), getData(psaModel))
View(round((cor(model.matrix(formula(psaModel), getData(psaModel))[,-1])),3))
kappa(cor(model.matrix(formula(psaModel), getData(psaModel))[,-1]))
