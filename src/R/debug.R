joint_psa = jointModelBayes(psaModel_spline, coxModel, timeVar = "visitTimeYears",
                            param = "td-both",
                            extraForm = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(1, 2, 4), Boundary.knots=c(0,15)),
                                             random=~0 + dns(visitTimeYears, knots=c(1), Boundary.knots=c(0,15)),
                                             indFixed = 4:7, indRandom=2:3))
