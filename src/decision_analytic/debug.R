jointModel_debug(lmeFit, jointModelData$survModel_simDs, timeVar = "visitTimeYears", 
           method = "spline-PH-aGH", parameterization = "both",
           derivForm = list(fixed = ~ 0 + dns(visitTimeYears, knots=c(0.1, 0.5, 4), Boundary.knots=c(0, 7)),
                            random=~0 + dns(visitTimeYears, knots=c(0.1), Boundary.knots=c(0, 7)),
                            indFixed = 4:7, indRandom=2:3))