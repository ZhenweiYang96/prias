tt = apply(mvJoint_psa_spline_pt1pt54_pt1_tdboth$mcmc$b, MARGIN = c(1,2), FUN = mean)

times = generateLongtiudinalTimeBySchedule()[1:35]

fixedSlopeMatrix = c(dns(times, knots = c(0.1, 0.5, 4), Boundary.knots = c(0, 7)) %*% getBetas(mvJoint_psa_spline_pt1pt54_pt1_tdboth)[4:7])
ranSlopeMatrix= t(dns(times, knots = c(0.1), Boundary.knots = c(0, 7)) %*% t(tt[,2:3]))

finalVelocity = t(apply(ranSlopeMatrix, MARGIN = 1, FUN = function(x){x + fixedSlopeMatrix}))

quantile(c(finalVelocity))
