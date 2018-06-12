# DR SUGGESTION 14/04/2017

# things to check: 
# weight function in simulation: use non-standard (as per input)
# phi value : decreasing phi increases alpha, but reduces n

# DR code: suggested model
library("JM")

pbc2$status2 <- as.numeric(pbc2$status != "alive")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

lmeFit <- lme(log(serBilir) ~ ns(year, 3, B = c(0, 9.3)), 
              data = pbc2,
              random = list(id = pdDiag(form = ~ ns(year, 3, B = c(0, 9.3)))))

survFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

# weight function is standard normal with sigma = 1
wFun <- function (s, t) dnorm(t - s)

wiForm <- list(fixed = ~ 0 + I(pnorm(year) - 0.5) + ins(year, 3, B = c(0, 9.3), 
                                                        weight.fun = wFun), 
               indFixed = 1:4,
               random = ~ 0 + I(pnorm(year) - 0.5) + ins(year, 3, B = c(0, 9.3), 
                                                         weight.fun = wFun), 
               indRandom = 1:4)

jointFit <- jointModel(lmeFit, survFit, timeVar = "year", parameterization = "slope",
                       derivForm = wiForm, method = "weibull-PH-aGH")
summary(jointFit)

# SIMULATING DATA USING ABOVE VALUES----------------------------------------------------------------------------------

random.seed <- runif(1, 1, 10000)
set.seed(random.seed)
#set.seed(184.3527)

library(plyr)
library(matrixcalc)

################################################

n <- length(unique(pbc2.id$id)) + 200 # number of subjects (312 + 200)
K <- max(count(pbc2,c('id'))[,2])*3 # number of planned repeated measurements per subject, per outcome (16)
t.max <- max(pbc2.id$years) # maximum follow-up time (14.30566)

#K <- 45
#t.max <- 40

################################################

# parameters for the linear mixed effects model
betas <- c("Intercept" =  as.numeric(fixef(jointFit))[1], "ns1" = as.numeric(fixef(jointFit))[2], 
           "ns2" = as.numeric(fixef(jointFit))[3], "ns3" = as.numeric(fixef(jointFit))[4])

# measurement error standard deviation
sigma.y <- summary(jointFit)$sigma

# parameters for the survival model
gammas <- c("(Intercept)" = as.numeric(summary(jointFit)$"CoefTable-Event"[1]), 
            "drug" = as.numeric(summary(jointFit)$"CoefTable-Event"[2])) # coefficients for baseline covariates

alpha <- as.numeric(summary(jointFit)$"CoefTable-Event"[3]) # association parameter
#alpha <- 0
#alpha <- 0.6029
shp1 <- 1

phi <- 0.8223

#mean.Cens <- mean(pbc2.id$years[pbc2.id$status=="alive"])

D <- matrix(0, 4, 4)
D[indLowerTri <- lower.tri(D, TRUE)] <- c(summary(jointFit)$D[1], 0, 0, 0,
                                          summary(jointFit)$D[2], 0, 0,
                                          summary(jointFit)$D[3], 0, 
                                          summary(jointFit)$D[4])
D <- D + t(D)
diag(D) <- diag(D) * 0.5

is.positive.definite(D)

################################################

# design matrices for the longitudinal measurement model

times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max)))))

table(pbc2.id$drug)
drug <- rep(0:1, each = n/2)

DF <- data.frame(id = rep(1:n, each = K), year = times, 
                 drug = factor(rep(drug, each = K)))

Bkn <- c(0, 9.3)
#Bkn <- c(0, 19.8)
#Bkn <- c(0, 26.4)

spline <- ns(pbc2$year, 3, B = Bkn)
attr(spline, "knots")

kn <- c(as.numeric(attr(spline, "knots")[1]), as.numeric(attr(spline, "knots")[2]))
#kn <- c(5, 10)
#kn <- c(7, 14)

X <- model.matrix(~ ns(year, knots = kn, Boundary.knots = Bkn), data = DF)
Z <- model.matrix(~ ns(year, knots = kn, Boundary.knots = Bkn), data = DF)

# design matrix for the survival model
W <- model.matrix(~ drug, data = DF[!duplicated(DF$id), ])

################################################

# simulate random effects
b <- mvrnorm(n, rep(0, nrow(D)), D)

# simulate longitudinal responses
id <- rep(1:n, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ])) # linear predictor
y <- rnorm(n * K, eta.y, sigma.y)

# simulate event times
eta.t <- as.vector(W %*% gammas)

# WEIGHTED VALUE SIMULATION-----------------------------------------------------------------------------------------------

# using unstandardized as per input
fff <- function (ss) {
  dnorm(x = ss, sd = shp1)
}

# correct formulation? 
invS_weight1 <- function (t) {
  h <- function (s) {
    w <- function (u) {
      NS <- ns(u, knots = kn, Boundary.knots = Bkn)
      XXi <- cbind(1, NS)
      ZZi <- cbind(1, NS)
      
      mi <- as.vector(XXi %*% betas + rowSums(ZZi * b[rep(i, nrow(ZZi)), ]))
      
      ww <- fff(s - u)
      ww * mi
    }
    ff <- integrate(w, lower = 0, upper = s, subdivisions = 200L, 
                    rel.tol = 0.001, stop.on.error = FALSE)$value
    
    exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + ff * alpha)
  }
  integrate(h, lower = 0, upper = t)$value + log(u[i])
}

u <- runif(n)
trueTimes <- numeric(n)

for (i in 1:n) {
  Up <- 15
  tries <- 3
  Root <- try(uniroot(invS_weight1, interval = c(1e-04, Up))$root, FALSE)
  while(inherits(Root, "try-error") && tries > 0) {
    tries <- tries - 1
    Up <- Up * 2
    Root <- try(uniroot(invS_weight1, interval = c(1e-04, Up))$root, FALSE)
  }
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
}

# alternative version - as per Elrozy's area simulation
############################################################################################################################

#wFun <- function (s, t) dnorm(t - s)

#invS <- function (t, u, i) {
#  h <- function (s) {

#    INS <- ins(s, 3, B = c(0, 9.3), weight.fun = wFun)
#    XXi <- cbind(pnorm(s) - 0.5, INS)
#    ZZi <- cbind(pnorm(s) - 0.5, INS)

#    f1 <- as.vector(XXi %*% betas + rowSums(ZZi * b[rep(i, nrow(ZZi)), ]))
#    exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha)
#  }
#  integrate(h, lower = 0, upper = t)$value + log(u)
#}

#u <- runif(1)
#i <- 1

#> u
#[1] 0.1883626

#invS_weight1(9)
#[1] -1.126471
#invS(9, u, i)
#-1.578524

#u <- runif(n)
#trueTimes <- numeric(n)

#for (i in 1:n) {
#  Up <- 15
#  tries <- 5
#  Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
#  while(inherits(Root, "try-error") && tries > 0) {
#    tries <- tries - 1
#    Up <- Up + 5000
#    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
#  }
#  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
#}

############################################################################################################################

na.ind <- !is.na(trueTimes)
trueTimes <- trueTimes[na.ind] 
median_times <- median(trueTimes)
hist(trueTimes)
length(trueTimes) 
summary(trueTimes)

b <- b[na.ind, ]
W <- W[na.ind, , drop = FALSE]

long.na.ind <- rep(na.ind, each = K)
y <- y[long.na.ind]
X <- X[long.na.ind, , drop = FALSE]
Z <- Z[long.na.ind, , drop = FALSE]
DF <- DF[long.na.ind, ]

n <- length(trueTimes)

# simulate censoring times from a uniform distribution,
# and calculate the observed event times, i.e., min(true event times, censoring times)

mean.Cens <- mean(pbc2.id$years[pbc2.id$status=="alive"])*0.3

Ctimes <- runif(n, 0, 2*mean.Cens)  
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator

hist(Ctimes)
table(event)

################################################

ind <- times[long.na.ind] <= rep(Time, each = K)
y <- y[ind]
X <- X[ind, , drop = FALSE]
Z <- Z[ind, , drop = FALSE]
id <- id[long.na.ind][ind]
id <- match(id, unique(id))

dat <- DF[ind, ]
dat$id <- id
dat$y <- y
dat$Time <- Time[id]
dat$event <- event[id]

########################## censoring events past t.max
#dat$event[dat$event == 1 & dat$Time > t.max] <- 0
#dat$Time[dat$Time > t.max] <- t.max

dat.id <- dat[!duplicated(dat$id), ]

names(dat) <- names(dat.id) <- c("id", "year", "drug", "y", "Time", "event")
summary(tapply(id, id, length))

table(dat.id$event)
n

mean(dat.id$event)
summary(dat.id$Time)

# true values for parameters and random effects
trueValues <- list(betas = betas, tau = 1/sigma.y^2,  
                   alphas = alpha, sigma.t = phi, sigma = sigma.y, inv.D = solve(D), 
                   b = b, shape = shp1)

# running same model as before to check

lmeFit_sim  <- lme(y ~ ns(year, 3, B = c(0, 9.3)), 
                   #B = c(0, 19.8)),
                   #B = c(0, 26.4)),
                   data = dat,
                   random = list(id = pdDiag(form = ~ ns(year, 3, B = c(0, 9.3)))))
#B = c(0, 19.8)))))
#B = c(0, 26.4)))))

survFit_sim  <- coxph(Surv(Time, event) ~ drug, data = dat.id, x = TRUE)

# weight function is standard normal with sigma = 1
wFun <- function (s, t) dnorm(t - s)

wiForm <- list(fixed = ~ 0 + I(pnorm(year) - 0.5) + ins(year, 3, B = c(0, 9.3), 
                                                        #B = c(0, 19.8),
                                                        #B = c(0, 26.4),
                                                        weight.fun = wFun), 
               indFixed = 1:4,
               random = ~ 0 + I(pnorm(year) - 0.5) + ins(year, 3, B = c(0, 9.3),
                                                         #B = c(0, 19.8),
                                                         #B = c(0, 26.4),
                                                         weight.fun = wFun), 
               indRandom = 1:4)

jointFit_sim <- jointModel(lmeFit_sim, survFit_sim, timeVar = "year", parameterization = "slope",
                           derivForm = wiForm, method = "weibull-PH-aGH")
summary(jointFit_sim)

# bayes version
# current value formulation
library(JMbayes)

jointFit_sim_bayes <- jointModelBayes(lmeFit_sim, survFit_sim,
                                      timeVar = "year")
summary(jointFit_sim_bayes)
#plot(jointFit_sim_bayes)

# normal (un-standardized)

wf <- function (u, parms, t.max) {
  num <- dnorm(x = u, sd = parms)
}

jointFit_sim_bayesw <- update(jointFit_sim_bayes, estimateWeightFun = TRUE,
                              weightFun = wf,
                              priorShapes = list(shape1 = dunif),
                              priors = list(priorshape1 = c(0, 5)))
summary(jointFit_sim_bayesw)
plot(jointFit_sim_bayesw)

#coef <- vector("list", 10) 

#modelnum <- 1

#coef[[modelnum]] <- jointFit_sim$coefficients

#modelnum <- modelnum + 1

#coef[[modelnum]] <- jointFit_sim$coefficients

rm(list=setdiff(ls(), c("jointFit", "coef", "modelnum")))

