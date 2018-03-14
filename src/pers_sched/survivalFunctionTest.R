set.seed(1000)

shape = 2
scale = 6
sample = rweibull(10000, shape, scale)

plot(density(sample))

trueMean = scale * gamma(1 + 1/shape)
sampleMean = mean(sample)

#Simple check with E(t)
integrate(function(x, shape, scale){
  x * dweibull(x, shape, scale)
}, lower = 0, upper = Inf, shape, scale)


integrate(function(x, shape, scale){
  exp(-(x/scale)^shape)
}, lower = 0, upper = Inf, shape, scale)

#Check for E(T|T>3)
t0 = 5
integrate(function(x, shape, scale){
  (x * dweibull(x, shape, scale)) / (exp(-(t0/scale)^shape))
}, lower = t0, upper = Inf, shape, scale)

t0 + integrate(function(x, shape, scale){
  exp(-(x/scale)^shape) / exp(-(t0/scale)^shape)
}, lower = t0, upper = Inf, shape, scale)$value

