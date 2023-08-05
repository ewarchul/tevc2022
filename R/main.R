source("cmaes.R")


x0 = rep(95, 2);

sphere_fn = function(x) { sum(x * x) }

ret = cma_es(x0, sphere_fn, lower = -100,  upper = 100)


