## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## ----setup--------------------------------------------------------------------
library(admix)

## -----------------------------------------------------------------------------
####### Under the null hypothesis H0.
## Parameters of the gaussian distribution to be tested:
list.comp <- list(f = 'norm', g = 'norm')
list.param <- list(f = c(mean = 2, sd = 0.5), g = c(mean = 0, sd = 1))
## Simulate data:
obs.data <- rsimmix(n = 700, unknownComp_weight = 0.8, comp.dist = list.comp,
                    comp.param = list.param)[['mixt.data']]
## Performs the test:
list.comp <- list(f = NULL, g = 'norm')
list.param <- list(f = NULL, g = c(mean = 0, sd = 1))
gaussianity_test(sample1 = obs.data, comp.dist = list.comp, comp.param = list.param,
                 K = 3, lambda = 0.1, support = 'Real')$rejection_rule

## -----------------------------------------------------------------------------
##### Under the null hypothesis H0.
## Simulate data:
list.comp <- list(f1 = 'norm', g1 = 'norm',
                  f2 = 'norm', g2 = 'norm')
list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
                   f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 6, sd = 1.2))
sample1 <- rsimmix(n=1000, unknownComp_weight=0.8, comp.dist = list(list.comp$f1,list.comp$g1),
                   comp.param = list(list.param$f1,list.param$g1))[['mixt.data']]
sample2 <- rsimmix(n=1100, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
                   comp.param = list(list.param$f2,list.param$g2))[['mixt.data']]
plot_mixt_density(samples = list(sample1,sample2), user.bounds = NULL, support='continuous')
##### Performs the test:
list.comp <- list(f1 = NULL, g1 = 'norm',
                  f2 = NULL, g2 = 'norm')
list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
                   f2 = NULL, g2 = list(mean = 6, sd = 1.2))
## Using expansion coefficients in orthonormal polynomial basis:
two_samples_test(samples = list(sample1,sample2), comp.dist = list.comp,
                 comp.param = list.param, method = 'Poly', K = 3, support = 'Real',
                 est.method = 'BVdk', s = 0.3, nb.ssEch=2, var.explicit = TRUE)$rejection_rule

## -----------------------------------------------------------------------------
## Simulate data under the null H0:
list.comp <- list(f1 = 'norm', g1 = 'norm',
                  f2 = 'norm', g2 = 'norm',
                  f3 = 'norm', g3 = 'norm')
list.param <- list(f1 = list(mean = 0, sd = 1), g1 = list(mean = 2, sd = 0.7),
                   f2 = list(mean = 0, sd = 1), g2 = list(mean = 4, sd = 1.1),
                   f3 = list(mean = 0, sd = 1), g3 = list(mean = 3, sd = 0.8))
sim1 <- rsimmix(n = 1500, unknownComp_weight = 0.8, comp.dist = list(list.comp$f1,list.comp$g1),
                comp.param = list(list.param$f1, list.param$g1))$mixt.data
sim2 <- rsimmix(n = 1700, unknownComp_weight = 0.7, comp.dist = list(list.comp$f2,list.comp$g2),
                comp.param = list(list.param$f2, list.param$g2))$mixt.data
sim3 <- rsimmix(n = 2000, unknownComp_weight = 0.6, comp.dist = list(list.comp$f3,list.comp$g3),
                comp.param = list(list.param$f3, list.param$g3))$mixt.data
## Perform the 3-samples test in a real-life setting:
list.comp <- list(f1 = NULL, g1 = 'norm',
                  f2 = NULL, g2 = 'norm',
                  f3 = NULL, g3 = 'norm')
list.param <- list(f1 = NULL, g1 = list(mean = 2, sd = 0.7),
                   f2 = NULL, g2 = list(mean = 4, sd = 1.1),
                   f3 = NULL, g3 = list(mean = 3, sd = 0.8))
obj <- IBM_k_samples_test(samples = list(sim1,sim2,sim3), sim_U = NULL, n_sim_tab = 8,
                          comp.dist = list.comp, comp.param = list.param, tune.penalty = TRUE,
                          parallel = FALSE, n_cpu = 2)
obj$rejection_rule

