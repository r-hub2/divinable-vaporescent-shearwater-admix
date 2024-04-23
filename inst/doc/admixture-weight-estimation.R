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
## Simulate data:
list.comp <- list(f = 'norm', g = 'norm')
list.param <- list(f = list(mean = 3, sd = 0.5), g = list(mean = 0, sd = 1))
data1 <- rsimmix(n = 800, unknownComp_weight = 0.7, list.comp, list.param)[['mixt.data']]
## Perform the estimation of parameters in real-life:
list.comp <- list(f = NULL, g = 'norm')
list.param <- list(f = NULL, g = list(mean = 0, sd = 1))
BVdk_estimParam(data1, method = 'L-BFGS-B', list.comp, list.param)

## -----------------------------------------------------------------------------
## Simulate data:
list.comp <- list(f = 'norm', g = 'norm')
list.param <- list(f = list(mean = 3, sd = 0.5), g = list(mean = 0, sd = 1))
data1 <- rsimmix(n = 1000, unknownComp_weight = 0.6, list.comp, list.param)[['mixt.data']]
## Transform the known component of the admixture model into a Uniform(O,1) distribution:
list.comp <- list(f = NULL, g = 'norm')
list.param <- list(f = NULL, g = list(mean = 0, sd = 1))
data1_transfo <- knownComp_to_uniform(data = data1, comp.dist=list.comp, comp.param=list.param)
PatraSen_est_mix_model(data = data1_transfo, method = 'fixed',
                        c.n = 0.1*log(log(length(data1_transfo))), gridsize = 1000)$alp.hat

## -----------------------------------------------------------------------------
## Simulate data:
list.comp <- list(f1 = 'norm', g1 = 'norm',
                  f2 = 'norm', g2 = 'norm')
list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
                   f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 5, sd = 2))
sample1 <- rsimmix(n=1700, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
                                                   comp.param=list(list.param$f1,list.param$g1))
sample2 <- rsimmix(n=1500, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
                                                   comp.param=list(list.param$f2,list.param$g2))
##### On a real-life example (unknown component densities, unknown mixture weights).
list.comp <- list(f1 = NULL, g1 = 'norm',
                  f2 = NULL, g2 = 'norm')
list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
                   f2 = NULL, g2 = list(mean = 5, sd = 2))
## Estimate the mixture weights of the two admixture models (provide only hat(theta)_n):
estim <- IBM_estimProp(sample1 = sample1[['mixt.data']], sample2 = sample2[['mixt.data']],
                       known.prop = NULL, comp.dist = list.comp, comp.param = list.param,
                       with.correction = FALSE, n.integ = 1000)
estim[['prop.estim']]

## -----------------------------------------------------------------------------
## Simulate data:
list.comp <- list(f1 = 'norm', g1 = 'norm',
                  f2 = 'norm', g2 = 'norm')
list.param <- list(f1 = list(mean = 1, sd = 0.5), g1 = list(mean = 0, sd = 1),
                   f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 5, sd = 2))
sample1 <- rsimmix(n=1700, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
                                                   comp.param=list(list.param$f1,list.param$g1))
sample2 <- rsimmix(n=1500, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
                                                   comp.param=list(list.param$f2,list.param$g2))
##### On a real-life example (unknown component densities, unknown mixture weights).
list.comp <- list(f1 = NULL, g1 = 'norm',
                  f2 = NULL, g2 = 'norm')
list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
                   f2 = NULL, g2 = list(mean = 5, sd = 2))
## Estimate the mixture weights of the two admixture models (provide only hat(theta)_n):
estim <- IBM_estimProp(sample1 = sample1[['mixt.data']], sample2 = sample2[['mixt.data']],
                       known.prop = NULL, comp.dist = list.comp, comp.param = list.param,
                       with.correction = FALSE, n.integ = 1000)
estim[['prop.estim']]

## -----------------------------------------------------------------------------
## Simulate data:
list.comp <- list(f1 = 'norm', g1 = 'norm',
                  f2 = 'norm', g2 = 'norm')
list.param <- list(f1 = list(mean = 3, sd = 0.5), g1 = list(mean = 0, sd = 1),
                   f2 = list(mean = 3, sd = 0.5), g2 = list(mean = 5, sd = 2))
sample1 <- rsimmix(n=1700, unknownComp_weight=0.5, comp.dist = list(list.comp$f1,list.comp$g1),
                   comp.param=list(list.param$f1,list.param$g1))
sample2 <- rsimmix(n=1500, unknownComp_weight=0.7, comp.dist = list(list.comp$f2,list.comp$g2),
                   comp.param=list(list.param$f2,list.param$g2))
## Estimate the mixture weight in each of the sample in real-life setting:
list.comp <- list(f1 = NULL, g1 = 'norm',
                  f2 = NULL, g2 = 'norm')
list.param <- list(f1 = NULL, g1 = list(mean = 0, sd = 1),
                   f2 = NULL, g2 = list(mean = 5, sd = 2))
estimate <- IBM_estimProp(sample1[['mixt.data']], sample2[['mixt.data']], comp.dist = list.comp,
                          comp.param = list.param, with.correction = FALSE, n.integ = 1000)
## Determine the decontaminated version of the unknown density by inversion:
res1 <- decontaminated_density(sample1 = sample1[['mixt.data']], comp.dist = list.comp[1:2], 
                               comp.param = list.param[1:2], estim.p = estimate$prop.estim[1])
res2 <- decontaminated_density(sample1 = sample2[['mixt.data']], comp.dist = list.comp[3:4], 
                               comp.param = list.param[3:4], estim.p = estimate$prop.estim[2])
plot(x = res1, type = "l", x_val = seq(from=-1,to=6,length.out=40), add_plot = FALSE)
plot(x = res2, type="l", col="red", x_val = seq(from=-1,to=6,length.out=40), add_plot = TRUE)

