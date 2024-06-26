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
## Simulate data (chosen parameters indicate 3 clusters (populations (1,3), (2,5) and 4)!):
list.comp <- list(f1 = "gamma", g1 = "exp",
                  f2 = "gamma", g2 = "exp",
                  f3 = "gamma", g3 = "gamma",
                  f4 = "exp", g4 = "exp",
                  f5 = "gamma", g5 = "exp")
list.param <- list(f1 = list(shape = 16, rate = 4), g1 = list(rate = 1/3.5),
                   f2 = list(shape = 14, rate = 2), g2 = list(rate = 1/5),
                   f3 = list(shape = 16, rate = 4), g3 = list(shape = 12, rate = 2),
                   f4 = list(rate = 1/2), g4 = list(rate = 1/7),
                   f5 = list(shape = 14, rate = 2), g5 = list(rate = 1/6))
A.sim <- rsimmix(n=6000, unknownComp_weight=0.8, comp.dist = list(list.comp$f1,list.comp$g1),
                 comp.param = list(list.param$f1, list.param$g1))$mixt.data
B.sim <- rsimmix(n=6000, unknownComp_weight=0.5, comp.dist = list(list.comp$f2,list.comp$g2),
                 comp.param = list(list.param$f2, list.param$g2))$mixt.data
C.sim <- rsimmix(n=6000, unknownComp_weight=0.65, comp.dist = list(list.comp$f3,list.comp$g3),
                 comp.param = list(list.param$f3, list.param$g3))$mixt.data
D.sim <- rsimmix(n=6000, unknownComp_weight=0.75, comp.dist = list(list.comp$f4,list.comp$g4),
                 comp.param = list(list.param$f4, list.param$g4))$mixt.data
E.sim <- rsimmix(n=6000, unknownComp_weight=0.55, comp.dist = list(list.comp$f5,list.comp$g5),
                 comp.param = list(list.param$f5, list.param$g5))$mixt.data
## Look for the clusters:
list.comp <- list(f1 = NULL, g1 = "exp",
                  f2 = NULL, g2 = "exp",
                  f3 = NULL, g3 = "gamma",
                  f4 = NULL, g4 = "exp",
                  f5 = NULL, g5 = "exp")
list.param <- list(f1 = NULL, g1 = list(rate = 1/3.5),
                   f2 = NULL, g2 = list(rate = 1/5),
                   f3 = NULL, g3 = list(shape = 12, rate = 2),
                   f4 = NULL, g4 = list(rate = 1/7),
                   f5 = NULL, g5 = list(rate = 1/6))
clusters <- admix_clustering(samples = list(A.sim,B.sim,C.sim,D.sim,E.sim), n_sim_tab = 10,
                             comp.dist=list.comp, comp.param=list.param, parallel=FALSE, n_cpu=2)
clusters

