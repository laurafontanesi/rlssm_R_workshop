rm(list = ls())
library(tidyverse)
library(dfoptim)
library(rtdists)
library(rstan)
library(bayesplot)
source('rwald code.r') # from https://osf.io/3sp9t/
source('functions_notebooks.R') # load  functions defined in notebooks

# DDM path
gen_drift = .3
gen_threshold = 1
gen_ndt = .23
sim_path <- dm_path(gen_drift, gen_threshold, gen_ndt)

ggplot(data = sim_path, aes(x = time, y = dv))+
  geom_line(size = .5) +
  geom_hline(yintercept=gen_threshold, size=1.5) +
  geom_hline(yintercept=0, size=1.5)

# RDM path
gen_drift1 = 1.7
gen_drift2 = 2
gen_threshold = 1
gen_ndt = .23
sim_path <- rdm_path(gen_drift1, gen_drift2, gen_threshold, gen_ndt)

ggplot(data = sim_path, aes(x = time, y = dv, color=factor(accumulator)))+
  geom_line(size = .5) +
  geom_hline(yintercept=gen_threshold, size=1.5)

# posterior predictives
sim_data <- rWaldRace(n=1000, v=c(.7, 1.2), B=1.5, A=0, t0=.23, gf = 0)
summary(sim_data)

sim_data_for_stan = list(
  N = dim(sim_data)[1],
  choices = sim_data$R,
  rt = sim_data$RT
)

fit1 <- stan(
  file = "stan_models/RDM.stan",                     # Stan program
  data = sim_data_for_stan,                          # named list of data
  chains = 2,                                        # number of Markov chains
  warmup = 1000,                                     # number of warmup iterations per chain
  iter = 3000,                                       # total number of iterations per chain
  cores = 2                                          # number of cores (could use one per chain)
)

print(fit1, pars = c("transf_drift1", "transf_drift2", "transf_threshold", "transf_ndt"))

# Extract single-trial drift-rates
extract(fit1, paste("drift1_t[", seq(1, dim(sim_data)[1], 1), "]", sep = ''))
