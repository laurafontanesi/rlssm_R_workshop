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

# fit hierarchical DM on real data
data <- read.csv('data/fontanesi2019.csv')
data <- select(data, -X) # drop pandas index column
data <- data[data$participant < 6,] # select only 5 participants
data$accuracy_recoded <- data$accuracy
data[data$accuracy==0, "accuracy_recoded"] <- -1

data[(data$inc_option == 1)&(data$cor_option == 2),"condition_label"] <- "difficult_low"
data[(data$inc_option == 1)&(data$cor_option == 3),"condition_label"] <- "easy_low"
data[(data$inc_option == 2)&(data$cor_option == 4),"condition_label"] <- "easy_high"
data[(data$inc_option == 3)&(data$cor_option == 4),"condition_label"] <- "difficult_high"

data[(data$inc_option == 1)&(data$cor_option == 2),"condition"] <- 1
data[(data$inc_option == 1)&(data$cor_option == 3),"condition"] <- 2
data[(data$inc_option == 2)&(data$cor_option == 4),"condition"] <- 3
data[(data$inc_option == 3)&(data$cor_option == 4),"condition"] <- 4

ggplot(data = data, aes(x = condition, y = accuracy, color=condition_label))+
  stat_summary(fun = "mean", geom="point",  size=5) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",  size=.5, width=.3)

ggplot(data = data, aes(x = condition, y = rt, color=condition_label))+
  stat_summary(fun = "mean", geom="point",  size=5) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",  size=.5, width=.3)

head(data)

unique(data$participant)

sim_data_for_stan <- list(
  N = dim(data)[1],
  L = length(unique(data$participant)),
  C = 4,
  condition = data$condition,
  participant = data$participant,
  accuracy = data$accuracy_recoded,
  rt = data$rt,
  starting_point = 0.5
)

fit1 <- stan(
  file = "stan_models/hierDDM_cond.stan",            # Stan program
  data = sim_data_for_stan,                          # named list of data
  chains = 2,                                        # number of Markov chains
  warmup = 1000,                                     # number of warmup iterations per chain
  iter = 3000,                                       # total number of iterations per chain
  cores = 2                                          # number of cores (could use one per chain)
)

print(fit1, pars = c("mu_drift[1]", "mu_drift[2]", "mu_drift[3]", "mu_drift[4]",
                     "sd_drift[1]", "sd_drift[2]", "sd_drift[3]", "sd_drift[4]",
                     "mu_threshold[1]", "mu_threshold[2]", "mu_threshold[3]", "mu_threshold[4]",
                     "sd_threshold[1]", "sd_threshold[2]", "sd_threshold[3]", "sd_threshold[4]",
                     "mu_ndt", "sd_ndt"))

posterior <- as.matrix(fit1)

plot_title <- ggtitle("Posterior distributions mu drift-rate per condition",
                      "with medians and 95% intervals")
mcmc_areas(posterior,
           pars = c("mu_drift[1]", "mu_drift[2]", "mu_drift[3]", "mu_drift[4]"),
           prob = 0.95) + plot_title

plot_title <- ggtitle("Posterior distributions mu threshold per condition",
                      "with medians and 95% intervals")

mcmc_areas(posterior,
           pars = c("mu_threshold[1]", "mu_threshold[2]", "mu_threshold[3]", "mu_threshold[4]"),
           prob = 0.95) + plot_title

plot_title <- ggtitle("Posterior distributions individual NDT parameters",
                      "with medians and 95% intervals")

mcmc_areas(posterior,
           pars = c("ndt_sbj[1]", "ndt_sbj[2]", "ndt_sbj[3]", "ndt_sbj[4]"),
           prob = 0.95) + plot_title
