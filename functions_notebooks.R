library(tidyverse)
library(dfoptim)
library(rtdists)
library(rstan)
library(bayesplot)
source('rwald code.r') # from https://osf.io/3sp9t/

### DM functions
dm_path <- function(drift, threshold, ndt, rel_sp=.5, noise_constant=1, dt=0.001, max_rt=10) {
  max_tsteps <- max_rt/dt
  
  # initialize the diffusion process
  tstep <- 0
  x <- c(rel_sp*threshold) # vector of accumulated evidence at t=tstep
  time <- c(ndt)
  
  # start accumulating
  while (0 < x[tstep+1] & x[tstep+1] < threshold & tstep < max_tsteps) {
    x <- c(x, x[tstep+1] + rnorm(mean=drift*dt, sd=noise_constant*sqrt(dt), n=1))
    time <- c(time, dt*tstep + ndt)
    tstep <- tstep + 1
  }
  return (data.frame(time=time, dv=x))
}

random_dm <- function(drift, threshold, ndt, rel_sp=.5, noise_constant=1, dt=0.001, max_rt=10) {
  
  if (length(drift) != length(threshold) | length(threshold)!=length(ndt)) {
    stop("drift, threshold, and ndt should have the same length")
  }
  
  n_trials <- length(drift)
  acc <- rep(NA, n_trials)
  rt <- rep(NA, n_trials)
  max_tsteps <- max_rt/dt
  
  # initialize the diffusion process
  tstep <- 0
  x <- rel_sp*threshold # vector of accumulated evidence at t=tstep
  ongoing <- rep(TRUE, n_trials) # have the accumulators reached the bound?
  
  # start accumulating
  while (sum(ongoing) > 0 & tstep < max_tsteps) {
    x[ongoing] <- x[ongoing] + rnorm(mean=drift[ongoing]*dt,
                                     sd=noise_constant*sqrt(dt),
                                     n=sum(ongoing))
    tstep <- tstep + 1
    
    # ended trials
    ended_correct <- (x >= threshold)
    ended_incorrect <- (x <= 0)
    
    # store results and filter out ended trials
    if(sum(ended_correct) > 0) {
      acc[ended_correct & ongoing] <- 1
      rt[ended_correct & ongoing] <- dt*tstep + ndt[ended_correct & ongoing]
      ongoing[ended_correct] <- FALSE
    }
    
    if(sum(ended_incorrect) > 0) {
      acc[ended_incorrect & ongoing] <- 0
      rt[ended_incorrect & ongoing] <- dt*tstep + ndt[ended_incorrect & ongoing]
      ongoing[ended_incorrect] <- FALSE
    }
  }
  return (data.frame(trial=seq(1, n_trials), accuracy=acc, rt=rt))
}

log_likelihood_dm <- function(par, data, ll_threshold=1e-10) {
  
  density <- ddiffusion(rt=data$rt, response=data$response, a=par[1], v=par[2], t0=par[3])
  
  density[density <= ll_threshold] = ll_threshold # put a threhsold on very low likelihoods for computability
  
  return(sum(log(density)))
}

### RDM functions
rdm_path <- function(drift1, drift2, threshold, ndt, sp1=0, sp2=0, noise_constant=1, dt=0.001, max_rt=10) {
  max_tsteps <- max_rt/dt
  
  # initialize the diffusion process
  tstep <- 0
  x1 <- c(sp1*threshold) # vector of accumulated evidence at t=tstep
  x2 <- c(sp2*threshold) # vector of accumulated evidence at t=tstep
  time <- c(ndt)
  
  # start accumulating
  while (x1[tstep+1] < threshold & x2[tstep+1] < threshold & tstep < max_tsteps) {
    x1 <- c(x1, x1[tstep+1] + rnorm(mean=drift1*dt, sd=noise_constant*sqrt(dt), n=1))
    x2 <- c(x2, x2[tstep+1] + rnorm(mean=drift2*dt, sd=noise_constant*sqrt(dt), n=1))
    time <- c(time, dt*tstep + ndt)
    tstep <- tstep + 1
  }
  return (data.frame(time=rep(time, 2), dv=c(x1, x2), accumulator=c(rep(1, length(x1)), rep(2, length(x2)))))
}

random_rdm <- function(n_trials, drift1, drift2, threshold, ndt, sp1=0, sp2=0, noise_constant=1, dt=0.001, max_rt=10) {
  
  choice <- rep(NA, n_trials)
  rt <- rep(NA, n_trials)
  max_tsteps <- max_rt/dt
  
  # initialize the diffusion process
  tstep <- 0
  x1 <- rep(sp1*threshold, n_trials) # vector of accumulated evidence at t=tstep
  x2 <- rep(sp2*threshold, n_trials) # vector of accumulated evidence at t=tstep
  ongoing <- rep(TRUE, n_trials) # have the accumulators reached the bound?
  
  # start accumulating
  while (sum(ongoing) > 0 & tstep < max_tsteps) {
    x1[ongoing] <- x1[ongoing] + rnorm(mean=drift1*dt, sd=noise_constant*sqrt(dt), n=sum(ongoing))
    x2[ongoing] <- x2[ongoing] + rnorm(mean=drift2*dt, sd=noise_constant*sqrt(dt), n=sum(ongoing))
    tstep <- tstep + 1
    
    # ended trials
    ended1 <- (x1 >= threshold)
    ended2 <- (x2 >= threshold)
    
    # store results and filter out ended trials
    if(sum(ended1) > 0) {
      choice[ended1 & ongoing] <- 1
      rt[ended1 & ongoing] <- dt*tstep + ndt
      ongoing[ended1] <- FALSE
    }
    
    if(sum(ended2) > 0) {
      choice[ended2 & ongoing] <- 2
      rt[ended2 & ongoing] <- dt*tstep + ndt
      ongoing[ended2] <- FALSE
    }
  }
  return (data.frame(trial=seq(1, n_trials), choice=choice, rt=rt))
}

log_likelihood_rdm <- function(par, data, ll_threshold=1e-10) {
  # par order: drift1, drift2, threshold, ndt
  density <- rep(NA, dim(data)[1])
  data$RT <- (data$RT - par[4]) # shift the distribution by the NDT
  # dWald: density for single accumulator
  # pWald: cumulative density for single accumulator
  density[data$R==1] <- dWald(data$RT[data$R==1], v=par[1], B=par[3], A=0)*(1 - pWald(data$RT[data$R==1], v=par[2], B=par[3], A=0))
  density[data$R==2] <- dWald(data$RT[data$R==2], v=par[2], B=par[3], A=0)*(1 - pWald(data$RT[data$R==2], v=par[1], B=par[3], A=0))
  
  density[density <= ll_threshold] = ll_threshold # put a threhsold on very low likelihoods for computability
  
  return(sum(log(density)))
}

### RLDM functions
generate_RL_stimuli <- function (n_trials, trial_types, mean_options, sd_options) {
  # n_trials_block : Number of trials per learning session.
  # n_blocks : Number of learning session per participant.
  # n_participants : Number of participants.
  
  # trial_types : List containing possible pairs of options.
  # E.g., in the original experiment: c('1-2', '1-3', '2-4', '3-4').
  # It is important that they are separated by a '-',
  # and that they are numbered from 1 to n_options (4 in the example).
  # Also, the "incorrect" option of the couple should go first in each pair.
  
  # mean_options : Mean reward for each option.
  # The length should correspond to n_options.
  
  # sd_options : SD reward for each option.
  # The length should correspond to n_options.
  
  generate_trial_type_sequence <- function (n_trials, trial_types) {
    n_trial_types = length(trial_types)
    sequence = rep(trial_types, n_trials/n_trial_types)
    sample(sequence)
    
    count = 3
    while (count < length(sequence)) {
      if (sequence[count]==sequence[count-1] & sequence[count-1]==sequence[count-2]) {
        sample(sequence)
        count = 2
      }
      count = count + 1
    }
    return(sequence)
  }
  
  task_design = data.frame(
    trial = seq(1, n_trials),
    trial_type = generate_trial_type_sequence(n_trials, trial_types)
  )
  
  task_design <- separate(task_design, 
                          col=trial_type,
                          into=c('inc_option', 'cor_option'),
                          sep='-',
                          remove=FALSE,
                          convert=TRUE)
  
  options = unique(c(unique(task_design$inc_option), unique(task_design$cor_option)))
  n_options = length(options)
  print("The task will be created with the following options:")
  print(data.frame(options=options, mean=mean_options, sd=sd_options))
  
  task_design$f_inc = round(rnorm(
    mean = mean_options[task_design$inc_option], 
    sd = sd_options[task_design$inc_option], 
    n=dim(task_design)[1]))
  task_design$f_cor = round(rnorm(
    mean = mean_options[task_design$cor_option], 
    sd = sd_options[task_design$cor_option], 
    n=dim(task_design)[1]))
  
  return(task_design)
}

# simulate data from simple RLDM
random_rldm <- function(stimuli, alpha, drift_scaling, threshold, ndt, rel_sp=.5, noise_constant=1, dt=0.001, max_rt=10, initial_value_learning=20) {
  n_trials = dim(stimuli)[1]
  options = unique(c(unique(stimuli$inc_option), unique(stimuli$cor_option)))
  n_options = length(options)
  
  # Rescorla Wagner learning rule
  stimuli$Q_cor = NA
  stimuli$Q_inc = NA
  Q = rep(initial_value_learning, n_options)
  
  for (n in seq(1, n_trials, 1)) {
    index_cor = stimuli[n, "cor_option"]
    index_inc = stimuli[n, "inc_option"]
    
    # current expectations
    stimuli[n, "Q_cor"] = Q[index_cor]
    stimuli[n, "Q_inc"] = Q[index_inc]
    
    # update expectations
    Q[index_cor] = Q[index_cor] + alpha*(stimuli[n, "f_cor"] - Q[index_cor])
    Q[index_inc] = Q[index_inc] + alpha*(stimuli[n, "f_inc"] - Q[index_inc])
  }
  
  # Diffusion Model as decision rule
  stimuli$drift = drift_scaling*(stimuli$Q_cor - stimuli$Q_inc)
  stimuli$threshold = threshold
  stimuli$ndt = ndt
  
  performance = random_dm(stimuli$drift, stimuli$threshold, stimuli$ndt)
  
  return(cbind(stimuli, performance[,c("accuracy", "rt")]))
}

log_likelihood_rldm <- function(par, data, initial_value_learning=20, ll_threshold=1e-10) {
  # par order: alpha, drift_scaling, threshold, ndt
  
  n_trials = dim(data)[1]
  options = unique(c(unique(data$inc_option), unique(data$cor_option)))
  n_options = length(options)
  
  # Rescorla Wagner learning rule
  data$Q_cor = NA
  data$Q_inc = NA
  Q = rep(initial_value_learning, n_options)
  
  for (n in seq(1, n_trials, 1)) {
    index_cor = data[n, "cor_option"]
    index_inc = data[n, "inc_option"]
    
    # current expectations
    data[n, "Q_cor"] = Q[index_cor]
    data[n, "Q_inc"] = Q[index_inc]
    
    # update expectations
    Q[index_cor] = Q[index_cor] + par[1]*(data[n, "f_cor"] - Q[index_cor])
    Q[index_inc] = Q[index_inc] + par[1]*(data[n, "f_inc"] - Q[index_inc])
  }
  
  # Diffusion Model as decision rule
  drift <- par[2]*(data$Q_cor - data$Q_inc)
  density <- ddiffusion(rt=data$rt, response=data$response, a=par[3], v=drift, t0=par[4])
  
  density[density <= ll_threshold] = ll_threshold # put a threhsold on very low likelihoods for computability
  
  return(sum(log(density)))
}


