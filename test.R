rm(list = ls())
library(tidyverse)
library(dfoptim)
library(rtdists)
library(rstan)
library(bayesplot)
source('rwald code.r') # from https://osf.io/3sp9t/
source('functions_notebooks.R') # load  functions defined in notebooks


data <- read.csv('data/fontanesi2019.csv')
data <- select(data, -X) # drop pandas index column
head(data)
