library(deSolve) # package for solving differential equations.
library(minpack.lm) # package for non-linear regression #rks to laz: I think we probably can just use nls() in stats, not nlsLM from linpack. Please check in R documentation if there is any functional difference at all
library(mvtnorm) # package for calculating confidence intervals by Monte Carlo simulation based on variance-covariance matrices #rks to laz: I added to comment.Please check that my addition is OK.
library(Hmisc)
library(dplyr)
rm(list=ls())
#Create dataframes that store the fibroblast WGE simple CA data used in 16Cacao 

#Oxygen = data.frame(d = c(0, .0125, .02, .025, .05, .075, .1, .2, .4),CA = c(.24, 1.66, 2.43, 2.37, 1.16, 2.85, 2.58, 6.94, 6.91), ion = "O", Z =8, L = 75, Z.b = 595) #remove 0 doses right from the start RKS 7/15/2018
# Some GCR components are high-speed Oxygen nuclei that are almost fully ionized. d=dose; CA are per hundred cells.
#putting it in one big data frame. #rks: the data frame incorporates a correction to Fe600 at dose 0.06, near line 34

modified_df = read.csv(file = 'CSV.csv')
