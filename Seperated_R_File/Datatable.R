library(deSolve) # package for solving differential equations.
library(minpack.lm) # package for non-linear regression #rks to laz: I think we probably can just use nls() in stats, not nlsLM from linpack. Please check in R documentation if there is any functional difference at all
library(mvtnorm) # package for calculating confidence intervals by Monte Carlo simulation based on variance-covariance matrices #rks to laz: I added to comment.Please check that my addition is OK.
library(Hmisc)
library(dplyr)
rm(list=ls())
#Create dataframes that store the fibroblast WGE simple CA data used in 16Cacao 

#Oxygen = data.frame(d = c(0, .0125, .02, .025, .05, .075, .1, .2, .4),CA = c(.24, 1.66, 2.43, 2.37, 1.16, 2.85, 2.58, 6.94, 6.91), ion = "O", Z =8, L = 75, Z.b = 595) #remove 0 doses right from the start RKS 7/15/2018
# Some GCR components are high-speed Oxygen nuclei that are almost fully ionized. d=dose; CA are per hundred cells.
Oxygen = read.csv(file = '/Users/zhaoliyang/Desktop/URAP\ FALL2018/Oxygen.csv')[,-1]
Si = read.csv(file = '/Users/zhaoliyang/Desktop/URAP\ FALL2018/Si.csv')[,-1]
Ti = read.csv(file = '/Users/zhaoliyang/Desktop/URAP\ FALL2018/Ti.csv')[,-1]
Fe600 = read.csv(file = '/Users/zhaoliyang/Desktop/URAP\ FALL2018/Fe600.csv')[,-1]
#600 refers to the energy in MeV per atomic mass unit in this Iron beam
Fe450 = read.csv(file = '/Users/zhaoliyang/Desktop/URAP\ FALL2018/Fe450.csv')[,-1]
Fe300 = read.csv(file = '/Users/zhaoliyang/Desktop/URAP\ FALL2018/Fe300.csv')[,-1]

#putting it in one big data frame. #rks: the data frame incorporates a correction to Fe600 at dose 0.06, near line 34
modified_df = rbind(Oxygen, Si, Ti, Fe600, Fe450, Fe300)

#Next modify the data frame to get rid of the zero dose points. Background CA frequency was determined seperately.
#modified_df = big_df[big_df$d != 0, ]
modified_df$CA = modified_df$CA*0.01 
modified_df$error = modified_df$error*0.01
modified_df$errorbar_lower = modified_df$CA - modified_df$error
modified_df$errorbar_upper = modified_df$CA + modified_df$error
