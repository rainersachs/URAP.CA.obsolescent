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

main_df = read.csv(file = 'CSV.csv')

#Comments:
comment = c(
  "There are some LET discrepancies. We decided to use the values given by Hada in all cases",
  "For now, we changed the LET measurements for two rows: Row28 from 17 to 20.93, Row91 from 190 to 175",
  "For now, we changed the Z.b value for row42 from 504.9148 to 504.918"
)


#Column Descriptions
KE_label = "Specific kinetic energy (SKE in MeV/u)"
CA_label = "WGE apparently simple aberrations"
X.T_label = "Total apparently simple (AS) abberrations observed with 3-color painting"
X_label = "Exchanges"
X.1_label = "Incomplete Exchanges"
X.2_label = "Dicentrics"
X.3_label = "Incomplete Dicentrics"
nn_label = "Number of nucleons (protons + neutrons) in the atomic nucleus"

#Output text file
writeLines(c(
  "Comments:",
  "",
  comment,
  "",
  "###################################################################################################################",
  "",
  "Column Descriptions:",
  "",
  paste("KE:",KE_label),
  paste("CA:", CA_label),
  paste("X.T:", X.T_label),
  paste("X:", X_label),
  paste("X.1:", X.1_label),
  paste("X.2:", X.2_label),
  paste("X.3:", X.3_label),
  paste("nn:", nn_label)
), "DataDescription.txt")







#Setting attributes
attr(main_df$KE, 'label') <- KE_label
attr(main_df$CA, 'label') <- CA_label
attr(main_df$X.T, 'label') <- X.T_label
attr(main_df$X, 'label') <- X_label 
attr(main_df$X.1, 'label') <- X.1_label
attr(main_df$X.2, 'label') <- X.2_label
attr(main_df$X.3, 'label') <- X.3_label
attr(main_df$nn, 'label') <- nn_label
comment(df) = comment



# Calculate the background

dose0 = main_df[main_df$d == 0.0000,]
num = sum(dose0$X.T)
risk = sum(dose0$At.Risk)
p = num/risk
BG_CA = p
BG_error = (sqrt((p*(1-p))/(risk)))
