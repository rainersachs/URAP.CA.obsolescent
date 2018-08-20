library(deSolve) # package for solving differential equations
library(minpack.lm) # package for non-linear regression #rks to laz: I think we probably can just use nls() in stats, not nlsLM from linpack. Please check in R documentation if there is any functional difference at all
library(mvtnorm) # package for calculating confidence intervals by Monte Carlo simulation based on variance-covariance matrices #rks to laz: I added to comment.Please check that my addition is OK.
library(Hmisc)
rm(list=ls())

#Create dataframes that store the fibroblast WGE simple CA data used in 16Cacao 

#Oxygen = data.frame(d = c(0, .0125, .02, .025, .05, .075, .1, .2, .4),CA = c(.24, 1.66, 2.43, 2.37, 1.16, 2.85, 2.58, 6.94, 6.91), ion = "O", Z =8, L = 75, Z.b = 595) #remove 0 doses right from the start RKS 7/15/2018
# Some GCR components are high-speed Oxygen nuclei that are almost fully ionized. d=dose; CA are per hundred cells.
Oxygen = data.frame(d = c( .0125, .02, .025, .05, .075, .1, .2, .4),CA = c( 1.66, 2.43, 2.37, 1.16, 2.85, 2.58, 6.94, 6.91),error= c(0.63, 0.77, 0.75, 0.52, 0.82, 0.78, 1.31, 1.59), ion = "O", Z =8, L = 75, Z.b = 595)
Si = data.frame(d = c( .02, .04, .06, .08, .1, .12, .2, .4, .8, 1.2), CA = c(1.26, 1.14, 1.58, 1.22, 1.89, 3.47, 4.6, 9.79, 27.01, 38.84), error =c(0.05, 0.07, 0.56, 0.18, 0.60, 1.23, 1.60, 1.55, 4.27, 7.21), ion = "Si", Z = 14, L = 100, Z.b = 690)
Ti = data.frame(d = c(0.02, .04, .06, .08, .1, .15, .3, .6), CA = c(1.99, 1.88, 1.44, 2.67, 2.57, 2.50, 5.64, 11.19), error = c(0.70, 0.66, 0.59, 0.80, 0.78, 0.48, 1.15, 2.39), ion = "Ti", Z = 22, L = 125, Z.b = 770)
Fe600 = data.frame(d = c(.01, .02, .04, .06, .08, .1, .12, .2, .4, .8), CA = c( .76, .99, 1.2, 1.74, 1.28, 1.2, 1.7, 3.02, 5.52, 12.42), error = c(0.38, 0.24, 0.21, 0.43, 0.37, 0.54, 0.17, 0.55, 1.75, 2.59),ion = "Fe600", Z = 26, L = 175, Z.b = 1075) 
#600 refers to the energy in MeV per atomic mass unit in this Iron beam
Fe450 = data.frame(d = c(.02, .04, .06, .08, .1, .2, .4), CA = c(.86, .6, .8, 1.22, 2.02, 2.3, 4.77), error = c(0.43, 0.34, 0.40, 0.50, 0.64, 0.73, 1.09), ion = "Fe450", Z = 26, L = 195, Z.b = 1245)
Fe300 = data.frame(d = c(.005, .01,  0.02, .04, .07, .1, .2, .4, .8), CA = c( 1.23, 1.47, 1.22, .97, 1.46, 1.21, 4.38, 6.22, 13.6), error = c(0.55, 0.60, 0.55, 0.49, 0.60, 0.54, 1.03, 1.22, 3.62), ion = "Fe300", Z = 26, L = 240, Z.b = 1585)

#putting it in one big data frame. #rks: the data frame incorporates a correction to Fe600 at dose 0.06, near line 34
modified_df = rbind(Oxygen, Si, Ti, Fe600, Fe450, Fe300)

#Next modify the data frame to get rid of the zero dose points. Background CA frequency was determined seperately.
#modified_df = big_df[big_df$d != 0, ]
modified_df$CA = modified_df$CA*0.01 
modified_df$error = modified_df$error*0.01
modified_df$errorbar_lower = modified_df$CA - modified_df$error
modified_df$errorbar_upper = modified_df$CA + modified_df$error

#Calibrating parameters and Model Selection (3 Models: NTE1, NTE2 and IDER)
NTE1_function = function(d, L, Z.b, eta0 = 0.00011, eta1 = 0.007, sig0 = 6.12, kap = 796) {
  0.0017 + eta0*L*exp(-eta1*L)*(d != 0) + 
    (6.242*(d/L))*(sig0*(1-exp(-Z.b/kap))^2 + 0.041/6.24*L*(1 - (1-exp(-Z.b/kap))^2))
} 
NTE2_function = function(d, L, Z.b, eta0 = 0.00047, eta1 = 0.011, sig0 = 6.75, kap = 590) {
  0.0017 + eta0*L*exp(-eta1*L)*exp(-(1012*(d/L)))*(d != 0) + 
    (6.242*(d/L))*(1-exp(-(1012*(d/L))))*
    (sig0*(1-exp(-Z.b/kap))^2 + 0.041/6.24*L*(1 - (1-exp(-Z.b/kap))^2))
}
#Our IDERs (Individual Dose Effect Relations). Applicable to the 1-ion components of a mixed simulated GCR beam 
#Modifying NTE1 and NTE2 by insisting they be twice continuously differentiable and monotonic increasing. Double check NTE1, NTE2, Our model
IDER = function(d, L, Z.b, eta00, eta10, sig00, kap0) {
  P = (1-exp(-Z.b/kap0))^2
  sig = sig00*P + 0.041/6.24*L*(1-P) # 0.041 +- 0.0051 comes from 16Cacao
  eta = eta00*L*exp(-eta10*L)
  0.00071 + sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d))  #0.00071 + sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^3*d))#don't use
} 
Parsimonious=function(d,L,eta_0,sigm_0){
  0.00071 + sigm_0*6.24*d/L*(1-exp(-1024*d/L)) + eta_0*(1-exp(-10^5*d))  
}
#nls (non-linear least square) method to get the parameters needed (4 parameter estimation) #rks to laz: see note under issues
IDER_model = nlsLM(CA ~ IDER(d, L, Z.b, eta00, eta10, sig00, kap0), data = modified_df, start = list(eta00 = 0.001, eta10 = 0.01, sig00 = 5, kap0 = 500), 
                   weights = (1/(modified_df$error)^2))
#Getting coefficients of the IDER model from nlsLM
coefs <- coef(IDER_model)
eta0_hat <- as.numeric(coefs[1])
eta1_hat <- as.numeric(coefs[2])
sig0_hat <- as.numeric(coefs[3])
kap_hat <- as.numeric(coefs[4])

Parsimonious_model = nlsLM(CA ~Parsimonious (d,L,eta_0,sigm_0) , data = modified_df, start = list(eta_0 = 0.001, sigm_0 = 5), 
                  weights = (1/(modified_df$error)^2))
ccoef = coef(Parsimonious_model)
eta_0_hat = as.numeric(ccoef[1])
sigm_0_hat = as.numeric(ccoef[2])

L_function = function(func, eta0, eta1, sig0, kap) { #A function trhat returns the residual squared
  a = vector(length = 0)
  for (i in 1:length(modified_df[, 1])) {
    a = c(a, modified_df$CA[i] - func(d = modified_df$d[i], L = modified_df$L[i], Z.b = modified_df$Z.b[i], eta0 = eta0, eta1 = eta1, sig0 = sig0, kap = kap))
  }
  return(a^2)
}
L_NTE1 = L_function(NTE1_function, eta0 = 0.00011, eta1 = 0.007, sig0 = 6.12, kap = 796)
L_NTE2 = L_function(NTE2_function, eta0 = 0.00047, eta1 = 0.011, sig0 = 6.75, kap = 590)
L_IDER = L_function(IDER, eta0 = eta0_hat, eta1 = eta1_hat, sig0 = sig0_hat, kap = kap_hat)
WRSS_NTE1 = sum((1/modified_df$error^2)*L_NTE1)
WRSS_NTE2 = sum((1/modified_df$error^2)*L_NTE2)
WRSS_IDER = sum((1/modified_df$error^2)*L_IDER)

#functions for AIC and BIC calculation for Weighted Least Square regression (using WRSS calculated above)
AIC_function = function(RSS, k = 4, n = length(modified_df[ , 1])) {
  n + n*log(2*pi) + n*log(RSS/n) + 2*(k+1)
}

BIC_function = function(n = length(modified_df[, 1]), k = 4, RSS) {
  n + n*log(2*pi) + n*log(RSS/n) + log(n)*(k+1)
}

NTE1_AIC = AIC_function(RSS = WRSS_NTE1)
NTE2_AIC = AIC_function(RSS = WRSS_NTE2)
IDER_AIC = AIC_function(RSS = WRSS_IDER)
NTE1_BIC = BIC_function(RSS = WRSS_NTE1)
NTE2_BIC = BIC_function(RSS = WRSS_NTE2)
IDER_BIC = BIC_function(RSS = WRSS_IDER)
information_critera_df = data.frame(AIC = c(NTE1_AIC, NTE2_AIC, IDER_AIC), BIC = c(NTE1_BIC, NTE2_BIC, IDER_BIC), row.names = c("NTE1 model", "NTE2 model", "IDER model"))
information_critera_df #IDER has the lowest score (performs better) in both criteria.


##MIXDER_function and graphs
MIXDER_function = function(r, L, Z.b, d = seq(0, 0.2, by = 0.001), eta0 = eta0_hat, eta1 = eta1_hat, sig0 = sig0_hat ,kap = kap_hat) {
  dE=function(yini,State,Pars){
    eta0 = eta0; eta1 = eta1; sig0 = sig0; kap = kap
    with(as.list(c(State, Pars)), {
      P = vector(length = length(L))
      sig = vector(length = length(L))
      etaa = vector(length = length(L))
      u = vector(length = length(L))
      for (i in 1:length(L)) {
        P[i] = (1-exp(-Z.b[i]/kap))^2
        sig[i] = sig0*P[i] + 0.041/6.24*L[i]*(1-P[i])
        etaa[i] = eta0*L[i]*exp(-eta1*L[i])
        u[i] = uniroot(function(d) sig[i]*6.24*d/L[i]*(1-exp(-1024*d/L[i])) + etaa[i]*(1-exp(-10^5*d)) - I, lower = 0, upper = 1, extendInt = "yes", tol = 10^-10)$root
      }
      dI = vector(length = length(L))
      for (i in 1:length(L)) {
        dI[i] = r[i]*(sig[i]*6.24/L[i]*exp(-1024*u[i]/L[i])*(exp(1024*u[i]/L[i]) + 1024*u[i]/L[i] - 1) + etaa[i]*10^5*exp(-10^5*u[i]))
      }
      dI = sum(dI)
      return(list(c(dI)))
    })
  }
  pars = NULL; yini = c(I= 0); d = d
  out = ode(yini,times = d, dE, pars, method = "radau")
  return(out)
}

#====== an example of a plot: 6 HZE, each 1/6 of total dose ======#
#Simple Effect Additivity (black line)
SEA=function(d){
  IDER(d/6,L[1],Z.b[1])+IDER(d/6,L[2],Z.b[2])+IDER(d/6,L[3],Z.b[3])+IDER(d/6,L[4],Z.b[4])+IDER(d/6,L[5],Z.b[5])+IDER(d/6,L[6],Z.b[6])
}

r=rep(1/6,6);L = c(75, 100, 125, 175, 195, 240); Z.b = c(595, 690, 770, 1075, 1245, 1585); dose=seq(0, 0.4, by=0.001)

MX=MIXDER_function(r, L, Z.b, d=dose) #for this mixture can't go much above 0.6#

IDER = function(d, L, Z.b, eta0 = eta0_hat, eta1 = eta1_hat, sig0 = sig0_hat, kap = kap_hat) {
  P = (1-exp(-Z.b/kap_hat))^2
  sig = sig0*P + 0.041/6.24*L*(1-P)
  eta = eta0*L*exp(-eta1*L)
  sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d))
}
#Plotting
plot(MX[,1],MX[,2],type='l',bty='l',col='red',ann='F',ylim=c(0,.10), main="CA vs Dose", xlab="Dose", ylab="CA")
for (ii in 1:length(L)){lines(dose, IDER(dose,L[ii],Z.b[ii]),col='green')}
lines(dose,SEA(dose),lty=2)
# ===== end plot example of 6 HZE sharing dose equally ====== #

# ===== example of using errbarr in a plot ========#
#========== Figure for Chang's new proton data point 5/22/2018=============#
ddose <- c(.00001*0:9,.0001*1:9,.001*1:9,.005*1:80)
plot(ddose,Parsimonious(ddose,L=75,eta_0_hat,sigm_0_hat),
     type='l',lwd=2, ann=FALSE, xlim=c(0,.4),ylim=c(0,.12), bty='u',col="black")#need model!
     errbar(modified_df[1:8, "d"], modified_df[1:8, "CA"],yplus=modified_df[1:8, "CA"]+1.96*modified_df[1:8, "error"],
       yminus=modified_df[1:8, "CA"]-1.96*modified_df[1:8, "error"], pch = 19,cap=0.03, add=TRUE, col='black',
       errbar.col = 'black',lwd=2) #  RKS: oxygen data points 95%CI
#======= end errorbar example ========#

#Confidence Intervals (Monte Carlo)
set.seed(19970101)
MM <- 500
#Sample parameters from their distributions with covariances
sig = vcov(IDER_model)
monte_carlo_parameters = rmvnorm(n = MM, mean = c(eta0 = eta0_hat, eta1 = eta1_hat, sig0 = sig0_hat, kap = kap_hat), sigma = sig)
eta0_MC = monte_carlo_parameters[, 1]
eta1_MC = monte_carlo_parameters[, 2]
sig0_MC = monte_carlo_parameters[, 3]
kap_MC = monte_carlo_parameters[, 4]
kap_MC[kap_MC <= 1e-6] = 1e-5 #10 of them were fixed

#Sample parameters from their distributions without covariances
eta0_MC_var = rnorm(MM, mean = eta0_hat, sd = 4.631e-05)
eta1_MC_var = rnorm(MM, mean = eta1_hat, sd = 1.617e-04)
sig0_MC_var = rnorm(MM, mean = sig0_hat, sd = 1.966e+00)
kap_MC_var = rnorm(MM, mean = kap_hat, sd = 2.628e+02)
kap_MC_var[kap_MC_var <= 1e-6] = 1e-5 

# #Monte Carlo functions (without and with covariance)
# CI_function_MIXDER_var = function(d, d_interested, r, interval = 0.95, L, Z.beta) { #Outputs the CI without using the vcov matrix of the parameters
#   MIXDER_curve = list(0)
#   for (i in 1:MM) {
#     MIXDER_curve[[i]] = MIXDER_function(r = r, d = d, L = L, Z.b = Z.b, eta0 = eta0_MC_var[i], eta1 = eta1_MC_var[i], sig0 = sig0_MC_var[i], kap = kap_MC_var[i])
#   }
#   info = vector(length = 0)
#   for (i in 1:MM) {
#     info = c(info, MIXDER_curve[[i]][, 2][2])
#   }
#   info = sort(info)
#   lower_bound = info[(1-interval)/2*MM]
#   upper_bound = info[(interval + (1-interval)/2)*MM]
#   CI = c(lower_bound, upper_bound)
#   return(CI)
# }

################### # EGH 08.20.18 REFACTORING

# Monte Carlo functions (without and with covariance)
refactor_CI_function_MIXDER_var = function(d, d_interested, r, interval = 0.95, L, Z.beta) { #Outputs the CI without using the vcov matrix of the parameters
  MIXDER_curve = list(0)
  for (i in 1:MM) {
    MIXDER_curve[[i]] = MIXDER_function(r = r, d = d, L = L, Z.b = Z.b, eta0 = eta0_MC_var[i], eta1 = eta1_MC_var[i], sig0 = sig0_MC_var[i], kap = kap_MC_var[i])
    cat(paste("  Currently at Monte Carlo step:", toString(i), "of", 
              toString(MM)), sprintf('\r'))
  }
  return(MIXDER_curve)
}

################### END REFACTORING

start_time <- Sys.time()
#MIXDER of 2-ion (Silicon and Fe600) with 60 different dosage points and Monte Carlo Ribbon Plot using covariances
d1 = c(seq(0, 0.009, 0.001), seq(0.01, 0.5, by = 0.01))
out = MIXDER_function(r = rep(1/2, times = 2), L = c(100, 175), Z.b = c(690, 1075), d = d1, eta0 = eta0_hat, eta1 = eta1_hat, sig0 = sig0_hat, kap = kap_hat)
two_ion_MIXDER = data.frame(d = out[, 1], CA = out[, 2] + 0.00071)
ninty_five_CI_lower = vector(length = length(d1))
ninty_five_CI_upper = vector(length = length(d1))

DERs <- matrix(nrow = MM, ncol = length(d1))
# for (i in 1:length(d1)){ # EGH 07.22.18 I strongly believe this is the area of the slowdown; we are calculating 29500 more curves then necessary.
  # info = vector(length = 0)
for (j in 1:MM) {
  DERs[j, ] <- MIXDER_function(r = rep(1/2, times = 2), 
                                 d = two_ion_MIXDER$d[1:length(d1)], L = c(100, 175), 
                                 Z.b = c(690, 1075), eta0 = eta0_MC[j], 
                                 eta1 = eta1_MC[j], sig0 = sig0_MC[j], 
                                 kap = kap_MC[j])[, 2]
  cat(paste("  Currently at Monte Carlo step:", toString(j), "of", 
            toString(MM)), sprintf('\r'))
}
for (i in 1:length(d1)) {
  sample_values <- sort(DERs[, i])
  # Returning resulting CI
  ninty_five_CI_lower[i] <- sample_values[ceiling((1 - 0.95) / 2 * MM)]
  ninty_five_CI_upper[i] <- sample_values[(0.95 + (1 - 0.95) / 2) * MM]
  # ninty_five_CI_lower[i] = quantile(info,(1 - 0.95) / 2) + 0.00071
  # ninty_five_CI_upper[i] = quantile(info,1 - (1 - 0.95) / 2) + 0.00071
}
Sys.time() - start_time

two_ion_MIXDER$CI_lower = ninty_five_CI_lower
two_ion_MIXDER$CI_upper = ninty_five_CI_upper
two_ion_MIXDER$simpleeffect = IDER(d = 0.5*d1, L = 100, Z.b = 690) + IDER(d = 0.5*d1, L = 175, Z.b = 1075) + 0.00071
two_ion_MIXDER$silicon = IDER(d = d1, L = 100, Z.b = 690) + 0.00071
two_ion_MIXDER$ironsix = IDER(d = d1, L = 175, Z.b = 1075) + 0.00071
#The graphing part
d1 <- two_ion_MIXDER$d
CA1 <- two_ion_MIXDER$CA
simpleeffect1 <- two_ion_MIXDER$simpleeffect
silicon1 <- two_ion_MIXDER$silicon
ironsix1 <- two_ion_MIXDER$ironsix
plot(x = d1 * 100, y = CA1 * 100, type = "l", col = "red")
lines(x = d1 * 100, y = simpleeffect1 * 100, col = "black", lty = 2, lwd = 0.5)
lines(x = d1 * 100, y = silicon1* 100, col = "green")  
lines(x = d1 * 100, y = ironsix1* 100, col = "green")
lines(x= d1*100 , y = two_ion_MIXDER$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d1*100 , y = two_ion_MIXDER$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d1*100,rev(d1*100)),c(two_ion_MIXDER$CI_lower * 100, rev(two_ion_MIXDER$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)

#2-ion 60 doses without using covariances
# two_ion_MIXDER_var = data.frame(d = out[, 1], CA = out[, 2] + 0.00071)
# ninty_five_CI_lower = vector(length = 0)
# ninty_five_CI_upper = vector(length = 0)
# a = vector(length = 0)
# for (i in 2:length(d1)) {
#   print(d1[i])
#   a = CI_function_MIXDER_var(d = c(0, two_ion_MIXDER_var$d[i]), r = rep(1/2, times = 2), L = c(100, 175), Z.b = c(690, 1075))
#   ninty_five_CI_lower = c(ninty_five_CI_lower, a[1])
#   ninty_five_CI_upper = c(ninty_five_CI_upper, a[2])
# }

################## 08.20.18 EGH REFACTORING
start_time <- Sys.time()
curves <- refactor_CI_function_MIXDER_var(d = d1, r = rep(1/2, times = 2), 
                                          L = c(100, 175), Z.b = c(690, 1075))
ninty_five_CI_lower = vector(length = 0)
ninty_five_CI_upper = vector(length = 0)
for (i in 2:length(d1)) {
  info = vector(length = 0)
  for (j in 1:MM) {
    info = c(info, curves[[j]][, 2][i])
  }
  info = sort(info)
  lower = info[(1 - 0.95) / 2 * MM]
  upper = info[(0.95 + (1 - 0.95) / 2) * MM]
  ninty_five_CI_lower = c(ninty_five_CI_lower, lower)
  ninty_five_CI_upper = c(ninty_five_CI_upper, upper)
  cat(paste("  Currently at Monte Carlo step:", toString(i), "of", 
            toString(length(d1))), sprintf('\r'))
}
Sys.time() - start_time
################## END REFACTORING

two_ion_MIXDER_var$CI_lower = c(0, ninty_five_CI_lower) + 0.00071
two_ion_MIXDER_var$CI_upper = c(0, ninty_five_CI_upper) + 0.00071
two_ion_MIXDER_var$simpleeffect = IDER(d = 0.5*d1, L = 100, Z.b = 690) + IDER(d = 0.5*d1, L = 175, Z.b = 1075) + 0.00071
two_ion_MIXDER_var$silicon = IDER(d = d1, L = 100, Z.b = 690) + 0.00071
two_ion_MIXDER_var$ironsix = IDER(d = d1, L = 175, Z.b = 1075) + 0.00071

d1 <- two_ion_MIXDER_var$d
CA1 <- two_ion_MIXDER_var$CA
simpleeffect1 <- two_ion_MIXDER_var$simpleeffect
silicon1 <- two_ion_MIXDER_var$silicon
ironsix1 <- two_ion_MIXDER_var$ironsix
plot(x = d1 * 100, y = CA1 * 100, type = "l", col = "red")
lines(x = d1 * 100, y = simpleeffect1 * 100, col = "black", lty = 2, lwd = 0.5)
lines(x = d1 * 100, y = silicon1* 100, col = "green")  
lines(x = d1 * 100, y = ironsix1* 100, col = "green")
lines(x= d1*100 , y = two_ion_MIXDER_var$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d1*100 , y = two_ion_MIXDER_var$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d1*100,rev(d1*100)),c(two_ion_MIXDER_var$CI_lower * 100, rev(two_ion_MIXDER_var$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)

#MIXDER of 6-ion with 50 different dosage points and Monte Carlo Ribbon Plot using covariances
d3 = c(seq(0, 0.009, 0.001), seq(0.01, 0.4, by = 0.01))
out = MIXDER_function(r = rep(1/6, times = 6), L = c(75, 100, 125, 175, 195, 240), Z.b = c(595, 690, 770, 1075, 1245, 1585), d = c(seq(0, 0.009, 0.001), seq(0.01, 0.4, by = 0.01)))
six_ion_MIXDER = data.frame(d = out[, 1], CA = out[, 2] + 0.00071)
ninty_five_CI_lower = vector(length = length(d3))
ninty_five_CI_upper = vector(length = length(d3))

# for (i in 1:length(d3)){
#   info = vector(length = 0)
#   for (j in 1:MM){
#     info = c(info, MIXDER_function(r = rep(1/6, times = 6), c(0,  six_ion_MIXDER$d[i]),L = c(75, 100, 125, 175, 195, 240), Z.b = c(595, 690, 770, 1075, 1245, 1585), eta0 = eta0_MC[j], eta1 = eta1_MC[j], sig0 = sig0_MC[j], kap = kap_MC[j])[,2][2])
#   }
#   info = sort(info)
#   ninty_five_CI_lower[i] = quantile(info,(1-0.95)/2) + 0.00071
#   ninty_five_CI_upper[i] = quantile(info,1-(1-0.95)/2) + 0.00071
# }
# six_ion_MIXDER$CI_lower = ninty_five_CI_lower
# six_ion_MIXDER$CI_upper = ninty_five_CI_upper
# six_ion_MIXDER$simpleeffect = IDER(d = 1/6*d3, L = 125, Z.b = 770) + IDER(d = 1/6*d3, L = 175, Z.b = 1075) + IDER(d = 1/6*d3, L = 75, Z.b = 595) + IDER(d = 1/6*d3, L = 100, Z.b = 690) + IDER(d = 1/6*d3, L = 195, Z.b = 1245) + IDER(d = 1/6*d3, L = 240, Z.b = 1585) + 0.00071

######################## 08.20.18 EGH REFACTORING START
start_time <- Sys.time()
DERs <- matrix(nrow = MM, ncol = length(d3))
# for (i in 1:length(d3)){ # EGH 08.20.18 I strongly believe this is the area of the slowdown; we are calculating 29500 more curves then necessary.
# info = vector(length = 0)
for (j in 1:MM) {
  DERs[j, ] <- MIXDER_function(r = rep(1/6, times = 6), 
                               d = six_ion_MIXDER$d[1:length(d3)], 
                               L = c(75, 100, 125, 175, 195, 240), 
                               Z.b = c(595, 690, 770, 1075, 1245, 1585),
                               eta0 = eta0_MC[j], eta1 = eta1_MC[j], 
                               sig0 = sig0_MC[j], kap = kap_MC[j])[, 2]
  cat(paste("  Currently at Monte Carlo step:", toString(j), "of", 
            toString(MM)), sprintf('\r'))
}
for (i in 1:length(d3)) {
  sample_values <- sort(DERs[, i])
  # Returning resulting CI
  ninty_five_CI_lower[i] <- sample_values[ceiling((1 - 0.95) / 2 * MM)]
  ninty_five_CI_upper[i] <- sample_values[(0.95 + (1 - 0.95) / 2) * MM]
  # ninty_five_CI_lower[i] = quantile(info,(1 - 0.95) / 2) + 0.00071
  # ninty_five_CI_upper[i] = quantile(info,1 - (1 - 0.95) / 2) + 0.00071
}
start_time <- Sys.time()

################### REFACTORING END

six_ion_MIXDER$CI_lower = ninty_five_CI_lower
six_ion_MIXDER$CI_upper = ninty_five_CI_upper
six_ion_MIXDER$simpleeffect = IDER(d = 1/6*d3, L = 125, Z.b = 770) + IDER(d = 1/6*d3, L = 175, Z.b = 1075) + IDER(d = 1/6*d3, L = 75, Z.b = 595) + IDER(d = 1/6*d3, L = 100, Z.b = 690) + IDER(d = 1/6*d3, L = 195, Z.b = 1245) + IDER(d = 1/6*d3, L = 240, Z.b = 1585) + 0.00071


#Get the individual IDERS
six_ion_MIXDER$oxygen = IDER(d = d3, L = 75, Z.b = 595) + 0.00071
six_ion_MIXDER$silicon = IDER(d = d3, L = 100, Z.b = 690) + 0.00071
six_ion_MIXDER$titanium = IDER(d = d3, L = 125, Z.b = 770) + 0.00071
six_ion_MIXDER$ironsix = IDER(d = d3, L = 175, Z.b = 1075) + 0.00071
six_ion_MIXDER$ironfour = IDER(d = d3, L = 195, Z.b = 1245) + 0.00071
six_ion_MIXDER$ironthree = IDER(d = d3, L = 240, Z.b = 1585) + 0.00071

#The graphing part
d3 <- six_ion_MIXDER$d
CA3 <- six_ion_MIXDER$CA
simpleeffect3 <- six_ion_MIXDER$simpleeffect
silicon3 <- six_ion_MIXDER$silicon
titanium3 <- six_ion_MIXDER$titanium
ironthree3 <- six_ion_MIXDER$ironthree
ironfour3 <- six_ion_MIXDER$ironfour
ironsix3 <- six_ion_MIXDER$ironsix
oxygen3 <- six_ion_MIXDER$oxygen
plot(x = d3 * 100, y = CA3 * 100, type = "l", col = "red")
lines(x = d3 * 100, y = simpleeffect3 * 100, col = "black", lty = 2, lwd = 0.5)
lines(x = d3 * 100, y = silicon3* 100, col = "green")  
lines(x = d3 * 100, y = titanium3* 100, col = "green")
lines(x = d3 * 100, y = ironthree3* 100, col = "green")  
lines(x = d3 * 100, y = ironfour3* 100, col = "green")
lines(x = d3 * 100, y = ironsix3* 100, col = "green")  
lines(x = d3 * 100, y = oxygen3* 100, col = "green")
lines(x= d3*100 , y = six_ion_MIXDER$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d3*100 , y = six_ion_MIXDER$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d3*100,rev(d3*100)),c(six_ion_MIXDER$CI_lower * 100, rev(six_ion_MIXDER$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)

#6-ion 50 doses without using covariances
CI_function_MIXDER_var = function(d, d_interested, r, interval = 0.95, L, Z.beta) {
  MIXDER_curve = list(0)
  for (i in 1:MM) {
    MIXDER_curve[[i]] = MIXDER_function(r = r, d = d, L = L, Z.b = Z.b, eta0 = eta0_MC_var[i], eta1 = eta1_MC_var[i], sig0 = sig0_MC_var[i], kap = kap_MC_var[i])
  }
  info = vector(length = 0)
  for (i in 1:MM) {
    info = c(info, MIXDER_curve[[i]][, 2][2])
  }
  info = sort(info)
  lower_bound = info[(1-interval)/2*MM]
  upper_bound = info[(interval + (1-interval)/2)*MM]
  CI = c(lower_bound, upper_bound)
  return(CI)
}

six_ion_MIXDER_var = six_ion_MIXDER
ninty_five_CI_lower = vector(length = 0)
ninty_five_CI_upper = vector(length = 0)
a = vector(length = 0)

for (i in 2:length(d3)) {
  a = CI_function_MIXDER_var(d = c(0, six_ion_MIXDER_var$d[i]), r = rep(1/6, times = 6), L = c(75, 100, 125, 175, 195, 240), Z.b = c(595, 690, 770, 1075, 1245, 1585))
  ninty_five_CI_lower = c(ninty_five_CI_lower, a[1])
  ninty_five_CI_upper = c(ninty_five_CI_upper, a[2])
}
six_ion_MIXDER_var$CI_lower = c(0, ninty_five_CI_lower + 0.00071)
six_ion_MIXDER_var$CI_upper = c(0, ninty_five_CI_upper + 0.00071)
six_ion_MIXDER_var$simpleeffect = IDER(d = 1/6*d3, L = 125, Z.b = 770) + IDER(d = 1/6*d3, L = 175, Z.b = 1075) + IDER(d = 1/6*d3, L = 75, Z.b = 595) + IDER(d = 1/6*d3, L = 100, Z.b = 690) + IDER(d = 1/6*d3, L = 195, Z.b = 1245) + IDER(d = 1/6*d3, L = 240, Z.b = 1585) + 0.00071

#Get the individual IDERS
six_ion_MIXDER_var$oxygen = IDER(d = d3, L = 75, Z.b = 595) + 0.00071
six_ion_MIXDER_var$silicon = IDER(d = d3, L = 100, Z.b = 690) + 0.00071
six_ion_MIXDER_var$titanium = IDER(d = d3, L = 125, Z.b = 770) + 0.00071
six_ion_MIXDER_var$ironsix = IDER(d = d3, L = 175, Z.b = 1075) + 0.00071
six_ion_MIXDER_var$ironfour = IDER(d = d3, L = 195, Z.b = 1245) + 0.00071
six_ion_MIXDER_var$ironthree = IDER(d = d3, L = 240, Z.b = 1585) + 0.00071