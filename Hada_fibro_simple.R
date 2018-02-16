#This is free, open-source software under the (lenient) GNU GPLv3. It comes with no warranty. 
#.Rmd version was written by Dae Woong Ham between Sept. 2016 and May 2017. Some quality control by Rainer (Ray) K. Sachs (rks) May-August 2017.
# This R version was written by Liyang (Andy) Zhao (laz) UCB semester fall 2017. Some quality control by rks and Julien Yu.
#Script concerns synergy analysis of WGE simple chromosome aberrations (CA) induced in 82-6 fibroblast cells by simulated GCR (Galactic Cosmic Radiation) mixed fields. It is an R version of parts of GCRfibroCA2GH.Rmd
#The script uses various mixture component IDERs (Individual Dose Effect Relations), summarized in "16Cacao" = 
# = Cacao, Hada, Saganti, George and Cucinotta. "Relative Biological Effectiveness of HZE Particles for Chromosomal Exchanges and Other Surrogate Cancer Risk Endpoints." PLoS One 11(4): e0153998. (2016)].
#The libraries needed for this script
library(deSolve) # package for solving differential equations
library(minpack.lm) # package for non-linear regression #rks to laz: I think we probably can just use nls() in stats, not nlsLM from linpack. Please check in R documentation if there is any functional difference at all
library(mvtnorm) # package for calculating confidence intervals by Monte Carlo simulation based on variance-covariance matrices #rks to laz: I added to comment.Please check that my addition is OK.
rm(list=ls())

#Create dataframes that store the fibroblast WGE simple CA data used in 16Cacao 

Oxygen = data.frame(d = c(0, .0125, .02, .025, .05, .075, .1, .2, .4),CA = c(.24, 1.66, 2.43, 2.37, 1.16, 2.85, 2.58, 6.94, 6.91))
# Some GCR components are high-speed Oxygen nuclei that are almost fully ionized. d=dose; CA are per hundred cells.

Si = data.frame(d = c(0, .02, .04, .06, .08, .1, .12, .2, .4, .8, 1.2), CA = c(.11, 1.26, 1.14, 1.58, 1.22, 1.89, 3.47, 4.6, 9.79, 27.01, 38.84))

Fe600 = data.frame(d = c(0, .01, .02, .04, .06, .08, .1, .12, .2, .4, .8), CA = c(.13, .76, .99, 1.2, 1.74, 1.28, 1.2, 1.7, 3.02, 5.52, 12.42)) 
#600 refers to the energy in MeV per atomic mass unit in this Iron beam

Fe450 = data.frame(d = c(0, .02, .04, .06, .08, .1, .2, .4), CA = c(0, .86, .6, .8, 1.22, 2.02, 2.3, 4.77))

Fe300 = data.frame(d = c(0, .005, .01,  0.02, .04, .07, .1, .2, .4, .8), CA = c(0.41, 1.23, 1.47, 1.22, .97, 1.46, 1.21, 4.38, 6.22, 13.6))

Ti = data.frame(d = c(0,  0.02, .04, .06, .08, .1, .15, .3, .6), CA = c(0, 1.99, 1.88, 1.44, 2.67, 2.57, 2.50, 5.64, 11.19))

param = data.frame(ion = c("O", "Si", "Ti", "Fe600", "Fe450", "Fe300"),
                   Z = c(8, 14, 22, 26, 26, 26), L = c(75, 100, 125, 175, 195, 240), #Z=atomic charge; L=LET
                   Z.b = c(595, 690, 770, 1075, 1245, 1585))#Z^2/beta*^2

#putting it in one big data frame. # the data frame incorporates a correction to Fe600 at dose 0.06
big_df = rbind(Oxygen, Si, Ti, Fe600, Fe450, Fe300)
big_df$Z = rep(param$Z, times = c(9, 11, 9, 11, 8, 10))
big_df$Z.b = rep(param$Z.b, times = c(9, 11, 9, 11, 8, 10))
big_df$L = rep(param$L, times = c(9, 11, 9, 11, 8, 10))
big_df$error = c(0.24, 0.63, 0.77, 0.75, 0.52, 0.82, 0.78, 1.31, 1.59, 0.12, 0.05, 0.07, 0.56, 0.18, 0.60, 1.23, 1.60, 1.55, 4.27, 7.21, 0, 0.70, 0.66, 0.59, 0.80, 0.78, 0.48, 1.15, 2.39, 0.16, 0.38, 0.24, 0.21, 0.43, 0.37, 0.54, 0.17, 0.55, 1.75, 2.59, 0, 0.43, 0.34, 0.40, 0.50, 0.64, 0.73, 1.09, 0.29, 0.55, 0.60, 0.55, 0.49, 0.60, 0.54, 1.03, 1.22, 3.62) #error bars for the CA measurements
big_df$ion = rep(param$ion, times = c(9, 11, 9, 11, 8, 10))

#Next modify the data frame to get rid of the zero dose points. Background CA frequency was determined seperately.
modified_df = big_df[big_df$d != 0, ]
modified_df$CA = modified_df$CA*0.01 
modified_df$error = modified_df$error*0.01
big_df$CA = big_df$CA * 0.01
big_df$error = big_df$error * 0.01
big_df$errorbar_lower = big_df$CA - big_df$error
big_df$errorbar_upper = big_df$CA + big_df$error

#NTE1 and NTE2 models in 16Cacao using their parameters. NTE is used in 16Cacao to signify that non-targeted effects are included in the model; Conventional targeted effects (TE) are included in all models.
#NTE1 function
NTE1_function = function(d, L, Z.b, eta0 = 0.00011, eta1 = 0.007, sig0 = 6.12, kap = 796) {
  0.0017 + eta0*L*exp(-eta1*L)*(d != 0) + 
    (6.242*(d/L))*(sig0*(1-exp(-Z.b/kap))^2 + 0.041/6.24*L*(1 - (1-exp(-Z.b/kap))^2))
} 

#NTE2 function
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

#nls (non-linear least square) method to get the parameters needed (4 parameter estimation) #rks to laz: see note under issues
IDER_model = nlsLM(CA ~ IDER(d, L, Z.b, eta00, eta10, sig00, kap0), data = modified_df, start = list(eta00 = 0.001, eta10 = 0.01, sig00 = 5, kap0 = 500), 
                   weights = (1/(modified_df$error)^2))
coef(IDER_model)
#e.g. eta0 = 1.484687e-04, eta1 = 3.513927e-03, sig0 = 4.149660e+00, kap = 4.688287e+02
vcov(IDER_model)# variance-covariance matrix, needed later for analyzing 95% confidence limits in baseline no-synergy/no-antagonism MIXDER (Mixture dose effect relation)
summary(IDER_model, cor = TRUE) #model parameters and their correlation matrix

XX <- coef(IDER_model)
eta00 <- as.numeric(XX[1])
eta10 <- as.numeric(XX[2])
sig00 <- as.numeric(XX[3])
kap0 <- as.numeric(XX[4])


#R function to give Information criteria (AIC and BIC) #rks to laz: R has AIC functions. Please check that they give the same ordering (not necessarily the same values) as the ones constructed below by Dae or you.
#L_function gives the residuals squared
L_function = function(func, eta0, eta1, sig0, kap) {
  a = vector(length = 0)
  for (i in 1:length(modified_df[, 1])) {
    a = c(a, modified_df$CA[i] - func(d = modified_df$d[i], L = modified_df$L[i], Z.b = modified_df$Z.b[i], eta0 = eta0, eta1 = eta1, sig0 = sig0, kap = kap))
  }
  return(a^2)
}
L_NTE1 = L_function(NTE1_function, eta0 = 0.00011, eta1 = 0.007, sig0 = 6.12, kap = 796)
L_NTE2 = L_function(NTE2_function, eta0 = 0.00047, eta1 = 0.011, sig0 = 6.75, kap = 590)
L_IDER = L_function(IDER, eta0 = eta00, eta1 = eta10, sig0 = sig00, kap = kap0)

#Since all models are weighted least square regression, will weight our model with our weights to get the WRSS (weighted residual squared sum)
WRSS_NTE1 = sum((1/modified_df$error^2)*L_NTE1)
WRSS_NTE2 = sum((1/modified_df$error^2)*L_NTE2)
WRSS_IDER = sum((1/modified_df$error^2)*L_IDER)

#functions for AIC and BIC calculation for Weighted Least Square regression
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
information_critera_df ## RKS runs well up to here. Lots of checks made 8/12/2017

# #Baseline no-synergy/no-antagonism MIXDER  based on any input
MIXDER_function = function(r, L, Z.b, d = seq(0, 0.2, by = 0.001), eta0 = eta00, eta1 = eta10, sig0 = sig00 ,kap = kap0) {
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


#Graphs for NTE1 and our IDER models: 12_ion figure
#When creating any sort of MIXDER figures we use this IDER now without the background effect then add the background effect at the end.
IDER = function(d, L, Z.b, eta0 = eta00, eta1 = eta10, sig0 = sig00, kap = kap0) {
  P = (1-exp(-Z.b/kap0))^2
  sig = sig0*P + 0.041/6.24*L*(1-P)
  eta = eta0*L*exp(-eta1*L)
  sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d))#sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-.5*10^3*d))#don't use
} 
r=rep(1/6,6);L = c(75, 100, 125, 175, 195, 240); Z.b = c(595, 690, 770, 1075, 1245, 1585);dose=seq(0,.4,by=0.001)
MX=MIXDER_function(r,L,Z.b,d=dose)#for this mixture can't go much above 0.6#
plot(MX[,1],MX[,2],type='l',bty='l',col='red',ann='F',ylim=c(0,.10))
for (ii in 1:length(L)){lines(dose,IDER (dose,L[ii],Z.b[ii]),col='green')}
#SEA
SEA=function(d){
  IDER(d/6,L[1],Z.b[1])+IDER(d/6,L[2],Z.b[2])+IDER(d/6,L[3],Z.b[3])+IDER(d/6,L[4],Z.b[4])+IDER(d/6,L[5],Z.b[5])+IDER(d/6,L[6],Z.b[6])
}
lines(dose,SEA(dose),lty=2)

#lines(dose,sum(SEA))
#Comparing the NTE1 model in 16Cacao with our NASA SNTE model

#Use these dose points to generate the curves necessary. Had to include a lot more dose points at lower doses where it is very sensitive to change
# d_O = c(0, 1e-6*(1:9), 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:9), 0.01*(10:40))
# d_1 = c(0, 1e-6*(1:9), 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:9), 1e-2*(1:9), 0.01*(10:15))
# d_2 = c(0, 1e-6*(1:9), 1e-5*(1:9), 1e-4*(1:9), 1e-3*(1:300))
# plot(c(0,.3),c(0,0.04), bty='l', ann='F')
# lines(d_2,IDER(d_2,240,1585), bty='l', ann='F')
# postscript("NTE1_model test.eps", width = 5, height = 7)
# par(mfrow=c(2,1))
# 
# plot(x = big_df$d[1:9]*100, y = big_df$CA[1:9]*100, main = "Oxygen", ylim = c(0, 10),pch=19, ylab = "CA", xlab = "Dose (cGy)")
# arrows(big_df$d[1:9]*100, big_df$errorbar_lower[1:9]*100, big_df$d[1:9]*100, big_df$errorbar_upper[1:9]*100, length=0.02, angle=90, code=3)
# lines(x = d_O*100, y = 100*(0.00071 + IDER(d = d_O, L = 75, Z.b = 595)), col = "red")
# lines(x = d_O*100, y = 100 * NTE1_function(d = d_O, L= 75, Z.b = 595), col = "blue", lty = 2)
# 
# plot(x = big_df$d[1:7]*100, y = big_df$CA[1:7]*100, main = "Oxygen",pch=19, ylab = "CA", xlab = "Dose (cGy)", ylim = c(0, 4), xlim = c(0, 15))
# arrows(big_df$d[1:7]*100, big_df$errorbar_lower[1:7]*100, big_df$d[1:7]*100, big_df$errorbar_upper[1:7]*100, length=0.02, angle=90, code=3)
# lines(x = d_1*100, y = 100*(0.00071 + IDER(d = d_1, L = 75, Z.b = 595)), col = "red")
# lines(x = d_1*100, y = 100 * NTE1_function(d = d_1, L= 75, Z.b = 595), col = "blue", lty = 2)




#Monte Carlo
sig = vcov(IDER_model)
set.seed(10)
monte_carlo_parameters = rmvnorm(n = 500, mean = c(eta0 = eta00, eta1 = eta10, sig0 = sig00, kap = kap0), sigma = sig)
eta0_MC = monte_carlo_parameters[, 1]
eta1_MC = monte_carlo_parameters[, 2]
sig0_MC = monte_carlo_parameters[, 3]
kap_MC = monte_carlo_parameters[, 4]
kap_MC[kap_MC <= 1e-6] = 1e-5 

#note because our last parameter were less significant when inputting our variance covariance matrix the parameters had to actually be readjusted to fix for negative value. This is important because if kap parameter went negative our model is nonsensical. (Luckily the fix was only for exactly 15 out of 500 MC samples or only roughly 3%)

#This is a general CI_function where you just get a CI for a specific dose point for any synergy analysis using the monte carlo simulations using vcov().
CI_function_MIXDER = function(d, d_interested, r, interval = 0.95, L, Z.beta) {
  MIXDER_curve = list(0)
  for (i in 1:500) {
    MIXDER_curve[[i]] = MIXDER_function(r = r, d = d, L = L, Z.b = Z.b, eta0 = eta0_MC[i], eta1 = eta1_MC[i], sig0 = sig0_MC[i], kap = kap_MC[i])
  }
  info = vector(length = 0)
  for (i in 1:500) {
    info = c(info, MIXDER_curve[[i]][, 2][2])
  }
  info = sort(info)
  lower_bound = info[(1-interval)/2*500]
  upper_bound = info[(interval + (1-interval)/2)*500]
  CI = c(lower_bound, upper_bound)
  return(CI)
}

#Monte Carlo using the error bars from summary()
set.seed(11)
eta0_MC_2 = rnorm(1000, mean = eta00, sd = 4.631e-05)
eta1_MC_2 = rnorm(1000, mean = eta10, sd = 1.617e-04)
sig0_MC_2 = rnorm(1000, mean = sig00, sd = 1.966e+00)
kap_MC_2 = rnorm(1000, mean = kap0, sd = 2.628e+02)
kap_MC_2[kap_MC_2 <= 1e-6] = 1e-5 

CI_function_MIXDER_MC_2 = function(d, d_interested, r, interval = 0.95, L, Z.beta) {
  MIXDER_curve = list(0)
  for (i in 1:500) {
    MIXDER_curve[[i]] = MIXDER_function(r = r, d = d, L = L, Z.b = Z.b, eta0 = eta0_MC_2[i], eta1 = eta1_MC_2[i], sig0 = sig0_MC_2[i], kap = kap_MC_2[i])
  }
  info = vector(length = 0)
  for (i in 1:500) {
    info = c(info, MIXDER_curve[[i]][, 2][2])
  }
  info = sort(info)
  lower_bound = info[(1-interval)/2*500]
  upper_bound = info[(interval + (1-interval)/2)*500]
  CI = c(lower_bound, upper_bound)
  return(CI)
}
out = MIXDER_function(r = rep(1/2, times = 2), L = c(100, 175), Z.b = c(690, 1075), d = c(seq(0, 0.01, 0.001), seq(0.01, 0.5, by = 0.01)))
two_ion_MIXDER = data.frame(d = out[, 1], CA = out[, 2] + 0.00071)
d = two_ion_MIXDER$d


#Then the 95% CI multivariate MC simulated ribbon
ninty_five_CI_lower = vector(length = 0)
ninty_five_CI_upper = vector(length = 0)
a = vector(length = 0)
for (i in 2:length(d)) {
  a = CI_function_MIXDER(d = c(0, two_ion_MIXDER$d[i]), d_interested = two_ion_MIXDER$d[i], r = rep(1/2, times = 2),  L = c(100, 175), Z.b = c(690, 1075))
  ninty_five_CI_lower = c(ninty_five_CI_lower, a[1])
  ninty_five_CI_upper = c(ninty_five_CI_upper, a[2])
}
two_ion_MIXDER$CI_lower = c(0, ninty_five_CI_lower + 0.00071)
two_ion_MIXDER$CI_upper = c(0, ninty_five_CI_upper + 0.00071)

#Get the simple effect additivity MIXDER
two_ion_MIXDER$simpleeffect = IDER(d = 0.5*d, L = 100, Z.b = 690) + IDER(d = 0.5*d, L = 175, Z.b = 1075) + 0.00071
#Get the individual IDERS
two_ion_MIXDER$silicon = IDER(d = d, L = 100, Z.b = 690) + 0.00071
two_ion_MIXDER$ironsix = IDER(d = d, L = 175, Z.b = 1075) + 0.00071
save(two_ion_MIXDER, file = "two_ion_95%.Rda")

d1 <- two_ion_MIXDER$d
CA1 <- two_ion_MIXDER$CA
simpleeffect1 <- two_ion_MIXDER$simpleeffect
silicon1 <- two_ion_MIXDER$silicon
ironsix1 <- two_ion_MIXDER$ironsix
plot(x = d1 * 100, y = CA1 * 100, type = "l", col = "red")
lines(x = d1 * 100, y = simpleeffect1 * 100, col = "black", lty = 2, lwd = 0.5) + 
  lines(x = d1 * 100, y = silicon1* 100, col = "green")  
lines(x = d1 * 100, y = ironsix1* 100, col = "green")
lines(x= d1*100 , y = two_ion_MIXDER$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d1*100 , y = two_ion_MIXDER$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d1*100,rev(d1*100)),c(two_ion_MIXDER$CI_lower * 100, rev(two_ion_MIXDER$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)

#This is for the six_ion synergy analysis incorporating our monte carlo simulations with the yellow ribbon as the 95% confidence interval using our monte carlo simulated parameters with two panels for two methods of the MC simulation. Similarly, blues curves are individual IDERs and red curve the mixed IDER and the black the simple effect addtivity curve.

#First the MIXDER of the true curve 
out = MIXDER_function(r = rep(1/6, times = 6), L = c(75, 100, 125, 175, 195, 240), Z.b = c(595, 690, 770, 1075, 1245, 1585), d = c(seq(0, 0.01, 0.001), seq(0.01, 0.4, by = 0.01)))
six_ion_MIXDER = data.frame(d = out[, 1], CA = out[, 2] + 0.00071)
d = six_ion_MIXDER$d
#Then the 95% CI multivariate MC simulated ribbon
ninty_five_CI_lower = vector(length = 0)
ninty_five_CI_upper = vector(length = 0)
a = vector(length = 0)
for (i in 2:length(d)) {
  a = CI_function_MIXDER(d = c(0, six_ion_MIXDER$d[i]), d_interested = six_ion_MIXDER$d[i], r = rep(1/6, times = 6), L = c(75, 100, 125, 175, 195, 240), Z.b = c(595, 690, 770, 1075, 1245, 1585))
  ninty_five_CI_lower = c(ninty_five_CI_lower, a[1])
  ninty_five_CI_upper = c(ninty_five_CI_upper, a[2])
}
six_ion_MIXDER$CI_lower = c(0, ninty_five_CI_lower + 0.00071)
six_ion_MIXDER$CI_upper = c(0, ninty_five_CI_upper + 0.00071)

#Then the 95% CI for the MC simulated independent gaussians
ninty_five_CI_lower_two = vector(length = 0)
ninty_five_CI_upper_two = vector(length = 0)
a = vector(length = 0)
for (i in 2:length(d)) {
  a = CI_function_MIXDER_MC_2(d = c(0, six_ion_MIXDER$d[i]), d_interested = six_ion_MIXDER$d[i], r = rep(1/6, times = 6), L = c(75, 100, 125, 175, 195, 240), Z.b = c(595, 690, 770, 1075, 1245, 1585))
  ninty_five_CI_lower_two = c(ninty_five_CI_lower_two, a[1])
  ninty_five_CI_upper_two = c(ninty_five_CI_upper_two, a[2])
}
six_ion_MIXDER$CI_lower_two = c(0, ninty_five_CI_lower_two + 0.00071)
six_ion_MIXDER$CI_upper_two = c(0, ninty_five_CI_upper_two + 0.00071)

#Get the simple effect additivity MIXDER 
six_ion_MIXDER$simpleeffect = IDER(d = 1/6*d, L = 125, Z.b = 770) + IDER(d = 1/6*d, L = 175, Z.b = 1075) + IDER(d = 1/6*d, L = 75, Z.b = 595) + IDER(d = 1/6*d, L = 100, Z.b = 690) + IDER(d = 1/6*d, L = 195, Z.b = 1245) + IDER(d = 1/6*d, L = 240, Z.b = 1585) + 0.00071

#Get the individual IDERS
six_ion_MIXDER$oxygen = IDER(d = d, L = 75, Z.b = 595) + 0.00071
six_ion_MIXDER$silicon = IDER(d = d, L = 100, Z.b = 690) + 0.00071
six_ion_MIXDER$titanium = IDER(d = d, L = 125, Z.b = 770) + 0.00071
six_ion_MIXDER$ironsix = IDER(d = d, L = 175, Z.b = 1075) + 0.00071
six_ion_MIXDER$ironfour = IDER(d = d, L = 195, Z.b = 1245) + 0.00071
six_ion_MIXDER$ironthree = IDER(d = d, L = 240, Z.b = 1585) + 0.00071
save(six_ion_MIXDER, file = "six_ion_95%.Rda")

d2 <- six_ion_MIXDER$d
CA2 <- six_ion_MIXDER$CA
simpleeffect2 <- six_ion_MIXDER$simpleeffect
silicon2 <- six_ion_MIXDER$silicon
titanium2 <- six_ion_MIXDER$titanium
ironthree2 <- six_ion_MIXDER$ironthree
ironfour2 <- six_ion_MIXDER$ironfour
ironsix2 <- six_ion_MIXDER$ironsix
oxygen2 <- six_ion_MIXDER$oxygen
plot(x = d2 * 100, y = CA2 * 100, type = "l", col = "red")
lines(x = d2 * 100, y = simpleeffect2 * 100, col = "black", lty = 2, lwd = 0.5) + 
lines(x = d2 * 100, y = silicon2* 100, col = "green")  
lines(x = d2 * 100, y = titanium2* 100, col = "green")
lines(x = d2 * 100, y = ironthree2* 100, col = "green")  
lines(x = d2 * 100, y = ironfour2* 100, col = "green")
lines(x = d2 * 100, y = ironsix2* 100, col = "green")  
lines(x = d2 * 100, y = oxygen2* 100, col = "green")
lines(x= d2*100 , y = six_ion_MIXDER$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d2*100 , y = six_ion_MIXDER$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d2*100,rev(d2*100)),c(six_ion_MIXDER$CI_lower * 100, rev(six_ion_MIXDER$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)


plot(x = d2 * 100, y = CA2 * 100, type = "l", col = "red")
lines(x = d2 * 100, y = simpleeffect2 * 100, col = "black", lty = 2, lwd = 0.5) + 
lines(x = d2 * 100, y = silicon2* 100, col = "green")  
lines(x = d2 * 100, y = titanium2* 100, col = "green")
lines(x = d2 * 100, y = ironthree2* 100, col = "green")  
lines(x = d2 * 100, y = ironfour2* 100, col = "green")
lines(x = d2 * 100, y = ironsix2* 100, col = "green")  
lines(x = d2 * 100, y = oxygen2* 100, col = "green")
lines(x= d2*100 , y = six_ion_MIXDER$CI_upper_two * 100, lty = 'dashed', col = 'red')
lines(x= d2*100 , y = six_ion_MIXDER$CI_lower_two * 100, lty = 'dashed', col = 'red')
polygon(c(d2*100,rev(d2*100)),c(six_ion_MIXDER$CI_lower_two * 100, rev(six_ion_MIXDER$CI_upper_two * 100)),col = rgb(1, 0, 0,0.5), border = NA)





