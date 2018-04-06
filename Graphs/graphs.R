
library(deSolve) # package for solving differential equations
library(minpack.lm) # package for non-linear regression #rks to laz: I think we probably can just use nls() in stats, not nlsLM from linpack. Please check in R documentation if there is any functional difference at all
library(mvtnorm) # package for calculating confidence intervals by Monte Carlo simulation based on variance-covariance matrices #rks to laz: I added to comment.Please check that my addition is OK.
rm(list=ls())
set.seed(19970101)

####2-ion 60 dose
d = c(seq(0, 0.009, 0.001), seq(0.01, 0.5, by = 0.01))
#Assuming the working directory is the location of this R file (On Github the csv file is in the same location as this R file)
two_ion_MIXDER = read.csv("2ionMonteCarlo.csv") #This is the Monte Carlo results with seed 19970101

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





####2-ion 12 dose
d_small = d[seq(1,60, by = 5)]
two_ion_MIXDER_small <- two_ion_MIXDER[seq(1,60, by = 5), ]

#The graphing part
d2 <- two_ion_MIXDER_small$d
CA2 <- two_ion_MIXDER_small$CA
simpleeffect2 <- two_ion_MIXDER_small$simpleeffect
silicon2 <- two_ion_MIXDER_small$silicon
ironsix2 <- two_ion_MIXDER_small$ironsix
plot(x = d2 * 100, y = CA2 * 100, type = "l", col = "red")
lines(x = d2 * 100, y = simpleeffect2 * 100, col = "black", lty = 2, lwd = 0.5)
lines(x = d2 * 100, y = silicon2* 100, col = "green")  
lines(x = d2 * 100, y = ironsix2* 100, col = "green")
lines(x= d2*100 , y = two_ion_MIXDER_small$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d2*100 , y = two_ion_MIXDER_small$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d2*100,rev(d2*100)),c(two_ion_MIXDER_small$CI_lower * 100, rev(two_ion_MIXDER_small$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)



####2-ion 12 dose without covariances (larger CI's)
#Assuming the working directory is the location of this R file (On Github the csv file is in the same location as this R file)
two_ion_MIXDER_small_var <- read.csv("2ionMonteCarlovar.csv")

#The graphing part
plot(x = d2 * 100, y = CA2 * 100, type = "l", col = "red")
lines(x = d2 * 100, y = simpleeffect2 * 100, col = "black", lty = 2, lwd = 0.5)
lines(x = d2 * 100, y = silicon2* 100, col = "green")  
lines(x = d2 * 100, y = ironsix2* 100, col = "green")
lines(x= d2*100 , y = two_ion_MIXDER_small_var$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d2*100 , y = two_ion_MIXDER_small_var$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d2*100,rev(d2*100)),c(two_ion_MIXDER_small_var$CI_lower * 100, rev(two_ion_MIXDER_small_var$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)


####6-ion 60 dose 
six_ion_MIXDER <- read.csv("6ionMonteCarlo.csv")

#the graphing part
d11 <- six_ion_MIXDER$d
CA11 <- six_ion_MIXDER$CA
simpleeffect11 <- six_ion_MIXDER$simpleeffect
silicon11 <- six_ion_MIXDER$silicon
titanium11 <- six_ion_MIXDER$titanium
ironthree11 <- six_ion_MIXDER$ironthree
ironfour11 <- six_ion_MIXDER$ironfour
ironsix11 <- six_ion_MIXDER$ironsix
oxygen11 <- six_ion_MIXDER$oxygen
plot(x = d11 * 100, y = CA11 * 100, type = "l", col = "red")
lines(x = d11 * 100, y = simpleeffect11 * 100, col = "black", lty = 2, lwd = 0.5)
lines(x = d11 * 100, y = silicon11* 100, col = "green")  
lines(x = d11 * 100, y = titanium11* 100, col = "green")
lines(x = d11 * 100, y = ironthree11* 100, col = "green")  
lines(x = d11 * 100, y = ironfour11* 100, col = "green")
lines(x = d11 * 100, y = ironsix11* 100, col = "green")  
lines(x = d11 * 100, y = oxygen11* 100, col = "green")
lines(x= d11*100 , y = six_ion_MIXDER$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d11*100 , y = six_ion_MIXDER$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d11*100,rev(d11*100)),c(six_ion_MIXDER$CI_lower * 100, rev(six_ion_MIXDER$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)