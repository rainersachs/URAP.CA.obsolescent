
library(deSolve) # package for solving differential equations
library(minpack.lm) # package for non-linear regression #rks to laz: I think we probably can just use nls() in stats, not nlsLM from linpack. Please check in R documentation if there is any functional difference at all
library(mvtnorm) # package for calculating confidence intervals by Monte Carlo simulation based on variance-covariance matrices #rks to laz: I added to comment.Please check that my addition is OK.
#rm(list=ls())
set.seed(19970101)

####2-ion 60 dose
d_2ion = c(seq(0, 0.009, 0.001), seq(0.01, 0.5, by = 0.01))
#Assuming the working directory is the location of this R file (On Github the csv file is in the same location as this R file)
two_ion_MIXDER = read.csv("2ionMonteCarlo.csv") #This is the Monte Carlo results with seed 19970101

#The graphing part
d1 <- two_ion_MIXDER$d
CA1 <- two_ion_MIXDER$CA
simpleeffect1 <- two_ion_MIXDER$simpleeffect
silicon1 <- two_ion_MIXDER$silicon
ironsix1 <- two_ion_MIXDER$ironsix
plot(x = d1 * 100, y = CA1 * 100, type = "l", col = "red", main="2-ion Model", sub="Full 60 Doses", 
     xlab="Dose * 100", ylab="CA * 100")
lines(x = d1 * 100, y = simpleeffect1 * 100, col = "black", lty = 2, lwd = 0.5)
lines(x = d1 * 100, y = silicon1* 100, col = "green")  
lines(x = d1 * 100, y = ironsix1* 100, col = "green")
lines(x= d1*100 , y = two_ion_MIXDER$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d1*100 , y = two_ion_MIXDER$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d1*100,rev(d1*100)),c(two_ion_MIXDER$CI_lower * 100, rev(two_ion_MIXDER$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)

###2-ion 60 dose without covariances (larger CI's)
#Assuming the working directory is the location of this R file (On Github the csv file is in the same location as this R file)
two_ion_MIXDER_var <- read.csv("2ionMonteCarlovar.csv")

#The graphing part
plot(x = d1 * 100, y = CA1 * 100, type = "l", col = "red", main="2-ion Model", sub="60 Doses without Covariances", 
     xlab="Dose * 100", ylab="CA * 100")
lines(x = d1 * 100, y = simpleeffect1 * 100, col = "black", lty = 2, lwd = 0.5)
lines(x = d1 * 100, y = silicon1* 100, col = "green")  
lines(x = d1 * 100, y = ironsix1* 100, col = "green")
lines(x= d1*100 , y = two_ion_MIXDER_var$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d1*100 , y = two_ion_MIXDER_var$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d1*100,rev(d1*100)),c(two_ion_MIXDER_var$CI_lower * 100, rev(two_ion_MIXDER_var$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)


####6-ion 50 dose
d_6ion <- c(seq(0, 0.009, 0.001), seq(0.01, 0.4, by = 0.01))
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
plot(x = d11 * 100, y = CA11 * 100, type = "l", col = "red", main="6-ion Model", sub="Full 50 Doses", 
     xlab="Dose * 100", ylab="CA * 100")
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

####6-ion 10 dose
d_6ion_small <- d_6ion[seq(1, 50, by = 5)] 
six_ion_MIXDER_small <- six_ion_MIXDER[seq(1, 50, by = 5),]

#the graphing part
d22 <- six_ion_MIXDER_small$d
CA22 <- six_ion_MIXDER_small$CA
simpleeffect22 <- six_ion_MIXDER_small$simpleeffect
silicon22 <- six_ion_MIXDER_small$silicon
titanium22 <- six_ion_MIXDER_small$titanium
ironthree22 <- six_ion_MIXDER_small$ironthree
ironfour22 <- six_ion_MIXDER_small$ironfour
ironsix22 <- six_ion_MIXDER_small$ironsix
oxygen22 <- six_ion_MIXDER_small$oxygen
plot(x = d22 * 100, y = CA22 * 100, type = "l", col = "red", main="6-ion Model", sub="10 Doses", 
     xlab="Dose * 100", ylab="CA * 100")
lines(x = d22 * 100, y = simpleeffect22 * 100, col = "black", lty = 2, lwd = 0.5)
lines(x = d22 * 100, y = silicon22* 100, col = "green")  
lines(x = d22 * 100, y = titanium22* 100, col = "green")
lines(x = d22 * 100, y = ironthree22* 100, col = "green")  
lines(x = d22 * 100, y = ironfour22* 100, col = "green")
lines(x = d22 * 100, y = ironsix22* 100, col = "green")  
lines(x = d22 * 100, y = oxygen22* 100, col = "green")
lines(x= d22*100 , y = six_ion_MIXDER_small$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d22*100 , y = six_ion_MIXDER_small$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d22*100,rev(d22*100)),c(six_ion_MIXDER_small$CI_lower * 100, rev(six_ion_MIXDER_small$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)

####6-ion 50 dose without covariances (larger Ci's)
six_ion_MIXDER_var <- read.csv("6ionMonteCarlovar.csv")

#the graphing part
d22 <- six_ion_MIXDER_var$d
CA22 <- six_ion_MIXDER_var$CA
simpleeffect22 <- six_ion_MIXDER_var$simpleeffect
silicon22 <- six_ion_MIXDER_var$silicon
titanium22 <- six_ion_MIXDER_var$titanium
ironthree22 <- six_ion_MIXDER_var$ironthree
ironfour22 <- six_ion_MIXDER_var$ironfour
ironsix22 <- six_ion_MIXDER_var$ironsix
oxygen22 <- six_ion_MIXDER_var$oxygen
plot(x = d22 * 100, y = CA22 * 100, type = "l", col = "red", main="6-ion Model", sub="50 Doses without Covariances", 
     xlab="Dose * 100", ylab="CA * 100")
lines(x = d22 * 100, y = simpleeffect22 * 100, col = "black", lty = 2, lwd = 0.5)
lines(x = d22 * 100, y = silicon22* 100, col = "green")  
lines(x = d22 * 100, y = titanium22* 100, col = "green")
lines(x = d22 * 100, y = ironthree22* 100, col = "green")  
lines(x = d22 * 100, y = ironfour22* 100, col = "green")
lines(x = d22 * 100, y = ironsix22* 100, col = "green")  
lines(x = d22 * 100, y = oxygen22* 100, col = "green")
lines(x= d22*100 , y = six_ion_MIXDER_var$CI_upper * 100, lty = 'dashed', col = 'red')
lines(x= d22*100 , y = six_ion_MIXDER_var$CI_lower * 100, lty = 'dashed', col = 'red')
polygon(c(d22*100,rev(d22*100)),c(six_ion_MIXDER_var$CI_lower * 100, rev(six_ion_MIXDER_var$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)
