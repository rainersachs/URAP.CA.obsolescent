source('/Users/zhaoliyang/Desktop/URAP\ FALL2018/MonteCarlo.R')
#====== an example of a plot: 6 HZE, each 1/6 of total dose ======#
#Simple Effect Additivity (black line)
SEA=function(d){
  IDER(d/6,L[1],Z.b[1])+IDER(d/6,L[2],Z.b[2])+IDER(d/6,L[3],Z.b[3])+IDER(d/6,L[4],Z.b[4])+IDER(d/6,L[5],Z.b[5])+IDER(d/6,L[6],Z.b[6])
}

r=rep(1/6,6);L = c(75, 100, 125, 175, 195, 240); Z.b = c(595, 690, 770, 1075, 1245, 1585); dose=seq(0, 0.4, by=0.001)

MX=MIXDER_function(r, L, Z.b, d=dose) #for this mixture can't go much above 0.6#

IDER = function(d, L = NULL, Z.b = NULL, ions = NULL, eta0 = eta0_hat, eta1 = eta1_hat, sig0 = sig0_hat, kap = kap_hat) {
  #SCW: Added an option for inputting ion names instead of L, Z.b to decrease chance of transcription error and keep track of which ions are being mixed. (Old code should still work for now, but ideally we will change everything to this format in the future)
  #ions is a vector of strings of the names of the ions
  if (is.null(ions)){
    P = (1-exp(-Z.b/kap_hat))^2
    sig = sig0*P + 0.041/6.24*L*(1-P)
    eta = eta0*L*exp(-eta1*L)
    return(sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d)))
  }
  else{#After testing, this results in the same output as the old one besides rounding error that cannot be seen in the available digits
    r = 1/length(ions)
    info_table = modified_df %>% group_by(ion, L, Z.b) %>% summarise()
    info_table = suppressWarnings(left_join(data.frame(ions), info_table, by = c("ions" = "ion")))
    output = 0
    for (i in 1: nrow(info_table)){
      output = output + IDER(d*r, L = info_table$L[i], Z.b = info_table$Z.b[i])
    }
    return(output)
  }
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

#========== Figure for MonteCarlo=============#
monte_carlo(c("Si", "Fe600"))
monte_carlo(c("Si", "Fe600"), cov = F)
monte_carlo(c("O", "Si", "Ti", " Fe600", "Fe450", "Fe300"))

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
two_ion_MIXDER$simpleeffect = IDER(d = 0.5*d1, L = 100, Z.b = 690) + IDER(d = 0.5*d1, L = 175, Z.b = 1075) + 0.00071 #SCW: I think here we should add background prevalence to every IDER instead (although this won't change the result much)
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


