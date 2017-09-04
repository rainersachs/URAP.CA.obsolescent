##This is free, open-source software under the (lenient) GNU GPLv3. It comes with no warranty. 
#It concerns synergy analysis of mixtures whose component IDER are LQ with generalized Lea-Catcheside dose protraction functional.  
#Written by Jian Gao; tested and modified by Ray Sachs and Julien Yu.
########################
library(deSolve)
rm(list=ls())
R=0.5#total mixture dose-rate here taken as 50 cGy per week
lambda=0.05#repair rate is assumed determined by target biology, and is the same for both ions
#Next line sets LQ-Lea parameters for both mixture components
alpha1=0;beta1=8.4;lambda1=lambda;alpha2 =0;beta2=4.2;lambda2=lambda;r1=.5
r2=1-r1# ion 1 delivers fraction r1 of the total dose; ion 2 delivers the rest.
t = seq(0, 1, by = 0.01)# irradiate one week at constant dose-rate
d = R * t# so total dose is 50 cGy
d1 = r1*d; d2 = r2*d# d is split between two ions
E1t <- alpha1*R*t + 2*R^2*beta1/lambda1^2*(lambda1*t -1+exp(-lambda1*t))#integrating the Lea-Catcheside functional explicitly
E2t <- alpha2*R*t + 2*R^2*beta2/lambda2^2*(lambda2*t -1+exp(-lambda2*t))#with constant dose rate gives these IDER

dI1=function(yini,State,Pars){ #incremental effect additivity ordinary differential equation (ODE) initial value problem
  with(as.list(c(State, Pars)), {#uniroot in the next lines determines the inverse functions numerically
    u1=uniroot(function(x) alpha1*R*x + 2*R^2*beta1/lambda1^2*(lambda1*x -1+exp(-lambda1*x)) - I1, lower = 0, upper = 10)$root
    u2=uniroot(function(x) alpha2*R*x + 2*R^2*beta2/lambda2^2*(lambda2*x -1+exp(-lambda2*x)) - I1, lower = 0, upper = 10)$root
    s1 = alpha1*R + 2*R^2*beta1/lambda1*(1-exp(-lambda1*u1))#contribution of ion 1 to the total rate
    s2 = alpha2*R + 2*R^2*beta2/lambda2*(1-exp(-lambda2*u2))
    dI1=r1*s1 + r2*s2# total slope defined as just the weighted sum of individual slopes at each effect value.
    return(list(c(dI1)))
  })
}
pars=c(alpha1,beta1,lambda1,r1, alpha2, beta2,lambda2, r2); yini=c(I1=1e-6); by=.01; times=t#ode is finicky about grammar
out1=ode(yini,times,dI1,pars)#performs the integration
#####simple effect additivity
SEA=alpha1*d1 + 2*R^2*beta1/lambda1^2*(lambda1*d1/R - 1+exp(-lambda1*d1/R))+ #each component contributes the effect due to 
  alpha2*d2 + 2*R^2*beta2/lambda2^2*(lambda2*d2/R - 1+exp(-lambda2*d2/R))    #its own dose
########plot
plot(t, E1t,  col = "blue",type='l', lty=3, lwd=2, ann='F', bty='u')
lines(out1, col="red", ann='F', bty='u')# ,ylim=c(0,alpha2/3))
#lines(d, E1,  col = "blue", lty=3, lwd=2)#sometimes plot out1 and add E1 with lines
lines(t, E2t, col = "blue",lty=3, lwd=2)
lines (t, SEA)
### make the above neater. Should generalize to mixtures with more than 2 ions but not for the first paper where we only ever use 2