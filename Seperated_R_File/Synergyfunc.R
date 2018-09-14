#First import all the dataframes generated from csv files.
source('/Users/zhaoliyang/Desktop/URAP\ FALL2018/Datatable.R')

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