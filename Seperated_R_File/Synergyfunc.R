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


#Parameter Calibration (Andy is in charge here)

###########################


#Individual Dose-Effect Relationship
#Models include: 4para, 3para, 2para, 2paraTE, lowLET
IDER = function(d, L = NULL, Z.b = NULL, ions = NULL, eta0 = eta0_hat, eta1 = eta1_hat, sig0 = sig0_hat, kap = kap_hat, sigm_0 = sigm_0_hat, eta_0 = eta_0_hat, model = "4para", alpha = alpha_hat) {
  #Main function that outputs the Individual dose-effect relationship. If "ions" are specified, take the average of them instead. 
  #Right now this function requires that the calibrated parameters to be defined in the exact names as above. Need to work on this part.
  if (is.null(ions)){
    if (model == "4para"){
      P = (1-exp(-Z.b/kap_hat))^2
      sig = sig0*P + 0.041/6.24*L*(1-P)
      eta = eta0*L*exp(-eta1*L)
      return(sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d)))
    }
    if (model == "3para"){
      eta = eta0*L*exp(-eta1*L)
      return(6.24*sig0*d/L*(1 - exp(-1024*d/L)) + eta*(1-exp(-10^5*d)))
    }
    if (model == "2para"){
      return(sigm_0*6.24*d/L*(1-exp(-1024*d/L)) + eta_0*(1-exp(-10^5*d)))
    }
    if (model == "2paraTE"){
      P = (1 - exp(-Z^2/(Z.b^2)))^2
      sigma = te_sig0*P + (0.00041 * L/6.242)*(1-P)
      return((sigma*6.24*(d/L))*(1 - exp(-1024*(d/L))))
    }
    if (model == "lowLET")
      return(0.00001+ alpha * d)
  }
  r = 1/length(ions)
  info_table = modified_df %>% group_by(ion, L, Z.b) %>% summarise()
  info_table = suppressWarnings(left_join(data.frame(ions), info_table, by = c("ions" = "ion")))
  output = 0
  for (i in 1: nrow(info_table)){
    output = output + IDER(d*r, L = info_table$L[i], Z.b = info_table$Z.b[i], model = model)
  }
  return(output)
  stop("Model Not Recognized")
}

#Mixture
MIXDER_function = function(r, L, Z.b, d = seq(0, 0.2, by = 0.001), eta0 = eta0_hat, eta1 = eta1_hat, sig0 = sig0_hat ,kap = kap_hat, sigm_0 = sigm_0_hat, eta_0 = eta_0_hat, model = "4para") {
  if (model == "4para"){
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
  }
  else if (model == "3para"){
    dE=function(yini,State,Pars){
      eta0 = eta0; eta1 = eta1; sig0 = sig0
      with(as.list(c(State, Pars)), {
        sig = vector(length = length(L))
        etaa = vector(length = length(L))
        u = vector(length = length(L))
        for (i in 1:length(L)) {
          etaa[i] = eta0*L[i]*exp(-eta1*L[i])
          u[i] = uniroot(function(d) sig0*6.24*d/L[i]*(1-exp(-1024*d/L[i])) + etaa[i]*(1-exp(-10^5*d)) - I, lower = 0, upper = 1, extendInt = "yes", tol = 10^-10)$root
        }
        dI = vector(length = length(L))
        for (i in 1:length(L)) {
          dI[i] = r[i]*(sig0*6.24/L[i]*exp(-1024*u[i]/L[i])*(exp(1024*u[i]/L[i]) + 1024*u[i]/L[i] - 1) + etaa[i]*10^5*exp(-10^5*u[i]))
        }
        dI = sum(dI)
        return(list(c(dI)))
      })
    }
  }
  else if (model == "2para"){
    dE = function(yini, State, Pars){
      sigm_0 = sigm_0; eta_0 = eta_0
      with(as.list(c(State, Pars)), {
        u = vector(length = length(L))
        for (i in 1:length(L)) {
          u[i] = uniroot(function(d) sigm_0*6.24*d/L[i]*(1-exp(-1024*d/L[i])) + eta_0*(1-exp(-10^5*d)) - I, lower = 0, upper = 1, extendInt = "yes", tol = 10^-10)$root
        }
        dI = vector(length = length(L))
        for (i in 1:length(L)) {
          dI[i] = r[i]*(sigm_0*6.24/L[i]*exp(-1024*u[i]/L[i])*(exp(1024*u[i]/L[i]) + 1024*u[i]/L[i] - 1) + eta_0*10^5*exp(-10^5*u[i]))
        }
        dI = sum(dI)
        return(list(c(dI)))
      })
    }
  }
  else{
    stop("Model Not Recognized")
  }
  pars = NULL; yini = c(I= 0); d = d
  out = ode(yini,times = d, dE, pars, method = "radau")
  return(out)
}
