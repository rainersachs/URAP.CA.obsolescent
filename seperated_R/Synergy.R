source('/Users/zhaoliyang/Desktop/Seperated_R_File/Datatable.R')

#Our IDERs (Individual Dose Effect Relations). Applicable to the 1-ion components of a mixed simulated GCR beam 
#Modifying NTE1 and NTE2 by insisting they be twice continuously differentiable and monotonic increasing. Double check NTE1, NTE2, Our model

#4-para model
func_4para = function(d, L, Z.b, eta0, eta1, sig0, kap0) {
  P = (1-exp(-Z.b/kap0))^2
  sig = sig0*P + 0.041/6.24*L*(1-P) # 0.041 +- 0.0051 comes from 16Cacao
  eta = eta0*L*exp(-eta1*L)
  0.00071 + sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d))  #0.00071 + sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^3*d))#don't use
} 

#3-para model
func_3para = function(d,L,eta0,eta1,sig0){
  eta = eta0*L*exp(-eta1*L)
  0.00071 + 6.24*sig0*d/L*(1 - exp(-1024*d/L)) + eta*(1-exp(-10^5*d))
}


#2-para model
func_2para=function(d,L,eta0,sig0){
  0.00071+ sig0*6.24*d/L*(1-exp(-1024*d/L)) + eta0*(1-exp(-10^5*d))  
}

#2-paraTE
func_2paraTE = function(d,L,Z.b,kap0,sig0){
  P = (1 - exp(-Z.b/kap0))^2
  sigma = sig0*P + (0.00041 * L/6.242)*(1-P)
  return(0.00071 + (sigma*6.24*(d/L))*(1 - exp(-1024*(d/L))))
}
  

#Calibration for Parsimonious
model_2para = nlsLM(CA ~func_2para(d,L,eta0,sig0) , data = modified_df, start = list(eta0 = 0.001, sig0 = 5), 
                               weights = (1/(modified_df$error)^2))
ccoef = coef(model_2para)
eta0_parsimonious = as.numeric(ccoef[1])
sigm0_parsimonious = as.numeric(ccoef[2])
sig_2para= vcov(model_2para)
#Calibration for 3-parameter model
model_3para = nlsLM(CA~func_3para(d,L,eta0,eta1,sig0), data=modified_df,start = list(eta0 = 0.001, eta1 = 0.01, sig0 = 5),
                         weights = (1/(modified_df$error)^2))
ccoef = coef(model_3para)
eta0_Threepara = as.numeric(ccoef[1])
sigm0_Threepara = as.numeric(ccoef[2])
eta1_Threepara = as.numeric(ccoef[3])
sig_3para = vcov(model_3para)
#Calibration for 4-parameter model
model_4para = nlsLM(CA ~ func_4para(d, L, Z.b, eta0, eta1, sig0, kap0), data = modified_df, start = list(eta0 = 0.001, eta1 = 0.01, sig0 = 5, kap0 = 500), 
                   weights = (1/(modified_df$error)^2))
#Getting coefficients of the IDER model from nlsLM
coefs <- coef(model_4para)
eta0_IDER <- as.numeric(coefs[1])
eta1_IDER <- as.numeric(coefs[2])
sigm0_IDER <- as.numeric(coefs[3])
kap_IDER <- as.numeric(coefs[4])
sig_4para = vcov(model_4para)

##Calibration for 2-parameter TE Only Model
model_2paraTE = nlsLM(CA ~ func_2paraTE(d,L,Z,Z.b,beta,kap0,sig0), data = modified_df, start = list(sig0 = 5, kap0 = 500), 
                   weights = (1/(modified_df$error)^2))
coefs <- coef(model_2paraTE)
sig0_TEonly <- as.numeric(coefs[1])
kap0_TEonly <- as.numeric(coefs[2])
sig_2paraTE = vcov(model_2paraTE)



Data_parameter = data.frame(value=c(eta0_parsimonious,sigm0_parsimonious,eta0_Threepara,sigm0_Threepara,eta1_Threepara,eta0_IDER,
                                    eta1_IDER,sigm0_IDER,kap_IDER,sig0_TEonly,kap0_TEonly),model=c('2para','2para','3para','3para','3para','4para','4para','4para','4para','2paraTE','2paraTE'),
                            parameter = c('eta0','sig0','eta0','sig0','eta1','eta0','eta1','sig0','kap0','sig0','kap0'))
                            
sig_2para
sig_3para
sig_4para
sig_2paraTE

#AIC and BIC
L_function = function(func, eta0=0, eta1=0, sig0=0, kap=0,model){
  #A function trhat returns the residual squared
  if (model == "4para"){
    a = vector(length = 0)
    for (i in 1:length(modified_df[, 1])) {
      a = c(a, modified_df$CA[i] - func(d = modified_df$d[i], L = modified_df$L[i], Z.b = modified_df$Z.b[i], eta0 = eta0, eta1 = eta1, sig0 = sig0, kap = kap))}
  }
  if (model =="3para"){
    a = vector(length = 0)
    for (i in 1:length(modified_df[, 1])) {
      a = c(a, modified_df$CA[i] - func(d = modified_df$d[i], L = modified_df$L[i], Z.b = modified_df$Z.b[i], eta0 = eta0, eta1 = eta1, sig0 = sig0))}
  }
  if(model == "2para"){
    a = vector(length = 0)
    for (i in 1:length(modified_df[, 1])) {
      a = c(a, modified_df$CA[i] - func(d = modified_df$d[i], L = modified_df$L[i], Z.b = modified_df$Z.b[i], eta0 = eta0,sig0 = sig0))}
  }
  if(model =="2paraTE"){
    a = vector(length = 0)
    for (i in 1:length(modified_df[, 1])) {
      a = c(a, modified_df$CA[i] - func(d = modified_df$d[i], L = modified_df$L[i], Z.b = modified_df$Z.b[i], sig0 = sig0, kap = kap))}
  }
  return(a^2)
}

L_2para = L_function(func_2para, eta0 = Data_parameter$value[1],sig0 = Data_parameter$value[2],model ="2para")
L_3para = L_function(func_3para, eta0 = Data_parameter$value[3], eta1 = Data_parameter$value[5], sig0 = Data_parameter$value[4],model="3para")
L_4para  = L_function(func_4para, eta0 = Data_parameter$value[6], eta1 = Data_parameter$value[7], sig0 = Data_parameter$value[8], kap = Data_parameter$value[9],model = "4para")
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
#Confidence Intervals (Monte Carlo)
set.seed(19970101)
MM <- 500
#Sample parameters from their distributions with covariances for 4para
monte_carlo_parameters = rmvnorm(n = MM, mean = c(eta0 = Data_parameter[6,1], eta1 = Data_parameter[7,1], sig0 = Data_parameter[8,1], kap = Data_parameter[9,1]), sigma = sig_4para)
eta0_MC = monte_carlo_parameters[, 1]
eta1_MC = monte_carlo_parameters[, 2]
sig0_MC = monte_carlo_parameters[, 3]
kap_MC = monte_carlo_parameters[, 4]
kap_MC[kap_MC <= 1e-6] = 1e-5 #10 of them were fixed

#Sample parameters from their distributions without covariances
eta0_MC_var = rnorm(MM, mean = Data_parameter[6,1], sd = sqrt(sig_4para[1,1]))
eta1_MC_var = rnorm(MM, mean = Data_parameter[7,1], sd = sqrt(sig_4para[2,2]))
sig0_MC_var = rnorm(MM, mean = Data_parameter[8,1], sd = sqrt(sig_4para[3,3]))
kap_MC_var = rnorm(MM, mean = Data_parameter[9,1], sqrt(sig_4para[4,4]))
kap_MC_var[kap_MC_var <= 1e-6] = 1e-5 


para_4para_nocov =  c(list(eta0_MC),list(eta1_MC),list(sig0_MC),list(kap_MC))
para_4para_cov = c(list(eta0_MC_var),list(eta1_MC_var),list(sig0_MC_var),list(kap_MC_var))

#Sample parameters from their distributions with covariances for Three parameter model.
monte_carlo_parameters = rmvnorm(n = MM, mean = c(eta0 = Data_parameter[3,1], sig0 = Data_parameter[4,1], eta1 = Data_parameter[5,1]), sigma = sig_3para)
eta0_MC = monte_carlo_parameters[, 1]
sig0_MC = monte_carlo_parameters[, 2]
eta1_MC = monte_carlo_parameters[, 3]

#Sample parameters from their distributions without covariances for Three parameter model.
eta0_MC_var = rnorm(MM, mean = Data_parameter[3,1], sd = sqrt(sig_3para[1,1]))
sig0_MC_var = rnorm(MM, mean = Data_parameter[4,1], sd = sqrt(sig_3para[2,2]))
eta1_MC_var = rnorm(MM, mean = Data_parameter[5,1], sd = sqrt(sig_3para[3,3]))
para_3para_nocov =  c(list(eta0_MC),list(sig0_MC),list(eta1_MC))
para_3para_cov = c(list(eta0_MC_var),list(sig0_MC_var),list(eta1_MC_var))

#Sample parameters from their distributions with covariances for 2para
monte_carlo_parameters = rmvnorm(n = MM, mean = c(eta0 = Data_parameter[1,1], sig0 = Data_parameter[2,1], sigma = sig_2para))
eta0_MC = monte_carlo_parameters[, 1]
sig0_MC = monte_carlo_parameters[, 2]


#Sample parameters from their distributions without covariances for 2para

eta0_MC_var = rnorm(MM, mean = Data_parameter[1,1], sd = sqrt(sig_2para[1,1]))
sig0_MC_var = rnorm(MM, mean = Data_parameter[2,1], sd = sqrt(sig_2para[2,2]))
para_2para_nocov =  c(list(eta0_MC),list(sig0_MC))
para_2para_cov= c(list(eta0_MC_var),list(sig0_MC_var))

#Sample parameters from their distributions with covariances for 2para TE only
monte_carlo_parameters = rmvnorm(n = MM, mean = c(eta0 = Data_parameter[9,1], sig0 = Data_parameter[10,1], sigma = sig_2paraTE))
sig0_MC = monte_carlo_parameters[, 1]
kap0_MC = monte_carlo_parameters[, 2]
#Sample parameters from their distributions without covariances for 2para TE only
sig0_MC_var = rnorm(MM, mean = Data_parameter[9,1], sd = sqrt(sig_2paraTE[1,1]))
kap0_MC_var = rnorm(MM, mean = Data_parameter[10,1], sd = sqrt(sig_2paraTE[2,2]))
para_2paraTE_nocov =  c(list(sig0_MC),list(kap0_MC))
para_2paraTE_cov= c(list(sig0_MC_var),list(kap0_MC_var))


#Useful output
Data_parameter
sig_2para
sig_3para
sig_4para
sig_2paraTE
#Monte Carlo
para_4para_nocov
para_4para_cov
para_3para_nocov
para_3para_cov
para_2para_nocov
para_2para_cov
para_2paraTE_nocov
para_2paraTE_cov