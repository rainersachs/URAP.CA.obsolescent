source('Datatable.R')

#Our IDERs (Individual Dose Effect Relations). Applicable to the 1-ion components of a mixed simulated GCR beam 
#Modifying NTE1 and NTE2 by insisting they be twice continuously differentiable and monotonic increasing. Double check NTE1, NTE2, Our model

#4-para model
func_4para = function(d, L, Z.b, eta0, eta1, sig0, kap0) {
  P = (1-exp(-Z.b/kap0))^2
  sig = sig0*P + 0.041/6.24*L*(1-P) # 0.041 +- 0.0051 comes from 16Cacao
  eta = eta0*L*exp(-eta1*L)
  return(0.00071 + sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d)))  #0.00071 + sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^3*d))#don't use
} 

#3-para model
func_3para = function(d,L,eta0,eta1,sig0){
  eta = eta0*L*exp(-eta1*L)
  return(0.00071 + 6.24*sig0*d/L*(1 - exp(-1024*d/L)) + eta*(1-exp(-10^5*d)))
}


#2-para model
func_2para=function(d,L,eta0,sig0){
  return(0.00071 + sig0*6.24*d/L*(1-exp(-1024*d/L)) + eta0*(1-exp(-10^5*d)))
}

#2-paraTE
func_2paraTE = function(d,L,Z.b,kap0,sig0){
  P = (1 - exp(-Z.b/kap0))^2
  sigma = sig0*P + (0.041 * L/6.24)*(1-P)
  return(0.00071 + (sigma*6.24*(d/L))*(1 - exp(-1024*(d/L))))
}
  
modified_df = modified_df %>% filter(Z > 3)
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
model_2paraTE = nlsLM(CA ~ func_2paraTE(d,L,Z.b,kap0,sig0), data = modified_df, start = list(sig0 = 5, kap0 = 500), 
                   weights = (1/(modified_df$error)^2))
coefs <- coef(model_2paraTE)
sig0_TEonly <- as.numeric(coefs[1])
kap0_TEonly <- as.numeric(coefs[2])
sig_2paraTE = vcov(model_2paraTE)



Data_parameter = data.frame(value=c(eta0_parsimonious,sigm0_parsimonious,eta0_Threepara,sigm0_Threepara,eta1_Threepara,eta0_IDER,
                                    eta1_IDER,sigm0_IDER,kap_IDER,sig0_TEonly,kap0_TEonly),model=c('2para','2para','3para','3para','3para','4para','4para','4para','4para','2paraTE','2paraTE'),
                            parameter = c('eta0','sig0','eta0','sig0','eta1','eta0','eta1','sig0','kap0','sig0','kap0'))


#Useful output
Data_parameter
sig_2para
sig_3para
sig_4para
sig_2paraTE


########################################
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
      a = c(a, modified_df$CA[i] - func(d = modified_df$d[i], L = modified_df$L[i], eta0 = eta0, eta1 = eta1, sig0 = sig0))}
  }
  if(model == "2para"){
    a = vector(length = 0)
    for (i in 1:length(modified_df[, 1])) {
      a = c(a, modified_df$CA[i] - func(d = modified_df$d[i], L = modified_df$L[i],eta0 = eta0,sig0 = sig0))}
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
L_2paraTE = L_function(func_2paraTE,sig0 = Data_parameter$value[10], kap = Data_parameter$value[11],model = "2paraTE")

WRSS_2para = sum((1/modified_df$error^2)*L_2para)
WRSS_3para = sum((1/modified_df$error^2)*L_3para)
WRSS_4para = sum((1/modified_df$error^2)*L_4para)
WRSS_2paraTE = sum((1/modified_df$error^2)*L_2paraTE)

#functions for AIC and BIC calculation for Weighted Least Square regression (using WRSS calculated above)
AIC_function = function(RSS, k, n = length(modified_df[ , 1])) {
  n + n*log(2*pi) + n*log(RSS/n) + 2*(k+1)
}

BIC_function = function(n = length(modified_df[, 1]), k, RSS) {
  n + n*log(2*pi) + n*log(RSS/n) + log(n)*(k+1)
}

AIC_2para = AIC_function(k=2, RSS = WRSS_2para)
BIC_2para = BIC_function(k=2, RSS = WRSS_2para)
AIC_3para = AIC_function(k=3, RSS = WRSS_3para)
BIC_3para = BIC_function(k=3, RSS = WRSS_3para)
AIC_4para = AIC_function(k=4, RSS = WRSS_4para)
BIC_4para = BIC_function(k=4, RSS = WRSS_4para)
AIC_2paraTE = AIC_function(k=2, RSS = WRSS_2paraTE)
BIC_2paraTE = BIC_function(k=2, RSS = WRSS_2paraTE)
information_critera_df = data.frame(AIC = c(AIC_2para,AIC_3para,AIC_4para,AIC_2paraTE), BIC= c(BIC_2para,BIC_3para,BIC_4para,BIC_2paraTE), row.names = c("2para model", "3para model", "4para model", "2paraTE model"))
information_critera_df #IDER has the lowest score (performs better) in both criteria.

#Leave-one-out cross validation
LOO = function(model,modified_df,func){
  value <- data.frame(matrix(nrow=length(modified_df[, 1]), ncol=2)) 
  x <- c("CA_value", "CA_predict")
  colnames(value) <- x
  if (model == "2para"){
    for (i in 1:length(modified_df[, 1])){
      train = modified_df[-i,]
      test = modified_df[i,]
      model_2para = nlsLM(CA ~func(d,L,eta0,sig0) , data = train, start = list(eta0 = 0.001, sig0 = 5), 
                          weights = (1/(train$error)^2))
      ccoef = coef(model_2para)
      eta0_2para = as.numeric(ccoef[1])
      sigm0_2para = as.numeric(ccoef[2])
      value$CA_value[i] = test$CA
      value$CA_predict[i] = func(test$d,test$L,eta0_2para,sigm0_2para)}
  }
  if (model == "3para"){
    for (i in 1:length(modified_df[, 1])){
      train = modified_df[-i,]
      test = modified_df[i,]
      model_3para = nlsLM(CA~func_3para(d,L,eta0,eta1,sig0), data=modified_df,start = list(eta0 = 0.001, eta1 = 0.01, sig0 = 5),
                          weights = (1/(modified_df$error)^2))
      ccoef = coef(model_3para)
      eta0_Threepara = as.numeric(ccoef[1])
      sigm0_Threepara = as.numeric(ccoef[2])
      eta1_Threepara = as.numeric(ccoef[3])
      value$CA_value[i] = test$CA
      value$CA_predict[i] = func(test$d,test$L,eta0_Threepara,sigm0_Threepara,eta1_Threepara)}
  }
  if (model == "4para"){
    for (i in 1:length(modified_df[, 1])){
      train = modified_df[-i,]
      test = modified_df[i,]
      model_4para = nlsLM(CA ~ func_4para(d, L, Z.b, eta0, eta1, sig0, kap0), data = modified_df, start = list(eta0 = 0.001, eta1 = 0.01, sig0 = 5, kap0 = 500), 
                          weights = (1/(modified_df$error)^2))
      coefs <- coef(model_4para)
      eta0_IDER <- as.numeric(coefs[1])
      eta1_IDER <- as.numeric(coefs[2])
      sigm0_IDER <- as.numeric(coefs[3])
      kap_IDER <- as.numeric(coefs[4])
      value$CA_value[i] = test$CA
      value$CA_predict[i] = func(test$d,test$L,test$Z.b,eta0_IDER,eta1_IDER,sigm0_IDER,kap_IDER)}
  }
  if (model == "2paraTE"){
    for (i in 1:length(modified_df[, 1])){
      train = modified_df[-i,]
      test = modified_df[i,]
      model_2paraTE = nlsLM(CA ~ func_2paraTE(d,L,Z.b,kap0,sig0), data = modified_df, start = list(sig0 = 5, kap0 = 500), 
                            weights = (1/(modified_df$error)^2))
      coefs <- coef(model_2paraTE)
      sig0_TEonly <- as.numeric(coefs[1])
      kap0_TEonly <- as.numeric(coefs[2])
      value$CA_value[i] = test$CA
      value$CA_predict[i] = func(test$d,test$L,test$Z.b,sig0_TEonly,kap0_TEonly)}
  }
  return(value)
}

LOO_2para_df = LOO("2para",modified_df,func_2para)
LOO_2para = (1/length(modified_df[, 1]))*sum((LOO_2para_df$CA_value - LOO_2para_df$CA_predict)^2)
LOO_3para_df = LOO("3para",modified_df,func_3para)
LOO_3para = (1/length(modified_df[, 1]))*sum((LOO_3para_df$CA_value - LOO_3para_df$CA_predict)^2)
LOO_4para_df = LOO("4para",modified_df,func_4para)
LOO_4para = (1/length(modified_df[, 1]))*sum((LOO_4para_df$CA_value - LOO_4para_df$CA_predict)^2)
LOO_2paraTE_df = LOO("2paraTE",modified_df,func_2paraTE)
LOO_2paraTE = (1/length(modified_df[, 1]))*sum((LOO_2paraTE_df$CA_value - LOO_2paraTE_df$CA_predict)^2)
LOO_CV_df = data.frame(CV_value = c(LOO_2para,LOO_3para,LOO_4para,LOO_2paraTE), row.names = c("2para model", "3para model", "4para model", "2paraTE model"))
LOO_CV_df
########################################

IDER = function(d, L = NULL, Z.b = NULL, ions = NULL, r = NULL, parameters = Data_parameter, model = "4para") {
  #Main function that outputs the Individual dose-effect relationship. If "ions" are specified, take the average of them instead. 
  #Right now this function requires that the calibrated parameters to be defined in the exact names as above. Need to work on this part.
  
  name = model #To avoid "model" being both a string (model name) and a column name
  parameters = parameters %>% filter(model == name) #Only the parameters we want
  
  if (is.null(ions)){
    if (model == "4para"){ #Parameters are kap, sig0, eta0, eta1
      kap = (parameters %>% filter(parameter == "kap0"))$value
      sig0 = (parameters %>% filter(parameter == "sig0"))$value
      eta0 = (parameters %>% filter(parameter == "eta0"))$value
      eta1 = (parameters %>% filter(parameter == "eta1"))$value
      
      P = (1-exp(-Z.b/kap))^2
      sig = sig0*P + 0.041/6.24*L*(1-P)
      eta = eta0*L*exp(-eta1*L)
      return(sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d)))
    }
    if (model == "3para"){ #Parameters are sig0, eta0, eta1
      sig0 = (parameters %>% filter(parameter == "sig0"))$value
      eta0 = (parameters %>% filter(parameter == "eta0"))$value
      eta1 = (parameters %>% filter(parameter == "eta1"))$value
      
      eta = eta0*L*exp(-eta1*L)
      return(6.24*sig0*d/L*(1 - exp(-1024*d/L)) + eta*(1-exp(-10^5*d)))
    }
    if (model == "2para"){ #Parameters are sig0, eta0
      sig0 = (parameters %>% filter(parameter == "sig0"))$value
      eta0 = (parameters %>% filter(parameter == "eta0"))$value
      
      return(sig0*6.24*d/L*(1-exp(-1024*d/L)) + eta0*(1-exp(-10^5*d)))
    }
    if (model == "2paraTE"){ #Parameters are kap, sig0
      kap = (parameters %>% filter(parameter == "kap0"))$value
      sig0 = (parameters %>% filter(parameter == "sig0"))$value
      
      P = (1-exp(-Z.b/kap))^2
      sigma = sig0*P + (0.00041 * L/6.242)*(1-P)
      return((sigma*6.24*(d/L))*(1 - exp(-1024*(d/L))))
    }
    # if (model == "lowLET"){
    #   return(0.00001+ alpha * d)
    # }
  }
  
  
  ########### SEA 
  if (is.null(r)){#Default r is average
    r = 1/length(ions)
  }
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
MIXDER_function = function(r, L, Z.b, d = seq(0, 0.2, by = 0.001), parameters = Data_parameter, model = "4para") {
  name = model
  parameters = parameters %>% filter(model == name)
  if (model == "4para"){
    kap = (parameters %>% filter(parameter == "kap0"))$value
    sig0 = (parameters %>% filter(parameter == "sig0"))$value
    eta0 = (parameters %>% filter(parameter == "eta0"))$value
    eta1 = (parameters %>% filter(parameter == "eta1"))$value
    
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
    sig0 = (parameters %>% filter(parameter == "sig0"))$value
    eta0 = (parameters %>% filter(parameter == "eta0"))$value
    eta1 = (parameters %>% filter(parameter == "eta1"))$value
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
    sig0 = (parameters %>% filter(parameter == "sig0"))$value
    eta0 = (parameters %>% filter(parameter == "eta0"))$value
    dE = function(yini, State, Pars){
      sig0 = sig0; eta0 = eta0
      with(as.list(c(State, Pars)), {
        u = vector(length = length(L))
        for (i in 1:length(L)) {
          u[i] = uniroot(function(d) sig0*6.24*d/L[i]*(1-exp(-1024*d/L[i])) + eta0*(1-exp(-10^5*d)) - I, lower = 0, upper = 1, extendInt = "yes", tol = 10^-10)$root
        }
        dI = vector(length = length(L))
        for (i in 1:length(L)) {
          dI[i] = r[i]*(sig0*6.24/L[i]*exp(-1024*u[i]/L[i])*(exp(1024*u[i]/L[i]) + 1024*u[i]/L[i] - 1) + eta0*10^5*exp(-10^5*u[i]))
        }
        dI = sum(dI)
        return(list(c(dI)))
      })
    }
  }
  
  else if (model == "2paraTE"){
    kap = (parameters %>% filter(parameter == "kap0"))$value
    sig0 = (parameters %>% filter(parameter == "sig0"))$value
    dE=function(yini,State,Pars){
      kap = kap; sig0 = sig0
      with(as.list(c(State, Pars)), {
        P = vector(length = length(L))
        sig = vector(length = length(L))
        u = vector(length = length(L))
        for (i in 1:length(L)) {
          P[i] = (1-exp(-Z.b[i]/kap))^2
          sig[i] = sig0*P[i] + 0.041/6.24*L[i]*(1-P[i])
          u[i] = uniroot(function(d) sig[i]*6.24*d/L[i]*(1-exp(-1024*d/L[i])) - I, lower = 0, upper = 1, extendInt = "yes", tol = 10^-10)$root
        }
        dI = vector(length = length(L))
        for (i in 1:length(L)) {
          dI[i] = r[i]*(sig[i]*6.24/L[i]*(1 - exp(-1024*u[i]/L[i]) + 1024*u[i]/L[i]*(exp(-1024*u[i]/L[i]))) + etaa[i]*10^5*exp(-10^5*u[i]))
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
