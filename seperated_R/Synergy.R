source('Datatable.R')

#Our IDERs (Individual Dose Effect Relations). Applicable to the 1-ion components of a mixed simulated GCR beam 
#Modifying NTE1 and NTE2 by insisting they be twice continuously differentiable and monotonic increasing. Double check NTE1, NTE2, Our model
library("forecast")

#4-para model
func_4para = function(d, L, Z.b, eta0, eta1, sig0, kap0) {
  P = (1-exp(-Z.b/kap0))^2
  sig = sig0*P + 0.041/6.24*L*(1-P) # 0.041 +- 0.0051 comes from 16Cacao
  eta = eta0*L*exp(-eta1*L)
  return(BG_CA + sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d)))  #0.00071 + sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^3*d))#don't use
} 

#3-para model (##PW: obsolete)
func_3para = function(d,L,eta0,eta1,sig0){
  eta = eta0*L*exp(-eta1*L)
  return(BG_CA + 6.24*sig0*d/L*(1 - exp(-1024*d/L)) + eta*(1-exp(-10^5*d)))
}


#2-para model (##PW: obsolete)
func_2para=function(d,L,eta0,sig0){
  return(BG_CA + sig0*6.24*d/L*(1-exp(-1024*d/L)) + eta0*(1-exp(-10^5*d)))
}

#2-paraTE (##PW: Is this the "straw man TE-only 2 parameter model" Sachs referred to?)
func_2paraTE = function(d,L,Z.b,kap0,sig0){
  P = (1 - exp(-Z.b/kap0))^2
  sigma = sig0*P + (0.041 * L/6.24)*(1-P)
  return(BG_CA + (sigma*6.24*(d/L))*(1 - exp(-1024*(d/L))))
}

#Swift Light Ion 2 parameter (less parsimonious) 
func_sli_exp2 <- function(a, C, d) {
  return(BG_CA + C*(exp(a*d)-1))
}

#Swift Light Ion 3 paramter (less parsimonious)
func_sli_exp3 <- function(a, b, C, d) {
  return(BG_CA + C*(exp(a*d+b*d^2)-1))
}


swift_light_df <- filter(main_df, Z <= 2 & d>0)
main_df = main_df %>% filter(Z > 3)
main_df = main_df %>% filter(d > 0)

#Calibration for Parsimonious
model_2para = nlsLM(CA ~func_2para(d,L,eta0,sig0) , data = main_df, start = list(eta0 = 0.001, sig0 = 5), 
                    weights = (1/(main_df$error)^2))
ccoef = coef(model_2para)
parameters_2para = data.frame(value = as.numeric(ccoef), model = "2para", parameter = names(ccoef))
sig_2para= vcov(model_2para)
#Calibration for 3-parameter model
model_3para = nlsLM(CA~func_3para(d,L,eta0,eta1,sig0), data=main_df,start = list(eta0 = 0.001, eta1 = 0.01, sig0 = 5),
                    weights = (1/(main_df$error)^2))
ccoef = coef(model_3para)
parameters_3para = data.frame(value = as.numeric(ccoef), model = "3para", parameter = names(ccoef))

sig_3para = vcov(model_3para)
#Calibration for 4-parameter model
model_4para = nlsLM(CA ~ func_4para(d, L, Z.b, eta0, eta1, sig0, kap0), data = main_df, start = list(eta0 = 0.001, eta1 = 0.01, sig0 = 5, kap0 = 500), 
                    weights = (1/(main_df$error)^2))
#Getting coefficients of the IDER model from nlsLM
ccoef <- coef(model_4para)
parameters_4para = data.frame(value = as.numeric(ccoef), model = "4para", parameter = names(ccoef))

sig_4para = vcov(model_4para)

##Calibration for 2-parameter TE Only Model
model_2paraTE = nlsLM(CA ~ func_2paraTE(d,L,Z.b,kap0,sig0), data = main_df, start = list(sig0 = 5, kap0 = 500), 
                      weights = (1/(main_df$error)^2))
ccoef <- coef(model_2paraTE)
parameters_2paraTE = data.frame(value = as.numeric(ccoef), model = "2paraTE", parameter = names(ccoef))

sig_2paraTE = vcov(model_2paraTE)

swift_light_df <- swift_light_df %>% mutate(dd = d^2)

#Calibration for SLI Linear 1 parameter
model_sli_linear1 <- lm(I(CA - BG_CA) ~ 0 + d, data = swift_light_df, weights = (1/(swift_light_df$error)^2))
ccoef <- coef(model_sli_linear1)
parameters_sli_lin1 <- data.frame(value = as.numeric(ccoef), model = "sli_lin1", parameter = names(ccoef))
sig_lin1 <- vcov(model_sli_linear1)

#Calibration for SLI Linear 2 parameter
model_sli_linear2 <- lm(I(CA - BG_CA) ~ 0 + d + I(d^2), data = swift_light_df, weights = (1/(swift_light_df$error)^2))
ccoef <- coef(model_sli_linear2)
parameters_sli_lin2 <- data.frame(value = as.numeric(ccoef), model = "sli_lin2", parameter = names(ccoef))
sig_lin2 <- vcov(model_sli_linear2)

#lets take out values below bg
edit_sli <- filter(swift_light_df, CA >= BG_CA)
#calibration exponential
model_sli_exp2 <- nlsLM(CA ~ func_sli_exp2(a, C, d), data = swift_light_df, start = list(a = -0.01, C = -1), weights = (1/(swift_light_df$error)^2))
ccoef <- coef(model_sli_exp2)
parameters_sli_exp2 <- data.frame(value = as.numeric(ccoef), model = "sli_exp2", parameter = names(ccoef))
sig_exp2 <- vcov(model_sli_exp2)


model_sli_exp3 <-  nlsLM(CA ~ func_sli_exp3(a,b, C, d), data = swift_light_df, start = list(a = -0.5, C = -1, b = -10), weights = (1/(swift_light_df$error)^2))
ccoef <- coef(model_sli_exp3)
parameters_sli_exp3 <- data.frame(value = as.numeric(ccoef), model = "sli_exp3", parameter = names(ccoef))
sig_exp3 <- vcov(model_sli_exp3)
#ggplot(swift_light_df) + geom_point(aes(x = d, y = CA)) + geom_line(aes(x = d, y = predict(model_sli_exp2)), col = "blue") + geom_line(aes(x = d, y = predict(model_sli_exp3)), col = "red")


#Combining the calibrated parameters into one dataframe.
# Data_parameter = data.frame(value=c(eta0_parsimonious,sigm0_parsimonious,eta0_Threepara,sigm0_Threepara,eta1_Threepara,eta0_IDER,
#                                     eta1_IDER,sigm0_IDER,kap_IDER,sig0_TEonly,kap0_TEonly),model=c('2para','2para','3para','3para','3para','4para','4para','4para','4para','2paraTE','2paraTE'),
#                             parameter = c('eta0','sig0','eta0','sig0','eta1','eta0','eta1','sig0','kap0','sig0','kap0'))
#The code above is commented out because the values are typed in pretty much manually. Imrpoved the code to make it not sensible to number of parameters.
Data_parameter = rbind(parameters_2para, 
                       parameters_3para, 
                       parameters_4para, 
                       parameters_2paraTE, 
                       parameters_sli_lin1, 
                       parameters_sli_lin2, 
                       parameters_sli_exp2, 
                       parameters_sli_exp3)

#Useful output
Data_parameter
sig_2para
sig_3para
sig_4para
sig_2paraTE
sig_lin1
sig_lin2


########################################
#AIC and BIC
L_function = function(func, eta0=0, eta1=0, sig0=0, kap=0,model){
  #A function trhat returns the residual squared
  if (model == "4para"){
    a = vector(length = 0)
    for (i in 1:length(main_df[, 1])) {
      a = c(a, main_df$CA[i] - func(d = main_df$d[i], L = main_df$L[i], Z.b = main_df$Z.b[i], eta0 = eta0, eta1 = eta1, sig0 = sig0, kap = kap))}
  }
  if (model =="3para"){
    a = vector(length = 0)
    for (i in 1:length(main_df[, 1])) {
      a = c(a, main_df$CA[i] - func(d = main_df$d[i], L = main_df$L[i], eta0 = eta0, eta1 = eta1, sig0 = sig0))}
  }
  if(model == "2para"){
    a = vector(length = 0)
    for (i in 1:length(main_df[, 1])) {
      a = c(a, main_df$CA[i] - func(d = main_df$d[i], L = main_df$L[i],eta0 = eta0,sig0 = sig0))}
  }
  if(model =="2paraTE"){
    a = vector(length = 0)
    for (i in 1:length(main_df[, 1])) {
      a = c(a, main_df$CA[i] - func(d = main_df$d[i], L = main_df$L[i], Z.b = main_df$Z.b[i], sig0 = sig0, kap = kap))}
  }
  return(a^2)
}

L_2para = L_function(func_2para, eta0 = filter(parameters_2para, parameter == "eta0")$value ,sig0 = filter(parameters_2para, parameter == "sig0")$value,model ="2para")
L_3para = L_function(func_3para, eta0 = filter(parameters_3para, parameter == "eta0")$value, eta1 = filter(parameters_3para, parameter == "eta1")$value, sig0 = filter(parameters_3para, parameter == "sig0")$value, model = "3para")
L_4para  = L_function(func_4para, eta0 = filter(parameters_4para, parameter == "eta0")$value, eta1 = filter(parameters_4para, parameter == "eta1")$value, sig0 = filter(parameters_4para, parameter == "sig0")$value, kap = filter(parameters_4para, parameter == "kap0")$value,model = "4para")
L_2paraTE = L_function(func_2paraTE,sig0 = filter(parameters_2paraTE, parameter == "sig0")$value, kap = filter(parameters_2paraTE, parameter == "kap0")$value,model = "2paraTE")

WRSS_2para = sum((1/main_df$error^2)*L_2para)
WRSS_3para = sum((1/main_df$error^2)*L_3para)
WRSS_4para = sum((1/main_df$error^2)*L_4para)
WRSS_2paraTE = sum((1/main_df$error^2)*L_2paraTE)

#functions for AIC and BIC calculation for Weighted Least Square regression (using WRSS calculated above)
AIC_function = function(RSS, k, n = nrow(main_df)) {
  n + n*log(2*pi) + n*log(RSS/n) + 2*(k+1)
}

BIC_function = function(n = nrow(main_df), k, RSS) {
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
information_criteria_df = data.frame(AIC = c(AIC_2para,AIC_3para,AIC_4para,AIC_2paraTE), BIC= c(BIC_2para,BIC_3para,BIC_4para,BIC_2paraTE), row.names = c("2para model", "3para model", "4para model", "2paraTE model"))
information_criteria_df #IDER has the lowest score (performs better) in both criteria.

#Leave-one-out cross validation
LOO = function(model,main_df,func){
  value <- data.frame(matrix(nrow=nrow(main_df), ncol=2)) 
  x <- c("CA_value", "CA_predict")
  colnames(value) <- x
  if (model == "2para"){
    for (i in 1:nrow(main_df)){
      train = main_df[-i,]
      test = main_df[i,]
      model_2para = nlsLM(CA ~func(d,L,eta0,sig0) , data = train, start = list(eta0 = 0.001, sig0 = 5), 
                          weights = (1/(train$error)^2))
      ccoef = coef(model_2para)
      eta0_2para = as.numeric(ccoef[1])
      sigm0_2para = as.numeric(ccoef[2])
      value$CA_value[i] = test$CA
      value$CA_predict[i] = func(test$d,test$L,eta0_2para,sigm0_2para)}
  }
  if (model == "3para"){
    for (i in 1:nrow(main_df)){
      train = main_df[-i,]
      test = main_df[i,]
      model_3para = nlsLM(CA~func_3para(d,L,eta0,eta1,sig0), data=train,start = list(eta0 = 0.001, eta1 = 0.01, sig0 = 5),
                          weights = (1/(train$error)^2))
      ccoef = coef(model_3para)
      eta0_Threepara = as.numeric(ccoef[1])
      sigm0_Threepara = as.numeric(ccoef[3])
      eta1_Threepara = as.numeric(ccoef[2])
      value$CA_value[i] = test$CA
      value$CA_predict[i] = func(test$d,test$L,eta0_Threepara,eta1_Threepara, sigm0_Threepara)}
  }
  if (model == "4para"){
    for (i in 1:nrow(main_df)){
      train = main_df[-i,]
      test = main_df[i,]
      model_4para = nlsLM(CA ~ func_4para(d, L, Z.b, eta0, eta1, sig0, kap0), data = train, start = list(eta0 = 0.0001, eta1 = 0.001, sig0 = 5, kap0 = 500), 
                          weights = (1/(train$error)^2))
      coefs <- coef(model_4para)
      eta0_IDER <- as.numeric(coefs[1])
      eta1_IDER <- as.numeric(coefs[2])
      sigm0_IDER <- as.numeric(coefs[3])
      kap_IDER <- as.numeric(coefs[4])
      value$CA_value[i] = test$CA
      value$CA_predict[i] = func(test$d,test$L,test$Z.b,eta0_IDER,eta1_IDER,sigm0_IDER,kap_IDER)}
  }
  if (model == "2paraTE"){
    for (i in 1:nrow(main_df)){
      train = main_df[-i,]
      test = main_df[i,]
      model_2paraTE = nlsLM(CA ~ func_2paraTE(d,L,Z.b,kap0,sig0), data = train, start = list(sig0 = 5, kap0 = 500), 
                            weights = (1/(train$error)^2))
      coefs <- coef(model_2paraTE)
      sig0_TEonly <- as.numeric(coefs[1])
      kap0_TEonly <- as.numeric(coefs[2])
      value$CA_value[i] = test$CA
      value$CA_predict[i] = func(test$d,test$L,test$Z.b,kap0_TEonly, sig0_TEonly)}
  }
  return(value)
}

LOO_2para_df = LOO("2para",main_df,func_2para)
LOO_2para = (1/nrow(main_df))*sum((LOO_2para_df$CA_value - LOO_2para_df$CA_predict)^2)
LOO_3para_df = LOO("3para",main_df,func_3para)
LOO_3para = (1/nrow(main_df))*sum((LOO_3para_df$CA_value - LOO_3para_df$CA_predict)^2)
LOO_4para_df = LOO("4para",main_df,func_4para)
LOO_4para = (1/nrow(main_df))*sum((LOO_4para_df$CA_value - LOO_4para_df$CA_predict)^2)
LOO_2paraTE_df = LOO("2paraTE",main_df,func_2paraTE)
LOO_2paraTE = (1/nrow(main_df))*sum((LOO_2paraTE_df$CA_value - LOO_2paraTE_df$CA_predict)^2)
LOO_CV_df = data.frame(CV_value = c(LOO_2para,LOO_3para,LOO_4para,LOO_2paraTE), row.names = c("2para model", "3para model", "4para model", "2paraTE model"))
LOO_CV_df


#Leave On Out Cross Validation for linear
LOOCV_linear <- function(model) {
  h = lm.influence(model)$h
  return(mean((residuals(model)/(1-h))^2))
}
#Above function is not actually needed, we can use CV function from forecast package to get LOOCV, AIC and BIC
LOO_sli_linear1 <- CV(model_sli_linear1)[[1]]
AIC_sli_linear1 <- CV(model_sli_linear1)[[2]]
BIC_sli_linear1 <- CV(model_sli_linear1)[[4]]

LOO_sli_linear2 <- CV(model_sli_linear2)[[1]]
AIC_sli_linear2 <- CV(model_sli_linear2)[[2]]
BIC_sli_linear2 <- CV(model_sli_linear2)[[4]]
SLI_AICs = c(AIC_sli_linear1, AIC_sli_linear2)
SLI_BICs = c(BIC_sli_linear1, BIC_sli_linear2)
SLI_LOOs = c(LOO_sli_linear1, LOO_sli_linear2)
SLI_summary = data.frame(model = c("Linear", "Quadratic"), AIC = SLI_AICs, BIC = SLI_BICs, LOO = SLI_LOOs)
SLI_summary
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
  info_table = main_df %>% group_by(ion, L, Z.b) %>% summarise()
  info_table = suppressWarnings(left_join(data.frame(ions), info_table, by = c("ions" = "ion")))
  output = 0
  for (i in 1: nrow(info_table)){
    output = output + IDER(d*r, L = info_table$L[i], Z.b = info_table$Z.b[i], model = model, parameters = parameters)
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
          dI[i] = r[i]*(sig[i]*6.24/L[i]*(1 - exp(-1024*u[i]/L[i]) + 1024/L[i]*(exp(-1024*u[i]/L[i]))))
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
