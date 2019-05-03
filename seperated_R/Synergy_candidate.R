source('Synergy.R')

#Our IDERs (Individual Dose Effect Relations). Applicable to the 1-ion components of a mixed simulated GCR beam 
#Modifying NTE1 and NTE2 by insisting they be twice continuously differentiable and monotonic increasing. Double check NTE1, NTE2, Our model
library("forecast")
swift_light_df 
main_df = main_df %>% filter(Z > 3)
main_df = main_df %>% filter(d > 0)
##########################################################################################
#Model Equation:
##Example: 4-para model
func_4para = function(d, L, Z.b, eta0, eta1, sig0, kap0) {
  P = (1-exp(-Z.b/kap0))^2
  sig = sig0*P + 0.041/6.24*L*(1-P) # 0.041 +- 0.0051 comes from 16Cacao
  eta = eta0*L*exp(-eta1*L)
  return(BG_CA + sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^5*d)))  #0.00071 + sig*6.24*d/L*(1-exp(-1024*d/L)) + eta*(1-exp(-10^3*d))#don't use
} 

##New candidate model 
func_NAME= function(...)

##########################################################################################
#Model Calibration:

##Example: Calibration for 4-parameter model
model_4para = nlsLM(CA ~ func_4para(d, L, Z.b, eta0, eta1, sig0, kap0), 
                    data = main_df, 
                    start = list(eta0 = 0.001, eta1 = 0.01, sig0 = 5, kap0 = 500), 
                    weights = (1/(main_df$error)^2))
ccoef <- coef(model_4para)
parameters_4para = data.frame(value = as.numeric(ccoef), model = "4para", parameter = names(ccoef))
sig_4para = vcov(model_4para)

##Calibration for candidate
model_NAME = nlsLM(...) #need function name, parameters, etc. See example above
ccoef <- coef(model_NAME)
parameters_NAME = data.frame(value = as.numeric(ccoef), model = "NAME", parameter = names(ccoef))
sig_NAME = vcov(model_NAME) #storing the covariance matrix



Data_parameter = rbind(#This was already calculated in Synergy.R already. Just add the nlsLM fit for the new candidate here
                       parameters_2para,
                       #parameters_NewCandidate,
                       parameters_3para, 
                       parameters_4para, 
                       parameters_2paraTE, 
                       parameters_sli_lin1, 
                       parameters_sli_lin2, 
                       parameters_sli_exp2, 
                       parameters_sli_exp3)

#####################################################################################################
##AIC and BIC
L_function_candiate = function(func, eta0=0, eta1=0, sig0=0, kap=0, model){ #Note that if the new candidates have more than 4 parameters, this needs to be adjusted
  #A function trhat returns the residual squared
  if (model == "Candidate1"){
    a = vector(length = 0)
    for (i in 1:length(main_df[, 1])) {
      a = c(a, main_df$CA[i] - func(d = main_df$d[i], 
                                    L = main_df$L[i], 
                                    Z.b = main_df$Z.b[i], 
                                    eta0 = eta0, eta1 = eta1, sig0 = sig0, kap = kap))}
    return(a^2)
  }
  if (model == "Candidate2"){
    a = vector(length = 0)
    for (i in 1:length(main_df[, 1])) {
      a = c(a, main_df$CA[i] - func(d = main_df$d[i], 
                                    L = main_df$L[i], 
                                    Z.b = main_df$Z.b[i], 
                                    eta0 = eta0, eta1 = eta1, sig0 = sig0, kap = kap))}
    return(a^2)
  }
  else{stop("Use L_function instead. This function is only used for new candidates.")}
}

L_NEWCANDIDATE = L_function_candiate(func_NEWCANDIDATE, eta0 = filter(parameters_NEWCANDIDATE, parameter == "eta0")$value,
                                     sig0 = filter(parameters_NEWCANDIDATE, parameter == "sig0")$value, 
                                     model ="candidate1")
WRSS_NEWCANDIDATE = sum((1/main_df$error^2)*L_NEWCANDIDATE)


AIC_4para = AIC_function(k=4, RSS = WRSS_4para) #Example
BIC_4para = BIC_function(k=4, RSS = WRSS_4para)


AIC_NEWCANDIDATE = AIC_function(k=3, RSS = WRSS_NEWCANDIDATE) #k is the number of parameters
BIC_NEWCANDIDATE = BIC_function(k=3, RSS = WRSS_NEWCANDIDATE)

information_criteria = data.frame(AIC = c(#AIC_NEWCANDIDATE,
                                          AIC_2para,AIC_3para,AIC_4para,AIC_2paraTE), 
                                  BIC = c(#BIC_NEWCANDIDATE,
                                          BIC_2para,BIC_3para,BIC_4para,BIC_2paraTE), 
                                    row.names = c(#"candidate1",
                                                  "2para model", "3para model", "4para model", "2paraTE model"
                                                  ))
information_criteria 
#This dataframe for only the original models is available as well, stored in object name: information_criteria_df

####################################################################################################################
##Leave-one-out cross validation
LOO_candidate = function(model,main_df,func){
  df <- data.frame(matrix(nrow=nrow(main_df), ncol=2)) 
  x <- c("CA_value", "CA_predict")
  colnames(value) <- x
  if (model == "2paraTE"){ #Example
    for (i in 1:nrow(main_df)){
      train = main_df[-i,]
      test = main_df[i,]
      model_fit = nlsLM(CA ~ func_2paraTE(d,L,Z.b,kap0,sig0), #Put in the corresponding candidate model function and parameters
                            data = train, start = list(sig0 = 5, kap0 = 500), #Corresponding Starting Values
                            weights = (1/(train$error)^2))
      df$CA_value[i] = test$CA
      df$CA_predict[i] = predict(model_fit, test)}
    
    return((1/nrow(main_df))*sum((df$CA_value - df$CA_predict)^2))
  }
  if(model == "candidate1"){
    
  }
  else{stop("Use L_function instead. This function is only used for new candidates.")}
}

LOO_NEWCANDIDATE = LOO_candidate("candidate1",main_df,func_NEWCANDIDATE)
LOO_CV_df = data.frame(CV_value = c(#LOO_NEWCANDIDATE,
                                    LOO_2para,LOO_3para,LOO_4para,LOO_2paraTE), 
                       row.names = c(#"candidate1",
                                      "2para model", "3para model", "4para model", "2paraTE model"))
LOO_CV_df
