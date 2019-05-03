source('Synergy.R')

#In this section we are estimating the confidence intervals of the individual and mixture DERs using Monte Carlo

set.seed(19970101)
N <- 500

#Monte-Carlo Parameters from Multivariate Normal Distributions
###Note: We are conditioning on all parameters being positive (by dropping them negative ones and redraw the sample). 
sample_parameters = function(model, n = N, cov = T){ #This function samples n sets of parameters from their joint distribution defined by the calibration process.
  #Note: The calibrated parameter dataframes from Synergy.R file is required
  if (model == "4para"){
    calibrated = parameters_4para
    covariance = sig_4para
  }
  else if (model == "3para"){
    calibrated = parameters_3para
    covariance = sig_3para
  }
  else if (model == "2para"){
    calibrated = parameters_2para
    covariance = sig_2para
  }
  else if (model == "2paraTE"){
    calibrated = parameters_2paraTE
    covariance = sig_2paraTE
  }
  if (!cov){ #If we are not using the covariance matrix, then the parameters will be independent
    n_para = nrow(calibrated)
    zeros = matrix(0, n_para, n_para)
    variances = as.numeric(diag(covariance))
    diag(zeros) = variances
    covariance = zeros
  }
  mc_parameters = rmvnorm(n = n, mean = calibrated$value, sigma = covariance)
  negatives = which(apply(mc_parameters, 1, function(row) any(row < 1e-6)))
  n_negatives = length(negatives)
  print(paste("There are ", n_negatives, " rows containing negative values out of ", n, " draws", sep = ""))
  if(n_negatives > 0){
    mc_parameters = mc_parameters[-negatives,]
  }
  while (nrow(mc_parameters) < n){ #Running a while loop to get enough draws from the multivariate Normal
    new_row = rmvnorm(n = 1, mean = calibrated$value, sigma = covariance)
    if (!any(new_row < 1e-6)){
      mc_parameters = rbind(mc_parameters, new_row)
    }
  }
  colnames(mc_parameters) = calibrated$parameter
  return(mc_parameters)
}

#A function that makes monte_carlo parameters into N dataframes
make_datapara <- function(l, model, n = N){
  #l is a vector of lists, each list containing N sampled values
  name = model
  df_list = list()
  if (model == "4para"){
    paras = c("eta0", "eta1", "sig0", "kap0")
    for (i in 1:n){
      value = c((l$eta0)[i], (l$eta1)[i], (l$sig0)[i], (l$kap0)[i])
      df = data.frame(value, model = name, parameter = paras)
      df_list[[i]] = df
    }
  }
  else if (model == "3para"){
    paras = c("eta0", "eta1", "sig0")
    for (i in 1:n){
      value = c((l$eta0)[i], (l$eta1)[i], (l$sig0)[i])
      df = data.frame(value, model = name, parameter = paras)
      df_list[[i]] = df
    }
  }
  else if (model == "2para"){
    paras = c("eta0", "sig0")
    for (i in 1:n){
      value = c((l$eta0)[i], (l$sig0)[i])
      df = data.frame(value, model = name, parameter = paras)
      df_list[[i]] = df
    }
  }
  else if (model == "2paraTE"){
    paras = c("sig0", "kap0")
    for (i in 1:n){
      value = c((l$sig0)[i], (l$kap0)[i])
      df = data.frame(value, model = name, parameter = paras)
      df_list[[i]] = df
    }
  }
  return(df_list)
}

#The main function that calculates the errors for each of the N = 500 sets of parameters.
monte_carlo <- function(ions, r = rep(1/length(ions), length(ions)), para = MC_4para, d = c(seq(0, 0.009, 0.001), seq(0.01, 0.5, by = 0.01)), n = N, background = 0, cov = T, model = NULL, IDER = F){
  ##SCW: This is a generic function that takes names of the ions in the desired mixture as input and either plots out the ribbon graph or outputs the dataset
  #r is a vector that sums up to 1
  #ions is a vector of strings of the names of the ions in the mixture
  #d is a vector of different dosage
  #para is a list (MC_2para, MC_3para, or MC_4para), and it will tell the function which model to use
  #If cov = F plots out the result without using the cov matrix
  #the maximum number for n is N (in this case N = 500)
  #If IDER = T then this outputs ribbons for IDER instead of MIXDER
  if (is.null(model)){
    model = para$name
  }  
  if (IDER){ #This section is for IDER ribbon
    if (length(ions) > 1) {stop("Only one ion allowed for IDER")}
    out = IDER(d, ions = ions, model = model, parameters = Data_parameter)
    IDER_table = data.frame(d = d, IDER = out + background)
    if (cov == F){
      MC_parameter = para[[2]]
    }
    else{
      MC_parameter = para[[1]]
    }
    bad_j = c()
    DERs <- matrix(nrow = n, ncol = length(d))
    for (j in 1:n) {
      der = IDER(d, ions = ions, model = model, parameters = MC_parameter[[j]])
      if (class(der) == "try-error"){ #Does not hault the 
        bad_j = c(bad_j, j)
        DERs[j,] <- out
      }
      else{
        DERs[j,] = der
      }
      cat(paste("  Currently at Monte Carlo step:", toString(j), "of", 
                toString(n)), sprintf('\r'))
    }
    if (length(bad_j) > 0){
      print(paste("Non-Convergence at i = ", bad_j))
    }
    ninty_five_CI_lower = numeric(length(d))
    ninty_five_CI_upper = numeric(length(d))
    for (i in 1:length(d)) {
      sample_values <- sort(DERs[, i])
      # Returning resulting CI
      ninty_five_CI_lower[i] <- as.numeric(quantile(sample_values, 0.025)) + background  #Background default is 0 (for DER)
      ninty_five_CI_upper[i] <- as.numeric(quantile(sample_values, 0.975)) + background
    }
    IDER_table$CI_lower = ninty_five_CI_lower
    IDER_table$CI_upper = ninty_five_CI_upper
    return(list(IDER_table, bad_j))
  }
  
  else{ #This section is for MIXDER ribbon
    info_table = main_df %>% group_by(ion, L, Z.b) %>% summarise() %>% filter(ion %in% ions)
    info_table = suppressWarnings(left_join(data.frame(ions), info_table, by = c("ions" = "ion")))
    L = info_table$L
    Z.b = info_table$Z.b
    out = MIXDER_function(r = r, L = L, Z.b = Z.b, d = d, model = model, parameters = Data_parameter)
    MIXDER = data.frame(d = d, CA = out[, 2] + background)
    ninty_five_CI_lower = vector(length = length(d))
    ninty_five_CI_upper = vector(length = length(d))
    DERs <- matrix(nrow = n, ncol = length(d))
    if (cov == F){
      MC_parameter = para[[2]]
    }
    else{
      MC_parameter = para[[1]]
    }
    bad_j = c()
    for (j in 1:n) {
      der = try(MIXDER_function(r = r, d = d, L = L, Z.b = Z.b, model = model, parameters = MC_parameter[[j]])[, 2])
      if (class(der) == "try-error"){ #Does not hault the 
        bad_j = c(bad_j, j)
        DERs[j,] <- out[,2]
      }
      else{
        DERs[j,] = der
      }
      cat(paste("  Currently at Monte Carlo step:", toString(j), "of", 
                toString(n)), sprintf('\r'))
    }
    if (length(bad_j) > 0){
      print(paste("Non-Convergence at i = ", bad_j))
    }
    for (i in 1:length(d)) {
      sample_values <- sort(DERs[, i])
      # Returning resulting CI
      ninty_five_CI_lower[i] <- as.numeric(quantile(sample_values, 0.025)) + background  #Background default is 0 (for DER)
      ninty_five_CI_upper[i] <- as.numeric(quantile(sample_values, 0.975)) + background
    }
    MIXDER$CI_lower = ninty_five_CI_lower
    MIXDER$CI_upper = ninty_five_CI_upper
    MIXDER$simpleeffect = IDER(d, ions = ions, model = model, parameters = Data_parameter) + background #SCW: I believe we should add background prevalence to every IDER
    names = colnames(MIXDER)
    for (k in ions){
      ider = IDER(d, ions = k, model = model, parameters = Data_parameter)
      MIXDER <- cbind(MIXDER, data.frame(ider))
    }
    colnames(MIXDER) <- c(names, ions)
    return(list(MIXDER, bad_j))
  }
}


#################################################### Begin Sampling


#4-Parameter Model Monte Carlo Sampling
MC_4para_cov_samples = sample_parameters(model = "4para", n = N, cov = T)
para_4para_cov =  c(eta0 = list(MC_4para_cov_samples[,"eta0"]), 
                    eta1 = list(MC_4para_cov_samples[,"eta1"]),
                    sig0 = list(MC_4para_cov_samples[,"sig0"]), 
                    kap0 = list(MC_4para_cov_samples[,"kap0"]))
MC_4para_var_samples = sample_parameters(model = "4para", n = N, cov = F)
para_4para_var =  c(eta0 = list(MC_4para_var_samples[,"eta0"]), 
                    eta1 = list(MC_4para_var_samples[,"eta1"]),
                    sig0 = list(MC_4para_var_samples[,"sig0"]), 
                    kap0 = list(MC_4para_var_samples[,"kap0"]))
MC_4para_cov = make_datapara(para_4para_cov, model = "4para", n = N)
MC_4para_var = make_datapara(para_4para_var, model = "4para", n = N)

#3-Parameter Model Monte Carlo Sampling

MC_3para_cov_samples = sample_parameters(model = "3para", n = N, cov = T)
para_3para_cov =  c(eta0 = list(MC_3para_cov_samples[,"eta0"]), 
                    eta1 = list(MC_3para_cov_samples[,"eta1"]),
                    sig0 = list(MC_3para_cov_samples[,"sig0"]))

MC_3para_var_samples = sample_parameters(model = "3para", n = N, cov = F)
para_3para_var =  c(eta0 = list(MC_3para_var_samples[,"eta0"]), 
                    eta1 = list(MC_3para_var_samples[,"eta1"]),
                    sig0 = list(MC_3para_var_samples[,"sig0"]))

MC_3para_cov = make_datapara(para_3para_cov, model = "3para",n = N)
MC_3para_var = make_datapara(para_3para_var, model = "3para",n = N)

#2-Parameter Model Monte Carlo Sampling

MC_2para_cov_samples = sample_parameters(model = "2para", n = N, cov = T)
para_2para_cov =  c(eta0 = list(MC_2para_cov_samples[,"eta0"]), 
                    sig0 = list(MC_2para_cov_samples[,"sig0"]))

MC_2para_var_samples = sample_parameters(model = "2para", n = N, cov = F)
para_2para_var =  c(eta0 = list(MC_2para_var_samples[,"eta0"]), 
                    sig0 = list(MC_2para_var_samples[,"sig0"]))

MC_2para_cov = make_datapara(para_2para_cov, model = "2para",n = N)
MC_2para_var = make_datapara(para_2para_var, model = "2para",n = N)

#2-Parameter Model Monte Carlo Sampling

MC_2paraTE_cov_samples = sample_parameters(model = "2paraTE", n = N, cov = T)
para_2paraTE_cov =  c(sig0 = list(MC_2paraTE_cov_samples[,"sig0"]), 
                    kap0 = list(MC_2paraTE_cov_samples[,"kap0"]))

MC_2paraTE_var_samples = sample_parameters(model = "2paraTE", n = N, cov = F)
para_2paraTE_var =  c(sig0 = list(MC_2paraTE_var_samples[,"sig0"]), 
                      kap0 = list(MC_2paraTE_var_samples[,"kap0"]))

MC_2paraTE_cov = make_datapara(para_2paraTE_cov, model = "2paraTE",n = N)
MC_2paraTE_var = make_datapara(para_2paraTE_var, model = "2paraTE",n = N)



#Combining them into lists for monte_carlo function
MC_2para = list(MC_2para_cov, MC_2para_var, name = "2para")
MC_3para = list(MC_3para_cov, MC_3para_var, name = "3para")
MC_4para = list(MC_4para_cov, MC_4para_var, name = "4para")
MC_2paraTE =  list(MC_2paraTE_cov, MC_2paraTE_var, name = "2paraTE")
###############################################



