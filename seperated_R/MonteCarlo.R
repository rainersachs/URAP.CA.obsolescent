source('Synergy.R')

#Confidence Intervals (Monte Carlo)
set.seed(19970101)
MM <- 500

#A function that makes monte_carlo parameters into MM dataframes
make_datapara <- function(l, model, n = MM){
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
  return(df_list)
}

####################################################
#Monte-Carlo Parameters
#In this section, we are conditioning on kappa > 0 (by dropping them negative ones and redraw the sample). 
#Also, we are making sig0 > 0 (by replacing the negative values with a small positive number 1e-6)

#Sample parameters from their distributions with covariances for 4para
monte_carlo_parameters = rmvnorm(n = MM, mean = c(eta0 = filter(parameters_4para, parameter == "eta0")$value, eta1 = filter(parameters_4para, parameter == "eta1")$value, sig0 = filter(parameters_4para, parameter == "sig0")$value, kap = filter(parameters_4para, parameter == "kap0")$value), sigma = sig_4para)
monte_carlo_parameters = monte_carlo_parameters[monte_carlo_parameters[, "kap"] > 1e-6, ] #Dropping negative kappas
monte_carlo_parameters = monte_carlo_parameters[monte_carlo_parameters[, "eta0"] > 1e-6, ]
monte_carlo_parameters = monte_carlo_parameters[monte_carlo_parameters[, "eta1"] > 1e-6, ]
monte_carlo_parameters = monte_carlo_parameters[monte_carlo_parameters[, "sig0"] > 1e-6, ]

negative_4para_cov = MM - nrow(monte_carlo_parameters) #The number of kap < 0


############THIS
while (nrow(monte_carlo_parameters) < MM){ #Running a while loop to get enough draws from the multivariate Normal
  new_row = rmvnorm(n = 1,mean = c(eta0 = filter(parameters_4para, parameter == "eta0")$value, eta1 = filter(parameters_4para, parameter == "eta1")$value, sig0 = filter(parameters_4para, parameter == "sig0")$value, kap = filter(parameters_4para, parameter == "kap0")$value), sigma = sig_4para)
  if (as.numeric(new_row[, "kap"]) > 1e-6 & as.numeric(new_row[, "eta0"]) > 1e-6 & as.numeric(new_row[, "eta1"]) > 1e-6 & as.numeric(new_row[, "sig0"]) > 1e-6){
    monte_carlo_parameters = rbind(monte_carlo_parameters, new_row)
  }
}

eta0_MC = monte_carlo_parameters[, "eta0"]
eta1_MC = monte_carlo_parameters[, "eta1"]
sig0_MC = monte_carlo_parameters[, "sig0"]
kap_MC = monte_carlo_parameters[, "kap"]


#Sample parameters from their distributions without covariances
eta0_MC_var = rnorm(MM, mean = filter(parameters_4para, parameter == "eta0")$value, sd = sqrt(sig_4para[1,1]))
eta1_MC_var = rnorm(MM, mean = filter(parameters_4para, parameter == "eta1")$value, sd = sqrt(sig_4para[2,2]))
sig0_MC_var = rnorm(MM, mean = filter(parameters_4para, parameter == "sig0")$value, sd = sqrt(sig_4para[3,3]))
kap_MC_var = rnorm(MM, mean = filter(parameters_4para, parameter == "kap0")$value, sd = sqrt(sig_4para[4,4]))

negative_kap_4para_var = MM - length(kap_MC_var) #The number of kap < 0

while (length(kap_MC_var) < MM){  #Running a while loop to get enough draws from the normal 
  new_value = rnorm(1, mean =filter(parameters_4para, parameter == "kap0")$value, sqrt(sig_4para[4,4]))
  if (new_value > 1e-6){
    kap_MC_var = c(kap_MC_var, new_value)
  }
}

para_4para_cov =  c(eta0 = list(eta0_MC), eta1 = list(eta1_MC),sig0 = list(sig0_MC), kap0 = list(kap_MC))
para_4para_nocov = c(eta0 = list(eta0_MC_var), eta1 = list(eta1_MC_var),sig0 = list(sig0_MC_var),kap0 =list(kap_MC_var))

#Sample parameters from their distributions with covariances for Three parameter model.
monte_carlo_parameters = rmvnorm(n = MM, mean = c(eta0 = filter(parameters_3para, parameter == "eta0")$value, eta1 = filter(parameters_3para, parameter == "eta1")$value, sig0 = filter(parameters_3para, parameter == "sig0")$value), sigma = sig_3para)

negative_sig0_3para_cov = sum(monte_carlo_parameters[,"sig0"] <= 1e-6) #The number of sig0 < 0

#monte_carlo_parameters[monte_carlo_parameters[,"sig0"] <= 1e-6, ]["sig0"] = 1e-5
monte_carlo_parameters = monte_carlo_parameters[monte_carlo_parameters[, "eta0"] > 1e-6, ] # Drops
monte_carlo_parameters = monte_carlo_parameters[monte_carlo_parameters[, "eta1"] > 1e-6, ]
monte_carlo_parameters = monte_carlo_parameters[monte_carlo_parameters[, "sig0"] > 1e-6, ]

while (nrow(monte_carlo_parameters) < MM){ #Running a while loop to get enough draws from the multivariate Normal
  new_row = rmvnorm(n = 1, mean = c(eta0 = filter(parameters_3para, parameter == "eta0")$value, eta1 = filter(parameters_3para, parameter == "eta1")$value, sig0 = filter(parameters_3para, parameter == "sig0")$value), sigma = sig_3para)
  if (as.numeric(new_row[, "eta0"]) > 1e-6 & as.numeric(new_row[, "eta1"]) > 1e-6 & as.numeric(new_row[, "sig0"]) > 1e-6) {
    monte_carlo_parameters = rbind(monte_carlo_parameters, new_row)
  }
}



eta0_MC = monte_carlo_parameters[,"eta0"]
eta1_MC = monte_carlo_parameters[,"eta1"]
sig0_MC = monte_carlo_parameters[,"sig0"]

#Sample parameters from their distributions without covariances for Three parameter model.
eta0_MC_var = rnorm(MM, mean = filter(parameters_3para, parameter == "eta0")$value, sd = sqrt(sig_3para[1,1]))
eta1_MC_var = rnorm(MM, mean = filter(parameters_3para, parameter == "eta1")$value, sd = sqrt(sig_3para[2,2]))
sig0_MC_var = rnorm(MM, mean = filter(parameters_3para, parameter == "sig0")$value, sd = sqrt(sig_3para[3,3]))


negative_sig0_3para_var = sum(sig0_MC_var <= 1e-6) #The number of sig0 < 0
sig0_MC_var[sig0_MC_var <= 1e-6] = 1e-5

para_3para_cov =  c(eta0 = list(eta0_MC),eta1 = list(eta1_MC), sig0 = list(sig0_MC))
para_3para_nocov = c(eta0 = list(eta0_MC_var),eta1 = list(eta1_MC_var), sig0 = list(sig0_MC_var))

#Sample parameters from their distributions with covariances for 2para
monte_carlo_parameters = rmvnorm(n = MM, mean = c(eta0 = filter(parameters_2para, parameter == "eta0")$value ,sig0 = filter(parameters_2para, parameter == "sig0")$value), sigma = sig_2para)
negative_sig0_2para_cov = sum(monte_carlo_parameters[, "sig0"] <= 1e-6) #The number of sig0 < 0

#monte_carlo_parameters[monte_carlo_parameters[,"sig0"] <= 1e-6, ]["sig0"] = 1e-5
monte_carlo_parameters = monte_carlo_parameters[monte_carlo_parameters[, "eta0"] > 1e-6, ] # Drops
#monte_carlo_parameters = monte_carlo_parameters[monte_carlo_parameters[, "eta1"] > 1e-6, ]
monte_carlo_parameters = monte_carlo_parameters[monte_carlo_parameters[, "sig0"] > 1e-6, ]

while (nrow(monte_carlo_parameters) < MM){ #Running a while loop to get enough draws from the multivariate Normal
  new_row = rmvnorm(n = 1, mean = c(eta0 = filter(parameters_2para, parameter == "eta0")$value ,sig0 = filter(parameters_2para, parameter == "sig0")$value), sigma = sig_2para)
  if (as.numeric(new_row[, "eta0"]) > 1e-6 & as.numeric(new_row[, "sig0"]) > 1e-6) {
    monte_carlo_parameters = rbind(monte_carlo_parameters, new_row)
  }
}

eta0_MC = monte_carlo_parameters[, "eta0"]
sig0_MC = monte_carlo_parameters[, "sig0"]


#Sample parameters from their distributions without covariances for 2para

eta0_MC_var = rnorm(MM, mean = filter(parameters_2para, parameter == "eta0")$value, sd = sqrt(sig_2para[1,1]))
sig0_MC_var = rnorm(MM, mean = filter(parameters_2para, parameter == "sig0")$value, sd = sqrt(sig_2para[2,2]))

negative_sig0_2para_var = sum(sig0_MC_var <= 1e-6) #The number of sig0 < 0
sig0_MC_var[sig0_MC_var <= 1e-6] = 1e-5

para_2para_cov =  c(eta0 = list(eta0_MC),sig0 = list(sig0_MC))
para_2para_nocov= c(eta0 = list(eta0_MC_var),sig0 = list(sig0_MC_var))

# #Sample parameters from their distributions with covariances for 2para TE only
# monte_carlo_parameters = rmvnorm(n = MM, mean = c(eta0 = Data_parameter[9,1], sig0 = Data_parameter[10,1], sigma = sig_2paraTE))
# sig0_MC = monte_carlo_parameters[, 1]
# kap0_MC = monte_carlo_parameters[, 2]
# #Sample parameters from their distributions without covariances for 2para TE only
# sig0_MC_var = rnorm(MM, mean = Data_parameter[9,1], sd = sqrt(sig_2paraTE[1,1]))
# kap0_MC_var = rnorm(MM, mean = Data_parameter[10,1], sd = sqrt(sig_2paraTE[2,2]))
# para_2paraTE_nocov =  c(list(sig0_MC),list(kap0_MC))
# para_2paraTE_cov= c(list(sig0_MC_var),list(kap0_MC_var))

#Monte Carlo
MC_2para_cov = make_datapara(para_2para_cov, model = "2para",n = MM)
MC_2para_var = make_datapara(para_2para_nocov, model = "2para",n = MM)
MC_3para_cov = make_datapara(para_3para_cov, model = "3para",n = MM)
MC_3para_var = make_datapara(para_3para_nocov, model = "3para",n = MM)
MC_4para_cov = make_datapara(para_4para_cov, model = "4para", n = MM)
MC_4para_var = make_datapara(para_4para_nocov, model = "4para", n = MM)

MC_2para = list(MC_2para_cov, MC_2para_var, name = "2para")
MC_3para = list(MC_3para_cov, MC_3para_var, name = "3para")
MC_4para = list(MC_4para_cov, MC_4para_var, name = "4para")
###############################################

monte_carlo <- function(ions, r = rep(1/length(ions), length(ions)), para = MC_4para, d = c(seq(0, 0.009, 0.001), seq(0.01, 0.5, by = 0.01)), n = MM, background = 0, cov = T, model = NULL, IDER = F){
  ##SCW: This is a generic function that takes names of the ions in the desired mixture as input and either plots out the ribbon graph or outputs the dataset
  #r is a vector that sums up to 1
  #ions is a vector of strings of the names of the ions in the mixture
  #d is a vector of different dosage
  #para is a list (MC_2para, MC_3para, or MC_4para), and it will tell the function which model to use
  #If cov = F plots out the result without using the cov matrix
  #the maximum number for n is MM (in this case MM = 500)
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
    info_table = modified_df %>% group_by(ion, L, Z.b) %>% summarise() %>% filter(ion %in% ions)
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


