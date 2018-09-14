source('/Users/zhaoliyang/Desktop/URAP\ FALL2018/Synergyfunc.R')
#Confidence Intervals (Monte Carlo)
set.seed(19970101)
MM <- 500
#Sample parameters from their distributions with covariances
sig = vcov(IDER_model)
monte_carlo_parameters = rmvnorm(n = MM, mean = c(eta0 = eta0_hat, eta1 = eta1_hat, sig0 = sig0_hat, kap = kap_hat), sigma = sig)
eta0_MC = monte_carlo_parameters[, 1]
eta1_MC = monte_carlo_parameters[, 2]
sig0_MC = monte_carlo_parameters[, 3]
kap_MC = monte_carlo_parameters[, 4]
kap_MC[kap_MC <= 1e-6] = 1e-5 #10 of them were fixed

#Sample parameters from their distributions without covariances
eta0_MC_var = rnorm(MM, mean = eta0_hat, sd = 4.631e-05)
eta1_MC_var = rnorm(MM, mean = eta1_hat, sd = 1.617e-04)
sig0_MC_var = rnorm(MM, mean = sig0_hat, sd = 1.966e+00)
kap_MC_var = rnorm(MM, mean = kap_hat, sd = 2.628e+02)
kap_MC_var[kap_MC_var <= 1e-6] = 1e-5 

# #Monte Carlo functions (without and with covariance)
# CI_function_MIXDER_var = function(d, d_interested, r, interval = 0.95, L, Z.beta) { #Outputs the CI without using the vcov matrix of the parameters
#   MIXDER_curve = list(0)
#   for (i in 1:MM) {
#     MIXDER_curve[[i]] = MIXDER_function(r = r, d = d, L = L, Z.b = Z.b, eta0 = eta0_MC_var[i], eta1 = eta1_MC_var[i], sig0 = sig0_MC_var[i], kap = kap_MC_var[i])
#   }
#   info = vector(length = 0)
#   for (i in 1:MM) {
#     info = c(info, MIXDER_curve[[i]][, 2][2])
#   }
#   info = sort(info)
#   lower_bound = info[(1-interval)/2*MM]
#   upper_bound = info[(interval + (1-interval)/2)*MM]
#   CI = c(lower_bound, upper_bound)
#   return(CI)
# }

################### # EGH 08.20.18 REFACTORING

# Monte Carlo functions (without and with covariance)
refactor_CI_function_MIXDER_var = function(d, d_interested, r, interval = 0.95, L, Z.beta) { #Outputs the CI without using the vcov matrix of the parameters
  MIXDER_curve = list(0)
  for (i in 1:MM) {
    MIXDER_curve[[i]] = MIXDER_function(r = r, d = d, L = L, Z.b = Z.b, eta0 = eta0_MC_var[i], eta1 = eta1_MC_var[i], sig0 = sig0_MC_var[i], kap = kap_MC_var[i])
    cat(paste("  Currently at Monte Carlo step:", toString(i), "of", 
              toString(MM)), sprintf('\r'))
  }
  return(MIXDER_curve)
}

################### END REFACTORING

monte_carlo <- function(ions, r = rep(1/length(ions), length(ions)), d = c(seq(0, 0.009, 0.001), seq(0.01, 0.5, by = 0.01)), eta0 = eta0_hat, eta1 = eta1_hat, sig0 = sig0_hat, kap = kap_hat, n = 500, graph = T, background = 0.00071, cov = T){
  ##SCW: This is a generic function that takes names of the ions in the desired mixture as input and either plots out the ribbon graph or outputs the dataset
  #r is a vector that sums up to 1
  #ions is a vector of strings of the names of the ions in the mixture
  #d is a vector of different dosage
  #If cov = F plots out the result without using the cov matrix
  info_table = modified_df %>% group_by(ion, L, Z.b) %>% summarise() %>% filter(ion %in% ions)
  info_table = suppressWarnings(left_join(data.frame(ions), info_table, by = c("ions" = "ion")))
  L = info_table$L
  Z.b = info_table$Z.b
  out = MIXDER_function(r = r, L = L, Z.b = Z.b, d = d, eta0 = eta0, eta1=eta1, sig0=sig0, kap = kap)
  MIXDER = data.frame(d = d, CA = out[, 2] + background)
  ninty_five_CI_lower = vector(length = length(d))
  ninty_five_CI_upper = vector(length = length(d))
  DERs <- matrix(nrow = n, ncol = length(d))
  if (cov == F){
    for (j in 1:n) {
      DERs[j,] <-  MIXDER_function(r = r, d = d, L = L, Z.b = Z.b, eta0 = eta0_MC_var[j], eta1 = eta1_MC_var[j], sig0 = sig0_MC_var[j], kap = kap_MC_var[j])[, 2]
      cat(paste("  Currently at Monte Carlo step:", toString(j), "of", 
                toString(n)), sprintf('\r'))
    }
  }
  else{
    for (j in 1:n) {
      DERs[j, ] <- MIXDER_function(r = r, L = L, Z.b = Z.b, d = d, eta0 = eta0_MC[j], eta1 = eta1_MC[j], sig0 = sig0_MC[j], kap = kap_MC[j])[, 2]
      cat(paste("  Currently at Monte Carlo step:", toString(j), "of", 
                toString(n)), sprintf('\r'))
    }
  }
  for (i in 1:length(d)) {
    sample_values <- sort(DERs[, i])
    # Returning resulting CI
    ninty_five_CI_lower[i] <- sample_values[ceiling((1 - 0.95) / 2 * n)] + background  #Need to add background prevalence of 0.00071 to every output
    ninty_five_CI_upper[i] <- sample_values[(0.95 + (1 - 0.95) / 2) * n] + background
  }
  
  MIXDER$CI_lower = ninty_five_CI_lower
  MIXDER$CI_upper = ninty_five_CI_upper
  MIXDER$simpleeffect = IDER(d, ions = ions) + length(ions) * background #SCW: I believe we should add background prevalence to every IDER
  names = colnames(MIXDER)
  for (k in ions){
    ider = IDER(d, ions = k)
    MIXDER <- cbind(MIXDER, data.frame(ider))
  }
  colnames(MIXDER) <- c(names, ions)
  if(graph == F){
    return(MIXDER)
  }
  else{
    CA <- MIXDER$CA
    simpleeffect <- MIXDER$simpleeffect
    plot(x = d * 100, y = CA * 100, type = "l", col = "red")
    lines(x = d * 100, y = simpleeffect * 100, col = "black", lty = 2, lwd = 0.5)
    for(i in 6:ncol(MIXDER)){
      lines(x = d * 100, y = 100*as.vector(MIXDER[,i]), col = "green")
    }
    lines(x= d*100 , y = MIXDER$CI_upper * 100, lty = 'dashed', col = 'red')
    lines(x= d*100 , y = MIXDER$CI_lower * 100, lty = 'dashed', col = 'red')
    polygon(c(d*100,rev(d*100)),c(MIXDER$CI_lower * 100, rev(MIXDER$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)
  }
}
