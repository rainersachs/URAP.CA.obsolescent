source('MonteCarlo.R')

MC_results_4para_cov = monte_carlo(ions = c("Fe600", "Si170", "O55", "O350"), para = MC_4para, n = 500) #This outputs a list of MIXDER dataframe and a vector of the indexes at which there were convergence issues. 
MIXDER_4para_cov = MC_results_4para_cov[[1]] #This is the dataframe that will be passed into the graph
non_convergence_4para_cov = MC_results_4para_cov[[2]] #This is a vector of indexes that did not converge. The length of this should be <= 3.

MC_results_4para_var = monte_carlo(ions = c("Fe600", "Si170", "O55", "O350"), para = MC_4para, n = 500, cov = F)
MIXDER_4para_var = MC_results_4para_var[[1]]
non_convergence_4para_var = MC_results_4para_var[[2]]

MC_results_2para_cov = monte_carlo(ions = c("Fe600", "Si170", "O55", "O350"), para = MC_2para, n = 500) 
MIXDER_2para_cov = MC_results_2para_cov[[1]] #This is the dataframe that will be passed into the graph
non_convergence_2para_cov = MC_results_2para_cov[[2]]

#############################################################################################
#Graphing part. 

monte_carlo_graph <- function(MIXDER, model_name = "Model Name"){
  #MIXDER is the output file of the monte carlo function. It could come from 6 different model * var combinations (3*2).
  d = MIXDER$d
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
  title(paste("Monte Carlo Plot for", model_name))
}

monte_carlo_graph(MIXDER_4para_cov, "4-Parameter Model with Covariances")
monte_carlo_graph(MIXDER_4para_var, "4-Parameter Model without Covariances")
monte_carlo_graph(MIXDER_2para_cov, "2-Parameter Model with Covariances")

