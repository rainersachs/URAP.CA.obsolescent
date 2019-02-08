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

MC_results_3para_cov = monte_carlo(ions = c("Fe600", "Si170", "O55", "O350"), para = MC_3para, n = 500) #This outputs a list of MIXDER dataframe and a vector of the indexes at which there were convergence issues. 
MIXDER_3para_cov = MC_results_3para_cov[[1]] #This is the dataframe that will be passed into the graph
non_convergence_3para_cov = MC_results_3para_cov[[2]]

MC_results_3para_var = monte_carlo(ions = c("Fe600", "Si170", "O55", "O350"), para = MC_3para, n = 500, cov = F) #This outputs a list of MIXDER dataframe and a vector of the indexes at which there were convergence issues. 
MIXDER_3para_var = MC_results_3para_var[[1]] #This is the dataframe that will be passed into the graph
non_convergence_3para_var = MC_results_3para_var[[2]]
#############################################################################################
#Functions 

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

IDER_graph <- function(ions = "O350", model = "4para", data = modified_df, point = F, d = c(seq(0, 0.009, 0.001), seq(0.01, 0.5, by = 0.01))){
  CA = IDER(d = d, ions = ions, model = model)
  data_filtered = data[(data$ion %in% ions) & (data$d <= 0.5),]
  data_d = 100*data_filtered$d
  data_CA = data_filtered$CA
  n = length(ions)
  names = ions[1]
  if (n > 1){
    for (i in 2:n){
      names = paste(names, ",", ions[i])
    }
  }
  print(max(data_CA))
  plot(x = d * 100, y = CA, type = "l", col = "green", ylim = c(0, 1.5*max(data_CA)))
  title(paste(model, "IDER Plot for", names))
  if(point){#Option to add the scatter plot from the true data onto the IDER plot
    points(data_d, data_CA)
  }
}


##################################################################
#Plots

monte_carlo_graph(MIXDER_4para_cov, "4-Parameter Model with Covariances")
monte_carlo_graph(MIXDER_4para_var, "4-Parameter Model without Covariances")
monte_carlo_graph(MIXDER_3para_cov, "3-Parameter Model with Covariances")
monte_carlo_graph(MIXDER_3para_var, "3-Parameter Model without Covariances")

monte_carlo_graph(MIXDER_2para_cov, "2-Parameter Model with Covariances")

IDER_graph(ions = "O55", model = "4para", point = T)
IDER_graph(ions = "O350", model = "4para", point = T)
IDER_graph(ions = "O55", model = "3para", point = T)
IDER_graph(ions = "O350", model = "3para", point = T)

IDER_graph(ions = "O55", model = "2para", point = T)

IDER_graph(ions = c("O350", "O55"), model = "3para", point = T) #If the ions argument has more than 1 ion, the output line is an average (r = 1/n)
IDER_graph(model = "2para", point = F)