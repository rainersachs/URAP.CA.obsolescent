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

IDER_graph <- function(data = modified_df, ions = "O350", model = "4para", point = T, d = c(seq(0, 0.009, 0.001), seq(0.01, 0.5, by = 0.01))){
  if(nrow(data) == nrow(modified_df)){
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
  else{
    ions = unique(data$ion)
    CA = IDER(d = d, ions = ions, model = model)
    data_filtered = data
    data_d = 100*data_filtered$d
    data_CA = data_filtered$CA
    n = length(ions)
    names = ions[1]
    if (n > 1){
      for (i in 2:n){
        names = paste(names, ",", ions[i])
      }
    }
    plot(x = d * 100, y = CA, type = "l", col = "green", ylim = c(0, 1.5*max(data_CA)))
    title(paste(model, "IDER Plot for", names))
    if(point){#Option to add the scatter plot from the true data onto the IDER plot
      points(data_d, data_CA)
    }
  }
}

##################################################################
#Selecting Specific Data by Using Function subset().
#Example: Subset: d >= 0.05, CA > 0.01, Z = 22 and Comment section includes "2018"
subset(modified_df, d >= 0.05 & CA > 0.01 & Z == 22 & grepl("2018", Comment))
#This returns a sub dataset with all the rows that satisfy the restrictions put in, which you can put in to the function "IDER_graph" as the only input. See example below.
IDER_graph(subset(modified_df, d >= 0.05 & CA > 0.01 & Z == 22 & grepl("2018", Comment)))
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


###################################################################
#Plotting Ideas

library(ggplot2)
library(tidyr)
theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #This centers plot titles for ggplots
DER_df = modified_df
bg = mean(modified_df[modified_df$d ==0,]$CA)
DER_df$CA = DER_df$CA - bg
DER_df = filter(DER_df, d > 0)

#Monte Carlo with Actual Data
example = subset(DER_df, ion %in% c("Fe600", "Si170", "O55", "O350") & d <= 0.5)
MIXDER = MIXDER_4para_cov
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
points(example$d * 100, example$CA * 100)
polygon(c(d*100,rev(d*100)),c(MIXDER$CI_lower * 100, rev(MIXDER$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)

#IDER Plot with ribbons (Si170 as example since it has the most data points)
Si170 = subset(DER_df, ion == "Si170" & d <= 0.5)
MC_results_4para_cov_Si170 = monte_carlo(ions = c("Si170"), para = MC_4para, n = 500, IDER = T) 
IDER_4para_cov_Si170 = MC_results_4para_cov_Si170[[1]] #This is the dataframe that will be passed into the graph
non_convergence_4para_cov_Si170 = MC_results_4para_cov_Si170[[2]] #This is a vector of indexes that did not converge. The length of this should be <= 3.
MC_results_4para_var_Si170 = monte_carlo(ions = c("Si170"), para = MC_4para, n = 500, IDER = T, cov = F) 
IDER_4para_var_Si170 = MC_results_4para_var_Si170[[1]] #This is the dataframe that will be passed into the graph
non_convergence_4para_var_Si170 = MC_results_4para_var_Si170[[2]] #This is a vector of indexes that did not converge. The length of this should be <= 3.

IDER_Si170 = IDER_4para_cov_Si170
ggplot(data = IDER_Si170, aes(x = d * 100)) + #This example consists of one ion only: Si170
  
  geom_ribbon(data = IDER_4para_var_Si170, aes(ymin = 100 * CI_lower, ymax = 100 * CI_upper, fill = "Variance"), alpha = 0.3) + #MonteCarlo Ribbons without cov (blue ribbons)
  geom_ribbon(aes(ymin = 100 * CI_lower, ymax = 100 * CI_upper, fill = "Covariance"), alpha = 0.5) + #MonteCarlo Ribbons with cov (red ribbons)
  geom_line(aes(y = 100 * IDER), size = 1) +  #IDER (colored solid lines)

  geom_point(data = Si170, aes(x= 100*d, y = 100*CA), size = 2) + #Actual Data for these for types of ions (colored points)
  geom_errorbar(data = Si170, aes(x = 100*d, ymin = 100*(CA - error), ymax = 100*(CA + error)), width = 0.5) + #Errorbars for each point from original data
  
  scale_fill_manual("Ions", breaks = c("Covariance", "Variance"), values=c("blue", "red")) + #legend for cov vs var
  
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),  #Some random tweeks of the ggplot theme. Need someone more artistic to do this part :)
        panel.grid.major.y = element_line(colour = "grey80"),
        panel.grid.minor.y = element_line(linetype = "dashed", colour = "grey80"),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_line(),
        plot.title = element_text(size = 20, face = "bold"), 
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")
  ) +
  
  labs(x = "100*Dosage", y = "100*Chromosomal Aberrations", title = "Individual Dose Effect Relationship", subtitle = "Si170") #axis labels and plot title


#MIXDER Plot, using the same 4 ions.
names(MIXDER)[2] <- "mixDER"
ggplot(data = MIXDER, aes(x = d * 100)) + #This example consists a mixture of 4 ions: "Fe600", "Si170", "O55", "O350". 
  geom_line(aes(y = 100 * Fe600, color = "Fe600"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100 * O350, color = "O350"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100 * O55, color = "O55"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100 * Si170, color = "Si170"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100*mixDER, color = "MIXDER"), size = 1.2) + #MIXDER (black solid line)
    #geom_ribbon(aes(ymin = 100 * CI_lower, ymax = 100 * CI_upper), fill = "red", alpha = 0.3) + #MonteCarlo Ribbons (red ribbons)
  geom_point(data = example, aes(x= 100*d, y = 100*CA, color = ion, size = 100*error)) + #Actual Data for these for types of ions (colored points)
  #geom_line(aes(y = 100*simpleeffect, linetype = "SEA")) + #Simple Effect Additivity (black dashed line)
  
  scale_color_manual("Ions", breaks = c("Fe600", "O350", "O55", "Si170", "MIXDER"), values=c("blue", "red", "orange", "darkgreen", "black")) + #legend for Ions
  scale_size_continuous("Error", range = c(1,5)) +
  
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),  #Some random tweeks of the ggplot theme. Need someone more artistic to do this part :)
        panel.grid.major.y = element_line(colour = "grey80"),
        panel.grid.minor.y = element_line(linetype = "dashed", colour = "grey80"),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_line(),
        plot.title = element_text(size = 20, face = "bold"), 
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")
        ) +
  
  labs(x = "100*Dosage", y = "100*Chromosomal Aberrations", title = "Dose Effect Relationships", subtitle = "Fe600, Si170, O55, O350") #axis labels and plot title


###########################################################################
