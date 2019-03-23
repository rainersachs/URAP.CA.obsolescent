source('MonteCarlo.R')
library(ggplot2)
library(tidyr)
theme_update(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #This centers plot titles for ggplots

MC_results_4para_cov_full = monte_carlo(ions = c("Fe300", "Fe450","Fe600", "Ti300", "Si170", "Si260", "O55", "O350"), para = MC_4para, n = 500) #This outputs a list of MIXDER dataframe and a vector of the indexes at which there were convergence issues. 
MIXDER_4para_cov_full = MC_results_4para_cov_full[[1]] #This is the dataframe that will be passed into the graph


##################################################################
#Selecting Specific Data by Using Function subset().
#Example: Subset: d >= 0.05, CA > 0.01, Z = 22 and Comment section includes "2018"
subset(main_df, d >= 0.05 & CA > 0.01 & Z == 22 & grepl("2018", Comment))
#This returns a sub dataset with all the rows that satisfy the restrictions put in, which you can put in to the function "IDER_graph" as the only input. See example below.
IDER_graph(subset(main_df, d >= 0.05 & CA > 0.01 & Z == 22 & grepl("2018", Comment)))
##################################################################

#Plots

#Model Free Plot in Intro
p_intro_all = ggplot(main_df %>% filter(d <= 0.5), aes(x = 100*d, y = 100*CA)) + 
  geom_point(aes(col = ion)) +
  geom_smooth(method = "lm", se = F) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 100*BG_CA), size = 2, color = "red") +
  geom_errorbar(aes(x = 0, ymin = 0, ymax = 100*(BG_CA + BG_error)), width = 0.5, color = "red") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),  #Some random tweeks of the ggplot theme. Need someone more artistic to do this part :)
        panel.grid.major.y = element_line(colour = "grey80"),
        panel.grid.minor.y = element_line(linetype = "dashed", colour = "grey80"),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_line(),
        plot.title = element_text(size = 20, face = "bold"), 
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")) +
  labs(x = "100*Dosage", y = "100*Chromosomal Aberrations", title = "Visual Evidence for NTE", subtitle = "Raw Data")  #axis labels and plot title
p_intro_all

ggsave("Intro_Scatter.png", width = 10, height = 7)

squeezed_df = main_df %>% filter(d <= 0.5) %>% group_by(d) %>% summarise(p = sum(X.T)/sum(At.Risk), CA = p * 2.478, error = (sqrt((p*(1-p))/(sum(At.Risk))))*2.478) 
p_intro_squeezed = ggplot(squeezed_df, aes(x = 100*d, y = 100*CA)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  geom_errorbar(aes(ymin = 100*(CA - error), ymax = 100*(CA + error)), width = 0.5) + #Errorbars for each point from original data
  geom_errorbar(aes(x = 0, ymin = 0, ymax = 100*(BG_CA + BG_error)), width = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 100*BG_CA, color = "Value = 0.38"), size = 2.5) +
  
  scale_color_manual("Background", breaks = "Value = 0.38", values= "red2") + #legend for Ions
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),  #Some random tweeks of the ggplot theme. Need someone more artistic to do this part :)
      panel.grid.major.y = element_line(colour = "grey80"),
      panel.grid.minor.y = element_line(linetype = "dashed", colour = "grey80"),
      panel.grid.major.x = element_blank(),
      axis.ticks = element_line(),
      plot.title = element_text(size = 20, face = "bold"), 
      axis.title = element_text(size = 15, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = "100*Dosage", y = "100*Chromosomal Aberrations", title = "Visual Evidence for NTE", subtitle = "Raw Data Grouped by Dosage") #axis labels and plot title

p_intro_squeezed
ggsave("Intro_Scatter.png", width = 10, height = 7)

####################################

DER_df = main_df %>% filter(ion %in% c("Fe300", "Fe450","Fe600", "Ti300", "Si170", "Si260","O55", "O350"))
DER_df$CA = DER_df$CA - BG_CA #Taking out background
DER_df = filter(DER_df, (d > 0) & (d <= 0.5)) #Taking out 0 dosage points and cut at 0.5
MIXDER = MIXDER_4para_cov_full



#IDER Plot with ribbons (Si170 as example since it has the most data points)
Si170 = subset(DER_df, ion == "Si170" & d <= 0.5)
MC_results_4para_cov_Si170 = monte_carlo(ions = c("Si170"), para = MC_4para, n = 500, IDER = T) 
IDER_4para_cov_Si170 = MC_results_4para_cov_Si170[[1]] #This is the dataframe that will be passed into the graph
non_convergence_4para_cov_Si170 = MC_results_4para_cov_Si170[[2]] #This is a vector of indexes that did not converge. The length of this should be <= 3.
MC_results_4para_var_Si170 = monte_carlo(ions = c("Si170"), para = MC_4para, n = 500, IDER = T, cov = F) 
IDER_4para_var_Si170 = MC_results_4para_var_Si170[[1]] #This is the dataframe that will be passed into the graph
non_convergence_4para_var_Si170 = MC_results_4para_var_Si170[[2]] #This is a vector of indexes that did not converge. The length of this should be <= 3.

IDER_Si170 = IDER_4para_cov_Si170
p_Si170 = ggplot(data = IDER_Si170, aes(x = d * 100)) + #This example consists of one ion only: Si170
  
  geom_ribbon(data = IDER_4para_var_Si170, aes(ymin = 100 * CI_lower, ymax = 100 * CI_upper, fill = "Variance"), alpha = 0.3) + #MonteCarlo Ribbons without cov (blue ribbons)
  geom_ribbon(aes(ymin = 100 * CI_lower, ymax = 100 * CI_upper, fill = "Covariance")) + #MonteCarlo Ribbons with cov (red ribbons)
  geom_line(aes(y = 100 * IDER), size = 1) +  #IDER (colored solid lines)

  geom_point(data = Si170, aes(x= 100*d, y = 100*CA), size = 2) + #Actual Data for these for types of ions (colored points)
  geom_errorbar(data = Si170, aes(x = 100*d, ymin = 100*(CA - error), ymax = 100*(CA + error)), width = 0.5) + #Errorbars for each point from original data
  
  scale_fill_manual("MC CI", breaks = c("Covariance", "Variance"), values=c("yellow", "blue")) + #legend for cov vs var
  
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

p_Si170
ggsave("Si170IDER.png", width = 10, height = 7)
#IDER Plot with ribbons Fe600
Fe600 = subset(DER_df, ion == "Fe600" & d <= 0.5)
MC_results_4para_cov_Fe600 = monte_carlo(ions = c("Fe600"), para = MC_4para, n = 500, IDER = T) 
IDER_4para_cov_Fe600 = MC_results_4para_cov_Fe600[[1]] #This is the dataframe that will be passed into the graph
non_convergence_4para_cov_Fe600 = MC_results_4para_cov_Fe600[[2]] #This is a vector of indexes that did not converge. The length of this should be <= 3.
MC_results_4para_var_Fe600 = monte_carlo(ions = c("Fe600"), para = MC_4para, n = 500, IDER = T, cov = F) 
IDER_4para_var_Fe600 = MC_results_4para_var_Fe600[[1]] #This is the dataframe that will be passed into the graph
non_convergence_4para_var_Fe600 = MC_results_4para_var_Fe600[[2]] #This is a vector of indexes that did not converge. The length of this should be <= 3.

IDER_Fe600 = IDER_4para_cov_Fe600
p_Fe600 = ggplot(data = IDER_Fe600, aes(x = d * 100)) + #This example consists of one ion only: Fe600
  
  geom_ribbon(data = IDER_4para_var_Fe600, aes(ymin = 100 * CI_lower, ymax = 100 * CI_upper, fill = "Variance"), alpha = 0.3) + #MonteCarlo Ribbons without cov (blue ribbons)
  geom_ribbon(aes(ymin = 100 * CI_lower, ymax = 100 * CI_upper, fill = "Covariance")) + #MonteCarlo Ribbons with cov (red ribbons)
  geom_line(aes(y = 100 * IDER), size = 1) +  #IDER (colored solid lines)
  
  geom_point(data = Fe600, aes(x= 100*d, y = 100*CA), size = 2) + #Actual Data for these for types of ions (colored points)
  geom_errorbar(data = Fe600, aes(x = 100*d, ymin = 100*(CA - error), ymax = 100*(CA + error)), width = 0.5) + #Errorbars for each point from original data
  
  scale_fill_manual("MC CI", breaks = c("Covariance", "Variance"), values=c("yellow", "blue")) + #legend for cov vs var
  
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),  #Some random tweeks of the ggplot theme. Need someone more artistic to do this part :)
        panel.grid.major.y = element_line(colour = "grey80"),
        panel.grid.minor.y = element_line(linetype = "dashed", colour = "grey80"),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_line(),
        plot.title = element_text(size = 20, face = "bold"), 
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")
  ) +
  
  labs(x = "100*Dosage", y = "100*Chromosomal Aberrations", title = "Individual Dose Effect Relationship", subtitle = "Fe600") #axis labels and plot title
p_Fe600
ggsave("Fe600IDER.png", width = 10, height = 7)

#MIXDER Plot, using the same 4 ions.
names(MIXDER)[2] <- "mixDER"
p_mix = ggplot(data = MIXDER, aes(x = d * 100)) + #This example consists a mixture of 4 ions: "Fe600", "Si170", "O55", "O350". 
  geom_line(aes(y = 100 * Fe300, color = "Fe300"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100 * Fe450, color = "Fe450"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100 * Fe600, color = "Fe600"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100 * Ti300, color = "Ti300"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100 * O350, color = "O350"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100 * O55, color = "O55"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100 * Si170, color = "Si170"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100 * Si170, color = "Si260"), linetype = "dashed", size = 1) +
  geom_line(aes(y = 100*mixDER, color = "MIXDER"), linetype = "solid", size = 1.2) + #MIXDER (black solid line)
    #geom_ribbon(aes(ymin = 100 * CI_lower, ymax = 100 * CI_upper), fill = "red", alpha = 0.3) + #MonteCarlo Ribbons (red ribbons)
  geom_point(data = DER_df, aes(x= 100*d, y = 100*CA, color = ion, size = 1/sqrt(error), alpha = 1/sqrt(error))) + #Actual Data for these for types of ions (colored points)
  #geom_line(aes(y = 100*simpleeffect, linetype = "SEA")) + #Simple Effect Additivity (black dashed line)
  scale_alpha_continuous(range = c(0.3, 1)) +
  scale_size_continuous("Accuracy",range = c(1.5,4)) +
  
  scale_color_manual("Ions", breaks = c("MIXDER", "Fe300", "Fe450","Fe600",  "Ti300", "Si170", "Si260", "O55", "O350"), values=c("cyan2",  "cornflowerblue", "blue3", "black", "purple", "red", "orange","tan4", "darkgreen")) + #legend for Ions
  guides(alpha=FALSE) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),  #Some random tweeks of the ggplot theme. Need someone more artistic to do this part :)
        panel.grid.major.y = element_line(colour = "grey80"),
        panel.grid.minor.y = element_line(linetype = "dashed", colour = "grey80"),
        panel.grid.major.x = element_blank(),
        axis.ticks = element_line(),
        plot.title = element_text(size = 20, face = "bold"), 
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 12, face = "bold")
        ) +
  
  labs(x = "100*Dosage", y = "100*Chromosomal Aberrations", title = "Dose Effect Relationships", subtitle = "4 Parameter Model") #axis labels and plot title
p_mix
ggsave("DoseEffectRelationships.png", width = 10, height = 7)

###########################################################################
p_mix_15 = p_mix + xlim(c(0, 15)) + ylim(c(0, 3.5))
p_mix_15
ggsave("DoseEffectRelationships15.png", width = 10, height = 7)


#############################################################################################
#Functions 

# monte_carlo_graph <- function(MIXDER, model_name = "Model Name"){
#   #MIXDER is the output file of the monte carlo function. It could come from 6 different model * var combinations (3*2).
#   d = MIXDER$d
#   CA <- MIXDER$CA
#   simpleeffect <- MIXDER$simpleeffect
#   plot(x = d * 100, y = CA * 100, type = "l", col = "red")
#   lines(x = d * 100, y = simpleeffect * 100, col = "black", lty = 2, lwd = 0.5)
#   for(i in 6:ncol(MIXDER)){
#     lines(x = d * 100, y = 100*as.vector(MIXDER[,i]), col = "green")
#   }
#   lines(x= d*100 , y = MIXDER$CI_upper * 100, lty = 'dashed', col = 'red')
#   lines(x= d*100 , y = MIXDER$CI_lower * 100, lty = 'dashed', col = 'red')
#   polygon(c(d*100,rev(d*100)),c(MIXDER$CI_lower * 100, rev(MIXDER$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)
#   title(paste("Monte Carlo Plot for", model_name))
# }
# 
# IDER_graph <- function(data = main_df, ions = "O350", model = "4para", point = T, d = c(seq(0, 0.009, 0.001), seq(0.01, 0.5, by = 0.01))){
#   if(nrow(data) == nrow(main_df)){
#     CA = IDER(d = d, ions = ions, model = model)
#     data_filtered = data[(data$ion %in% ions) & (data$d <= 0.5),]
#     data_d = 100*data_filtered$d
#     data_CA = data_filtered$CA
#     n = length(ions)
#     names = ions[1]
#     if (n > 1){
#       for (i in 2:n){
#         names = paste(names, ",", ions[i])
#       }
#     }
#     print(max(data_CA))
#     plot(x = d * 100, y = CA, type = "l", col = "green", ylim = c(0, 1.5*max(data_CA)))
#     title(paste(model, "IDER Plot for", names))
#     if(point){#Option to add the scatter plot from the true data onto the IDER plot
#       points(data_d, data_CA)
#     }
#   }
#   else{
#     ions = unique(data$ion)
#     CA = IDER(d = d, ions = ions, model = model)
#     data_filtered = data
#     data_d = 100*data_filtered$d
#     data_CA = data_filtered$CA
#     n = length(ions)
#     names = ions[1]
#     if (n > 1){
#       for (i in 2:n){
#         names = paste(names, ",", ions[i])
#       }
#     }
#     plot(x = d * 100, y = CA, type = "l", col = "green", ylim = c(0, 1.5*max(data_CA)))
#     title(paste(model, "IDER Plot for", names))
#     if(point){#Option to add the scatter plot from the true data onto the IDER plot
#       points(data_d, data_CA)
#     }
#   }
# }

# Old Plot
# d = MIXDER$d
# CA <- MIXDER$CA
# simpleeffect <- MIXDER$simpleeffect
# plot(x = d * 100, y = CA * 100, type = "l", col = "red")
# lines(x = d * 100, y = simpleeffect * 100, col = "black", lty = 2, lwd = 0.5)
# for(i in 6:ncol(MIXDER)){
#   lines(x = d * 100, y = 100*as.vector(MIXDER[,i]), col = "green")
# }
# lines(x= d*100 , y = MIXDER$CI_upper * 100, lty = 'dashed', col = 'red')
# lines(x= d*100 , y = MIXDER$CI_lower * 100, lty = 'dashed', col = 'red')
# points(DER_df$d * 100, DER_df$CA * 100)
# polygon(c(d*100,rev(d*100)),c(MIXDER$CI_lower * 100, rev(MIXDER$CI_upper * 100)),col = rgb(1, 0, 0,0.5), border = NA)