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

##Model Free Plot in Intro
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
ggsave("Intro_Squeezed.png", width = 10, height = 7)

###################################################################################################################################################################

##Dose Effect Relationship Plots
DER_df = main_df %>% filter(ion %in% c("Fe300", "Fe450","Fe600", "Ti300", "Si170", "Si260","O55", "O350"))
DER_df$CA = DER_df$CA - BG_CA #Taking out background
DER_df = filter(DER_df, (d > 0) & (d <= 0.5)) #Taking out 0 dosage points and cut at 0.5
MIXDER = MIXDER_4para_cov_full



###Individual DER Plot with ribbons (Si170 as example since it has the most data points)
IDER_plot = function(data = DER_df, ion, d_cap = 0.5, var_ribbon = T, para = MC_4para){
  name = ion
  ion_data = subset(x = data, ion == name & d <= d_cap)
  MC_results_cov = monte_carlo(ions = ion, para = para, n = 500, IDER = T) 
  IDER_cov = MC_results_cov[[1]]
  
  p = ggplot(data = IDER_cov, aes(x = d * 100))
  if (var_ribbon){
    MC_results_var = monte_carlo(ions = ion, para = para, n = 500, IDER = T, cov = F) 
    IDER_var = MC_results_var[[1]]
    p = p + geom_ribbon(data = IDER_var, aes(ymin = 100 * CI_lower, ymax = 100 * CI_upper, fill = "Variance"), alpha = 0.3) #MonteCarlo Ribbons without cov (blue ribbons)
  }
  p = p + 
    geom_ribbon(aes(ymin = 100 * CI_lower, ymax = 100 * CI_upper, fill = "Covariance")) +
    geom_line(aes(y = 100 * IDER), size = 1) +  #IDER (colored solid lines)
    geom_point(data = ion_data, aes(x= 100*d, y = 100*CA), size = 2) + #Actual Data for these for types of ions (colored points)
    geom_errorbar(data = ion_data, aes(x = 100*d, ymin = 100*(CA - error), ymax = 100*(CA + error)), width = 0.5) + #Errorbars for each point from original data
    scale_fill_manual("MonteCarlo CI", breaks = c("Covariance", "Variance"), values=c("yellow", "blue")) + #legend for cov vs var
    
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),  #Some random tweeks of the ggplot theme. Need someone more artistic to do this part :)
          panel.grid.major.y = element_line(colour = "grey80"),
          panel.grid.minor.y = element_line(linetype = "dashed", colour = "grey80"),
          panel.grid.major.x = element_blank(),
          axis.ticks = element_line(),
          plot.title = element_text(size = 20, face = "bold"), 
          axis.title = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 12, face = "bold")
    ) +
    labs(x = "100*Dosage", y = "100*Chromosomal Aberrations", title = "Individual Dose Effect Relationship", subtitle = ion) #axis labels and plot title
  return(p)
}

p_Si170 = IDER_plot(ion = "Si170")
p_Si170
ggsave("Si170IDER.png", width = 10, height = 7)

p_Fe600 = IDER_plot(ion = "Fe600")
p_Fe600
ggsave("Fe600IDER.png", width = 10, height = 7)

p_Fe300 = IDER_plot(ion = "Fe300")
p_Fe300
ggsave("Fe300IDER.png", width = 10, height = 7)



##MIXDER Plot, using the same ions.
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

##Zoom in
p_mix_15 = p_mix + xlim(c(0, 15)) + ylim(c(0, 3.5))
p_mix_15
ggsave("DoseEffectRelationships15.png", width = 10, height = 7)

###########################################################################
#The SLI Plots
##Linear Model
DER_SLI_df = swift_light_df %>% mutate(CA = CA - BG_CA) %>% select(d, CA, ion, error)
slope = coef(model_sli_linear1) %>% as.numeric()

ggplot(DER_SLI_df, aes(x = 100*d)) +
  
  geom_abline(slope = slope) + 
  geom_point(aes(y = 100*CA*(CA > 0), color = ion, size = 1/sqrt(error), alpha = 1/sqrt(error))) + #Making negative values 0
  
  scale_color_manual("Ions", breaks = c("H250", "He250"), values = c("blue3", "red")) +
  scale_size_continuous("Accuracy",range = c(1.5,4)) + 
  scale_alpha_continuous(range = c(0.3, 1)) +
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
  lims(x = c(0,60)) +
  labs(x = "100*Dosage", y = "100*Chromosomal Aberrations", title = "Dose Effect Relationships", subtitle = "Swift Light Ions Linear Model")
ggsave("SLILinearDER.png", width = 10, height = 7)




#############################################################################################
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