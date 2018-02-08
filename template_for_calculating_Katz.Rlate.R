#   Filename: HGsynergyMain_merge2.R 
#   Purpose: Concerns radiogenic mouse HG tumorigenesis.

#    "16Srn" = Siranart et al."Mixed Beam Murine Harderian Gland Tumorigenesis: Predicted Dose-Effect Relationships if neither Synergism nor Antagonism Occurs." Radiat Res 186(6): 577-591 (2016).  

rm(list=ls())
#=========================== DATA START ===========================#
dfr <- data.frame( #  data used in 16Chang; includes data analyzed in .93Alp and .94Alp) # RKS -> Peter This is a completely different data set but the method of calculating Katz from L and Z and MeV/u is the same. Here a data frame is first set up and then the Katz calculation is added later. There are more systematic and elegant ways to program here but this way works.  
  dose.1 = c(0.2,0.4,0.6,1.2,2.4,3.2,5.1,7,0.05,0.1,0.15,0.2,0.4,0.8,1.6,0.05,0.1,0.2,0.4,0,0.1,0.2,0.4,0.8,1.6,0.4,0.8,1.6,3.2,0.05,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.8,0.1,0.2,0.4,0.04,0.08,0.16,0.32,0.033,0.066,0.13,0.26,0.52,.2, .4, .6),
  HG = c(0.091,0.045,0.101,0.169,0.347,0.431,0.667,0.623,0.156,0.215,0.232,0.307,0.325,0.554,0.649,0.123,0.145,0.207,0.31,0.026,0.083,0.25,0.39,0.438,0.424,0.093,0.195,0.302,0.292,0.109,0.054,0.066,0.128,0.286,0.183,0.167,0.396,0.536,0.192,0.234,0.317,0.092,0.131,0.124,0.297,0.082,0.088,0.146,0.236,0.371,.154,.132,.333), #  HG prevalence as defined in 16Chang
  NWeight = c(520,2048,1145,584,313,232,293,221,1162,877,455,409,374,223,320,742,661,347,131,6081,1091,251,244,191,131,645,255,199,111,649,378,973,833,201,468,381,197,109,496,257,185,1902,1063,884,350,1767,1408,874,299,261,322,206,67), #  nominal weight for weighted least squaresregression; see .93Alp. 
  index=c(rep(1,8),rep(0,17), rep(1,4),  rep(0,24)), #  index=0 for Z>3 ions, 1 otherwise. Not needed in some models
  L = c(rep(1.6,8), rep(193, 7), rep(250, 4), rep(195, 6), rep(0.4, 4), rep(25, 5), rep(464, 4), rep(193, 3),rep(70, 4), rep(100, 5), rep(953, 3)), #  L = LET = LET_infinity = stopping power (keV/micron)
  Z = c(rep(2, 8), rep(26, 17), rep(1, 4), rep(10, 5), rep(43, 4), rep(26, 3), rep(14, 4), rep(22, 5), rep(57, 3)), #  atomic number, charge in units of proton charge on fully ionized atomic nucleus, e.g. 2 for 2He4
   beta = c(rep("TBD", 53)), #  ion speed, relative to speed of light, calculated below
  MeVperu = c(rep(228, 8), rep(600, 7), rep(300, 4), rep(600, 6), rep(250, 4), rep(670, 5), rep(600, 4), rep(600, 3), rep(260, 4), rep(1000, 5), rep(593, 3)), #  Kinetic energy in MeV, divided by atomic mass, e.g. divided by 4u=4x931.5 MeV/c^2 for 2He4
  Katz = c(rep("TBD", 53)), #  for fully ionized nuclei, Katz's Z^2/beta^2, Calculated below. It is part of the Bethe Barkas Bloch equation for stopping power. Our calculations don't use Katz, but various similar calculations do.
  ion = c(rep("He4", 8), rep("Fe56", 17), rep("p", 4), rep("Ne20", 5), rep("Nb93", 4), rep("Fe56", 3), rep("Si28", 4), rep("Ti48", 5), rep("La139", 3)),
  comments = c(".93AlpLooksOK", rep("", 7), ".93AlplooksOK", rep("", 11), ".93Alp.no.iso", "not in 17Cuc (or 16Chang?)", rep("", 3), "16Chang all OK?", rep('', 24), ".94Alp","From graphs",'e.g. in 17Cuc')
) 

# Data for HG induced by photons from Cs-137 or Co-60 beta decay; from 16Chang (and calibration of LQ model)
ddd <- data.frame(dose.1 = c(0, 0.4, 0.8, 1.6, 3.2, 7, 0, .4, .8, .12, 1.6),
                  HG = c(.026, .048, .093, .137, .322, .462, .0497, .054, .067, .128, .202),
                  NWeight = c(6081.2, 4989.5, 1896.8, 981.1, 522.2, 205.2, 7474.1, 2877.6, 1423.7, 689.9, 514.9),
                  Nucleus = c(rep("Cobalt-60", 6), rep("Cesium-137", 5)),
                  Comments = c(rep("TBD", 11))
)
GeVu <- 0.001 * dfr[, "MeVperu"] #  convert to GeV/u for convenience in a calculation
dfr[, "Katz"] <- round(dfr[, "Z"] ^2 * (2.57 * GeVu ^2 + 4.781 * GeVu + 2.233) / (2.57 * GeVu ^2 + 4.781 * GeVu), 2) #  special relativistic calculation of Z^2/beta^2. The numerics include conversion from GeV to joules and from u to kg.
dfr[, "beta"] <- round(dfr[, "Z"] * sqrt(1 / dfr[, "Katz"]), 3) #  i.e. Z*sqrt(beta^2/Z^2) 
