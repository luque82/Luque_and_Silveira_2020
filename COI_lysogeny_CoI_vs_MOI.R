## This R script generates figures to visualize the relationship between COI and MOI for different parameters involved

### Working directory
setwd("/Volumes/GoogleDrive/My Drive/1_contributions/1_articles/1_in_progress/6_COI_and_lysogeny_low_cell_densities/wip/code_development/1_COI_vs_MOI_code")

### Packages
library("ggplot2")
library("svglite")

# Color-blind-friendly palettes, one with gray, and one with black (palettes source: https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/).
cbp1 = c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

####################
### CoI vs MOI
####################
### Environmental limits
dmin = 7.22E-10 ### ml/hr
dmax = 1.18E-6 ### ml/hr
taumin = 0.2*2.74 ### hr
taumax = 0.2*log(2)/(0.00412/24) ### hr
Bmin = 3.78E4 ## cells/ml 
Bmax = 7.65E9 ## cells/ml


### E. coli -- Lambda (default values)
d0 = 5.58E-7 ## infection rate constant (ml/h). (median from meta-analysis)
tau0 = 0.2*(30/60) ## lysogenic commitment period (h) (mean value from lab)
B0 = 5E8 ## bacterial concentration (cells/ml)
MOI0 = 2 ## Default MOI
COI0 = tau0*d0*B0*MOI0 ## average number of coinfections 

### Varying bacterial concentration
B1 = 1E5
B2 = 1E6
B3 = 1E7
B4 = 1E8
B5 = 1E9
B6 = 1E10
MOI = seq(0.01,200, by =0.01)
COIB1 = tau0*d0*B1*MOI 
COIB2 = tau0*d0*B2*MOI 
COIB3 = tau0*d0*B3*MOI 
COIB4 = tau0*d0*B4*MOI 
COIB5 = tau0*d0*B5*MOI 
COIB6 = tau0*d0*B6*MOI 
#HorLn = seq(2,2, length.out = length(MOI)) 
dfB = data.frame("MOI"=MOI, "COI1" = COIB1, "COI2" = COIB2, "COI3" = COIB3, "COI4" = COIB4, "COI5" = COIB5, "COI6" = COIB6)

## plot with lines
options(scipen = 999)
COIB = ggplot(dfB, aes(x=MOI, y=COI1)) + geom_smooth(method = "lm", se = FALSE, color= cbp1[1], linetype = "dashed", size = 0.5)
COIB = COIB + scale_x_log10() + scale_y_log10()
COIB = COIB + coord_cartesian(xlim = c(0.01,200), ylim = c(0.01,10))
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI2), method = "lm", se = FALSE, color= cbp1[2], linetype = "dashed", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI3), method = "lm", se = FALSE, color= cbp1[3], linetype = "dashed", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI4), method = "lm", se = FALSE, color= cbp1[4], linetype = "dashed", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI5), method = "lm", se = FALSE, color= cbp1[5], linetype = "solid", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI6), method = "lm", se = FALSE, color= cbp1[6], linetype = "dashed", size = 0.5)
COIB = COIB + geom_hline(yintercept = 2, linetype = "solid", color = "black", size = 0.5)
COIB = COIB + xlab("MOI") + ylab("CoI") + labs(title = "Bacterial concentration impact on CoI")
plot(COIB)
ggsave("COI_lysogeny_CoI_vs_MOI_for_B_lines.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("COI_lysogeny_CoI_vs_MOI_for_B_lines.png", device ="png", width = 8, height = 8, units = "cm")

#### Contour plot
df_MOI_COI_B.log10 = data.frame("MOI"=double(),"COI"=double(),"B"=double())
x = seq(log10(0.01),log10(100), length.out=100)
y = seq(log10(0.01),log10(100), length.out=100)
for(valx in x){
  for(valy in y){
    valz = (10^valy)/( (10^valx)*tau0*d0)
    nRow = data.frame("MOI" = valx,"COI" = valy,"B" = log10(valz))
    df_MOI_COI_B.log10 <- rbind(df_MOI_COI_B.log10,nRow)
  }
}
gpcont = ggplot(df_MOI_COI_B.log10, aes(MOI,COI,z =B)) +
geom_raster(aes(fill = B)) +
scale_fill_gradient(low="white", high=cbp1[8], guide="colorbar", limits=c(log10(Bmin),log10(Bmax))) +
geom_contour(colour = "white", linetype = "dashed") +
theme_classic() +
#theme(legend.position="none") +
geom_abline(intercept=log10(tau0*d0*B0),slope=1, colour = "white") +
geom_hline(yintercept = log10(2), linetype = "solid", color = "black", size = 0.5)
plot(gpcont)
ggsave("COI_lysogeny_CoI_vs_MOI_for_B_gradient.svg", device ="svg", width = 10, height = 8, units = "cm")
ggsave("COI_lysogeny_CoI_vs_MOI_for_B_gradient.png", device ="png", width = 10, height = 8, units = "cm")

#plot(MOI,COIB1, type = 'n', log='xy', xlim = c(0.01,200), ylim = c(1E-2,1E1))
#lines(MOI,COIB1)
#lines(MOI,COIB2)
#lines(MOI,COIB3)
#lines(MOI,COIB4)
#lines(MOI,COIB5)
#lines(MOI,COIB6)
#lines(MOI,HorLn)

### Varying commitment time
tau1 = 1E-3
tau2 = 1E-2
tau3 = 1E-1
tau4 = 1E0
tau5 = 1E1
tau6 = 1E2
MOI = seq(0.01,200, by =0.01)
COItau1 = tau1*d0*B0*MOI 
COItau2 = tau2*d0*B0*MOI 
COItau3 = tau3*d0*B0*MOI 
COItau4 = tau4*d0*B0*MOI 
COItau5 = tau5*d0*B0*MOI 
COItau6 = tau6*d0*B0*MOI 
HorLn = seq(2,2, length.out = length(MOI)) 
dfB = data.frame("MOI"=MOI, "COI1" = COItau1, "COI2" = COItau2, "COI3" = COItau3, "COI4" = COItau4, "COI5" = COItau5, "COI6" = COItau6)

## plot with lines
options(scipen = 999)
COIB = ggplot(dfB, aes(x=MOI, y=COI1)) + geom_smooth(method = "lm", se = FALSE, color= cbp1[1], linetype = "dashed", size = 0.5)
COIB = COIB + scale_x_log10() + scale_y_log10()
COIB = COIB + coord_cartesian(xlim = c(0.01,200), ylim = c(0.01,10))
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI2), method = "lm", se = FALSE, color= cbp1[2], linetype = "dashed", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI3), method = "lm", se = FALSE, color= cbp1[3], linetype = "solid", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI4), method = "lm", se = FALSE, color= cbp1[4], linetype = "dashed", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI5), method = "lm", se = FALSE, color= cbp1[5], linetype = "dashed", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI6), method = "lm", se = FALSE, color= cbp1[6], linetype = "dashed", size = 0.5)
COIB = COIB + geom_hline(yintercept = 2, linetype = "solid", color = "black", size = 0.5)
COIB = COIB + xlab("MOI") + ylab("CoI") + labs(title = "Commitment time impact on CoI")
plot(COIB)
ggsave("COI_lysogeny_CoI_vs_MOI_for_tau_lines.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("COI_lysogeny_CoI_vs_MOI_for_tau_lines.png", device ="png", width = 8, height = 8, units = "cm")

#### Contour plot
df_MOI_COI_tau.log10 = data.frame("MOI"=double(),"COI"=double(),"tau"=double())
x = seq(log10(0.01),log10(100), length.out=100)
y = seq(log10(0.01),log10(100), length.out=100)
for(valx in x){
  for(valy in y){
    valz = (10^valy)/( (10^valx)*B0*d0)
    nRow = data.frame("MOI" = valx,"COI" = valy,"tau" = log10(valz))
    df_MOI_COI_tau.log10 <- rbind(df_MOI_COI_tau.log10,nRow)
  }
}
gpcont = ggplot(df_MOI_COI_tau.log10, aes(MOI,COI,z = tau)) +
  geom_raster(aes(fill = tau)) +
  scale_fill_gradient(low="white", high=cbp1[7], guide="colorbar", limits=c(log10(taumin),log10(taumax))) +
  geom_contour(colour = "white", linetype = "dashed") +
  theme_classic() +
  #theme(legend.position="none") +
  geom_abline(intercept=log10(tau0*d0*B0),slope=1, colour = "white") +
  geom_hline(yintercept = log10(2), linetype = "solid", color = "black", size = 0.5)
plot(gpcont)
ggsave("COI_lysogeny_CoI_vs_MOI_for_tau_gradient.svg", device ="svg", width = 10, height = 8, units = "cm")
ggsave("COI_lysogeny_CoI_vs_MOI_for_tau_gradient.png", device ="png", width = 10, height = 8, units = "cm")


### Varying infection rate constant
d1 = 1E-11
d2 = 1E-10
d3 = 1E-9
d4 = 1E-8
d5 = 1E-7
d6 = 1E-6
MOI = seq(0.01,200, by =0.01)
COId1 = tau0*d1*B0*MOI 
COId2 = tau0*d2*B0*MOI 
COId3 = tau0*d3*B0*MOI 
COId4 = tau0*d4*B0*MOI 
COId5 = tau0*d5*B0*MOI 
COId6 = tau0*d6*B0*MOI 
#HorLn = seq(2,2, length.out = length(MOI)) 
dfB = data.frame("MOI"=MOI, "COI1" = COId1, "COI2" = COId2, "COI3" = COId3, "COI4" = COId4, "COI5" = COId5, "COI6" = COId6)

## plot with lines
options(scipen = 999)
COIB = ggplot(dfB, aes(x=MOI, y=COI1)) + geom_smooth(method = "lm", se = FALSE, color= cbp1[1], linetype = "dashed", size = 0.5)
COIB = COIB + scale_x_log10() + scale_y_log10()
COIB = COIB + coord_cartesian(xlim = c(0.01,200), ylim = c(0.01,10))
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI2), method = "lm", se = FALSE, color= cbp1[2], linetype = "dashed", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI3), method = "lm", se = FALSE, color= cbp1[3], linetype = "dashed", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI4), method = "lm", se = FALSE, color= cbp1[4], linetype = "solid", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI5), method = "lm", se = FALSE, color= cbp1[5], linetype = "dashed", size = 0.5)
COIB = COIB + geom_smooth(data = dfB, aes(x=MOI, y = COI6), method = "lm", se = FALSE, color= cbp1[6], linetype = "dashed", size = 0.5)
COIB = COIB + geom_hline(yintercept = 2, linetype = "solid", color = "black", size = 0.5)
COIB = COIB + xlab("MOI") + ylab("CoI") + labs(title = "Infection rate constant impact on CoI")
plot(COIB)
ggsave("COI_lysogeny_CoI_vs_MOI_for_d_lines.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("COI_lysogeny_CoI_vs_MOI_for_d_lines.png", device ="png", width = 8, height = 8, units = "cm")

#### Contour plot
df_MOI_COI_d.log10 = data.frame("MOI"=double(),"COI"=double(),"d"=double())
x = seq(log10(0.01),log10(100), length.out=100)
y = seq(log10(0.01),log10(100), length.out=100)
for(valx in x){
  for(valy in y){
    valz = (10^valy)/( (10^valx)*tau0*B0)
    nRow = data.frame("MOI" = valx,"COI" = valy,"d" = log10(valz))
    df_MOI_COI_d.log10 <- rbind(df_MOI_COI_d.log10,nRow)
  }
}
gpcont = ggplot(df_MOI_COI_d.log10, aes(MOI,COI,z =d)) +
  geom_raster(aes(fill = d)) +
  scale_fill_gradient(low="white", high=cbp1[4], guide="colorbar", limits=c(log10(dmin),log10(dmax))) +
  geom_contour(colour = "white", linetype = "dashed") +
  theme_classic() +
  #theme(legend.position="none") +
  geom_abline(intercept=log10(tau0*d0*B0),slope=1, colour = "white") +
  geom_hline(yintercept = log10(2), linetype = "solid", color = "black", size = 0.5)
plot(gpcont)
ggsave("COI_lysogeny_CoI_vs_MOI_for_d_gradient.svg", device ="svg", width = 10, height = 8, units = "cm")
ggsave("COI_lysogeny_CoI_vs_MOI_for_d_gradient.png", device ="png", width = 10, height = 8, units = "cm")


#plot(MOI,COId1, type = 'n', log='xy', xlim = c(0.01,200), ylim = c(1E-2,1E1))
#lines(MOI,COId1)
#lines(MOI,COId2)
#lines(MOI,COId3)
#lines(MOI,COId4)
#lines(MOI,COId5)
#lines(MOI,COId6)
#lines(MOI,HorLn)
