## This R script plots and analyzes MOI and COI for lambda

### Working directory
setwd("/Volumes/GoogleDrive/My Drive/1_contributions/1_articles/1_in_progress/6_COI_and_lysogeny_low_cell_densities/wip/code_development/2_lambda_analysis_MOI_and_COI")

### Packages
library("ggplot2")
library("svglite")

# Color-blind-friendly palettes, one with gray, and one with black (palettes source: http://jfly.iam.u-tokyo.ac.jp/color/).
cbp1 = c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Read MOI data
MOIdata = read.csv("literature_MOI_lysogeny_2020_01_15_mod.csv")

####################
### Lysogeny and MOI for lambda
####################
# Fitted model: Hill-Langmuir equation of order 2 (modified by prefactor)
x = MOIdata$MOI
y = MOIdata$Percent.Lysogeny
model_n1 = nls(y~a*x^1/(b+x^1), start = list(a=100, b=1))
model_n2 = nls(y~a*x^2/(b+x^2), start = list(a=100, b=1))
model_n3 = nls(y~a*x^3/(b+x^3), start = list(a=100, b=1))
predicted_df = data.frame(MOI=MOIdata$MOI, PL.predn1 = predict(model_n1,MOIdata$MOI), PL.predn2 = predict(model_n2,MOIdata$MOI),PL.predn3 = predict(model_n3,MOIdata$MOI))

# Output summary to file
sink('lambda_analysis_MOI_Hill_model.txt')
print('Probability of lysogeny vs. MOI: Hill-Langmuir model of order n=1')
summary(model_n1)
print('Probability of lysogeny vs. MOI: Hill-Langmuir model of order n=2')
summary(model_n2)
print('Probability of lysogeny vs. MOI: Hill-Langmuir model of order n=3')
summary(model_n3)
sink()

###################
### Intial COI vs initial MOI (lambda)
###################
### E. coli -- Lambda (default values)
d0 = 3E-7 ## infection rate constant (ml/h).
dmin = 5.93E-8 ## (ml/h)
dmax = 1.18E-6 ## (ml/h)
dmid = (dmin+dmax)/2
tau0 = 1/6 ## lysogenic commitment period (h)
taumin = 0.2*(20/60) ## lysogenic commitment period (h) (5 min)
taumax = 0.2*(40/60) ## lysogenic commitment period (h) (10 min)
taumid = (taumin+taumax)/2
B0 = 5E8 ## bacterial concentration (cells/ml)

MOI0 = seq(0.01, 100, by=0.01)
COI0 = tau0*d0*B0*MOI0
COImin = taumin*dmin*B0*MOI0
COImax = taumax*dmax*B0*MOI0
df = data.frame("MOI0"=MOI0, "COI0"=COI0, "COImin"=COImin, "COImax"=COImax)
df.pol = data.frame(x=c(0.01,0.01,100,100), y=c(taumin*d0*B0*0.01,taumax*d0*B0*0.01,taumin*d0*B0*100,taumax*d0*B0*100))

####################
### COI (general)
####################
COI = seq(0.01,100,by = 0.01)
PL = (1 - exp(-COI) - COI*exp(-COI)) ### probability of lysgoeny
PercL = 100*PL ### percentage of lysogeny

df.COI = data.frame("COI"=COI, "PercL"=PercL)
df.predicted.COI = data.frame(COI=df.COI$COI, PL.predn2 = predict(model_n2,data.frame(x=COI), interval = "prediction"))

####################
### MOI --> COI --> Probability (lambda)
####################
PLmin = (1 - exp(-COImin) - COImin*exp(-COImin)) ### probability of lysgoeny (tau min)
PLmax = (1 - exp(-COImax) - COImax*exp(-COImax)) ### probability of lysgoeny (tau max)
COImid = taumid*dmid*B0*MOI0
PLmid = (1 - exp(-COImid) - COImid*exp(-COImid)) ### probability of lysgoeny (tau mid)
PercLmin = 100*PLmin ### percentage of lysogeny
PercLmax = 100*PLmax
PercLmid = 100*PLmid

df.PercL = data.frame("MOI"=MOI0, "PercLmin"=PercLmin, "PercLmid"=PercLmid, "PercLmax"=PercLmax)

#######
# PLOTS
#######
# Plot Lysogeny and MOI for lambda
options(scipen = 999)
gglysogMOI = ggplot(MOIdata, aes(x=MOI,y=Percent.Lysogeny)) +
  geom_point() +
  #scale_x_log10()
  scale_x_log10() + scale_y_log10() +
  annotation_logticks() +
  coord_cartesian(xlim = c(0.01,100), ylim = c(0.01,100)) +
  #theme_classic() +
  theme_bw() +
  xlab("MOI") + ylab("Percentage of lysogeny (%)") + labs(title = "MOI impact on lysogeny in lambda") +
  geom_line(data = predicted_df, aes(x=MOI,y=PL.predn1), color = cbp1[1], linetype = "dashed") +
  geom_line(data = predicted_df, aes(x=MOI,y=PL.predn2), color = cbp1[1]) +
  geom_line(data = predicted_df, aes(x=MOI,y=PL.predn3), color = cbp1[1], linetype = "dotted")
  #geom_line(data=df.COI, aes(x=COI, y=PercL))
plot(gglysogMOI)
#ggsave("lambda_analysis_MOI.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lambda_analysis_MOI.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lambda_analysis_MOI.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lambda_analysis_MOI.png", device ="png", width = 8, height = 8, units = "cm")

### Plot: Intial COI vs initial MOI (lambda)
ggCOI0vsMOI0 = ggplot(df, aes(x=MOI0)) +
  #geom_line(aes(y=COI0)) +
  geom_line(aes(y=COImin), color = cbp1[1]) +
  geom_line(aes(y=COImax), color = cbp1[1]) +
  #geom_line(df, aes(x=MOI0, y = COImin)) +
  #geom_line(df, aes(x=MOI0, y = COImax)) +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks() +
  #theme_classic() +
  theme_bw() +
  geom_ribbon(aes(ymin=COImin,ymax=COImax), fill=cbp1[1], alpha = "0.2") +
  geom_hline(yintercept = 2, linetype = "solid", color = "black", size = 0.5) +
  xlab("Initial MOI") + ylab("Initial CoI") + labs(title = "Initial CoI in Lambda experiments")
plot(ggCOI0vsMOI0)
ggsave("lambda_analysis_MOI_and_COI.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lambda_analysis_MOI_and_COI.png", device ="png", width = 8, height = 8, units = "cm")

### PLOT COI (general)
ggCOI = ggplot(df.COI, aes(x=COI, y=PercL)) +
  geom_line() +
  geom_line(data = df.predicted.COI, aes(x=COI,y=PL.predn2), color = cbp1[1]) +
  #theme_classic()
  theme_bw() +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks() +
  xlab("COI and MOI") + ylab("Percentage of lysogeny (%)") + labs(title="Percentage for lysogeny for CoI")
plot(ggCOI)
#ggsave("lambda_analysis_COI.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lambda_analysis_COI.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lambda_analysis_COI.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lambda_analysis_COI.png", device ="png", width = 8, height = 8, units = "cm")

### Plot: MOI --> COI --> Probability (lambda)
ggPercL = ggplot(df.PercL, aes(x=MOI)) +
  geom_line(aes(y=PercLmid)) +
  geom_ribbon(aes(ymin=PercLmin,ymax=PercLmax), fill=cbp1[1], alpha = "0.2") +
  #theme_classic()
  theme_bw() +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks() +
  coord_cartesian(xlim = c(0.01,100), ylim = c(0.01,100)) +
  xlab("Constant MOI") + ylab("Percentage of lysogeny (%)") + labs(title="Predicted percentage of lysogeny (lambda)")
plot(ggPercL)
ggsave("lambda_analysis_PercLys_predicted.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lambda_analysis_PercLys_predicted.png", device ="png", width = 8, height = 8, units = "cm")


