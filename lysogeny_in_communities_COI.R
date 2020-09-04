## Probability of lysogeny thorugh phage coinfections in microbial communities

### Working directory
setwd("/Volumes/GoogleDrive/My Drive/1_contributions/1_articles/1_in_progress/6_COI_and_lysogeny_low_cell_densities/wip/code_development/4_lysogeny_in_communities")

### Packages
library("ggplot2")
library("svglite")
library(scales) # to access break formatting functions
library(tidyverse) # this includes the dplyr library for data frame subsets

# Color-blind-friendly palettes, one with gray, and one with black (palettes source: http://jfly.iam.u-tokyo.ac.jp/color/).
cbp1 = c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

### Read data
# Rank abundances
df.ranks.bacteria = read.csv("rank_abundance_analysis_bacteria.csv")
df.ranks.phage = read.csv("Viral_rank_abundance.csv")
nranks = nrow(df.ranks.bacteria)

### General paramters
nsamples = 100000 ### number of samples per community
PercL.breaks = c(0.1,1,10,25,50,75,100.01) ## Breaks for the percentage of lysogeny analysis in communities
PercLThr = 1 ## Threshold to analyze communities with percentage of lysogeny
PercLThr2 = 10 ## Threshold to analyze communities with percentage of lysogeny
PercLThr3 = 25 ## Threshold to analyze communities with percentage of lysogeny
PercL.ranks = 10 ## selected ranks for contribution analysis
ndata = 200 ## Number of points in scatter dataplots
vec.ndata.samples.rd = sample(seq(1,nsamples),ndata) ## random sampling of ndata within nsamples

#### Input parameters
## Gut
d.gut.min = 5.93E-8 ## adsorption rate ml/h
d.gut.max = 1.18E-6 ## adsorption rate ml/h
r.gut.min = log(2)/(7.27) ## minimum growth rate (1/hr)
r.gut.max = log(2)/(2.74) ## minimum growth rate (1/hr)
M.gut.min = 3.46E5 ## cells/ml
M.gut.max = 7.65E9 ## cells/ml
V.gut.min = 5.09E6 ## VLPs/ml
V.gut.max = 1.05E10 ## VLPs/ml
a.gut = 10^(5.35) ## V(M) power function relation: V = a*M^b
b.gut = 0.388 ## V(M) power function relation: V = a*M^b
rel.mean.gut = 1 ## relative mean to randomize virus-microbe relation
rel.sd.gut = 0.05 ## relative standard deviation to randomize virus-microbe relation
tau.frac.gut = 0.2 ## fraction of duplication time allocated to lysogenic commitment

### Marine
d.marine.min = 7.22E-10 ## adsorption rate ml/h
d.marine.max = 3.66E-7 ## adsorption rate ml/h
r.marine.min = (0.00412/24) ## minimum growth rate (1/hr)
r.marine.max = (0.3034/24) ## growth rate constant 1/h
M.marine.min = 3.78E4 ## cells/ml
M.marine.max = 6.75E6 ## cells/ml
V.marine.min = 1.4E5 ## VLPs/ml
V.marine.max = 3.86E7 ## VLPs/ml
a.marine = 10^(2.5) ## V(M) power function relation: V = a*M^b
b.marine = 0.712 ## V(M) power function relation: V = a*M^b
rel.mean.marine = 1 ## relative mean to randomize virus-microbe relation
rel.sd.marine = 0.05 ## relative standard deviation to randomize virus-microbe relation
tau.frac.marine = 0.2 ## fraction of duplication time allocated to lysogenic commitment

###########################
### GUT: sampling (direct network)
###########################
# Latin Hypercube sampling
# Hypercube
logd.vec.gut = seq(log10(d.gut.min),log10(d.gut.max), length.out = nsamples)
tau.gut.min = tau.frac.gut*log(2)/r.gut.max
tau.gut.max = tau.frac.gut*log(2)/r.gut.min
logtau.vec.gut = seq(log10(tau.gut.min),log10(tau.gut.max), length.out = nsamples)
logM.vec.gut = seq(log10(M.gut.min),log10(M.gut.max), length.out = nsamples)
#logV.vec.gut = log10(a.gut) + b.gut*logM.vec.gut ## Viral population based on power function
#logV.vec.gut = seq(log10(V.gut.min),log10(V.gut.max), length.out = nsamples)
V.frac.gut = df.ranks.phage$Gut
M.frac.gut = df.ranks.bacteria$gut
# Randomized vectors
logd.vec.gut.rd = sample(logd.vec.gut)
logtau.vec.gut.rd = sample(logtau.vec.gut)
logM.vec.gut.rd = logM.vec.gut
logV.vec.gut = log10(a.gut) + b.gut*logM.vec.gut.rd ## Viral population based on power function
#rel.vec.gut = abs(rnorm(nsamples,1,0.1)) ## relative random distribution vector. Avoid negative values
rel.vec.gut = abs(rnorm(nsamples,rel.mean.gut,rel.sd.gut)) ## relative random distribution vector. Avoid negative values
#logV.vec.gut.rd = logV.vec.gut
logV.vec.gut.rd = logV.vec.gut*rel.vec.gut
logV.vec.gut.rd = ifelse(logV.vec.gut.rd>log10(V.gut.max),log10(V.gut.max),logV.vec.gut.rd)
logV.vec.gut.rd = ifelse(logV.vec.gut.rd<log10(V.gut.min),log10(V.gut.min),logV.vec.gut.rd)
#logV.vec.gut.rd = sample(logV.vec.gut)


########################
### GUT: Community contribution to lysogeny
########################
COI.LHS.gut.com = data.frame(matrix("", ncol = nranks, nrow = nsamples))
PL.LHS.gut.com = data.frame(matrix("", ncol = nranks, nrow = nsamples))
Lys.LHS.gut.com = data.frame(matrix("", ncol = nranks, nrow = nsamples))
iiseq = seq(1,nranks)
for (ii in iiseq) {
  ## Rank per column. Rows are sampling
  COI.LHS.gut.com[ii] = (10^logd.vec.gut)*(10^logtau.vec.gut.rd)*(10^logV.vec.gut.rd)*V.frac.gut[ii]
  PL.LHS.gut.com[ii] = 1 - (1+COI.LHS.gut.com[ii])*exp(-COI.LHS.gut.com[ii])
  Lys.LHS.gut.com[ii] = (10^logM.vec.gut.rd)*M.frac.gut[ii]*PL.LHS.gut.com[ii]
}

Lys.LHS.gut.com.sum = rowSums(Lys.LHS.gut.com) ## Total concentration of lysogens
PercLys.LHS.gut.com = 100*Lys.LHS.gut.com.sum/(10^logM.vec.gut.rd) ## Percentage of lysogens formed

# Data frame
df.gut.com = data.frame("M"=10^logM.vec.gut.rd, "V"=10^logV.vec.gut.rd, "Lys"=Lys.LHS.gut.com.sum, "PercL"=PercLys.LHS.gut.com)

# Frequency of samples within lysogeny percentage breaks
iimax = length(PercL.breaks)
iiseq = seq(1,iimax)
nmax = nrow(df.gut.com)
cumul.breaks = matrix(0, ncol = 1, nrow = iimax)
for (ii in iiseq){
  br = PercL.breaks[ii]
  n = nrow(subset(df.gut.com, PercL < br))
  cumul.breaks[ii] = n
}
gut.freq.Perclys = matrix(0, ncol = 1, nrow = iimax) 
gut.freq.Perclys[1] = cumul.breaks[1]/nmax
iiseq = seq(2,iimax)
for (ii in iiseq){
  m = cumul.breaks[ii]
  n = cumul.breaks[ii-1]
  gut.freq.Perclys[ii] = (m-n)/nmax
}

# Data frame
df.gut.com.PercLys.breaks = data.frame("breaks"=PercL.breaks,"freq"=gut.freq.Perclys)

# Subset of communities above percentage threshold (ranks analysis)
index.lys = which(PercLys.LHS.gut.com>PercLThr)
subset.lys = Lys.LHS.gut.com.sum[index.lys] ### Total lysogeny for selected subset
subset.ranks.lys = Lys.LHS.gut.com[index.lys,1:PercL.ranks] ### Lysogeny concentration for the first ten ranks
frac.ranks.lys.gut.com = subset.ranks.lys/subset.lys
ave.lys.ranks.gut.com = matrix(0,nrow=PercL.ranks,ncol = 1)
sd.lys.ranks.gut.com = matrix(0,nrow=PercL.ranks,ncol = 1)
iiseq = seq(1,PercL.ranks)
matrix.tmp = data.matrix(frac.ranks.lys.gut.com)
for (ii in iiseq) {
  ave.lys.ranks.gut.com[ii] = mean(matrix.tmp[,ii])
  sd.lys.ranks.gut.com[ii] = sd(matrix.tmp[,ii])
}

logd.vec.gut.rd.lys = logd.vec.gut.rd[index.lys] ## subset of adsorption rates
logtau.vec.gut.rd.lys = logtau.vec.gut.rd[index.lys] ## subset of commitment times
logM.vec.gut.rd.lys = logM.vec.gut.rd[index.lys] ## subset of bacteria concentrations
logV.vec.gut.rd.lys = logV.vec.gut.rd[index.lys] ## subset of phage concentrations

# Subset of communities above percentage threshold (ranks analysis). Threshold 3
index.lys3 = which(PercLys.LHS.gut.com>PercLThr3)
subset.lys3 = Lys.LHS.gut.com.sum[index.lys3] ### Total lysogeny for selected subset
subset.ranks.lys3 = Lys.LHS.gut.com[index.lys3,1:PercL.ranks] ### Lysogeny concentration for the first ten ranks
frac.ranks.lys.gut.com3 = subset.ranks.lys3/subset.lys3
ave.lys.ranks.gut.com3 = matrix(0,nrow=PercL.ranks,ncol = 1)
sd.lys.ranks.gut.com3 = matrix(0,nrow=PercL.ranks,ncol = 1)
iiseq = seq(1,PercL.ranks)
matrix.tmp = data.matrix(frac.ranks.lys.gut.com3)
for (ii in iiseq) {
  ave.lys.ranks.gut.com3[ii] = mean(matrix.tmp[,ii])
  sd.lys.ranks.gut.com3[ii] = sd(matrix.tmp[,ii])
}

logd.vec.gut.rd.lys3 = logd.vec.gut.rd[index.lys3] ## subset of adsorption rates
logtau.vec.gut.rd.lys3 = logtau.vec.gut.rd[index.lys3] ## subset of commitment times
logM.vec.gut.rd.lys3 = logM.vec.gut.rd[index.lys3] ## subset of bacteria concentrations
logV.vec.gut.rd.lys3 = logV.vec.gut.rd[index.lys3] ## subset of phage concentrations


##########
## GUT: Contribution of dominant phages
#########
# COI
COI.LHS.gut.dom = (10^logd.vec.gut)*(10^logtau.vec.gut.rd)*(10^logV.vec.gut.rd)*V.frac.gut[1]
COI.LHS.gut.dom1 = COI.LHS.gut.dom
COI.LHS.gut.dom2 = (10^logd.vec.gut)*(10^logtau.vec.gut.rd)*(10^logV.vec.gut.rd)*V.frac.gut[2]
COI.LHS.gut.dom3 = (10^logd.vec.gut)*(10^logtau.vec.gut.rd)*(10^logV.vec.gut.rd)*V.frac.gut[3]

# Subset of communities above percentage threshold (COI)
index.lys = which(PercLys.LHS.gut.com>PercLThr)
COI.LHS.gut.dom1.lys = COI.LHS.gut.dom1[index.lys] ### rank 1
COI.LHS.gut.dom2.lys = COI.LHS.gut.dom2[index.lys] ### rank 2
COI.LHS.gut.dom3.lys = COI.LHS.gut.dom3[index.lys] ### rank 3

# Probability of lysogeny rank 1
PL.LHS.gut.dom1 = 1 - (1+COI.LHS.gut.dom1)*exp(-COI.LHS.gut.dom1)
# Concentration of lysogens rank 1
Lys.LHS.gut.dom1 = (10^logM.vec.gut.rd)*M.frac.gut[1]*PL.LHS.gut.dom1

# Data frame rank 1
df.gut.dom1 = data.frame("M"=10^logM.vec.gut.rd, "V"=10^logV.vec.gut.rd, "COI"=COI.LHS.gut.dom1, "PL"=PL.LHS.gut.dom1, "Lys"=Lys.LHS.gut.dom1)


###########################
### MARINE: sampling (direct network)
###########################
# Latin Hypercube sampling
# Hypercube
logd.vec.marine = seq(log10(d.marine.min),log10(d.marine.max), length.out = nsamples)
tau.marine.min = tau.frac.marine*log(2)/r.marine.max
tau.marine.max = tau.frac.marine*log(2)/r.marine.min
logtau.vec.marine = seq(log10(tau.marine.min),log10(tau.marine.max), length.out = nsamples)
logM.vec.marine = seq(log10(M.marine.min),log10(M.marine.max), length.out = nsamples)
#logV.vec.marine = log10(a.marine) + b.marine*logM.vec.marine ## Viral population based on power function
#logV.vec.marine = seq(log10(V.marine.min),log10(V.marine.max), length.out = nsamples)
V.frac.marine = df.ranks.phage$Marine
M.frac.marine = df.ranks.bacteria$marine
# Randomized vectors
logd.vec.marine.rd = logd.vec.marine
logtau.vec.marine.rd = sample(logtau.vec.marine)
logM.vec.marine.rd = sample(logM.vec.marine)
logV.vec.marine = log10(a.marine) + b.marine*logM.vec.marine.rd ## Viral population based on power function
rel.vec.marine = abs(rnorm(nsamples,rel.mean.marine,rel.sd.marine)) ## relative random distribution vector. Avoid negative vaues
logV.vec.marine.rd = logV.vec.marine*rel.vec.marine
logV.vec.marine.rd = ifelse(logV.vec.marine.rd>log10(V.marine.max),log10(V.marine.max),logV.vec.marine.rd)
logV.vec.marine.rd = ifelse(logV.vec.marine.rd<log10(V.marine.min),log10(V.marine.min),logV.vec.marine.rd)
#logV.vec.marine.rd = logV.vec.marine
#logV.vec.marine.rd = sample(logV.vec.marine)

########################
### MARINE: Community contribution to lysogeny
########################
COI.LHS.marine.com = data.frame(matrix("", ncol = nranks, nrow = nsamples))
PL.LHS.marine.com = data.frame(matrix("", ncol = nranks, nrow = nsamples))
Lys.LHS.marine.com = data.frame(matrix("", ncol = nranks, nrow = nsamples))
iiseq = seq(1,nranks)
for (ii in iiseq) {
  ## Rank per column. Rows are sampling
  COI.LHS.marine.com[ii] = (10^logd.vec.marine)*(10^logtau.vec.marine.rd)*(10^logV.vec.marine.rd)*V.frac.marine[ii]
  PL.LHS.marine.com[ii] = 1 - (1+COI.LHS.marine.com[ii])*exp(-COI.LHS.marine.com[ii])
  Lys.LHS.marine.com[ii] = (10^logM.vec.marine.rd)*M.frac.marine[ii]*PL.LHS.marine.com[ii]
}

Lys.LHS.marine.com.sum = rowSums(Lys.LHS.marine.com) ## Total concentration of lysogens
PercLys.LHS.marine.com = 100*Lys.LHS.marine.com.sum/(10^logM.vec.marine.rd) ## Percentage of lysogens formed

# Data frame
df.marine.com = data.frame("M"=10^logM.vec.marine.rd, "V"=10^logV.vec.marine.rd, "Lys"=Lys.LHS.marine.com.sum, "PercL"=PercLys.LHS.marine.com)

# Frequency of samples within lysogeny percentage breaks
iimax = length(PercL.breaks)
iiseq = seq(1,iimax)
nmax = nrow(df.marine.com)
cumul.breaks = matrix(0, ncol = 1, nrow = iimax)
for (ii in iiseq){
  br = PercL.breaks[ii]
  n = nrow(subset(df.marine.com, PercL < br))
  cumul.breaks[ii] = n
}
marine.freq.Perclys = matrix(0, ncol = 1, nrow = iimax) 
marine.freq.Perclys[1] = cumul.breaks[1]/nmax
iiseq = seq(2,iimax)
for (ii in iiseq){
  m = cumul.breaks[ii]
  n = cumul.breaks[ii-1]
  marine.freq.Perclys[ii] = (m-n)/nmax
}
# Data frame
df.marine.com.PercLys.breaks = data.frame("breaks"=PercL.breaks,"freq"=marine.freq.Perclys)

# Subset of communities above percentage threshold (ranks analysis)
index.lys = which(PercLys.LHS.marine.com>PercLThr)
subset.lys = Lys.LHS.marine.com.sum[index.lys] ### Total lysogeny for selected subset
subset.ranks.lys = Lys.LHS.marine.com[index.lys,1:PercL.ranks] ### Lysogeny concentration for the first ten ranks
frac.ranks.lys.marine.com = subset.ranks.lys/subset.lys
ave.lys.ranks.marine.com = matrix(0,nrow=PercL.ranks,ncol = 1)
sd.lys.ranks.marine.com = matrix(0,nrow=PercL.ranks,ncol = 1)
iiseq = seq(1,PercL.ranks)
matrix.tmp = data.matrix(frac.ranks.lys.marine.com)
for (ii in iiseq) {
  ave.lys.ranks.marine.com[ii] = mean(matrix.tmp[,ii])
  sd.lys.ranks.marine.com[ii] = sd(matrix.tmp[,ii])
}

logd.vec.marine.rd.lys = logd.vec.marine.rd[index.lys] ## subset of adsorption rates
logtau.vec.marine.rd.lys = logtau.vec.marine.rd[index.lys] ## subset of commitment times
logM.vec.marine.rd.lys = logM.vec.marine.rd[index.lys] ## subset of bacteria concentrations
logV.vec.marine.rd.lys = logV.vec.marine.rd[index.lys] ## subset of phage concentrations

# Subset of communities above percentage threshold (ranks analysis). Threshold 2
index.lys2 = which(PercLys.LHS.marine.com>PercLThr2)
subset.lys2 = Lys.LHS.marine.com.sum[index.lys2] ### Total lysogeny for selected subset
subset.ranks.lys2 = Lys.LHS.marine.com[index.lys2,1:PercL.ranks] ### Lysogeny concentration for the first ten ranks
frac.ranks.lys.marine.com2 = subset.ranks.lys2/subset.lys2
ave.lys.ranks.marine.com2 = matrix(0,nrow=PercL.ranks,ncol = 1)
sd.lys.ranks.marine.com2 = matrix(0,nrow=PercL.ranks,ncol = 1)
iiseq = seq(1,PercL.ranks)
matrix.tmp = data.matrix(frac.ranks.lys.marine.com2)
for (ii in iiseq) {
  ave.lys.ranks.marine.com2[ii] = mean(matrix.tmp[,ii])
  sd.lys.ranks.marine.com2[ii] = sd(matrix.tmp[,ii])
}

logd.vec.marine.rd.lys2 = logd.vec.marine.rd[index.lys2] ## subset of adsorption rates
logtau.vec.marine.rd.lys2 = logtau.vec.marine.rd[index.lys2] ## subset of commitment times
logM.vec.marine.rd.lys2 = logM.vec.marine.rd[index.lys2] ## subset of bacteria concentrations
logV.vec.marine.rd.lys2 = logV.vec.marine.rd[index.lys2] ## subset of phage concentrations

##########
## MARINE: Contribution of dominant phages
#########
# COI
COI.LHS.marine.dom = (10^logd.vec.marine)*(10^logtau.vec.marine.rd)*(10^logV.vec.marine.rd)*V.frac.marine[1]
COI.LHS.marine.dom1 = COI.LHS.marine.dom
COI.LHS.marine.dom2 = (10^logd.vec.marine)*(10^logtau.vec.marine.rd)*(10^logV.vec.marine.rd)*V.frac.marine[2]
COI.LHS.marine.dom3 = (10^logd.vec.marine)*(10^logtau.vec.marine.rd)*(10^logV.vec.marine.rd)*V.frac.marine[3]

# Subset of communities above percentage threshold (COI)
index.lys = which(PercLys.LHS.marine.com>PercLThr)
COI.LHS.marine.dom1.lys = COI.LHS.marine.dom1[index.lys] ### rank 1
COI.LHS.marine.dom2.lys = COI.LHS.marine.dom2[index.lys] ### rank 2
COI.LHS.marine.dom3.lys = COI.LHS.marine.dom3[index.lys] ### rank 3

# Probability of lysogeny rank 1
PL.LHS.marine.dom1 = 1 - (1+COI.LHS.marine.dom1)*exp(-COI.LHS.marine.dom1)
# Concentration of lysogens rank 1
Lys.LHS.marine.dom1 = (10^logM.vec.marine.rd)*M.frac.marine[1]*PL.LHS.marine.dom1

# Data frame rank 1
df.marine.dom1 = data.frame("M"=10^logM.vec.marine.rd, "V"=10^logV.vec.marine.rd, "COI"=COI.LHS.marine.dom1, "PL"=PL.LHS.marine.dom1, "Lys"=Lys.LHS.marine.dom1)

#########
## PLOTS 
#########

### Percentage of communities with different ranges of lysogeny
nbreaks = length(PercL.breaks)
brk = rep(c("0-0.1%","0.1-1%","1-10%","10-25%","25-50%","50-75%","75-100%"),2)
eco = c(rep("1.Marine",nbreaks),rep("2.Gut",nbreaks))
comperc = c(100*marine.freq.Perclys,100*gut.freq.Perclys)
df.breaks.lys = data.frame(brk,eco,comperc)

plot.PercL.breaks.bar = ggplot(df.breaks.lys, aes(fill=eco, y = comperc, x = brk)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("dodgerBlue4","darkgoldenrod1")) +
  theme_classic() +
  ylim(0,100) +
  xlab("Percentage of lysogens in the community (%)") +
  ylab("Percentage of communities (%)") +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  
plot(plot.PercL.breaks.bar)
ggsave("lysogeny_in_communities_percentage_breaks.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_percentage_breaks.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_percentage_breaks.png", device ="png", width = 8, height = 8, units = "cm")


## Contribution of top ranks to lysogeny
nbreaks = 10
brk = rep(c("1","2","3","4","5","6","7","8","9","9+1"),2)
eco = c(rep("1.Marine",nbreaks),rep("2.Gut",nbreaks))
rankfrac.ave = c(ave.lys.ranks.marine.com[1:nbreaks],ave.lys.ranks.gut.com[1:nbreaks])
rankfrac.sd = c(sd.lys.ranks.marine.com[1:nbreaks],sd.lys.ranks.gut.com[1:nbreaks])
df.ranks.lys = data.frame(brk,eco,rankfrac.ave,rankfrac.sd)
plot.ranks.lys.bar = ggplot(df.ranks.lys, aes(fill=eco, y = rankfrac.ave, x = brk)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("dodgerBlue4","darkgoldenrod1")) +
  geom_errorbar(aes(ymin=rankfrac.ave-rankfrac.sd, ymax=rankfrac.ave+rankfrac.sd), width=.2,
                position=position_dodge(.9)) +
  theme_classic() +
  ylim(0,1) +
  xlab("Phage-host pair rank") +
  ylab("Lysogenic contribution fraction") +
  theme(legend.position = "none")

plot(plot.ranks.lys.bar)
#ggsave("lysogeny_in_communities_lys_contr_ranks.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_lys_contr_ranks.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_lys_contr_ranks.png", device ="png", width = 8, height = 8, units = "cm")

## Density plot adsorption rates contributing to lysogeny
dvec = c(logd.vec.marine.rd.lys,logd.vec.gut.rd.lys)
nmax.marine = length(logd.vec.marine.rd.lys)
nmax.gut = length(logd.vec.gut.rd.lys)
eco = c(rep("1. Marine",nmax.marine),rep("2. Gut",nmax.gut))
df = data.frame("eco"=eco,d=dvec)
ggplot(df,aes(x=d, colour = eco)) +
  #geom_density() +
  stat_density(geom = "line", position = "identity") +
  scale_colour_manual(values=c("dodgerBlue4","darkgoldenrod1")) +
  geom_segment(aes(x=log10(d.marine.min),xend=log10(d.marine.max),y=-0.1,yend=-0.1),colour= "dodgerBlue4", size = 0.2, linetype = 2) +
  geom_segment(aes(x=log10(d.gut.min),xend=log10(d.gut.max),y=-0.15,yend=-0.15),colour= "darkgoldenrod1", size = 0.2, linetype = 2) +
  ylim(-0.2,1.5) +
  #expand_limits(y = -0.1)+
  theme_classic() +
  scale_x_continuous(breaks = log10(c(1E-10,1E-9,1E-8,1E-7,1E-6)),
                     labels = c("10^-10","10^-9","10^-8","10^-7","10^-6")) +
  xlab("Phage adsorption rate, d (ml/hr)") +
  ylab("Probability density distribution") +
  theme(legend.position = "none",text = element_text(size=6), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2))
#ggsave("lysogeny_in_communities_lys_adsorption_rates_dens.svg", device ="svg", width = 4.35, height = 4.35, units = "cm")
#ggsave("lysogeny_in_communities_lys_adsorption_rates_dens.pdf", device ="pdf", width = 4.35, height = 4.35, units = "cm")
#ggsave("lysogeny_in_communities_lys_adsorption_rates_dens.eps", device ="eps", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_lys_adsorption_rates_dens.png", device ="png", width = 4.35, height = 4.35, units = "cm")


## Density plot commitment times contributing to lysogeny
tauvec = c(logtau.vec.marine.rd.lys,logtau.vec.gut.rd.lys)
nmax.marine = length(logtau.vec.marine.rd.lys)
nmax.gut = length(logtau.vec.gut.rd.lys)
eco = c(rep("1. Marine",nmax.marine),rep("2. Gut",nmax.gut))
df = data.frame("eco"=eco,tau=tauvec)
ggplot(df,aes(x=tau, colour = eco)) +
  #geom_density() +
  stat_density(geom = "line", position = "identity") +
  scale_colour_manual(values=c("dodgerBlue4","darkgoldenrod1")) +
  geom_segment(aes(x=log10(tau.marine.min),xend=log10(tau.marine.max),y=-0.1,yend=-0.1),colour= "dodgerBlue4", size = 0.2, linetype = 2) +
  geom_segment(aes(x=log10(tau.gut.min),xend=log10(tau.gut.max),y=-0.15,yend=-0.15),colour= "darkgoldenrod1", size = 0.2, linetype = 2) +
  #xlim(log10(0.5),log10(50)) +
  ylim(-0.2,3) +
  theme_classic() +
  scale_x_continuous(breaks = log10(c(1E-1,1E0,1E1,1E2,1E3)),
                     labels = c("10^-1","10^0","10^1","10^2","10^3")) +
  xlab("Lysogenic commitment time, tau (hr)") +
  ylab("Probability density distribution") +
  theme(legend.position = "none",text = element_text(size=6), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2))
#ggsave("lysogeny_in_communities_lys_commitment_times_dens.svg", device ="svg", width = 4.35, height = 4.35, units = "cm")
#ggsave("lysogeny_in_communities_lys_commitment_times_dens.pdf", device ="pdf", width = 4.35, height = 4.35, units = "cm")
#ggsave("lysogeny_in_communities_lys_commitment_times_dens.eps", device ="eps", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_lys_commitment_times_dens.png", device ="png", width = 4.35, height = 4.35, units = "cm")

## Density plot bacteria concentrations contributing to lysogeny
Mvec = c(logM.vec.marine.rd.lys,logM.vec.gut.rd.lys)
nmax.marine = length(logM.vec.marine.rd.lys)
nmax.gut = length(logM.vec.gut.rd.lys)
eco = c(rep("1. Marine",nmax.marine),rep("2. Gut",nmax.gut))
df = data.frame("eco"=eco,M=Mvec)
ggplot(df,aes(x=M, colour = eco)) +
  #geom_density() +
  stat_density(geom = "line", position = "identity") +
  scale_colour_manual(values=c("dodgerBlue4","darkgoldenrod1")) +
  geom_segment(aes(x=log10(M.marine.min),xend=log10(M.marine.max),y=-0.1,yend=-0.1),colour= "dodgerBlue4", size = 0.2, linetype = 2) +
  geom_segment(aes(x=log10(M.gut.min),xend=log10(M.gut.max),y=-0.15,yend=-0.15),colour= "darkgoldenrod1", size = 0.2, linetype = 2) +
  ylim(-0.2,1.0) +
  theme_classic() +
  scale_x_continuous(limits = log10(c(1E4,2E10)),
                     breaks = log10(c(1E4,1E6,1E8,1E10)),
                     labels = c("10^4","10^6","10^8","10^10")) +
  xlab("Total bacterial concentration, B (cells/ml)") +
  ylab("Probability density distribution") +
  theme(legend.position = "none",text = element_text(size=6), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2))
#ggsave("lysogeny_in_communities_lys_bacteria_concentrations_dens.svg", device ="svg", width = 4.35, height = 4.35, units = "cm")
#ggsave("lysogeny_in_communities_lys_bacteria_concentrations_dens.pdf", device ="pdf", width = 4.35, height = 4.35, units = "cm")
#ggsave("lysogeny_in_communities_lys_bacteria_concentrations_dens.eps", device ="eps", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_lys_bacteria_concentrations_dens.png", device ="png", width = 4.35, height = 4.35, units = "cm")

## Density plot bacteria concentrations contributing to lysogeny
Vvec = c(logV.vec.marine.rd.lys,logV.vec.gut.rd.lys)
nmax.marine = length(logV.vec.marine.rd.lys)
nmax.gut = length(logV.vec.gut.rd.lys)
eco = c(rep("1. Marine",nmax.marine),rep("2. Gut",nmax.gut))
df = data.frame("eco"=eco,V=Vvec)
ggplot(df,aes(x=V, colour = eco)) +
  #geom_density() +
  stat_density(geom = "line", position = "identity") +
  scale_colour_manual(values=c("dodgerBlue4","darkgoldenrod1")) +
  geom_segment(aes(x=log10(V.marine.min),xend=log10(V.marine.max),y=-0.1,yend=-0.1),colour= "dodgerBlue4", size = 0.2, linetype = 2) +
  geom_segment(aes(x=log10(V.gut.min),xend=log10(V.gut.max),y=-0.15,yend=-0.15),colour= "darkgoldenrod1", size = 0.2, linetype = 2) +
  ylim(-0.2,1.5) +
  theme_classic() +
  scale_x_continuous(limits = log10(c(3E4,2E10)),
                     breaks = log10(c(1E5,1E6,1E7,1E8,1E9,1E10)),
                     labels = c("10^5","10^6","10^7","10^8","10^9","10^10")) +
  xlab("Total phage concentration, P (phages/ml)") +
  ylab("Probability density distribution") +
  theme(legend.position = "none",text = element_text(size=6), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2))
#ggsave("lysogeny_in_communities_lys_phage_concentrations_dens.svg", device ="svg", width = 4.35, height = 4.35, units = "cm")
#ggsave("lysogeny_in_communities_lys_phage_concentrations_dens.pdf", device ="pdf", width = 4.35, height = 4.35, units = "cm")
#ggsave("lysogeny_in_communities_lys_phage_concentrations_dens.eps", device ="eps", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_lys_phage_concentrations_dens.png", device ="png", width = 4.35, height = 4.35, units = "cm")

## Density plot COI top ranks gut
COIvec = log10(c(COI.LHS.gut.dom1.lys,COI.LHS.gut.dom2.lys,COI.LHS.gut.dom3.lys))
nmax = length(COI.LHS.gut.dom1.lys)
rank = c(rep("Rank 1",nmax),rep("Rank 2",nmax),rep("Rank 3",nmax))
df = data.frame("rank"=rank,COI=COIvec)
#ggplot(df,aes(x=COI, colour = rank)) +
ggplot(df,aes(x=COI, linetype = rank, color = "darkgoldenrod1")) +  
  #geom_density() +
  stat_density(geom = "line", position = "identity") +
  scale_colour_manual(values=c("darkgoldenrod1","darkgoldenrod1","darkgoldenrod1")) +
  scale_linetype_manual(values=c("solid","dashed","dotted")) +
  geom_segment(aes(x=log10(2),xend=log10(2),y=0,yend=2),colour= "black", size = 0.01) +
  ylim(0.0,2) +
  #xlim(-1,2) +
  theme_classic() +
  scale_x_continuous(limits = log10(c(1E-1,1E2)),
                     breaks = log10(c(1E-1,3E-1,1E0,3E0,1E1,3E1,1E2)),
                     labels = c("10^-1","3.10^-1","10^0","3.10^0","10^1","3.10^1","10^2")) +
  xlab("Average phage coinfections, COI") +
  ylab("Probability density distribution") +
  theme(legend.position = "none")
#ggsave("lysogeny_in_communities_lys_COI_ranks_gut_dens.svg", device ="svg", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_lys_COI_ranks_gut_dens.pdf", device ="pdf", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_lys_COI_ranks_gut_dens.png", device ="png", width = 4.35, height = 4.35, units = "cm")


## Density plot COI top ranks marine
COIvec = log10(c(COI.LHS.marine.dom1.lys,COI.LHS.marine.dom2.lys,COI.LHS.marine.dom3.lys))
nmax = length(COI.LHS.marine.dom1.lys)
rank = c(rep("Rank 1",nmax),rep("Rank 2",nmax),rep("Rank 3",nmax))
df = data.frame("rank"=rank,COI=COIvec)
#ggplot(df,aes(x=COI, colour = rank)) +
ggplot(df,aes(x=COI, linetype = rank, color = "dodgerBlue4")) +  
  #geom_density() +
  stat_density(geom = "line", position = "identity") +
  scale_colour_manual(values=c("dodgerBlue4","dodgerBlue4","dodgerBlue4")) +
  scale_linetype_manual(values=c("solid","dashed","dotted")) +
  geom_segment(aes(x=log10(2),xend=log10(2),y=0,yend=2),colour= "black", size = 0.01) +
  ylim(0,2) +
  #xlim(-1,2) +
  theme_classic() +
  scale_x_continuous(limits = log10(c(1E-1,1E2)),
                     breaks = log10(c(1E-1,3E-1,1E0,3E0,1E1,3E1,1E2)),
                     labels = c("10^-1","3.10^-1","10^0","3.10^0","10^1","3.10^1","10^2")) +
  xlab("Average phage coinfections, COI") +
  ylab("Probability density distribution") +
  theme(legend.position = "none")
#ggsave("lysogeny_in_communities_lys_COI_ranks_marine_dens.svg", device ="svg", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_lys_COI_ranks_marine_dens.pdf", device ="pdf", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_lys_COI_ranks_marine_dens.png", device ="png", width = 4.35, height = 4.35, units = "cm")

# Plot Percentage of Lysogens in the population vs M (reduced axes)
df.marine.com.ndata.rd = df.marine.com[vec.ndata.samples.rd,]
df.gut.com.ndata.rd = df.gut.com[vec.ndata.samples.rd,]
ggLysPercMmarinegutcomred = ggplot() +
  geom_point(data=df.marine.com.ndata.rd, aes(x=M, y=Lys/M*100), color="dodgerBlue4") +
  #geom_smooth(data=df.marine.com.ndata.rd,aes(x=M, y=Lys/M*100), color="dodgerBlue4") +
  geom_smooth(data=df.marine.com,aes(x=M, y=Lys/M*100), color="dodgerBlue4") +
  geom_point(data=df.gut.com.ndata.rd, aes(x=M, y=Lys/M*100), color = "darkgoldenrod1") +
  #geom_smooth(data=df.gut.com.ndata.rd,aes(x=M, y=Lys/M*100), color = "darkgoldenrod1") +
  geom_smooth(data=df.gut.com,aes(x=M, y=Lys/M*100), color = "darkgoldenrod1") +
  ylim(0,100) +
  #scale_x_log10() +
  #scale_y_log10(limits-c(0.01,100)) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits=c(5E4,1E10)) +
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x)), limits=c(0.01,100)) +
  annotation_logticks(sides = "b") +
  xlab("Bacteria (1/ml)") + ylab("Percentage of lysogens (%)") + labs(title="Percentage of lysogeny") + 
  #theme_classic()
  theme(legend.position = "none",text = element_text(size=6), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2))
plot(ggLysPercMmarinegutcomred)
ggsave("lysogeny_in_communities_LysPer_vs_M_comred.svg", device ="svg", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_LysPer_vs_M_comred.eps", device ="eps", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_LysPer_vs_M_comred.pdf", device ="pdf", width = 4.35, height = 4.35, units = "cm")
ggsave("lysogeny_in_communities_LysPer_vs_M_comred.png", device ="png", width = 4.35, height = 4.35, units = "cm")


# Plot Concentration of Lysogens in the population vs M
df.marine.com.ndata.rd = df.marine.com[vec.ndata.samples.rd,]
df.gut.com.ndata.rd = df.gut.com[vec.ndata.samples.rd,]
ggLysMmarinegutcom = ggplot() +
  geom_point(data=df.marine.com.ndata.rd, aes(x=M, y=Lys), color="dodgerBlue4", size = 0.5) +
  geom_smooth(data=df.marine.com,aes(x=M, y=Lys), color="dodgerBlue4") +
  geom_point(data=df.gut.com.ndata.rd, aes(x=M, y=Lys), color = "darkgoldenrod1", size = 0.5) +
  geom_smooth(data=df.gut.com,aes(x=M, y=Lys), color = "darkgoldenrod1") +
  #scale_x_log10() +
  #scale_y_log10() +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e4,1e10)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(size = 0.2) +
  xlab("Bacteria (cells/ml)") + ylab("Concentration of lysogens (1/ml)") + labs(title="Lysogeny for full community") + 
  #theme_classic()
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
plot(ggLysMmarinegutcom)
#ggsave("lysogeny_in_communities_Lys_vs_B_com.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_Lys_vs_B_com.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_Lys_vs_B_com.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_Lys_vs_B_com.png", device ="png", width = 8, height = 8, units = "cm")

# Plot Concentration of Lysogens in the population vs V
df.marine.com.ndata.rd = df.marine.com[vec.ndata.samples.rd,]
df.gut.com.ndata.rd = df.gut.com[vec.ndata.samples.rd,]
ggLysMmarinegutcom = ggplot() +
  geom_point(data=df.marine.com.ndata.rd, aes(x=V, y=Lys), color="dodgerBlue4", size = 0.5) +
  geom_smooth(data=df.marine.com,aes(x=V, y=Lys), color="dodgerBlue4") +
  geom_point(data=df.gut.com.ndata.rd, aes(x=V, y=Lys), color = "darkgoldenrod1", size = 0.5) +
  geom_smooth(data=df.gut.com,aes(x=V, y=Lys), color = "darkgoldenrod1") +
  #scale_x_log10() +
  #scale_y_log10() +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e4,1e10)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(size = 0.2) +
  xlab("Phages (phages/ml)") + ylab("Concentration of lysogens (1/ml)") + labs(title="Lysogeny for full community") + 
  #theme_classic()
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
plot(ggLysMmarinegutcom)
#ggsave("lysogeny_in_communities_Lys_vs_P_com.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_Lys_vs_P_com.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_Lys_vs_P_com.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_Lys_vs_P_com.png", device ="png", width = 8, height = 8, units = "cm")

# Plot Percentage of Lysogens in the population vs M
df.marine.com.ndata.rd = df.marine.com[vec.ndata.samples.rd,]
df.gut.com.ndata.rd = df.gut.com[vec.ndata.samples.rd,]
ggLysMmarinegutcom = ggplot() +
  geom_point(data=df.marine.com.ndata.rd, aes(x=M, y=Lys/M*100), color="dodgerBlue4", size = 0.5) +
  geom_smooth(data=df.marine.com,aes(x=M, y=Lys/M*100), color="dodgerBlue4") +
  geom_point(data=df.gut.com.ndata.rd, aes(x=M, y=Lys/M*100), color = "darkgoldenrod1", size = 0.5) +
  geom_smooth(data=df.gut.com,aes(x=M, y=Lys/M*100), color = "darkgoldenrod1") +
  #scale_x_log10() +
  #scale_y_log10() +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e4,1e10)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(size = 0.2) +
  xlab("Bacteria (cells/ml)") + ylab("Percentage of lysogens (%)") + labs(title="Lysogeny for full community") + 
  #theme_classic()
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
plot(ggLysMmarinegutcom)
#ggsave("lysogeny_in_communities_LysPerc_vs_B_com.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_B_com.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_B_com.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_B_com.png", device ="png", width = 8, height = 8, units = "cm")

# Plot Percentage of Lysogens in the population vs V
df.marine.com.ndata.rd = df.marine.com[vec.ndata.samples.rd,]
df.gut.com.ndata.rd = df.gut.com[vec.ndata.samples.rd,]
ggLysMmarinegutcom = ggplot() +
  geom_point(data=df.marine.com.ndata.rd, aes(x=V, y=Lys/M*100), color="dodgerBlue4", size = 0.5) +
  geom_smooth(data=df.marine.com,aes(x=V, y=Lys/M*100), color="dodgerBlue4") +
  geom_point(data=df.gut.com.ndata.rd, aes(x=V, y=Lys/M*100), color = "darkgoldenrod1", size = 0.5) +
  geom_smooth(data=df.gut.com,aes(x=V, y=Lys/M*100), color = "darkgoldenrod1") +
  #scale_x_log10() +
  #scale_y_log10() +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e4,1e10)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(size = 0.2) +
  xlab("Phages (phages/ml)") + ylab("Percentage of lysogens (%)") + labs(title="Lysogeny for full community") + 
  #theme_classic()
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
plot(ggLysMmarinegutcom)
#ggsave("lysogeny_in_communities_LysPerc_vs_P_com.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_P_com.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_P_com.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_P_com.png", device ="png", width = 8, height = 8, units = "cm")

# Plot COI vs M for rank1
df.marine.dom1.ndata.rd = df.marine.dom1[vec.ndata.samples.rd,]
df.gut.dom1.ndata.rd = df.gut.dom1[vec.ndata.samples.rd,]
ggplot() +
  geom_point(data=df.marine.dom1.ndata.rd, aes(x=M, y=COI), color="dodgerBlue4", size = 0.5) +
  geom_smooth(data=df.marine.dom1,aes(x=M, y=COI), color="dodgerBlue4") +
  geom_point(data=df.gut.dom1.ndata.rd, aes(x=M, y=COI), color = "darkgoldenrod1", size = 0.5) +
  geom_smooth(data=df.gut.dom1,aes(x=M, y=COI), color = "darkgoldenrod1") +
  #scale_x_log10() +
  #scale_y_log10() +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e4,1e10)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e-4,1e2)) +
  annotation_logticks(size = 0.2) +
  xlab("Bacteria (cells/ml)") + ylab("COI") + labs(title="Lysogeny for full community") + 
  #theme_classic()
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
#ggsave("lysogeny_in_communities_COI_vs_M_rank1.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_COI_vs_M_rank1.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_COI_vs_M_rank1.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_COI_vs_M_rank1.png", device ="png", width = 8, height = 8, units = "cm")

# Plot COI vs V for rank1
df.marine.dom1.ndata.rd = df.marine.dom1[vec.ndata.samples.rd,]
df.gut.dom1.ndata.rd = df.gut.dom1[vec.ndata.samples.rd,]
ggplot() +
  geom_point(data=df.marine.dom1.ndata.rd, aes(x=V, y=COI), color="dodgerBlue4", size = 0.5) +
  geom_smooth(data=df.marine.dom1,aes(x=V, y=COI), color="dodgerBlue4") +
  geom_point(data=df.gut.dom1.ndata.rd, aes(x=V, y=COI), color = "darkgoldenrod1", size = 0.5) +
  geom_smooth(data=df.gut.dom1,aes(x=V, y=COI), color = "darkgoldenrod1") +
  #scale_x_log10() +
  #scale_y_log10() +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e4,1e10)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e-4,1e2)) +
  annotation_logticks(size = 0.2) +
  xlab("Phages (phages/ml)") + ylab("COI") + labs(title="Lysogeny for full community") + 
  #theme_classic()
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
#ggsave("lysogeny_in_communities_COI_vs_V_rank1.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_COI_vs_V_rank1.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_COI_vs_V_rank1.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_COI_vs_V_rank1.png", device ="png", width = 8, height = 8, units = "cm")

# Plot Percentage of Lysogens from Rank vs M
df.marine.dom1.ndata.rd = df.marine.dom1[vec.ndata.samples.rd,]
df.gut.dom1.ndata.rd = df.gut.dom1[vec.ndata.samples.rd,]
ggLysMmarinegutcom = ggplot() +
  geom_point(data=df.marine.dom1.ndata.rd, aes(x=M, y=Lys/M*100), color="dodgerBlue4", size = 0.5) +
  geom_smooth(data=df.marine.dom1,aes(x=M, y=Lys/M*100), color="dodgerBlue4") +
  geom_point(data=df.gut.dom1.ndata.rd, aes(x=M, y=Lys/M*100), color = "darkgoldenrod1", size = 0.5) +
  geom_smooth(data=df.gut.dom1,aes(x=M, y=Lys/M*100), color = "darkgoldenrod1") +
  #scale_x_log10() +
  #scale_y_log10() +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e4,1e10)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(size = 0.2) +
  xlab("Bacteria (cells/ml)") + ylab("Percentage of lysogens from rank 1 (%)") + labs(title="Lysogeny for Rank 1") + 
  #theme_classic()
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
plot(ggLysMmarinegutcom)
#ggsave("lysogeny_in_communities_LysPerc_vs_B_rank1.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_B_rank1.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_B_rank1.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_B_rank1.png", device ="png", width = 8, height = 8, units = "cm")


# Plot Percentage of Lysogens from Rank vs V
df.marine.dom1.ndata.rd = df.marine.dom1[vec.ndata.samples.rd,]
df.gut.dom1.ndata.rd = df.gut.dom1[vec.ndata.samples.rd,]
ggLysMmarinegutcom = ggplot() +
  geom_point(data=df.marine.dom1.ndata.rd, aes(x=V, y=Lys/M*100), color="dodgerBlue4", size = 0.5) +
  geom_smooth(data=df.marine.dom1,aes(x=V, y=Lys/M*100), color="dodgerBlue4") +
  geom_point(data=df.gut.dom1.ndata.rd, aes(x=V, y=Lys/M*100), color = "darkgoldenrod1", size = 0.5) +
  geom_smooth(data=df.gut.dom1,aes(x=V, y=Lys/M*100), color = "darkgoldenrod1") +
  #scale_x_log10() +
  #scale_y_log10() +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e4,1e10)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(size = 0.2) +
  xlab("Phages (phages/ml)") + ylab("Percentage of lysogens from rank 1 (%)") + labs(title="Lysogeny for Rank 1") + 
  #theme_classic()
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
plot(ggLysMmarinegutcom)
#ggsave("lysogeny_in_communities_LysPerc_vs_V_rank1.svg", device ="svg", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_V_rank1.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_V_rank1.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_LysPerc_vs_V_rank1.png", device ="png", width = 8, height = 8, units = "cm")

# Plot Percentage of lysogeny and VMR
df.filt <- df.marine.com %>% filter(PercL > 1)
ggLysVMRmar = ggplot()+
  geom_point(data=df.filt, aes(x=V/M, y=PercL), color="dodgerBlue4", size = 0.05) +
  geom_vline(xintercept=1) +
  ylim(0,100) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e-1,5e2)) +
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e0,1e2)) +
  annotation_logticks(side = 'b',size = 0.2) +
  xlab("virus-to-microbe ratio (VMR)") + ylab("Percentage of lysogeny (>1%)") + labs(title="Lysogney and VMR (marine)") +
  theme_classic() +
  #theme_bw()
  theme(text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
plot(ggLysVMRmar)
ggsave("lysogeny_in_communities_marine_VMR.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_marine_VMR.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_marine_VMR.png", device ="png", width = 8, height = 8, units = "cm")

df.filt <- df.gut.com %>% filter(PercL > 1)
ggLysVMRgut = ggplot()+
  geom_point(data=df.filt, aes(x=V/M, y=PercL), color = "darkgoldenrod1", size = 0.01) +
  geom_vline(xintercept=1) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  ylim(0,100) +
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(1e0,1e2)) +
  annotation_logticks(side = 'b',size = 0.2) +
  xlab("virus-to-microbe ratio (VMR)") + ylab("Percentage of lysogeny (>1%)") + labs(title="Lysogney and VMR (marine)") +
  theme_classic() +
  #theme_bw()
  theme(text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
plot(ggLysVMRgut)
ggsave("lysogeny_in_communities_gut_VMR.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_gut_VMR.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_gut_VMR.png", device ="png", width = 8, height = 8, units = "cm")

### Percentage of communities with VMR<1 (among communities with lysogeny > 1%)
df = data.frame('VMR'=(df.gut.com$V)/df.gut.com$M,'PercL'=df.gut.com$PercL)
nVMRl1 = nrow(df %>% filter(PercL>1 & VMR<1))
nPercm1 = nrow(df %>% filter(PercL>1))
PercVMRl1 = nVMRl1/nPercm1*100

# Plot correlation abundance and commitment time
df = data.frame("M"=df.marine.com$M, "tau"=10^logtau.vec.marine.rd, "PercL"=df.marine.com$PercL)
df.filt <- df %>% filter(PercL > 1)
ggLysDensTauMar = ggplot()+
  geom_point(data=df.filt, aes(x=M, y=tau, colour = PercL), size = 0.5) +
  scale_colour_gradient(limits=c(0,100),breaks=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(5e4,1e7)) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(size = 0.2) +
  xlab("Bacteria (1/ml)") + ylab("Commitment time, tau (h)") + labs(title="Combined effect %L > 1 %(marine)") +
  theme_classic() +
  #theme_bw()
  theme(text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
plot(ggLysDensTauMar)

ggsave("lysogeny_in_communities_marine_B_and_tau.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_marine_B_and_tau.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_marine_B_and_tau.png", device ="png", width = 8, height = 8, units = "cm")


df = data.frame("M"=df.gut.com$M, "tau"=10^logtau.vec.gut.rd, "PercL"=df.gut.com$PercL)
df.filt <- df %>% filter(PercL > 1)
ggLysDensTauGut = ggplot()+
  geom_point(data=df.filt, aes(x=M, y=tau, colour = PercL), size = 0.25) +
  scale_colour_gradient(low="darkred", high="yellow",
                        limits=c(0,100),breaks=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + #, limits = c(5e4,1e7)) +
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x)), limits = c(4e-1,2e0)) +
  ylim(0.5,1.5) +
  annotation_logticks(sides = 'b', size = 0.2) +
  xlab("Bacteria (1/ml)") + ylab("Commitment time, tau (h)") + labs(title="Combined effect %L > 1 %(marine)") +
  theme_classic() +
  #theme_bw()
  theme(text = element_text(size=6), axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
plot(ggLysDensTauGut)

ggsave("lysogeny_in_communities_gut_B_and_tau.eps", device ="eps", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_gut_B_and_tau.pdf", device ="pdf", width = 8, height = 8, units = "cm")
ggsave("lysogeny_in_communities_gut_B_and_tau.png", device ="png", width = 8, height = 8, units = "cm")


df = data.frame("tau"=10^logtau.vec.gut.rd)

########
# OUTPUT data
########

### Output community data
df = data.frame("Bacteria per ml"=df.gut.com$M,"Viruses per ml"=df.gut.com$V,"adsorption.rate in ml per h"=10^logd.vec.gut.rd,"commitment time in h"=10^logtau.vec.gut.rd,"percentage of lysogens"=df.gut.com$PercL)
write.csv(df,"sampling_communities_gut.csv")

df = data.frame("Bacteria per ml"=df.marine.com$M,"Viruses per ml"=df.marine.com$V,"adsorption.rate in ml per h"=10^logd.vec.marine.rd,"commitment time in h"=10^logtau.vec.marine.rd,"percentage of lysogens"=df.marine.com$PercL)
write.csv(df,"sampling_communities_marine.csv")

### Output rank data
write.csv(COI.LHS.gut.com,"sampling_ranks_gut_COI.csv")
write.csv(COI.LHS.marine.com,"sampling_ranks_marine_COI.csv")
write.csv(PL.LHS.gut.com,"sampling_ranks_gut_ProbabLysog.csv")
write.csv(PL.LHS.marine.com,"sampling_ranks_marine_ProbabLysog.csv")
write.csv(Lys.LHS.gut.com,"sampling_ranks_gut_Lysogens.csv")
write.csv(Lys.LHS.marine.com,"sampling_ranks_marine_Lysogens.csv")

### Percentage of communities with different ranges of lysogeny
sink("percentage_lysogens_breaks.txt")
df.breaks.lys
sink()

## Contribution of top ranks to lysogeny
sink("contribution_lysogeney_ranks.txt")
df.ranks.lys
sink()

## Summary statitics of main properties for communities above lysogenic threshold
sink("lysogeny_in_communities_above_threshold_summary_statistics.txt")
print("Bacterial concentrations (cells/ml), gut")
summary(10^logM.vec.gut.rd.lys)
print("Bacterial concentrations (cells/ml), marine")
summary(10^logM.vec.marine.rd.lys)
print("")
print("Phage concentrations (phages/ml), gut")
summary(10^logV.vec.gut.rd.lys)
print("Phage concentrations (phages/ml), marine")
summary(10^logV.vec.marine.rd.lys)
print("")
print("Adsorption rates (ml/hr), gut")
summary(10^logd.vec.gut.rd.lys)
print("Adsorption rates (ml/hr), marine")
summary(10^logd.vec.marine.rd.lys)
print("")
print("Commitment times (hr), gut")
summary(10^logtau.vec.gut.rd.lys)
print("Commitment times (hr), marine")
summary(10^logtau.vec.marine.rd.lys)
print("")
print("COI top ranks (1,2,3), gut")
summary(COI.LHS.gut.dom1.lys)
summary(COI.LHS.gut.dom2.lys)
summary(COI.LHS.gut.dom3.lys)
print("COI top ranks (1,2,3), marine")
summary(COI.LHS.marine.dom1.lys)
summary(COI.LHS.marine.dom2.lys)
summary(COI.LHS.marine.dom3.lys)
sink()

## Summary statitics of main properties for communities above lysogenic threshold
sink("lysogeny_in_communities_above_threshold_summary_statistics_gut_threshold_25_or_above.txt")
print("Bacterial concentrations (cells/ml), gut")
summary(10^logM.vec.gut.rd.lys3)
print("")
print("Phage concentrations (phages/ml), gut")
summary(10^logV.vec.gut.rd.lys3)
print("")
print("Adsorption rates (ml/hr), gut")
summary(10^logd.vec.gut.rd.lys3)
print("")
print("Commitment times (hr), gut")
summary(10^logtau.vec.gut.rd.lys3)
print("")
sink()

## Summary statitics of main properties for communities above lysogenic threshold
sink("lysogeny_in_communities_above_threshold_summary_statistics_marine_threshold_10_or_above.txt")
print("Bacterial concentrations (cells/ml), marine")
summary(10^logM.vec.marine.rd.lys2)
print("")
print("Phage concentrations (phages/ml), marine")
summary(10^logV.vec.marine.rd.lys2)
print("")
print("Adsorption rates (ml/hr), marine")
summary(10^logd.vec.marine.rd.lys2)
print("")
print("Commitment times (hr), marine")
summary(10^logtau.vec.marine.rd.lys2)
print("")
sink()
