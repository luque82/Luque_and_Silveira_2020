## This R script analyzes the rank abundance for marine and gut viruses and microbes

### Working directory
setwd("/Volumes/GoogleDrive/My Drive/1_contributions/1_articles/1_in_progress/6_COI_and_lysogeny_low_cell_densities/wip/code_development/3_rank_abundance_code")

### Packages
library("ggplot2")
library("svglite")

# Color-blind-friendly palettes, one with gray, and one with black (palettes source: http://jfly.iam.u-tokyo.ac.jp/color/).
cbp1 = c("#999999", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Read MOI data
df.marine.bacteria = read.table(file = 'tara_mOTU.linkage-groups.relab.release.tsv', sep = '\t', header = TRUE)
df.gut.bacteria = read.csv(file="gut_bwa_counts-unique.csv")

####################
### Marine bacteria
####################
# Check normalization of columns
df.marine.bacteria.sum = colSums(df.marine.bacteria[,-1])
# Output file
sink('rank_abundance_analysis_marine_bacteria_sum.txt')
print('Dimension of the data frame rows x columns')
dim(df.marine.bacteria)
print('Sum of each column in the marine bacteria data set. Sorted output: min to max')
sort(df.marine.bacteria.sum)
sink()

# Rank abunance in each sample
rows = nrow(df.marine.bacteria)
cols = ncol(df.marine.bacteria) -1
df.marine.bacteria.ranked = data.frame(matrix("", ncol = cols, nrow = rows))
iiseq = seq(1,cols)
for(ii in iiseq){
  yranked = sort(df.marine.bacteria[[ii+1]], decreasing = TRUE)
  df.marine.bacteria.ranked[[ii]] = yranked
}

# Medians for the top 100 across samples
top.matrix = as.matrix(df.marine.bacteria.ranked[1:100,])
rows = nrow(top.matrix)
marine.bacteria.medians = matrix(0 ,ncol = 1, nrow = rows)
iiseq = seq(1,rows)
for(ii in iiseq){
  marine.bacteria.medians[ii] = median(top.matrix[ii,])  
}

# plot
df = as.data.frame(marine.bacteria.medians)
ggmarbr = ggplot(df) +
  geom_point(aes(x = as.numeric(rownames(df)), y = V1)) +
  scale_x_log10() + scale_y_log10() +
  theme_classic()
plot(ggmarbr)

####################
### Gut bacteria (not normlized)
####################

# Normalized rank abunance in each sample
rows = nrow(df.gut.bacteria)
cols = ncol(df.gut.bacteria) -1
df.gut.bacteria.ranked = data.frame(matrix("", ncol = cols, nrow = rows))
iiseq = seq(1,cols)
for(ii in iiseq){
  yranked = sort(df.gut.bacteria[[ii+1]], decreasing = TRUE)
  df.gut.bacteria.ranked[[ii]] = yranked/sum(yranked)
  print(sum(yranked))
}

# Medians for the top 100 across samples
top.matrix = as.matrix(df.gut.bacteria.ranked[1:100,])
xNA = which(is.na(top.matrix[1, ])) ## identify NA values
top.matrix.cure = top.matrix[, -xNA] ## curate NA values
rows = nrow(top.matrix.cure)
gut.bacteria.medians = matrix(0 ,ncol = 1, nrow = rows)
iiseq = seq(1,rows)
for(ii in iiseq){
  gut.bacteria.medians[ii] = median(top.matrix.cure[ii,])  
}

# plot
df = as.data.frame(gut.bacteria.medians)
gggutbr = ggplot(df) +
  geom_point(aes(x = as.numeric(rownames(df)), y = V1)) +
  scale_x_log10() + scale_y_log10() +
  theme_classic()
plot(gggutbr)

#### Combined marine and gut
df = data.frame("gut"=gut.bacteria.medians,"marine"=marine.bacteria.medians)
ggrankbr = ggplot(df) +
  geom_point(aes(x = as.numeric(rownames(df)), y = gut), color=cbp1[1]) +
  geom_point(aes(x = as.numeric(rownames(df)), y = marine), color=cbp1[2]) +
  scale_x_log10() + scale_y_log10() +
  theme_classic()
plot(ggrankbr)

### Output data
df.bacteria.ranks = data.frame("rank"=as.numeric(rownames(df)),"gut"=gut.bacteria.medians,"marine"=marine.bacteria.medians)
write.csv(df.bacteria.ranks,"rank_abundance_analysis_bacteria.csv")
