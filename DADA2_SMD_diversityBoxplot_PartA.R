# making diversity boxplot comparing pops and nests

setwd("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/DADA2_SMD/")

################### IMPORT #####################################################

# OTU(or ASV) data has samples as columns, and is subsampled to the 1st quartile
data.otu.norm <- read.csv("SMD_DADA_Subsampl3527.OTUtable.csv", row.names=1)
# Taxonomy data has OTUs (or ASVs) as rows, and column names names = "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus" 
data.tax.norm <- read.csv("SMD_DADA_Subsampl3527.TAXtable.csv", row.names=1)
# meta data (or sample data) has samples as rows, and information about samples in columns, including lat and lon
data.meta.norm <- read.csv("SMD_DADA_Subsampl3527.SAMPLEtable.csv", row.names=1)

################### NO LAB + SINGLE SPECIES ####################################

# take away lab data
data.meta.norm <- data.meta.norm[which(data.meta.norm$population != "lab"), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

# species can be changed here to run script for different species.
#species <- "mimosarum"
species <- "dumicola"
#species <- "sarasinorum"
data.meta.norm <- data.meta.norm[which(data.meta.norm$species == species), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

# also change color for plot
# species colors:
#species_color = '#7570b3' # mim
species_color = '#d95f02' # dum
#species_color = '#1b9e77' # sar

# check that things line up
rownames(data.meta.norm) == colnames(data.otu.norm)

################### GET DIVERSITIES ############################################
# transpose
data.otu.norm.t <- as.data.frame(t(data.otu.norm))

# get Sorensen diversity
library(vegan)
beta.dist <- vegdist(data.otu.norm.t, method="bray", binary=TRUE)
# melt
library(reshape2)
beta.dist.df <- melt(as.matrix(beta.dist), varnames = c("row", "col"))

# remove self comparison
beta.dist.df <- beta.dist.df[beta.dist.df$row != beta.dist.df$col, ]

# remove double comparison
hej <- beta.dist.df[, 1:2]
hej.sort <- t(apply(hej, 1, sort))
beta.dist.df.unique <- beta.dist.df[!duplicated(hej.sort), ]

################# SORT BY COMPARISON TYPE ######################################

# insert population
sample.to.pop <- as.data.frame(cbind(rownames(data.meta.norm), as.character(data.meta.norm$population)))
colnames(sample.to.pop) <- c("sample", "population")
# row pops
new.frame <- as.data.frame(beta.dist.df.unique[, "row"])
new.frame[] <- lapply(new.frame, function(x) sample.to.pop$population[match(x, sample.to.pop$sample)])
colnames(new.frame) <- "pop1"
beta.dist.df.unique$pop1 <- new.frame$pop1
# col pops
new.frame <- as.data.frame(beta.dist.df.unique[, "col"])
new.frame[] <- lapply(new.frame, function(x) sample.to.pop$population[match(x, sample.to.pop$sample)])
colnames(new.frame) <- "pop2"
beta.dist.df.unique$pop2 <- new.frame$pop2

# insert nest
sample.to.nest <- as.data.frame(cbind(rownames(data.meta.norm), as.character(data.meta.norm$colony)))
colnames(sample.to.nest) <- c("sample", "nest")
# row pops
new.frame <- as.data.frame(beta.dist.df.unique[, "row"])
new.frame[] <- lapply(new.frame, function(x) sample.to.nest$nest[match(x, sample.to.nest$sample)])
colnames(new.frame) <- "nest1"
beta.dist.df.unique$nest1 <- new.frame$nest1
# row pops
new.frame <- as.data.frame(beta.dist.df.unique[, "col"])
new.frame[] <- lapply(new.frame, function(x) sample.to.nest$nest[match(x, sample.to.nest$sample)])
colnames(new.frame) <- "nest2"
beta.dist.df.unique$nest2 <- new.frame$nest2

# renaming
bddu <- beta.dist.df.unique

# different populations
bd.diff.pop <- bddu[bddu$pop1 != bddu$pop2, ]
# same population
bd.same.pop <- bddu[bddu$pop1 == bddu$pop2, ]
# same nest
bd.same.nest <- bddu[bddu$nest1 == bddu$nest2, ]

################# PLOT #########################################################

par(font.main = 3)
#jpeg(file=sprintf("Bdiv_Boxplot_%s.jpg", species), width = 17, height = 10, units="cm", res = 300) # this prints
boxplot(bd.same.nest$value, bd.same.pop$value, bd.diff.pop$value,
        border = species_color,
        ylab = "Sorensen's dissimilarity", 
        ylim = c(0, 1),
        main = sprintf("S.%s beta diversity", species), 
        names = c('within nests', 'within populations', 'between populations'))
#dev.off()


################# ANOVA TO COMPARE GROUPS ######################################
# method from: http://www.sthda.com/english/wiki/one-way-anova-test-in-r
same_nest <- cbind(as.numeric(as.character(bd.same.nest$value)), rep("s_nest", nrow(bd.same.nest)))
same_pop <- cbind(as.numeric(as.character(bd.same.pop$value)), rep("s_pop", nrow(bd.same.pop)))
diff_pop <- cbind(as.numeric(as.character(bd.diff.pop$value)), rep("d_pop", nrow(bd.diff.pop)))
anova.data <- as.data.frame(rbind(same_nest, same_pop, diff_pop))   
colnames(anova.data) <- c("Bdiv", "group")
anova.data$Bdiv <- as.numeric(as.character(anova.data$Bdiv))

levels(anova.data$group)
anova.data$group <- ordered(anova.data$group, levels = c("d_pop", "s_nest", "s_pop" ))
library(dplyr)
group_by(anova.data, group) %>%
  summarise(
    count = n(),
    mean = mean(Bdiv, na.rm = TRUE),
    sd = sd(Bdiv, na.rm = TRUE)
  )

bdiv.anova <- aov(Bdiv ~ group, data = anova.data)
summary(bdiv.anova)
TukeyHSD(bdiv.anova)

pairwise.t.test(anova.data$Bdiv, anova.data$group,
                p.adjust.method = "BH")
#### sar
# Df Sum Sq Mean Sq F value Pr(>F)    
# group          2  41.29  20.645    1148 <2e-16 ***
#   Residuals   4789  86.12   0.018                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > TukeyHSD(bdiv.anova)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Bdiv ~ group, data = anova.data)
# 
# $group
# diff        lwr        upr p adj
# s_nest-d_pop -0.4101678 -0.4375954 -0.3827402     0
# s_pop-d_pop  -0.1851372 -0.1974623 -0.1728122     0
# s_pop-s_nest  0.2250305  0.1958217  0.2542393     0
# 
# > pairwise.t.test(anova.data$Bdiv, anova.data$group,
#                   +                 p.adjust.method = "BH")
# 
# Pairwise comparisons using t tests with pooled SD 
# 
# data:  anova.data$Bdiv and anova.data$group 
# 
# d_pop  s_nest
# s_nest <2e-16 -     
#   s_pop  <2e-16 <2e-16
# 
# P value adjustment method: BH 

### mim
# > summary(bdiv.anova)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# group          2  0.556 0.27804   24.87 2.27e-11 ***
#   Residuals   1698 18.985 0.01118                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > TukeyHSD(bdiv.anova)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Bdiv ~ group, data = anova.data)
# 
# $group
# diff         lwr         upr     p adj
# s_nest-d_pop -0.08865232 -0.12509170 -0.05221294 0.0000000
# s_pop-d_pop  -0.03029242 -0.04575132 -0.01483352 0.0000137
# s_pop-s_nest  0.05835990  0.01995978  0.09676001 0.0010883
# 
# > pairwise.t.test(anova.data$Bdiv, anova.data$group,
#                   +                 p.adjust.method = "BH")
# 
# Pairwise comparisons using t tests with pooled SD 
# 
# data:  anova.data$Bdiv and anova.data$group 
# 
# d_pop   s_nest 
# s_nest 4.1e-08 -      
#   s_pop  6.9e-06 0.00037
# 
# P value adjustment method: BH 

### dum
# > summary(bdiv.anova)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# group          2   1.89  0.9457   29.71 2.04e-13 ***
#   Residuals   1748  55.63  0.0318                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > TukeyHSD(bdiv.anova)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = Bdiv ~ group, data = anova.data)
# 
# $group
# diff         lwr         upr     p adj
# s_nest-d_pop -0.20374626 -0.27086818 -0.13662433 0.0000000
# s_pop-d_pop  -0.03762379 -0.06297604 -0.01227155 0.0014838
# s_pop-s_nest  0.16612247  0.09617500  0.23606993 0.0000001
# 
# > pairwise.t.test(anova.data$Bdiv, anova.data$group,
#                   +                 p.adjust.method = "BH")
# 
# Pairwise comparisons using t tests with pooled SD 
# 
# data:  anova.data$Bdiv and anova.data$group 
# 
# d_pop   s_nest 
# s_nest 4.7e-12 -      
#   s_pop  0.00051 4.4e-08
# 
# P value adjustment method: BH 