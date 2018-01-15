# gen_dist vs beta_dist

setwd("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/DADA2_SMD/")

################### DATA IMPORT ################################################
dum.names <- read.csv("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/rrSpidersFeb17analysis/3rdtry/dum.nameconversion.csv", header=FALSE, sep=";")
rownames(dum.names) <- dum.names$V1
#mim.names <- read.csv("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/rrSpidersFeb17analysis/3rdtry/mim.nameconversion.csv", header=FALSE, sep=";")
#rownames(mim.names) <- mim.names$V1

dum.gen.dist <- read.csv("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/rrSpidersFeb17analysis/3rdtry/dumicola_gendist.csv", sep=";")
dum.gen.dist <- dum.gen.dist[1:1128, ]
#mim.gen.dist <- read.csv("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/rrSpidersFeb17analysis/3rdtry/mimosarum_gendist.csv", sep=";")

# # I dont want to rewrite, so just do this
dum.gen.dist <- mim.gen.dist
dum.names <- mim.names

# species can be changed here to run scripts for different species.
# dont' have gen data for sar :(
#species <- "mimosarum"
species <- "dumicola"

# also change color for plot
# species colors:
#species_color = '#7570b3' # mim
species_color = '#d95f02' # dum

################### GET REAL NAMES FOR GEN DIST DATA ###########################
ind.name.1 <- vector()
for (ind in dum.gen.dist$ind_1) {
  if (ind %in% rownames(dum.names)) {
    hej <- as.vector(dum.names[ind, "V2"])
  }
  else {
    hej <- NA
  }
  
  ind.name.1 <- c(ind.name.1, hej)
}

ind.name.2 <- vector()
for (ind in dum.gen.dist$ind_2) {
  if (ind %in% rownames(dum.names)) {
    hej <- as.vector(dum.names[ind, "V2"])
  }
    else {
      hej <- NA
  }
  
  ind.name.2 <- c(ind.name.2, hej)
}

dum.gen.dist$pair <- paste(ind.name.1, ind.name.2) 

################### IMPORT #####################################################

# this should be normalized data in some way
data.otu.norm <- read.csv("SMD_DADA_Subsampl3527.OTUtable.csv", row.names=1)
data.tax.norm <- read.csv("SMD_DADA_Subsampl3527.TAXtable.csv", row.names=1)
data.meta.norm <- read.csv("SMD_DADA_Subsampl3527.SAMPLEtable.csv", row.names=1)

# take away lab data
data.meta.norm <- data.meta.norm[which(data.meta.norm$population != "lab"), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

################### NO LAB + JUST DUMICOLA #####################################

# take away lab data
data.meta.norm <- data.meta.norm[which(data.meta.norm$population != "lab"), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

data.meta.norm <- data.meta.norm[which(data.meta.norm$species== species), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

rownames(data.meta.norm) == colnames(data.otu.norm)

##################### SUM FOR COLONIES #########################################

# transpose
data.otu.norm.t <- as.data.frame(t(data.otu.norm))

# add column with col/pop/species you want summed
spider <- 'colony'
data.otu.norm.t$spider <- data.meta.norm[, spider]
sum.me <- data.otu.norm.t

#add together data from same spider level
require(data.table)
SUM.ME <- data.table(sum.me)
im.summed <- SUM.ME[, lapply(.SD, sum ), by = spider]

# set rownames to spider group level and order by bact abundance 
data.summed <- as.data.frame(im.summed)
rownames(data.summed) <- data.summed$spider
data.summed$spider <- NULL

# not just summed, but average of sum for each colony
hej <- table(droplevels(data.otu.norm.t$spider))
data.summed <- data.summed/hej

# reset variable
data.otu.colony <- as.data.frame(data.summed)

##################### GET B-DIVERSITY ##########################################
# If we denote the number of species shared between two sites as a and the 
# numbers of unique species (not shared) as b and c, then S = a + b + c and 
# ?? = (2 a + b + c)/2 so that ??_w = (b+c)/(2 a + b + c). 
# This is the Sorensen dissimilarity 
# as defined in vegan function vegdist with argument binary = TRUE
library(vegan)

method <- 'Sorensen'

col.B.dist <- list()
col1.list <- list()
col2.list <- list()
for (colony in rownames(data.otu.colony)) {
  for (colony2 in rownames(data.otu.colony)) {
    B <- vegdist(rbind(data.otu.colony[colony, ],data.otu.colony[colony2,]), method="bray", binary=TRUE)
    col.B.dist <- c(col.B.dist, B)
    col1.list <- c(col1.list, colony)
    col2.list <- c(col2.list, colony2)
  }
}
col.B.dist <- unlist(col.B.dist)
col2.list <- unlist(col2.list)
col1.list <- unlist(col1.list)

listed.beta.dist <- as.data.frame(cbind(col1.list, col2.list, col.B.dist))
colnames(listed.beta.dist) <- c("Var1", "Var2", "value")

######################### COMPARISON ##########################################

# set pairs as rownames 
listed.beta.dist$pair <- paste(listed.beta.dist$Var1, listed.beta.dist$Var2)
rownames(listed.beta.dist) <- paste(listed.beta.dist$Var1, listed.beta.dist$Var2)
rownames(dum.gen.dist) <- dum.gen.dist$pair

# only keep common pairs + sort in same order
beta.dist.common <- listed.beta.dist[rownames(listed.beta.dist) %in% rownames(dum.gen.dist), ]
gen.dist.common <- dum.gen.dist[rownames(dum.gen.dist) %in% rownames(beta.dist.common), ]
gen.dist.common.sort <- gen.dist.common[match(rownames(beta.dist.common), rownames(gen.dist.common)), ]

# test it
rownames(gen.dist.common.sort) == rownames(beta.dist.common)

# make one df and export for later use
dum.gendist.betadist <- as.data.frame(cbind(gen.dist.common.sort$pair, as.numeric(as.character(gen.dist.common.sort$gen_distance_)), as.numeric(as.character(beta.dist.common$value))))
colnames(dum.gendist.betadist) <- c("pair", "gendist", "betadiv")
write.csv(dum.gendist.betadist, sprintf("%s_gendist_betadist.csv", species))

# import
chart.data <- read.csv(sprintf("%s_gendist_betadist.csv", species), row.names = 1)

######################### PLOT AND SAVE ########################################

# ggplot version
library(ggplot2)
x <- as.vector(chart.data$gendist)
y <- as.vector(chart.data$betadiv)
plot.data <- data.frame(geo = x, beta = as.numeric(y))
ggplot(plot.data, aes(x=geo, y=beta)) +
  geom_point(shape = 1, color = sprintf("%s", species_color)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid lines
  ggtitle(sprintf("Beta diversity vs Genetic distance, S.%s nests", species)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Genetic distance (p-distance)") +
  ylab("Sorensen's dissimilarity") +
  theme(axis.text = element_text(size = 12)) + 
  geom_smooth(method = lm, se = FALSE, color = 1, lty = "dashed", size=0.5) + # add linear regression line, don't add shaded confidence region
  scale_y_continuous(limits = c(0.3, 1.0))

ggsave(sprintf("DADA2_SMD_PARTA_betaDiv_vs_GenDist_%s.png", species), 
       scale = 1, width = 15, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE)

#  cor test. null hypothesis: slope=0
cor.test(x,as.numeric(y))
