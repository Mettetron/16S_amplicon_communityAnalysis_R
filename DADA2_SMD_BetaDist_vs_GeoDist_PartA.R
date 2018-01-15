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

# species can be changed here to run scripts for different species.
#species <- "mimosarum"
species <- "dumicola"
#species <- "sarasinorum"
data.meta.norm <- data.meta.norm[which(data.meta.norm$species == species), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

# also change color for plot
# species colors:
#'#1b9e77', '#7570b3', '#d95f02'
#"S.sarasinorum","S.mimosarum", "S.dumicola"
species_color = '#d95f02'

# check that things line up
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

################### GET GEO DISTANCES ##########################################
library(sp)
# take only col, lat and lon
col.coords <- data.meta.norm[, c('colony', 'lat','lon')]
# get unique cases - so one lien for each colony
col.coords <- unique(col.coords)
#set colony as rowname
rownames(col.coords) <- col.coords$colony
# remove colony column
col.coords$colony <- NULL

# make the geo distance matrix. dist in km
library(sp)
col.geo.dist <- apply(col.coords, 1, function(eachPoint) spDistsN1(as.matrix(col.coords), eachPoint, longlat=TRUE))
rownames(col.geo.dist) <- colnames(col.geo.dist)

listed.geo.dist <- melt(col.geo.dist)

##################### GET B-DIVERSITY ##########################################
# If we denote the number of species shared between two sites as a and the 
# numbers of unique species (not shared) as b and c, then S = a + b + c and 
# ?? = (2 a + b + c)/2 so that ??_w = (b+c)/(2 a + b + c). 
# This is the Sorensen dissimilarity 
# as defined in vegan function vegdist with arguments method="bray" and  binary=TRUE
library(vegan)

method <- 'bray_pa'
#hm. lets do it the hard way
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

##################### PLOT #####################################################

# make table with both variables
listed.beta.dist$Var2 == listed.geo.dist$Var1 # checking that things line up
listed.beta.dist$Var1 == listed.geo.dist$Var2 # checking that things line up
listed.beta.geo.dist <- listed.beta.dist
listed.beta.geo.dist$geo <- listed.geo.dist$value
colnames(listed.beta.geo.dist) <- c("colony1", "colony2", "beta", "geo")

# remove self comparisons
listed.beta.geo.dist <- listed.beta.geo.dist[listed.beta.geo.dist$beta != 0, ]

# remove double comparisons
rownames(listed.beta.geo.dist) <- seq(1:nrow(listed.beta.geo.dist))
for (i in rownames(listed.beta.geo.dist)){
  hej <- paste(listed.beta.geo.dist[i, ]$colony1, listed.beta.geo.dist[i, ]$colony2) 
  for (j in rownames(listed.beta.geo.dist)) {
    hej2 <- paste(listed.beta.geo.dist[j, ]$colony2, listed.beta.geo.dist[j, ]$colony1) 
    if (hej == hej2) {
      if (!is.na(listed.beta.geo.dist[i, "beta"])) {
        listed.beta.geo.dist[j, "beta"] <- NA
      }
    }
  }
}


listed.beta.geo.dist.fixed <- listed.beta.geo.dist[complete.cases(listed.beta.geo.dist), ]

listed.beta.geo.dist <- listed.beta.geo.dist.fixed

# export for future plot tune-up etc
#write.csv(listed.beta.geo.dist, sprintf("%s_listed_beta_geo_dist.csv", species))


# species can be changed here to run script for different species.
species <- "mimosarum"
#species <- "dumicola"
#species <- "sarasinorum"

# also change color for plot
# species colors:
species_color = '#7570b3' # mim
#species_color = '#d95f02' # dum
#species_color = '#1b9e77' # sar

# import ready made data
listed.beta.geo.dist <- read.csv(sprintf("%s_listed_beta_geo_dist.csv", species), row.names = 1)


# # plot
# x <- as.vector(listed.beta.geo.dist$geo)
# y <- as.vector(listed.beta.geo.dist$beta)
# plot(x, y, col = '#d95f02', main = sprintf("%s", species),  xlab = 'Geographical distance (km)', ylab ="Beta diversity (Sorensen's dissimilarity)")
# abline(lm(y ~ x), lty = 2:6)

library(ggplot2)
# ggplot version
x <- as.vector(listed.beta.geo.dist$geo)
y <- as.vector(listed.beta.geo.dist$beta)
plot.data <- data.frame(geo = x, beta = as.numeric(y))
ggplot(plot.data, aes(x=geo, y=beta)) +
  geom_point(shape = 1, color = sprintf("%s", species_color)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid lines
  ggtitle(sprintf("Beta diversity vs Geo distance, S.%s nests", species)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Geographical distance (km)") +
  ylab("Sorensen's dissimilarity") +
  theme(axis.text = element_text(size = 12)) + 
  geom_smooth(method = lm, se = FALSE, color = 1, lty = "dashed", size=0.5) + # add linear regression line, don't add shaded confidence region
  scale_y_continuous(limits = c(0.3, 1.0))

ggsave(sprintf("DADA2_SMD_PARTA_betaDiv_vs_geoDist_%s.png", species), 
       scale = 1, width = 15, height = 8, units = "cm",
       dpi = 300, limitsize = TRUE)

#  cor test. null hypothesis: slope=0
cor.test(x,as.numeric(y))

##################### LOOK AT  INDIVIDUAL POPULATIONS ##########################
# species can be changed here to run scripts for different species.
# dont' have accurate geo data for sar :(
species <- "mimosarum"
#species <- "dumicola"

# also change color for plot
# species colors:
species_color = '#7570b3' # mim
#species_color = '#d95f02' # dum

# import ready made data
listed.beta.geo.dist <- read.csv(sprintf("%s_listed_beta_geo_dist.csv", species), row.names = 1)

pops <- strsplit(as.character(listed.beta.geo.dist$colony1),'_') 
pops <- as.data.frame(do.call(rbind, pops))
listed.beta.geo.dist$pop1 <- pops$V1

pops <- strsplit(as.character(listed.beta.geo.dist$colony2),'_') 
pops <- as.data.frame(do.call(rbind, pops))
listed.beta.geo.dist$pop2 <- pops$V1

lbgd <- listed.beta.geo.dist

for (pop in unique(lbgd$pop2)){
  population = pop
  lbgd.pop <- lbgd[lbgd$pop1 == population & lbgd$pop2 == population, ]
  
  x <- as.vector(lbgd.pop$geo)
  y <- as.vector(lbgd.pop$beta)
  
  # simple plot
  jpeg(file=sprintf("Bdiv_vs_Geodist_%s_%s.jpg", species, population), width = 529, height = 359, units="px") # this prints
  plot(x, y, col=species_color, xlab = 'Geographical distance (km)', ylab ="Sorensen's dissimilarity", 
       main=sprintf("%s", population),
       ylim=c(0.3, 1))
  abline(lm(y ~ x), lty = 2:6)
  hej <- cor.test(x,as.numeric(y))
  legend("bottomright", inset=0.02,
         paste("cor = ", round(hej$estimate,3),"\np-val = ", round(hej$p.value, 3), sep=""))
  dev.off()
}
