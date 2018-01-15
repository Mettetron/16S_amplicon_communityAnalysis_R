# ordination of spider microbiome data

setwd("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/DADA2_SMD/")

library(vegan)

################### IMPORT #####################################################

# subsampl
data.otu <- read.csv("SMD_DADA_Subsampl3527.OTUtable.csv", row.names=1)
data.tax <- read.csv("SMD_DADA_Subsampl3527.TAXtable.csv", row.names=1)
data.meta <- read.csv("SMD_DADA_Subsampl3527.SAMPLEtable.csv", row.names=1)

data.norm <- "subsampl3527"
data.type <- "abundance"

# take away lab data
data.meta <- data.meta[which(data.meta$population != "lab"), ]
data.otu <- data.otu[, which(colnames(data.otu) %in% rownames(data.meta))]

# transpose otu data
data.otu.t <- as.data.frame(t(data.otu))
# make color vector
rownames(data.otu.t) == rownames(data.meta)
sars <- 2*(data.meta$species == "sarasinorum")
dums <- 3*(data.meta$species == "dumicola")
mims <- 4*(data.meta$species == "mimosarum")
col.vector <- sars+mims+dums
# replace with other colors
# species colors:
#'#1b9e77', '#7570b3', '#d95f02'
#"S.sarasinorum","S.mimosarum", "S.dumicola"
col.vector <- replace(col.vector, col.vector == 4, '#7570b3')
col.vector <- replace(col.vector, col.vector == 2, '#1b9e77')
col.vector <- replace(col.vector, col.vector == 3, '#d95f02')

### NMDS
nmds.otu.t <- metaMDS(data.otu.t, k=2, distance = "bray") # from vegan package
# stress > 20!

# make an empty plot
plot(nmds.otu.t, display = "sites", type = "n", main = "Stegodyphus Microbiome Ordination")
# and plot our NMDS species with color vector
orditorp(nmds.otu.t, display = "sites", labels= "n", col = col.vector, pch = 16) 
# try to fit OTUs
# colnames(data.otu.t) <- data.tax$Genus
# fit <- envfit(nmds.otu.t ~ Mycoplasma+NA_Chlamydiales+Borrelia+Rickettsia+Entomoplasma+Bergeyella+Brevibacterium+Diplorickettsia+Acaricomes, data=data.otu.t)
# plot(fit, p.max = 0.001, col=1, cex=0.8, font=4)
legend(1.5, 0, c(expression(italic("S.sarasinorum")),expression(italic("S.mimosarum")),expression(italic("S.dumicola"))),
       pch = 16, col = c('#1b9e77', '#7570b3', '#d95f02'), cex = 0.8)

