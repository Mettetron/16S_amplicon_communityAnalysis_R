setwd("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/DADA2_SMD/")

################### IMPORT #####################################################

# this should be normalized data in some way
data.otu.norm <- read.csv("SMD_DADA_Fraction.OTUtable.csv", row.names=1)
data.tax.norm <- read.csv("SMD_DADA_Fraction.TAXtable.csv", row.names=1)
data.meta.norm <- read.csv("SMD_DADA_Fraction.SAMPLEtable.csv", row.names=1)

# merge tax and otu
data.tax.otu.norm <- cbind(data.tax.norm, data.otu.norm)

################### NO LAB #####################################################

# take away lab data
data.meta.norm <- data.meta.norm[which(data.meta.norm$population != "lab"), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

# merge tax and otu
data.tax.otu.norm <- cbind(data.tax.norm, data.otu.norm)

################### FIND CORE ##################################################

# make presence/absence data
library(vegan)
data.otu.pa <- data.otu.norm
threshold <- 0.001  # for fraction data, this is 0.1%
data.otu.pa[data.otu.pa < threshold] <- 0 
data.otu.pa <- decostand(data.otu.pa, method="pa")

# now sum up all the rows 
otu.sums <- rowSums(data.otu.pa)

# just get all otus ordered by most core-like
data.tax.otu.corenum <- data.tax.otu.norm
data.tax.otu.corenum$fiveplus <- otu.sums
data.tax.otu.corenum.sorted <- data.tax.otu.corenum[order(-data.tax.otu.corenum$fiveplus), ]

# set minimum number to be called core:
# must be present in 20%
core.num <- ncol(data.otu.norm)/5

# only keep OTUs above core num
data.tax.otu.corenum.sorted.core <- data.tax.otu.corenum.sorted[data.tax.otu.corenum.sorted$fiveplus > core.num, ]

# get shorter name
data.core <- data.tax.otu.corenum.sorted.core

################### GET PERCENTAGES ############################################
# make population names unique for sepcies
species.letter <- substr(data.meta.norm$species, 1, 3)
data.meta.norm$population <- paste("S.", species.letter, " ", data.meta.norm$population, sep = "")
# make a list of unique population names
pop.list <- unique(data.meta.norm$population)
pop.list

# custom pop list for order
pop.list <- c("S.dum PAA","S.dum ADDO","S.dum WEE","S.dum PON", "S.dum KRU",
              "S.mim WEE", "S.mim PON", "S.mim SAK", "S.mim MAH", "S.mim TANA", 
              "S.sar SriLanka", "S.sar E", "S.sar D", "S.sar A", "S.sar B","S.sar C",
              "S.sar RAO", "S.sar Himalaya")

# # alternative pop list for lab inclusion
# pop.list <- c("S.dum PAA","S.dum ADDO","S.dum WEE","S.dum PON", "S.dum KRU", "S.dum lab",
#               "S.mim WEE", "S.mim PON", "S.mim SAK", "S.mim MAH", "S.mim TANA", "S.mim lab", 
#               "S.sar SriLanka", "S.sar E", "S.sar D", "S.sar A", "S.sar B","S.sar C",
#               "S.sar RAO", "S.sar Himalaya")

# separate otu and tax and put otu column in data.core.tax
data.core.otu <- data.core[, 7:(ncol(data.core)-1)] # getting rid of the fiveplus column
data.core.tax <- data.core[, 1:6]
data.core.tax$otu <- rownames(data.core.tax)

# test that otu and meta data matches
colnames(data.core.otu) == rownames(data.meta.norm)

# make p/a again
data.core.otu.t <- as.data.frame(t(data.core.otu))
#data.core.otu.t[data.core.otu.t < 5] <- 0 # threshold for presence = 5
data.core.otu.t[data.core.otu.t < threshold] <- 0 # threshold for presence = 0.1%
data.core.otu.t.pa <- decostand(data.core.otu.t, method="pa")

# initialize result data frame
percent.per.pop <- data.core.otu[, 1:2]

# loop over pops and take % of indv in pop or col containing the otu
for (pop in pop.list){
  pop.matrix <- data.core.otu.t.pa[data.meta.norm$population == pop, ] # change here if you're doing colony or population!
  pop.otus <- colSums(pop.matrix)/nrow(pop.matrix)*100
  percent.per.pop[, pop] <- pop.otus
}

# remove those first two columns
percent.per.pop[, 1] <- NULL
percent.per.pop[, 1] <- NULL

################### HEATMAP DATA PREP + EXPORT FOR TREE ########################
# insert genus names
rownames(percent.per.pop) <- paste(data.core.tax$Genus, data.core.tax$otu)


##########TREE STUFF VVV ###############################################
# 
# ### make fasta to export for making tree
# library(Biostrings)
# # import seqs
# tax.seqs <- read.csv2("SMD_DADA2_seqTAXtable.csv")
# tax.seqs$names <- paste(tax.seqs$Genus, tax.seqs$numbers)
# # get just the ones used in heatmap
# tax.seqs.map <- tax.seqs[tax.seqs$names %in% rownames(percent.per.pop), ]
# # make biostrings object and export as fasta
# seq <- as.character(tax.seqs.map$seq)
# hej <- as.character(tax.seqs.map$names)
# names(seq) <- hej
# dna <- DNAStringSet(seq)
# writeXStringSet(dna, "prevalenceHeatmapSeqs.fasta")
# 
# # import Tree
# library(ape)
# tree <- read.tree("HeatmapTree20percent.nwk")
# # deal with branch lengths of zero + line up straight
# tree$edge.length[which(tree$edge.length == 0)] <- 0.00001
# tree <- chronopl(tree, lambda = 0.1, tol = 0)
# tree <- as.dendrogram(as.hclust.phylo(tree))
# # plot
# plot(tree)
# 
# # make sure same genomes in both and get same order in heatmap rows as in tree
# clade_order <- order.dendrogram(tree)
# clade_name <- labels(tree)
# # insert _ in rownames of heatmap data
# rownames(percent.per.pop) <- gsub(' ', '_', rownames(percent.per.pop))
# 
# # check that names match
# rownames(percent.per.pop) %in% clade_name
# 
# # order
# clade_position <- data.frame(clade_name, clade_order)
# clade_position <- clade_position[order(clade_position$clade_order),]
# new_order <- match(clade_position$clade_name, row.names(percent.per.pop))
# percent.per.pop.ordered <- percent.per.pop[new_order,]
# # check
# rownames(percent.per.pop.ordered) == clade_name
# 
# # reset variable
# percent.per.pop <- percent.per.pop.ordered
# 
# # order columns
# #cluster.matrix.ordered <- cluster.matrix.ordered[, order(colSums(cluster.matrix.ordered))]
# 
# # make matrix
# percentpop.matrix <- as.matrix(percent.per.pop)
# 
# # make plot that gives the tree structure (cot this out and paste on next one)
# jpeg("TREEforPrevalence_heatmap_PartA_20percent.jpg", width = 15, height = 10, units="cm", res=300) # this prints
# heatmap(percentpop.matrix, Rowv = tree, Colv = NA)
# dev.off()
############ TREE STUFF ^^^ ####################################################

# make matrix
percentpop.matrix <- as.matrix(percent.per.pop)

#transpose
percentpop.matrix.t <- t(percentpop.matrix)

################### PLOT AND SAVE ##############################################

library(reshape2)
library(ggplot2)
dat <- percentpop.matrix.t
dat2 <- melt(dat)

# preparing for facets
dat2$species <-  dat2$Var1
species.names <- strsplit(as.character(dat2$species),' ') 
species.names <- as.data.frame(do.call(rbind, species.names))
dat2$species <- species.names$V1
dat2$species <- gsub('S.dum', 'S. dumicola', dat2$species)
dat2$species <- gsub('S.mim', 'S. mimosarum', dat2$species)
dat2$species <- gsub('S.sar', 'S. sarasinorum', dat2$species)

dat2$Var1 <- species.names$V2

## plot
#jpeg("Prevalence_heatmap_PartA_20percent.jpg", width = 15, height = 10, units="cm", res=300) # this prints
ggplot(dat2, aes(as.factor(Var1), Var2, group=Var2)) +
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme(axis.title.y = element_blank()) +
  ylab("Symbiont OTUs") + 
  xlab("Spider populations") +
  facet_grid(.~species, scales="free_x", space="free_x") +
  theme(strip.text=element_text(face = "italic")) + 
  theme(axis.title.x = element_text(size = 9)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = 'black')) +
  theme(axis.text.y = element_text(size = 8, color = 'black', face="italic"))
#dev.off()
  