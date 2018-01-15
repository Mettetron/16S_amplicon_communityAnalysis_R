# using DADA2 output which has been through phyloseq

setwd("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/DADA2_SMD/")

# import raw data
data.otu.raw <- read.csv2("SMD_DADA2_OTUtable.csv", row.names=1)
data.tax.raw <- read.csv2("SMD_DADA2_TAXtable.csv", row.names=1)
data.meta.raw <- read.csv2("SMD_DADA2_SAMPLEtable.csv", row.names=1)

# checking that everything matches up
rownames(data.meta.raw) == rownames(data.otu.raw)
rownames(data.tax.raw) == colnames(data.otu.raw)

# filter away non-bact
data.tax.bact <- data.tax.raw[data.tax.raw$Kingdom == "Bacteria", ]
data.otu.bact <- data.otu.raw[, colnames(data.otu.raw) %in% rownames(data.tax.bact)]

# reset variables
data.tax <- data.tax.bact
data.otu <- data.otu.bact
data.meta <- data.meta.raw
data.meta$sampleID <- rownames(data.meta)

################### EXPORT "RAW" ###############################################
# export non-normaized data - only bact
write.csv(data.tax, "SMD_DADA2_NoBact_TAXtable.csv")
write.csv(data.otu, "SMD_DADA2_NoBact_OTUtable.csv")
write.csv(data.meta, "SMD_DADA2_NoBact_METAtable.csv")

################### SUBSAMPLING ################################################

# method A - subsampling
# this method randomly samples x reads from each sampleID with x+ reads

library(vegan)

# set a minimun number of reads e.g. 1st Quartile of the otu data summed per sample
min_reads <- summary(rowSums(data.otu))[2]
# remove samples-ids below min_reads
data.otu.norm <- data.otu[rowSums(data.otu)>min_reads, ]
# randomly sample the same number of reads from all the remaining sample_ids
set.seed(42)
data.otu.norm <- as.data.frame(rrarefy(data.otu.norm, min_reads))
# remove otus which are no longer represented in the data
data.otu.norm <- data.otu.norm[, colSums(data.otu.norm)>0]

# make table of sample IDs that were lost in subsampling
lost.sampleIDs <- data.meta[!(data.meta$sampleID %in% rownames(data.otu.norm)), ]
lost.sampleIDs$reads <- rowSums(data.otu[rownames(lost.sampleIDs), ])

# remove samples from meta data and tax data
data.meta.norm <- data.meta[(rownames(data.meta) %in% rownames(data.otu.norm)), ]
data.tax.norm <- data.tax[rownames(data.tax) %in% colnames(data.otu.norm), ]

# make tax-otu overview table
data.otu.norm.t <- as.data.frame(t(data.otu.norm))
data.tax.otu.norm <- cbind(data.tax.norm, data.otu.norm.t)
sample_nest_names <- paste(data.meta.norm$colony, data.meta.norm$sampleID)
colnames(data.tax.otu.norm)[7:ncol(data.tax.otu.norm)] <- sample_nest_names

# export
write.csv(data.otu.norm.t, file=sprintf("SMD_DADA_Subsampl%s.OTUtable.csv", min_reads))
write.csv(data.tax.norm, file=sprintf("SMD_DADA_Subsampl%s.TAXtable.csv", min_reads))
write.csv(data.meta.norm, file=sprintf("SMD_DADA_Subsampl%s.SAMPLEtable.csv", min_reads))
write.csv2(data.tax.otu.norm, file=sprintf("SMD_DADA_Subsampl%s.TAXandOTUtable.csv", min_reads))
################### FRACTION OF TOTAL READS #####################################

# set a minimun number of reads e.g. 1st Quartile of the otu data summed per sample
min_reads <- summary(rowSums(data.otu))[2]
# transform otu table to fit with old code
data.otu.t <- as.data.frame(t(data.otu))
# discard samples with less than you minimum of reads reads
data.otu.t.minreads <- data.otu.t[, which(colSums(data.otu.t) >= min_reads)]
data.meta.norm <- data.meta[which(rownames(data.meta) %in% colnames(data.otu.t.minreads)), ]

# fractionify

data.otu.t.norm <- data.otu.t.minreads
for (i in 1:ncol(data.otu.t.norm)){
  k <- sum(data.otu.t.norm[, i])
  for (j in 1:nrow(data.otu.t.norm)){
    data.otu.t.norm[j,i] <- data.otu.t.norm[j,i]/k
  }
}

# just checking that everything sums to 1
colSums(data.otu.t.norm)


# make tax-otu overview table
rownames(data.tax) == rownames(data.otu.t.norm)
data.tax.otu.norm <- cbind(data.tax, round(data.otu.t.norm,3))
sample_nest_names <- paste(data.meta.norm$colony, data.meta.norm$sampleID)
colnames(data.tax.otu.norm)[7:ncol(data.tax.otu.norm)] <- sample_nest_names

# export
write.csv(data.otu.t.norm, "SMD_DADA_Fraction.OTUtable.csv")
write.csv(data.tax, "SMD_DADA_Fraction.TAXtable.csv")
write.csv(data.meta.norm, "SMD_DADA_Fraction.SAMPLEtable.csv")
write.csv(data.tax.otu.norm, "SMD_DADA_Fraction.TAXandOTUtable.csv")

