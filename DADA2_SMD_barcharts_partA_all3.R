# make barcharts to compare relative abundance of bacteria in spiders at different grouping levels

setwd("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/DADA2_SMD/")
################### IMPORT #####################################################
data.otu.norm <- read.csv("SMD_DADA_Fraction.OTUtable.csv", row.names=1)
data.tax.norm <- read.csv("SMD_DADA_Fraction.TAXtable.csv", row.names=1)
data.meta.norm <- read.csv("SMD_DADA_Fraction.SAMPLEtable.csv", row.names=1)

################### NO LAB #####################################################

# take away lab data
data.meta.norm <- data.meta.norm[which(data.meta.norm$population != "lab"), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

################### SUM FOR SPIDER SPECIES #####################################

# transpose
data.otu.t <- as.data.frame(t(data.otu.norm))

library(plyr)
# add column with col/pop/species you want summed
spider <- 'species'
data.otu.t$spider <- data.meta.norm[, spider]
sum.me <- data.otu.t

#add together data from same spider level
require(data.table)
SUM.ME <- data.table(sum.me)
im.summed <- SUM.ME[, lapply(.SD, sum ), by = spider]

# set rownames to spider group level and order by bact abundance 
data.summed <- as.data.frame(im.summed)
rownames(data.summed) <- data.summed$spider
data.summed$spider <- NULL

################### COLLECT OTHER ##############################################

# order by bact abundance and sum other
data.summed <- data.summed[, order(colSums(data.summed))]
data.summed[, (ncol(data.summed)-9):ncol(data.summed)]
other <- rowSums(data.summed[, 1:(ncol(data.summed)-10)])
data.summed.other <- data.summed[, (ncol(data.summed)-9):ncol(data.summed)]
data.summed.other <- cbind(other, data.summed.other)

# insert rownames as column
data.summed.other$subject <- rownames(data.summed.other)

# insert genus as columnnames
colnames(data.summed.other) <- c("other", 
                                 paste(as.vector(data.tax.norm[colnames(data.summed.other[, 2:11]), ]$Genus), 
                                       as.vector(rownames(data.tax.norm[colnames(data.summed.other[, 2:11]), ]))), "subject")
# order for custom colors:
column.order <- c("other", "Bergeyella ASV_11", "NA_Rickettsiaceae ASV_13", "Brevibacterium ASV_15",  "Brevibacterium ASV_6", "Acaricomes ASV_8", "NA_Chlamydiales ASV_1", "Diplorickettsia ASV_5", "Borrelia ASV_2", "Mycoplasma ASV_14", "Mycoplasma ASV_3", "subject")
data.summed.other <- data.summed.other[column.order]
color.order <- c("#9E0142",	"#F46D43","#FDAE61","#FEE08B","#FFFFBF","#E6F598","#ABDDA4","#66C2A5", "#3288BD", "#796CB2","#5E4FA2")
################### PLOT AND EXPORT ############################################
# melt
library(reshape2)
chart.data <- melt(data.summed.other, id = 'subject') 
# change columnnames
colnames(chart.data) <- c('sample', 'ASV', 'count')

library(RColorBrewer)
library(ggplot2)
require(scales)

# plot 
ggplot(chart.data, aes(x = sample)) + geom_bar(aes(weight=count, fill = ASV), position = 'fill') + 
  scale_y_continuous(labels = percent) + # percentages on y-axis
  scale_fill_manual(values = color.order) + # choses color palette
  theme_bw() + # white background
  theme(axis.title.x = element_blank()) +  # removing x-axis label
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = 'black', face = 'italic')) + # changes to x-axis text
  ggtitle("Microbiome of Social Stegodyphus") +
  ylab("Relative Abundance")  # adding y-axis label

ggsave("PartA_3speciesBarchartColorSort.png", 
       scale = 1, width = 17.5, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE)
