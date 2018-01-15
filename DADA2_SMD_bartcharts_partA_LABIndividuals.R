# make barcharts to compare relative abundance of bacteria in spiders at different grouping levels

setwd("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/DADA2_SMD/")
################### IMPORT #####################################################
data.otu.norm <- read.csv("SMD_DADA_Fraction.OTUtable.csv", row.names=1)
data.tax.norm <- read.csv("SMD_DADA_Fraction.TAXtable.csv", row.names=1)
data.meta.norm <- read.csv("SMD_DADA_Fraction.SAMPLEtable.csv", row.names=1)

################### NO LAB + JUST DUMICOLA #####################################

# take away lab data
data.meta.norm <- data.meta.norm[which(data.meta.norm$population == "lab"), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

# # only Dumicola
# species <- "dumicola"
# data.meta.norm <- data.meta.norm[which(data.meta.norm$species == species), ]
# data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

################### COLLECT OTHER ##############################################
data.summed <- as.data.frame(t(data.otu.norm))

# order by bact abundance and sum other
data.summed <- data.summed[, order(colSums(data.summed))]
data.summed[, (ncol(data.summed)-9):ncol(data.summed)]
other <- rowSums(data.summed[, 1:(ncol(data.summed)-10)])
data.summed.other <- data.summed[, (ncol(data.summed)-9):ncol(data.summed)]
data.summed.other <- cbind(other, data.summed.other)

# insert genus as columnnames
colnames(data.summed.other) <- c("other", 
                                 paste(as.vector(data.tax.norm[colnames(data.summed.other[, 2:11]), ]$Genus), 
                                       as.vector(rownames(data.tax.norm[colnames(data.summed.other[, 2:11]), ]))))
# insert nest names in rownames
rownames(data.summed.other) == rownames(data.meta.norm)
rownames(data.summed.other) <- paste(data.meta.norm$colony, rownames(data.summed.other))

# insert rownames as column
data.summed.other$subject <- rownames(data.summed.other)

# manualorder for colors
column.order <- c("other", "Bergeyella ASV_9", "Acaricomes ASV_248", "NA_Chlamydiales ASV_125", "NA_Chlamydiales ASV_1", "Diplorickettsia ASV_81", "Diplorickettsia ASV_10", "Borrelia ASV_58", "Borrelia ASV_2","Mycoplasma ASV_147", "Mycoplasma ASV_25", "subject")
color.order <- c("#9E0142","#F68562","#EAF7A9","#C7E8C2","#ABDDA4","#B3E1D2","#99D6C3","#99C4DE","#3288BD","#AFA7D1","#948AC1")
data.summed.other <- data.summed.other[column.order]
################### PLOT AND EXPORT ############################################
# melt
library(reshape2)
chart.data <- melt(data.summed.other, id = 'subject') 

# insert nest name for plot split
nests <- strsplit(as.character(chart.data$subject),' ') 
nests <- as.data.frame(do.call(rbind, nests))
chart.data$nest <- nests$V1
# remove nest names from sample name
chart.data$subject <- nests$V2
# change columnnames
colnames(chart.data) <- c('sample', 'ASV', 'count', 'nest')

library(RColorBrewer)
library(ggplot2)
require(scales)

chart.data$nest <- gsub("lab_PON", "mim Lab 1", chart.data$nest)
chart.data$nest <- gsub("lab_SKRU", "mim Lab 2", chart.data$nest)
chart.data$nest <- gsub("lab_KZN16", "dum Lab 1", chart.data$nest)
chart.data$nest <- gsub("lab_KZN11", "dum Lab 2", chart.data$nest)
ggplot(chart.data, aes(x = sample)) + geom_bar(aes(weight=count, fill = ASV), position = 'fill') +
  scale_y_continuous(labels = percent) + # percentages on y-axis
  scale_fill_manual(values = color.order) + # choses color palette
  theme_bw() + # white background
  theme(axis.title.x = element_blank()) +  # removing x-axis label
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = 'black')) + # changes to x-axis text
  facet_grid(.~nest, scales="free_x", space="free_x") + # splits
  theme(strip.text = element_text(size = 7)) +
  xlab("Individual spiders") +
  ylab("Relative Abundance")  # adding y-axis label

ggsave("PartA_LABIndividualsBarchart_colorOrederd.png",
       scale = 1, width = 15, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE)

