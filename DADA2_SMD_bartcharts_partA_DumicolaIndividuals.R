# make barcharts to compare relative abundance of bacteria in spiders at different grouping levels

setwd("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/DADA2_SMD/")
################### IMPORT #####################################################
data.otu.norm <- read.csv("SMD_DADA_Fraction.OTUtable.csv", row.names=1)
data.tax.norm <- read.csv("SMD_DADA_Fraction.TAXtable.csv", row.names=1)
data.meta.norm <- read.csv("SMD_DADA_Fraction.SAMPLEtable.csv", row.names=1)

################### NO LAB + JUST DUMICOLA #####################################

# take away lab data
data.meta.norm <- data.meta.norm[which(data.meta.norm$population != "lab"), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

# only Dumicola
species <- "dumicola"
data.meta.norm <- data.meta.norm[which(data.meta.norm$species == species), ]
data.otu.norm <- data.otu.norm[, which(colnames(data.otu.norm) %in% rownames(data.meta.norm))]

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

# manually fix order to normalize colors with the 3 species barchart
column.order <- c("other", "Bergeyella ASV_11", "Rickettsia ASV_4", "NA_Rickettsiaceae ASV_13", "Acaricomes ASV_8", "Diplorickettsia ASV_10", "Diplorickettsia ASV_5", "Borrelia ASV_41", "Borrelia ASV_2", "Mycoplasma ASV_14",  "Mycoplasma ASV_3", "subject")
color.order <- c("#9E0142","#F46D43","#FEC996","#FDAE61","#E6F598","#99D6C3","#66C2A5","#549CC8","#3288BD","#796CB2","#5E4FA2")

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

## plots used for partA - see below for other options

pops <- strsplit(as.character(chart.data$nest),'_')
pops <- as.data.frame(do.call(rbind, pops))
chart.data$pop <- pops$V1

# manually split into PAA + ADDO
chart.data.pop <- chart.data[chart.data$pop == 'PAA' | chart.data$pop == 'ADDO' , ]
ggplot(chart.data.pop, aes(x = sample)) + geom_bar(aes(weight=count, fill = ASV), position = 'fill') + 
  scale_y_continuous(labels = percent) + # percentages on y-axis
  scale_fill_manual(values = color.order ) + # choses color palette
  theme_bw() + # white background
  theme(axis.title.x = element_blank()) +  # removing x-axis label
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = 'black')) + # changes to x-axis text
  facet_grid(.~nest, scales="free_x", space="free_x") + # splits
  theme(strip.text = element_text(size=5)) + 
  theme(legend.text=element_text(size=6.8)) +
  theme(panel.spacing= unit(0.1, "cm")) +
  ylab("Relative Abundance")  # adding y-axis label
ggsave("PartA_DumicolaIndividualsBarchart_ADDO_PAA_colorOrdered.png",
       scale = 1, width = 20, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE)


# manually split into PON + KRU + WEE
chart.data.pop <- chart.data[chart.data$pop == 'PON' | chart.data$pop == 'KRU' | chart.data$pop == 'WEE' , ]
ggplot(chart.data.pop, aes(x=sample)) + geom_bar(aes(weight=count, fill=ASV), position='fill') + 
  scale_y_continuous(labels=percent) + # percentages on y-axis
  scale_fill_manual(values=color.order) + # choses color palette
  theme_bw() + # white background
  theme(axis.title.x=element_blank()) +  # removing x-axis label
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8, color='black')) + # changes to x-axis text
  facet_grid(.~nest, scales="free_x", space="free_x") + # splits
  theme(strip.text=element_text(size=5)) + 
  theme(panel.spacing=unit(0.1, "cm")) +
  theme(legend.position="none") + #removing legend
  ylab("Relative Abundance")  # adding y-axis label
ggsave("PartA_DumicolaIndividualsBarchart_PON_KRU_WEE_colorOrdered.png",
       scale = 1, width = 20, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE)

##### other options #######
# # plot 
# ggplot(chart.data, aes(x = sample)) + geom_bar(aes(weight=count, fill = ASV), position = 'fill') + 
#   scale_y_continuous(labels = percent) + # percentages on y-axis
#   scale_fill_manual(values = brewer.pal(11, "Spectral")) + # choses color palette
#   theme_bw() + # white background
#   theme(axis.title.x = element_blank()) +  # removing x-axis label
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = 'black')) + # changes to x-axis text
#   facet_grid(.~nest, scales="free_x", space="free_x") + # splits
#   theme(strip.text = element_text(size = 5)) + 
#   xlab("Individual spiders") +
#   ylab("Relative Abundance")  # adding y-axis label
# 
# ggsave("PartA_DumicolaIndividualsBarchart_sketch.png",
#        scale = 1, width = 20, height = 10, units = "cm",
#        dpi = 300, limitsize = TRUE)
# 
# # this is a nice and informative plot, but too big. 
# # split into one plot for each population
# # insert population name for plot split
# pops <- strsplit(as.character(chart.data$nest),'_') 
# pops <- as.data.frame(do.call(rbind, pops))
# chart.data$pop <- pops$V1
# 
# for (pop in unique(chart.data$pop)) {
#   chart.data.pop <- chart.data[chart.data$pop == pop, ]
#   nests <- strsplit(as.character(chart.data.pop$nest),'_') 
#   nests <- as.data.frame(do.call(rbind, nests))
#   chart.data.pop$nest <- nests$V2
#   ggplot(chart.data.pop, aes(x = sample)) + geom_bar(aes(weight=count, fill = ASV), position = 'fill') + 
#     scale_y_continuous(labels = percent) + # percentages on y-axis
#     scale_fill_manual(values = brewer.pal(11, "Spectral")) + # choses color palette
#     theme_bw() + # white background
#     theme(axis.title.x = element_blank()) +  # removing x-axis label
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = 'black')) + # changes to x-axis text
#     facet_grid(.~nest, scales="free_x", space="free_x") + # splits
#     theme(strip.text = element_text(size = 9)) + 
#     ggtitle(sprintf("S.%s, population %s", species, pop)) +
#     xlab("Individual spiders") +
#     theme(panel.spacing= unit(0.1, "cm")) +
#     theme(legend.position="none") + #removing legend
#     ylab("Relative Abundance")  # adding y-axis label
#   ggsave(sprintf("PartA_DumicolaIndividualsBarchart_%s.png", pop),
#          scale = 1, width = (1+nrow(chart.data.pop)/11*1.2), height = 10, units = "cm",
#          dpi = 300, limitsize = TRUE)
# }
# 
# 
