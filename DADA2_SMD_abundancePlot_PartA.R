# making plot with sex ratio distorters for sex ratio article

setwd("~/Dropbox/Bioinf_R/stegodyphus_microbiome_diversity/DADA2_SMD/")

################### IMPORT #####################################################

# this should be normalized data in some way
data.otu.norm <- read.csv("SMD_DADA_Fraction.OTUtable.csv", row.names=1)
data.tax.norm <- read.csv("SMD_DADA_Fraction.TAXtable.csv", row.names=1)
data.meta.norm <- read.csv("SMD_DADA_Fraction.SAMPLEtable.csv", row.names=1)

colnames(data.otu.norm) == rownames(data.meta.norm)

# make names = colony_sample species
colnames(data.otu.norm) <- paste0(data.meta.norm$colony,"_", rownames(data.meta.norm), " ", data.meta.norm$species)

###### chose your data to plot: multiple ASVs (3) #######
# asv1 <- "ASV_2" # borrelia1
# asv2 <- "ASV_41" # borrelia2
# asv3 <- "ASV_58" # # borrelia3
# 
# df <- rbind(colnames(data.otu.norm), data.otu.norm[c(asv1, asv2, asv3), ])
# rownames(df) <- c("sample", paste(asv1, data.tax.norm[asv1, "Genus"]),
#                   paste(asv2, data.tax.norm[asv1, "Genus"]), paste(asv3, data.tax.norm[asv1, "Genus"]))

#### chose, only 1 ASV ########
asv1 <- "ASV_34" # sciuri

df <- rbind(colnames(data.otu.norm), data.otu.norm[asv1, ])
rownames(df) <- c("sample", paste(asv1, data.tax.norm[asv1, "Genus"]))


df <- as.data.frame(t(df))
library(reshape2)
test_dfm <- melt(df,id.vars = 1)
test_dfm$value <- as.numeric(test_dfm$value)

# split 
hej <- strsplit(as.character(test_dfm$sample), " ")
hej2 <- do.call(rbind, hej)
dfm <- cbind(hej2[, 1], as.character(test_dfm$variable), test_dfm$value, hej2[, 2])
colnames(dfm) <- c("sample", "variable", "value", "species")
dfm <- as.data.frame(dfm)
dfm$value <- as.numeric(as.character(dfm$value))

# publication ready figure attempt
# from: https://rpubs.com/Koundy/71792

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =0, size = 10),
            axis.title.x = element_text(vjust = -0.2, size = 10),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            #panel.grid.minor = element_line(colour="#f0f0f0"),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(2, "cm"),
            legend.title = element_text(size = 8),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

library(ggplot2)
library(gridExtra)


ggplot(dfm, aes(x = sample, y = value*100)) + # data x100 to get percent
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge") + 
  scale_fill_brewer(palette="Set1") +
  ylab("Percentage of total reads") + 
  xlab("Individual spiders") + 
  labs(fill = "Symbiont ASV:") +
  theme_Publication() +
  theme(legend.text = element_text(size = 8, face = "italic")) + 
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust= 0.5, size = 4)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(panel.grid.major.x= element_blank()) + # gets rid of grid lines
  scale_y_continuous(minor_breaks = seq(0, 100, 25), breaks = seq(0, 100, 25)) + # places the grid lines that are left
  facet_grid(.~species, scales="free_x") + # this makes it so the x values do not have to be the same
  theme(strip.text = element_text(face = "bold.italic")) # italic spider species labels


# ggsave("ASVDistribution.pdf", 
#        scale = 1, width = 21, height = 10.5, units = "cm",
#        dpi = 300, limitsize = TRUE)
