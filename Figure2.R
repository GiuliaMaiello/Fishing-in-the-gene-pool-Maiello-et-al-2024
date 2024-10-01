### Load libraries 
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("openxlsx")
library("pracma")
library("reshape2")
library("vegan")
library("VennDiagram")

### Load metabarcoding data
##Teleo
Teleo = read.csv("TeleoDataset.csv", stringsAsFactors = FALSE, sep = ",")
Teleo$count <- apply(Teleo[,-c(1:13)], 1, FUN = sum)
Teleo <- Teleo[-which(Teleo$count == 0),-56]
Teleo <- rbind(Teleo[which(Teleo$Class %in% c('Actinopterygii', 'Chondrichthyes')),]) ##Clean dataset for selected taxa
Teleo$count <- apply(Teleo[,14:55], 1, sum)
TeleProp <- Teleo[,14:55] 
for (i in 1:ncol(TeleProp)) {
  TeleProp[which(TeleProp[,i] <= 1),i] <- 0
}
Teleo <- data.frame(Teleo[1:13], TeleProp)
##Elasmo
Elasmo = read.csv("ElasmoDataset.csv", stringsAsFactors = FALSE, sep = ",")
Elasmo$count <- apply(Elasmo[,-c(1:13)], 1, FUN = sum)
Elasmo <- Elasmo[-which(Elasmo$count == 0),-56]
Elasmo <- rbind(Elasmo[which(Elasmo$Class %in% c('Actinopterygii', 'Chondrichthyes')),]) ##Clean dataset for selected taxa
Elasmo$count <- apply(Elasmo[,14:55], 1, sum)
ElasmoProp <- Elasmo[,14:55] 
for (i in 1:ncol(ElasmoProp)) {
  ElasmoProp[which((ElasmoProp[,i] <= 1)), i] <- 0
}
Elasmo <- data.frame(Elasmo[1:13], ElasmoProp)

### Load environmental data 
Coo = read.xlsx("Coordinates_Medits_2022.xlsx", sep = ",")

### Elasmo/Teleo detection
##Teleosts
TeleoTotSp <- Teleo[,c(2,7)]
TeleoTotSp$reads <- apply(Teleo[,-c(1:13)], 1, FUN = sum)
TeleoTotSp$primer <- 'Tele02'
TeleoTotSp$species <- rep(1)
NrSpTeleo <- sum(TeleoTotSp$species)
NrReadsTeleo <- sum(TeleoTotSp$reads)
TeleoTotSp$species <- TeleoTotSp$species/NrSpTeleo
TeleoTotSp$reads <- TeleoTotSp$reads/NrReadsTeleo
TeleoAct_sp <- length(which(TeleoTotSp$Class=='Actinopterygii'))/NrSpTeleo
TeleoCondr_sp <- length(which(TeleoTotSp$Class=='Chondrichthyes'))/NrSpTeleo
TeleoAct_reads <- TeleoTotSp[which(TeleoTotSp$Class=='Actinopterygii'),]
TeleoAct_reads <- sum(TeleoAct_reads$reads)
TeleoCondr_reads <- TeleoTotSp[which(TeleoTotSp$Class=='Chondrichthyes'),]
TeleoCondr_reads <- sum(TeleoCondr_reads$reads)
TeleoCondr_reads_tot <- TeleoCondr_reads * NrReadsTeleo
##Elasmobranchs
ElasmoTotSp <- Elasmo[,c(2,7)]
ElasmoTotSp$reads <- apply(Elasmo[,-c(1:13)], 1, FUN = sum)
ElasmoTotSp$primer <- 'Elasm02'
ElasmoTotSp$species <- rep(1)
NrSpElasmo <- sum(ElasmoTotSp$species)
NrReadsElasmo <- sum(ElasmoTotSp$reads)
ElasmoTotSp$species <- ElasmoTotSp$species/NrSpElasmo
ElasmoTotSp$reads <- ElasmoTotSp$reads/NrReadsElasmo
ElasmoAct_sp <- length(which(ElasmoTotSp$Class=='Actinopterygii'))/NrSpElasmo
ElasmoCondr_sp <- length(which(ElasmoTotSp$Class=='Chondrichthyes'))/NrSpElasmo
ElasmoAct_reads <- ElasmoTotSp[which(ElasmoTotSp$Class=='Actinopterygii'),]
ElasmoAct_reads <- sum(ElasmoAct_reads$reads)
ElasmoCondr_reads <- ElasmoTotSp[which(ElasmoTotSp$Class=='Chondrichthyes'),]
ElasmoCondr_reads <- sum(ElasmoCondr_reads$reads)
ElasmoCondr_reads_tot <- ElasmoCondr_reads * NrReadsElasmo

AllTotSp <- rbind(TeleoTotSp, ElasmoTotSp)
AllTotSp$primer <- as.character(AllTotSp$primer)
AllTotSp$primer <- factor(AllTotSp$primer, levels=c("Tele02", "Elasm02"))

##Proportion test
test_reads <- prop.test(x = c(93965, 17777), n = c(1471036, 1878487), correct = TRUE)
test_species <- prop.test(x = c(12, 20), n = c(14, 58), correct = TRUE)

##BarPlot
SpBarPlot <- ggplot(data=AllTotSp, aes(y=species, x=primer, fill=factor(Class))) + 
  geom_bar(stat="identity") +
  scale_color_manual("black") +
  scale_fill_manual(values=c(alpha("darkgoldenrod1", 0.8), alpha("#339900", 0.8))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.title=element_text(size = 12)) +
  labs(fill="Class") +
  theme_bw() 
ReadsBarPlot <- ggplot(data=AllTotSp, aes(y=reads, x=primer, fill=factor(Class))) + 
  geom_bar(stat="identity") +
  scale_color_manual("black") +
  scale_fill_manual(values=c(alpha("darkgoldenrod1", 0.8), alpha("#339900", 0.8))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.title=element_text(size = 12)) +
  labs(fill="Class") +
  theme_bw() 
TotBarPlot = ggarrange(SpBarPlot, ReadsBarPlot,  
                       nrow = 1, ncol = 2, common.legend = TRUE, legend="bottom")


##VennDiagram
ElasmoSpecies <- Elasmo$scientific_name
TeleoSpecies <- Teleo$scientific_name

VennTeleEla <- draw.pairwise.venn(area1 = length(ElasmoSpecies), 
                                  area2 = length(TeleoSpecies), 
                                  cross.area = length(intersect(ElasmoSpecies, TeleoSpecies)), 
                                  category = c("Elasm02", "Tele02"),
                                  col = c("#006666", "#FF3333"), lwd=3, font.face=rep("bold",3),
                                  fill = c("#006666", "#FF3333"), ext.line.lty= "dashed", font.face=rep("bold",3),
                                  fontfamily = rep("Arial", 3), cat.fontfamily = "Arial")
species_onlyElasmo <- setdiff(ElasmoSpecies, TeleoSpecies)
species_onlyTeleo <- setdiff(TeleoSpecies, ElasmoSpecies)
species_ElaTele <- intersect(ElasmoSpecies, TeleoSpecies)

##NMDS
TeleoData <- Teleo[,c(2, 14:55)]
TeleoData <- TeleoData[which(TeleoData$scientific_name %in% species_ElaTele),]
TeleoData <- as.data.frame(t(TeleoData))
colnames(TeleoData) <- TeleoData[1,]
TeleoData <- TeleoData[-1,]
TeleoData$Staz <- rownames(TeleoData)
rownames(TeleoData) <- NULL
TeleoData <- melt(TeleoData, id.vars = "Staz")
TeleoData$value <- as.numeric(TeleoData$value)
colnames(TeleoData) <- c("Staz", "Species", "presence")
TeleoData$REP <- rep(c("EtOH","Si"), (nrow(TeleoData)/2))
Staz = character(length(TeleoData$Staz))
und_pos = numeric(length(TeleoData$Staz))
for(i in 1:length(TeleoData$Staz)){
  und_pos = unlist(gregexpr(pattern = "_", text = TeleoData$Staz[i]))[2]
  Staz[i] =  substr(TeleoData$Staz[i], 1, und_pos-1)
}
TeleoData$Staz = Staz
TeleoData$Source <- "Tele02"
TeleoData$reads <- nthroot(TeleoData$presence, 2)
TeleoData$presence <- 1*(TeleoData$presence>0)

ElasmoData <- Elasmo[,c(2, 14:55)]
ElasmoData <- ElasmoData[which(ElasmoData$scientific_name %in% species_ElaTele),]
ElasmoData <- as.data.frame(t(ElasmoData))
colnames(ElasmoData) <- ElasmoData[1,]
ElasmoData <- ElasmoData[-1,]
ElasmoData$Staz <- rownames(ElasmoData)
rownames(ElasmoData) <- NULL
ElasmoData <- melt(ElasmoData, id.vars = "Staz")
ElasmoData$value <- as.numeric(ElasmoData$value)
colnames(ElasmoData) <- c("Staz", "Species", "presence")
ElasmoData$REP <- rep(c("EtOH","Si"), (nrow(ElasmoData)/2))
Staz = character(length(ElasmoData$Staz))
und_pos = numeric(length(ElasmoData$Staz))
for(i in 1:length(ElasmoData$Staz)){
  und_pos = unlist(gregexpr(pattern = "_", text = ElasmoData$Staz[i]))[2]
  Staz[i] =  substr(ElasmoData$Staz[i], 1, und_pos-1)
}
ElasmoData$Staz = Staz
ElasmoData$Source <- "Elasmo02"
ElasmoData$presence <- as.numeric(ElasmoData$presence)
ElasmoData$reads <- nthroot(ElasmoData$presence, 2)
ElasmoData$presence <- 1*(ElasmoData$presence>0)

ElaTeleData <- rbind(TeleoData, ElasmoData)
ElaTeleData <- dcast(ElaTeleData[,-3], Staz + Source + REP ~ Species, value.var="reads", fun.aggregate = sum)

ElaTeleData <- ElaTeleData[,-3]
ElaTeleData$StSource <- paste0(ElaTeleData$Staz, "_", ElaTeleData$Source)
ElaTeleData <- aggregate(.~  StSource, ElaTeleData[,-c(1,2)], FUN=mean)
ElaTeleData <- melt(ElaTeleData, id.vars = "StSource")
Staz = character(length(ElaTeleData$StSource))
und_pos = numeric(length(ElaTeleData$StSource))
for(i in 1:length(ElaTeleData$StSource)){
  und_pos = unlist(gregexpr(pattern = "_", text = ElaTeleData$StSource[i]))[2]
  Staz[i] =  substr(ElaTeleData$StSource[i], 1, und_pos-1)
}
ElaTeleData$Staz = Staz
ElaTeleData$Source <- rep(c("Elas02","Tele02"), (nrow(ElaTeleData)/2))
ElaTeleData$value <- as.numeric(ElaTeleData$value)

ElaTeleData <- dcast(ElaTeleData[,-c(1,6)], Staz + Source ~ variable, value.var="value", fun.aggregate = sum)

##MRPP
grouping = ElaTeleData[,1]
distmatrix = vegdist(ElaTeleData[,3:ncol(ElaTeleData)], method = "bray")
mrpp_ElaTele = mrpp(dat = distmatrix, grouping = grouping, permutations = 999)

##nMDS calucation and plotting 
NMDS_ElaTele = metaMDS(ElaTeleData[,-c(1,2)], distance="bray", k = 3, trymax = 200, autotransform = F)

NMDSspecies = data.frame(Species=rownames(NMDS_ElaTele$species),
                         x=as.numeric(NMDS_ElaTele$species[,1]),
                         y=as.numeric(NMDS_ElaTele$species[,2]))

NMDSstaz = data.frame(Staz=ElaTeleData$Staz, Source=ElaTeleData$Source,
                      x=as.numeric(NMDS_ElaTele$points[,1]), 
                      y=as.numeric(NMDS_ElaTele$points[,2]))

NMDSstaz$Source = as.character(NMDSstaz$Source)
NMDSstaz_center = aggregate(.~  Staz, NMDSstaz[,c("Staz","x","y")], FUN=mean)

NMDSstazd = merge(NMDSstaz, Coo, by="Staz")
NMDSstazd_labels =  NMDSstazd[,c("x", "y", "Staz")]
NMDSstazd_labels =  aggregate(data = NMDSstazd_labels,
                              .~  Staz, FUN = "mean")
NMDSstazd_labels$type = "label"
NMDSstazd_points = unique(NMDSstazd[,c("Staz", "x", "y")])
NMDSstazd_points$type = "point"
NMDSstazd_lp = rbind(NMDSstazd_labels, NMDSstazd_points)
NMDSstazd_lp = left_join(NMDSstazd_lp, NMDSstazd[,c("Staz", "x", "y","Source")])
NMDSstazd_Teleo = NMDSstazd_lp[which((NMDSstazd_lp$Source =="Tele02")&(NMDSstazd_lp$type == "point")),c("x", "y")]
NMDSstazd_Elasmo = NMDSstazd_lp[which((NMDSstazd_lp$Source =="Elas02")&(NMDSstazd_lp$type == "point")),c("x", "y")]
colnames(NMDSstazd_Elasmo) = c("x1", "y1")
NMDSstazd_seg = cbind(NMDSstazd_Teleo, NMDSstazd_Elasmo)

NMDS_points_ElaTele =  ggplot(data = NMDSstazd_lp, aes(x=x, y=y, fill = Source)) +
  geom_hline(yintercept = 0, color = "darkgrey", size=0.5) + 
  geom_vline(xintercept = 0, color = "darkgrey", size=0.5) + 
  xlab("NMDS1") + ylab("NMDS2") + 
  geom_point(data=NMDSstazd_lp[which(NMDSstazd_lp$type == "point"),],
             aes(x=x,y=y,shape=Source,colour=Source),size=3, inherit.aes = F) + 
  scale_colour_manual(values=c("#006666", "#FF3333")) +
  geom_segment(data = NMDSstazd_seg, linetype = 2,
               aes(x = x, y = y, xend = x1, yend = y1), inherit.aes = F) + 
  geom_text_repel(data = NMDSstazd_lp[which(NMDSstazd_lp$type == "label"),], 
                  aes(x=x, y=y, label=Staz), size = 3.5, force = 10,
                  fontface = "bold", inherit.aes = F, max.overlaps = 20) + 
  theme(legend.position = "bottom") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size=12)) +
  annotate("text", x = 0.7, y = -0.5, size = 4, fontface = "bold",
           label = paste0("Stress: ",round(NMDS_ElaTele$stress,2)))
NMDS_points_ElaTele