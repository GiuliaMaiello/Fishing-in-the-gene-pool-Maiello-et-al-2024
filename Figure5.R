### Load libraries 
library("dplyr")
library("ggplot2")
library("ggrepel")
library("indicspecies")
library("iNEXT")
library("openxlsx")
library("pairwiseAdonis")
library("plyr")
library("pracma")
library("RColorBrewer")
library("reshape2")
library("vegan")

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

##Merge
gen_data = rbind(Teleo[,-1], Elasmo[,-1])
SpAggBest <- aggregate(data = gen_data[,c(1,2)], .~ scientific_name, FUN = max)
SpAgg <- gen_data[, c(1,12:54)]
SpAgg <- aggregate(. ~ scientific_name, data = SpAgg, mean)
TaxaAgg <- left_join(SpAggBest, unique(gen_data[,c(1,3:11)]))
AllTaxa <- left_join(TaxaAgg, SpAgg)

MetaData <- AllTaxa[,c(1, 13:54)]
MetaData <- as.data.frame(t(MetaData))
colnames(MetaData) <- MetaData[1,]
MetaData <- MetaData[-1,]
MetaData$Staz <- rownames(MetaData)
rownames(MetaData) <- NULL
MetaData <- melt(MetaData, id.vars = "Staz")
colnames(MetaData) <- c("Staz", "Species", "presence")
MetaData$REP <- rep(c("EtOH","Si"), (nrow(MetaData)/2))
Staz = character(length(MetaData$Staz))
und_pos = numeric(length(MetaData$Staz))
for(i in 1:length(MetaData$Staz)){
  und_pos = unlist(gregexpr(pattern = "_", text = MetaData$Staz[i]))[2]
  Staz[i] =  substr(MetaData$Staz[i], 1, und_pos-1)
}
MetaData$Staz = Staz
MetaData$Source <- "Metabarcoding"
MetaData$presence <- as.numeric(MetaData$presence)
MetaData$reads <- nthroot(MetaData$presence, 2)
MetaData$presence <- 1*(MetaData$presence>0)
BinMetaData <- dcast(MetaData[,-6], Staz + Source + REP ~ Species, value.var="presence", fun.aggregate = sum)
AbundMetaData <- dcast(MetaData[,-3], Staz + Source + REP ~ Species, value.var="reads", fun.aggregate = sum)

### Load environmental data 
Coo = read.xlsx("Coordinates_Medits_2022.xlsx", sep = ";")

##nMDS
DataNMDS <- AbundMetaData
NMDS_MetaData = metaMDS(DataNMDS[,-c(1:3)], 
                        distance="bray", k = 3, trymax = 200, 
                        autotransform = F)
NMDSspecies = data.frame(Species=rownames(NMDS_MetaData$species),
                         x=as.numeric(NMDS_MetaData$species[,1]),
                         y=as.numeric(NMDS_MetaData$species[,2]))
NMDSstaz = data.frame(Staz=DataNMDS$Staz, Source = AbundMetaData$REP[which(AbundMetaData$REP != "H20")], 
                      x=as.numeric(NMDS_MetaData$points[,1]), 
                      y=as.numeric(NMDS_MetaData$points[,2]))
NMDSstazd = merge(NMDSstaz, Coo, by="Staz")
NMDSstazd$Depth_range = c("[0-100m)", "[100-200m)", "[200-400m)", ">400m")[findInterval(NMDSstazd$Depth, c(0, 100, 200, 400, Inf))]
NMDSstazd$Depth_range <- factor(NMDSstazd$Depth_range,
                                levels = c("[0-100m)", "[100-200m)", "[200-400m)", ">400m"),ordered = TRUE)

NMDSstazd_labels =  NMDSstazd[,c("x", "y", "Staz")]
NMDSstazd_labels =  aggregate(data = NMDSstazd_labels,
                              .~  Staz, FUN = "mean")
NMDSstazd_labels$type = "label"
NMDSstazd_points = unique(NMDSstazd[,c("Staz", "x", "y")])
NMDSstazd_points$type = "point"
NMDSstazd_lp = rbind(NMDSstazd_labels, NMDSstazd_points)
NMDSstaz$Source = as.character(NMDSstazd$Source)
NMDSstaz_center = aggregate(.~  Staz, NMDSstazd[,c("Staz","x","y")], FUN=mean)
find_hull <- function(df) df[chull(df$x, df$y), ]
hulls_GSA <- ddply(NMDSstazd, .(GSA), find_hull)
hulls_GSA$GSA = factor(hulls_GSA$GSA,
                       levels = sort(unique(hulls_GSA$GSA)))
color_GSA = brewer.pal(n = length(unique(hulls_GSA$GSA)), name = "Dark2")
hulls_depth <- ddply(NMDSstazd, .(Depth_range), find_hull)
hulls_depth$Depth_range = factor(hulls_depth$Depth_range,
                                 levels = sort(unique(hulls_depth$Depth_range)))
blues_depth = brewer.pal(n = length(unique(hulls_depth$Depth_range))+1 , name = "Blues")

NMDSstazd_lp = left_join(NMDSstazd_lp, NMDSstazd[,c("Staz", "x", "y","Source")])
NMDSstazd_Si = NMDSstazd_lp[which((NMDSstazd_lp$Source =="Si")&(NMDSstazd_lp$type == "point")),c("x", "y")]
NMDSstazd_EtOH = NMDSstazd_lp[which((NMDSstazd_lp$Source =="EtOH")&(NMDSstazd_lp$type == "point")),c("x", "y")]
colnames(NMDSstazd_EtOH) = c("x1", "y1")
NMDSstazd_seg = cbind(NMDSstazd_Si, NMDSstazd_EtOH)


NMDS_GSA =  ggplot(data = NMDSstazd_lp, aes(x=x, y=y)) +
  geom_hline(yintercept = 0, color = "darkgrey", size=0.5) + 
  geom_vline(xintercept = 0, color = "darkgrey", size=0.5) + 
  xlab("NMDS1") + ylab("NMDS2") + 
  geom_point(data = NMDSstazd_lp[which(NMDSstazd_lp$type == "point"),]
             , size = 2) + 
  geom_polygon(data = hulls_GSA, aes(x=x, y=y, group=GSA, fill=GSA), 
               colour="NA", alpha = 0.7) +
  scale_fill_manual(values=color_GSA)+ 
  geom_segment(data = NMDSstazd_seg, linetype = 2,
               aes(x = x, y = y, xend = x1, yend = y1), inherit.aes = F) + 
  geom_text_repel(data = NMDSstazd_lp[which(NMDSstazd_lp$type == "label"),], 
                  aes(x=x, y=y, label=Staz), size = 3.5, force = 10,
                  inherit.aes = F, max.overlaps = 20) + 
  theme(legend.position = "bottom") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size=12)) +
  annotate("text", x = 1.0, y = -0.5, size = 4, fontface = "bold",
           label = paste0("Stress: ",round(NMDS_MetaData$stress,2)))

NMDS_depth =  ggplot(data = NMDSstazd_lp, aes(x=x, y=y)) +
  geom_hline(yintercept = 0, color = "darkgrey", size=0.5) + 
  geom_vline(xintercept = 0, color = "darkgrey", size=0.5) + 
  xlab("NMDS1") + ylab("NMDS2") + 
  geom_point(data = NMDSstazd_lp[which(NMDSstazd_lp$type == "point"),]
             , size = 2) + 
  geom_polygon(data = hulls_depth, aes(x=x, y=y, group=Depth_range, fill=Depth_range), 
               colour="NA", alpha = 0.8) +
  scale_fill_manual(values=blues_depth[-1]) + 
  geom_segment(data = NMDSstazd_seg, linetype = 2,
               aes(x = x, y = y, xend = x1, yend = y1), inherit.aes = F) + 
  geom_text_repel(data = NMDSstazd_lp[which(NMDSstazd_lp$type == "label"),], 
                  aes(x=x, y=y, label=Staz), size = 3.5, force = 10,
                  inherit.aes = F, max.overlaps = 20) + 
  theme(legend.position = "bottom") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text=element_text(size=12)) +
  annotate("text", x = 1.0, y = -0.5, size = 4, fontface = "bold",
           label = paste0("Stress: ",round(NMDS_MetaData$stress,2)))

Xs = DataNMDS[,c(1,3)]
Xs$GSA <- NMDSstazd$GSA
Xs$Depth_range <- NMDSstazd$Depth_range
Ys = DataNMDS[,-c(1:3)]
permanova = adonis2(Ys ~ Xs$Staz + Xs$REP, permutations = 9999, method = "bray")
permanova_env = adonis2(Ys ~ Xs$GSA + Xs$Depth_range, permutations = 9999, method = "bray")

pairwise.adonis(Ys, factors = Xs$GSA)
pairwise.adonis(Ys, factors = Xs$Depth_range)

##Indicator Species Analysis
IndicSp <- AbundMetaData[,c(1,4:111)]

IndicSp <- aggregate(. ~ Staz, data = IndicSp, sum)
IndicSp <- left_join(IndicSp, Coo)
IndicSp$Depth_range = c("[0-100m)", "[100-200m)", "[200-400m)", ">400m")[findInterval(IndicSp$Depth, c(0, 100, 200, 400, Inf))]
IndicSp$Depth_range <- factor(IndicSp$Depth_range,
                              levels = c("[0-100m)", "[100-200m)", "[200-400m)", ">400m"),ordered = TRUE)
abund = IndicSp[,2:109]
GSA<- IndicSp$GSA
Depth <- IndicSp$Depth_range
inv = multipatt(abund, GSA, func = "r.g", control = how(nperm=999))
summary(inv)

##Species accumulation curve
AccCurve <- BinMetaData[,c(1,4:111)]

AccAdr <- AccCurve[1:10,]
AccAdr = aggregate(.~  Staz, AccAdr, FUN=max)
AccAdr <- as.data.frame(t(AccAdr))
colnames(AccAdr) <- AccAdr[1,]
AccAdr <- AccAdr[-1,]
AccAdr$species <- rownames(AccAdr)
rownames(AccAdr) <- NULL

AccSard <- AccCurve[11:30,]
AccSard = aggregate(.~  Staz, AccSard, FUN=max)
AccSard <- as.data.frame(t(AccSard))
colnames(AccSard) <- AccSard[1,]
AccSard <- AccSard[-1,]
AccSard$species <- rownames(AccSard)
rownames(AccSard) <- NULL

AccTyrr <- AccCurve[31:42,]
AccTyrr = aggregate(.~  Staz, AccTyrr, FUN=max)
AccTyrr <- as.data.frame(t(AccTyrr))
colnames(AccTyrr) <- AccTyrr[1,]
AccTyrr <- AccTyrr[-1,]
AccTyrr$species <- rownames(AccTyrr)
rownames(AccTyrr) <- NULL

Adr <- as.matrix(apply(AccAdr[,-6],2,as.integer))
row.names(Adr) <- AccAdr[,6]
Sard <- as.matrix(apply(AccSard[,-11],2,as.integer))
row.names(Sard) <- AccSard[,11]
Tyrr <- as.matrix(apply(AccTyrr[,-7],2,as.integer))
row.names(Tyrr) <- AccTyrr[,7]

AccCurveGSA = list(Adriatic_sea = Adr, Sardinian_sea = Sard, Tyrrhenian_sea = Tyrr)

AccGSA <- iNEXT(AccCurveGSA, datatype="incidence_raw", endpoint=20)
AccGSA <- ggiNEXT(AccGSA) +
  scale_fill_manual(values=color_GSA) +
  scale_color_manual(values = color_GSA) +
  xlab("Number of sampling units") + ylab("Species richness") +
  theme_bw()
AccGSA <- AccGSA + theme(legend.position = "bottom")


AccDepth <- left_join(AccCurve, Coo)
AccDepth <- AccDepth[, c(2:109, 114)]
AccDepth = aggregate(.~ Depth, AccDepth, FUN=max)
AccDepth$Depth_range = c("[0-100m)", "[100-200m)", "[200-400m)", ">400m")[findInterval(AccDepth$Depth, c(0, 100, 200, 400, Inf))]
AccDepth$Depth_range <- factor(AccDepth$Depth_range,
                               levels = c("[0-100m)", "[100-200m)", "[200-400m)", ">400m"),ordered = TRUE)
AccDepth <- AccDepth[,-1]

AccDepth1 <- AccDepth[which(AccDepth$Depth_range == '[0-100m)'),]
AccDepth1 <- as.data.frame(t(AccDepth1))
colnames(AccDepth1) <- AccDepth1[109,]
AccDepth1 <- AccDepth1[-109,]
AccDepth1$species <- rownames(AccDepth1)
rownames(AccDepth1) <- NULL

AccDepth2 <- AccDepth[which(AccDepth$Depth_range == '[100-200m)'),]
AccDepth2 <- as.data.frame(t(AccDepth2))
colnames(AccDepth2) <- AccDepth2[109,]
AccDepth2 <- AccDepth2[-109,]
AccDepth2$species <- rownames(AccDepth2)
rownames(AccDepth2) <- NULL

AccDepth3 <- AccDepth[which(AccDepth$Depth_range == '[200-400m)'),]
AccDepth3 <- as.data.frame(t(AccDepth3))
colnames(AccDepth3) <- AccDepth3[109,]
AccDepth3 <- AccDepth3[-109,]
AccDepth3$species <- rownames(AccDepth3)
rownames(AccDepth3) <- NULL

AccDepth4 <- AccDepth[which(AccDepth$Depth_range == '>400m'),]
AccDepth4 <- as.data.frame(t(AccDepth4))
colnames(AccDepth4) <- AccDepth4[109,]
AccDepth4 <- AccDepth4[-109,]
AccDepth4$species <- rownames(AccDepth4)
rownames(AccDepth4) <- NULL

Depth1 <- as.matrix(apply(AccDepth1[,-12],2,as.integer))
row.names(Depth1) <- AccDepth1[,12]
Depth2 <- as.matrix(apply(AccDepth2[,-4],2,as.integer))
row.names(Depth2) <- AccDepth2[,4]
Depth3 <- as.matrix(apply(AccDepth3[,-4],2,as.integer))
row.names(Depth3) <- AccDepth3[,4]
Depth4 <- as.matrix(apply(AccDepth4[,-5],2,as.integer))
row.names(Depth4) <- AccDepth4[,5]

AccCurveDepth = list(Depth1 = Depth1, Depth2 = Depth2, Depth3 = Depth3, Depth4 = Depth4)

AccDepth <- iNEXT(AccCurveDepth, datatype="incidence_raw", endpoint=20)
AccDepth <- ggiNEXT(AccDepth) +
  scale_fill_manual(values=blues_depth[-1]) +
  scale_color_manual(values = blues_depth[-1]) +
  xlab("Number of sampling units") + ylab("Species richness") +
  theme_bw() 

AccDepth <- AccDepth + theme(legend.position = "bottom")
