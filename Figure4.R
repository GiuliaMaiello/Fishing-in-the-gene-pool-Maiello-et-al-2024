### Load libraries 
library("ggplot2")
library("gridExtra")
library("stringr")

### Load abundance data
MetaAbund = read.csv("Metabarcoding.csv", stringsAsFactors = FALSE, sep = ",")
colnames(MetaAbund) = str_replace_all(colnames(MetaAbund), pattern = "[.]", replacement = " ")
CatchAbund <- read.csv("Species_Medits_2022_vertebrate.csv", stringsAsFactors = FALSE, sep = ",")
tot_abundance12S <- MetaAbund[,c(4:112)]
tot_abundance12S <- data.frame(species = colnames(tot_abundance12S),
                               totreads = apply(tot_abundance12S, 2, sum))
rownames(tot_abundance12S) <- NULL

##histograms by sampling site
abundance_12S <- MetaAbund[,c(1,4:112)]
abund_agg_12S <- aggregate(. ~ Staz, data = abundance_12S, mean)

#GSA11_69
abund_12S_69 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[1,-1], 2, sum))
rownames(abund_12S_69) <- NULL
abund_12S_69 <- abund_12S_69[-which(abund_12S_69$totreads == 0),]
abund_12S_69 <- abund_12S_69[order(abund_12S_69$totreads, decreasing = TRUE),]
abund_12S_69$species =factor(abund_12S_69$species,levels=abund_12S_69$species)
abund_12S_69$group <- NA
bonus_species <- setdiff(unique(abund_12S_69$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_69)) {
  if(abund_12S_69$species[i] %in% bonus_species)
  {abund_12S_69$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_69$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_69 <- sum(abund_12S_69$totreads[1:2])/sum(abund_12S_69$totreads)

HistoAbund12S_69 <- ggplot(abund_12S_69,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Sard_69") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA11_70
abund_12S_70 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[2,-1], 2, sum))
rownames(abund_12S_70) <- NULL
abund_12S_70 <- abund_12S_70[-which(abund_12S_70$totreads == 0),]
abund_12S_70 <- abund_12S_70[order(abund_12S_70$totreads, decreasing = TRUE),]
abund_12S_70$species =factor(abund_12S_70$species,levels=abund_12S_70$species)
abund_12S_70$group <- NA
bonus_species <- setdiff(unique(abund_12S_70$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_70)) {
  if(abund_12S_70$species[i] %in% bonus_species)
  {abund_12S_70$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_70$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_70 <- sum(abund_12S_70$totreads[1:5])/sum(abund_12S_70$totreads)

HistoAbund12S_70 <- ggplot(abund_12S_70,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Sard_70") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA11_71
abund_12S_71 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[3,-1], 2, sum))
rownames(abund_12S_71) <- NULL
abund_12S_71 <- abund_12S_71[-which(abund_12S_71$totreads == 0),]
abund_12S_71 <- abund_12S_71[order(abund_12S_71$totreads, decreasing = TRUE),]
abund_12S_71$species =factor(abund_12S_71$species,levels=abund_12S_71$species)
abund_12S_71$group <- NA
bonus_species <- setdiff(unique(abund_12S_71$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_71)) {
  if(abund_12S_71$species[i] %in% bonus_species)
  {abund_12S_71$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_71$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_71 <- sum(abund_12S_71$totreads[1:12])/sum(abund_12S_71$totreads)

HistoAbund12S_71 <- ggplot(abund_12S_71,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Sard_71") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA11_72
abund_12S_72 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[4,-1], 2, sum))
rownames(abund_12S_72) <- NULL
abund_12S_72 <- abund_12S_72[-which(abund_12S_72$totreads == 0),]
abund_12S_72 <- abund_12S_72[order(abund_12S_72$totreads, decreasing = TRUE),]
abund_12S_72$species =factor(abund_12S_72$species,levels=abund_12S_72$species)
abund_12S_72$group <- NA
bonus_species <- setdiff(unique(abund_12S_72$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_72)) {
  if(abund_12S_72$species[i] %in% bonus_species)
  {abund_12S_72$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_72$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_72 <- sum(abund_12S_72$totreads[1:8])/sum(abund_12S_72$totreads)

HistoAbund12S_72 <- ggplot(abund_12S_72,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Sard_72") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA11_73
abund_12S_73 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[5,-1], 2, sum))
rownames(abund_12S_73) <- NULL
abund_12S_73 <- abund_12S_73[-which(abund_12S_73$totreads == 0),]
abund_12S_73 <- abund_12S_73[order(abund_12S_73$totreads, decreasing = TRUE),]
abund_12S_73$species =factor(abund_12S_73$species,levels=abund_12S_73$species)
abund_12S_73$group <- NA
bonus_species <- setdiff(unique(abund_12S_73$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_73)) {
  if(abund_12S_73$species[i] %in% bonus_species)
  {abund_12S_73$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_73$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_73 <- sum(abund_12S_73$totreads[1:6])/sum(abund_12S_73$totreads)

HistoAbund12S_73 <- ggplot(abund_12S_73,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Sard_73") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA11_74
abund_12S_74 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[6,-1], 2, sum))
rownames(abund_12S_74) <- NULL
abund_12S_74 <- abund_12S_74[-which(abund_12S_74$totreads == 0),]
abund_12S_74 <- abund_12S_74[order(abund_12S_74$totreads, decreasing = TRUE),]
abund_12S_74$species =factor(abund_12S_74$species,levels=abund_12S_74$species)
abund_12S_74$group <- NA
bonus_species <- setdiff(unique(abund_12S_74$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_74)) {
  if(abund_12S_74$species[i] %in% bonus_species)
  {abund_12S_74$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_74$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_74 <- sum(abund_12S_74$totreads[1:14])/sum(abund_12S_74$totreads)

HistoAbund12S_74 <- ggplot(abund_12S_74,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Sard_74") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA11_75
abund_12S_75 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[7,-1], 2, sum))
rownames(abund_12S_75) <- NULL
abund_12S_75 <- abund_12S_75[-which(abund_12S_75$totreads == 0),]
abund_12S_75 <- abund_12S_75[order(abund_12S_75$totreads, decreasing = TRUE),]
abund_12S_75$species =factor(abund_12S_75$species,levels=abund_12S_75$species)
abund_12S_75$group <- NA
bonus_species <- setdiff(unique(abund_12S_75$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_75)) {
  if(abund_12S_75$species[i] %in% bonus_species)
  {abund_12S_75$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_75$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_75 <- sum(abund_12S_75$totreads[1:13])/sum(abund_12S_75$totreads)

HistoAbund12S_75 <- ggplot(abund_12S_75,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Sard_75") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA11_77
abund_12S_77 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[8,-1], 2, sum))
rownames(abund_12S_77) <- NULL
abund_12S_77 <- abund_12S_77[-which(abund_12S_77$totreads == 0),]
abund_12S_77 <- abund_12S_77[order(abund_12S_77$totreads, decreasing = TRUE),]
abund_12S_77$species =factor(abund_12S_77$species,levels=abund_12S_77$species)
abund_12S_77$group <- NA
bonus_species <- setdiff(unique(abund_12S_77$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_77)) {
  if(abund_12S_77$species[i] %in% bonus_species)
  {abund_12S_77$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_77$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_77 <- sum(abund_12S_77$totreads[1:12])/sum(abund_12S_77$totreads)

HistoAbund12S_77 <- ggplot(abund_12S_77,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Sard_77") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA11_80
abund_12S_80 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[9,-1], 2, sum))
rownames(abund_12S_80) <- NULL
abund_12S_80 <- abund_12S_80[-which(abund_12S_80$totreads == 0),]
abund_12S_80 <- abund_12S_80[order(abund_12S_80$totreads, decreasing = TRUE),]
abund_12S_80$species =factor(abund_12S_80$species,levels=abund_12S_80$species)
abund_12S_80$group <- NA
bonus_species <- setdiff(unique(abund_12S_80$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_80)) {
  if(abund_12S_80$species[i] %in% bonus_species)
  {abund_12S_80$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_80$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_80 <- sum(abund_12S_80$totreads[1:8])/sum(abund_12S_80$totreads)

HistoAbund12S_80 <- ggplot(abund_12S_80,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text( size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Sard_80") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA11_84
abund_12S_84 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[10,-1], 2, sum))
rownames(abund_12S_84) <- NULL
abund_12S_84 <- abund_12S_84[-which(abund_12S_84$totreads == 0),]
abund_12S_84 <- abund_12S_84[order(abund_12S_84$totreads, decreasing = TRUE),]
abund_12S_84$species =factor(abund_12S_84$species,levels=abund_12S_84$species)
abund_12S_84$group <- NA
bonus_species <- setdiff(unique(abund_12S_84$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_84)) {
  if(abund_12S_84$species[i] %in% bonus_species)
  {abund_12S_84$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_84$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_84 <- sum(abund_12S_84$totreads[1:12])/sum(abund_12S_84$totreads)

HistoAbund12S_84 <- ggplot(abund_12S_84,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Sard_84") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA17_110
abund_12S_110 <- data.frame(species = colnames(abund_agg_12S[-1]),
                            totreads = apply(abund_agg_12S[11,-1], 2, sum))
rownames(abund_12S_110) <- NULL
abund_12S_110 <- abund_12S_110[-which(abund_12S_110$totreads == 0),]
abund_12S_110 <- abund_12S_110[order(abund_12S_110$totreads, decreasing = TRUE),]
abund_12S_110$species =factor(abund_12S_110$species,levels=abund_12S_110$species)
abund_12S_110$group <- NA
bonus_species <- setdiff(unique(abund_12S_110$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_110)) {
  if(abund_12S_110$species[i] %in% bonus_species)
  {abund_12S_110$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_110$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_110 <- sum(abund_12S_110$totreads[1:4])/sum(abund_12S_110$totreads)

HistoAbund12S_110 <- ggplot(abund_12S_110,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text( size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Adr_110") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA17_115
abund_12S_115 <- data.frame(species = colnames(abund_agg_12S[-1]),
                            totreads = apply(abund_agg_12S[12,-1], 2, sum))
rownames(abund_12S_115) <- NULL
abund_12S_115 <- abund_12S_115[-which(abund_12S_115$totreads == 0),]
abund_12S_115 <- abund_12S_115[order(abund_12S_115$totreads, decreasing = TRUE),]
abund_12S_115$species =factor(abund_12S_115$species,levels=abund_12S_115$species)
abund_12S_115$group <- NA
bonus_species <- setdiff(unique(abund_12S_115$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_115)) {
  if(abund_12S_115$species[i] %in% bonus_species)
  {abund_12S_115$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_115$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_115 <- sum(abund_12S_115$totreads[1:16])/sum(abund_12S_115$totreads)

HistoAbund12S_115 <- ggplot(abund_12S_115,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Adr_115") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA17_118
abund_12S_118 <- data.frame(species = colnames(abund_agg_12S[-1]),
                            totreads = apply(abund_agg_12S[13,-1], 2, sum))
rownames(abund_12S_118) <- NULL
abund_12S_118 <- abund_12S_118[-which(abund_12S_118$totreads == 0),]
abund_12S_118 <- abund_12S_118[order(abund_12S_118$totreads, decreasing = TRUE),]
abund_12S_118$species =factor(abund_12S_118$species,levels=abund_12S_118$species)
abund_12S_118$group <- NA
bonus_species <- setdiff(unique(abund_12S_118$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_118)) {
  if(abund_12S_118$species[i] %in% bonus_species)
  {abund_12S_118$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_118$group[i] <- "Metabarcoding&Catch"
  }
}

abund_12S_118$group[10] <- "Metabarcoding&Catch"

Prop_118 <- sum(abund_12S_118$totreads[1:20])/sum(abund_12S_118$totreads)

HistoAbund12S_118 <- ggplot(abund_12S_118,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Adr_118") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA17_87
abund_12S_87 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[14,-1], 2, sum))
rownames(abund_12S_87) <- NULL
abund_12S_87 <- abund_12S_87[-which(abund_12S_87$totreads == 0),]
abund_12S_87 <- abund_12S_87[order(abund_12S_87$totreads, decreasing = TRUE),]
abund_12S_87$species =factor(abund_12S_87$species,levels=abund_12S_87$species)
abund_12S_87$group <- NA
bonus_species <- setdiff(unique(abund_12S_87$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_87)) {
  if(abund_12S_87$species[i] %in% bonus_species)
  {abund_12S_87$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_87$group[i] <- "Metabarcoding&Catch"
  }
}

abund_12S_87$group[1] <- "Metabarcoding&Catch"

Prop_87 <- sum(abund_12S_87$totreads[1:16])/sum(abund_12S_87$totreads)

HistoAbund12S_87 <- ggplot(abund_12S_87,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Adr_87") +
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA17_89
abund_12S_89 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[15,-1], 2, sum))
rownames(abund_12S_89) <- NULL
abund_12S_89 <- abund_12S_89[-which(abund_12S_89$totreads == 0),]
abund_12S_89 <- abund_12S_89[order(abund_12S_89$totreads, decreasing = TRUE),]
abund_12S_89$species =factor(abund_12S_89$species,levels=abund_12S_89$species)
abund_12S_89$group <- NA
bonus_species <- setdiff(unique(abund_12S_89$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_89)) {
  if(abund_12S_89$species[i] %in% bonus_species)
  {abund_12S_89$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_89$group[i] <- "Metabarcoding&Catch"
  }
}

abund_12S_89$group[4] <- "Metabarcoding&Catch"

Prop_89 <- sum(abund_12S_89$totreads[1:6])/sum(abund_12S_89$totreads)

HistoAbund12S_89 <- ggplot(abund_12S_89,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Adr_89") + 
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA9_10
abund_12S_10 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[16,-1], 2, sum))
rownames(abund_12S_10) <- NULL
abund_12S_10 <- abund_12S_10[-which(abund_12S_10$totreads == 0),]
abund_12S_10 <- abund_12S_10[order(abund_12S_10$totreads, decreasing = TRUE),]
abund_12S_10$species =factor(abund_12S_10$species,levels=abund_12S_10$species)
abund_12S_10$group <- NA
bonus_species <- setdiff(unique(abund_12S_10$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_10)) {
  if(abund_12S_10$species[i] %in% bonus_species)
  {abund_12S_10$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_10$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_10 <- sum(abund_12S_10$totreads[1:4])/sum(abund_12S_10$totreads)

HistoAbund12S_10 <- ggplot(abund_12S_10,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Tyrr_10") + 
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA9_48
abund_12S_48 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[17,-1], 2, sum))
rownames(abund_12S_48) <- NULL
abund_12S_48 <- abund_12S_48[-which(abund_12S_48$totreads == 0),]
abund_12S_48 <- abund_12S_48[order(abund_12S_48$totreads, decreasing = TRUE),]
abund_12S_48$species =factor(abund_12S_48$species,levels=abund_12S_48$species)
abund_12S_48$group <- NA
bonus_species <- setdiff(unique(abund_12S_48$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_48)) {
  if(abund_12S_48$species[i] %in% bonus_species)
  {abund_12S_48$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_48$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_48 <- sum(abund_12S_48$totreads[1:8])/sum(abund_12S_48$totreads)

HistoAbund12S_48 <- ggplot(abund_12S_48,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Tyrr_48") + 
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA9_59
abund_12S_59 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[18,-1], 2, sum))
rownames(abund_12S_59) <- NULL
abund_12S_59 <- abund_12S_59[-which(abund_12S_59$totreads == 0),]
abund_12S_59 <- abund_12S_59[order(abund_12S_59$totreads, decreasing = TRUE),]
abund_12S_59$species =factor(abund_12S_59$species,levels=abund_12S_59$species)
abund_12S_59$group <- NA
bonus_species <- setdiff(unique(abund_12S_59$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_59)) {
  if(abund_12S_59$species[i] %in% bonus_species)
  {abund_12S_59$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_59$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_59 <- sum(abund_12S_59$totreads[1:16])/sum(abund_12S_59$totreads)

HistoAbund12S_59 <- ggplot(abund_12S_59,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Tyrr_59") + 
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA9_71
abund_12S_71_9 <- data.frame(species = colnames(abund_agg_12S[-1]),
                             totreads = apply(abund_agg_12S[19,-1], 2, sum))
rownames(abund_12S_71_9) <- NULL
abund_12S_71_9 <- abund_12S_71_9[-which(abund_12S_71_9$totreads == 0),]
abund_12S_71_9 <- abund_12S_71_9[order(abund_12S_71_9$totreads, decreasing = TRUE),]
abund_12S_71_9$species =factor(abund_12S_71_9$species,levels=abund_12S_71_9$species)
abund_12S_71_9$group <- NA
bonus_species <- setdiff(unique(abund_12S_71_9$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_71_9)) {
  if(abund_12S_71_9$species[i] %in% bonus_species)
  {abund_12S_71_9$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_71_9$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_71_9 <- sum(abund_12S_71_9$totreads[1:10])/sum(abund_12S_71_9$totreads)

HistoAbund12S_71_9 <- ggplot(abund_12S_71_9,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Tyrr_71_9") + 
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA9_83
abund_12S_83 <- data.frame(species = colnames(abund_agg_12S[-1]),
                           totreads = apply(abund_agg_12S[20,-1], 2, sum))
rownames(abund_12S_83) <- NULL
abund_12S_83 <- abund_12S_83[-which(abund_12S_83$totreads == 0),]
abund_12S_83 <- abund_12S_83[order(abund_12S_83$totreads, decreasing = TRUE),]
abund_12S_83$species =factor(abund_12S_83$species,levels=abund_12S_83$species)
abund_12S_83$group <- NA
bonus_species <- setdiff(unique(abund_12S_83$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_83)) {
  if(abund_12S_83$species[i] %in% bonus_species)
  {abund_12S_83$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_83$group[i] <- "Metabarcoding&Catch"
  }
}

Prop_83 <- sum(abund_12S_83$totreads[1:17])/sum(abund_12S_83$totreads)

HistoAbund12S_83 <- ggplot(abund_12S_83,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Tyrr_83") + 
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

#GSA9_87_9
abund_12S_87_9 <- data.frame(species = colnames(abund_agg_12S[-1]),
                             totreads = apply(abund_agg_12S[21,-1], 2, sum))
rownames(abund_12S_87_9) <- NULL
abund_12S_87_9 <- abund_12S_87_9[-which(abund_12S_87_9$totreads == 0),]
abund_12S_87_9 <- abund_12S_87_9[order(abund_12S_87_9$totreads, decreasing = TRUE),]
abund_12S_87_9$species =factor(abund_12S_87_9$species,levels=abund_12S_87_9$species)
abund_12S_87_9$group <- NA
bonus_species <- setdiff(unique(abund_12S_87_9$species), unique(CatchAbund$scientific_name))
for (i in 1:nrow(abund_12S_87_9)) {
  if(abund_12S_87_9$species[i] %in% bonus_species)
  {abund_12S_87_9$group[i] <- "MetabarcodingBonus"
  } else {
    abund_12S_87_9$group[i] <- "Metabarcoding&Catch"
  }
}

abund_12S_87_9$group[12] <- "Metabarcoding&Catch"

Prop_87_9 <- sum(abund_12S_87_9$totreads[1:16])/sum(abund_12S_87_9$totreads)

HistoAbund12S_87_9 <- ggplot(abund_12S_87_9,aes(x=species,y=totreads,fill=factor(group))) + 
  geom_bar(stat='identity', width = 1) + 
  scale_fill_manual(values = c("seagreen", "darkorange")) + 
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(colour = "black", size = 8, face = "bold"),
        panel.background = element_blank()) +
  guides(x = "none") +
  ggtitle("Tyrr_89") + 
  labs(x = "taxa", y = "sqrt(reads)", col = "group")

mean_prop <- sum(Prop_10, Prop_110, Prop_115, Prop_118, Prop_48,
                 Prop_59, Prop_69, Prop_70, Prop_71, Prop_71_9, Prop_72,
                 Prop_73, Prop_74, Prop_75, Prop_77, Prop_80, Prop_83,
                 Prop_84, Prop_87, Prop_87_9, Prop_89)/21

St_dev_prop <- sd(c(Prop_10, Prop_110, Prop_115, Prop_118, Prop_48,
                    Prop_59, Prop_69, Prop_70, Prop_71, Prop_71_9, Prop_72,
                    Prop_73, Prop_74, Prop_75, Prop_77, Prop_80, Prop_83,
                    Prop_84, Prop_87, Prop_87_9, Prop_89))


allHisto <- grid.arrange(HistoAbund12S_69, HistoAbund12S_70, HistoAbund12S_71, HistoAbund12S_72, HistoAbund12S_73,
                         HistoAbund12S_74, HistoAbund12S_75, HistoAbund12S_77, HistoAbund12S_80, HistoAbund12S_84,
                         HistoAbund12S_110, HistoAbund12S_115, HistoAbund12S_118, HistoAbund12S_87, HistoAbund12S_89,
                         HistoAbund12S_10, HistoAbund12S_48, HistoAbund12S_59, HistoAbund12S_71_9, HistoAbund12S_83, 
                         HistoAbund12S_87_9)