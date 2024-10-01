### Load libraries 
library("dplyr")
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

##Merge
gen_data = rbind(Teleo[,-1], Elasmo[,-1])
SpAggBest <- aggregate(data = gen_data[,c(1,2)], .~ scientific_name, FUN = max)
SpAgg <- gen_data[, c(1,12:54)]
SpAgg <- aggregate(. ~ scientific_name, data = SpAgg, mean)
TaxaAgg <- left_join(SpAggBest, unique(gen_data[,c(1,3:11)]))
AllTaxa <- left_join(TaxaAgg, SpAgg)

### Load Catch data
Catch <- read.csv("Species_Medits_2022_vertebrate.csv", stringsAsFactors = FALSE, sep = ";")

##comparison metabarcoding/catch
VennTot <- draw.pairwise.venn(area1 = length(unique(Catch$scientific_name)), 
                              area2 = length(unique(AllTaxa$scientific_name)), 
                              cross.area = length(intersect(unique(Catch$scientific_name), unique(AllTaxa$scientific_name))), 
                              category = c("Catch", "Metaprobe"),
                              col = c("black", "black"), lwd=3, font.face=rep("bold",3),
                              fontfamily = rep("Arial", 3), cat.fontfamily = "Arial")

VennTotGenus <- draw.pairwise.venn(area1 = length(unique(Catch$Genus)), 
                                   area2 = length(unique(AllTaxa$Genus)), 
                                   cross.area = length(intersect(unique(Catch$Genus), unique(AllTaxa$Genus))), 
                                   category = c("Catch", "Metaprobe"),
                                   col = c("black", "black"), lwd=3, font.face=rep("bold",3),
                                   fontfamily = rep("Arial", 3), cat.fontfamily = "Arial")

species_overlap <- intersect(unique(Catch$scientific_name), unique(AllTaxa$scientific_name))
species_onlymeta <- setdiff(unique(AllTaxa$scientific_name), unique(Catch$scientific_name))
species_onlycatch <- setdiff(unique(Catch$scientific_name), unique(AllTaxa$scientific_name))
