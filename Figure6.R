### Libraries
library("data.table")
library("dplyr")
library("ggplot2")
library("grid")
library("gridExtra")
library("irr")
library("openxlsx")
library("plyr")
library("randomForest")
library("reshape2")
library("stringr")
library("yardstick")

### Load environmental data 
Coo = read.xlsx("Coordinates_Medits_2022.xlsx", sep = ";")

### load metabarcoding data
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

data.metabarcoding = AllTaxa[, -c(2:12)]
rownames(data.metabarcoding) = data.metabarcoding$scientific_name
data.metabarcodingw = transpose(data.metabarcoding, 
                               make.names = "scientific_name")
data.metabarcodingw$Staz = colnames(data.metabarcoding)[-1]
data.metabarcodingw$Staz = str_remove_all(data.metabarcodingw$Staz, "_EtOH")
data.metabarcodingw$Staz = str_remove_all(data.metabarcodingw$Staz, "_Si")

data.metabarcodingw = left_join(data.metabarcodingw, Coo)
data.metabarcodingw$GSA = str_replace_all(data.metabarcodingw$GSA,
                                          "Sard", "GSA11")
data.metabarcodingw$GSA = str_replace_all(data.metabarcodingw$GSA,
                                          "Adr", "GSA17")
data.metabarcodingw$GSA = str_replace_all(data.metabarcodingw$GSA,
                                          "Tyrr", "GSA09")
data.metabarcodingw$Depth = factor(c("0-100", "100-200", "200-400", "400-800")[
  findInterval(data.metabarcodingw$Depth, c(0, 100, 200, 400, 800))],
  levels = c("0-100", "100-200", "200-400", "400-800"))
data.metabarcodingw = data.metabarcodingw[,-which(colnames(data.metabarcodingw) %in% c("Staz",
                                                                                       "Haul",
                                                                                       "Latitude",
                                                                                       "Longitude",
                                                                                       "Distance"))]

### Load catch data
data.catch = read.csv("Species_Medits_2022_vertebrate.csv", sep = ";")[, c(1, 2, 6)]
data.catch = dcast(data = data.catch, Haul ~ scientific_name)
data.catch[which(is.na(data.catch), arr.ind = T)] = 0
colnames(data.catch)[1] = "Staz"
data.catchw = left_join(data.catch, Coo)
data.catchw$GSA = str_replace_all(data.catchw$GSA,
                                          "Sard", "GSA11")
data.catchw$GSA = str_replace_all(data.catchw$GSA,
                                          "Adr", "GSA17")
data.catchw$GSA = str_replace_all(data.catchw$GSA,
                                          "Tyrr", "GSA09")
data.catchw$Depth = factor(c("0-100", "100-200", "200-400", "400-800")[
  findInterval(data.catchw$Depth, c(0, 100, 200, 400, 800))],
  levels = c("0-100", "100-200", "200-400", "400-800"))
data.catchw = data.catchw[,-which(colnames(data.catchw) %in% c("Staz",
                                                               "Haul",
                                                               "Latitude",
                                                               "Longitude",
                                                               "Distance"))]


### Combine catch and metabarcoding data 
##merge meta and catch data
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
MetaData$presence <- 1*(MetaData$presence>0)
CatchTaxa <- melt(data.catch, id.vars = "Staz")
CatchTaxa$value = 1*(CatchTaxa$value>0)
colnames(CatchTaxa) <- c("Staz", "Species", "presence")
CatchTaxa$Source <- "Catch"
CatchTaxa$REP <- rep(0, nrow(CatchTaxa))
MetaCatch <- rbind(MetaData[,-6], CatchTaxa)
MetaCatch <- dcast(MetaCatch, Staz + Source + REP ~ Species, value.var="presence", fun.aggregate = sum)

data.metabarcoding.catch = left_join(MetaCatch, Coo)
data.metabarcoding.catch$GSA = str_replace_all(data.metabarcoding.catch$GSA,
                                          "Sard", "GSA11")
data.metabarcoding.catch$GSA = str_replace_all(data.metabarcoding.catch$GSA,
                                          "Adr", "GSA17")
data.metabarcoding.catch$GSA = str_replace_all(data.metabarcoding.catch$GSA,
                                          "Tyrr", "GSA09")
data.metabarcoding.catch$Depth = factor(c("0-100", "100-200", "200-400", "400-800")[
  findInterval(data.metabarcoding.catch$Depth, c(0, 100, 200, 400, 800))],
  levels = c("0-100", "100-200", "200-400", "400-800"))
data.metabarcoding.catch = data.metabarcoding.catch[,-which(colnames(data.metabarcoding.catch) %in% c("Staz", "Source", "Rep", 
                                                                                       "Haul",
                                                                                       "Latitude",
                                                                                       "Longitude",
                                                                                       "Distance"))]

### GSA
n = 100
train_size = 0.8

data.catchw_gsa = data.catchw[,-which(colnames(data.catchw) == "Depth")]
colnames(data.catchw_gsa) = str_replace_all(colnames(data.catchw_gsa), " ", "_")
data.catchw_gsa$GSA = factor(data.catchw_gsa$GSA)
predictions_catch = observations_catch = character(0)

for(i in 1:n){
  training_set = sample(1:nrow(data.catchw), round(nrow(data.catchw) * train_size, 0))
  test_set = setdiff(1:nrow(data.catchw), training_set)
  rf_prod <- randomForest(data.catchw_gsa[training_set, - which(colnames(data.catchw) == "GSA")],
                          data.catchw_gsa[training_set, "GSA"], 
                          mtry = 7, 
                          nodesize = 5, 
                          ntree = 1000, 
                          importance = TRUE,
                          keep.forest = TRUE)
  
  preds = as.character(predict(rf_prod, data.catchw_gsa[test_set, - which(colnames(data.catchw) == "GSA")]))
  obs = as.character(data.catchw_gsa[test_set, "GSA"])
  predictions_catch = c(predictions_catch, preds)
  observations_catch = c(observations_catch, obs)
}

k_catch = kappa2(data.frame(observations_catch, predictions_catch), weight = "equal")

data.metabarcodingw_gsa = data.metabarcodingw[,-which(colnames(data.metabarcodingw) == "Depth")]
colnames(data.metabarcodingw_gsa) = str_replace_all(colnames(data.metabarcodingw_gsa), " ", "_")
data.metabarcodingw_gsa$GSA = factor(data.metabarcodingw_gsa$GSA)
predictions_meta = observations_meta = character(0)

for(i in 1:n){
  training_set = sample(1:nrow(data.metabarcodingw), round(nrow(data.metabarcodingw) * train_size, 0))
  test_set = setdiff(1:nrow(data.metabarcodingw), training_set)
  rf_prod <- randomForest(data.metabarcodingw_gsa[training_set, - which(colnames(data.metabarcodingw) == "GSA")],
                          data.metabarcodingw_gsa[training_set, "GSA"], 
                          mtry = 7, 
                          nodesize = 5, 
                          ntree = 1000, 
                          importance = TRUE,
                          keep.forest = TRUE)
  preds = as.character(predict(rf_prod, data.metabarcodingw_gsa[test_set, - which(colnames(data.metabarcodingw) == "GSA")]))
  obs = as.character(data.metabarcodingw_gsa[test_set, "GSA"])
  predictions_meta = c(predictions_meta, preds)
  observations_meta = c(observations_meta, obs)
}

k_meta <- kappa2(data.frame(observations_meta, predictions_meta), weight = "equal")

data.meta.catch_gsa = data.metabarcoding.catch[,-which(colnames(data.metabarcoding.catch) == "Depth")]
colnames(data.meta.catch_gsa) = str_replace_all(colnames(data.meta.catch_gsa), " ", "_")
data.meta.catch_gsa$GSA = factor(data.meta.catch_gsa$GSA)
predictions_meta_catch = observations_meta_catch = character(0)

for(i in 1:n){
  training_set = sample(1:nrow(data.metabarcoding.catch), round(nrow(data.metabarcoding.catch) * train_size, 0))
  test_set = setdiff(1:nrow(data.metabarcoding.catch), training_set)
  rf_prod <- randomForest(data.meta.catch_gsa[training_set, - which(colnames(data.metabarcoding.catch) == "GSA")],
                          data.meta.catch_gsa[training_set, "GSA"], 
                          mtry = 7, 
                          nodesize = 5, 
                          ntree = 1000, 
                          importance = TRUE,
                          keep.forest = TRUE)
  
  preds = as.character(predict(rf_prod, data.meta.catch_gsa[test_set, - which(colnames(data.metabarcoding.catch) == "GSA")]))
  obs = as.character(data.meta.catch_gsa[test_set, "GSA"])
  predictions_meta_catch = c(predictions_meta_catch, preds)
  observations_meta_catch = c(observations_meta_catch, obs)
}

k_meta_catch = kappa2(data.frame(observations_meta_catch, predictions_meta_catch), weight = "equal")

df.meta.gsa = data.frame(observation = factor(observations_meta), 
                    prediction = factor(predictions_meta))
df.catch.gsa = data.frame(observation = factor(observations_catch), 
                     prediction = factor(predictions_catch))
df.meta.catch.gsa = data.frame(observation = factor(observations_meta_catch), 
                         prediction = factor(predictions_meta_catch))

# The confusion matrix from a single assessment set
cm_meta.gsa <- conf_mat(data = df.meta.gsa,
                   truth = "observation", estimate = "prediction")
cm_catch.gsa <- conf_mat(data = df.catch.gsa,
                   truth = "observation", estimate = "prediction")
cm_meta_catch.gsa <- conf_mat(data = df.meta.catch.gsa,
                        truth = "observation", estimate = "prediction")

cm_meta.gsa$table <- prop.table(yardstick::conf_mat(data = df.meta.gsa,
                                              truth = "observation", estimate = "prediction")$table)
cm_catch.gsa$table <- prop.table(yardstick::conf_mat(data = df.catch.gsa,
                                              truth = "observation", estimate = "prediction")$table)
cm_meta_catch.gsa$table <- prop.table(yardstick::conf_mat(data = df.meta.catch.gsa,
                                                    truth = "observation", estimate = "prediction")$table)

gsa_catch = autoplot(cm_catch.gsa, type = "heatmap") +
            scale_fill_gradient(low="#ecfae8",high = "#1d700c") + 
            ggtitle(paste0("Catch data - ", " Cohen k = ", round(k_catch$value, 2)))
gsa_meta = autoplot(cm_meta.gsa, type = "heatmap") +
            scale_fill_gradient(low="#ecfae8",high = "#1d700c") + 
            ggtitle(paste0("Metabarcoding data - ", " Cohen k = ", round(k_meta$value, 2)))
gsa_meta_catch = autoplot(cm_meta_catch.gsa, type = "heatmap") +
            scale_fill_gradient(low="#ecfae8",high = "#1d700c") + 
            ggtitle(paste0("Catch + Metabarcoding data - ", " Cohen k = ", round(k_meta_catch$value, 2)))

g.gsa = grid.arrange(gsa_catch, gsa_meta, gsa_meta_catch, nrow = 1, top=textGrob("Confusion matrix GSA", gp=gpar(fontsize=20,font=2)))

### Depth
n = 100
train_size = 0.8

data.catchw_Depth = data.catchw[,-which(colnames(data.catchw) == "GSA")]
colnames(data.catchw_Depth) = str_replace_all(colnames(data.catchw_Depth), " ", "_")
data.catchw_Depth$Depth = factor(data.catchw_Depth$Depth, levels = c("0-100", "100-200", "200-400", "400-800"))
predictions_catch = observations_catch = character(0)

for(i in 1:n){
  training_set = sample(1:nrow(data.catchw), round(nrow(data.catchw) * train_size, 0))
  test_set = setdiff(1:nrow(data.catchw), training_set)
  rf_prod <- randomForest(data.catchw_Depth[training_set, - which(colnames(data.catchw) == "Depth")],
                          data.catchw_Depth[training_set, "Depth"], 
                          mtry = 7, 
                          nodesize = 5, 
                          ntree = 1000, 
                          importance = TRUE,
                          keep.forest = TRUE)
  
  preds = as.character(predict(rf_prod, data.catchw_Depth[test_set, - which(colnames(data.catchw) == "Depth")]))
  obs = as.character(data.catchw_Depth[test_set, "Depth"])
  predictions_catch = c(predictions_catch, preds)
  observations_catch = c(observations_catch, obs)
}
  
k_catch = kappa2(data.frame(observations_catch, predictions_catch), weight = "equal")

data.metabarcodingw_Depth = data.metabarcodingw[,-which(colnames(data.metabarcodingw) == "GSA")]
colnames(data.metabarcodingw_Depth) = str_replace_all(colnames(data.metabarcodingw_Depth)," ", "_")
data.metabarcodingw_Depth$Depth = factor(data.metabarcodingw_Depth$Depth,
                                         levels = c("0-100", "100-200", "200-400", "400-800"))
predictions_meta = observations_meta = character(0)

for(i in 1:n){
  training_set = sample(1:nrow(data.metabarcodingw), round(nrow(data.metabarcodingw) * train_size, 0))
  test_set = setdiff(1:nrow(data.metabarcodingw), training_set)
  
  rf_prod <- randomForest(data.metabarcodingw_Depth[training_set, - which(colnames(data.metabarcodingw) == "Depth")],
                          data.metabarcodingw_Depth[training_set, "Depth"], 
                          mtry = 7, 
                          nodesize = 5, 
                          ntree = 1000, 
                          importance = TRUE,
                          keep.forest = TRUE)
  
  preds = as.character(predict(rf_prod, data.metabarcodingw_Depth[test_set, - which(colnames(data.metabarcodingw) == "Depth")]))
  obs = as.character(data.metabarcodingw_Depth[test_set, "Depth"])
  
  predictions_meta = c(predictions_meta, preds)
  observations_meta = c(observations_meta, obs)
}


K_meta <- kappa2(data.frame(observations_meta, predictions_meta), weight = "equal")

data.meta.catch_Depth = data.metabarcoding.catch[,-which(colnames(data.metabarcoding.catch) == "GSA")]
colnames(data.meta.catch_Depth) = str_replace_all(colnames(data.meta.catch_Depth), " ", "_")
data.meta.catch_Depth$Depth = factor(data.meta.catch_Depth$Depth, levels = c("0-100", "100-200", "200-400", "400-800"))
predictions_meta_catch = observations_meta_catch = character(0)

for(i in 1:n){
  training_set = sample(1:nrow(data.metabarcoding.catch), round(nrow(data.metabarcoding.catch) * train_size, 0))
  test_set = setdiff(1:nrow(data.metabarcoding.catch), training_set)
  rf_prod <- randomForest(data.meta.catch_Depth[training_set, - which(colnames(data.metabarcoding.catch) == "Depth")],
                          data.meta.catch_Depth[training_set, "Depth"], 
                          mtry = 7, 
                          nodesize = 5, 
                          ntree = 1000, 
                          importance = TRUE,
                          keep.forest = TRUE)
  
  preds = as.character(predict(rf_prod, data.meta.catch_Depth[test_set, - which(colnames(data.metabarcoding.catch) == "Depth")]))
  obs = as.character(data.meta.catch_Depth[test_set, "Depth"])
  predictions_meta_catch = c(predictions_meta_catch, preds)
  observations_meta_catch = c(observations_meta_catch, obs)
}

k_meta_catch = kappa2(data.frame(observations_meta_catch, predictions_meta_catch), weight = "equal")

df.meta.depth = data.frame(observation = factor(observations_meta, 
                                              levels = c("0-100", "100-200", "200-400", "400-800")), 
                         prediction = factor(predictions_meta, 
                                             levels = c("0-100", "100-200", "200-400", "400-800")))
df.catch.depth = data.frame(observation = factor(observations_catch, 
                                               levels = c("0-100", "100-200", "200-400", "400-800")), 
                          prediction = factor(predictions_catch, 
                                              levels = c("0-100", "100-200", "200-400", "400-800")))
df.meta.catch.depth = data.frame(observation = factor(observations_meta_catch, 
                                                levels = c("0-100", "100-200", "200-400", "400-800")), 
                           prediction = factor(predictions_meta_catch, 
                                               levels = c("0-100", "100-200", "200-400", "400-800")))

cm_meta.depth <- conf_mat(data = df.meta.depth,
                        truth = "observation", estimate = "prediction")
cm_catch.depth <- conf_mat(data = df.catch.depth,
                         truth = "observation", estimate = "prediction")
cm_meta_catch.depth <- conf_mat(data = df.meta.catch.depth,
                          truth = "observation", estimate = "prediction")

cm_meta.depth$table <- prop.table(yardstick::conf_mat(data = df.meta.depth,
                                                    truth = "observation", estimate = "prediction")$table)
cm_catch.depth$table <- prop.table(yardstick::conf_mat(data = df.catch.depth,
                                                     truth = "observation", estimate = "prediction")$table)
cm_meta_catch.depth$table <- prop.table(yardstick::conf_mat(data = df.meta.catch.depth,
                                                      truth = "observation", estimate = "prediction")$table)


depth_catch = autoplot(cm_catch.depth, type = "heatmap") +
  scale_fill_gradient(low="#D6EAF8",high = "#2E86C1") + 
  ggtitle(paste0("Catch data - ", " Cohen k = ", round(k_catch$value, 2)))
depth_meta = autoplot(cm_meta.depth, type = "heatmap") +
  scale_fill_gradient(low="#D6EAF8",high = "#2E86C1") + 
  ggtitle(paste0("Metabarcoding data - ", " Cohen k = ", round(k_meta$value, 2)))
depth_meta_catch = autoplot(cm_meta_catch.depth, type = "heatmap") +
  scale_fill_gradient(low="#D6EAF8",high = "#2E86C1") + 
  ggtitle(paste0("Catch + Metabarcoding data - ", " Cohen k = ", round(k_meta_catch$value, 2)))

g.depth = grid.arrange(depth_catch, depth_meta, depth_meta_catch, nrow = 1, top=textGrob("Confusion matrix depth", gp=gpar(fontsize=20,font=2)))


g.all = grid.arrange(g.gsa, g.depth, nrow = 2)