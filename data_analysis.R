source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("oligo")
library(oligo)
biocLite("pd.mouse430.2")
library(pd.mouse430.2)

muscle_celfiles <- list.files("./data/CEL",pattern = "^M", full = TRUE)
liver_celfiles <- list.files("./data/CEL",pattern = "^L", full = TRUE)
muscle_rawData <- read.celfiles(muscle_celfiles)
liver_rawData <- read.celfiles(liver_celfiles)
#rawData
liver_sampleNames <- sampleNames(liver_rawData)
muscle_sampleNames <- sampleNames(muscle_rawData)

#sampleNames
getClass("GeneFeatureSet")
# exprs(rawData)
# exprs(rawData)[1:4, 1:3]
# boxplot(rawData)
# max(exprs(rawData))
# dim(exprs(rawData))
# pData(rawData)
muscle_normData <- rma(muscle_rawData)
liver_normData <- rma(liver_rawData)
# normData
# boxplot(normData)
# exprs(normData)[1:4, 1:3]
# pData(rawData)
# dim(pData(rawData))
pData(liver_normData)$group <- ifelse(grepl("Kt*", liver_sampleNames), "Keto", "Control")
pData(muscle_normData)$group <- ifelse(grepl("Kt*", muscle_sampleNames), "Keto", "Control")
#pData(normData)
pData(muscle_normData)$group <- factor(pData(muscle_normData)$group)
pData(liver_normData)$group <- factor(pData(liver_normData)$group)
#pData(normDaa)$group
biocLite("limma")
library(limma)
liver_design <- model.matrix(~pData(liver_normData)$group)
muscle_design <- model.matrix(~pData(muscle_normData)$group)
# dim(des,,,,,ign_diet)
# head(design_diet)
lm_fit_muscle <- lmFit(exprs(muscle_normData), muscle_design)
lm_fit_liver <- lmFit(exprs(liver_normData), liver_design)
#lm_#fit_diet
lm_fit_liver <- eBayes(lm_fit_liver)
lm_fit_muscle <- eBayes(lm_fit_muscle)
muscle_tbl <- toptable(lm_fit_muscle, coef = 2, number = 500)
liver_tbl <- toptable(lm_fit_liver, coef = 2, number = 500)
?toptable
#unlist(mget(tbl$ID, pd.mouse430.2SYMBOL))
muscle_geneList <- rownames(muscle_tbl)
liver_geneList <- rownames(liver_tbl)
# head(geneList)
# ?write.table
write.table(liver_geneList, file = "liver_geneList.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(muscle_geneList, file = "muscle_geneList.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# geneName <- rownames(toptable(lm_fit_diet, coef = 2, n = 2)[2, ])
# geneName
# mean_types <- tapply(exprs(normData)[geneName, ], pData(rawData)$group, mean)
# mean_types
