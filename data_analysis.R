source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("oligo")
library(oligo)

#Our experimental data (CEL files) are stored as a data matrix of expression values (each row corresponds to a porticul gene
#and each column correspoonds to a sample).

#In order to quanitfy the defferential expression of samples (control samples and keto samples)
#according to diet we fit a linear regression model to each row of the data matrix (each row corresponds to a porticul gene)
#this allowes us to test the hypotheses that this gene is differentially expressed in a partituarl sample.  

#get the CEL files 
muscle_celfiles <- list.files("./data/CEL",pattern = "^M", full = TRUE)
liver_celfiles <- list.files("./data/CEL",pattern = "^L", full = TRUE)

#get the raw data
muscle_rawData <- read.celfiles(muscle_celfiles)
liver_rawData <- read.celfiles(liver_celfiles)
#muscle_rawData
#exprs(muscle_rawData)

#convert an AffyBatch object into an ExpressionSet object using the robust multi-array average (RMA) expression measure.
#rma() does background correction, normalization and calculation of expression
muscle_normData <- rma(muscle_rawData)
liver_normData <- rma(liver_rawData)
#muscle_normData
#exprs(muscle_normData)

#sampleNames
muscle_sampleNames <- sampleNames(muscle_rawData)
liver_sampleNames <- sampleNames(liver_rawData)
#liver_sampleNames


#add group to samples pData
pData(liver_normData)$group <- ifelse(grepl("Kt*", liver_sampleNames), "Keto", "Control")
pData(muscle_normData)$group <- ifelse(grepl("Kt*", muscle_sampleNames), "Keto", "Control")

biocLite("limma")
library(limma)
#In order to quanitfy the defferential expression of samples (control samples and keto samples)
#according to diet (group variable) we fit a linear regression model to each row of the data matrix (each row corresponds to a porticul gene)
#this allowes us to test the hypotheses that this gene is differentially expressed in a partituarl sample.  
liver_design <- model.matrix(~pData(liver_normData)$group)
muscle_design <- model.matrix(~pData(muscle_normData)$group)
lm_fit_muscle <- lmFit(exprs(muscle_normData), muscle_design)
lm_fit_liver <- lmFit(exprs(liver_normData), liver_design)
lm_fit_muscle
lm_fit_liver <- eBayes(lm_fit_liver)
lm_fit_muscle <- eBayes(lm_fit_muscle)

#Extract a table of the top-ranked genes from a linear model fit.
muscle_tbl <- toptable(lm_fit_muscle, coef = 2, number = 500)
liver_tbl <- toptable(lm_fit_liver, coef = 2, number = 500)
head(muscle_tbl)
muscle_geneList <- rownames(muscle_tbl)
liver_geneList <- rownames(liver_tbl)

#output to files
write.table(liver_geneList, file = "liver_geneList.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(muscle_geneList, file = "muscle_geneList.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

