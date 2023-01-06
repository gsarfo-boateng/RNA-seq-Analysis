getwd()
countData <- read.table("Bombus.txt", header = TRUE)
head(countData)
rownames(countData) = countData$Geneid
countData = countData[, -1]
PhenoData <- read.table("SamT.txt", header = TRUE)

columns = PhenoData$Geneid 
colnames(countData) = columns
#BeeCol = PhenoData$Beecolony

#treatment = (PhenoData$Chill)
colData = as.data.frame(cbind(colnames(countData), PhenoData)) 
colData = colData[, -2] 
#PhenoData$Beecolony, PhenoData$Beecast, PhenoData$Chill, PhenoData$Mate, PhenoData$Tissue
library(DESeq2)
DDS <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~Beecolony+Beecast+Chill+Mate+Tissue)
DDS = DESeq(DDS)
nrow(DDS)

#Remove zeros from count less than 