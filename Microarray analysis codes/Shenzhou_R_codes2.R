#rm(list = ls())

install.packages("rmarkdown") 
"library("Package Name")":
  source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("affy","limma","oligo","affyPLM","annotate","GO.db"))
install.packages(c("dplyr","ggplot2","ggrepel"))

##Activate packageS lot

library("affy")
library("limma")
library("oligo")
library("affyPLM")
library("annotate")
library("dplyr")
library("parallel")
library("GO.db")
library("ggplot2")
library("ggrepel")
source("https://bioconductor.org/biocLite.R")
biocLite("ath1121501.db")
library("ath1121501.db")
library("org.At.tair.db")
library("GO.db")

### Shenzhou microarray location
setwd("C:/Users/richardbarker/Google Drive/R_&_Stats/Raw_microarray_data/E-MTAB-2518_Shenzhou8_CEL/E-MTAB-2518.raw.1")
getwd()
data <- ReadAffy()
sampleNames(data)
AEsetnorm <- affy::rma(data)

#####Create MA plots to evaluate the normalization:

for (i in 1:length(sampleNames(data)))
{name = paste("MAplotnorm",i,".jpg",sep="")
jpeg(name)
oligo::MAplot(AEsetnorm,which=i)
dev.off()}


## import file info from targets file
targets <- read.csv(file="targets.csv", stringsAsFactors=FALSE)

## 12.	Fit a robust linear model to the probe level data of the raw microarray data (used to analyze overall quality of the microarrays):
##	Subset arrays of interest:
##  data.Tissue <- data[,c(7,8,9,10,11,12,19,20,21,22,23,24)]
## sampleNames(data.Tissue) <- targets$Name[c(7,8,9,10,11,12,19,20,21,22,23,24)]

#### double check the order of the centrafuge and flight samples!!
data.ground <- data[,c(5,6)]
sampleNames(data.ground) = targets$Name[c(5,6)] #####
data.flight<- data[,c(1,2)]
sampleNames(data.flight) = targets$Name[c(1,2)]#####
data.1G_flight<- data[,c(3,4)]
sampleNames(data.1G_flight) = targets$Name[c(3,4)]######

data.ground.PLM <- fitPLM(data.ground)
data.flight.PLM <- fitPLM(data.flight)
data.1G_flight.PLM <- fitPLM(data.1G_flight)

#sampleNames(data) = targets$Name
#data.PLM <- fitPLM(data)

### 13.	Assess RNA degradation (differentially degraded RNA will produce less accurate results):
RNAdeg <- AffyRNAdeg(data.ground)
colors <- palette(rainbow(12))
plotAffyRNAdeg(RNAdeg, col=colors)
legend("topleft", sampleNames(data.ground),lty=1,col=colors,lwd=2, cex=0.5, ncol = 4)

RNAdeg <- AffyRNAdeg(data.flight)
colors <- palette(rainbow(12))
plotAffyRNAdeg(RNAdeg, col=colors)
legend("topleft", sampleNames(data.flight),lty=1,col=colors,lwd=2, cex=0.5, ncol = 4)

RNAdeg <- AffyRNAdeg(data.1G_flight)
colors <- palette(rainbow(12))
plotAffyRNAdeg(RNAdeg, col=colors)
legend("topleft", sampleNames(data.1G_flight),lty=1,col=colors,lwd=2, cex=0.5, ncol = 4)

RNAdeg <- AffyRNAdeg(data)
colors <- palette(rainbow(12))
plotAffyRNAdeg(RNAdeg, col=colors)
legend("topleft", sampleNames(data),lty=1,col=colors,lwd=2, cex=0.5, ncol = 4)


### 14.	Save the RNA degradation measurements to a readable file:
RNADegSummary <- summaryAffyRNAdeg(RNAdeg)
write.csv(RNADegSummary, file="RNADegSummary.csv")

## 15.	Analyze the quality of the Relative Log Expression (RLE) and Normalized Unscaled Standard Error (NUSE) plots:

RLE(data.PLM, main="RLE Plot", col=colors)
legend("bottomleft", sampleNames(data),lty=1,col=colors,lwd=2, cex=0.5, ncol = 4)

NUSE(data.PLM, main="NUSE Plot",col=colors)
legend("topleft", sampleNames(data),lty=1,col=c(2,3,4,5,6,7),lwd=2, cex=0.5, ncol = 4)

# 16.	Extract the experimental factors used in the study:
facs <- targets[,c(1,3)]
facs = paste(facs[,2], sep="")
f = factor(facs)

# 17.	Create a design matrix, which matches .CEL files to their corresponding conditions:
design = model.matrix(~0+f)
colnames(design) = levels(f)

# 18.	Create a contrast matrix for comparison of Spaceflight and Control microarrays:

# contrast matrix for 1G_Flight vs Ground Control
cont.matrix.1 = makeContrasts(GroundvsCenttrifuge = Ground_Control-Inflight_Centrifuge, levels=design) 

# contrast matrix for Flight vs Ground Control
cont.matrix.2 = makeContrasts(GroundvsMicrogravity = Ground_Control-Space_Flight_Microgravity, levels=design) 

# contrast matrix for 1G_Flight vs Flight
cont.matrix.3 = makeContrasts(CentrifugevsMicrogravity = Inflight_Centrifuge-Space_Flight_Microgravity, levels=design) 

# Compute estimated coefficients and standard errors
fit = lmFit(AEsetnorm, design)
fit2.1 = contrasts.fit(fit, cont.matrix.1)
fit2.1 = eBayes(fit2.1)
res1 = topTable(fit2.1, coef = "GroundvsCenttrifuge", adjust = "BH", number = Inf)

fit = lmFit(AEsetnorm, design)
fit2.2 = contrasts.fit(fit, cont.matrix.2)
fit2.2 = eBayes(fit2.2)
res2 = topTable(fit2.2, coef = "GroundvsMicrogravity", adjust = "BH", number = Inf)

fit = lmFit(AEsetnorm, design)
fit2.3 = contrasts.fit(fit, cont.matrix.3)
fit2.3 = eBayes(fit2.3)
res3 = topTable(fit2.3, coef = "CentrifugevsMicrogravity", adjust = "BH", number = Inf)

#Create an optional threshold, to highlight genes with a log fold-change > 2 and p-value < 0.05:

#For 1G_Flight vs Ground Control
res1$threshold = as.factor(abs(res1$logFC) > 2 & res1$adj.P.Val < 0.05)

#For Flight vs Ground Control
res2$threshold = as.factor(abs(res2$logFC) > 2 & res2$adj.P.Val < 0.05)

#For 1G_Flight vs Flight
res3$threshold = as.factor(abs(res3$logFC) > 2 & res3$adj.P.Val < 0.05)

#21.	Download and load the correct annotation database for the microarrays used 

#22.	Link each probe to its corresponding gene symbol and Entrez ID number:
TAIR_Link <- read.table(
  "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_NCBI_mapping_files/TAIR10_NCBI_GENEID_mapping",
  sep="\t", header=FALSE, col.names= c("EntrezID", "LocusID"))

#res1

gene.ID <- mget(rownames(res1),ath1121501ACCNUM, ifnotfound = NA)
ID <- sapply(gene.ID, paste, collapse=",")
res1$LocusID <- ID

gene.Symbol <- mget(rownames(res1),ath1121501SYMBOL, ifnotfound = NA)
Symbol <- sapply(gene.Symbol, paste, collapse=",")
res1$Symbol <- Symbol

Locus1 <- as.character(res1$LocusID)
LocusLink <- TAIR_Link$LocusID
EntrezLink <- TAIR_Link$EntrezID
x <- match(Locus1, LocusLink)
y <- 1:length(rownames(res1))
for (i in 1:length(x)){
  y[i] <- EntrezLink[x[i]]
}
res1$EntrezID <- y


#res2

gene.ID <- mget(rownames(res2),ath1121501ACCNUM, ifnotfound = NA)
ID <- sapply(gene.ID, paste, collapse=",")
res2$LocusID <- ID

gene.Symbol <- mget(rownames(res2),ath1121501SYMBOL, ifnotfound = NA)
Symbol <- sapply(gene.Symbol, paste, collapse=",")
res2$Symbol <- Symbol

Locus1 <- as.character(res2$LocusID)
LocusLink <- TAIR_Link$LocusID
EntrezLink <- TAIR_Link$EntrezID
x <- match(Locus1, LocusLink)
y <- 1:length(rownames(res2))
for (i in 1:length(x)){
  y[i] <- EntrezLink[x[i]]
}
res2$EntrezID <- y

#res3
gene.ID.3 <- mget(rownames(res3),ath1121501ACCNUM, ifnotfound = NA)
ID3 <- sapply(gene.ID.3, paste, collapse=",")
res3$LocusID <- ID3

gene.Symbol3 <- mget(rownames(res3),ath1121501SYMBOL, ifnotfound = NA)
Symbol3 <- sapply(gene.Symbol3, paste, collapse=",")
res3$Symbol <- Symbol3

Locus3 <- as.character(res3$LocusID)
LocusLink <- TAIR_Link$LocusID

EntrezLink <- TAIR_Link$EntrezID
x <- match(Locus3, LocusLink)
y <- 1:length(rownames(res3))
for (i in 1:length(x)){
  y[i] <- EntrezLink[x[i]]
}
res3$EntrezID <- y

###################################################################################

#24: Link genes to correct homology group via homologene 68
homologene <- read.table(
  "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data",
  sep="\t", quote= "", header=FALSE, col.names= c("HomologyGroup", "taxonomy_id", "EntrezID", "symbol", "protein_gi", "protein_accession"))

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

#Subset Homologene table to only contain genes from Arabidopsis
homologene1 <- subset(homologene, taxonomy_id == 3702)

##############################################################################
#Res1
# Send variables and non-base functions to the cluster
clusterExport(cl=cl, varlist=c("homologene", "res1", "homologene1"), envir=environment())

### res2
clusterExport(cl=cl, varlist=c("homologene", "res2", "homologene1"), envir=environment())

### res3
clusterExport(cl=cl, varlist=c("homologene", "res3", "homologene1"), envir=environment())
clusterEvalQ(cl, {library(dplyr)})

#############################################################################
#Res1
# Parallel functions to build groups table
Link <- parLapply(cl, res1$EntrezID, function(x) {
  subtable <- subset(homologene1, EntrezID == x)
  ifelse(nrow(subtable) == 0, return(NA), return((subtable$HomologyGroup)))})

res1$HomologyGroup <- Link

#######################################################################
Link2 <- parLapply(cl, res2$EntrezID, function(x) {
  subtable <- subset(homologene1, EntrezID == x)
  ifelse(nrow(subtable) == 0, return(NA), return((subtable$HomologyGroup)))})

res2$HomologyGroup <- Link2

#########################################################################
Link3 <- parLapply(cl, res3$EntrezID, function(x) {
  subtable <- subset(homologene1, EntrezID == x)
  ifelse(nrow(subtable) == 0, return(NA), return((subtable$HomologyGroup)))})

res3$HomologyGroup <- Link3

#########################################################################
#Stop the current cluster
stopCluster(cl)

#########################################################################
#res1
#un-list column
res1$HomologyGroup <- sapply(res1$HomologyGroup, paste, collapse="")
res1$HomologyGroup <- as.numeric(res1$HomologyGroup)

res2$HomologyGroup <- sapply(res2$HomologyGroup, paste, collapse="")
res2$HomologyGroup <- as.numeric(res2$HomologyGroup)

res3$HomologyGroup <- sapply(res3$HomologyGroup, paste, collapse="")
res3$HomologyGroup <- as.numeric(res3$HomologyGroup)

#########################################################################
# res1
#27. Save top Table results:
write.table(res1,file="E-MTAB-2518_Shenzhou_res1.txt", sep = "\t")
write.csv(res1,file="E-MTAB-2518_Shenzhou_res1.csv")
# res2
write.table(res2,file="E-MTAB-2518_Shenzhou_res2.txt", sep = "\t")
write.csv(res2,file="E-MTAB-2518_Shenzhou_res2.csv")
# res3
write.table(res3,file="E-MTAB-2518_Shenzhou_res3.txt", sep = "\t")
write.csv(res3,file="E-MTAB-2518_Shenzhou_res3.csv")

