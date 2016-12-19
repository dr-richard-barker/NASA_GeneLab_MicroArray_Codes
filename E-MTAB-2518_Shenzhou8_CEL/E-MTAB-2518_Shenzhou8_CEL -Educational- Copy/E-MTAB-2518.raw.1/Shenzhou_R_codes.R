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
library("GO.db")
library("ggplot2")
library("ggrepel")


### Chitin microarray location
setwd("C:/Users/richardbarker/Google Drive/R & Stats/Raw microarray data/E-MTAB-2518_Shenzhou8_CEL/E-MTAB-2518.raw.1")
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
data.flight<- data[,c(1,2)]
data.1G_flight<- data[,c(3,4)]

data.ground.PLM <- fitPLM(data.ground)
data.flight.PLM <- fitPLM(data.flight)
data.1G_flight.PLM <- fitPLM(data.1G_flight)
data.PLM <- fitPLM(data)

### 13.	Assess RNA degradation (differentially degraded RNA will produce less accurate results):
RNAdeg <- AffyRNAdeg(data)
colors <- palette(rainbow(12))
plotAffyRNAdeg(RNAdeg, col=colors)
legend("topleft", sampleNames(data),lty=1,col=colors,lwd=2, cex=0.5, ncol = 4)

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


### 14.	Save the RNA degradation measurements to a readable file:
RNADegSummary <- summaryAffyRNAdeg(RNAdeg)
write.csv(RNADegSummary, file="RNADegSummary.csv")

## 15.	Analyze the quality of the Relative Log Expression (RLE) and Normalized Unscaled Standard Error (NUSE) plots:

RLE(data.PLM, main="RLE Plot", col=colors)
legend("bottomleft", sampleNames(data),lty=1,col=colors,lwd=2, cex=0.5, ncol = 4)

NUSE(data.PLM, main="NUSE Plot",col=colors)
legend("topleft", sampleNames(data),lty=1,col=c(2,3,4,5,6,7),lwd=2, cex=0.5, ncol = 4)

#????# 16.	Extract the experimental factors used in the study:
facs <- targets[,c(1,3)]
facs = paste(facs[,2], sep="")
f = factor(facs)

#????# 17.	Create a design matrix, which matches .CEL files to their corresponding conditions:
design = model.matrix(~0+f)
colnames(design) = levels(f)

#????# 18.	Create a contrast matrix for comparison of Spaceflight and Control microarrays:
##cont.matrix = makeContrasts(Tissue_SpaceflightvsNorm = Spaceflight_Tissue-Ground_Tissue, levels=design)

#????# contrast matrix for Flight vs Ground Control
cont.matrix = makeContrasts(Flight = Ground, levels=design)

#????# contrast matrix for 1G_Flight vs Ground Control
cont.matrix = makeContrasts(1G_Flight = Ground, levels=design)

#????# contrast matrix for 1G_Flight vs Flight
cont.matrix = makeContrasts(1G_Flight = Flight, levels=design)




