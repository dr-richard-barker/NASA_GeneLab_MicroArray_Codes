### R code from vignette source '/Library/Frameworks/R.framework/Versions/3.3/Resources/library/HTqPCR/doc/HTqPCR.Rnw'

###################################################
### code chunk number 1: Prepare parameters
###################################################
options(width=65)
set.seed(123)


###################################################
### code chunk number 2: Load package
###################################################
library("HTqPCR")


###################################################
### code chunk number 3: Extract R code (eval = FALSE)
###################################################
## all.R.commands <- system.file("doc", "HTqPCR.Rnw", package = "HTqPCR")
## Stangle(all.R.commands)


###################################################
### code chunk number 4: All functions
###################################################
ls("package:HTqPCR")


###################################################
### code chunk number 5: Load example data
###################################################
data(qPCRraw)
data(qPCRpros)
class(qPCRraw)


###################################################
### code chunk number 6: Information contained in qPCRsets
###################################################
slotNames(qPCRraw)
phenoData(qPCRraw)
pData(qPCRraw)
pData(qPCRraw) <- data.frame(Genotype=rep(c("A", "B"), each=3), Replicate=rep(1:3, 2))
pData(qPCRraw)
featureData(qPCRraw)
head(fData(qPCRraw))


###################################################
### code chunk number 7: Example input files
###################################################
path <- system.file("exData", package="HTqPCR")
head(read.delim(file.path(path, "files.txt")))


###################################################
### code chunk number 8: Read raw data
###################################################
files <- read.delim(file.path(path, "files.txt"))
raw <- readCtData(files=files$File, path=path)


###################################################
### code chunk number 9: Show qPCRset data object
###################################################
show(raw)


###################################################
### code chunk number 10: Ct overview ex 1
###################################################
g <- featureNames(raw)[1:10]
plotCtOverview(raw, genes=g, xlim=c(0,50), groups=files$Treatment, conf.int=TRUE, 
ylim=c(0,55))


###################################################
### code chunk number 11: Ct overview ex 2
###################################################
plotCtOverview(raw, genes=g, xlim=c(0,50), groups=files$Treatment, 
calibrator="Control")


###################################################
### code chunk number 12: HTqPCR.Rnw:242-245
###################################################
par(mfrow=c(2,1))
g <- featureNames(raw)[1:10]
plotCtOverview(raw, genes=g, xlim=c(0,50), groups=files$Treatment, conf.int=TRUE, 
ylim=c(0,55))
plotCtOverview(raw, genes=g, xlim=c(0,50), groups=files$Treatment, 
calibrator="Control")


###################################################
### code chunk number 13: Ct card ex 1
###################################################
plotCtCard(raw, col.range=c(10,35), well.size=2.6)


###################################################
### code chunk number 14: Ct card ex 2
###################################################
featureClass(raw) <- factor(c("Marker", "TF", "Kinase")[sample(c(1,1,2,2,1,3), 
384, replace=TRUE)])
plotCtCard(raw, plot="class", well.size=2.6)


###################################################
### code chunk number 15: HTqPCR.Rnw:267-268
###################################################
plotCtCard(raw, col.range=c(10,35), well.size=2.6)


###################################################
### code chunk number 16: HTqPCR.Rnw:270-271
###################################################
featureClass(raw) <- factor(c("Marker", "TF", "Kinase")[sample(c(1,1,2,2,1,3), 
384, replace=TRUE)])
plotCtCard(raw, plot="class", well.size=2.6)


###################################################
### code chunk number 17: Ct replicates
###################################################
plotCtReps(qPCRraw, card=2, percent=20)


###################################################
### code chunk number 18: HTqPCR.Rnw:288-289
###################################################
plotCtReps(qPCRraw, card=2, percent=20)


###################################################
### code chunk number 19: Ct variation ex 1
###################################################
raw.mix	<- raw
exprs(raw.mix)[,6]	<- sample(exprs(raw[,6]))
plotCtVariation(raw.mix, variation="sd", log=TRUE, main="SD of replicated features", col="lightgrey")


###################################################
### code chunk number 20: Ct variation ex 2
###################################################
raw.variation	<- plotCtVariation(raw.mix, type="detail", add.featurenames=TRUE, pch=" ", cex=1.2)


###################################################
### code chunk number 21: Ct variation ex 2 in detail
###################################################
names(raw.variation)
head(raw.variation[["Var"]][,1:4])
head(raw.variation[["Mean"]][,1:4])
apply(raw.variation[["Var"]][,3:7], 2, summary)
colSums(raw.variation[["Var"]][,3:7]>20)


###################################################
### code chunk number 22: HTqPCR.Rnw:327-328
###################################################
raw.mix	<- raw
exprs(raw.mix)[,6]	<- sample(exprs(raw[,6]))
plotCtVariation(raw.mix, variation="sd", log=TRUE, main="SD of replicated features", col="lightgrey")


###################################################
### code chunk number 23: HTqPCR.Rnw:330-331
###################################################
raw.variation	<- plotCtVariation(raw.mix, type="detail", add.featurenames=TRUE, pch=" ", cex=1.2)


###################################################
### code chunk number 24: Set Ct categories
###################################################
raw.cat <- setCategory(raw, groups=files$Treatment, quantile=0.8)


###################################################
### code chunk number 25: Plot Ct categories ex 1
###################################################
plotCtCategory(raw.cat)


###################################################
### code chunk number 26: Plot Ct categories ex 2
###################################################
plotCtCategory(raw.cat, stratify="class")


###################################################
### code chunk number 27: HTqPCR.Rnw:371-374
###################################################
par(mfrow=c(2,1))
plotCtCategory(raw.cat)
plotCtCategory(raw.cat, stratify="class")


###################################################
### code chunk number 28: Plot Ct categories ex 3
###################################################
plotCtCategory(raw.cat, by.feature=TRUE, cexRow=0.1)


###################################################
### code chunk number 29: HTqPCR.Rnw:389-390
###################################################
plotCtCategory(raw.cat, by.feature=TRUE, cexRow=0.1)


###################################################
### code chunk number 30: Normalise data
###################################################
q.norm <- normalizeCtData(raw.cat, norm="quantile")
sr.norm <- normalizeCtData(raw.cat, norm="scale.rank")
nr.norm <- normalizeCtData(raw.cat, norm="norm.rank")
d.norm <- normalizeCtData(raw.cat, norm="deltaCt", deltaCt.genes=c("Gene1", "Gene60"))
g.norm <- normalizeCtData(raw.cat, norm="geometric.mean")


###################################################
### code chunk number 31: Normalisation comparison
###################################################
plot(exprs(raw), exprs(q.norm), pch=20, main="Quantile normalisation", col=rep(brewer.pal(6, "Spectral"), each=384))


###################################################
### code chunk number 32: HTqPCR.Rnw:442-458
###################################################
col <- rep(brewer.pal(6, "Spectral"), each=384)
col2 <- brewer.pal(5, "Dark2")
par(mfrow=c(3,2), mar=c(2,2,2,2))
# All methods individually
plot(exprs(raw), exprs(q.norm), pch=20, main="Quantile normalisation", col=col)
plot(exprs(raw), exprs(sr.norm), pch=20, main="Rank invariant scaling", col=col)
plot(exprs(raw), exprs(nr.norm), pch=20, main="Rank invariant normalisation", col=col)
plot(exprs(raw), exprs(d.norm), pch=20, main="deltaCt normalisation", col=col)
plot(exprs(raw), exprs(g.norm), pch=20, main="Geometric mean normalisation", col=col)
# Just a single sample, across methods
plot(exprs(raw)[,3], exprs(q.norm)[,3], pch=20, col=col2[1], main="Comparison of methods for sample 3", ylim=c(-10,40))
points(exprs(raw)[,3], exprs(sr.norm)[,3], pch=20, col=col2[2])
points(exprs(raw)[,3], exprs(nr.norm)[,3], pch=20, col=col2[3])
points(exprs(raw)[,3], exprs(d.norm)[,3], pch=20, col=col2[4])
points(exprs(raw)[,3], exprs(g.norm)[,3], pch=20, col=col2[5])
legend(8, 40, legend=c("Quantile", "Rank.invariant scaling", "Rank.invariant normalization", "deltaCt", "Geometric.mean"), col=col2, lwd=2, bty="n")


###################################################
### code chunk number 33: Subset data
###################################################
nr.norm[1:10,]
nr.norm[,c(1,3,5)]


###################################################
### code chunk number 34: Filter data 1
###################################################
qFilt   <- filterCtData(nr.norm, remove.type="Endogenous Control")
qFilt   <- filterCtData(nr.norm, remove.name=c("Gene1", "Gene20", "Gene30"))
qFilt   <- filterCtData(nr.norm, remove.class="Kinase")
qFilt   <- filterCtData(nr.norm, remove.type=c("Endogenous Control"), remove.name=c("Gene1", "Gene20", "Gene30"))


###################################################
### code chunk number 35: Filter data 2
###################################################
qFilt   <- filterCtData(nr.norm, remove.category="Undetermined")
qFilt   <- filterCtData(nr.norm, remove.category="Undetermined", n.category=5)


###################################################
### code chunk number 36: IQR plot
###################################################
iqr.values <- apply(exprs(nr.norm), 1, IQR)
hist(iqr.values, n=20, main="", xlab="IQR across samples")
abline(v=1.5, col=2)


###################################################
### code chunk number 37: Filter data 3
###################################################
qFilt   <- filterCtData(nr.norm, remove.IQR=1.5)


###################################################
### code chunk number 38: HTqPCR.Rnw:511-512
###################################################
iqr.values <- apply(exprs(nr.norm), 1, IQR)
hist(iqr.values, n=20, main="", xlab="IQR across samples")
abline(v=1.5, col=2)


###################################################
### code chunk number 39: Ct correlations
###################################################
plotCtCor(raw, main="Ct correlation")


###################################################
### code chunk number 40: HTqPCR.Rnw:538-539
###################################################
plotCtCor(raw, main="Ct correlation")


###################################################
### code chunk number 41: Summary of Ct values
###################################################
summary(raw)


###################################################
### code chunk number 42: Ct density
###################################################
plotCtDensity(sr.norm)


###################################################
### code chunk number 43: Ct histogram
###################################################
plotCtHistogram(sr.norm)


###################################################
### code chunk number 44: HTqPCR.Rnw:569-572
###################################################
par(mfrow=c(1,2), mar=c(3,3,2,1))
plotCtDensity(sr.norm)
plotCtHistogram(sr.norm)


###################################################
### code chunk number 45: HTqPCR.Rnw:583-590
###################################################
par(mfrow=c(3,2), mar=c(2,2,2,1))
plotCtDensity(qPCRraw, main="Raw Ct values")
plotCtDensity(q.norm, main="quantile")
plotCtDensity(sr.norm, main="scale.rankinvariant")
plotCtDensity(nr.norm, main="norm.rankinvariant")
plotCtDensity(d.norm, main="deltaCt")
plotCtDensity(g.norm, main="geometric.mean")


###################################################
### code chunk number 46: Ct boxes
###################################################
plotCtBoxes(sr.norm, stratify="class")


###################################################
### code chunk number 47: HTqPCR.Rnw:605-606
###################################################
plotCtBoxes(sr.norm, stratify="class")


###################################################
### code chunk number 48: Ct scatter ex 1
###################################################
plotCtScatter(sr.norm, cards=c(1,2), col="type", diag=TRUE)


###################################################
### code chunk number 49: Ct scatter ex 2
###################################################
plotCtScatter(sr.norm, cards=c(1,4), col="class", diag=TRUE)


###################################################
### code chunk number 50: HTqPCR.Rnw:626-629
###################################################
par(mfrow=c(1,2), mar=c(3,3,2,1))
plotCtScatter(sr.norm, cards=c(1,2), col="type", diag=TRUE)
plotCtScatter(sr.norm, cards=c(1,4), col="class", diag=TRUE)


###################################################
### code chunk number 51: Ct pairs
###################################################
plotCtPairs(sr.norm, col="type", diag=TRUE)


###################################################
### code chunk number 52: HTqPCR.Rnw:647-648
###################################################
plotCtPairs(sr.norm, col="type", diag=TRUE)


###################################################
### code chunk number 53: Ct heatmap
###################################################
plotCtHeatmap(raw, gene.names="", dist="euclidean")


###################################################
### code chunk number 54: HTqPCR.Rnw:667-668
###################################################
plotCtHeatmap(raw, gene.names="", dist="euclidean")


###################################################
### code chunk number 55: CV across samples
###################################################
plotCVBoxes(qPCRraw, stratify="class")
plotCVBoxes(qPCRraw, stratify="type")


###################################################
### code chunk number 56: HTqPCR.Rnw:686-689
###################################################
par(mfrow=c(1,2), mar=c(2,2,2,1))
plotCVBoxes(qPCRraw, stratify="class")
plotCVBoxes(qPCRraw, stratify="type")


###################################################
### code chunk number 57: Cluster Ct
###################################################
clusterCt(sr.norm, type="samples")


###################################################
### code chunk number 58: HTqPCR.Rnw:715-716
###################################################
clusterCt(sr.norm, type="samples")


###################################################
### code chunk number 59: Plot subclusters
###################################################
cluster.list <- clusterCt(sr.norm, type="genes", n.cluster=6, cex=0.5)


###################################################
### code chunk number 60: Show subcluster
###################################################
c6 <- cluster.list[[6]]
print(c6)
show(sr.norm[c6,])


###################################################
### code chunk number 61: HTqPCR.Rnw:735-736
###################################################
cluster.list <- clusterCt(sr.norm, type="genes", n.cluster=6, cex=0.5)


###################################################
### code chunk number 62: Principal components analysis
###################################################
plotCtPCA(qPCRraw)
plotCtPCA(qPCRraw, features=FALSE)


###################################################
### code chunk number 63: HTqPCR.Rnw:754-757
###################################################
par(mfrow=c(1,2), mar=c(2,2,2,1))
plotCtPCA(qPCRraw)
plotCtPCA(qPCRraw, features=FALSE)


###################################################
### code chunk number 64: Object history
###################################################
getCtHistory(sr.norm)
getCtHistory(qFilt)


###################################################
### code chunk number 65: Perform standard t-test
###################################################
qDE.ttest <- ttestCtData(sr.norm[,1:4], groups=files$Treatment[1:4], 
calibrator="Control")
head(qDE.ttest)


###################################################
### code chunk number 66: Perform Mann-Whitney test
###################################################
qDE.mwtest <- mannwhitneyCtData(sr.norm[,1:4], groups=files$Treatment[1:4], 
calibrator="Control")
head(qDE.mwtest)


###################################################
### code chunk number 67: Perform limma test
###################################################
# Preparing experiment design
design <- model.matrix(~0+files$Treatment)
colnames(design) <- c("Control", "LongStarve", "Starve")
print(design)
contrasts <- makeContrasts(LongStarve-Control, LongStarve-Starve, 
Starve-Control, (Starve+LongStarve)/2-Control, levels=design)
colnames(contrasts) <- c("LS-C", "LS-S", "S-C", "bothS-C")
print(contrasts)
# Reorder data to get the genes in consecutive rows
sr.norm2    <- sr.norm[order(featureNames(sr.norm)),] 
qDE.limma       <- limmaCtData(sr.norm2, design=design, contrasts=contrasts, 
ndups=2, spacing=1)


###################################################
### code chunk number 68: limma test output
###################################################
class(qDE.limma)
names(qDE.limma)
head(qDE.limma[["LS-C"]])


###################################################
### code chunk number 69: limma summary output
###################################################
qDE.limma[["Summary"]][21:30,]


###################################################
### code chunk number 70: Relative quantification ex 1
###################################################
plotCtRQ(qDE.ttest, genes=1:15)


###################################################
### code chunk number 71: Relative quantification ex 2
###################################################
plotCtRQ(qDE.limma, p.val=0.085, transform="log10", col="#9E0142")


###################################################
### code chunk number 72: HTqPCR.Rnw:867-868
###################################################
plotCtRQ(qDE.ttest, genes=1:15)


###################################################
### code chunk number 73: HTqPCR.Rnw:877-878
###################################################
plotCtRQ(qDE.limma, p.val=0.085, transform="log10", col="#9E0142")


###################################################
### code chunk number 74: Significant Ct
###################################################
plotCtSignificance(qDE.limma, q=sr.norm, 
groups=files$Treatment, target="LongStarve", 
calibrator="Control", genes=featureNames(sr.norm)[11:20], un.col="#3288BD", jitter=0.2)


###################################################
### code chunk number 75: HTqPCR.Rnw:897-898
###################################################
plotCtSignificance(qDE.limma, q=sr.norm, 
groups=files$Treatment, target="LongStarve", 
calibrator="Control", genes=featureNames(sr.norm)[11:20], un.col="#3288BD", jitter=0.2)


###################################################
### code chunk number 76: Heatmap significant Ct
###################################################
heatmapSig(qDE.limma, dist="euclidean")


###################################################
### code chunk number 77: HTqPCR.Rnw:915-916
###################################################
heatmapSig(qDE.limma, dist="euclidean")


###################################################
### code chunk number 78: Multiple samples per card
###################################################
# Example with 2  or 4 samples per 384 well card.
sample2.order	<- rep(c("subSampleA", "subSampleB"), each=192)
sample4.order	<- rep(c("subA", "subB", "subC", "subD"), each=96)
# Splitting the data into all individual samples
qPCRnew2 <- changeCtLayout(sr.norm, sample.order=sample2.order)
show(qPCRnew2)
qPCRnew4 <- changeCtLayout(sr.norm, sample.order=sample4.order)
show(qPCRnew4)


###################################################
### code chunk number 79: Card history
###################################################
getCtHistory(qPCRnew4)


###################################################
### code chunk number 80: Combine qPCRset objects
###################################################
q.comb	<- cbind(q.norm[,1:3], sr.norm[,4], nr.norm[,c(1,5,6)])
q.comb
q.comb2	<- rbind(q.norm, sr.norm[1:4,], nr.norm)
q.comb2


###################################################
### code chunk number 81: Combined card history
###################################################
getCtHistory(q.comb)


###################################################
### code chunk number 82: Example SDS data
###################################################
path <- system.file("exData", package="HTqPCR")
cat(paste(readLines(file.path(path, "SDS_sample.txt"), n=19), "\n"))


###################################################
### code chunk number 83: Example SDS data 2
###################################################
readLines(file.path(path, "SDS_sample.txt"), n=20)


###################################################
### code chunk number 84: SDS format
###################################################
path <- system.file("exData", package = "HTqPCR")
raw <- readCtData(files="SDS_sample.txt", path=path, format="SDS")
show(raw)


###################################################
### code chunk number 85: SDS format
###################################################
path <- system.file("exData", package = "HTqPCR")
raw <- readCtData(files="LightCycler_sample.txt", path=path, format="LightCycler")
show(raw)


###################################################
### code chunk number 86: CFX format
###################################################
path <- system.file("exData", package = "HTqPCR")
raw <- readCtData(files="CFX_sample.txt", path=path, format="CFX", n.features=330)
show(raw)


###################################################
### code chunk number 87: BioMark format
###################################################
exPath <- system.file("exData", package="HTqPCR")
raw1 <- readCtData(files="BioMark_sample.csv", path=exPath, format="BioMark", n.features=48, n.data=48) 
dim(raw1)
raw2 <- readCtData(files="BioMark_sample.csv", path=exPath, format="BioMark", n.features=48*48, n.data=1) 
dim(raw2)


###################################################
### code chunk number 88: Ct microfluidic array
###################################################
plotCtArray(raw1)


###################################################
### code chunk number 89: HTqPCR.Rnw:1072-1073
###################################################
plotCtArray(raw1)


###################################################
### code chunk number 90: OpenArray format
###################################################
exPath <- system.file("exData", package="HTqPCR")
raw1 <- readCtData(files="OpenArray_sample.csv", path=exPath, format="OpenArray", n.features=846, n.data=6) 
dim(raw1)
raw2 <- readCtData(files="OpenArray_sample.csv", path=exPath, format="OpenArray", n.features=846*6, n.data=1) 
dim(raw2)


###################################################
### code chunk number 91: Create Fluidigm set 1
###################################################
# Get example data
exPath <- system.file("exData", package="HTqPCR")
exFiles <- "BioMark_sample.csv"
# Reading data into a data frame
temp	<- read.delim(file.path(exPath, exFiles), skip=11, sep=",", colClasses="character")
n	<- 48
# Turn into matrix
mat	<- matrix(as.numeric(temp$Value), ncol=n, nrow=n, byrow=FALSE)
mat[mat>40]	<- NA
# Create qPCRset
raw <- new("qPCRset", exprs=mat, featureCategory=as.data.frame(array("OK", c(n,n))))
sampleNames(raw) <- paste("S", 1:n, sep="") 
featureNames(raw)	<- paste("A", 1:n, sep="")


###################################################
### code chunk number 92: Create Fluidigm set 2
###################################################
# Create qPCRset object
temp	<- readCtData(exFiles, path=exPath, n.features=48*48, column.info=list(flag=9, feature=5, type=6, Ct=7, position=1), skip=12, sep=",")
# Re-format from 1x2304 samples in input file into 48x48 as on array
raw	<- changeCtLayout(temp, sample.order=rep(1:48, each=48))


###################################################
### code chunk number 93: Check HTqPCR news
###################################################
news(Version>1.7, package="HTqPCR")


###################################################
### code chunk number 94: sessionInfo
###################################################
toLatex(sessionInfo())


