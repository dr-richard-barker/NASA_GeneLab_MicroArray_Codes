##
# generic template for loading the raw count data
##

PROJECT_DIR <- "C:/Users/Robot/Documents/BRIC19_DESeq"
PROJECT_NAME <- basename(PROJECT_DIR)

setwd(PROJECT_DIR)

source("source/lib/info.R")

raw <- read.table(file=paste("detailed_results/",PROJECT_NAME, "_raw_countstable.txt", sep=""), header=T, row.names=1, sep="\t")
raw <- raw[, order(colnames(raw))]

# make sure there are no missing values (NAs) in any column
raw[is.na(raw)] <- 0


groups <- as.factor( 
    as.factor( c(rep("Col-0FL", 4), rep("Col-0GC", 4), rep("Cvi-0FL", 3), rep("Cvi-0GC", 3), rep("Ler-0FL", 3), rep("Ler-0GC", 3), rep("WS-2FL", 4), rep("WS-2GC", 4), rep("cml242-FL", 3), rep("cml242-GC", 3), rep("cml244-FL", 3), rep("cml244-GC", 3)))
)


## template for DESeq-based analysis
library(DESeq)

source("source/lib/ellipse.R")
source("source/lib/robinVennDiagram.R")
source("source/lib/malowess.R")

###################################
# settings to be completed by GUI #
###################################

P_VAL_ADJUST_METHOD <- "BH"
P_VAL_CUTOFF        <- 0.05
MIN_LFC2        	<- ifelse(1, 1, 0)

cds <- newCountDataSet( raw, groups)
cds <- estimateSizeFactors( cds )
replicated <- names(which(tapply(conditions(cds), conditions(cds), length) > 1))
if (length(replicated) < 1) {
	print("no replication within conditions")
	cds <- estimateDispersions( cds, method="pooled" )
} else {
	test <- try( estimateDispersions( cds, method="pooled", fitType="parametric" )) # there are more options availbale - include in next update
	
	if (class(test) == "try-error") {
		print("parametric fit type failed - falling back to locFit and blind")
		cds <- estimateDispersions( cds, method="blind", fitType="local", sharingMode="fit-only" )
	} else {
		cds <- test
	}
}

plotDispEsts <- function( cds ) {
	plot(
	rowMeans( counts( cds, normalized=TRUE ) ),
	fitInfo(cds)$perGeneDispEsts,
	pch = '.', log="xy" )
	xg <- 10^seq( -.5, 5, length.out=300 )
	lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}

#png(file="qualitychecks/dispersion_fit_diag.png",
#        width=800,
#        height=400)
#        par(mfrow=c(1,2))
#plotDispEsts(cds)
#dev.off()


## end checking

# test for diff exp between two conditions
# always does BH adjusting

## define pairwise contrasts
contrast.table <- matrix(c("Ler-0GC", "Ler-0FL",
"Cvi-0GC", "Cvi-0FL",
"WS-2FL", "WS-2GC",
"Col-0GC", "Col-0FL",
"Col-0GC", "cml242-FL",
"Col-0GC", "cml244-FL",
"cml242-GC", "cml242-FL",
"cml244-GC", "cml244-FL"
),byrow=T, ncol=2)

# compute them
result <- c()
res.cols <- character(0)

for (i in 1:nrow(contrast.table)) {
	cond1 <- contrast.table[i, 1]
	cond2 <- contrast.table[i, 2]
	contrast <- paste(cond1, "-", cond2, sep="")
	
	print(paste("computing contrast", contrast))
	detail.result <- nbinomTest(cds, cond1, cond2)
	detail.result$padj <- p.adjust(detail.result $pval, method=P_VAL_ADJUST_METHOD)
	
	if (FALSE %in% (rownames(result) == detail.result$id)) {
		stop("data table rownames inconsistent")		
	}
	
	## write the detailed result for each contrast to file
	write.table(detail.result, file=paste("detailed_results/full_table_", contrast, ".txt", sep=""),
				sep="\t",
				quote=F,
				row.names=F)
	
	valuecol <- detail.result$log2FoldChange
	signcol <- rep(0, length(valuecol))
	
	signcol[ which( (detail.result$padj < P_VAL_CUTOFF) & (detail.result$log2FoldChange >= MIN_LFC2) ) ] <- 1
	signcol[ which( (detail.result$padj < P_VAL_CUTOFF) & (detail.result$log2FoldChange <= MIN_LFC2 * -1) ) ] <- -1
	
	result <- cbind(result, valuecol, signcol)
	res.cols <- c(res.cols, contrast, paste("d_", contrast, sep="") )
	
	rownames(result) <- detail.result$id
	colnames(result) <- res.cols	
	print(head(result))
	png(file=paste(sep="","plots/MAplot_", contrast, ".png"),
        width=600,
        height=600)
	plotMAlowess(	M=detail.result$log2FoldChange, 
					A=detail.result$baseMean, 
					log="x", 
					stats=T, 
					array="none", 
					sig=signcol, 
					main=contrast)
	dev.off()
}

# some overview plots on the variance stabilized version of the data
cds.blind <- estimateDispersions( cds, method="blind", fitType="local", sharingMode="fit-only" )
vsd <- getVarianceStabilizedData( cds.blind ) 


## groups PCA
pca <- prcomp(t(vsd))

# plot principal components 1 and 2
png(file=paste(sep="","plots/PCAplot_", ".png"),
        width=800,
        height=600)

layout(matrix(1:2, nc = 2), c(3,1))
pca.info <- summary(pca)

plot(pca$x[,1:2],
	main="Principal component analysis",
	col="white",
	xlab=paste("PC1:", round(pca.info$imp[2,1]*100, digits=2), "%"),
	ylab=paste("PC2:", round(pca.info$imp[2,2]*100, digits=2), "%") )


for (i in 1:length(levels(groups))) {

	points(pca$x[groups == levels(groups)[i], 1], pca$x[groups == levels(groups)[i], 2], col=i, pch=i )

	center.x <- mean(pca$x[groups ==levels(groups)[i], 1])
	center.y <- mean(pca$x[groups ==levels(groups)[i], 2])

	a <- (max(pca$x[groups ==levels(groups)[i], 1]) - min(pca$x[groups ==levels(groups)[i], 1]))
	b <- (max(pca$x[groups ==levels(groups)[i], 2]) - min(pca$x[groups ==levels(groups)[i], 2]))

	ellipse(center=c(center.x, center.y), radius=c(a,b), rotate=0, add=T, col=i, lty="dotted")
}
par(mai = c(0,0,1.01,0))
plot(1:10, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
legend("topleft", legend=levels(groups), col=c(1:length(levels(groups))), pch=c(1: length(levels(groups))), bty="n", cex=0.75)
dev.off()

## plot again as PDF
pdf(file=paste(sep="","plots/PCAplot", ".pdf"), width=8, height=6)
layout(matrix(1:2, nc = 2), c(3,1))
pca.info <- summary(pca)

plot(pca$x[,1:2],
	main="Principal component analysis",
	col="white",
	xlab=paste("PC1:", round(pca.info$imp[2,1]*100, digits=2), "%"),
	ylab=paste("PC2:", round(pca.info$imp[2,2]*100, digits=2), "%") )


for (i in 1:length(levels(groups))) {

	points(pca$x[groups == levels(groups)[i], 1], pca$x[groups == levels(groups)[i], 2], col=i, pch=i )

	center.x <- mean(pca$x[groups ==levels(groups)[i], 1])
	center.y <- mean(pca$x[groups ==levels(groups)[i], 2])

	a <- (max(pca$x[groups ==levels(groups)[i], 1]) - min(pca$x[groups ==levels(groups)[i], 1]))
	b <- (max(pca$x[groups ==levels(groups)[i], 2]) - min(pca$x[groups ==levels(groups)[i], 2]))

	ellipse(center=c(center.x, center.y), radius=c(a,b), rotate=0, add=T, col=i, lty="dotted")
}
par(mai = c(0,0,1.01,0))
plot(1:10, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
legend("topleft", legend=levels(groups), col=c(1:length(levels(groups))), pch=c(1: length(levels(groups))), bty="n", cex=0.75)
dev.off()

## hierarchical clustering based on person correlation
png(file=paste(sep="","plots/hclust", ".png"),
        width=600,
        height=600)        
plot( hclust(as.dist(1-cor(vsd))) )
dev.off()

## Venn diagrams
vennres <- as.matrix(result[, seq(from=2, to=ncol(result), by=2)])
colnames(vennres) <- sub("^d_", "", colnames(result)[seq(from=2, to=ncol(result), by=2)])

if (nrow(contrast.table) <= 4) {
    png(filename="plots/vennDiagram_total.png", height=6, width=10, units="in", res=150)
    par(cex = 0.75)
    robinVennDiagram(vennres, main="Significantly regulated genes")
    dev.off()

    png(filename="plots/vennDiagram_up.png", height=6, width=10, units="in", res=150)
    par(cex = 0.75)
    robinVennDiagram(vennres, main="Upregulated genes", include="up")
    dev.off()

    png(filename="plots/vennDiagram_down.png", height=6, width=10, units="in", res=150)
    par(cex = 0.75)
    robinVennDiagram(vennres, main="Downregulated genes", include="down")
    dev.off()
}

# now write the results for each comparison to file
# the resulting table will have a pair of columns for
# each constrast. The first coumns contains the normalized
# log fold change and the second the result of the significance
# test for each gene: 1=significantly upregulated, 0=not sig.
# -1=sig. downregulated.

header <- c("Identifier", colnames(result))
result.file.name <- paste(PROJECT_NAME, "_results.txt", sep="")

write.table(as.list(header), file=result.file.name, sep="\t", quote=F, row.names=F, col.names=F)
write.table(result, file=result.file.name, sep="\t", quote=F, col.names=F, append=T)

write.short.sessionInfo("source/R.session.info.txt")
