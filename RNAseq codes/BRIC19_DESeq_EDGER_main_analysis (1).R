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


 ## template for DESeq-based analysis - GLM version for more than 2 factors
library(edgeR)

source("source/lib/ellipse.R")
source("source/lib/robinVennDiagram.R")
source("source/lib/malowess.R")
source("source/lib/robinPlotMDS.dge.R")

###################################
# settings to be completed by GUI #
###################################

P_VAL_ADJUST_METHOD <- "BH"
P_VAL_CUTOFF        <- 0.05
MIN_LFC2_ONE        <- 1
DISPERSION          <- "auto"

design <- model.matrix(~ -1+groups)
colnames(design) <- as.character(unique(groups))

d <- DGEList(counts=raw, group=groups)
d <- calcNormFactors(d)

if (all(tabulate(groups) <= 1)) {
	## there is no replication whatsoever. We should actually
	## stop the script here and ask the user what he/she wants to
	## do. To keep the analysis running we set the common disp.
	## to a rel. high value
	d$common.dispersion <- 0.4	
} else {
	d <- estimateGLMCommonDisp(d, design)
}

if (DISPERSION == "tagwise") {
	d <- estimateGLMTagwiseDisp(d, design) 
} else if (DISPERSION == "trended") {	
	d <- estimateGLMTrendedDisp(d, design)
}

fit <- glmFit(d, design)

## define pairwise contrasts
contrast.table <- matrix(c("Col-0GC", "Col-0FL",
"cml242-GC", "cml242-FL",
"cml244-GC", "cml244-FL",
"Col-0GC", "cml244-FL",
"Col-0GC", "cml242-FL",
"WS-2GC", "WS-2FL",
"Cvi-0GC", "Cvi-0FL",
"Ler-0GC", "Ler-0FL"
),byrow=T, ncol=2)

result <- c()
res.cols <- character(0)

for (i in 1:nrow(contrast.table)) {
	comp <- c(rep(0, ncol(design)))
	comp[which(colnames(design)==contrast.table[i,1])] <- 1
	comp[which(colnames(design)==contrast.table[i,2])] <- -1
		
	cond1 <- contrast.table[i, 1]
	cond2 <- contrast.table[i, 2]
	contrast <- paste(cond1, "-", cond2, sep="")
	
	print(comp)	
	print(paste("computing contrast", contrast))	

	lrt <- glmLRT(d, fit, contrast=comp)
	top <- topTags(lrt, n=nrow(lrt), adjust.method=P_VAL_ADJUST_METHOD)	
	
	if (FALSE %in% (rownames(result) == rownames(lrt$table))) {
		stop("data table rownames inconsistent")		
	}
	
	print(top$comparison)
	print(head(top$table))
	
	## write the complete detailed result for each contrast to file
	write.table(as.list(c("Identifier", colnames(top$table))), file=paste("detailed_results/full_table_", contrast, ".txt", sep=""),
				sep="\t",
				quote=F,
				row.names=F,
				col.names=F)
	write.table(top$table, file=paste("detailed_results/full_table_", contrast, ".txt", sep=""),
				sep="\t",
				quote=F,
				row.names=T,
				col.names=F,
				append = T)
	
	## write the significantly changing tags only
	write.table(as.list(c("Identifier", colnames(top$table))), file=paste("detailed_results/significant_", contrast, ".txt", sep=""),
				sep="\t",
				quote=F,
				row.names=F,
				col.names=F)
	write.table(top$table[which(top$table$FDR < P_VAL_CUTOFF),], file=paste("detailed_results/significant_", contrast, ".txt", sep=""),
				sep="\t",
				quote=F,
				row.names=T,
				col.names=F,
				append = T)

	
	dt <- decideTestsDGE(lrt, adjust.method=P_VAL_ADJUST_METHOD, p.value=P_VAL_CUTOFF)
	
	result <- cbind(result, lrt$table$logFC, dt)
	res.cols <- c(res.cols, contrast, paste("d_", contrast, sep="") )
	colnames(result) <- res.cols
	rownames(result) <- rownames(lrt$table)
	
	sig.tags <- rownames(top$table)[which(top$table$FDR < P_VAL_CUTOFF)]
	png(file=paste(sep="","plots/MAplot_", contrast, ".png"),
        width=600,
        height=600)

        if (nrow(top$table) <= 50) {
            plotSmear(lrt, de.tags=sig.tags,	main=paste("MA plot of contrast", contrast), lowess=F ) 
        } else {
            plotSmear(lrt, de.tags=sig.tags,	main=paste("MA plot of contrast", contrast), lowess=T ) 
        }
	dev.off()
	
}


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

header <- c("Identifier")
header <- c(colnames(result))

result.file.name <- paste(PROJECT_NAME, "_results.txt", sep="")

write.table(as.list(header), file=result.file.name, sep="\t", quote=F, row.names=F, col.names=F)
write.table(result, file=result.file.name, sep="\t", quote=F, col.names=F, append=T)

write.short.sessionInfo("source/R.session.info.txt")

tryCatch (
	{
		## do an overview MDS plot
		## but only if there is replication - otherwise the function breaks
		if (length(groups) > length(levels(groups)) ) {
		    png(file="plots/MDSplot.png",
		        width=800,
		        height=600)
		    robinPlotMDS.dge(d, top=1000, main="MDS plot", groups=groups)
			dev.off()
		}
	},
	error = function (e) {
		print(str(e))
		print(traceback())
		if (length(groups) > length(levels(groups)) ) {
		    png(file="plots/MDSplot.png",
		        width=800,
		        height=600)
		    plotMDS(d, top=1000, main="MDS plot")
			dev.off()
		}
		# if that worked exit with 0
		quit(save="no", status=0)
	} 
)


