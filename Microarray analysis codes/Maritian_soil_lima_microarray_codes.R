

# Load Limma Package 
library("limma")

# Load targets file
targets <- readTargets("targets.txt")
## setwd("F:/NASA/Barker/Data Sets/Arabidopsis/GLDS-22/GLDS-22_microarray_E-GEOD-20109.raw.1")
setwd("C:/Users/richardbarker/Google Drive/R_&_Stats/Raw_microarray_data/GPR and text formats/E-GEOD-20109 Martian Soil_text_format/E-GEOD-20109.raw.1")

# Read in Agilent txt files
x <- read.maimages(targets, path=getwd(), source="agilent",green.only=TRUE)
#setwd(C:/Users/richardbarker/Google Drive/R_&_Stats/Raw_microarray_data/GPR and text formats/E-GEOD-20109 Martian Soil_text_format/E-GEOD-20109.raw.1)

# Create Design Matrix
f <- factor(targets$Condition, levels = unique(targets$Condition))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)

# Create contrast matrix of interest
contrast.matrix <- makeContrasts(WT_mgvscontrol_45 = WT_mg_45 - WT_control_45, caxvsWT_mg_180 = cax_mg_180 - WT_mg_180, levels=design)

# Determine best Offset value by analyzing the Bayesian prior value
calc_prior <- function(offset_val) {
  y0 <- limma::backgroundCorrect(x, method="normexp", normexp.method="rma", offset = offset_val)
  y0 <- normalizeBetweenArrays(y0, method="quantile")
  y0.ave <- avereps(y0, ID=y0$genes$ProbeName) 
  MA.y0 <- y0.ave[y0.ave$genes$ControlType==0,] 
  fit0 <- lmFit(MA.y0, design)
  fit0 <- contrasts.fit(fit0, contrast.matrix)
  eBayes(fit0)
}

fit_offset_0 <- calc_prior(0)
fit_offset_10 <- calc_prior(10)
fit_offset_25 <- calc_prior(25)

# Analyze the Bayesian prior value for each offset value
fit_offset_0$df.prior
## [1] 2.430236
fit_offset_10$df.prior
## [1] 2.294831
fit_offset_25$df.prior
## [1] 2.017256
# For most accurate results, use the offset which gave the highest prior value
diffexp <- function(offset_val, designmat, contmat, coefname) {
  y0 <- limma::backgroundCorrect(x, method="normexp", normexp.method="rma", offset = offset_val)
  y0 <- normalizeBetweenArrays(y0, method="quantile")
  y0.ave <- avereps(y0, ID=y0$genes$ProbeName) 
  MA.y0 <- y0.ave[y0.ave$genes$ControlType==0,] 
  fit0 <- lmFit(MA.y0, designmat)
  fit0 <- contrasts.fit(fit0, contmat)
  fit0 <- eBayes(fit0)
  topTable(fit0, adjust="BH", coef=coefname, genelist=MA.y0$genes, number=Inf)
}

res_WT_mgvscontrol_45_0 <- diffexp(0,design,contrast.matrix, "WT_mgvscontrol_45")
res_caxvsWT_mg_180_0 <- diffexp(0,design,contrast.matrix, "caxvsWT_mg_180")
# highlight genes with a log fold-change > 2 and p-value < 0.05 
newthreshold <- function(FC, Pval) {
  as.factor(abs(FC) > 2 & Pval < 0.05)
} 

res_WT_mgvscontrol_45_0$threshold <- newthreshold(res_WT_mgvscontrol_45_0$logFC, res_WT_mgvscontrol_45_0$adj.P.Val)
res_caxvsWT_mg_180_0$threshold <- newthreshold(res_caxvsWT_mg_180_0$logFC, res_caxvsWT_mg_180_0$adj.P.Val)

write.csv(res_WT_mgvscontrol_45_0,file="res_WT_mgvscontrol_45_0.csv")
write.csv(res_caxvsWT_mg_180_0,file="res_caxvsWT_mg_180_0.csv")
