setwd("~/Dropbox/1_misc/Richards_plots/outlines_edited")
library(Momocs)
### ### ### ### ### ### ### 
### First we read in the outlines
### ### ### ### ### ### ### 

jpg.list <- list.files(pattern = ".jpg")
import_jpg(jpg.list)->returns_leaves
Out(returns_leaves)-> leaves

### look at the outlines
panel(leaves, col=1)

### Lets read in the grouping variable which is in the folder with the outlines
read.csv("groups.csv")-> groups
groups_leaves <-data.frame(groups[,1:2])
leaves$fac <- groups_leaves

par(mfrow=c(1,3))
stack(leaves)
coo_center(leaves) -> centered_leaves
stack(centered_leaves)
coo_template(centered_leaves) -> template_leaves
stack(template_leaves)

### Now lets run the elliptic Fourier analysis
efourier(template_leaves, norm=T, nb.h=32, smooth.it=1)-> efou_leaves

### And now the PCA of the eFa harmonic coefficients
pca_leaves <- PCA(efou_leaves,  center=T)

### Here we can see the contribution of each component
summary(pca_leaves)

### Lets create a vector of colors for coloring
cols<-c("red","blue")

plant_col <- character(length(groups$plant)) 
plant_col[groups$plant == "control"] <- "blue"
plant_col[groups$plant == "treatment"] <- "red"


### Now lets save a pdf with the results
pdf("leaves_Richard.pdf")
panel(leaves, col=plant_col, names="plant", cex=0.3)
plot(pca_leaves, "plant",xax = 1, yax = 2, ellipsesax=T, density=F, contour=F,ellipses=T, col=cols, eigen=F, rug=F, delaunay=F, loadings=F, pch=16,cex=0.4, points=T)
plot(pca_leaves, "plant", pos="li",xax = 1, yax = 2, ellipsesax=T, density=F, size.shp=0.5,contour=F,ellipses=T, col=cols, eigen=F, rug=F, delaunay=F, loadings=F, cex=0.75, points=F)
PCcontrib(pca_leaves)
boxplot(pca_leaves, "plant",cex.legend = 0.75)
mshapes(efou_leaves, "plant") -> mshapes_pca_leaves
mshapes_pca_leaves$shp -> mshapes
tps_iso(mshapes$control, mshapes$treatment)

par(mfrow=c(1,2))                    
coo_plot(mshapes$control, col="blue", border=1, lwd=1, centroid=F, xy.axis=F, first.point=F)
title('Mean shape control')
coo_plot(mshapes$treatment, col="red", border=1, lwd=1, centroid=F, xy.axis=F, first.point=F)
title('Mean shape treatment')


dev.off()
