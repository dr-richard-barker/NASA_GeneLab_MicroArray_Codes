#➢  source("http://bioconductor.org/biocLite.R")

)
a#➢	> browseVignettes("cummeRbund")
#o	This provides weblink to more information

#This is the start possition if you restart R

citation("cummeRbund")

#Find your cuffdiff output folder
#/Users/RichardBarker/Desktop/cuffdiff_out GeneLab_Sheldon/
#  Note: press / &theTAB buttonget

#######################
setwd("/Users/RichardBarker/Desktop/cuffdiff_out_GreenLine_Col-0")
setwd("/Users/RichardBarker/Desktop/Cuffdiff2_Duck_Weed_2-2015-11-20-16-19-51.1/cuffdiff_out")
setwd("/Users/RichardBarker/Desktop/Cuffdiff2_BRIC20_Col-0_3dpg-2015-11-19-16-14-28.7/cuffdiff_out")
setwd("/Users/RichardBarker/Desktop/Cuffdiff2_BRIC20_Col-0_3dpg_RepsCombined-2015-11-21-19-46-38.1/cuffdiff_out")
setwd("/Users/RichardBarker/Documents/RNAseq_data/Cuffdiff2_BRIC20_Col-0_Novel_GTF-2015-11-22-21-15-03.7/cuffdiff_out")
setwd("/Users/RichardBarker/Documents/RNAseq_data/Cuffdiff2_Col-0_BRIC20_anf_BRIC19-2015-11-22-03-40-51.9/cuffdiff_out")
setwd("/Users/RichardBarker/Documents/RNAseq_data/Cuffdiff2_Masson_0G,0M,5G,5M,20G,20M_FDR0.05_UQN_MHC-2015-07-09-19-47-56.8/cuffdiff_out/")
setwd("/Users/RichardBarker/Documents/RNAseq_data/Cuffdiff2_MS_time_series-2015-11-24-15-50-21.1/cuffdiff_out")
setwd("/Users/RichardBarker/Documents/RNAseq_data/Cuffdiff2_BRIC19_3Col_3cml24-2_3cml24-4-2015-11-27-05-44-41.9/cuffdiff_out")
setwd("/Users/RichardBarker/Documents/RNAseq_data/tch2_cuffdiff_out")
setwd("/Users/RichardBarker/Documents/RNAseq_data/cuffdiff2_sheldon_output_BRIC19")
setwd("/Users/RichardBarker/downloads/cuffdiff_out")



getwd
setwd("/Users/RichardBarker/Desktop/cuffdiff2_Sheldon_output_actual/")
setwd("/Users/RichardBarker/Desktop/tch2_cuffdiff_out/")
setwd("/Users/RichardBarker/Desktop/BRIC17_FDR0.-5_UQN_CuffDiff2/cuffdiff_out/")
setwd("/Users/RichardBarker/Desktop/mRNA_cuffdiff_3reps_PlusPairedEnd_NewAlignment-2013-02-15-11-04-05.968/cuffdiff_out/")
setwd("/Users/RichardBarker/Desktop/Cuffdiff2_Masson_0G,0M,5G,5M,20G,20M_FDR0.05_UQN_MHC-2015-07-09-19-47-56.8/cuffdiff_out/")
setwd("/Users/RichardBarker/Desktop/Cuffdiff2_Snail_slime-2015-08-29-09-29-24_4/cuffdiff_out_snail_slime")

##New Harddrive
setwd("Volumes/Barker_2016/RNAseq_data/Cuffdiff2_BRIC19_3Col_3cml24-2_3cml24-4-2015-11-27-05-44-41.9/cuffdiff_out")


getwd()

dir()
 
slime <- readCufflinks(rebuild=T)
tch2 <- readCufflinks(rebuild=T)
BRIC17 <- readCufflinks(rebuild=T)
Mason <- readCufflinks(rebuild=T)
genelab <- readCufflinkcufs(rebuild=T)

cuf
tch2
#(note: cuf is the variable name. but this could be anything, such as edgeR)
#(Note: half type a command eg disper & then press tab to see the options)
#######################
gene.features<-annotation(genes(cuf)) 
##############################
dispersionPlot(genes(cuf))
dispersionPlot(isoforms(cuf))

#Note: this shows the difference between the replicates
csBoxplot(genes(cuf))
csBoxplot(genes(cuf),replicates=T)

#Note: this shows the difference between the replicates
csDensity(genes(cuf))
csDensity(genes(cuf),replicates=T)

#Note: this shows the difference between the replicates
csDendro(genes(cuf))
csDendro(genes(cuf),replicates=T)

############################# Scatter plots
csScatterMatrix(genes(cuf))
csScatterMatrix(genes(cuf),replicates=T)
csScatter(genes(cuf),"Col_0_GC","Col_0_FL",smooth=T)
csScatter(genes(cuf),"Cvi_0_GC","Cvi_0_FL",smooth=T)
csScatter(genes(cuf),"Ler_0_GC","Ler_0_FL",smooth=T)
csScatter(genes(cuf),"WS_GC","WS_2_FL",smooth=T)


############################# Volcano plots
csVolcanoMatrix(genes(cuf))
csVolcanoMatrix(genes(cuf),"Col_0_GC","Col_0_FL",smooth=T)
csVolcano(genes(cuf),"Col_0_GC","Col_0_FL",smooth=T)

###############################What happened to this command?##########################
fpkmSCVPlot(genes(cuf))
fpkmSCVPlot(genes(genelab))

#############################
mySigMat<-sigMatrix(tch2,level='genes',alpha=0.05)
mySigMat
sigMatrix(genelab,level='genes',alpha=0.05)

MAplot(genes(cuf),"Col_0_GC", "Col_0_FL")

########### Siginificant filter
mysiggenes5GS <-getSig(Mason,"X0_GS","X5_GS")
myGeneset0GS5GS<-getGenes(Mason,mysiggenes,sampleIdList=c("X0_GS","X5_GS"))
mysiggenes5GS <-getSig(Mason,alpha=0.01,level="isoforms","X0_GS","X5_GS")
csHeatmap(myGeneset0GS5GS)
csHeatmap(mysiggenes5GS)
csHeatmap((mySigGenes),cluster='both')

csHeatmap(object, rescaling='none', clustering='none', labCol=T, labRow=T, logMode=T, pseudocount=1.0, border=FALSE, heatscale= c(low='darkred',mid='orange',high='white'))

##########
mysiggenes <-getSig(Mason)
sigGeneIds<-getSig(tch2,alpha=0.05,level="isoforms")

head(sigGeneIds)
length(sigGeneIds)
mySigGenes<-getGenes(tch2,sigGeneIds)

#################Make heat map
csHeatmap(mySigGenes)
####################
csHeatmap(myheatmapgenes)
csHeatmap((mySigGenes),cluster='both')

###MDS PLotss 
genes.MDS<-MDSplot(genes(slime))
genes.MDS
genes.MDS<-MDSplot(genes(slime),replicates=T)
genes.MDS

genes.PCA<-PCAplot(genes(cuff),"PC1","PC2")
genes.MDS<-MDSplot(genes(cuff))
genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
genes.MDS.rep<-MDSplot(genes(cuff),replicates=T)

###### See bottom for extra cool Heat map notes

######################## K-means clustering
ic<-csCluster(mySigGenes,k=4)
head(ic$cluster)
icp<-csClusterPlot(ic)
icp

ic6<-csCluster(mySigGenes,k=6)
head(ic6$cluster)
icp6<-csClusterPlot(ic6)
icp6
######################## Hand picking genes
myGeneIds<-c("AT1G74310", "AT1G64370", "AT2G32120", "AT1G05260", "AT5G51440", "AT5G44417", "AT3G12580", "AT2G38380", "AT5G52640", "AT1G73480", "AT4G08300", "AT5G12030", "AT1G07400", "AT2G46240", "AT3G32980", "AT4G01037", "AT5G10580", "AT1G61590", "AT5G09220", "AT4G11320", "AT1G54050", "AT2G32120", "AT3G16050", "AT1G14160", "AT4G11290")

myGenes<-getGenes(cuf,myGeneIds)
head(fpkm(myGenes))
head(fpkm(isoforms(myGenes)))
expressionBarplot(myGenes)



#### 
myGeneIds <- read.table(file = '~/Desktop/FLG22_flight_list.txt',
                        header = FALSE,
                        sep = 't')

db <- fpkmMatrix(genes(cuf))
db.sigGenes <- db[myGeneIds$V1,]

db.col <- logFC(db.sigGenes, 
               mutants = 'Col_0_FL',
               WT = 'Col_0_GC',
               logBase = 2,
               pseudo = 1)


db.cml.2 <- logFC(db.sigGenes, 
             mutants = 'cml24_2_FL',
             WT = 'cml24_2_GC',
             logBase = 2,
             pseudo = 1)

db.cml.4 <- logFC(db.sigGenes, 
                  mutants = 'cml24_4_FL',
                  WT = 'cml24_4_GC',
                  logBase = 2,
                  pseudo = 1)

final.table <- cbind(db.col, db.cml.2, db.cml.4)

colnames(final.table) <- c('col-0', 'cml24-2', 'cml24-4')

#### makes not log, note: if you change the log base 10 create the DB then you'll need to replace the 2 with 10
2^final.table


#Cullens
myGeneIds<-c( "AT3G43810", "AT3G51920", "AT5G37770", "AT1G18210", "AT3G10300", "AT5G06320", "AT5G54490", "AT1G66410", "AT2G04050", "AT2G39650", "AT2G39920", "AT3G25600", "AT3G01830", "AT1G21550", "AT4G24570",  "AT1G01560", "AT3G47450", "AT3G61080",  "AT2G18193", "AT3G53620", "AT4G12490", "AT5G39670", "AT5G60250")
##really low expression in root tip# "AT2G04050", "AT2G39650", "AT2G39920", "AT3G25600", "AT3G01830", "AT1G21550", "AT4G24570",  "AT1G01560", "AT3G47450", "AT3G61080",  "AT2G18193", "AT3G53620", "AT4G12490", "AT5G39670", "AT5G60250",  
myGeneIds<-c( "AT3G51920", "AT5G37770")
              
#myGeneIds<-c( "AT5G37770", "AT2G41100", "AT1G17420", "AT5G42650", "AT5G59820")
myGenes<-getGenes(cuf,myGeneIds)      
expressionBarplot(myGenes)
csHeatmap((myGenes),cluster='both')

myGeneIds<-c("AT5G37770","AT5G60250")
##CMLfamily members from cullen
myGeneIds<-c("AT5G37770", "AT5G39670", "At3g43810", "At3g51920", "At1g18210", "At1g66410", "At1g21550", "At3g01830", "At3g10300", "At3g25600", "At5g39670")

####Cullens target genes filtering for interesting space results
myGeneIds<-c( "AT2G41100", "AT2G04050", "AT2G39650", "AT2G39920", "AT5G54490")
### Col-0 immune response 
myGeneIds<-c( "AT3G46530", "AT5G44070", "AT1G29340", "AT3G52450", "AT2G35930", "AT5G18370", "AT1G66090", "AT3G44400", "AT4G12720", "AT5G22690", "AT2G42010", "AT4G23100", "AT5G41550", "AT1G28380", "AT3G45640", "AT5G06320", "AT4G19520", "AT4G11850", "AT1G56540", "AT1G29690")
#### Col-0 response to chitin
myGeneIds<-c('AT3G53200', 'AT3G52450', 'AT2G35930', 'AT4G01250', 'AT1G64380', 'AT2G13790', 'AT3G19580', 'AT5G59450', 'AT4G01350', 'AT3G45640', 'AT4G17230', 'AT4G26120', 'AT5G46910', 'AT5G03680')
             AT5G57220, AT2G30750, AT3G55800, AT1G06680, AT3G26830, AT3G04720, AT3G11630, AT1G53240, AT3G44300, AT1G02930, AT1G70690, AT1G25220, AT4G39950, AT1G19250, AT5G26920, AT4G03280, AT4G09650, AT3G49110, AT4G01050, AT5G03760, AT3G54640, AT1G02860, AT5G46050, AT1G72520, AT1G24100, AT2G30860
   
### GeneLab                       
myTopCoreGeneLabUpGeneIds<-c("AT3G12580", "AT2G32120", "AT5G52640", "AT1G74310", "AT5G51440")
myTopCoreGeneLabDownGeneIds<-c("AT1G64370", "AT1G05260", "AT5G44417", "AT2G38380", "AT4G11320", "AT3G32980")
myUpGeneIds<-c("AT1G74310", "AT2G32120", "AT5G51440", "AT3G12580", "AT5G52640", "AT1G73480", "AT5G12030", "AT1G07400", "AT2G46240", "AT4G01037", "AT1G54050", "AT2G32120", "AT3G16050", "AT1G14160", "AT4G11290")
myDownGeneIds<-c("AT1G64370", "AT1G05260", "AT5G44417", "AT2G38380", "AT4G08300", "AT3G32980", "AT5G10580", "AT1G61590", "AT5G09220", "AT4G11320")
mycore11geneIds<-c("AT3G12580", "AT2G32120", "AT5G52640", "AT1G74310", "AT5G51440","AT1G64370", "AT1G05260", "AT5G44417", "AT2G38380", "AT4G11320", "AT3G32980")
  
  
#### Irani et al., Radition genes found in all 4 genelab varieties ####
myRadiationGeneIds<-c("AT5G51440", "AT5G53250", "AT1G58340", "AT5G10580", "AT1G74310") 
myRadiationGenes<-getGenes(genelab,myRadiationGeneIds)

myRadiationGeneIds<-c("AT1G60190", "AT5G39050", "AT4G10040", "AT4G27410", "AT4G19030", "AT5G59520", "AT5G40390", "AT1G32450", "AT3G16050", "AT3G09350", "AT3G62270", "AT4G35250")
myRadiationGenes<-getGenes(genelab,myRadiationGeneIds)

####Plot
expressionBarplot(myRadiationGenes)

####################### How to target a specific gene

myGeneId<-"AT4G13590"
myGene<-getGene(cuff,myGeneId)
myGene
head(fpkm(myGene))
head(fpkm(isoforms(myGene)))

#######################Gene Line plots (check for plural "s")
gl<-expressionPlot(myGene)
gl
gl.rep<-expressionPlot(myGene,replicates=TRUE)
gl.rep
gl.iso.rep<-expressionPlot(isoforms(myGene),replicates=T)
gl.iso.rep
gl.cds.rep<-expressionPlot(CDS(myGene),replicates=T)
gl.cds.rep

########################Gene bar plots (check for plural "s")
gb<-expressionBarplot(myGene)
gb
gb.rep<-expressionBarplot(myGene,replicates=T)
gb.rep
igb<-expressionBarplot(isoforms(myGene),replicates=T)
igb
gp<-csPie(myGene,level="isoforms")
gp

########################This might make gene isoform models

gene.features<-annotation(genes(cuf)) 
myGeneIds<-"AT5G10580"
myGene<-getGene(cuff,myGeneId)
myGene

genetrack<-makeGeneRegionTrack(myGene)
plotTracks(genetrack)

######################## Find similar genes
mySimilar<-findSimilar(cuff,"PINK1",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)

myProfile<-c(500,0,400)
mySimilar2<-findSimilar(cuff,myProfile,n=10)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)



#This tool can search for co-expression patterns
source("http://bioconductor.org/biocLite.R")
biocLite("EBcoexpress")


############################################################
############### another way to make heat maps.... nearly works.....
http://sebastianraschka.com/Articles/heatmaps_in_r.html#installing
################
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

data <- read.csv("../datasets/heatmaps_in_r.csv", comment.char="#")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names 


#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0,0.8,length=100),              # for yellow
               seq(0.8,1,length=100))              # for green

# creates a 5 x 5 inch image
png("../images/heatmaps_in_r.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device


#####Prepareing for better heat map...
library(cummeRbund)
library(gplots)

cuf <- readCufflinks("cuffdiff_out",rebuild=T)
db <- fpkmMatrix(genes(cuf))
sigGenes <- getSig(cuf,
                   x = c('Col_0_GC'), #,'cml24_2_GC', 'cml24_4_GC'),
                   y = c('Col_0_FL'), #, 'cml24_2_FL', 'cml24_4_FL'),
                   alpha=0.05, 
                   level='genes')

db <- db[myGeneIds,]	


WT <- "BY4741" # Set the column name of the WT sample
mut <- c("ino80","arp5","ies6") # Set column names of test conditions you want in the analysis

logFC <- function(db,mutants,WT,logBase=2,pseudo=1) {
  if (length(WT) !=1 ) {
    stop('WT must refer to a single gene/column')
  }
  if (is.numeric(logBase)==FALSE) {
    stop('logBase must be a numeric value')
  }
  
  dbb <- (db[,mutants]+pseudo)/(db[,WT]+pseudo)
  dbb <- log(dbb,logBase)
  
  names(dbb) <- rownames(db)

  return(dbb)
  
}

logbase <- 2
ps <- 1

db1 <- logFC(db, 
             mutants = 'Col_0_FL',
             WT = 'Col_0_GC',
             logBase = logbase,
             pseudo = ps) # This does the log transformation

db2 <- logFC(db, 
             mutants = 'cml24_2_FL',
             WT = 'cml24_2_GC',
             logBase = logbase,
             pseudo = ps) # This does the log transformation

db3 <- logFC(db, 
             mutants = 'cml24_4_FL',
             WT = 'cml24_4_GC',
             logBase = logbase,
             pseudo = ps) # This does the log transformation

DB <- cbind(Col_0 = db1,
            cml24_2 = db2,
            cml24_4 = db3)

heatmap.2(
  DB,
  Rowv=TRUE,
  Colv=FALSE,
  dendrogram="row",
  trace="none",
  labRow= rownames(DB),
  density.info=c("none"),
  main="Log(Fold-Change) \nExpression Profiles"
)
