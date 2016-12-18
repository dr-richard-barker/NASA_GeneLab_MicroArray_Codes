source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")

#G Yu, LG Wang, Y Han, QY He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology 2012, 16(5):284-287.

#To view documentation for the version of this package installed in your system, start R and enter: 
browseVignettes("clusterProfiler")

mycore11genelab<-c("AT3G12580", "AT2G32120", "AT5G52640", "AT1G74310", "AT5G51440","AT1G64370", "AT1G05260", "AT5G44417", "AT2G38380", "AT4G11320", "AT3G32980")
mycoreGeneLabUp<-c("AT3G12580", "AT2G32120", "AT5G52640", "AT1G74310", "AT5G51440")
mycoreGeneLabDownGene<-c("AT1G64370", "AT1G05260", "AT5G44417", "AT2G38380", "AT4G11320", "AT3G32980")

myGenelab_2out3<-c("AT1G74310", "AT1G64370", "AT2G32120", "AT1G05260", "AT5G51440", "AT5G44417", "AT3G12580", "AT2G38380", "AT5G52640", "AT1G73480", "AT4G08300", "AT5G12030", "AT1G07400", "AT2G46240", "AT3G32980", "AT4G01037", "AT5G10580", "AT1G61590", "AT5G09220", "AT4G11320", "AT1G54050", "AT2G32120", "AT3G16050", "AT1G14160", "AT4G11290")
myUpGenelab_2out3<-c("AT1G74310", "AT2G32120", "AT5G51440", "AT3G12580", "AT5G52640", "AT1G73480", "AT5G12030", "AT1G07400", "AT2G46240", "AT4G01037", "AT1G54050", "AT2G32120", "AT3G16050", "AT1G14160", "AT4G11290")
myDownGenelab_2out3<-c("AT1G64370", "AT1G05260", "AT5G44417", "AT2G38380", "AT4G08300", "AT3G32980", "AT5G10580", "AT1G61590", "AT5G09220", "AT4G11320")

#bitr: Biological Id TranslatoR (Note: Arabidopsis is supported)
x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2", 
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", annoDb="org.Hs.eg.db")
head(eg)

#User should provides an annotation package, both fromType and toType can accept any types that supported.
#User can use idType to list all supporting types.
idType("org.Hs.eg.db")
#We can translate from one type to other types.
ids <- bitr(x, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), annoDb="org.Hs.eg.db")
head(ids)
library("DOSE")
data(geneList)
View(geneList)
gene <- names(mycoreGeneLabUp)[abs(mycoreGeneLabUp) > 2]
head(gene)

#  Gene Ontology Classification
ggo <- groupGO(gene     = gene,
               organism = "Arabidopsis",
               ont      = "BP",
               level    = 3,
               readable = TRUE)
head(summary(ggo))

# GO over-representation test
ego2 <- gseGO(geneList     = geneList,
              organism     = "Arabidopsis",
              ont          = "CC",
              nPerm        = 1000,
              minGSSize    = 120,
              pvalueCutoff = 0.01,
              verbose      = FALSE)

#KEGG analysis KEGG = "Kyoto Encyclopedia of genes and Genomes"

kk <- enrichKEGG(gene         = gene,
                 organism     = "Arabidopsis",
                 pvalueCutoff = 0.05, 
                 readable     = TRUE,
                 use_internal_data = TRUE)
head(summary(kk))

#KEGG Gene Set Enrichment Analysis

kk2 <- gseKEGG(geneList     = geneList,
               organism     = "Arabidopsis",
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.01,
               verbose      = FALSE,
               use_internal_data = TRUE)
head(summary(kk2))

#Reactome pathway analysis
#ReactomePA uses Reactome as a source of pathway data. The function call of enrichPathway and gsePathway in ReactomePA is consistent with enrichKEGG and gseKEGG.

#DAVID functional analysis
david <- enrichDAVID(gene = gene,
                     idType = "ENTREZ_GENE_ID",
                     listType = "Gene",
                     annotation = "KEGG_PATHWAY",
                     david.user = "clusterProfiler@hku.hk")

#Universal enrichment analysis

#Functional analysis of NGS data
http://bioconductor.org/packages/release/bioc/html/ChIPseeker.html

#Visualization section
#BarPlot 
barplot(ggo, drop=TRUE, showCategory=12)
barplot(ego, showCategory=8)

#Dotplot
dotplot(ego)

#enrichMap
enrichMap(ego)

#cnetplot
cnetplot(ego, categorySize="pvalue", foldChange=geneList)

cnetplot(kk, categorySize="geneNum", foldChange=geneList)

# gseaplot
gseaplot(kk2, geneSetID = "hsa04145")

# plotGOgraph
plotGOgraph(ego)

#pathview from pathview package
#clusterProfiler users can also use pathview from the pathview9 to visualize KEGG pathway.
library("pathview")
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

#Biological theme comparison
#clusterProfiler was developed for biological theme comparison2, and it provides a function, compareCluster, to automatically calculate enriched functional categories of each gene clusters.
data(gcSample)
lapply(gcSample, head)
#The input for geneCluster parameter should be a named list of gene IDs. To speed up the compilation of this document, we set use_internal_data = TRUE.
ck <- compareCluster(geneCluster = gcSample, fun = "enrichKEGG", use_internal_data = TRUE)
head(summary(ck))

