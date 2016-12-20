#################################################################
#################################################################
## HEATMAPS
##
## This script will ask for:
##    - a file containing at least two rows (names [must be first row] 
##      and Genes) and x columns. Each column must contain name for heatmap
##      and a commaseparated list of genes 
##    - a folder containing the output of a CuffDiff run
##    - a path to a folder in which the heatmaps will be saved
##
## It will then do the following for each list of genes (i.e. each column):
##    - convert the log changes to fold changes
##    - create heatmap
##    - save heatmap in specified folder
##
#################################################################
#################################################################
## Loading packages and define functions
##

## Load packages
library(cummeRbund, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(ggdendro, quietly = TRUE)
library(RGraphics, quietly = TRUE)
library(gridExtra, quietly = TRUE)
library(scales, quietly = TRUE)


source('/Users/richardbarker/Google Drive/genelab_BRIC19/heatmaps_functions.R')

#################################################################
#################################################################
## NOTE: THIS PART IS INTERACTIVE!!!!!!!!
## 
## Choose files and folders
##

cat('Choose file with lists:\n')
lists.file <- file.choose()
## '/Users/richardbarker/Google\ Drive/genelab_BRIC19/lists_for_heatmaps.txt'
cat(paste(lists.file, '\n', sep = ' '))

cat('Choose folder containing CuffDiff output:\n')
cuffdiff_out <- choose.dir()
cat(paste(cuffdiff_out, '\n', sep = ' '))

cat('Choose folder to save heatmaps to:\n')
path <- choose.dir()
cat(paste(path, '\n', sep = ' '))

##
##
#################################################################
#################################################################
## Set variables for later use
##

cat('\nDefault values: 
logbase = 2.
pseudo number = 1.
Colors chosen are red, white and yellow.
Non-significant fold changes will be colored grey95.
Dendogram will be included.
Line size is 1.
Database will NOT be rebuild.'
)
all.default <- readline('Use all default? y/n ')
if(all.default == 'y'){
  logbase <- 2
  ps <- 1
  cols <- c('red','white','yellow')
  check.for.significance <- 'y'
  na.col <- 'grey95'
  Dend <- 'y'
  LS <- 1
  rebuild <- 'n'
  p.value.cutoff <- 0.05
  
} else {
  message('To use default values, hit enter!')
  
  ## Choose logarithmic base for log transformation of fold changes 
  ## (default: 2):
  lb <- as.numeric(readline('What logbase would you like to use? (default: 2) '))
  logbase <- ifelse(is.na(lb),
                    2,
                    lb)
  cat(paste('You chose', logbase, '\n', sep = ' '))
  
  ## Choose pseudo value (default: 1):
  Ps <- as.numeric(readline('What pseudo value would you like to use? (default: 1) '))
  ps <- ifelse(is.na(Ps),
               1,
               Ps)
  cat(paste('You chose', ps, '\n', sep = ' '))
  
  ## Choose colors
  COL <- readline('What colors would you like to use? 
If white is used, at least three colors must be provided, with white not being the first nor last. 
White will then be used to indicate zero change. 
(Default: red, white, yellow) ')
  
  cols <- if(nchar(COL) < 1){
    c('red','white','yellow')
  } else {
    unlist(strsplit(COL, split = ',', fixed = FALSE))
  }
  cat(paste('You chose', paste(cols, collapse = ', '), '\n', sep = ' '))
  
  ## Color all non-significant fold changes white?
  check.for.significance <- readline('Color-code non-significant values? y/n (default: y) \n')
  if(nchar(check.for.significance) < 1){
    check.for.significance <- 'y'
  }
  
  if(check.for.significance == 'y'){
    na.col <- readline('What color do you want for non-significant values? (default: grey95) ')
    if(nchar(na.col) < 1){
      na.col <- 'grey95'
    }
    p.value.cutoff <- as.numeric(readline('p-value cutoff: (default: 0.05) '))
    if(is.na(p.value.cutoff)){
      p.value.cutoff <- 0.05
    }
    
    cat(paste('Non-significant values will be colored', na.col, '\n',sep = ' '))
  } else {
    cat('Non-significant values will not be color-coded. \n')
  }
  
  
  ## Choose border line size
  Ls <- readline('Choose line size for black borders: (default: 1) ')
  LS <- if(nchar(Ls) < 1){
    1
  } else { as.numeric(Ls) }
  cat(paste('You chose', LS, '\n', sep = ' '))
  
  
  ## Dendro/no dendro
  Dend <- readline('Do you want to include a dendrogram? y/n (default: y) ')
  if(nchar(Dend) < 1){
    Dend <- 'y'
  }
  
  Dend <- c('TRUE', 'FALSE')[match(Dend, c('y','n'))]
  
  if(Dend){
    cat('A dendrogram will be included. \n')
  } else {
    cat('No dendrogram will be included. \n')
  }
  
  
  #################################################################
  #################################################################
  ## Read CuffLinks and file containing lists before starting the 
  ## loop
  
  ## Should the database be rebuild?
  rebuild <- readline('Do you want to rebuild the database? y/n (default: n) ')
  if(nchar(rebuild) < 1){
    rebuild <- 'n'
  }
}

rebuild <- c(TRUE, FALSE)[match(rebuild, c('y','n'))]

## Read cuff links and create fpkm matrix
cuf <- readCufflinks(cuffdiff_out, 
                     rebuild = rebuild)
fpkmMat <- fpkmMatrix(genes(cuf))

## Load file containing lists
lists <- read.table(file = lists.file,
                    header = TRUE,
                    sep = '\t',
                    row.names = 1)


#################################################################
#################################################################
## LOOP
##

## Start loop over all lists
for (no in 1:length(lists)){
  
  message(paste('Number', no, 'out of', length(lists),
                sep = ' '))
  
  ## Select current list
  cur.list <- lists[no]
  
  ## Separate GO term and name from current list
  cur.name <- unlist(strsplit(names(cur.list), 
                              split = '.', 
                              fixed = TRUE))
  
  go.term <- paste(cur.name[1:2], collapse = ':')  
  name <- cur.name[-c(1:2)]
  
  ## Get p-value from current list
  p.value <- cur.list['PValue',]
  
  sigGenes <- as.character(cur.list['Genes',])
  sigGenes <- unlist(strsplit(sigGenes, split = ', ', fixed = TRUE))
  
  genesDiff <- getGenes(cuf, sigGenes)@diff
  
  fpkm.sigGenes <- fpkmMat[sigGenes,]	
  
  dbs <- list()
  dbs.na <- list()
  
  name.matrix <- matrix(names(fpkmMat), ncol = 2, byrow = TRUE)
  
  for(i in 1:nrow(name.matrix)){
    dbs[[i]] <- dbs.na[[i]] <- logFC(fpkm.sigGenes,
                                     mutants = name.matrix[i,1],
                                     WT = name.matrix[i,2],
                                     logBase = logbase,
                                     pseudo = ps)
    
    if(check.for.significance == 'y'){
      for(gene.id in sigGenes){
        print(subset(genesDiff, gene_id == gene.id &
                       sample_1 == name.matrix[i,1] & 
                       sample_2 == name.matrix[i,2])$significant)
        if(subset(genesDiff, gene_id == gene.id &
                  sample_1 == name.matrix[i,1] & 
                  sample_2 == name.matrix[i,2])$p_value > p.value.cutoff){
          
          dbs.na[[i]][gene.id] <- NA
        }
      }
    }
    
    names(dbs)[i] <- names(dbs.na)[i] <- paste(name.matrix[i,], collapse = '/')
    
  }
  
  DB <- do.call(cbind, dbs)
  DB.na <- do.call(cbind, dbs.na) 
  
  Filename <- paste(paste(name, collapse = '_'),
                    '.pdf',
                    sep = '')
  
  pdf(file = paste(path, Filename, sep = ''),
      height = unit(18,'cm'),
      width = unit(12,'cm'))
  HEATMAPS(db = DB,
           db.na = DB.na,
           cols = cols,
           na.col = na.col,
           go.term = go.term,
           name = name,
           p.value = p.value,
           line.size = LS,
           p.value.cutoff = p.value.cutoff)
  dev.off()
}
