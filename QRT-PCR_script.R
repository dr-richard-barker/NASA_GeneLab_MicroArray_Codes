#### This is my attempt at processing QRT-PCR data using R

### install packages and setup library's
## install.packages("xlsx")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("lattice")

## library(lattice) # for function 'barchart'
## library(xlsx)
library(dplyr)
library (ggplot2)

### This installs specific packages for QRT-PCR but i'm not sure how to use it...
source("https://bioconductor.org/biocLite.R")
biocLite("HTqPCR")
all.R.commands <- system.file("doc", "HTqPCR.Rnw", package = "HTqPCR")
Stangle(all.R.commands)
### this was the next command on the HTqPCR read me PDF but it didn't work https://www.bioconductor.org/packages/release/bioc/vignettes/HTqPCR/inst/doc/HTqPCR.pdf
## ls(package:HTqPCR)

### import data NOTE: replace т with T before importing
## mydata <- read.xlsx("/Users/richardbarker/Desktop/flg22_set_2_plate_1.xls", 1, startRow=7, header=TRUE, col.names=TRUE)

### import data NOTE: replace т with T before importing, replace space with _ and save as .csv ####
mydata <- read.csv ("/Users/richardbarker/Desktop/flg22_set_2_plate_1.csv", 1,)

### Have a quick look at the data
View(mydata)
str(mydata)

### count the number of rows so you can cut down the table to the correct size
nrow(mydata)
head(mydata,2)
names(mydata2) ## this tells you the head names

### remove the first rows that don't include data
##mydata_remove_first6<-mydata %>% slice (6:106)
##View(mydata_remove_first6)
##str(mydata_remove_first6)

### plot graphs

hist(mydata$CT.Mean)
hist(mydata$CT.SD)

ggplot( data = mydata, aes(y=CT.Mean,x=Target.Name), fill = "white") + geom_bar(stat="identity")
       
       
