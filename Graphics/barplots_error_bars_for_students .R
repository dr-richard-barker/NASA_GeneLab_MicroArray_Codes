setwd("/Users/RichardBarker/Dropbox")

# Before reading the table into R, do the following:
# 1. Remove quotations 
# 2. Do a global find and replace space for underscore.

library(reshape) # for the 'colsplit' function
library(lattice) # for function 'barchart'

# look at the files in your directory
list.files()

# read the appropriate file and look at the first lines

read.table("/Users/RichardBarker/DropboxYaz_Mike_skew", header=1) -> data
head(data)

# split the first column by underscores and name them. I dummy named the here
colsplit(data$Tag, split="_", names=c("Genotype", "Age",
                                      "Replicate"
                                )) -> Data

View(Data)

# the colsplit function got rid of the data columns, lets add them again with cbind
head(Data) 

#This attaches the original data/measurments to the new tag column structure and renameing the new system all_data
cbind(Data, data[2:7]) -> all_data

# not essential but you can view your data frame with these commands
View(all_data)
str(all_data)
head(all_data) 

# lets look at the groups
levels(all_data$Phytogel) 
levels(all_data$Sucrose)
levels(all_data$Genotype)

###################################
# lets make a 'barchart' with lattice
# I can't get it just right but you can play with it
# below that ones another one with error bars
###################################

#barchart(Total_Length~as.character(Fourth_column), data=all_data, groups=as.character(First_column), 
 #        scales=list(x=list(rot=50,cex=0.8)))

#barchart(Total_Length~as.character(First_column), data=all_data, groups=as.character(Fourth_column), 
 #        scales=list(x=list(rot=50,cex=0.8)))

################## barplot using the third column as group and with error bars

#all_data[all_data$Genotype %in% "Col-0",] -> group_1
#all_data[all_data$Genotype %in% "Cvi-0",] -> group_2
#all_data[all_data$Genotype %in% "Ler-0",] -> group_3
#all_data[all_data$Genotype %in% "WS-2",] -> group_4

#colnames(group_1)

#group_1[9] -> total_length_group_1
#group_2[9] -> total_length_group_2
#group_3[9] -> total_length_group_3
#group_4[9] -> total_length_group_4

#c(mean(total_length_group_1$Total_Length), mean(total_length_group_2$Total_Length),
  #mean(total_length_group_3$Total_Length),mean(total_length_group_4$Total_Length)) -> means

################################################################
#####Total length bar plot
################################################################

# the means have been separated so different media and ecotypes can be visusalised
means <- aggregate(Total.Length ~ Genotype, all_data, mean)
means1 <- means[1:10,] 
means2 <- means[11:20,]

#View(means)
#str(means)
#MeanTotalLength <- as.numeric(means$Total.Length)
MeanTotalLength1 <- as.numeric(means1$Total.Length)
MeanTotalLength2 <- as.numeric(means2$Total.Length)

#stDevs <- matrix(c(sd(total_length_group_1$Total_Length), sd(total_length_group_2$Total_Length), sd(total_length_group_3$Total_Length),sd(total_length_group_4$Total_Length)), 4)

#this calculates the standard error for the error bar
stDevs <- aggregate(Total.Length ~ Genotype, all_data, function(x) SE=sd(x)/sqrt(length(x)))
stDevs1 <- stDevs[1:10,]
stDevs2 <- stDevs[11:20,]
SDTotalLength1 <- as.numeric(stDevs1$Total.Length)
SDTotalLength2 <- as.numeric(stDevs2$Total.Length)

## Bar Plot with "OLIVE GREEN" bars
barplot(means1$Total.Length, axes=FALSE, axisnames=FALSE, ylim=c(0, 350),
              col=c("cadetblue1", "darkolivegreen4", "cadetblue2", "darkolivegreen4", "cadetblue3", "darkolivegreen4", "cadetblue4", "darkolivegreen4", "cadetblue", "darkolivegreen4"), main="Total length", xlab="Treatment", ylab="Total length (mm)", beside = true) -> barplot

axis(1, labels=c("M10.1", "M10.1(WT)", "M22.2", "M22.2(WT)", "M23.5", "M23.5(WT)", "M25.2", "M25.2(WT)", "M28.3", "M28.3(WT)"), at=barplot, cex=1)
axis(2, at=seq(0 , 500, by=50))
#legend("topright", legend=c("Col-0", "Cvi-0", "Ler-0", "WS-2"), col=c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4"), pch=15, bty = "n", cex= 1.5)
box()
segments(barplot, MeanTotalLength1 - SDTotalLength1, barplot, MeanTotalLength1 + SDTotalLength1, lwd=2)
segments(barplot - 0.1, MeanTotalLength1 - SDTotalLength1, barplot + 0.1, MeanTotalLength1 - SDTotalLength1, lwd=2)
segments(barplot - 0.1, MeanTotalLength1 + SDTotalLength1, barplot + 0.1, MeanTotalLength1 + SDTotalLength1, lwd=2)
# 15 by 9 inch pdf dimensions


barplot(means2$Total.Length, axes=FALSE, axisnames=FALSE, ylim=c(0, 350),
        col=c("aquamarine", "darkolivegreen4", "aquamarine1", "darkolivegreen4", "aquamarine2", "darkolivegreen4", "aquamarine3", "darkolivegreen4", "aquamarine4", "darkolivegreen4"), main="Total length", xlab="Treatment", ylab="Total length (mm)", beside = true) -> barplot

axis(1, labels=c("M29", "M29(WT)", "M31", "M31(WT)", "M33.5", "M33.5(WT)","M8.2", "M8.2(WT)", "M9.3.9", "M9.3.9(WT)"), at=barplot, cex=1)
axis(2, at=seq(0 , 500, by=50))
#legend("topright", legend=c("Col-0", "Cvi-0", "Ler-0", "WS-2"), col=c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4"), pch=15, bty = "n", cex= 1.5)
box()
segments(barplot, MeanTotalLength2 - SDTotalLength2, barplot, MeanTotalLength2 + SDTotalLength2, lwd=2)
segments(barplot - 0.1, MeanTotalLength2 - SDTotalLength2, barplot + 0.1, MeanTotalLength2 - SDTotalLength2, lwd=2)
segments(barplot - 0.1, MeanTotalLength2 + SDTotalLength2, barplot + 0.1, MeanTotalLength2 + SDTotalLength2, lwd=2)
######################



####### convex hull


# the means have been separated so different media and ecotypes can be visusalised
means <- aggregate(Convex.Hull ~ Genotype, all_data, mean)
means1 <- means[1:10,] 
means2 <- means[11:20,]

#View(means)
#str(means)
#MeanTotalLength <- as.numeric(means$Total.Length)
Meanhull1 <- as.numeric(means1$Convex.Hull)
Meanhull2 <- as.numeric(means2$Convex.Hull)

#stDevs <- matrix(c(sd(total_length_group_1$Total_Length), sd(total_length_group_2$Total_Length), sd(total_length_group_3$Total_Length),sd(total_length_group_4$Total_Length)), 4)

#this calculates the standard error for the error bar
stDevs <- aggregate(Convex.Hull ~ Genotype, all_data, function(x) SE=sd(x)/sqrt(length(x)))
stDevs1 <- stDevs[1:10,]
stDevs2 <- stDevs[11:20,]
SDhull1 <- as.numeric(stDevs1$Convex.Hull)
SDhull2 <- as.numeric(stDevs2$Convex.Hull)

## Bar Plot with "OLIVE GREEN" bars
barplot(means1$Convex.Hull, axes=FALSE, axisnames=FALSE, ylim=c(0, 3500),
        col=c("cadetblue1", "darkolivegreen4", "cadetblue2", "darkolivegreen4", "cadetblue3", "darkolivegreen4", "cadetblue4", "darkolivegreen4", "cadetblue", "darkolivegreen4"), main="Convex Hull", xlab="Treatment", ylab="Area (mm2)", beside = true) -> barplot

axis(1, labels=c("M10.1", "M10.1(WT)", "M22.2", "M22.2(WT)", "M23.5", "M23.5(WT)", "M25.2", "M25.2(WT)", "M28.3", "M28.3(WT)"), at=barplot, cex=1)
axis(2, at=seq(0 , 3500, by=100))
#legend("topright", legend=c("Col-0", "Cvi-0", "Ler-0", "WS-2"), col=c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4"), pch=15, bty = "n", cex= 1.5)
box()
segments(barplot, Meanhull1 - SDhull1, barplot, Meanhull1 + SDhull1, lwd=2)
segments(barplot - 0.1, Meanhull1 - SDhull1, barplot + 0.1, Meanhull1 - SDhull1, lwd=2)
segments(barplot - 0.1, Meanhull1 + SDhull1, barplot + 0.1, Meanhull1 + SDhull1, lwd=2)
# 15 by 9 inch pdf dimensions


barplot(means2$Convex.Hull, axes=FALSE, axisnames=FALSE, ylim=c(0, 3500),
        col=c("cadetblue1", "darkolivegreen4", "cadetblue2", "darkolivegreen4", "cadetblue3", "darkolivegreen4", "cadetblue4", "darkolivegreen4", "cadetblue", "darkolivegreen4"), main="Convex Hull", xlab="Treatment", ylab="Area (mm2)", beside = true) -> barplot

axis(1, labels=c("M29", "M29(WT)", "M31", "M31(WT)", "M33.5", "M33.5(WT)","M8.2", "M8.2(WT)", "M9.3.9", "M9.3.9(WT)"), at=barplot, cex=1)
axis(2, at=seq(0 , 3500, by=100))
#legend("topright", legend=c("Col-0", "Cvi-0", "Ler-0", "WS-2"), col=c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4"), pch=15, bty = "n", cex= 1.5)
box()
segments(barplot, Meanhull2 - SDhull2, barplot, Meanhull2 + SDhull2, lwd=2)
segments(barplot - 0.1, Meanhull2 - SDhull2, barplot + 0.1, Meanhull2 - SDhull2, lwd=2)
segments(barplot - 0.1, Meanhull2 + SDhull2, barplot + 0.1, Meanhull2 + SDhull2, lwd=2)




