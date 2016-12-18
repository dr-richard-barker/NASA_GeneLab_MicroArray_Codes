# Grass Grub data analyses
library(car)

GG <- read.csv("~/Desktop/GG BioassayRawData.csv")
str(GG)
GG["Block.."] <- as.factor(GG$Block..)


hist(GG$Corrected.Change.in.Root.Wet.Mass)
hist(GG$Larva_weight.change)


#### main model for change in root mass
GG1 <- aov(GG$Corrected.Change.in.Root.Wet.Mass ~ GG$Block.. + GG$larva_Initial.Wet + GG$Endophyte.Isolate + GG$Treatment + GG$Endophyte.Isolate*GG$Treatment + GG$Endophyte.Isolate*GG$Treatment*GG$larva_Initial.Wet)
plot(GG1)
summary(GG1) #none of the interactions with the covariate (larva initial wet weight) were significant, so they are taken out

GG2 <- aov(GG$Corrected.Change.in.Root.Wet.Mass ~ GG$Block.. + GG$larva_Initial.Wet + GG$Endophyte.Isolate + GG$Treatment + GG$Endophyte.Isolate*GG$Treatment)
plot(GG2)
summary(GG2)
coefficients(GG2)
outlierTest(GG2) # observation 124 is an outlier
GG = GG[-124,]
outlierTest(GG2) # 106 is an outlier
GG = GG[-106,]
outlierTest(GG2) # 31 is an outlier
GG = GG[-31,]
outlierTest(GG2) # 116 is an outlier
GG = GG[-116,]
outlierTest(GG2) # 117 is an outlier
GG = GG[-117,]
outlierTest(GG2) # 155 is an outlier
GG = GG[-155,]
outlierTest(GG2) # 126 is an outlier
GG = GG[-126,]
outlierTest(GG2) # 81 is an outlier
GG = GG[-81,]
outlierTest(GG2) # 146 is an outlier
GG = GG[-146,]
outlierTest(GG2) # 146 is an outlier
GG = GG[-146,]
outlierTest(GG2) # 16 is an outlier
GG = GG[-16,]
outlierTest(GG2) # 32 is an outlier
GG = GG[-32,]
outlierTest(GG2) # 39 is an outlier
GG = GG[-39,]
outlierTest(GG2) # 39 is an outlier
GG = GG[-38,]
outlierTest(GG2) # 102 is an outlier
GG = GG[-102,]
outlierTest(GG2) # 117 is an outlier
GG = GG[-117,]
outlierTest(GG2) # 87 is an outlier
GG = GG[-87,]
outlierTest(GG2) # 24 is an outlier
GG = GG[-24,]


# comparing models
anova(GG1, GG2) # the interactions with the covariate do not significantly affect the model, so the most parsimonius model is GG2

# contrasts
library(multcomp)
K <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,-1,0,0,0,0,0,0), 1)
t <- glht(GG2, linfct = K)
summary(t)
contrast1 <- read.csv("~/Downloads/contrast 1.csv")
C1 <- aov(contrast1$Corrected.Change.in.Root.Wet.Mass ~ contrast1$larva_Initial.Wet + contrast1$Treatment)
summary(C1)


K1 <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0,0,0,0,0,-1,0,0,0), 1)
t1 <- glht(GG2, linfct = K1)
summary(t1)
contrast2 <- read.csv("~/Downloads/contrast2.csv")
C2 <- aov(contrast2$Corrected.Change.in.Root.Wet.Mass ~ contrast2$larva_Initial.Wet + contrast2$Treatment)
summary(C2)



K2 <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0,-1,0,0,0,0,0,0,0,0), 1)
t2 <- glht(GG2, linfct = K2)
summary(t1)
contrast3 <- read.csv("~/Downloads/contrast3.csv")
C3 <- aov(contrast3$Corrected.Change.in.Root.Wet.Mass ~ contrast3$larva_Initial.Wet + contrast3$Treatment)
summary(C3)



contrast4 <- read.csv("~/Downloads/contrast4.csv")
C4 <- aov(contrast4$Corrected.Change.in.Root.Wet.Mass ~ contrast4$larva_Initial.Wet + contrast4$Treatment)
summary(C4)


#### main model for change in root mass
GG3 <- aov(GG$Larva_weight.change ~ GG$Block.. + GG$larva_Initial.Wet + GG$Endophyte.Isolate + GG$Treatment + GG$Endophyte.Isolate*GG$Treatment + GG$Endophyte.Isolate*GG$Treatment*GG$larva_Initial.Wet)
plot(GG3)
summary(GG3) #none of the interactions with the covariate (larva initial wet weight) were significant, so they are taken out

GG4 <- aov(GG$Larva_weight.change ~ GG$Block.. + GG$larva_Initial.Wet + GG$Endophyte.Isolate + GG$Treatment + GG$Endophyte.Isolate*GG$Treatment)
plot(GG4)
summary(GG4)
coefficients(GG4)

# comparing models
anova(GG3, GG4) # the interactions with the covariate do not significantly affect the model, so the most parsimonius model is GG2

outlierTest(GG3) # observation 89 is an outlier
GG = GG[-89,]
outlierTest(GG3) # observation 130 is an outlier
GG = GG[-130,]
outlierTest(GG3) # observation 116 is an outlier
GG = GG[-116,]
outlierTest(GG3) # observation 129 is an outlier
GG = GG[-129,]
outlierTest(GG3) # observation 85 is an outlier
GG = GG[-85,]
outlierTest(GG3) # observation 31 is an outlier
GG = GG[-31,]
outlierTest(GG3) # observation 57 is an outlier
GG = GG[-57,]
outlierTest(GG3) # observation 112 is an outlier
GG = GG[-112,]
outlierTest(GG3) # observation 91 is an outlier
GG = GG[-91,]
outlierTest(GG3) # observation 3 is an outlier
GG = GG[-3,]
outlierTest(GG3) # observation 3 is an outlier
GG = GG[-3,]

C1a <- aov(contrast1$Larva_weight.change ~ contrast1$larva_Initial.Wet + contrast1$Treatment)
summary(C1a)

C3a <- aov(contrast3$Larva_weight.change ~ contrast3$larva_Initial.Wet + contrast3$Treatment)
summary(C3a)
