# WeGas Analyses

# upload data
weGAS <- read.csv("~/Desktop/weGAS.csv")
View(weGAS)
str(weGAS)

Shoot <- weGAS[1:511,]
View(Shoot)

Root <- weGAS[512:1238,]
View(Root)

# viewing distributions
hist(weGAS$Emergence.Angle)
hist(1/(weGAS$Emergence.Angle))
hist(weGAS$Tip.Angle)
hist(weGAS$Total.Length)
hist(weGAS$Total.Primary.Angle)
hist(weGAS$Tortuosity)

#hist(Shoot$Emergence.Angle)
#hist(1/(Shoot$Emergence.Angle))
hist(Shoot$Tip.Angle)
hist(Shoot$Total.Length) #almost normal
hist((Shoot$Total.Length)) 
#hist(Shoot$Total.Primary.Angle)
#hist(Shoot$Tortuosity)

#hist(Root$Emergence.Angle)
#hist(1/(Root$Emergence.Angle))
hist(Root$Tip.Angle) ##IS THIS ACTUALLY TOTAL LENGTH?????
hist(Root$Total.Length) #normal IS THIS ACTUALLY TIP ANGLE
RTPABIN <- hist(Root$Total.Primary.Angle, breaks = c(-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110)) #normal
RTPABIN$counts
plot(RTPABIN)
#hist(Root$Tortuosity) #normal

# simple linear models (first round)
library("MASS")
library("stats")

STL <- lm(Total.Length ~ X02conc*Genotype, data=Shoot)
boxcox(STL) #0 or 0.5 log or sqrt
STL1 <- lm(log(Total.Length) ~ X02conc*Genotype, data=Shoot)
summary(STL1)

drop1(STL1,~.,test="F")
STL1A <- aov(log(Total.Length) ~ X02conc*Genotype, data=Shoot)
TukeyHSD(STL1A)

#STL2 <- lm(sqrt(Total.Length) ~ X02conc*Genotype, data=Shoot)
#summary(STL2)
#drop1(STL2,~.,test="F")

RTA <- lm(Tip.Angle ~ X02conc*Genotype, data=Root)
boxcox(RTA) #0 log
RTA1 <- lm(log(Tip.Angle) ~ X02conc*Genotype, data=Root)
drop1(RTA1,~.,test="F")
RTA1A <- aov(log(Tip.Angle) ~ X02conc*Genotype, data=Root)
TukeyHSD(RTA1A)

RTL <- lm(Total.Length ~ X02conc*Genotype, data=Root)
#boxcox(RTL) #0 
drop1(RTL,~.,test="F")
RTLA <- aov(Total.Length ~ X02conc*Genotype, data=Root)
TukeyHSD(RTLA)

RTPA <- lm(Total.Primary.Angle ~ X02conc*Genotype, data=Root)
#boxcox(RTPA) #0 
drop1(RTPA,~.,test="F")
RTPAA <- aov(Total.Primary.Angle ~ X02conc*Genotype, data=Root)
TukeyHSD(RTPAA)


#histogram counts

str(Root)
ColRTPABIN <- Root[ which(Root$Genotype =='Col-0'),]
CviRTPABIN <- Root[ which(Root$Genotype =='Cvi-0'),]
LerRTPABIN <- Root[ which(Root$Genotype =='Ler-0'),]
WSRTPABIN <- Root[ which(Root$Genotype =='WS-2'),]

ColRTPABINH <- hist(ColRTPABIN$Total.Primary.Angle, breaks = c(-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110)) #normal
ColRTPABINH$counts

CviRTPABINH <- hist(CviRTPABIN$Total.Primary.Angle, breaks = c(-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110)) #normal
CviRTPABINH$counts

LerRTPABINH <- hist(LerRTPABIN$Total.Primary.Angle, breaks = c(-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110)) #normal
LerRTPABINH$counts

WSRTPABINH <- hist(WSRTPABIN$Total.Primary.Angle, breaks = c(-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110)) #normal
WSRTPABINH$counts

