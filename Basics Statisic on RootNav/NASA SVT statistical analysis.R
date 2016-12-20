# NASA SVT Data Analyses

Root <- read.csv("~/Desktop/Root measurements_KSC_SVT_processed.csv")
str(Root)

# Checking distributions
hist(Root$Total.Primary.Angle) #Normal
hist(Root$Total.Length) #Normal

# Anova
library(car)

## Skew angle vs Genotype
angle <- lm(Root$Total.Primary.Angle ~ Root$Position + Root$Genotype)
anova(angle)
#Anova(angle, type=c("III"))
summary(angle)

## Skew angle
angle_geno <- lm(Root$Total.Primary.Angle ~ Root$Genotype)
anova(angle_geno)
summary(angle_geno)

## Root Length vs Possition
length_posi <- aov(Root$Total.Length ~ Root$Position + Root$Genotype)
anova(length_posi)
summary(length_posi)

## Root Length vs Genotype
length_geno <- lm(Root$Total.Length ~ Root$Genotype)
anova(length_geno)
summary(length_geno)

### Convex Hull vs position
curvature_posi <- lm(Root$Convex.Hull ~ Root$Position + Root$Genotype)
anova(curvature_posi)
#Anova(curvature, type=c("III"))
summary(curvature_posi)

##### Convex Hull vs genotype
curvature_geno <- lm(Root$Convex.Hull ~ Root$Genotype)
anova(curvature_geno)
#Anova(curvature, type=c("III"))
summary(curvature_geno)
