## comments are useful 
## like basic arithmetic
4+6
4*6
4^6 ### 4 to the power of 6
4/6
sqrt(46)
log(46)
log10(46)
log(46, 10)

## this sign <- causes creates a new object
age <- 2016 - 1969
## then you can use the object for further claculation
sqrt_age <- sqrt(age)
sqrt_age 
sqrt(age)

##Vectors
## Note: don't use C as a variable name as it also functions as a "combine" or concatinate functions
weights <- c(30, 100, 4000, 8000)
length(weights)
animals <-c("rat", "cats")
## you can add on to a vetor using <- c
animals <-c ("mouse", animals)

## what type or class of information is contained within a vector
class(animals)
class(weights)

## animals = text
## weights = numeric
## note: you can't mix vector types like this...
x<- c (1,2,3, 'a')
y<- c (1,2,3, TRUE)
z<- c ('a', TRUE, 'b', 'c')
tricky <- c(1, '2', 3, 4)

## this tells you about the data structure
str(weights)
str(animals) 

## this tells you what variables you have...
ls() ## list
rm() ## remove

## indexing / subsetting
weights
weights[1]
weights[2]
weights[c(2,3)]

animals [c(1,2,3,1,2,4)]
animals <- c("mouse", "rat", "dog", "cat")

## calculations to the whole vector
-weights
weights*2.2/1000
sqrt(weights)
log10(weights)
weights
save_these <- c(TRUE,FALSE, FALSE, FALSE, FALSE)
weights[save_these]

weights [weights <50]
animals [weights <1000]
animals [weights >1000]
animals [weights >50 & weights <5000]
weights[animals %in% c("rats", "cats")]

## missing values
heights <- c(2, 4, 4, NA, 6)
## note: missing values cause basic calculations to fail
mean(heights)
max(heights)
## thus
mean(heights, na.rm=TRUE)
max(heights, na.rm=TRUE)
na.omit(heights)

##Note: if you put NA in "NA" then it creates a different type of vector
v <- c(2, 4, 4, "NA", 6)
v
mean (v)
mean (v, na.rm=TRUE)
na.omit (v)
# make sure you don't have "NA" in your codes?...

## Note: ! means do the opposit

## load the data from the url "http://kbroman.org/datacarp/portal_data_joined.csv"
url <- "http://kbroman.org/datacarp/portal_data_joined.csv"
filename <- "portal_data_joined.csv"
download.file(url, filename)
surveys <- read.csv ("portal_data_joined.csv")

## or 
download.file("http://kbroman.org/datacarp/portal_data_joined.csv", "portal_data_joined.csv")
surveys <- read.csv ("portal_data_joined.csv")E
## note: read.csv makes every column as a factor

## this makes strings of numbers instead of factors
surveys <- read.csv ("portal_data_joined.csv" stringsAsFactors=HELLNO) 

## new data... start with STR and head
head(surveys)
head(surveys, 2)
tail (surveys)
tail (surveys, 3)
ncol(surveys) # columns
nrow(surveys) # rows
summary(surveys) # this might be a really interesting summary! 
names(surveys) # column names
rownames(surveys) # row names
dim(surveys)
sex <- factor( c("male", "female"))
sex
class(sex)
levels (sex)
nlevels(sex)
sex_char <- as.character(sex)
sex_char

numbers <- factor(c("1","5","8","4"))
levels (numbers)
as.numeric(as.character(numbers))
expt <- c("treat1", "treat2", "treat1", "treat1", "treat3", "treat1", "control", "treat1","treat2", "treat3")
## note: in cases like this just write t1, t2, t3, control

expt <- factor (expt)
table(expt)
levels(expt)

## Indexing of data frames
surveys[1,1]
surveys[4,7]

###These all do the same thing...
sex <- surveys[,"sex"]
sex <- surveys$sex
sex

## sequences, : and seq
5:10
1:10
seq (5, 10, by=2)
seq (5, 10, by=0.5)
seq (5, 10, length.out=26)

## slices of a data frame
surveys [10:20, 2:4]
surveys [10:20, 2:5]

male_weights <-surveys[surveys$sex =="M", "weights"]
male_weights <-surveys$weight[surveys$sex=="M"]
male_weights
mean(male_weights, na.rm=TRUE) ### mean male weight

###Useful packages 
install.packages("dplyr")
install.packages("ggpolot2")
library(dplyr)

Dplyr aims to provide a function for each basic verb of data manipulation:
  
##  https://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html  
filter() (and slice())
arrange()
select() (and rename())
distinct()
mutate() (and transmute())
summarise()
sample_n() and sample_frac()


selected_col <- select(surveys, species_id, sex)
filter (surveys, sex=="M")

###pipe
== ## means compare
library("dplyr")

#dplyr   SQL
#select  SELECT
#filter  WHERE
#mutate 
#group_by GROUP BY
#summarize
#arrange  ORDER BY  
#tally

#%>% (cmd+shift+M)
#AND
surveys %>% filter(sex=="M",weight < 10) %>%  select(species_id,sex) 
#OR
surveys %>% filter(sex=="M" | weight < 10) %>%  select(species_id,sex)

#Challenge answer
surveys %>% filter(year < 1995) %>% select(year, sex, weight)
surveys %>% mutate(sqrt_weight=sqrt(weight)) %>% head
surveys %>% filter (!is.na(weight)) %>% mutate(sqrt_weight=sqrt(weight), weight=weight*2.2)

###These both failed
head(select(filter(surveys, weights<6), species_id, weights))
filtered <- filter(surveys, weights<6)

## Failed...
challenge_results <- surveys %>% mutate(hindfoot_sqrt = sqrt(hindfoot_length)) %>% filter(!is.na(hindfoot_sqrt), hindfoot_sqrt < 3) %>% select (species_id, hindfoot_sqrt)

## arrange 
challenge_results %>% arrange (hindfoot_sqrt)
challenge_results %>% arrange (desc(hindfoot_sqrt),species_id)
## desc = decending

## group by
surveys %>% group_by(sex) %>% tally
surveys %>% group_by(sex, year) %>% tally

surveys %>% group_by(sex) %>% summarize(mean_weight=mean(weight, na.rm=TRUE))
##alternatively
surveys %>% group_by(sex) %>% filter(!is.na(weights)) %>% summarize(mean_weights=mean(weights)) ,min_weight=min(weights),max_weight=max(weights))

## some data cleaning
surveys_complete <- surveys %>% filter(!is.na(hindfoot_length))%>% filter (!is.na(weight)) %>% filter (sex != "")
# count species
count_species <- surveys_complete %>% group_by (species_id) %>% tally
# species with >= 10 count
frequent_species <- count_species %>% filter (n>=10)

#save just rows with species_id count >= 10
surveys_complete <-surveys_complete %>% filter (species_id %>% frequent_species$speciess_id)

###ggplot graph creation
##Note: this changes the graphs from a grey theme to a black and white one...
XXXYYYZZZ code + theme_bw()
### aes = asethetics = tells the software which column or vector should be on the X or Y axis... 
## Note the + operator allows you to add extra functions to the graph
surveys_complete %>% ggplot(aes(x=weight, y=hindfoot_length)) + geom_point ()
## this can be brocken down...
p1 <-surveys_complete %>% ggplot(aes(x=weight, y=hindfoot_length))
p2 <- p1 + geom_point() 
p2
p2 + scale_x_log10()
p2 + scale_x_sqrt()

p3 <-surveys_complete %>% ggplot(aes(x=species_id, y=hindfoot_length))
p4 <- p3 + geom_point()
p4

###Contains filter and error
surveys_complete %>% filter (species_id=="DM") %>%
  ggplot (aes((x=weight, y=hindfoot_length)) 
          + geom_point()
   
### other aesthetics : size, color, alpha, shape
genom_point(alpha=0.1, color="slateblue", size=20)

## this assigns colors to specific species
genom_point(aes(color=species_id))
 
## This one works                                                        
surveys_complete %>% 
filter(species_id %in% c("DM", "DS", "DO")) %>%
ggplot (aes(x=weight, y=hindfoot_length)) +
geom_point(aes(color=species_id, shape=species_id))                                                           

## this one works...!!
surveys_complete %>% group_by(year) %>% tally %>%
               ggplot (aes(x=year, y=n)) +
               geom_line(color="lightblue") +
               geom_point(color="violetred", size = 2)  
                                                                                                                            

### 
surveys_complete %>% filter (species_id %in% c("DM", "DS")) %>%
            group_by(year, species_id) %>%
            tally %>%
            ggplot (aes(x=year, y=n)) +
            geom_point(aes(color=species_id)) +
            geom_line(aes(color=species_id))

## histogram
surveys_complete %>% ggplot (aes(x=weight)) +
  geom_histogram(bins=100) 

### box plot
surveys_complete %>% 
  ggplot (aes(x=species_id, y=weight)) +
  geom_boxplot () 

### box plot with dots on top.
surveys_complete %>% 
  ggplot (aes(x=species_id, y=weight)) +
  geom_boxplot () + geom_jitter(col="Orchid", alpha=0.1)

#### faceting, this makes each "species" its own graphical pannel
surveys_complete %>%
  group_by(year, species_id) %>%
  tally %>% ggplot(aes(x=year, y=n)) +
  geom_line() + facet_wrap (~ species_id)

surveys_complete %>%
  ggplot (aes(x=weight, y=hindfoot_length)) +
  geom_point (aes(color=species_id))

surveys_complete %>%
  ggplot (aes(x=weight, y=hindfoot_length)) +
  geom_point (aes(color=species_id)) + facet_wrap(~ year)

surveys_complete %>% filter (species_id %in% c("DM", "DS")) %>%
  ggplot (aes(x=weight, y=hindfoot_length)) +
  geom_point (aes(color=species_id)) + facet_wrap(~ species_id)

surveys_complete %>% filter (species_id %in% c("DM", "DS")) %>%
  ggplot (aes(x=weight, y=hindfoot_length)) +
  geom_point (aes(color=species_id)) + facet_wrap(species_id ~year)

#### R mark down = Literate programming = Go to File -> new file 
