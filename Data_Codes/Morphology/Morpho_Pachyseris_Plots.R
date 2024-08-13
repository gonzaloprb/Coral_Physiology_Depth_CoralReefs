# Script for quick plot visualisations and to generate databases for All_features Script "Script0_Opening.R", and further analysis

# Pachyseris speciosa
# Graphs and database
require (ggplot2)
require (dplyr)
require (plyr) 
require (ggstatsplot)
require (tidyverse)

rm (list = ls ())

morpho_pachy <- read.csv ("Data_Codes/Morphology/Morpho_Pachy_Final.csv", header = T, dec = ",", sep = ";")
### Variables = "Septa.Distance", "Width.Valley","Height.Septa"

colnames (morpho_pachy) <- c("Island_Site", "Depth","ID_LAB","Septa.Distance", "Width.Valley","Height.Septa")

str (morpho_pachy)
quan.var <- c("Septa.Distance", "Width.Valley","Height.Septa")
morpho_pachy [quan.var] <- sapply (morpho_pachy [quan.var], as.numeric)

summary(morpho_pachy)

# There are some NAs because I could not do the 10 replicate reading.


###### STEP 1 -> mean and sd of the variables for EACH SAMPLE ----
resume_pachy <-ddply(morpho_pachy,~ID_LAB+Depth+Island_Site,function(x){c(Septa.Distance=mean(na.omit(x$Septa.Distance)), Width.Valley=mean(na.omit(x$Width.Valley)), Height.Septa=mean(na.omit(x$Height.Septa)), septa_sd=sd(na.rm = FALSE,(x$Septa.Distance)), height_sd=sd(na.rm = FALSE,(x$Height.Septa)), valley_sd=sd(na.rm = FALSE,(x$Width.Valley)))})
summary(resume_pachy)


###### STEP 2 -> mean and sd of the variables for EACH DEPTH AND ISLAND ----
resume_depth_island_pachy <-ddply(morpho_pachy,~Depth+Island_Site,function(x){c(Septa.Distance=mean(na.omit(x$Septa.Distance)), Width.Valley=mean(na.omit(x$Width.Valley)), Height.Septa=mean(na.omit(x$Height.Septa)), septa_sd=sd(na.omit(x$Septa.Distance)), height_sd=sd(na.omit(x$Height.Septa)), valley_sd=sd(na.rm = T,(x$Width.Valley)))})
summary(resume_depth_island_pachy)

##### STEP 3 -> Plot the three variables
resume_pachy$Depth <- as.factor (resume_pachy$Depth)
resume_depth_island_pachy$Depth <- as.factor (resume_depth_island_pachy$Depth)

#### 1) Plot of SEPTA SPACING ----
ggplot(resume_pachy, aes(x=Depth, y=Septa.Distance)) +
  geom_boxplot() + theme_bw() + stat_summary(fun=mean, geom="point", shape=18, color="red", size=2, position=position_dodge(width=0.75)) +
  facet_grid(. ~ Island_Site, scales = "free_x") + 
  ggtitle("Septa spacing") + ylab("septa spacing (mm)") +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.text.x = element_text(size=10, colour="black"),
        axis.text.y = element_text(size=9, face="bold", colour="black"),
        axis.title.y = element_text(size=11, face="bold", colour="black"),
        axis.title.x = element_text(size=11, face="bold", colour="black")) 



#### 2) VALLEY WIDTH ----
ggplot(resume_pachy, aes(x=Depth, y=Width.Valley)) +
  geom_boxplot() + theme_bw() + stat_summary(fun=mean, geom="point", shape=18, color="red", size=2, position=position_dodge(width=0.75)) +
  facet_grid(. ~ Island_Site, scales = "free_x") + 
  ggtitle("Valley width") + ylab("valley width (mm)") +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.text.x = element_text(size=10, colour="black"),
        axis.text.y = element_text(size=9, face="bold", colour="black"),
        axis.title.y = element_text(size=11, face="bold", colour="black"),
        axis.title.x = element_text(size=11, face="bold", colour="black")) 

 

#### 3) HEIGHT SEPTA ----
ggplot(resume_pachy, aes(x=Depth, y=Height.Septa)) +
  geom_boxplot() + theme_bw() + stat_summary(fun=mean, geom="point", shape=18, color="red", size=2, position=position_dodge(width=0.75)) +
  facet_grid(. ~ Island_Site, scales = "free_x") + 
  ggtitle("Height.Septa") + ylab("Height.Septa (mm)") +
  theme(plot.title = element_text(hjust=0.5, size=14, face="bold"),
        axis.text.x = element_text(size=10, colour="black"),
        axis.text.y = element_text(size=9, face="bold", colour="black"),
        axis.title.y = element_text(size=11, face="bold", colour="black"),
        axis.title.x = element_text(size=11, face="bold", colour="black")) 





# Generate output databases:
Pachyseris_morpho_features <- resume_pachy %>%
  group_by(ID_LAB, Depth, Island_Site) %>% 
  summarise_each(funs(mean))
# Pachyseris_morpho_features <- aggregate (.~ ID + Depth + Site, resume_pachy, mean)

write_csv(Pachyseris_morpho_features, "Data_Codes/Morphology/Pachyseris_Morphology.csv")


# Make a super summary for modelling
resume_Pachyseris_Morphology <-ddply(morpho_pachy,~Depth+Island_Site,function(x){c(Septa.Distance=mean(na.omit(x$Septa.Distance)), Width.Valley=mean(na.omit(x$Width.Valley)), Height.Septa=mean(na.omit(x$Height.Septa)), Septa.Distance_sd=sd(na.rm = FALSE,(x$Septa.Distance)), Width.Valley_sd=sd(na.rm = FALSE,(x$Width.Valley)), Height.Septa_sd=sd(na.rm = FALSE,(x$Height.Septa)))})
summary(resume_Pachyseris_Morphology)


