# Script for quick plot visualisations and to generate databases for All_features Script "Script0_Opening.R", and further analysis

# Pocillopora verrucosa
# Graphs and database
require (ggplot2)
require (dplyr)
require (plyr)
require (ggstatsplot)
require (tidyverse)
require (ggstatsplot) 
 

morpho_poc <- read.csv ("Data_Codes/Morphology/Morpho_Pocillopora_Final.csv", header = T, dec = ",", sep = ";")

colnames (morpho_poc) <- c("Island_Site", "Depth","ID_LAB","Size", "Distance","Part")


str (morpho_poc)
quan.var <- c( "Size", "Distance")
morpho_poc [quan.var] <- sapply (morpho_poc [quan.var], as.numeric)

# Here not the problem of the 62m. Already considered as 60
summary (as.factor(morpho_poc$Depth))

# summary(morpho_poc)



###### STEP 1 -> mean and sd of the variables for EACH SAMPLE ----
resume_poc <-ddply(morpho_poc,~ID_LAB+Depth+Island_Site+Part,function(x){c(Size=mean(na.omit(x$Size)), Distance=mean(na.omit(x$Distance)), Size_sd=sd(na.omit(x$Size)), Distance_sd=sd(na.omit(x$Distance)))})
summary(resume_poc)


##### STEP 2 -> mean and sd of the variables for EACH DEPTH AND ISLAND ----
resume_depth_island_poc <-ddply(morpho_poc,~Depth+Island_Site+Part,function(x){c(Size=mean(na.omit(x$Size)), Distance=mean(na.omit(x$Distance)), Size_sd=sd(na.omit(x$Size)), Distance_sd=sd(na.omit(x$Distance)))})
summary(resume_depth_island_poc)


##### STEP 3 -> PLOT distance and size of corallites  ----
resume_poc$Depth <- as.factor (resume_poc$Depth)
resume_depth_island_poc$Depth <- as.factor (resume_depth_island_poc$Depth)



#### 1) Plot of distance


## a) without distinction between basal and apical regions. 
# Considering boxplot with all sd from ID
ggplot(resume_poc, aes(x=Depth, y=Distance)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Corallite spacing") + ylab ("Distance (mm)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# Considering boxplot with no sd of ID, just from Island & Depth
ggplot(resume_depth_island_poc, aes(x=Depth, y=Distance)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) +
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2, position=position_dodge(width=0.75)) + 
  theme_bw() + ggtitle("Corallite spacing") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

# Barplot
ggplot(resume_depth_island_poc, aes(x = Depth, y = Distance, fill=Part)) + 
  geom_bar(stat="identity",position=position_dodge()) +
  facet_grid(. ~ Island_Site) +
  scale_fill_manual(values=c("gray55","gray5"))+
  theme_bw() +
  geom_errorbar(aes(ymin=Distance-Distance_sd,ymax=Distance+Distance_sd), width=.4, position=position_dodge(0.9))


## b) with distinction between basal and apical regions

ggplot(resume_poc, aes(x=Depth, y=Distance, color=Part)) +
  geom_boxplot() + theme_bw() + facet_grid(Part ~ Island_Site) +
  ggtitle("Corallite spacing") + ylab("Distance (mm)") + xlab ("Depth (m)")
  theme(plot.title = element_text(hjust=0.5, size=10, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))


  
ggplot(resume_depth_island_poc, aes(x=Depth, y=Distance, color=Part)) +
  geom_boxplot() + theme_bw() + facet_grid(Part ~ Island_Site) +
  ggtitle("Corallite spacing") + ylab("Distance (mm)") + xlab ("Depth (m)")
  theme(plot.title = element_text(hjust=0.5, size=10, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))



### 2) SIZE of corallites ----

## c) without distinction between basal and apical regions

# Considering boxplot with all sd from ID
ggplot(resume_poc, aes(x=Depth, y=Size)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Corallite size") + ylab ("Size (mm)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

  
# Considering boxplot with no sd of ID, just from Island & Depth
ggplot(resume_depth_island_poc, aes(x=Depth, y=Size)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Corallite size") + ylab ("Size (mm)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

  
## d) with distinction between basal and apical regions
  
ggplot(resume_poc, aes(x=Depth, y=Distance, color=Part)) +
  geom_boxplot() + theme_bw() + facet_grid(Part ~ Island_Site) +
  ggtitle("Corallite size") + ylab("Size (mm)") + xlab ("Depth (m)")
  theme(plot.title = element_text(hjust=0.5, size=10, face="bold"),
      axis.text = element_text(size=10, colour="black"),
      axis.title = element_text(size=11, face="bold", colour="black"))
  
  
ggplot(resume_depth_island_poc, aes(x=Depth, y=Distance, color=Part)) +
  geom_boxplot() + theme_bw() + facet_grid(Part ~ Island_Site) +
  ggtitle("Corallite size") + ylab("Size (mm)") + xlab ("Depth (m)")
  theme(plot.title = element_text(hjust=0.5, size=10, face="bold"),
      axis.text = element_text(size=10, colour="black"),
      axis.title = element_text(size=11, face="bold", colour="black"))
  

# There are some values that very likely are outliers...


# Genereate database with ID and the fourth measures differentiating between Apical and Basal in colnames)
# Separate apical and basal to put both together after. I am sure there's a shorter way. I was not inspired today...
Apical <- morpho_poc %>% filter (Part == "apical") 
colnames(Apical) <- c ("Island_Site","Depth","ID_LAB","Size_Ap","Distance_Ap","Part")
Apical <- subset(Apical, select = -c(Part))
Basal <- morpho_poc %>% filter (Part == "basal") 
colnames(Basal) <- c ("Island_Site","Depth","ID_LAB","Size_Ba","Distance_Ba","Part")
Basal <- subset(Basal, select = -c(Part))
morpho_pocillopora <- cbind (Basal, Apical)
morpho_pocillopora <- morpho_pocillopora [,c(1:5,9,10)]
# 

Pocillopora_morphology <-ddply(morpho_pocillopora,~ID_LAB+Depth+Island_Site,function(x){c(Size_Ba=mean(na.omit(x$Size_Ba)), Distance_Ba=mean(na.omit(x$Distance_Ba)),Size_Ap=mean(na.omit(x$Size_Ap)), Distance_Ap=mean(na.omit(x$Distance_Ap)),
                                                                                      Size_Ba_sd=sd(na.omit(x$Size_Ba)), Distance_Ba_sd=sd(na.omit(x$Distance_Ba)),Size_Ap_sd=sd(na.omit(x$Size_Ap)), Distance_Ap_sd=sd(na.omit(x$Distance_Ap)))})

write_csv(Pocillopora_morphology, "Data_Codes/Morphology/Pocillopora_Morphology.csv")



# Make a super summary per Depth and Site
resume_Pocillopora_Morphology <-ddply(morpho_pocillopora,~Depth+Island_Site,function(x){c(Size_Ba=mean(na.omit(x$Size_Ba)), Distance_Ba=mean(na.omit(x$Distance_Ba)), Size_Ap=mean(na.omit(x$Size_Ap)),Distance_Ap=mean(na.omit(x$Distance_Ap)), Size_Ba_sd=sd(na.rm = FALSE,(x$Size_Ba)), Distance_Ba_sd=sd(na.rm = FALSE,(x$Distance_Ba)), Size_Ap_sd=sd(na.rm = FALSE,(x$Size_Ap)),Distance_Ap_sd=sd(na.rm = FALSE,(x$Distance_Ap)))})
summary(resume_Pocillopora_Morphology)





 