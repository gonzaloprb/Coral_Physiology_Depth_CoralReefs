# Pocillopora analysis between Apical and basal parts of the colony.
# Sup. Fig. 2 and 3 and Sup. Table 1.

# Read Morphology database
Physio_Pocillopora_morpho <- read.csv ("Data_Codes/Morphology/Pocillopora_Morphology.csv", header = T, dec = ",", sep = ",")

rownames (Physio_Pocillopora_morpho) <- Physio_Pocillopora_morpho[,1]

# Select Morphology variables 
Physio_Pocillopora_morpho <- Physio_Pocillopora_morpho [,c(3,2,4,5,6,7)]


# Set quantitative variables as numeric
quan.var <- c("Depth","Size_Ba","Distance_Ba","Size_Ap","Distance_Ap")
Physio_Pocillopora_morpho [quan.var] <- sapply (Physio_Pocillopora_morpho [quan.var], as.numeric)

# Change the names of the Sites
Physio_Pocillopora_morpho$Island_Site <- gsub('Bora_2', 'Bora Bora', Physio_Pocillopora_morpho$Island_Site)
Physio_Pocillopora_morpho$Island_Site <- gsub('Moorea_2', 'Moorea', Physio_Pocillopora_morpho$Island_Site)
Physio_Pocillopora_morpho$Island_Site <- gsub('Rangiroa_3', 'Rangiroa', Physio_Pocillopora_morpho$Island_Site)
Physio_Pocillopora_morpho$Island_Site <- gsub('Tikehau_2', 'Tikehau', Physio_Pocillopora_morpho$Island_Site)

Physio_Pocillopora_morpho = Physio_Pocillopora_morpho %>% mutate(Island_Site = factor(Island_Site, levels = c("Moorea","Bora Bora","Tikehau","Rangiroa")))


# Statistical test to compare between the basal and apical (compare two numerical variables in function of depth)
# Distance corallites
distance_morpho <- Physio_Pocillopora_morpho [,c(1,2,4,6)]
# Extract island
colnames (distance_morpho) [1] <- "Island"

distance_morpho <- melt (distance_morpho, id.vars = c ("Depth", "Island"),na.rm = F, value.name = c("Distance"))
colnames (distance_morpho) <- c ("Depth", "Island", "Part", "Distance")


distance_morpho$Part <- gsub('Distance_Ba', 'Basal', distance_morpho$Part)
distance_morpho$Part <- gsub('Distance_Ap', 'Apical', distance_morpho$Part)


ggplot(distance_morpho, aes(x=factor(Part), y=Distance)) +
  geom_boxplot() + facet_grid(cols = vars(Island), vars (Depth)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Distance (mm)") + xlab ("Part of the colony") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                               axis.text = element_text(size=10, colour="black"),
                                                               axis.title = element_text(size=11, face="bold", colour="black")) 

Anova(lm(Distance ~ Depth*Part, data=distance_morpho, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Distance ~ Part*Depth, data=distance_morpho, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3) # Same result
# Statistical difference in the distance between corallites between the two parts

# Size corallites
size_morpho <- Physio_Pocillopora_morpho [,c(1,2,3,5)]
colnames (size_morpho) [1] <- "Island"

size_morpho <- melt (size_morpho, id.vars = c ("Depth", "Island"),na.rm = F, value.name = c("Size"))
colnames (size_morpho) <- c ("Depth", "Island", "Part", "Size")


size_morpho$Part <- gsub('Size_Ba', 'Basal', size_morpho$Part)
size_morpho$Part <- gsub('Size_Ap', 'Apical', size_morpho$Part)

ggplot(size_morpho, aes(x=factor(Part), y=Size)) +
  geom_boxplot() + facet_grid(cols = vars(Island), vars (Depth)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Size (mm)") + xlab ("Part of the colony") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                           axis.text = element_text(size=10, colour="black"),
                                                           axis.title = element_text(size=11, face="bold", colour="black")) 

Anova(lm(Size ~ Depth*Part, data=size_morpho, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Size ~ Part*Depth, data=size_morpho, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
# Statistical difference in the size of corallites between the two parts

# This Anova tests are used to generate Supplementary Table 1. 
# Visual photos are also available in Supplementary Figure 2 and 3
