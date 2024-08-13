# Script to open all together. PhotoAutotropy, Morphology, Isotopes and Adatation

rm(list =ls())

require (reshape2); require (ggplot2); require (plyr); require (tidyverse); require (dplyr); library(randomcoloR); library (RColorBrewer)


# Open Photoautotrophy
setwd("~/Documents/AAASea_Science/AAA_PhD_Thesis/Physio_Lab/Photoautotrophy")
Pocillopora_Photoautotrophy <- read.csv ("Pocillopora_Photoautotrophy.csv", header = T, dec = ",", sep = ",")


# Open Symbiodiniacea clean

setwd("~/Documents/AAASea_Science/AAA_PhD_Thesis/Physio_Lab/Symbiodiniacea")
Poc_Symbio <- read.csv ("Pocillopora_Symbio_Clean.csv", header = T, dec = ",", sep = ",")



# Measure relative values of Symbio Profiles - must be an easier way but I am tired!

Poc_Symbio$Total <- rowSums(Poc_Symbio[,5:35])

# Create new dataframe with relative values
Relative_Poc_Symbio <- Poc_Symbio[,c("ID_LAB","Depth","Island_Site","Colbar")]
Relative_Poc_Symbio2 <- Poc_Symbio [,5:35] / Poc_Symbio$Total
# Remove columns (Symbio Profiles) where all are 0
Relative_Poc_Symbio2 <- Relative_Poc_Symbio2[, colSums(Relative_Poc_Symbio2) != 0]

# Combine the initial values with the relatives
Relative_Poc_Symbio <- cbind (Relative_Poc_Symbio,Relative_Poc_Symbio2)

# Make columns of all 1s
Relative_Poc_Symbio$Total <- rowSums(Relative_Poc_Symbio[,5:16])



Percentatges_Profiles <- melt(Relative_Poc_Symbio, id = c("ID_LAB","Island_Site","Depth","Colbar", "Total"))

colnames (Percentatges_Profiles) <- c( "ID_LAB","Island_Site","Depth","Colbar","Total","Profile_ITS2","Proportion")


# As you won't be able to keep the unique sample ID, you need to work with averages for site and depth!
Percentatges_Profiles_mean <- ddply(Percentatges_Profiles, ~ Island_Site + Depth + Profile_ITS2, function(x){c(Proportion = mean(x$Proportion), Sd = sd(x$Proportion), Se=sd(x$Proportion) / sqrt(length(x$Proportion)))}) 

colours <- distinctColorPalette(12)

ggplot(data = Percentatges_Profiles_mean, aes(x = factor(Island_Site), y = Proportion))+
  geom_bar(aes(fill=factor(Profile_ITS2)), stat="identity", position = position_stack(reverse = TRUE)) + facet_wrap(~  Depth, ncol = 4) +
  scale_y_continuous(position = "right",breaks = c(0,0.25,0.5,0.75,1))+ scale_x_discrete ()+ ylab ("Proportion (%)") + xlab ("") + 
  theme_classic() + theme(legend.position = "bottom") 


ggplot(Percentatges_Profiles_mean, aes(x=Depth, y=Proportion)) +
  geom_bar(aes(fill=factor(Profile_ITS2)), stat="identity") + 
  scale_fill_manual(values= colours) +
  facet_wrap(~Island_Site,ncol = 4) + 
  scale_x_continuous(name ="Depth (m)", breaks = c(6,20,40,60)) +
  scale_y_continuous(name ="Proportion (%)", limits=c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=8, colour="black"),
                          axis.title = element_text(size=12, face="bold", colour="black"), 
                          strip.text = element_text(size=10),
                          legend.title = element_blank(), legend.position = "bottom") 


ggplot(Percentatges_Profiles_mean, aes(x=Island_Site, y=Proportion)) +
  geom_bar(aes(fill=factor(Profile_ITS2)), stat="identity", width = 0.3) +  facet_grid (rows = vars(Depth), switch = "y") +
  scale_fill_manual(values= colours) +
  scale_y_continuous(position = "right",breaks = c(0,0.25,0.5,0.75,1))+ scale_x_discrete ()+ ylab ("Proportion (%)") + xlab ("") + 
  ggtitle("Relative Symbiotic Commuinities per depth")+
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=16, face="bold"),
                          axis.text.y = element_text(size=10, colour="black"),
                          axis.text.x = element_text(size=10, colour="black", angle = 0),
                          strip.text.x = element_text(angle = 0),
                          strip.text.y = element_text(size=16, colour="black"),
                          axis.title = element_text(size=14, face="bold", colour="black"), legend.position = "bottom") 



# Without the mean and keeping all sample ID to make similar plot to Heloise
# Summary of the number of samples I have for each Island_Site and Depth
ddply(Percentatges_Profiles, ~ Island_Site + Depth, function(x){c(Nsamples = length(unique(x$ID_LAB)))})

# Delete all columns where I have 0
# Not sure if it does really make a difference!
Percentatges_Profiles <- Percentatges_Profiles[Percentatges_Profiles$Proportion != 0, ]

ggplot(Percentatges_Profiles, aes(x=Island_Site, y=Proportion)) +
  geom_bar(aes(fill=factor(Profile_ITS2)), stat="identity", position = "dodge") +  facet_grid (rows = vars(Depth), switch = "y") +
  scale_fill_manual(values= colours) +
  scale_y_continuous(position = "right",breaks = c(0,0.25,0.5,0.75,1))+ scale_x_discrete ()+ ylab ("Proportion (%)") + xlab ("") + 
  ggtitle("Relative Symbiotic Commuinities per depth")+
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=16, face="bold"),
                          axis.text.y = element_text(size=10, colour="black"),
                          axis.text.x = element_text(size=10, colour="black", angle = 0),
                          strip.text.x = element_text(angle = 0),
                          strip.text.y = element_text(size=16, colour="black"),
                          axis.title = element_text(size=14, face="bold", colour="black"), legend.position = "bottom") 

ggplot(Percentatges_Profiles, aes(x=as.factor (ID_LAB), y=Proportion)) +
  geom_bar(aes(fill=factor(Profile_ITS2)),colour="black", stat="identity", width = 1) +  
  facet_grid (Depth ~ Island_Site, switch = "y", scales="free_x") +
  scale_fill_manual(values= colours) +
  scale_y_continuous(position = "right",breaks = c(0,0.25,0.5,0.75,1))+ scale_x_discrete ()+ ylab ("Proportion (%)") + xlab ("") + 
  ggtitle("Relative Symbiotic Commuinities per depth")+
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=16, face="bold"),
                          axis.text.y = element_text(size=10, colour="black"),
                          axis.text.x = element_text(size=6, colour="black", angle = 90),
                          strip.text.x = element_text(angle = 0),
                          strip.text.y = element_text(size=16, colour="black"),
                          axis.title = element_text(size=14, face="bold", colour="black"), legend.position = "bottom") 





ggplot(Percentatges_Profiles, aes(x=as.factor (ID_LAB), y=Proportion)) +
  geom_bar(aes(fill=factor(Profile_ITS2)), stat="identity", width = 1, position = "dodge") +  
  facet_grid (rows = vars(Depth), switch = "y") +
  scale_fill_manual(values= colours) +
  scale_y_continuous(position = "right",breaks = c(0,0.25,0.5,0.75,1))+ scale_x_discrete ()+ ylab ("Proportion (%)") + xlab ("") + 
  ggtitle("Relative Symbiotic Commuinities per depth")+
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=16, face="bold"),
                          axis.text.y = element_text(size=10, colour="black"),
                          axis.text.x = element_text(size=10, colour="black", angle = 90),
                          strip.text.x = element_text(angle = 0),
                          strip.text.y = element_text(size=16, colour="black"),
                          axis.title = element_text(size=14, face="bold", colour="black"), legend.position = "bottom") 

 