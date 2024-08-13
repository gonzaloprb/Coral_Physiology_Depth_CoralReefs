

# This script uses data on the Photoautotrophy Zoox, Chl and genreates the first preliminary plots

# Database generation. Stats. Plots vs light. 


# Extra info: There are outliers that need to be deleted or fixed such as for Bora_2 a Pachy speciosa at 90m that has 0 Chl2 (which is impossible) and others. 
# All these samples (with few exceptions) have the associated data of Morphological features, Isotopes and qualitative ITS2 profile types of Symbiodiniacea.


# Open database, prepare and separate between Pocillopora verrucosa and Pachyseris speciosa

rm (list = ls())

require (reshape2); require (ggplot2); require (plyr); require (tidyverse); require (dplyr); 

physio <- read.csv(file = "Data_Codes/Photoautotrophy/Physio_Export_CSV_02032020.csv", header = T, dec = ".", sep = ";")

colnames (physio) <- c("ID_LAB", "Species", "Island_Site", "Depth", "Zoox", "Chl_A", "Chl_C2")
summary (physio)
summary (as.factor(physio$Species)) 

# N = 234
# 4 sites with more than 50 samples
# For each species, site and depth, number of replicate
resume <- ddply(physio,~ Species + Island_Site +Depth,function(x){c(nbindiv=nrow(x))})
# Intermediate depths there are less than 5 replicates. 

quan.var <- c("Depth", "Zoox", "Chl_A", "Chl_C2")
physio [quan.var] <- sapply (physio [quan.var], as.numeric)

### Add the light do the databases ####
Light_loggers <- read.csv(file =  "Data_Codes/Light/Light_loggers.csv", header = T, dec = ".", sep = ",")

# Necessary to change the names:
Light_loggers$Island_Site <- gsub('Bora Bora','Bora_2', Light_loggers$Island_Site)
Light_loggers$Island_Site <- gsub('Moorea','Moorea_2', Light_loggers$Island_Site)
Light_loggers$Island_Site <- gsub('Rangiroa', 'Rangiroa_3', Light_loggers$Island_Site)
Light_loggers$Island_Site <- gsub('Tikehau','Tikehau_2', Light_loggers$Island_Site)

# physio_light <- merge (physio,Light_loggers, by = c("Island_Site","Depth"), all=T) # If you want to keep non-target depths for Pachyseris
physio_light <- merge (physio,Light_loggers, by = c("Island_Site","Depth"))

# The three measures we want are: Zoox, Chla/Zoox and ChlC2/Zoox (values of Xe-5 and Xe-6)
physio_light$RatioChlA <- physio_light$Chl_A/physio_light$Zoox
physio_light$RatioChlC2 <- physio_light$Chl_C2/physio_light$Zoox


# Necessary to think what do with intermediate depths
physio_Pocillopora <- filter (physio_light, Species == "Pocillopora verrucosa")
summary (as.factor(physio_Pocillopora$Depth))
# Make the replacement
physio_Pocillopora$Depth[physio_Pocillopora$Depth==62] <- 60

physio_Pachyseris <- filter (physio_light, Species == "Pachyseris speciosa")
summary (as.factor(physio_Pachyseris$Depth))
# Pachyseris has some spare depths in physio Pachyseris


### Graphs
################## Pocillopora verrucosa #####################

# We do not obtain clear differences for (Zoox, Chls) if we look at the variation in depth without considering sites. 

ggplot(physio_Pocillopora, aes(x = factor(Depth), y = Zoox)) + geom_boxplot() 
ggplot(physio_Pocillopora, aes(x = factor(Depth), y = Chl_A)) + geom_boxplot()  
ggplot(physio_Pocillopora, aes(x = factor(Depth), y = Chl_C2)) + geom_boxplot() 
ggplot(physio_Pocillopora, aes(x = factor(Depth), y = RatioChlA)) + geom_boxplot() 
ggplot(physio_Pocillopora, aes(x = factor(Depth), y = RatioChlC2)) + geom_boxplot() 


# I we separate the graphs by island sites we do have a variation in depth on Chl A and C2 

ggplot(physio_Pocillopora, aes(x = factor(Depth), y = Zoox)) + geom_boxplot() + facet_grid(~ Island_Site)
ggplot(physio_Pocillopora, aes(x = factor(Depth), y = Chl_A)) + geom_boxplot()  + facet_grid(~ Island_Site)
ggplot(physio_Pocillopora, aes(x = factor(Depth), y = Chl_C2)) + geom_boxplot()  + facet_grid(~ Island_Site)
ggplot(physio_Pocillopora, aes(x = factor(Depth), y = RatioChlA)) + geom_boxplot()  + facet_grid(~ Island_Site)
ggplot(physio_Pocillopora, aes(x = factor(Depth), y = RatioChlC2)) + geom_boxplot()  + facet_grid(~ Island_Site)



ggplot(physio_Pocillopora, aes(x=factor(Depth), y=Zoox)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, 
                                                                        color="red", size=2) + theme_bw()  + ylab ("Zooxanthellae (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


ggplot(physio_Pocillopora, aes(x=factor(Depth), y=Chl_A)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, 
                                                                        color="red", size=2) + theme_bw()  + ylab ("Chlorophill A (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pocillopora, aes(x=factor(Depth), y=Chl_C2)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, 
                                                                        color="red", size=2) + theme_bw()  + ylab ("Chlorophill C2 (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pocillopora, aes(x=factor(Depth), y=RatioChlA)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, 
                                                                        color="red", size=2) + theme_bw()  + ylab ("RatioChlA (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pocillopora, aes(x=factor(Depth), y=RatioChlC2)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, 
                                                                        color="red", size=2) + theme_bw()  + ylab ("RatioChlC2 (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 




# Ploting the points instead of boxplot to see the internal variability 

ggplot(physio_Pocillopora, aes(x=factor(Depth), y=Zoox)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw() + ylab ("Zooxanthellae (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pocillopora, aes(x=factor(Depth), y=Chl_A)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw() + ylab ("Chlorophill A (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pocillopora, aes(x=factor(Depth), y=Chl_C2)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw()  + ylab ("Chlorophill C2 (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pocillopora, aes(x=factor(Depth), y=RatioChlA)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw()  + ylab ("RatioChlA (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pocillopora, aes(x=factor(Depth), y=RatioChlC2)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw()  + ylab ("RatioChlC2 (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# Plot with light instead of depths
ggplot(physio_Pocillopora, aes(x=IndexLoss, y=Zoox)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw() + ylab ("Zooxanthellae (XX)") + xlab ("Index Loss Light") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pocillopora, aes(x=IndexLoss, y=RatioChlA)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw() + ylab ("RatioChlA (XX)") + xlab ("Index Loss Light") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pocillopora, aes(x=IndexLoss, y=RatioChlC2)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw() + ylab ("RatioChlC2 (XX)") + xlab ("Index Loss Light") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

# Important: From now on, we will work with "Depth" instead of "Light". Please see Methods in the manuscript.


# Ploting the three parameters together with secondary axis.
# Ploting them together is not necessary

ggplot(physio_Pocillopora) + 
  geom_boxplot(aes(x = factor (Depth), y = Zoox, color="Zoox"))+
  geom_boxplot(aes(x = factor (Depth), y = Chl_A*2e+05, color="Chl_A")) +
  scale_y_continuous("Zooxanthellae", sec.axis = sec_axis(~.*2e+05, name = "Chlorophyll A2",breaks=seq(0,20,40))) +
  geom_boxplot(aes(x = factor (Depth), y = Chl_C2*2e+05, color="Chl_C2"))+
  scale_y_continuous("Zooxanthellae", sec.axis = sec_axis(~.*2e+05, name = "Chlorophyll C2",breaks=seq(0,20,40)))+ 
  facet_wrap(Island_Site~., nrow = 1)


ggplot(physio_Pocillopora) + 
  geom_boxplot(aes(x = factor (Depth), y = Zoox, color="Zoox"))+
  geom_boxplot(aes(x = factor (Depth), y = RatioChlA, color="RatioChlA")) +
  geom_boxplot(aes(x = factor (Depth), y = RatioChlC2, color="RatioChlC2")) +
  scale_y_continuous("Photoautotrophy", sec.axis = sec_axis(~.*2e+05, name = "Chlorophyll A2",breaks=seq(0,20,40))) +
  facet_wrap(Island_Site~., nrow = 1)

# Resume and save dataframe
# Make a resume
ddply(physio_Pocillopora,~Depth+Island_Site,function(x){c(Zooxanthellae=mean(na.omit(x$Zoox)), Chlorophyll_a=mean(na.omit(x$Chl_A)), Chlorophyll_c2=mean(na.omit(x$Chl_C2)), zoox_sd=sd(na.rm = FALSE,(x$Zoox)), chl_a_sd=sd(na.rm = FALSE,(x$Chl_A)), chl_c2_sd=sd(na.rm = FALSE,(x$Chl_C2)))})


# Save the database of Zoox and Chl for modelling
write_csv(physio_Pocillopora, "Data_Codes/Photoautotrophy/Pocillopora_Photoautotrophy.csv")


################## Pachyseris verrucosa #####################

# We do not obtain clear differences for (Zoox, Chls) if we look at the variation in depth without considering sites. 

ggplot(physio_Pachyseris, aes(x = factor(Depth), y = Zoox)) + geom_boxplot() 
ggplot(physio_Pachyseris, aes(x = factor(Depth), y = Chl_A)) + geom_boxplot()  
ggplot(physio_Pachyseris, aes(x = factor(Depth), y = Chl_C2)) + geom_boxplot()  
ggplot(physio_Pachyseris, aes(x = factor(Depth), y = RatioChlA)) + geom_boxplot()  
ggplot(physio_Pachyseris, aes(x = factor(Depth), y = RatioChlC2)) + geom_boxplot()  



# I we separate the graphs by island sites we do have a variation in depth on Chl A and C2 

ggplot(physio_Pachyseris, aes(x = factor(Depth), y = Zoox)) + geom_boxplot() + facet_grid(~ Island_Site)
ggplot(physio_Pachyseris, aes(x = factor(Depth), y = Chl_A)) + geom_boxplot()  + facet_grid(~ Island_Site)
ggplot(physio_Pachyseris, aes(x = factor(Depth), y = Chl_C2)) + geom_boxplot()  + facet_grid(~ Island_Site)
ggplot(physio_Pachyseris, aes(x = factor(Depth), y = RatioChlA)) + geom_boxplot()  + facet_grid(~ Island_Site)
ggplot(physio_Pachyseris, aes(x = factor(Depth), y = RatioChlC2)) + geom_boxplot()  + facet_grid(~ Island_Site)

# Other preliminary version plots: 
ggplot(physio_Pachyseris, aes(x=factor(Depth), y=Zoox)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, 
                                                                        color="red", size=2) + theme_bw()  + ylab ("Zooxanthellae (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pachyseris, aes(x=factor(Depth), y=Chl_A)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, 
                                                                        color="red", size=2) + theme_bw()  + ylab ("Chlorophill A (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pachyseris, aes(x=factor(Depth), y=Chl_C2)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, 
                                                                        color="red", size=2) + theme_bw()  + ylab ("Chlorophill C2 (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pachyseris, aes(x=factor(Depth), y=RatioChlA)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, 
                                                                        color="red", size=2) + theme_bw()  + ylab ("RatioChlA (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pachyseris, aes(x=factor(Depth), y=RatioChlC2)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, 
                                                                        color="red", size=2) + theme_bw()  + ylab ("RatioChlC2 (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# Ploting the points instead of boxplot to see the internal variability 

ggplot(physio_Pachyseris, aes(x=factor(Depth), y=Zoox)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw() + ylab ("Zooxanthellae (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pachyseris, aes(x=factor(Depth), y=Chl_A)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw() + ylab ("Chlorophill A (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pachyseris, aes(x=factor(Depth), y=Chl_C2)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw()  + ylab ("Chlorophill C2 (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pachyseris, aes(x=factor(Depth), y=RatioChlA)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw()  + ylab ("RatioChlA (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(physio_Pachyseris, aes(x=factor(Depth), y=RatioChlC2)) +
  geom_point() + facet_grid(cols = vars(Island_Site)) + theme_bw()  + ylab ("RatioChlC2 (XX)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# Ploting the three parameters together with a plot with secondary axis to plot all together... I think not necessary
ggplot(physio_Pachyseris) + 
  geom_boxplot(aes(x = factor (Depth), y = Zoox, color="Zoox"))+
  geom_boxplot(aes(x = factor (Depth), y = Chl_A*2e+05, color="Chl_A")) +
  scale_y_continuous("Zooxanthellae", sec.axis = sec_axis(~.*2e+05, name = "Chlorophyll A2",breaks=seq(0,20,40))) +
  geom_boxplot(aes(x = factor (Depth), y = Chl_C2*2e+05, color="Chl_C2"))+
  scale_y_continuous("Zooxanthellae", sec.axis = sec_axis(~.*2e+05, name = "Chlorophyll C2",breaks=seq(0,20,40)))+ 
  facet_wrap(Island_Site~., nrow = 1)

# Make a resume
ddply(physio_Pachyseris,~Depth+Island_Site,function(x){c(Zooxanthellae=mean(na.omit(x$Zoox)), Chlorophyll_a=mean(na.omit(x$Chl_A)), Chlorophyll_c2=mean(na.omit(x$Chl_C2)), zoox_sd=sd(na.rm = FALSE,(x$Zoox)), chl_a_sd=sd(na.rm = FALSE,(x$Chl_A)), chl_c2_sd=sd(na.rm = FALSE,(x$Chl_C2)))})

# Save the database of Zoox and Chl for modelling
write_csv(physio_Pachyseris, "Data_Codes/Photoautotrophy/Pachyseris_Photoautotrophy.csv")

