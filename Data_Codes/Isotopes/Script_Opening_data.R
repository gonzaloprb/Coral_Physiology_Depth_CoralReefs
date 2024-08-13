# Script for opening the data and cleaning before analysis with Scripts: "Script_Pocillopora_Siber_Hotelling.R" and "Script_Pachyseris_Siber_Hotelling.R"

# Install necessary packages
require (ggplot2)
require (dplyr)
require (plyr)
require (tidyverse)

# Clean workspace
rm (list=ls())

### AAA Read Me
# Across the Manuscript codes and data:
# Host == Polyps
# Zoox == Symbionts
# In the Manuscript they are finally displayed as "Host" and "Symbionts"



# To be aware of: there are some outliers! - In orange in the excel file "Isotopes_runall_cleaned.xlsx"
# Open datatable 
Isotopes <- read.csv ("Data_Codes/Isotopes/Isotopes_runall_cleaned.csv", header = T, dec = ",", sep = ";") # With data cleaned. No negative and > 30  (δ15N)

# Replace -P and -Z by Polyps and Zoox. Separate them. 
Isotopes$Sample.ID <- gsub("-P", "-Polyps", Isotopes$Sample.ID)
Isotopes$Sample.ID <- gsub("-Z", "-Zoox", Isotopes$Sample.ID)


# Separate 
Isotopes <- separate(Isotopes, col = Sample.ID, into = c("ID_LAB","Extraction"), sep = "-")


# Combine dataframes to get information on site and depth for all Id_Labs
physio <- read.csv(file = "Data_Codes/Photoautotrophy/Physio_Export_CSV_02032020.csv", header = T, dec = ".", sep = ";")


colnames (physio) <- c("ID_LAB", "Species", "Island_Site", "Depth", "Zoox", "Chl_A", "Chl_C2")
names (physio)
physio <- subset(physio, select=-c(Zoox, Chl_A, Chl_C2))

# Combine dataframes 
Isotopes <- merge (Isotopes, physio, by = "ID_LAB")


# Add the surfaces and volumes into the database - necessary to correct the rawData values. 
surfaces_volumes <- read.csv(file = "Data_Codes/Isotopes/Isotopes_surfaces_volumes.csv", header = T, dec = ".", sep = ";")

# Change order of columns 
Isotopes <- Isotopes[c("ID_LAB", "Species", "Island_Site", "Depth" , "Extraction", "Weight..mg.","N2.Amp","X.N","δ15N.vs..At..Air","CO2.Amp","X.C","δ13C.vs..VPDB")]
# Change name of columns 
colnames (Isotopes) <- c("ID_LAB", "Species", "Island_Site", "Depth" , "Extraction", "Weight.(mg)","N2.Amp","N","δ15N","CO2.Amp","C","δ13C")



# Check database
# str (Isotopes)
quan.var <- c("Weight.(mg)","N2.Amp","N","δ15N","CO2.Amp","C","δ13C")
Isotopes [quan.var] <- sapply (Isotopes [quan.var], as.numeric)

# Need to measure the ratio of C/N
Isotopes$Ratio_C_N <- Isotopes$C / Isotopes$N

# Measure the deltas that Heloise call: Host - "I understand Polyps  and Zoox - " I understand symbionts"
df <- Isotopes


Deltas <- data.frame ()
for (i in unique (df$ID_LAB)) {
  ID_LAB <- unique (df$ID_LAB)
  Delta13C_Pol_Zoox  <-  df$δ13C[df$Extraction=="Polyps"] - df$δ13C[df$Extraction=="Zoox"]
  Delta15N_Pol_Zoox  <-  df$δ15N[df$Extraction=="Polyps"] - df$δ15N[df$Extraction=="Zoox"]
  Deltas <- cbind (ID_LAB,Delta13C_Pol_Zoox,Delta15N_Pol_Zoox)
}

Deltas <- as.data.frame (Deltas)
# As numeric variables
quan.var <- c("Delta13C_Pol_Zoox", "Delta15N_Pol_Zoox")
Deltas [quan.var] <- sapply (Deltas [quan.var], as.numeric)


# Merge dataframes
Isotopes <- merge (Isotopes, Deltas, by = "ID_LAB")


# Measure the ratios of Polyps/Zoox for each C and N
df <- Isotopes

Ratios <- data.frame ()
for (i in unique (df$ID_LAB)) {
  ID_LAB <- unique (df$ID_LAB)
  Ratio13C_Pol_Zoox  <-  df$δ13C[df$Extraction=="Polyps"] / df$δ13C[df$Extraction=="Zoox"]
  Ratio15N_Pol_Zoox  <-  df$δ15N[df$Extraction=="Polyps"] / df$δ15N[df$Extraction=="Zoox"]
  Ratios <- cbind (ID_LAB,Ratio13C_Pol_Zoox,Ratio15N_Pol_Zoox)
}

Ratios <- as.data.frame (Ratios)
# As numeric variables
quan.var <- c("Ratio13C_Pol_Zoox", "Ratio15N_Pol_Zoox")
Ratios [quan.var] <- sapply (Ratios [quan.var], as.numeric)


# Merge dataframes
Isotopes <- merge (Isotopes, Ratios, by = "ID_LAB")



# Separate between Pocillopora and Pachyseris
Isotopes_Pocillopora <- filter (Isotopes, Species == "Pocillopora verrucosa")
summary (as.factor(Isotopes_Pocillopora$Depth))

Isotopes_Pachyseris <- filter (Isotopes, Species == "Pachyseris speciosa")
summary (as.factor(Isotopes_Pachyseris$Depth))






# Quick plots 

##### Pocillopora #######

# Separatting Polyps and Zoox
ggplot(Isotopes_Pocillopora, aes(x=Depth, y=δ15N, shape = Extraction)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("δ15N P. verrucosa") + ylab ("δ15N") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(Isotopes_Pocillopora, aes(x=Depth, y=δ13C, shape = Extraction)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("δ13C P. verrucosa") + ylab ("δ13C") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))


# % of Nitrogen and Carbon
ggplot(Isotopes_Pocillopora, aes(x=Depth, y= N, shape = Extraction)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("%N P. verrucosa") + ylab ("%N") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

ggplot(Isotopes_Pocillopora, aes(x=Depth, y= C, shape = Extraction)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("%C P. verrucosa") + ylab ("%C") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

# Ratio C:N
ggplot(Isotopes_Pocillopora, aes(x=Depth, y= Ratio_C_N, shape = Extraction)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Ratio C & N") + ylab ("Ratio C & N") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

# Delta 13C and Delta 15N / It does not change according to extraction
ggplot(Isotopes_Pocillopora, aes(x=Depth, y= Delta13C_Pol_Zoox)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Delta13C_Pol_Zoox") + ylab ("Delta13C_Pol_Zoox") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

ggplot(Isotopes_Pocillopora, aes(x=Depth, y= Delta15N_Pol_Zoox)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Delta15N_Pol_Zoox") + ylab ("Delta15N_Pol_Zoox") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

# Ratios 13C and Ratio 15N (Polyps/Zoox)
ggplot(Isotopes_Pocillopora, aes(x=Depth, y= Ratio13C_Pol_Zoox)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Ratio13C_Pol_Zoox") + ylab ("Ratio13C_Pol_Zoox") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

ggplot(Isotopes_Pocillopora, aes(x=Depth, y= Ratio15N_Pol_Zoox)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Ratio15N_Pol_Zoox") + ylab ("Ratio15N_Pol_Zoox") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))


# Both Isotopes in the same graph (I can delete the ratio if I want) 
forms <- c("Polyps"=15,"Zoox"=17)
colV <- c("Carbon"="green","Nitrogen"="blue", "Ratio_C_N" = "black")
# In the same graph

ggplot(Isotopes_Pocillopora, aes(x = Depth,  shape = Extraction)) +
  geom_point(aes(y = C, color = "Carbon"), size = 1.5) +
  geom_point(aes(y = N, color = "Nitrogen"), size = 1.5) +
  geom_point(aes(y = Ratio_C_N, color = "Ratio_C_N"), size = 1.5) +
  scale_color_manual(values = colV) +
  labs(x = "Depth",y = "(%)",color = "Legend") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

# Both Isotopes deltas
forms <- c("Polyps"=15,"Zoox"=17)
colV <- c("Carbon"="green","Nitrogen"="blue")
# In the same graph

ggplot(Isotopes_Pocillopora, aes(x = Depth,  shape = Extraction)) +
  geom_point(aes(y = δ13C, color = "Carbon"), size = 1.5) +
  geom_point(aes(y = δ15N, color = "Nitrogen"), size = 1.5) +
  scale_color_manual(values = colV) +
  labs(x = "Depth",y = " δ (%)",color = "Legend") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# Plot with the points and deviations - Measure the mean and deviations for each (site depth and extraction)

# I have a serious doubt of how I should consider the site variability here...

# The best to avoid problems with the NA is to remove the rows
Isotopes_Pocillopora <- na.omit (Isotopes_Pocillopora)

Isotopes_Pocillopora_mean <- ddply(Isotopes_Pocillopora, ~ Island_Site + Depth + Extraction , function(x){c(Mean_δ13C = mean(x$δ13C), Sd_δ13C = sd(x$δ13C), Se_δ13C = sd(x$δ13C) / sqrt(length(x$δ13C)), 
                                                                                                          Mean_δ15N = mean(x$δ15N), Sd_δ15N = sd(x$δ15N), Se_δ15N = sd(x$δ15N) / sqrt(length(x$δ15N))) })

forms <- c("Polyps"=0,"Zoox"=1)
ggplot(Isotopes_Pocillopora_mean, aes(x = Mean_δ13C,  y = Mean_δ15N, shape = Extraction, colour = as.factor (Depth))) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin=Mean_δ15N-Se_δ15N, ymax=Mean_δ15N+Se_δ15N), width=.2,position=position_dodge(0.05)) +
  geom_errorbar(aes(xmin=Mean_δ13C-Se_δ13C, xmax=Mean_δ13C+Se_δ13C), width=.2,position=position_dodge(0.05)) +
  facet_grid(~Island_Site) +
  scale_shape_manual(values = forms) +
  labs(x = "δ13C",y = "δ15N") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 
  

# The interpretation is: More distance between the polyps and the zoox, more heterotrophy!

##### Pocillopora #######



##### Pachyseris #######
# Separatting Polyps and Zoox
ggplot(Isotopes_Pachyseris, aes(x=Depth, y=δ15N, shape = Extraction)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("δ15N P. speciosa") + ylab ("δ15N") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

ggplot(Isotopes_Pachyseris, aes(x=Depth, y=δ13C, shape = Extraction)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("δ13C P. speciosa") + ylab ("δ13C") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))


# % of Nitrogen and Carbon
ggplot(Isotopes_Pachyseris, aes(x=Depth, y= N, shape = Extraction)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("%N P. speciosa") + ylab ("%N") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

ggplot(Isotopes_Pachyseris, aes(x=Depth, y= C, shape = Extraction)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("%C P. speciosa") + ylab ("%C") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

# Ratio C:N
ggplot(Isotopes_Pachyseris, aes(x=Depth, y= Ratio_C_N, shape = Extraction)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Ratio C & N") + ylab ("Ratio C & N") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

# Delta 13C and Delta 15N / It does not change according to extraction
ggplot(Isotopes_Pachyseris, aes(x=Depth, y= Delta13C_Pol_Zoox)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Delta13C_Pol_Zoox") + ylab ("Delta13C_Pol_Zoox") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

ggplot(Isotopes_Pachyseris, aes(x=Depth, y= Delta15N_Pol_Zoox)) +
  geom_point(size = 3) + facet_grid(cols = vars(Island_Site)) + 
  stat_summary(fun=mean, geom="point", shape=18, color="red", size=2) + 
  theme_bw() + ggtitle("Delta15N_Pol_Zoox") + ylab ("Delta15N_Pol_Zoox") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))


# Both Isotopes in the same graph (I can delete the ratio if I want) 
forms <- c("Polyps"=15,"Zoox"=17)
colV <- c("Carbon"="green","Nitrogen"="blue", "Ratio_C_N" = "black")
# In the same graph

ggplot(Isotopes_Pachyseris, aes(x = Depth,  shape = Extraction)) +
  geom_point(aes(y = C, color = "Carbon"), size = 1.5) +
  geom_point(aes(y = N, color = "Nitrogen"), size = 1.5) +
  geom_point(aes(y = Ratio_C_N, color = "Ratio_C_N"), size = 1.5) +
  scale_color_manual(values = colV) +
  labs(x = "Depth",
       y = "(%)",
       color = "Legend") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

# Both Isotopes deltas
forms <- c("Polyps"=15,"Zoox"=17)
colV <- c("Carbon"="green","Nitrogen"="blue")
# In the same graph

ggplot(Isotopes_Pachyseris, aes(x = Depth,  shape = Extraction)) +
  geom_point(aes(y = δ13C, color = "Carbon"), size = 1.5) +
  geom_point(aes(y = δ15N, color = "Nitrogen"), size = 1.5) +
  scale_color_manual(values = colV) +
  labs(x = "Depth",
       y = " δ (%)",
       color = "Legend") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 


# Plot with the points and deviations - Measure the mean and deviations for each (site depth and extraction)

# I have a serious doubt of how I should consider the site variability here... (in the same plot it does not work!)

# The best to avoid problems with the NA is to remove the rows
Isotopes_Pachyseris <- na.omit (Isotopes_Pachyseris)

Isotopes_Pachyseris_mean <- ddply(Isotopes_Pachyseris, ~ Island_Site + Depth + Extraction , function(x){c(Mean_δ13C = mean(x$δ13C), Sd_δ13C = sd(x$δ13C), Se_δ13C = sd(x$δ13C) / sqrt(length(x$δ13C)), 
                                                                                                          Mean_δ15N = mean(x$δ15N), Sd_δ15N = sd(x$δ15N), Se_δ15N = sd(x$δ15N) / sqrt(length(x$δ15N))) })

forms <- c("Polyps"=0,"Zoox"=1)
ggplot(Isotopes_Pachyseris_mean, aes(x = Mean_δ13C,  y = Mean_δ15N, shape = Extraction, colour = as.factor (Depth))) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin=Mean_δ15N-Se_δ15N, ymax=Mean_δ15N+Se_δ15N), width=.2,position=position_dodge(0.05)) +
  geom_errorbar(aes(xmin=Mean_δ13C-Se_δ13C, xmax=Mean_δ13C+Se_δ13C), width=.2,position=position_dodge(0.05)) +
  facet_grid(~Island_Site) +
  scale_shape_manual(values = forms) +
  labs(x = "δ13C",y = "δ15N") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 

# The interpretation is: More distance between the polyps and the zoox, more heterotrophy!

##### Pachyseris #######


###############################################################

# There are some samples that bring problems: 163, 1615, 1861, 1149, 1840 (163, 1769) 
# These were also highlighted in yellow during the extraction process because they had problems except: 163, 1769

# Pocillopora and PAchy
Isotopes_Pocillopora <- Isotopes_Pocillopora[- grep(c("163|1615|1861|1149"), Isotopes_Pocillopora$ID_LAB),]

# Already deleted from the Excel or NA
# Isotopes_Pachyseris <- Isotopes_Pachyseris[- grep(c("1840"), Isotopes_Pachyseris$ID_LAB),]

# Some extra values that do not make sense
Isotopes_Pachyseris <- Isotopes_Pachyseris[- grep(c("2566|857"), Isotopes_Pachyseris$ID_LAB),]


# Also delete random depths, where very little replicates

# Already deleted from the Excel or NA
# Isotopes_Pocillopora <- Isotopes_Pocillopora[- grep(c("62"), Isotopes_Pocillopora$Depth),]

Isotopes_Pachyseris <- Isotopes_Pachyseris[- grep(c("18|47"), Isotopes_Pachyseris$Depth),]

# Make a resume
# Save the database of istopes and resume


# Pocillopora
resume_Pocillopora_Isotopes <-ddply(Isotopes_Pocillopora,~Depth+Island_Site+Extraction,
                                   function(x){c(Nitrogen=mean(na.omit(x$N)), Delta_Nitrogen15=mean(na.omit(x$δ15N)), 
                                                 Carbon=mean(na.omit(x$C)), Delta_Carbon13=mean(na.omit(x$δ13C)),
                                                 Ratio_C_N=mean(na.omit(x$Ratio_C_N)),
                                                 Delta_Nitrogen15_SD=sd(na.rm = T,(x$δ15N)),
                                                 Delta_Carbon13_SD=sd(na.rm = T,(x$δ13C)),
                                                 Ratio_C_N_SD=sd(na.rm = T,(x$Ratio_C_N)),
                                                 Delta13C_Pol_Zoox=mean(na.rm = T,(x$Delta13C_Pol_Zoox)),
                                                 Delta15N_Pol_Zoox=mean(na.rm = T,(x$Delta15N_Pol_Zoox)),
                                                 Delta13C_Pol_Zoox_SD=sd(na.rm = T,(x$Delta13C_Pol_Zoox)),
                                                 Delta15N_Pol_Zoox_SD=sd(na.rm = T,(x$Delta15N_Pol_Zoox)),
                                                 Ratio13C_Pol_Zoox=mean(na.rm = T,(x$Ratio13C_Pol_Zoox)),
                                                 Ratio15N_Pol_Zoox=mean(na.rm = T,(x$Ratio15N_Pol_Zoox)),
                                                 Ratio13C_Pol_Zoox_SD=sd(na.rm = T,(x$Ratio13C_Pol_Zoox)),
                                                 Ratio15N_Pol_Zoox_SD=sd(na.rm = T,(x$Ratio15N_Pol_Zoox)))})
summary(resume_Pocillopora_Isotopes)
# Multiple means
resume_Pocillopora_Isotopes %>%
  group_by(Extraction,Depth) %>% 
  summarise_at(vars("Nitrogen","Delta_Nitrogen15", "Carbon", "Delta_Carbon13","Ratio_C_N","Delta13C_Pol_Zoox","Delta15N_Pol_Zoox", "Ratio13C_Pol_Zoox", "Ratio15N_Pol_Zoox"), mean)


# Pachyseris
resume_Pachyseris_Isotopes <-ddply(Isotopes_Pachyseris,~Depth+Island_Site+Extraction,
                                   function(x){c(Nitrogen=mean(na.omit(x$N)), Delta_Nitrogen15=mean(na.omit(x$δ15N)), 
                                                Carbon=mean(na.omit(x$C)), Delta_Carbon13=mean(na.omit(x$δ13C)),
                                                Ratio_C_N=mean(na.omit(x$Ratio_C_N)),
                                                Delta_Nitrogen15_SD=sd(na.rm = T,(x$δ15N)),
                                                Delta_Carbon13_SD=sd(na.rm = T,(x$δ13C)),
                                                Ratio_C_N_SD=sd(na.rm = T,(x$Ratio_C_N)),
                                                Delta13C_Pol_Zoox=mean(na.rm = T,(x$Delta13C_Pol_Zoox)),
                                                Delta15N_Pol_Zoox=mean(na.rm = T,(x$Delta13C_Pol_Zoox)),
                                                Delta13C_Pol_Zoox_SD=sd(na.rm = T,(x$Delta15N_Pol_Zoox)),
                                                Delta15N_Pol_Zoox_SD=sd(na.rm = T,(x$Delta15N_Pol_Zoox)),
                                                Ratio13C_Pol_Zoox=mean(na.rm = T,(x$Ratio13C_Pol_Zoox)),
                                                Ratio15N_Pol_Zoox=mean(na.rm = T,(x$Ratio15N_Pol_Zoox)),
                                                Ratio13C_Pol_Zoox_SD=sd(na.rm = T,(x$Ratio13C_Pol_Zoox)),
                                                Ratio15N_Pol_Zoox_SD=sd(na.rm = T,(x$Ratio15N_Pol_Zoox)))})

summary(resume_Pachyseris_Isotopes)
resume_Pachyseris_Isotopes %>%
  group_by(Extraction,Depth) %>% 
  summarise_at(vars("Nitrogen","Delta_Nitrogen15", "Carbon", "Delta_Carbon13","Ratio_C_N","Delta13C_Pol_Zoox","Delta15N_Pol_Zoox","Ratio13C_Pol_Zoox", "Ratio15N_Pol_Zoox"), mean)



# Save resumed databases
# Pocillopora
write_csv(resume_Pocillopora_Isotopes, "Data_Codes/Isotopes/Pocillopora_Isotopes_resume.csv")
# Pachyseris
write_csv(resume_Pachyseris_Isotopes, "Data_Codes/Isotopes/Pachyseris_Isotopes_resume.csv")

# Keep only the variables we need for Siber_Hotelling tests and save datasets in the right format. 
Isotopes_Pocillopora <- select(Isotopes_Pocillopora, ID_LAB, Island_Site, Depth, Extraction, δ13C, δ15N,  Ratio_C_N, Delta13C_Pol_Zoox,Delta15N_Pol_Zoox,Ratio13C_Pol_Zoox, Ratio15N_Pol_Zoox)
Isotopes_Pocillopora <- Isotopes_Pocillopora[!duplicated(Isotopes_Pocillopora), ]
write_csv(Isotopes_Pocillopora, "Data_Codes/Isotopes/Pocillopora_Isotopes.csv")


Isotopes_Pachyseris <- select(Isotopes_Pachyseris, ID_LAB, Island_Site, Depth, Extraction,  δ13C, δ15N, Ratio_C_N, Delta13C_Pol_Zoox,Delta15N_Pol_Zoox,Ratio13C_Pol_Zoox, Ratio15N_Pol_Zoox)
Isotopes_Pachyseris <- Isotopes_Pachyseris[!duplicated(Isotopes_Pachyseris), ]
write_csv(Isotopes_Pachyseris, "Data_Codes/Isotopes/Pachyseris_Isotopes.csv")

# Prepare datasets for All_Features Script
# Pocillopora
# Make the Polyps and Zoox as new columns and not keep them as category column 
# Iso 13
C13 <- Isotopes_Pocillopora[,c( "ID_LAB","Island_Site","Depth","Extraction", "δ13C")]
C13 <- dcast (C13, ID_LAB + Island_Site + Depth ~ Extraction, fun=sum, value.var = c("δ13C"))
colnames (C13) <- c ("ID_LAB","Island_Site","Depth","iso13C_Polyps", "iso13C_Zoox")

N15 <- Isotopes_Pocillopora[,c( "ID_LAB","Island_Site","Depth","Extraction", "δ15N")]
N15 <- dcast (N15, ID_LAB + Island_Site + Depth ~ Extraction, fun=sum, value.var = c("δ15N"))
colnames (N15) <- c ("ID_LAB","Island_Site","Depth","iso15N_Polyps", "iso15N_Zoox")

Ratio_iso13C_iso15N <- Isotopes_Pocillopora[,c( "ID_LAB","Island_Site","Depth","Extraction", "Ratio_C_N")]
Ratio_iso13C_iso15N <- dcast (Ratio_iso13C_iso15N, ID_LAB + Island_Site + Depth ~ Extraction, fun=sum, value.var = c("Ratio_C_N"))
colnames (Ratio_iso13C_iso15N) <- c ("ID_LAB","Island_Site","Depth","Ratio_CN_Polyps", "Ratio_CN_Zoox")


Pocillopora_Polyps_Zoox <- merge (C13, N15, by = c("ID_LAB","Island_Site","Depth") , all = T)
Pocillopora_Polyps_Zoox <- merge (Pocillopora_Polyps_Zoox, Ratio_iso13C_iso15N, by = c("ID_LAB","Island_Site","Depth") , all = T)

# Delete columns and duplicate of Isotopes_Pocillopora - and delete duplicated
Isotopes_Pocillopora_all <- Isotopes_Pocillopora[,c( "ID_LAB","Island_Site","Depth","Delta13C_Pol_Zoox", "Delta15N_Pol_Zoox")]
Isotopes_Pocillopora_all <- Isotopes_Pocillopora_all[duplicated(Isotopes_Pocillopora_all), ]

# Merge
Pocillopora_Polyps_Zoox <- merge (Pocillopora_Polyps_Zoox, Isotopes_Pocillopora_all, by = c("ID_LAB","Island_Site","Depth") , all = T)
write_csv(Pocillopora_Polyps_Zoox, "Data_Codes/Isotopes/Pocillopora_Isotopes_Polyps_Zoox.csv")


# Make the same for Pachyseris
# Iso 13
C13 <- Isotopes_Pachyseris[,c( "ID_LAB","Island_Site","Depth","Extraction", "δ13C")]
C13 <- dcast (C13, ID_LAB + Island_Site + Depth ~ Extraction, fun=sum, value.var = c("δ13C"))
colnames (C13) <- c ("ID_LAB","Island_Site","Depth","iso13C_Polyps", "iso13C_Zoox")

N15 <- Isotopes_Pachyseris[,c( "ID_LAB","Island_Site","Depth","Extraction", "δ15N")]
N15 <- dcast (N15, ID_LAB + Island_Site + Depth ~ Extraction, fun=sum, value.var = c("δ15N"))
colnames (N15) <- c ("ID_LAB","Island_Site","Depth","iso15N_Polyps", "iso15N_Zoox")

Ratio_iso13C_iso15N <- Isotopes_Pachyseris[,c( "ID_LAB","Island_Site","Depth","Extraction", "Ratio_C_N")]
Ratio_iso13C_iso15N <- dcast (Ratio_iso13C_iso15N, ID_LAB + Island_Site + Depth ~ Extraction, fun=sum, value.var = c("Ratio_C_N"))
colnames (Ratio_iso13C_iso15N) <- c ("ID_LAB","Island_Site","Depth","Ratio_CN_Polyps", "Ratio_CN_Zoox")


Pachyseris_Polyps_Zoox <- merge (C13, N15, by = c("ID_LAB","Island_Site","Depth") , all = T)
Pachyseris_Polyps_Zoox <- merge (Pachyseris_Polyps_Zoox, Ratio_iso13C_iso15N, by = c("ID_LAB","Island_Site","Depth") , all = T)

# Delete columns and duplicate of Isotopes_Pachyseris - and delete duplicated
Isotopes_Pachyseris_all <- Isotopes_Pachyseris[,c( "ID_LAB","Island_Site","Depth","Delta13C_Pol_Zoox", "Delta15N_Pol_Zoox")]
Isotopes_Pachyseris_all <- Isotopes_Pachyseris_all[duplicated(Isotopes_Pachyseris_all), ]

# Merge
Pachyseris_Polyps_Zoox <- merge (Pachyseris_Polyps_Zoox, Isotopes_Pachyseris_all, by = c("ID_LAB","Island_Site","Depth") , all = T)

write_csv(Pachyseris_Polyps_Zoox, "Data_Codes/Isotopes/Pachyseris_Isotopes_Polyps_Zoox.csv")

