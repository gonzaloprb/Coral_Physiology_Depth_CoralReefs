# Script to open all together. PhotoAutotropy, Morphology, Isotopes and Adatation

rm(list =ls())

require (reshape2); require (ggplot2); require (plyr); require (tidyverse); require (dplyr); library(randomcoloR); library (RColorBrewer)

# Venn Diagramm
library(VennDiagram)

# Open database
setwd("~/Documents/AAASea_Science/PhD_Thesis/Physio_Lab/Symbiodiniacea")
SymPocillopora <- read.csv ("SymPocillopora.csv", header = T, dec = ",", sep = ";")



# Generate the different sets
# 6 m set
set1_six <- filter(SymPocillopora, Depth == 6)


# How many times each profile appears
set1 <- set1_six[,5:16] %>% summarize_if(is.numeric, sum, na.rm=TRUE)

set1[,1]

set1_six <- rep(colnames (set1)[1], set1[,1]) # first ITS2 profile, there are 11
set1_six2 <- rep(colnames (set1)[2], set1[,2]) #second profile
set1_six3 <- rep(colnames (set1)[3], set1[,3])
set1_six4 <- rep(colnames (set1)[4], set1[,4])
set1_six5 <- rep(colnames (set1)[5], set1[,5])
set1_six6 <- rep(colnames (set1)[6], set1[,6]) 
set1_six7 <- rep(colnames (set1)[7], set1[,7])
set1_six8 <- rep(colnames (set1)[8], set1[,8])
set1_six9 <- rep(colnames (set1)[9], set1[,9])
set1_six10 <- rep(colnames (set1)[10], set1[,10])
set1_six11 <- rep(colnames (set1)[11], set1[,11])


set_final_six <- c(set1_six,set1_six2,set1_six3,set1_six4,set1_six5,set1_six6,set1_six7,set1_six8,
                   set1_six9,set1_six10,set1_six11)



# 20 m set
set2_twenty <- filter(SymPocillopora, Depth == 20)


# How many times each profile appears
set2 <- set2_twenty[,5:16] %>% summarize_if(is.numeric, sum, na.rm=TRUE)


set2_twenty <- rep(colnames (set2)[1], set2[,1]) # first ITS2 profile, there are 11
set2_twenty2 <- rep(colnames (set2)[2], set2[,2]) #second profile
set2_twenty3 <- rep(colnames (set2)[3], set2[,3])
set2_twenty4 <- rep(colnames (set2)[4], set2[,4])
set2_twenty5 <- rep(colnames (set2)[5], set2[,5])
set2_twenty6 <- rep(colnames (set2)[6], set2[,6])
set2_twenty7 <- rep(colnames (set2)[7], set2[,7])
set2_twenty8 <- rep(colnames (set2)[8], set2[,8])
set2_twenty9 <- rep(colnames (set2)[9], set2[,9])
set2_twenty10 <- rep(colnames (set2)[10], set2[,10])
set2_twenty11 <- rep(colnames (set2)[11], set2[,11])



set_final_twenty <- c(set2_twenty,set2_twenty2,set2_twenty3,set2_twenty4,set2_twenty5,set2_twenty6,
                      set2_twenty7,set2_twenty8,set2_twenty9,set2_twenty10,set2_twenty11)



# 40 m set
set3_forty <- filter(SymPocillopora, Depth == 40)


# How many times each profile appears
set3 <- set3_forty[,5:16] %>% summarize_if(is.numeric, sum, na.rm=TRUE)


set3_forty <- rep(colnames (set3)[1], set3[,1]) # first ITS2 profile, there are 11
set3_forty2 <- rep(colnames (set3)[2], set3[,2]) #second profile
set3_forty3 <- rep(colnames (set3)[3], set3[,3])
set3_forty4 <- rep(colnames (set3)[4], set3[,4])
set3_forty5 <- rep(colnames (set3)[5], set3[,5])
set3_forty6 <- rep(colnames (set3)[6], set3[,6])
set3_forty7 <- rep(colnames (set3)[7], set3[,7])
set3_forty8 <- rep(colnames (set3)[8], set3[,8])
set3_forty9 <- rep(colnames (set3)[9], set3[,9])
set3_forty10 <- rep(colnames (set3)[10], set3[,10])
set3_forty11 <- rep(colnames (set3)[11], set3[,11])



set_final_forty <- c(set3_forty,set3_forty2,set3_forty3,set3_forty4,set3_forty5,set3_forty6,
                      set3_forty7,set3_forty8,set3_forty9,set3_forty10,set3_forty11)



# 60 m set
set4_sixty <- filter(SymPocillopora, Depth == 60)


# How many times each profile appears
set4 <- set4_sixty[,5:16] %>% summarize_if(is.numeric, sum, na.rm=TRUE)


set4_sixty <- rep(colnames (set4)[1], set4[,1]) # first ITS2 profile, there are 11
set4_sixty2 <- rep(colnames (set4)[2], set4[,2]) #second profile
set4_sixty3 <- rep(colnames (set4)[3], set4[,3])
set4_sixty4 <- rep(colnames (set4)[4], set4[,4])
set4_sixty5 <- rep(colnames (set4)[5], set4[,5])
set4_sixty6 <- rep(colnames (set4)[6], set4[,6])
set4_sixty7 <- rep(colnames (set4)[7], set4[,7])
set4_sixty8 <- rep(colnames (set4)[8], set4[,8])
set4_sixty9 <- rep(colnames (set4)[9], set4[,9])
set4_sixty10 <- rep(colnames (set4)[10], set4[,10])
set4_sixty11 <- rep(colnames (set4)[11], set4[,11])



set_final_sixty <- c(set4_sixty,set4_sixty2,set4_sixty3,set4_sixty4,set4_sixty5,set4_sixty6,
                     set4_sixty7,set4_sixty8,set4_sixty9,set4_sixty10,set4_sixty11)



# Simplest graph
venn.diagram(
  x = list(set_final_six, set_final_twenty, set_final_forty,set_final_sixty),
  category.names = c("6 m" , "20 m" , "40 m", "60 m"),
  filename = '#14_venn_diagramm.png',
  output=TRUE, 
)



# Set the colors
myCol <- cols <- c("6" = "#BDD7E7", "20" = "#6BAED6", "40" = "#3182BD", "60" = "#08519C")


# Chart
venn.diagram(
  x = list(set_final_six, set_final_twenty, set_final_forty,set_final_sixty),
  category.names = c("6m", "20m", "40m","60m"),
  filename = 'Depth_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, -135, 135),
  cat.dist = c(0.055, 0.055, 0.085, 0.085),
  cat.fontfamily = "sans"
  )




#Make the plot
venn.diagram(
  x = list(set_final_six, set_final_twenty, set_final_forty, set_final_sixty),
  category.names = c("6m" , "20m", "40m","60m"),
  filename = 'Pocillopora_Depth_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#BDD7E7", '#6BAED6', '#3182BD', '#08519C'),
  fill = c(alpha("#BDD7E7",0.3), alpha('#6BAED6',0.3), alpha('#3182BD',0.3),alpha('#08519C',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, -135, 135),
  cat.dist = c(0.055, 0.055, 0.085, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#BDD7E7", '#6BAED6', '#3182BD', '#08519C')
)






library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(set_final_six, set_final_twenty),
  category.names = c("6 m" , "20 m"),
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)





# Generate the sets
# resume_pachy <-ddply(morpho_pachy,~ID_LAB+Depth+Island_Site,function(x){c(Septa.Distance=mean(na.omit(x$Septa.Distance)), Width.Valley=mean(na.omit(x$Width.Valley)), Height.Septa=mean(na.omit(x$Height.Septa)), septa_sd=sd(na.rm = FALSE,(x$Septa.Distance)), height_sd=sd(na.rm = FALSE,(x$Height.Septa)), valley_sd=sd(na.rm = FALSE,(x$Width.Valley)))})



 
# ID_Labs where I have Symbio and not Autotrophy
setdiff(unique (Poc_Symbio$ID_LAB),unique (Pocillopora_Photoautotrophy$ID_LAB))

# ID_Labs where I have Autotrophy and not Symbio
setdiff(unique (Pocillopora_Photoautotrophy$ID_LAB),unique (Poc_Symbio$ID_LAB))

# ID_Labs where I have both - 87 samples
intersect(unique (Pocillopora_Photoautotrophy$ID_LAB),unique (Poc_Symbio$ID_LAB))

# Keep only those
keep <- intersect(unique (Pocillopora_Photoautotrophy$ID_LAB),unique (Poc_Symbio$ID_LAB))
Poc_Symbio <- Poc_Symbio [Poc_Symbio$ID_LAB %in% keep, ]

### Add here info of sites and depth
Poc_Symbio <- merge (Pocillopora_Photoautotrophy,Poc_Symbio,  by = c("ID_LAB"))
Poc_Symbio <- Poc_Symbio[, -(4:10)]

# This is what you needed

# Measure relative values of Symbio Profiles - must be an easier way but I am tired!

Poc_Symbio$Total <- rowSums(Poc_Symbio[,5:35])

# Create new dataframe with relative values
Relative_Poc_Symbio <- Poc_Symbio[,c("ID_LAB","Depth","Island_Site","Colbar")]
Relative_Poc_Symbio2 <- Poc_Symbio [,5:35] 
# Remove columns (Symbio Profiles) where all are 0
Relative_Poc_Symbio2 <- Relative_Poc_Symbio2[, colSums(Relative_Poc_Symbio2) != 0]

# Combine the initial values with the relatives
Relative_Poc_Symbio <- cbind (Relative_Poc_Symbio,Relative_Poc_Symbio2)

#write.table(Relative_Poc_Symbio, file="/Users/Heloise/Desktop/DEEPHOPE/BIOMOL/NGS/Symportal/20201127_rouze_pocillopora/stats/gonzalo Physio/SymP-GON.csv", sep=";",dec=".", row.names = F)
### ARRET LA

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
  #scale_fill_manual(values= colours) +
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
  #scale_fill_manual(values= colours) +
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

 