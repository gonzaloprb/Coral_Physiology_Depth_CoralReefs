# Script for the Physiology manuscript " Grouping of light across islands " 



#Necessary packages
library(ggplot2); library(RColorBrewer); library(hms); library (tidyverse); library (reshape2); library(ggcorrplot)

# Open all databases first
Moo_Index_Light_loss<- read.csv(file = "Data_Codes/Light/MOO2/MOO-2_Index_Light_loss.csv", header = T, dec = ".", sep = ",")
colnames(Moo_Index_Light_loss)[1] <- "Island_Site"
Moo_Index_Light_loss$Island_Site <- "Moorea"

Bor_Index_Light_loss<- read.csv(file = "Data_Codes/Light/BOR2/BOR-2_Index_Light_loss.csv", header = T, dec = ".", sep = ",")
colnames(Bor_Index_Light_loss)[1] <- "Island_Site"
Bor_Index_Light_loss$Island_Site <- "Bora Bora"

Ran_Index_Light_loss<- read.csv(file = "Data_Codes/Light/RAN3/RAN-3_Index_Light_loss.csv", header = T, dec = ".", sep = ",")
colnames(Ran_Index_Light_loss)[1] <- "Island_Site"
Ran_Index_Light_loss$Island_Site <- "Rangiroa"

Tik_Index_Light_loss<- read.csv(file = "Data_Codes/Light/TIK2/TIK-2_Index_Light_loss.csv", header = T, dec = ".", sep = ",")
colnames(Tik_Index_Light_loss)[1] <- "Island_Site"
Tik_Index_Light_loss$Island_Site <- "Tikehau"

Light_loggers <- rbind (Moo_Index_Light_loss,Bor_Index_Light_loss,Ran_Index_Light_loss,Tik_Index_Light_loss)

Light_loggers = Light_loggers %>% mutate(Island_Site = factor(Island_Site, levels = c("Moorea","Bora Bora","Tikehau","Rangiroa")))


plot_light <- ggplot(data=Light_loggers, aes(y= Depth , x =IndexLoss,color =Island_Site))+
  geom_point (size = 0.8) + 
  #geom_smooth(method = "loess", color="grey", fill="grey", size=0.2, alpha=0.2)+
  geom_path (aes(y= Depth,x =IndexLoss,color=Island_Site),size = 2,width = 1,alpha = 0.5) +
 #  geom_point (aes(y= Depth , x =IndexLoss, size = 0.8)) +
  scale_y_reverse(name ="Depth (m)", limits=c(120,0),breaks = c(120,90, 60, 40, 20, 6,0))+
  scale_x_continuous (position = "top") + 
  labs(x = "Relative_Index_Light (μmol m-2 s-1)", y="Depth (m)")+
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "bottom")
plot_light
ggsave ( "Outputs_R/Figures/Plot_Light_Depth_Sites.pdf", plot_light,width = 4.8, height = 4.5)
# This will become Sup Fig. 1

# Make correlation plot
Light_loggers2 <- dcast(Light_loggers, Depth  ~ Island_Site, mean, add.missing = T)

colnames (Light_loggers2) <- c("Depth (m)", "Moorea", "Bora Bora",  "Rangiroa", "Tikehau")

# Compute a matrix of correlation p-values
p_mat <- cor_pmat(Light_loggers2)
corr <- round(cor(Light_loggers2), 2)
# Visualize the correlation matrix
cor_Depth_Light <- ggcorrplot(corr, method = "square", lab = T, type = "lower", sig.level = 0.05) + scale_x_discrete(position = "top") + theme(axis.text.x = element_text(angle = 90, color = 'black'))
cor_Depth_Light 
ggsave ("Outputs_R/Figures/cor_Depth_Light.pdf", cor_Depth_Light,width = 5, height = 5)
# This will become Fig Sup Fig 1
cor.test (Light_loggers$Depth, Light_loggers$IndexLoss)


