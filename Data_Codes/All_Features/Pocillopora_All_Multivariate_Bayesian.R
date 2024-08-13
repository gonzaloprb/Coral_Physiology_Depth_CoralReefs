# Script to open and analyse all Pocillopora data together. PhotoAutotropy, Morphology and Isotopes
# PCA and Bayesian modelling
# AAA Read Me: "Zoox" and "Polyps" are later replaced with "Symbionts" and "Host" in the Manuscript.

# Necessary packages:
require (dplyr)
require (tidyverse)
library (tidyr); library (plyr); 
require (reshape2)
require (ggridges)
require (PerformanceAnalytics)
library (devtools)
require (FactoMineR); require (Factoshiny)
require (cowplot); require (patchwork)
require (vegan)
require (goeveg)
require (car)
library (ggcorrplot)
library(formattable)
# Bayesian modelling packages!
library (gbm); library (brms); 
library('rstan');library('parallel');library('rstanarm'); 

# Clean previous workspace
rm(list =ls())


############# Open Pocillopora datasets 
# Systematic, I applied the same protocol for every single sample!

Physio_Pocillopora_all <- read.csv ("Data_Codes/All_Features/Physio_Pocillopora_all.csv", header = T, dec = ",", sep = ",")

Pocillopora_Isotopes_Polyps_Zoox <- read.csv ("Data_Codes/Isotopes/Pocillopora_Isotopes_Polyps_Zoox.csv", header = T, dec = ",", sep = ",")

Physio_Pocillopora_all <- merge (Physio_Pocillopora_all, Pocillopora_Isotopes_Polyps_Zoox, by = c("ID_LAB","Island_Site","Depth") , all = T)

# Set quantitative variables as numeric
quan.var <- c("Depth","Relative_Index_Light","Zoox","Chl_A","Chl_C2","RatioChlA","RatioChlC2","Size_Ba","Distance_Ba","Size_Ap","Distance_Ap", 
              "iso13C_Polyps","iso13C_Zoox","iso15N_Polyps","iso15N_Zoox","Ratio_CN_Polyps","Ratio_CN_Zoox","Delta13C_Pol_Zoox","Delta15N_Pol_Zoox")
Physio_Pocillopora_all [quan.var] <- sapply (Physio_Pocillopora_all [quan.var], as.numeric)

# Chlorophyll reading followed Jeffrey & Humphrey (1975) equations
# Chl A and C2 are = (1st Extrac+2nd Extrac)*VolCulotRediss(=5)/VolKeptChl (=1.5)*1/Surface (This measure constant for all) // Chl per surface
# RatioChlA and RatioChlC2 is ChlA / Zoox and C2 / Zoox  // Chl per Zoox

# Add possible variables according to literature:
Physio_Pocillopora_all$Chlorophylls <- Physio_Pocillopora_all$Chl_A + Physio_Pocillopora_all$Chl_C2 # Lesser et al 2010 = Chl a + Chl c2
Physio_Pocillopora_all$Chlorophylls_Ratio <- Physio_Pocillopora_all$Chlorophylls/Physio_Pocillopora_all$Zoox # Lesser et al 2010    =  (Chl a + Chl c2)/Zoox
Physio_Pocillopora_all$Chlorophylls_Ratio_a_c2 <- Physio_Pocillopora_all$Chl_A/Physio_Pocillopora_all$Chl_C2 # Lesser et al 2010    = Chl a/Chl c2
Physio_Pocillopora_all$Chlorophylls_Ratio_c2_a <- Physio_Pocillopora_all$Chl_C2/Physio_Pocillopora_all$Chl_A # Padilla et al 2019   = Chl c2/Chl ca
Physio_Pocillopora_all$Chlorophylls_Ratio_Ratio_c2_a <- Physio_Pocillopora_all$RatioChlC2/Physio_Pocillopora_all$RatioChlA # Padilla et al 2019  = RatioChlC2/RatioChlA  (per number of zoox)

# Change the names of the Sites
Physio_Pocillopora_all$Island_Site <- gsub('Bora_2', 'Bora Bora', Physio_Pocillopora_all$Island_Site)
Physio_Pocillopora_all$Island_Site <- gsub('Moorea_2', 'Moorea', Physio_Pocillopora_all$Island_Site)
Physio_Pocillopora_all$Island_Site <- gsub('Rangiroa_3', 'Rangiroa', Physio_Pocillopora_all$Island_Site)
Physio_Pocillopora_all$Island_Site <- gsub('Tikehau_2', 'Tikehau', Physio_Pocillopora_all$Island_Site)

Physio_Pocillopora_all = Physio_Pocillopora_all %>% mutate(Island_Site = factor(Island_Site, levels = c("Moorea","Bora Bora","Tikehau","Rangiroa")))


## 1st. Quick visual impression with depth. (Figures finally not displayed in the Manuscript)
Zoox_surf <- ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Zoox)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Zoox / surface") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                        axis.text = element_text(size=10, colour="black"),
                                                        axis.title = element_text(size=11, face="bold", colour="black")) 
Chla_surf <- ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Chl_A)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Chl a / surface ") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                           axis.text = element_text(size=10, colour="black"),
                                                           axis.title = element_text(size=11, face="bold", colour="black")) 
Chlc2_surf <- ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Chl_C2)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Chl c2 / surface ") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                           axis.text = element_text(size=10, colour="black"),
                                                           axis.title = element_text(size=11, face="bold", colour="black")) 
Chla_zoox <- ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=RatioChlA)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Chl a / Zoox") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                           axis.text = element_text(size=10, colour="black"),
                                                           axis.title = element_text(size=11, face="bold", colour="black")) 
Chlc2_zoox <- ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=RatioChlC2)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Chl c2 / Zoox") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                           axis.text = element_text(size=10, colour="black"),
                                                           axis.title = element_text(size=11, face="bold", colour="black")) 
Chlac2_surf <- ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Chlorophylls)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Chls a + c2 / surface ") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                           axis.text = element_text(size=10, colour="black"),
                                                           axis.title = element_text(size=11, face="bold", colour="black")) 
Chlac2_Zoox <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Chlorophylls_Ratio)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Chls a + c2 / Zoox ") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                           axis.text = element_text(size=10, colour="black"),
                                                           axis.title = element_text(size=11, face="bold", colour="black")) 
Chl_ratio_a_c2 <- ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Chlorophylls_Ratio_a_c2)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  +
  ylab ("Chlorophylls_Ratio_a_c2_Surface (XX)") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                           axis.text = element_text(size=10, colour="black"),
                                                           axis.title = element_text(size=11, face="bold", colour="black"))
Chl_ratio_c2_a <- ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Chlorophylls_Ratio_c2_a)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Chl c2 / a ") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                                     axis.text = element_text(size=10, colour="black"),
                                                                     axis.title = element_text(size=11, face="bold", colour="black")) 
Chl_ratio_c2_a_zoox <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Chlorophylls_Ratio_Ratio_c2_a)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Chl a / c2 ") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                                     axis.text = element_text(size=10, colour="black"),
                                                                     axis.title = element_text(size=11, face="bold", colour="black")) 
Size_Ba <- ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Size_Ba)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Size basal (mm)") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                                     axis.text = element_text(size=10, colour="black"),
                                                                     axis.title = element_text(size=11, face="bold", colour="black")) 
Dist_Ba <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Distance_Ba)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Distance basal (mm)") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                                     axis.text = element_text(size=10, colour="black"),
                                                                     axis.title = element_text(size=11, face="bold", colour="black")) 
Size_Ap <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Size_Ap)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Size apical (mm)") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                     axis.text = element_text(size=10, colour="black"),
                                                     axis.title = element_text(size=11, face="bold", colour="black")) 
Dist_Ap <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Distance_Ap)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Distance apical (mm)") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                         axis.text = element_text(size=10, colour="black"),
                                                         axis.title = element_text(size=11, face="bold", colour="black")) 
Iso13C_Po <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=iso13C_Polyps)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("δ13C Polyps") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                         axis.text = element_text(size=10, colour="black"),
                                                         axis.title = element_text(size=11, face="bold", colour="black")) 
Iso13C_Zo <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=iso13C_Zoox)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("δ13C Zoox") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                         axis.text = element_text(size=10, colour="black"),
                                                         axis.title = element_text(size=11, face="bold", colour="black")) 
Iso15N_Zo <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=iso15N_Polyps)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("δ15N Polyps") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                         axis.text = element_text(size=10, colour="black"),
                                                         axis.title = element_text(size=11, face="bold", colour="black")) 
Iso15N_Po <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=iso15N_Zoox)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("δ15N Zoox") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                         axis.text = element_text(size=10, colour="black"),
                                                         axis.title = element_text(size=11, face="bold", colour="black")) 
Ratio_CN_Po <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Ratio_CN_Polyps)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("δ13C / δ15N Polyps") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                         axis.text = element_text(size=10, colour="black"),
                                                         axis.title = element_text(size=11, face="bold", colour="black")) 
Ratio_CN_Zo <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Ratio_CN_Zoox)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("δ13C / δ15N Polyps") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                         axis.text = element_text(size=10, colour="black"),
                                                         axis.title = element_text(size=11, face="bold", colour="black")) 
Delta_13C_Pol_Zoo <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Delta13C_Pol_Zoox)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Delta δ13C Polyps - Zoox") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                         axis.text = element_text(size=10, colour="black"),
                                                         axis.title = element_text(size=11, face="bold", colour="black")) 
Delta_15N_Pol_Zoo <-ggplot(Physio_Pocillopora_all, aes(x=factor(Depth), y=Delta15N_Pol_Zoox)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Delta δ15N Polyps - Zoox") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                         axis.text = element_text(size=10, colour="black"),
                                                         axis.title = element_text(size=11, face="bold", colour="black")) 

# Photophysiology
(Fig_Pocillopora_Photo = Zoox_surf + Chla_surf + Chlc2_surf + Chla_zoox + Chlc2_zoox + Chlac2_surf + Chlac2_Zoox + 
     Chl_ratio_c2_a + Chl_ratio_c2_a_zoox + 
    plot_layout(guides = 'collect', ncol = 2)  & theme(legend.position='none'))

# Morphology
(Fig_Pocillopora_Morpho = Size_Ba + Dist_Ba + Size_Ap + Dist_Ap + 
    plot_layout(guides = 'collect', ncol = 2)  & theme(legend.position='none'))

# Isotopes
(Fig_Pocillopora_Iso =  Iso13C_Po + Iso13C_Zo + Iso15N_Po + Iso15N_Zo + Ratio_CN_Po + Ratio_CN_Zo + Delta_13C_Pol_Zoo + Delta_15N_Pol_Zoo + 
    plot_layout(guides = 'collect', ncol = 2)  & theme(legend.position='none'))

# Clean the workspace
rm (Zoox_surf,Chla_surf,Chlc2_surf,Chla_zoox,Chlc2_zoox,Chlac2_surf,Chlac2_Zoox,
      Chl_ratio_c2_a,Chl_ratio_a_c2,Chl_ratio_c2_a_zoox,Size_Ba,Dist_Ba,Size_Ap,Dist_Ap,Iso13C_Po,Iso13C_Zo,Iso15N_Po,Iso15N_Zo,Ratio_CN_Po,Ratio_CN_Zo,Delta_13C_Pol_Zoo,Delta_15N_Pol_Zoo)

## 2nd. Make a huge table resume 

# Necessary to remove the NA
Physio_Pocillopora_2 <- na.omit (Physio_Pocillopora_all)

resume_Pocillopora <-ddply(Physio_Pocillopora_2,~Island_Site + Depth,function(x){c(Zoox=mean(na.omit(x$Zoox)),zoox_sd=sd(na.rm = FALSE,(x$Zoox)),zoox_se=sd(na.rm = FALSE,x$Zoox) / sqrt(length(x$Zoox)),
                                                                                   Chl_A=mean(na.omit(x$Chl_A)),Chl_A_sd=sd(na.rm = FALSE,(x$Chl_A)),Chl_A_se=sd(na.rm = FALSE,x$Chl_A) / sqrt(length(x$Chl_A)),
                                                                                   Chl_C2=mean(na.omit(x$Chl_C2)),Chl_C2_sd=sd(na.rm = FALSE,(x$Chl_C2)),Chl_C2_se=sd(na.rm = FALSE,x$Chl_C2) / sqrt(length(x$Chl_C2)),
                                                                                   RatioChlA=mean(na.omit(x$RatioChlA)),RatioChlA_sd=sd(na.rm = FALSE,(x$RatioChlA)),RatioChlA_se=sd(na.rm = FALSE,x$RatioChlA) / sqrt(length(x$RatioChlA)),
                                                                                   RatioChlC2=mean(na.omit(x$RatioChlC2)),RatioChlC2_sd=sd(na.rm = FALSE,(x$RatioChlC2)),RatioChlC2_se=sd(na.rm = FALSE,x$RatioChlC2) / sqrt(length(x$RatioChlC2)),
                                                                                   Chlorophylls=mean(na.omit(x$Chlorophylls)),Chlorophylls_sd=sd(na.rm = FALSE,(x$Chlorophylls)),Chlorophylls_se=sd(na.rm = FALSE,x$Chlorophylls) / sqrt(length(x$Chlorophylls)),
                                                                                   Chlorophylls_Ratio=mean(na.omit(x$Chlorophylls_Ratio)),Chlorophylls_Ratio_sd=sd(na.rm = FALSE,(x$Chlorophylls_Ratio)),Chlorophylls_Ratio_se=sd(na.rm = FALSE,x$Chlorophylls_Ratio) / sqrt(length(x$Chlorophylls_Ratio)),
                                                                                   Chlorophylls_Ratio_a_c2=mean(na.omit(x$Chlorophylls_Ratio_a_c2)),Chlorophylls_Ratio_a_c2_sd=sd(na.rm = FALSE,(x$Chlorophylls_Ratio_a_c2)),Chlorophylls_Ratio_a_c2_se=sd(na.rm = FALSE,x$Chlorophylls_Ratio_a_c2) / sqrt(length(x$Chlorophylls_Ratio_a_c2)),
                                                                                   Chlorophylls_Ratio_c2_a=mean(na.omit(x$Chlorophylls_Ratio_c2_a)),Chlorophylls_Ratio_c2_a_sd=sd(na.rm = FALSE,(x$Chlorophylls_Ratio_c2_a)),Chlorophylls_Ratio_c2_a_se=sd(na.rm = FALSE,x$Chlorophylls_Ratio_c2_a) / sqrt(length(x$Chlorophylls_Ratio_c2_a)),
                                                                                   Chlorophylls_Ratio_Ratio_c2_a=mean(na.omit(x$Chlorophylls_Ratio_Ratio_c2_a)),Chlorophylls_Ratio_Ratio_c2_a_sd=sd(na.rm = FALSE,(x$Chlorophylls_Ratio_Ratio_c2_a)),Chlorophylls_Ratio_Ratio_c2_a_se=sd(na.rm = FALSE,x$Chlorophylls_Ratio_Ratio_c2_a) / sqrt(length(x$Chlorophylls_Ratio_Ratio_c2_a)),
                                                                                   Size_Ba=mean(na.omit(x$Size_Ba)),Size_Ba_sd=sd(na.rm = FALSE,(x$Size_Ba)),Size_Ba_se=sd(na.rm = FALSE,x$Size_Ba) / sqrt(length(x$Size_Ba)),
                                                                                   Distance_Ba=mean(na.omit(x$Distance_Ba)),Distance_Ba_sd=sd(na.rm = FALSE,(x$Distance_Ba)),Distance_Ba_se=sd(na.rm = FALSE,x$Distance_Ba) / sqrt(length(x$Distance_Ba)),
                                                                                   Size_Ap=mean(na.omit(x$Size_Ap)),Size_Ap_sd=sd(na.rm = FALSE,(x$Size_Ap)),Size_Ap_se=sd(na.rm = FALSE,x$Size_Ap) / sqrt(length(x$Size_Ap)),
                                                                                   Distance_Ap=mean(na.omit(x$Distance_Ap)),Distance_Ap_sd=sd(na.rm = FALSE,(x$Distance_Ap)),Distance_Ap_se=sd(na.rm = FALSE,x$Distance_Ap) / sqrt(length(x$Distance_Ap)),
                                                                                   iso13C_Polyps=mean(na.omit(x$iso13C_Polyps)),iso13C_Polyps_sd=sd(na.rm = FALSE,(x$iso13C_Polyps)),iso13C_Polyps_se=sd(na.rm = FALSE,x$iso13C_Polyps) / sqrt(length(x$iso13C_Polyps)),
                                                                                   iso13C_Zoox=mean(na.omit(x$iso13C_Zoox)),iso13C_Zoox_sd=sd(na.rm = FALSE,(x$iso13C_Zoox)),iso13C_Zoox_se=sd(na.rm = FALSE,x$iso13C_Zoox) / sqrt(length(x$iso13C_Zoox)),
                                                                                   iso15N_Polyps=mean(na.omit(x$iso15N_Polyps)),iso15N_Polyps_sd=sd(na.rm = FALSE,(x$iso15N_Polyps)),iso15N_Polyps_se=sd(na.rm = FALSE,x$iso15N_Polyps) / sqrt(length(x$iso15N_Polyps)),
                                                                                   iso15N_Zoox=mean(na.omit(x$iso15N_Zoox)),iso15N_Zoox_sd=sd(na.rm = FALSE,(x$iso15N_Zoox)),iso15N_Zoox_se=sd(na.rm = FALSE,x$iso15N_Zoox) / sqrt(length(x$iso15N_Zoox)),
                                                                                   Ratio_CN_Polyps=mean(na.omit(x$Ratio_CN_Polyps)),Ratio_CN_Polyps_sd=sd(na.rm = FALSE,(x$Ratio_CN_Polyps)),Ratio_CN_Polyps_se=sd(na.rm = FALSE,x$Ratio_CN_Polyps) / sqrt(length(x$Ratio_CN_Polyps)),
                                                                                   Ratio_CN_Zoox=mean(na.omit(x$Ratio_CN_Zoox)),Ratio_CN_Zoox_sd=sd(na.rm = FALSE,(x$Ratio_CN_Zoox)),Ratio_CN_Zoox_se=sd(na.rm = FALSE,x$Ratio_CN_Zoox) / sqrt(length(x$Ratio_CN_Zoox)),
                                                                                   Delta13C_Pol_Zoox=mean(na.omit(x$Delta13C_Pol_Zoox)),Delta13C_Pol_Zoox_sd=sd(na.rm = FALSE,(x$Delta13C_Pol_Zoox)),Delta13C_Pol_Zoox_se=sd(na.rm = FALSE,x$Delta13C_Pol_Zoox) / sqrt(length(x$Delta13C_Pol_Zoox)),
                                                                                   Delta15N_Pol_Zoox=mean(na.omit(x$Delta15N_Pol_Zoox)),Delta15N_Pol_Zoox_sd=sd(na.rm = FALSE,(x$Delta15N_Pol_Zoox)),Delta15N_Pol_Zoox_se=sd(na.rm = FALSE,x$Delta15N_Pol_Zoox) / sqrt(length(x$Delta15N_Pol_Zoox))
                                                                                   )})


resume_Pocillopora <- resume_Pocillopora  %>% select(Island_Site, Depth, Zoox,Chl_A,Chl_C2,RatioChlA,RatioChlC2,Chlorophylls,Chlorophylls_Ratio,
                                                     Chlorophylls_Ratio_a_c2,Chlorophylls_Ratio_c2_a,Chlorophylls_Ratio_Ratio_c2_a,Size_Ba,Distance_Ba,
                                                     Size_Ap,Distance_Ap,iso13C_Polyps,iso13C_Zoox,iso15N_Polyps,iso15N_Zoox,Ratio_CN_Polyps,Ratio_CN_Zoox,
                                                     Delta13C_Pol_Zoox,Delta15N_Pol_Zoox)

resume_Pocillopora_table <- melt (resume_Pocillopora, id = c ("Island_Site", "Depth"), measure.vars = c(3:24), value.name = c("Value"))
resume_Pocillopora_table <- dcast(resume_Pocillopora_table, variable ~ Island_Site + Depth , mean, add.missing = T)

# View (resume_Pocillopora_table)
write.table(format(resume_Pocillopora_table,digits = 2), file = "Outputs_R/Tables/resume_Pocillopora_table.csv", sep = ",", col.names = NA,qmethod = "double")
# This is Supplementary Table 3



## 3rd. Standardize the numeric data to be able to compare ####
df_1 <- Physio_Pocillopora_all
df_1 <- na.omit(df_1)

df_1$Zoox <- (df_1$Zoox - mean(df_1$Zoox)) / sd(df_1$Zoox)
df_1$Chl_A <- (df_1$Chl_A - mean(df_1$Chl_A)) / sd(df_1$Chl_A)
df_1$Chl_C2 <- (df_1$Chl_C2 - mean(df_1$Chl_C2)) / sd(df_1$Chl_C2)

df_1$RatioChlA <- (df_1$RatioChlA - mean(df_1$RatioChlA)) / sd(df_1$RatioChlA)
df_1$RatioChlC2 <- (df_1$RatioChlC2 - mean(df_1$RatioChlC2)) / sd(df_1$RatioChlC2)

df_1$Chlorophylls <- (df_1$Chlorophylls - mean(df_1$Chlorophylls)) / sd(df_1$Chlorophylls)
df_1$Chlorophylls_Ratio <- (df_1$Chlorophylls_Ratio - mean(df_1$Chlorophylls_Ratio)) / sd(df_1$Chlorophylls_Ratio)
# 
df_1$Chlorophylls_Ratio_a_c2 <- (df_1$Chlorophylls_Ratio_a_c2 - mean(df_1$Chlorophylls_Ratio_a_c2)) / sd(df_1$Chlorophylls_Ratio_a_c2)
df_1$Chlorophylls_Ratio_c2_a <- (df_1$Chlorophylls_Ratio_c2_a - mean(df_1$Chlorophylls_Ratio_c2_a)) / sd(df_1$Chlorophylls_Ratio_c2_a)
df_1$Chlorophylls_Ratio_Ratio_c2_a <- (df_1$Chlorophylls_Ratio_Ratio_c2_a - mean(df_1$Chlorophylls_Ratio_Ratio_c2_a)) / sd(df_1$Chlorophylls_Ratio_Ratio_c2_a)


df_1$Size_Ba <- (df_1$Size_Ba - mean(df_1$Size_Ba)) / sd(df_1$Size_Ba)
df_1$Distance_Ba <- (df_1$Distance_Ba - mean(df_1$Distance_Ba)) / sd(df_1$Distance_Ba)
df_1$Size_Ap <- (df_1$Size_Ap - mean(df_1$Size_Ap)) / sd(df_1$Size_Ap)
df_1$Distance_Ap <- (df_1$Distance_Ap - mean(df_1$Distance_Ap)) / sd(df_1$Distance_Ap)

df_1$iso13C_Polyps <- (df_1$iso13C_Polyps - mean(df_1$iso13C_Polyps)) / sd(df_1$iso13C_Polyps)
df_1$iso13C_Zoox <- (df_1$iso13C_Zoox - mean(df_1$iso13C_Zoox)) / sd(df_1$iso13C_Zoox)
df_1$iso15N_Polyps <- (df_1$iso15N_Polyps - mean(df_1$iso15N_Polyps)) / sd(df_1$iso15N_Polyps)
df_1$iso15N_Zoox <- (df_1$iso15N_Zoox - mean(df_1$iso15N_Zoox)) / sd(df_1$iso15N_Zoox)

df_1$Ratio_CN_Polyps <- (df_1$Ratio_CN_Polyps - mean(df_1$Ratio_CN_Polyps)) / sd(df_1$Ratio_CN_Polyps)
df_1$Ratio_CN_Zoox <- (df_1$Ratio_CN_Zoox - mean(df_1$Ratio_CN_Zoox)) / sd(df_1$Ratio_CN_Zoox)
df_1$Delta13C_Pol_Zoox <- (df_1$Delta13C_Pol_Zoox - mean(df_1$Delta13C_Pol_Zoox)) / sd(df_1$Delta13C_Pol_Zoox)
df_1$Delta15N_Pol_Zoox <- (df_1$Delta15N_Pol_Zoox - mean(df_1$Delta15N_Pol_Zoox)) / sd(df_1$Delta15N_Pol_Zoox)



## 4th. Check correlations - colinearity
df_2 <- na.omit(df_1)

# Compute a matrix of correlation p-values
names (df_2)

# Drop non-numerical columns
df_2 <- df_2[,-c(1,2,4,15:18)]

names (df_2)
colnames (df_2) [1] <-  "Depth"
colnames (df_2) [2] <-  "Zoox / surface"
colnames (df_2) [3] <-  "Chl a / surface"
colnames (df_2) [4] <-"Chl c2 / surface"
colnames (df_2) [5] <- "Light"
colnames (df_2) [6] <-  "Chl a / Zoox"
colnames (df_2) [7] <-"Chl c2 / Zoox"
colnames (df_2) [8] <-  "Size basal (mm)"
colnames (df_2) [9] <-  "Distance basal (mm)"
colnames (df_2) [10] <-"Size apical (mm)"
colnames (df_2) [11] <-  "Distance apical (mm)"

colnames (df_2) [12] <-  "δ13C Polyps"
colnames (df_2) [13] <-  "δ13C Zoox"
colnames (df_2) [14] <-"δ15N Polyps"
colnames (df_2) [15] <-  "δ15N Zoox"
colnames (df_2) [16] <-  "δ13C / δ15N Polyps"
colnames (df_2) [17] <-"δ13C / δ15N Zoox"
colnames (df_2) [18] <-  "Delta δ13C Polyps - δ13C Zoox"
colnames (df_2) [19] <-  "Delta δ15N Polyps - δ15N Zoox"
colnames (df_2) [20] <-"Chl a + Chl c2 / surface"
colnames (df_2) [21] <-  "Chl a + Chl c2 / Zoox"
colnames (df_2) [22] <-  "Chl a / c2"
colnames (df_2) [23] <- "Chl c2 / a"
colnames (df_2) [24] <-  "Chl c2 / a / Zoox"


# Change order of columns in r
df_2 <- df_2[, c(1,5,2,3,4,6,7,20,21,22,23,24,8,9,10,11,12,13,14,15,16,17,18,19)] 

p_mat <- cor_pmat(df_2)
corr <- round(cor(df_2),2)
# Visualize the correlation matrix
cor_Poc <- ggcorrplot(corr, lab = T, type = "upper") + scale_x_discrete(position = "top") + theme(axis.text.x = element_text(angle = 90, color = 'black'))
cor_Poc
ggsave ("Outputs_R/Figures/cor_Poc_all.pdf", cor_Poc, width = 14, height = 14) 
# This will become Supplementary Figure 4. Changed the names from "Zoox" and "Polyps" to "Symbionts" and "Host"



## 5th. Principal Component Analysis

# PCA without dropping correlated variables
# 2 Qualitative variables: Island Site and Depth

df_1_2 <- df_1

# Unite multiple columns with tidyr
df_1_2 <- df_1_2 %>% unite("Unique", Island_Site:ID_LAB, sep = "_",remove = FALSE)
df_1_2$Unique <- with(df_1_2, paste0(Depth, sep = "_",Unique))
rownames (df_1_2) <- df_1_2$Unique
# Remove unnecessary variables
df_1_2 <- df_1_2[,-c(1:3,5,16:19)]

res.pca <- PCA(df_1_2, quanti.sup = c(1,5))
barplot(res.pca$eig[,1],main="Eigenvalues",names.arg=1:nrow(res.pca$eig))
summary(res.pca)

# Make plot of vars and ind
res.pca <- PCA(df_1_2, quanti.sup = c(1))

# Make depth as factor 
df_1_2$Depth_fact <- as.factor (df_1_2$Depth)
# result <- Factoshiny(df_1_2)

# Supplementary variables have no influence on the principal components of the analysis. 
# I display them in blue to help with the analysis
res.PCA<-PCA(df_1_2,quanti.sup = c(1),quali.sup=c(25),graph=FALSE)
plot.PCA(res.PCA,choix='var',title="PCA graph of variables - Pocillopora",col.quanti.sup='#0000FF')

plotellipses(res.PCA, keepvar=25,invisible=c('quali','ind.sup'),cex=0.75,cex.main=0.75,cex.axis=0.75,title="PCA graph of individuals - Pocillopora",label ='none')
plotellipses(res.PCA, keepvar=25,invisible=c('quali','ind.sup'),cex=0.75,cex.main=0.75,cex.axis=0.75,title="PCA graph of individuals - Pocillopora",label ='ind')


# plot with ggplot
PCA_All <- as.data.frame(res.PCA$ind$coord)
PCA_All$ID <- rownames (PCA_All)
# Extract depths
PCA_All$Depth <- sub("\\_.*", "", PCA_All$ID)
# Extract Island 
PCA_All$Island <- sapply(strsplit(PCA_All$ID, "_"), function(x) x[2])

PCA_All = PCA_All %>% mutate(Depth = factor(Depth, levels = c("6", "20", "40", "60")))
PCA_All = PCA_All %>% mutate(Island = factor(Island, levels = c("Moorea","Bora Bora","Tikehau","Rangiroa")))

PCA_All %>% ggplot() + 
  geom_point(aes(x = Dim.1, y = Dim.2, color = Depth, fill = Depth, shape = Island),  show.legend = T) +
  geom_vline(xintercept=0, linetype="dotted", color = "black") +
  geom_hline(yintercept=0, linetype="dotted", color = "black") +
 # stat_ellipse(aes(x = Dim.1, y = Dim.2, color = Depth, fill = Depth), type = "norm", linetype = 3, level = 0.95) + 
  stat_ellipse(aes(x = Dim.1, y = Dim.2, color = Depth, fill = Depth), type = "norm", linetype = 1, level = 0.5) + 
  #  stat_ellipse(aes(x = Dim.1, y = Dim.2, color = Depth, fill = Depth), type = "euclid", level = 3) +  
  # scale_x_continuous(name ="Dim 1 (74.93%)", limits=c(-3,3), breaks = c(-2,0,2)) +
  # scale_y_continuous(name ="Dim 2 (25.07%)", limits=c(-2,4), breaks = c(-2,0,2)) +
  labs(x = paste("Dim 1:",round(res.PCA$eig[1,2],digits = 2),"(%)"), y = paste("Dim 2:",round(res.PCA$eig[2,2],digits = 2),"(%)"), title = "") + theme_classic() + theme (legend.position = "right")
# This plot is not used in the MS, nor Sup Fig.


# Drop variables with correlations higher than 0.85

# Depth and index Light very correlated - keep only depth
# Chl A and Chl C2 very correlated -  Keep only Chlorophylls
# Ratio A and Ratio C2 very correlated - Keep only Chlorophylls ratio
# Chlorophylls ratio A_C2 and C2_A very correlated - keep only A_C2 - not even delete because -0.78
# Iso 13 C Zoox iso 13 C Polyps 0.95 correlated! But, in this case, keeping both!
# Delta15N Pol Zoox and Iso 15 N Polyps very correlated! But, in this case, keeping both!



# PCA without correlated variables
# Unite multiple columns with tidyr
df_2 <- df_1 %>% unite("Unique", Island_Site:ID_LAB, sep = "_",remove = FALSE)
df_2$Unique <- with(df_2, paste0(Depth, sep = "_",Unique))
rownames (df_2) <- df_2$Unique

# Remove unnecessary correlated variables 
keep_columns <- c("Depth","Chlorophylls", "Chlorophylls_Ratio", "Chlorophylls_Ratio_Ratio_c2_a","Size_Ba", "Distance_Ba", "Size_Ap","Distance_Ap",
                  "iso13C_Polyps","iso13C_Zoox","iso15N_Polyps","iso15N_Zoox","Ratio_CN_Polyps","Ratio_CN_Zoox","Delta13C_Pol_Zoox","Delta15N_Pol_Zoox") # The best is not keeping "Chlorophylls_Ratio"
df_2_2 <- df_2[,c(keep_columns)]


# PCA Without the data standardised 
df_3 <- Physio_Pocillopora_2 %>% unite("Unique", Island_Site:ID_LAB, sep = "_",remove = FALSE)
df_3$Unique <- with(df_3, paste0(Depth, sep = "_",Unique))
rownames (df_3) <- df_3$Unique

# Remove unnecessary columns - considering the ones that are highly correlated between them
keep_columns <- c("Depth","Chlorophylls", "Chlorophylls_Ratio", "Chlorophylls_Ratio_Ratio_c2_a","Size_Ba", "Distance_Ba", "Size_Ap","Distance_Ap",
                  "iso13C_Polyps","iso13C_Zoox","iso15N_Polyps","iso15N_Zoox","Ratio_CN_Polyps","Ratio_CN_Zoox","Delta13C_Pol_Zoox","Delta15N_Pol_Zoox") # The best is not keeping "Chlorophylls_Ratio"
df_3_2 <- df_3[,c(keep_columns)]



# Factoshinny - PCA - With not standardised data, necessary for PCA
res.pca <- PCA(df_3_2, quanti.sup = c(1))

# Make depth as factor 
df_3_2$Depth_fact <- as.factor (df_3_2$Depth)
# result <- Factoshiny(df_2_2)

# Supplementary variables have no influence on the principal components of the analysis. (So Depth is not inside)
# I display "Depth" in blue to help with the analysis
res.PCA<-PCA(df_3_2,quanti.sup = c(1),quali.sup=c(17),graph=FALSE)
Plot_vars_Poc <- plot.PCA(res.PCA,choix='var',title="",col.quanti.sup='#0000FF', cex=0.8,graph.type = "ggplot" ) + theme_void()
Plot_vars_Poc 
ggsave ("Outputs_R/Figures/Plot_vars_Poc.pdf", Plot_vars_Poc, width = 4, height = 4)
# This will become Fig 3a, PCA; Changed the names from "Zoox" and "Polyps" to "Symbionts" and "Host"

# Extra plot for visualisation
plotellipses(res.PCA, keepvar=17,invisible=c('quali','ind.sup'),cex=0.75,cex.main=0.75,cex.axis=0.75,title="PCA graph of individuals - Pocillopora",label ='ind')


# PCA plot with ggplot
PCA_all2 <- as.data.frame(res.PCA$ind$coord) # all2 because correlated vars are deleted
PCA_all2$ID <- rownames (PCA_all2)
# Extract depths
PCA_all2$Depth <- sub("\\_.*", "", PCA_all2$ID)
# Extract Island 
PCA_all2$Island <- sapply(strsplit(PCA_all2$ID, "_"), function(x) x[2])

PCA_all2 = PCA_all2 %>% mutate(Depth = factor(Depth, levels = c("6", "20", "40", "60")))
PCA_all2 = PCA_all2 %>% mutate(Island = factor(Island, levels = c("Moorea","Bora Bora","Tikehau","Rangiroa")))

cols <- c("6" = "#BDD7E7", "20" = "#6BAED6", "40" = "#3182BD", "60" = "#08519C", "90" = "darkblue")

PCA_Plot_all2 <- PCA_all2 %>% ggplot() + 
  geom_point(aes(x = Dim.1, y = Dim.2, color = Depth, fill = Depth, shape = Island),  show.legend = T) +
  geom_vline(xintercept=0, linetype="dotted", color = "black") +
  geom_hline(yintercept=0, linetype="dotted", color = "black") +
  scale_colour_manual(values = cols) +
 # stat_ellipse(aes(x = Dim.1, y = Dim.2, color = Depth, fill = Depth), type = "norm", linetype = 3, level = 0.95) + 
  stat_ellipse(aes(x = Dim.1, y = Dim.2, color = Depth, fill = Depth), type = "norm", linetype = 1, level = 0.5) + 
  coord_fixed(ratio = 1) +
#  stat_ellipse(aes(x = Dim.1, y = Dim.2, color = Depth, fill = Depth), type = "euclid", level = 3) +  
  # scale_x_continuous(name ="Dim 1 (74.93%)", limits=c(-3,3), breaks = c(-2,0,2)) +
  # scale_y_continuous(name ="Dim 2 (25.07%)", limits=c(-2,4), breaks = c(-2,0,2)) +
  labs(x = paste("Dim 1:",round(res.PCA$eig[1,2],digits = 2),"(%)"), y = paste("Dim 2:",round(res.PCA$eig[2,2],digits = 2),"(%)"), title = "") + theme_classic() + theme (legend.position = "right") + theme(aspect.ratio=4/4)
PCA_Plot_all2
ggsave ("Outputs_R/Figures/PCA_Plot_Poc.pdf", PCA_Plot_all2, width = 6, height = 6)
# This will become Fig 3a, left. 


## 6th. Pariwise multivariate analysis between depths for each island with standardised data 
df_2_2$Island <- rownames (df_2_2)
# Extract Island 
df_2_2$Island <- sapply(strsplit(df_2_2$Island, "_"), function(x) x[2])

# Two way with interaction anova. Type 3 because unbalanced data
Anova(lm(Chlorophylls ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Chlorophylls_Ratio ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Chlorophylls_Ratio_Ratio_c2_a ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)

Anova(lm(Distance_Ba ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Size_Ap ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Distance_Ap ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Size_Ba ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)

Anova(lm(iso13C_Polyps ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(iso13C_Zoox ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(iso15N_Polyps ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(iso15N_Zoox ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Ratio_CN_Polyps ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Ratio_CN_Zoox ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Delta13C_Pol_Zoox ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Delta15N_Pol_Zoox ~ Depth*Island, data=df_2_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)

# Pariwise multivariate analysis between depths for each island with NON standardised data 
df_3_2$Island <- rownames (df_3_2)
# Extract Island 
df_3_2$Island <- sapply(strsplit(df_3_2$Island, "_"), function(x) x[2])

# Two way with interaction anova. Type 3 because unbalanced data
Anova(lm(Chlorophylls ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Chlorophylls_Ratio_Ratio_c2_a ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)

Anova(lm(Distance_Ba ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Size_Ap ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Distance_Ap ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Size_Ba ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)

Anova(lm(iso13C_Polyps ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(iso13C_Zoox ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(iso15N_Polyps ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(iso15N_Zoox ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Ratio_CN_Polyps ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Ratio_CN_Zoox ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Delta13C_Pol_Zoox ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)
Anova(lm(Delta15N_Pol_Zoox ~ Depth*Island, data=df_3_2, contrasts=list(topic=contr.sum, sys=contr.sum)), type=3)



# Extra statistics considering composition of all measures at once, with standardised data
total <- as.data.frame (df_2_2[,c(2:16)]) 
# To avoid negative values ... 
total <- total + abs (total)

dis <- vegdist (total, method="bray")
total_NMDS <- metaMDS(total, k=2) 

# If interaction is not significant it means there is not community shift.
# Every Island in function of Depth
mod1 <- adonis2(dis ~ Island*Depth, data = df_2_2, permutations = 999)
mod1
summary (mod1)
# Every Depth in function of Island
mod2 <- adonis2(dis ~ Depth*Island, data = df_2_2, permutations = 999)  # To check if sensitive to the order
mod2
summary (mod2)
# Same results, independently of the order
adonis(dis ~ Depth / Island, strata = df_2_2$Depth, by = margin, data = df_2_2, permutations = 999)
adonis(dis ~ Island / Depth, strata = df_2_2$Island, by = margin, data = df_2_2, permutations = 999)
# Overall depth is playing a significant effect in the differentiation

# Extra comparisons between depths, Island and the combination of Depth & Island 
source ("Outputs_R/Files/pairwise.adonis.R")
pairwise.adonis(df_2_2[,2:16],paste(df_2_2$Depth,df_2_2$Island))
pairwise.adonis(df_2_2[,2:16],paste(df_2_2$Island,df_2_2$Depth))
pairwise.adonis(df_2_2[,2:16],paste(df_2_2$Depth))

# Separating by Islands:
df_2_2_Moorea<- subset (df_2_2, Island == c("Moorea"))
pairwise.adonis(df_2_2_Moorea[,2:16],paste(df_2_2_Moorea$Depth)) # very few replicates

df_2_2_Bora<- subset (df_2_2, Island == c("Bora Bora"))
pairwise.adonis(df_2_2_Bora[,2:16],paste(df_2_2_Bora$Depth))

df_2_2_Tikehau<- subset (df_2_2, Island == c("Tikehau"))
pairwise.adonis(df_2_2_Tikehau[,2:16],paste(df_2_2_Tikehau$Depth))

df_2_2_Rangiroa<- subset (df_2_2, Island == c("Rangiroa"))
pairwise.adonis(df_2_2_Rangiroa[,2:16],paste(df_2_2_Rangiroa$Depth))
# Pairwise multivariate analysis between depths for each island 
# Instead of making individual comparisons, let's study the variables response to increasing depth


## 7th. Bayesian modelling
df_2_2$Site <- word(rownames (df_2_2), 2, sep = "_")

hist (df_2_2$Chlorophylls)
hist (df_2_2$Chlorophylls_Ratio)
hist (df_2_2$Distance_Ba)
hist (df_2_2$iso13C_Zoox)

# Bayesian modelling. 
# I tried several different combinations, but I am just showing the one used and two more!


# All variables at once, with random slope - This one does not work. Some variables are too correlated to allow convergence!
# fit_Bayes_all_Random_Slope1_Pocillopora <- brm(
#      mvbind(Chlorophylls, Chlorophylls_Ratio,Chlorophylls_Ratio_Ratio_c2_a,Size_Ba,Distance_Ba,Size_Ap,Distance_Ap,iso13C_Polyps,iso13C_Zoox,iso15N_Polyps,
#             iso15N_Zoox,Ratio_CN_Polyps,Ratio_CN_Zoox,Delta13C_Pol_Zoox,Delta15N_Pol_Zoox) ~ 1 + Depth + (1 + Depth | Island),
#      data = df_2_2, family = student(), control = list(adapt_delta = 0.95, max_treedepth = 15), iter =4000, chains = 2, cores = 2)
# save(fit_Bayes_all_Random_Slope1_Pocillopora, file="Outputs_R/Files/fit_Bayes_all_Random_Slope1_Pocillopora.RData")
# load("Outputs_R/Files/fit_Bayes_all_Random_Slope1_Pocillopora.RData") 
# summary(fit_Bayes_all_Random_Slope1_Pocillopora)


# Removing some variables: "Delta15N_Pol_Zoox" and "iso13C_Zoox" But keeping both!
# fit_Bayes_all_Random_Slope_2 <- brm(
#   mvbind(Chlorophylls, Chlorophylls_Ratio,Chlorophylls_Ratio_Ratio_c2_a,Size_Ba,Distance_Ba,Size_Ap,Distance_Ap,iso13C_Polyps,iso15N_Polyps,
#          iso15N_Zoox,Ratio_CN_Polyps,Ratio_CN_Zoox,Delta13C_Pol_Zoox) ~ 1 + Depth + (1 + Depth | Island),
#   data = df_2_2, family = student(), control = list(adapt_delta = 0.95, max_treedepth = 15), iter =4000, chains = 2, cores = 2)
# 
# save(fit_Bayes_all_Random_Slope_2, file="Outputs_R/Files/fit_Bayes_all_Random_Slope2_Pocillopora.RData")
load("Outputs_R/Files/fit_Bayes_all_Random_Slope2_Pocillopora.RData") 
summary(fit_Bayes_all_Random_Slope_2) # This model converged!

# Without random slope to check - we prefer to keep the random slope for site also
# Delta15N Pol Zoox and Iso 15 N Polyps very correlated! But keeping both!
# fit_Bayes_all_Poc <- brm(
#   mvbind(Chlorophylls, Chlorophylls_Ratio,Chlorophylls_Ratio_Ratio_c2_a,Size_Ba,Distance_Ba,Size_Ap,Distance_Ap,iso13C_Polyps,iso15N_Polyps,
#          iso15N_Zoox,Ratio_CN_Polyps,Ratio_CN_Zoox,Delta13C_Pol_Zoox) ~  Depth + (1|Site),
#   data = df_2_2, family = student(), control = list(adapt_delta = 0.95, max_treedepth = 15), iter =4000, chains = 2, cores = 2)
# save(fit_Bayes_all_Poc, file="Outputs_R/Files/fit_Bayes_all_Poc.RData")
# load("Outputs_R/Files/fit_Bayes_all_Poc.RData") 
# fit_Bayes_Pocillopora <- fit_Bayes_all_Poc


### The good model is fit_Bayes_all_Random_Slope_2 ###
fit_Bayes_Pocillopora <- fit_Bayes_all_Random_Slope_2


fit_Bayes_Pocillopora <- add_criterion(fit_Bayes_Pocillopora, "loo")
summary(fit_Bayes_Pocillopora)
bayes_R2(fit_Bayes_Pocillopora)
# coef(fit_Bayes_Pocillopora)
# Variance of the Random effect
# ranef(fit_Bayes_Pocillopora)
# VarCorr(fit_Bayes_Pocillopora)

# pp_check(fit_Bayes_Pocillopora, resp = "Chlorophylls")
# pp_check(fit_Bayes_Pocillopora, resp = "ChlorophyllsRatio")
# pp_check(fit_Bayes_Pocillopora, resp = "iso13CPolyps")

# Plot the conditional effects
ce_fit_Bayes_Poc <- conditional_effects(fit_Bayes_Pocillopora, nsamples = 200, spaghetti = F, probs = c(0.025, 0.975))
plot(ce_fit_Bayes_Poc, ask = FALSE, points = F) # Probability scale!
# Some examples with nicer ggplots
plot(ce_fit_Bayes_Poc, plot = FALSE)[["Chlorophylls.Chlorophylls_Depth"]] + 
  labs(y = "Chlorophylls A & C2", x = "Depth (m)") + theme_classic()
plot(ce_fit_Bayes_Poc, plot = FALSE)[["ChlorophyllsRatio.ChlorophyllsRatio_Depth"]] + 
  labs(y = "Chlorophylls (A & C2) / Zoox", x = "Depth (m)") + theme_classic()
plot(ce_fit_Bayes_Poc, plot = FALSE)[["iso13CPolyps.iso13CPolyps_Depth"]] + 
  labs(y = "iso 13 C Polyps & Zoox", x = "Depth (m)") + theme_classic()



# Plot the slopes of the Fixed effects
Fixed_Bayes_Poc_FE <- as.data.frame (fixef(fit_Bayes_Pocillopora))
Fixed_Bayes_Poc_FE <- Fixed_Bayes_Poc_FE[c("Chlorophylls_Depth","ChlorophyllsRatio_Depth","ChlorophyllsRatioRatioc2a_Depth",
                                               "SizeAp_Depth","SizeBa_Depth", "DistanceAp_Depth","DistanceBa_Depth", 
                                           "iso13CPolyps_Depth", "iso15NPolyps_Depth","iso15NZoox_Depth","RatioCNPolyps_Depth","RatioCNZoox_Depth","Delta13CPolZoox_Depth"),]


Fixed_Bayes_Poc_FE$RespVar <- rownames(Fixed_Bayes_Poc_FE)
Fixed_Bayes_Poc_FE$RespVar <- gsub('Chlorophylls_Depth', 'Chlorophylls a+c2 surface', Fixed_Bayes_Poc_FE$RespVar)
Fixed_Bayes_Poc_FE$RespVar <- gsub('ChlorophyllsRatio_Depth', 'Chlorophyll ratio Zoox', Fixed_Bayes_Poc_FE$RespVar)
Fixed_Bayes_Poc_FE$RespVar <- gsub('ChlorophyllsRatioRatioc2a_Depth', 'Chlorophyll ratio c2:a', Fixed_Bayes_Poc_FE$RespVar)

Fixed_Bayes_Poc_FE$RespVar <- gsub('SizeAp_Depth', 'Size apical', Fixed_Bayes_Poc_FE$RespVar)
Fixed_Bayes_Poc_FE$RespVar <- gsub('SizeBa_Depth', 'Size basal', Fixed_Bayes_Poc_FE$RespVar)
Fixed_Bayes_Poc_FE$RespVar <- gsub('DistanceAp_Depth', 'Distance apical', Fixed_Bayes_Poc_FE$RespVar)
Fixed_Bayes_Poc_FE$RespVar <- gsub('DistanceBa_Depth', 'Distance basal', Fixed_Bayes_Poc_FE$RespVar)

Fixed_Bayes_Poc_FE$RespVar <- gsub('iso13CPolyps_Depth', 'iso13C Polyps & Zoox', Fixed_Bayes_Poc_FE$RespVar)
Fixed_Bayes_Poc_FE$RespVar <- gsub('iso15NPolyps_Depth', 'iso15N Polyps & Delta 15N Polyps - Zoox', Fixed_Bayes_Poc_FE$RespVar)
Fixed_Bayes_Poc_FE$RespVar <- gsub('iso15NZoox_Depth', 'iso15N Zoox', Fixed_Bayes_Poc_FE$RespVar)
Fixed_Bayes_Poc_FE$RespVar <- gsub('RatioCNPolyps_Depth', 'Ratio C:N Polyps', Fixed_Bayes_Poc_FE$RespVar)
Fixed_Bayes_Poc_FE$RespVar <- gsub('RatioCNZoox_Depth', 'Ratio C:N Zoox', Fixed_Bayes_Poc_FE$RespVar)
Fixed_Bayes_Poc_FE$RespVar <- gsub('Delta13CPolZoox_Depth', 'Delta 13C Polyps - Zoox', Fixed_Bayes_Poc_FE$RespVar)


Fixed_Bayes_Poc_FE$RespVar =  factor(Fixed_Bayes_Poc_FE$RespVar,levels = c ("iso13C Polyps & Zoox", "Ratio C:N Zoox", "Chlorophyll ratio c2:a","iso15N Zoox",
                                                                                "iso15N Polyps & Delta 15N Polyps - Zoox","Ratio C:N Polyps","Distance apical",
                                                                                "Chlorophyll ratio Zoox","Size basal","Delta 13C Polyps - Zoox","Size apical","Chlorophylls a+c2 surface","Distance basal"))

# Format the Quantiles to 3 decimals
Fixed_Bayes_Poc_FE$Q2.5 <- formattable(Fixed_Bayes_Poc_FE$Q2.5, digits = 3, format = "f")
Fixed_Bayes_Poc_FE$Q2.5[Fixed_Bayes_Poc_FE$Q2.5 == -0.0029085549] <- 0.00

# With two decimals they are significant!
Fixed_Bayes_Poc_FE$Q2.5 [Fixed_Bayes_Poc_FE$Q2.5  > -0.005 & Fixed_Bayes_Poc_FE$Q2.5  < 0.00] <- 0.00



slopes_Poc <- ggplot(Fixed_Bayes_Poc_FE, aes(x = Estimate, y = RespVar,xmin = Q2.5, xmax = Q97.5, colour = "gray")) +
  geom_point() + geom_errorbar(width = 0.1) + 
  geom_vline(xintercept=-0.0, linetype="dotted", color = "red") +
  scale_colour_manual(name = 'Stats', values = setNames(c('red'),c("Across 0.00 = No Sig & Pos or Neg = Sign"))) +
  labs(x = "β slopes with 95% IC", y = "", title = "") + theme_classic() + theme (legend.position = "bottom")
slopes_Poc
ggsave ("Outputs_R/Figures/Bayesian_slopes_Poc.pdf", slopes_Poc, width = 6, height = 8)
# This is figure 3A right. 
# Significant variables are the ones where slope confidence intervals do NOT go across 0


