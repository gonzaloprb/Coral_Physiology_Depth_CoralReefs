# Script for analysing the isotopes with the already prepared Pachyseris speciosa database
# This is a simplified version going straight to the analysis and figures of the Manuscript, plus some SIBER plots and statistics

# The tests used here are:      SIBER: https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html
#                               Hotelling tests: https://cran.r-project.org/web/packages/Hotelling/Hotelling.pdf
#                               Relative Degree of Heterotrophy: Williams et al 2018, which is the study of the variable "Delta13C_Pol_Zoox" with depth. Here, there is only a plot, but additional study can be found in Script: "Pocillopora_All_Multivar_Bayesian.R" and "Pachyseris_All_Multivar_Bayesian.R" 

# AAA Read Me: "Zoox" and "Polyps" are later replaced with "Symbionts" and "Host" in the Manuscript.


# Clear database
rm (list = ls ())

# Necessary packages are:
library (SIBER); library (bayesm); library (dplyr); library (stringr); library (plyr); library (ggplot2); library (tidyr); library (car) #library (siar); library (mvnormtest);
library (Hotelling) ; library (corpcor)


# Open the data
Isotopes_Pachyseris <- read.csv ("Data_Codes/Isotopes/Pachyseris_Isotopes.csv", header = T, dec = ".", sep = ",")
resume_Pachyseris_Isotopes <- read.csv ("Data_Codes/Isotopes/Pachyseris_Isotopes_resume.csv", header = T, dec = ".", sep = ",")



### 1st  - Make the plot of δ15N vs δ13C 

# Remove all extra depths rather than (20, 40, 60 and 90) for visual simplification
Depths_Pachyseris <- c(20,40,60,90)
resume_Pachyseris_Isotopes <- subset(resume_Pachyseris_Isotopes, Depth %in% Depths_Pachyseris)

resume_Pachyseris_Isotopes$Island_Site <- gsub('Bora_2', 'Bora Bora', resume_Pachyseris_Isotopes$Island_Site)
resume_Pachyseris_Isotopes$Island_Site <- gsub('Moorea_2', 'Moorea', resume_Pachyseris_Isotopes$Island_Site)
resume_Pachyseris_Isotopes$Island_Site <- gsub('Rangiroa_3', 'Rangiroa', resume_Pachyseris_Isotopes$Island_Site)
resume_Pachyseris_Isotopes$Island_Site <- gsub('Tikehau_2', 'Tikehau', resume_Pachyseris_Isotopes$Island_Site)


resume_Pachyseris_Isotopes = resume_Pachyseris_Isotopes %>% mutate(Island_Site = factor(Island_Site, levels = c("Moorea","Bora Bora","Tikehau","Rangiroa")))


cols <- c("6" = "#BDD7E7", "20" = "#6BAED6", "40" = "#3182BD", "60" = "#08519C", "90" = "darkblue")
forms <- c("Polyps"=0,"Zoox"=1)

# 1st bis - or straight from the "resume_Pachyseris_Isotopes
Isot_P.spec <- ggplot(resume_Pachyseris_Isotopes, aes(x = Delta_Carbon13,  y = Delta_Nitrogen15, shape = Extraction, colour = as.factor (Depth))) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin=Delta_Nitrogen15-Delta_Nitrogen15_SD, ymax=Delta_Nitrogen15+Delta_Nitrogen15_SD), width=.2,position=position_dodge(0.05)) +
  geom_errorbar(aes(xmin=Delta_Carbon13-Delta_Carbon13_SD, xmax=Delta_Carbon13+Delta_Carbon13_SD), width=.2,position=position_dodge(0.05)) +
  facet_wrap(~Island_Site, ncol = 4) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = forms) +
  labs(x = "δ13C",y = "δ15N", title = "Pachyseris speciosa spp.", color = "Depth (m)") + 
  theme_bw () + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 
Isot_P.spec
ggsave ("Outputs_R/Figures/Isot_P.spec.pdf", Isot_P.spec, width = 12, height = 6)
# This is Fig 4b

# (Fig_Isotopes = Isot_P.verr + Isot_P.spec +  
#     plot_layout(guides = 'collect', ncol = 1)  & theme(legend.position='right'))
# ggsave( "Outputs_R/Figures/Fig_Isotopes_P.verr_P.spec.pdf", Fig_Isotopes,  width = 10, height = 6)



# Further analyses to check Polyp (host) - symbiont (zoox) per depth
ddply (Isotopes_Pachyseris, ~ Extraction, function (x){c(Max_δ13C = max(x$δ13C),Min_δ13C = min(x$δ13C),Max_δ15N = max(x$δ15N),Min_δ15N = min(x$δ15N))})
ddply (Isotopes_Pachyseris, ~ Depth + Extraction, function (x){c(Max_δ13C = max(x$δ13C),Min_δ13C = min(x$δ13C),Max_δ15N = max(x$δ15N),Min_δ15N = min(x$δ15N))})
ddply (Isotopes_Pachyseris, ~ Depth, function (x){c(Max_delta_δ13C = max(x$Delta13C_Pol_Zoox),Min_delta_δ13C = min(x$Delta13C_Pol_Zoox),Max_delta15N = max(x$Delta15N_Pol_Zoox),Min_delta15N = min(x$Delta15N_Pol_Zoox))})
ddply (Isotopes_Pachyseris, ~ Depth, function (x){c(Min_delta_δ13C = min(x$Delta13C_Pol_Zoox),Max_delta_δ13C = max(x$Delta13C_Pol_Zoox),Min_delta15N = min(x$Delta15N_Pol_Zoox),Max_delta15N = max(x$Delta15N_Pol_Zoox))})




# Check the normal distribution or not
hist (Isotopes_Pachyseris$δ15N)

# Check own independent normal distribution for each Island_Site and depth group (=depth)
# Change them manually: Site and Depth
df <- Isotopes_Pachyseris
for(i in 1:nrow(Isotopes_Pachyseris)) {
  if(Isotopes_Pachyseris$Island_Site[i] == "Tikehau_2" & Isotopes_Pachyseris$Depth[i] == "60") {
    df <- rbind(df, Isotopes_Pachyseris[i, ])
    # Create a histogram plot of one of the filtered columns
    hist(df$δ13C, main="")
  }
}



# SIBER
# First, for all islands together. Although it will be difficult to see anything.
# Create and give format to the database siber_data
siber_data <- Isotopes_Pachyseris %>%
  select(δ13C,δ15N,Depth, Extraction)

# I need to change the names 
# Depth = group
# Extraction = Community 
# iso 1 = δ13C
# iso 2 = δ15N
colnames (siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number). IMPORTANT! OTHERWISE code CRUSH!
siber_data <- siber_data[as.numeric(ave(siber_data$group, siber_data$group, FUN=length)) >= 5, ]

# Put in the right order
siber_data$group = factor(siber_data$group ,levels = c ("10","20","25","40","60","72","90"))
siber_data <- siber_data[order(siber_data$group),]

# Create the siber object 
siber.example <- createSiberObject(siber_data)

# Create lists of plotting arguments to be passed for future plotting functions.
# Plotting the raw data - calculate metrics
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")

palette(c("wheat", "green" ,"seagreen", "red","blue","navyblue","black"))

# Plot
par(mfrow=c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = F, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5)
# Add  ellipses by directly calling plot.group.ellipses() with p.interval % 
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  lty = 1, lwd = 2)
# Wheat = 10; # green = 20; # seagreen = 25; red = 40; blue = 60; navyblue = 72 and black = 90
# Circles = Polyps ; Triangles = Zoox

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)
# TA area of convex hull
# SEA standard ellipse area
# SEAc distance between the centroid mean (corrects also the differences between sample sizes)

# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(siber.example) 
print(community.ML)


## Run SIBER separately by Islands
#### Bora_2 ####
Bora_Isotopes_Pachyseris <- Isotopes_Pachyseris %>% filter(str_detect(Island_Site, "Bora_2"))

# Create and give format to the database
Bora_siber_data <- Bora_Isotopes_Pachyseris %>%
  select(δ13C,δ15N,Depth, Extraction)

# I need to change the names: Depth = group; Extraction = Community; iso 1 = δ13C; iso 2 = δ15N
colnames (Bora_siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number)
Bora_siber_data <- Bora_siber_data[as.numeric(ave(Bora_siber_data$group, Bora_siber_data$group, FUN=length)) >= 3, ]

# Put in the right order
Bora_siber_data$group = factor(Bora_siber_data$group ,levels = c ("40","60","90"))
Bora_siber_data <- Bora_siber_data[order(Bora_siber_data$group),]

# Create the siber object 
Bora_siber.example <- createSiberObject(Bora_siber_data)


## Plot with ggplot to show GroupHulls with lty according to community and colour to Depth group...

colggplot <- palette(c("red" ,"blue","black")) # Make sure this is run
shapeggplot <- c(19,1)

Bora_siber_data %>% ggplot(aes(iso1, iso2, fill = group, colour = group, shape = community, linetype=community))+
  geom_point()+ 
  stat_ellipse(type = "norm", geom = "polygon", level = 0.8, alpha = 0.1)+ 
  scale_fill_manual(values= colggplot) + scale_colour_manual (values= colggplot) + scale_shape_manual(values = shapeggplot) +
  xlab(expression({delta}^13*C~'\u2030'))+
  ylab(expression({delta}^15*N~'\u2030')) +
  theme_bw() + theme(plot.title = element_text(hjust=0.5, size=10, face="bold"),
                     axis.text.x = element_text(size=10, colour="black"),
                     axis.text.y = element_text(size=10, colour="black"),
                     axis.title.y = element_text(size=11, face="bold", colour="black"),
                     axis.title.x = element_text(size=11, face="bold", colour="black")) 

# Calculate summary statistics for each group: TA, SEA and SEAc
Bora_group.ML <- groupMetricsML(Bora_siber.example)
print(Bora_group.ML)

# Calculate the various Layman metrics on each of the communities.
Bora_community.ML <- communityMetricsML(Bora_siber.example) 
print(Bora_community.ML)

#### Bora_2 ####

#### Moorea_2 ####
Moorea_Isotopes_Pachyseris <- Isotopes_Pachyseris %>% filter(str_detect(Island_Site, "Moorea_2"))

# Create and give format to the database
Moorea_siber_data <- Moorea_Isotopes_Pachyseris %>%
  select(δ13C,δ15N,Depth, Extraction)

# I need to change the names: Depth = group; Extraction = Community; iso 1 = δ13C; iso 2 = δ15N
colnames (Moorea_siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number)
Moorea_siber_data <- Moorea_siber_data[as.numeric(ave(Moorea_siber_data$group, Moorea_siber_data$group, FUN=length)) >= 3, ]

# Put in the right order
Moorea_siber_data$group = factor(Moorea_siber_data$group ,levels = c ("40","60","72", "77", "90"))
Moorea_siber_data <- Moorea_siber_data[order(Moorea_siber_data$group),]

# Create the siber object 
Moorea_siber.example <- createSiberObject(Moorea_siber_data)



## Plot with ggplot 

colggplot <- palette(c( "red" ,"blue","black" )) # Filtering below the 72 and 77 depth groups
shapeggplot <- c(19,1)

Moorea_siber_data <- subset(Moorea_siber_data,group %in% c("40" , "60", "90")) # Too few points to draw ellipse
Moorea_siber_data %>% ggplot(aes(iso1, iso2, fill = group, colour = group, shape = community, linetype=community))+
  geom_point()+ 
  stat_ellipse(type = "norm", geom = "polygon", level = 0.8, alpha = 0.1)+ 
  scale_fill_manual(values= colggplot) + scale_colour_manual (values= colggplot) + scale_shape_manual(values = shapeggplot) +
  xlab(expression({delta}^13*C~'\u2030'))+
  ylab(expression({delta}^15*N~'\u2030')) +
  theme_bw() + theme(plot.title = element_text(hjust=0.5, size=10, face="bold"),
                     axis.text.x = element_text(size=10, colour="black"),
                     axis.text.y = element_text(size=10, colour="black"),
                     axis.title.y = element_text(size=11, face="bold", colour="black"),
                     axis.title.x = element_text(size=11, face="bold", colour="black")) 


# Calculate summary statistics for each group: TA, SEA and SEAc
Moorea_group.ML <- groupMetricsML(Moorea_siber.example)
print(Moorea_group.ML)

# Calculate the various Layman metrics on each of the communities.
Moorea_community.ML <- communityMetricsML(Moorea_siber.example) 
print(Moorea_community.ML)

#### Moorea_2 ####


#### Tikehau_2 ####
Tikehau_Isotopes_Pachyseris <- Isotopes_Pachyseris %>% filter(str_detect(Island_Site, "Tikehau_2"))

# Create and give format to the database
Tikehau_siber_data <- Tikehau_Isotopes_Pachyseris %>%
  select(δ13C,δ15N,Depth, Extraction)

# I need to change the names: Depth = group; Extraction = Community; iso 1 = δ13C; iso 2 = δ15N
colnames (Tikehau_siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number)
Tikehau_siber_data <- Tikehau_siber_data[as.numeric(ave(Tikehau_siber_data$group, Tikehau_siber_data$group, FUN=length)) >= 3, ]

# Put in the right order
Tikehau_siber_data$group = factor(Tikehau_siber_data$group ,levels = c ("10","20","40","60","90"))
Tikehau_siber_data <- Tikehau_siber_data[order(Tikehau_siber_data$group),]

# Create the siber object 
Tikehau_siber.example <- createSiberObject(Tikehau_siber_data)




## Plot with ggplot 

colggplot <- palette(c("wheat","green" ,"red","blue", "black")) # Make sure this is run!
shapeggplot <- c(19,1)

Tikehau_siber_data %>% ggplot(aes(iso1, iso2, fill = group, colour = group, shape = community, linetype=community))+
  geom_point()+ 
  stat_ellipse(type = "norm", geom = "polygon", level = 0.8, alpha = 0.1)+ 
  scale_fill_manual(values= colggplot) + scale_colour_manual (values= colggplot) + scale_shape_manual(values = shapeggplot) +
  xlab(expression({delta}^13*C~'\u2030'))+
  ylab(expression({delta}^15*N~'\u2030')) +
  theme_bw() + theme(plot.title = element_text(hjust=0.5, size=10, face="bold"),
                     axis.text.x = element_text(size=10, colour="black"),
                     axis.text.y = element_text(size=10, colour="black"),
                     axis.title.y = element_text(size=11, face="bold", colour="black"),
                     axis.title.x = element_text(size=11, face="bold", colour="black")) 


# Calculate summary statistics for each group: TA, SEA and SEAc
Tikehau_group.ML <- groupMetricsML(Tikehau_siber.example)
print(Tikehau_group.ML)

# Calculate the various Layman metrics on each of the communities.
Tikehau_community.ML <- communityMetricsML(Tikehau_siber.example) 
print(Tikehau_community.ML)

#### Tikehau_2 ####

#### Rangiroa_3 ####
Rangiroa_Isotopes_Pachyseris <- Isotopes_Pachyseris %>% filter(str_detect(Island_Site, "Rangiroa_3"))

# Create and give format to the database
Rangiroa_siber_data <- Rangiroa_Isotopes_Pachyseris %>%
  select(δ13C,δ15N,Depth, Extraction)

# I need to change the names: Depth = group; Extraction = Community; iso 1 = δ13C; iso 2 = δ15N
colnames (Rangiroa_siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number)
Rangiroa_siber_data <- Rangiroa_siber_data[as.numeric(ave(Rangiroa_siber_data$group, Rangiroa_siber_data$group, FUN=length)) >= 3, ]

# Put in the right order
Rangiroa_siber_data$group = factor(Rangiroa_siber_data$group ,levels = c ("25","40","60","90"))
Rangiroa_siber_data <- Rangiroa_siber_data[order(Rangiroa_siber_data$group),]

# Create the siber object 
Rangiroa_siber.example <- createSiberObject(Rangiroa_siber_data)

## Plot with ggplot 

colggplot <- palette(c( "green" ,"red","dodgerblue", "black")) # Make sure this is run!
shapeggplot <- c(19,1)

Rangiroa_siber_data %>% ggplot(aes(iso1, iso2, fill = group, colour = group, shape = community, linetype=community))+
  geom_point()+ 
  stat_ellipse(type = "norm", geom = "polygon", level = 0.8, alpha = 0.1)+ 
  scale_fill_manual(values= colggplot) + scale_colour_manual (values= colggplot) + scale_shape_manual(values = shapeggplot) +
  xlab(expression({delta}^13*C~'\u2030'))+
  ylab(expression({delta}^15*N~'\u2030')) +
  theme_bw() + theme(plot.title = element_text(hjust=0.5, size=10, face="bold"),
                     axis.text.x = element_text(size=10, colour="black"),
                     axis.text.y = element_text(size=10, colour="black"),
                     axis.title.y = element_text(size=11, face="bold", colour="black"),
                     axis.title.x = element_text(size=11, face="bold", colour="black")) 


# Calculate summary statistics for each group: TA, SEA and SEAc
Rangiroa_group.ML <- groupMetricsML(Rangiroa_siber.example)
print(Rangiroa_group.ML)


# Calculate the various Layman metrics on each of the communities.
Rangiroa_community.ML <- communityMetricsML(Rangiroa_siber.example) 
print(Rangiroa_community.ML)
#### Rangiroa_3 ####




### Groups combining Depth + Island or (Island + Depth) in the same plot; it will be difficult to visualise
Isotopes_Pachyseris <- Isotopes_Pachyseris %>% unite(Depth_Island, Depth,Island_Site, remove = F)

# Create and give format to the database - keeping all combinations of Island + Depth
Combined_siber_data <- Isotopes_Pachyseris %>%
  select(δ13C,δ15N,Depth_Island, Extraction)

# I need to change the names: Island_Depth = group; Extraction = Community; iso 1 = δ13C; iso 2 = δ15N
colnames (Combined_siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number); here I set 5 to avoid error and warnings in the plots
Combined_siber_data <- Combined_siber_data[as.numeric(ave(Combined_siber_data$group, Combined_siber_data$group, FUN=length)) >= 5, ]

# Put in the right order, I am actually not sure if this is necessary or correct
Combined_siber_data$group = factor(Combined_siber_data$group ,levels = c ("10_Tikehau_2","20_Tikehau_2","25_Rangiroa_3","40_Bora_2",
                                                                          "40_Moorea_2","40_Tikehau_2", "40_Rangiroa_3", "60_Bora_2", 
                                                                          "60_Moorea_2","60_Tikehau_2" ,"60_Rangiroa_3","72_Moorea_2",
                                                                          "90_Bora_2","90_Tikehau_2","90_Rangiroa_3"))

Combined_siber_data <- Combined_siber_data[order(Combined_siber_data$group),]

# Create the siber object 
Combined_siber.example <- createSiberObject(Combined_siber_data)


## Plot with ggplot 

colggplot <- palette(c("green","green1","green2","green3", "red", "red1", "red2", "red3","dodgerblue","dodgerblue1","dodgerblue2","dodgerblue3","mediumblue","navy","black")) # Make sure this is run!
shapeggplot <- c(19,1)

Combined_siber_data %>% ggplot(aes(iso1, iso2, fill = group, colour = group, shape = community, linetype=community))+
  geom_point()+ 
  stat_ellipse(type = "norm", geom = "polygon", level = 0.8, alpha = 0.1)+ 
  scale_fill_manual(values= colggplot) + scale_colour_manual (values= colggplot) + scale_shape_manual(values = shapeggplot) +
  xlab(expression({delta}^13*C~'\u2030'))+
  ylab(expression({delta}^15*N~'\u2030')) +
  theme_bw() + theme(plot.title = element_text(hjust=0.5, size=10, face="bold"),
                     axis.text.x = element_text(size=10, colour="black"),
                     axis.text.y = element_text(size=10, colour="black"),
                     axis.title.y = element_text(size=11, face="bold", colour="black"),
                     axis.title.x = element_text(size=11, face="bold", colour="black")) 



# Calculate summary statistics for each group: TA, SEA and SEAc
Combined_group.ML <- groupMetricsML(Combined_siber.example)
print(Combined_group.ML) 

# Calculate the various Layman metrics on each of the communities.
Combined_community.ML <- communityMetricsML(Combined_siber.example) 
print(Combined_community.ML)




### Plot the Standard Ellipse Area (%) with depth like in Radice et al 2019 - either individually or straight for all together "Combined_group.ML"
Combined_group.ML <- as.data.frame(print(Combined_group.ML))
# Keep only SEAc
Matrix_SEAc <- Combined_group.ML[c("SEAc"),]

# Colnames to rownames
Matrix_SEAc <- data.frame(y=unlist(Matrix_SEAc))
colnames (Matrix_SEAc) <- "SEAc"
Matrix_SEAc$Island = sapply(strsplit(rownames(Matrix_SEAc), "_"), function(x) x[2])
Matrix_SEAc$Fraction_Depth = sapply(strsplit(rownames(Matrix_SEAc), "_"), function(x) x[1])

Matrix_SEAc$Fraction_Depth <- gsub(".", "_", Matrix_SEAc$Fraction_Depth , fixed = TRUE)

Matrix_SEAc$Depth = sapply(strsplit(Matrix_SEAc$Fraction_Depth, "_"), function(x) x[2])
Matrix_SEAc$Fraction = sapply(strsplit(Matrix_SEAc$Fraction_Depth, "_"), function(x) x[1])

Matrix_SEAc = Matrix_SEAc %>% mutate(Depth = factor(Depth, levels = c("10","20", "25","40", "60","72", "90")))

ggplot(Matrix_SEAc, aes(x=factor(Depth), y=SEAc)) + geom_point (aes(x=factor(Depth), y=SEAc, color = Island), size = 1) +
  geom_boxplot() + facet_grid(cols = vars(Fraction), rows = vars (Island)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=3) + theme_bw()  + 
  ylab ("Standard Ellipse Area (‰)") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                                  axis.text = element_text(size=10, colour="black"),
                                                                  axis.title = element_text(size=11, face="bold", colour="black")) 

# It does not seem to follow any clear pattern with depth
# Extra statistics across all islands or individually are available upon request. 


### Hotelling tests 

# It measures the differences/distances between fraction of Polyps and Zoox (Symbionts)! 
# Hotelling test is a multivariate analogue of the univariate t test that is suited for comparison of population mean vectors!
# More distance more heterotrophic
# Less distance more autotrophic
# P-value > 0.05 there IS NOT difference between the the two multivariate means
# P-value < 0.05 there IS difference between the the two multivariate means


# Here I run it for the whole community together (sites and by depths). However, it is Necessary to subset by sites + depths
# Considering the whole community (groups/depths) all together to see fraction-community
print(hotelling.test(.~community , data = siber_data [,c("iso1","iso2","community")]))

# This is considering all sites and depths together...


# Extract the separate depths make test for 10, 20...90 / but keeping all sites together
Ten_m_siber_data <- subset (siber_data, group == 10)
print(hotelling.test(.~community, data = Ten_m_siber_data [,c(1,2,4)]))


Twenty_m_siber_data <- subset (siber_data, group == 20)
print(hotelling.test(.~community, data = Twenty_m_siber_data [,c(1,2,4)]))


Twentyfive_m_siber_data <- subset (siber_data, group == 25)
print(hotelling.test(.~community, data = Twentyfive_m_siber_data [,c(1,2,4)]))


Forty_m_siber_data <- subset (siber_data, group == 40)
print(hotelling.test(.~community, data = Forty_m_siber_data [,c(1,2,4)]))


Sixty_m_siber_data <- subset (siber_data, group == 60)
print(hotelling.test(.~community, data = Sixty_m_siber_data [,c(1,2,4)]))


Seventytwo_m_siber_data <- subset (siber_data, group == 72)
print(hotelling.test(.~community, data = Seventytwo_m_siber_data [,c(1,2,4)]))


Ninety_m_siber_data <- subset (siber_data, group == 90)
print(hotelling.test(.~community, data = Ninety_m_siber_data [,c(1,2,4)]))



# Considering separated by islands, depths and community 
# This is what is answering our questions!
# We are also adding the perm = T option to increase statistical power!


#### Bora_2 ####
# Extract 40
Forty_m_Bora_siber_data <- subset (Bora_siber_data, group == 40)
# with permutations
fit_Bora_fortym_perm <- hotelling.test(iso1+iso2 ~community, data  = Forty_m_Bora_siber_data, perm =  TRUE)
fit_Bora_fortym_perm

Sixty_m_Bora_siber_data <- subset (Bora_siber_data, group == 60)
fit_Bora_sixtym_perm <- hotelling.test(iso1+iso2 ~community, data  = Sixty_m_Bora_siber_data, perm =  TRUE)
fit_Bora_sixtym_perm

Ninety_m_Bora_siber_data <- subset (Bora_siber_data, group == 90)
fit_Bora_ninetym_perm <- hotelling.test(iso1+iso2 ~community, data  = Ninety_m_Bora_siber_data, perm =  TRUE)
fit_Bora_ninetym_perm
#### Bora_2 ####



#### Moorea_2 ####
Forty_m_Moorea_siber_data <- subset (Moorea_siber_data, group == 40)
fit_Moorea_fortym_perm <- hotelling.test(iso1+iso2 ~community, data  = Forty_m_Moorea_siber_data, perm =  TRUE)
fit_Moorea_fortym_perm

Sixty_m_Moorea_siber_data <- subset (Moorea_siber_data, group == 60)
fit_Moorea_sixtym_perm <- hotelling.test(iso1+iso2 ~community, data  = Sixty_m_Moorea_siber_data, perm =  TRUE)
fit_Moorea_sixtym_perm

Seventytwo_m_Moorea_siber_data <- subset (Moorea_siber_data, group == 72)
fit_Moorea_seventytwom_perm <- hotelling.test(iso1+iso2 ~community, data  = Seventytwo_m_Moorea_siber_data, perm =  TRUE)
fit_Moorea_seventytwom_perm
#### Moorea_2 ####


#### Tikehau_2 ####
Ten_m_Tikehau_siber_data <- subset (Tikehau_siber_data, group == 10)
fit_Tikehau_tenm_perm <- hotelling.test(iso1+iso2 ~community, data  = Ten_m_Tikehau_siber_data, perm =  TRUE)
fit_Tikehau_tenm_perm

# Extract 20 
Twenty_m_Tikehau_siber_data <- subset (Tikehau_siber_data, group == 20)
fit_Tikehau_twentym_perm <- hotelling.test(iso1+iso2 ~community, data  = Twenty_m_Tikehau_siber_data, perm =  TRUE)
fit_Tikehau_twentym_perm

Forty_m_Tikehau_siber_data <- subset (Tikehau_siber_data, group == 40)
fit_Tikehau_fortym_perm <- hotelling.test(iso1+iso2 ~community, data  = Forty_m_Tikehau_siber_data, perm =  TRUE)
fit_Tikehau_fortym_perm

Sixty_m_Tikehau_siber_data <- subset (Tikehau_siber_data, group == 60)
fit_Tikehau_sixtym_perm <- hotelling.test(iso1+iso2 ~community, data  = Sixty_m_Tikehau_siber_data, perm =  TRUE)
fit_Tikehau_sixtym_perm

Ninety_m_Tikehau_siber_data <- subset (Tikehau_siber_data, group == 90)
fit_Tikehau_ninetym_perm <- hotelling.test(iso1+iso2 ~community, data  = Ninety_m_Tikehau_siber_data, perm =  TRUE)
fit_Tikehau_ninetym_perm
#### Tikehau_2 ####

#### Rangiroa_3 ####
# Extract 25
Twentyfive_m_Rangiroa_siber_data <- subset (Rangiroa_siber_data, group == 25)
fit_Rangiroa_twentyfivem_perm <- hotelling.test(iso1+iso2 ~community, data  = Twentyfive_m_Rangiroa_siber_data, perm =  TRUE)
fit_Rangiroa_twentyfivem_perm


Forty_m_Rangiroa_siber_data <- subset (Rangiroa_siber_data, group == 40)
fit_Rangiroa_fortym_perm <- hotelling.test(iso1+iso2 ~community, data  = Forty_m_Rangiroa_siber_data, perm =  TRUE)
fit_Rangiroa_fortym_perm

Sixty_m_Rangiroa_siber_data <- subset (Rangiroa_siber_data, group == 60)
fit_Rangiroa_sixtym_perm <- hotelling.test(iso1+iso2 ~community, data  = Sixty_m_Rangiroa_siber_data, perm =  TRUE)
fit_Rangiroa_sixtym_perm

Ninety_m_Rangiroa_siber_data <- subset (Rangiroa_siber_data, group == 90)
fit_Rangiroa_ninetym_perm <- hotelling.test(iso1+iso2 ~community, data  = Ninety_m_Rangiroa_siber_data, perm =  TRUE)
fit_Rangiroa_ninetym_perm
#### Rangiroa_3 ####

### End of Hotelling tests


### Finally, Relative Degree of Heterotrophy (RDH) Williams et al 2018, quick visualization. 
ggplot(Isotopes_Pachyseris, aes(x=factor(Depth), y=Delta13C_Pol_Zoox)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Delta δ13C Polyps - Zoox") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                                 axis.text = element_text(size=10, colour="black"),
                                                                 axis.title = element_text(size=11, face="bold", colour="black")) 

# RDH = "Delta13C_Pol_Zoox" with depth. Additional tests on this RDH can be found in Script: "Pachyseris_All_Multivar_Bayesian.R" 



# Many other tests have been done, but for ease use of the data and interpretation, these were excluded from the manuscript and I have not included them here. 
# These extra tests are: 
# Additional tests from SIBER package: https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html and Jackson et al 2011 and Conti-Jerpe et al 2020
# Additional tests as suggested by Radice et al 2019
# Bayesian modelling of the isotope variables with depth. Also available in Script: "Pachyseris_All_Multivar_Bayesian.R" 

# The whole script with all different tests can be shared upon request. 





