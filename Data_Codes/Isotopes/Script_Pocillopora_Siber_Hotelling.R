# Script for analysing the isotopes with the already prepared Pocillopora verrucosa database
# This is a simplified version going straight to the analysis and figures of the Manuscript, plus some SIBER plots and statistics

# The tests used here are:      SIBER: https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html
#                               Hotelling tests: https://cran.r-project.org/web/packages/Hotelling/Hotelling.pdf
#                               Relative Degree of Heterotrophy: Williams et al 2018, which is the study of the variable "Delta13C_Pol_Zoox" with depth. Here, there is only a plot, but additional study can be found in Script: "Pocillopora_All_Multivar_Bayesian.R" and "Pachyseris_All_Multivar_Bayesian.R" 

# AAA Read Me: "Zoox" and "Polyps" are later replaced with "Symbionts" and "Host" in the Manuscript.

# Clear database
rm (list = ls ())

# Necessary packages are:
library (SIBER); library (bayesm); library (dplyr); library (stringr); library (plyr); library (ggplot2); library (tidyr) #library (siar); library (mvnormtest);
library (Hotelling) ; library (corpcor)


# Open the data:
Isotopes_Pocillopora <- read.csv ("Data_Codes/Isotopes/Pocillopora_Isotopes.csv", header = T, dec = ".", sep = ",")
resume_Pocillopora_Isotopes <- read.csv ("Data_Codes/Isotopes/Pocillopora_Isotopes_resume.csv", header = T, dec = ".", sep = ",")


### 1st - Make the plot of δ15N vs δ13C 
resume_Pocillopora_Isotopes$Island_Site <- gsub('Bora_2', 'Bora Bora', resume_Pocillopora_Isotopes$Island_Site)
resume_Pocillopora_Isotopes$Island_Site <- gsub('Moorea_2', 'Moorea', resume_Pocillopora_Isotopes$Island_Site)
resume_Pocillopora_Isotopes$Island_Site <- gsub('Rangiroa_3', 'Rangiroa', resume_Pocillopora_Isotopes$Island_Site)
resume_Pocillopora_Isotopes$Island_Site <- gsub('Tikehau_2', 'Tikehau', resume_Pocillopora_Isotopes$Island_Site)


resume_Pocillopora_Isotopes = resume_Pocillopora_Isotopes %>% mutate(Island_Site = factor(Island_Site, levels = c("Moorea","Bora Bora","Tikehau","Rangiroa")))

cols <- c("6" = "#BDD7E7", "20" = "#6BAED6", "40" = "#3182BD", "60" = "#08519C", "90" = "darkblue")
forms <- c("Polyps"=0,"Zoox"=1)

Isot_P.verr <- ggplot(resume_Pocillopora_Isotopes, aes(x = Delta_Carbon13,  y = Delta_Nitrogen15, shape = Extraction, colour = as.factor (Depth))) +
  geom_point(size = 3.5) +  
  geom_errorbar(aes(ymin=Delta_Nitrogen15-Delta_Nitrogen15_SD, ymax=Delta_Nitrogen15+Delta_Nitrogen15_SD), width=.2,position=position_dodge(0.05)) +
  geom_errorbar(aes(xmin=Delta_Carbon13-Delta_Carbon13_SD, xmax=Delta_Carbon13+Delta_Carbon13_SD), width=.2,position=position_dodge(0.05)) +
  facet_wrap(~Island_Site, ncol = 4) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = forms) +
  labs(x = "δ13C",y = "δ15N", title = "Pocillopora cf. verrucosa", color = "Depth (m)") + 
  theme_bw () + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black")) 
Isot_P.verr
ggsave ("Outputs_R/Figures/Isot_P.verr.pdf", Isot_P.verr, width = 12, height = 6)
# This is Fig. 4a


# Further analysis. Check Max and Min values Polyp (host) - symbiont (zoox) per depth
ddply (Isotopes_Pocillopora, ~ Extraction, function (x){c(Max_δ13C = max(x$δ13C),Min_δ13C = min(x$δ13C),Max_δ15N = max(x$δ15N),Min_δ15N = min(x$δ15N))})
ddply (Isotopes_Pocillopora, ~ Depth + Extraction, function (x){c(Max_δ13C = max(x$δ13C),Min_δ13C = min(x$δ13C),Max_δ15N = max(x$δ15N),Min_δ15N = min(x$δ15N))})
ddply (Isotopes_Pocillopora, ~ Depth, function (x){c(Max_delta_δ13C = max(x$Delta13C_Pol_Zoox),Min_delta_δ13C = min(x$Delta13C_Pol_Zoox),Max_delta15N = max(x$Delta15N_Pol_Zoox),Min_delta15N = min(x$Delta15N_Pol_Zoox))})
ddply (Isotopes_Pocillopora, ~ Depth, function (x){c(Min_delta_δ13C = min(x$Delta13C_Pol_Zoox),Max_delta_δ13C = max(x$Delta13C_Pol_Zoox),Min_delta15N = min(x$Delta15N_Pol_Zoox),Max_delta15N = max(x$Delta15N_Pol_Zoox))})


# Check if normal distribution 
hist (Isotopes_Pocillopora$δ15N)

# Check own independent normal distribution for each Site and depth 
# Change them manually: Site and Depth
df <- Isotopes_Pocillopora
for(i in 1:nrow(Isotopes_Pocillopora)) {
  if(Isotopes_Pocillopora$Island_Site[i] == "Moorea_2" & Isotopes_Pocillopora$Depth[i] == "20") {
    df <- rbind(df, Isotopes_Pocillopora[i, ])
    # Create a histogram plot of one of the filtered columns
    hist(df$δ15N, main="")
  }
}



# SIBER
# First, for all islands together. Although it will be difficult to see anything.
# Create and give format to the database siber_data
siber_data <- Isotopes_Pocillopora %>%
  select(δ13C,δ15N,Depth, Extraction)

# I need to change the names 
# Depth = group
# Extraction = Community 
# iso 1 = δ13C
# iso 2 = δ15N
colnames (siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations/replicates (absolute minimum number). IMPORTANT! OTHERWISE code CRUSH!
siber_data <- siber_data[as.numeric(ave(siber_data$group, siber_data$group, FUN=length)) >= 3, ]

# Put in the right order
siber_data$group = factor(siber_data$group ,levels = c ("6","20","40","60"))
siber_data <- siber_data[order(siber_data$group),]


# Create the Siber object 
siber.example <- createSiberObject(siber_data)

# Create lists of plotting arguments to be passed for future plotting functions.
# Plotting the raw data - calculate metrics
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")
palette(c("green" ,"red","blue","black"))


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

# Green = 6m ; # red = 20m ; # Blue = 40m ; # black = 60m 
# Circles = Polyps ; Triangles = Zoox

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)
# TA area of convex hull
# SEA standard ellipse area
# SEAc distance between the centroid mean (corrects also the differences between sample sizes)





## Run SIBER separately by Islands

#### Bora_2 ####
Bora_Isotopes_Pocillopora <- Isotopes_Pocillopora %>% filter(str_detect(Island_Site, "Bora_2"))

# Create and give format to the database
Bora_siber_data <- Bora_Isotopes_Pocillopora %>%
  select(δ13C,δ15N,Depth, Extraction)

# I need to change the names: Depth = group; Extraction = Community; iso 1 = δ13C; iso 2 = δ15N
colnames (Bora_siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number)
Bora_siber_data <- Bora_siber_data[as.numeric(ave(Bora_siber_data$group, Bora_siber_data$group, FUN=length)) >= 3, ]

# Put in the right order
Bora_siber_data$group = factor(Bora_siber_data$group ,levels = c ("6","20","40","60"))
Bora_siber_data <- Bora_siber_data[order(Bora_siber_data$group),]

# Create the siber object 
Bora_siber.example <- createSiberObject(Bora_siber_data)


## Plot with ggplot to show GroupHulls with lty according to community and colour to Depth group

colggplot <- palette(c( "green" ,"red","dodgerblue","black"))
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
Moorea_Isotopes_Pocillopora <- Isotopes_Pocillopora %>% filter(str_detect(Island_Site, "Moorea_2"))

# Create and give format to the database
Moorea_siber_data <- Moorea_Isotopes_Pocillopora %>%
  select(δ13C,δ15N,Depth, Extraction)

# I need to change the names: Depth = group; Extraction = Community; iso 1 = δ13C; iso 2 = δ15N
colnames (Moorea_siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number)
Moorea_siber_data <- Moorea_siber_data[as.numeric(ave(Moorea_siber_data$group, Moorea_siber_data$group, FUN=length)) >= 3, ]

# Put in the right order
Moorea_siber_data$group = factor(Moorea_siber_data$group ,levels = c ("6","20","40","60"))
Moorea_siber_data <- Moorea_siber_data[order(Moorea_siber_data$group),]

# Create the siber object 
Moorea_siber.example <- createSiberObject(Moorea_siber_data)


## Plot with ggplot to show GroupHulls with lty according to community and colour to Depth group

colggplot <- palette(c( "green" ,"red","dodgerblue","black"))
shapeggplot <- c(19,1)

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
print(Moorea_group.ML) # here it increases, in Bora Bora decreases

# Calculate the various Layman metrics on each of the communities.
Moorea_community.ML <- communityMetricsML(Moorea_siber.example) 
print(Moorea_community.ML)
#### Moorea_2 ####


#### Tikehau_2 ####
Tikehau_Isotopes_Pocillopora <- Isotopes_Pocillopora %>% filter(str_detect(Island_Site, "Tikehau_2"))

# Create and give format to the database
Tikehau_siber_data <- Tikehau_Isotopes_Pocillopora %>%
  select(δ13C,δ15N,Depth, Extraction)

# I need to change the names: Depth = group; Extraction = Community; iso 1 = δ13C; iso 2 = δ15N
colnames (Tikehau_siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number)
Tikehau_siber_data <- Tikehau_siber_data[as.numeric(ave(Tikehau_siber_data$group, Tikehau_siber_data$group, FUN=length)) >= 3, ]

# Put in the right order
Tikehau_siber_data$group = factor(Tikehau_siber_data$group ,levels = c ("6","20","40"))
Tikehau_siber_data <- Tikehau_siber_data[order(Tikehau_siber_data$group),]

# Create the siber object 
Tikehau_siber.example <- createSiberObject(Tikehau_siber_data)




## Plot with ggplot 

colggplot <- palette(c( "green" ,"red","dodgerblue")) # Make sure this is run!
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
print(Tikehau_group.ML) # Here it decreases but not for Zoox

# Calculate the various Layman metrics on each of the communities.
Tikehau_community.ML <- communityMetricsML(Tikehau_siber.example) 
print(Tikehau_community.ML)
#### Tikehau_2 ####


#### Rangiroa_3 ####
Rangiroa_Isotopes_Pocillopora <- Isotopes_Pocillopora %>% filter(str_detect(Island_Site, "Rangiroa_3"))

# Create and give format to the database
Rangiroa_siber_data <- Rangiroa_Isotopes_Pocillopora %>%
  select(δ13C,δ15N,Depth, Extraction)

# I need to change the names: Depth = group; Extraction = Community; iso 1 = δ13C; iso 2 = δ15N
colnames (Rangiroa_siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number)
Rangiroa_siber_data <- Rangiroa_siber_data[as.numeric(ave(Rangiroa_siber_data$group, Rangiroa_siber_data$group, FUN=length)) >= 3, ]

# Put in the right order
Rangiroa_siber_data$group = factor(Rangiroa_siber_data$group ,levels = c ("6","20","40"))
Rangiroa_siber_data <- Rangiroa_siber_data[order(Rangiroa_siber_data$group),]

# Create the siber object 
Rangiroa_siber.example <- createSiberObject(Rangiroa_siber_data)



## Plot with ggplot 

colggplot <- palette(c( "green" ,"red","dodgerblue"))
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
print(Rangiroa_group.ML) # Non-sense humped relationship with depth

# Calculate the various Layman metrics on each of the communities.
Rangiroa_community.ML <- communityMetricsML(Rangiroa_siber.example) 
print(Rangiroa_community.ML)
#### Rangiroa_3 ####



### Groups combining Depth + Island or (Island + Depth) in the same plot
Isotopes_Pocillopora <- Isotopes_Pocillopora %>% unite(Depth_Island, Depth,Island_Site, remove = F)

# Create and give format to the database - keeping all combinations of Island + Depth
Combined_siber_data <- Isotopes_Pocillopora %>%
  select(δ13C,δ15N,Depth_Island, Extraction)

# I need to change the names: Island_Depth = group; Extraction = Community; iso 1 = δ13C; iso 2 = δ15N
colnames (Combined_siber_data) <- c("iso1", "iso2", "group", "community")

# Necessary to delete all groups that have less than 3 observations (absolute minimum number); here I set 5 to avoid error and warnings in the plots
Combined_siber_data <- Combined_siber_data[as.numeric(ave(Combined_siber_data$group, Combined_siber_data$group, FUN=length)) >= 5, ]

# Put in the right order, I am actually not sure if this is necessary or correct
Combined_siber_data$group = factor(Combined_siber_data$group ,levels = c ("6_Bora_2","6_Moorea_2","6_Tikehau_2","6_Rangiroa_3",
                                                                                        "20_Bora_2","20_Moorea_2", "20_Tikehau_2", "20_Rangiroa_3", 
                                                                                        "40_Bora_2","40_Moorea_2" ,"40_Tikehau_2","40_Rangiroa_3",
                                                                                        "60_Bora_2","60_Moorea_2"))

Combined_siber_data <- Combined_siber_data[order(Combined_siber_data$group),]

# Create the siber object 
Combined_siber.example <- createSiberObject(Combined_siber_data)



## Plot with ggplot 

colggplot <- palette(c("green","green1","green2","green3", "red", "red1", "red2", "red3","dodgerblue","dodgerblue1","dodgerblue2","dodgerblue3","navy","black")) # Make sure this is well run!
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

Matrix_SEAc = Matrix_SEAc %>% mutate(Depth = factor(Depth, levels = c("6","20", "40", "60")))


ggplot(Matrix_SEAc, aes(x=factor(Depth), y=SEAc)) + geom_point (aes(x=factor(Depth), y=SEAc, color = Island), size = 1) +
  geom_boxplot() + facet_grid(cols = vars(Fraction), rows = vars (Island)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=3) + theme_bw()  + 
  ylab ("Standard Ellipse Area (‰)") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                                  axis.text = element_text(size=10, colour="black"),
                                                                  axis.title = element_text(size=11, face="bold", colour="black")) 
# You can see that in some cases it seems to decrease with depth but not always. Extra statistics across all islands or individually are available upon request. 




### Hotelling tests 

# It measures the differences/distances between fraction of Polyps and Zoox (Symbionts)! 
# Hotelling test is a multivariate analogue of the univariate t test that is suited for comparison of population mean vectors!
# More distance more heterotrophic
# Less distance more autotrophic
# P-value > 0.05 there IS NOT difference between the the two multivariate means
# P-value < 0.05 there IS difference between the the two multivariate means


# Here I run it for the whole community together (sites and by depths). However, it is Necessary to subset by sites + depths
# Considering the whole community (groups/depths) all together to see fraction-community
print (hotelling.test(.~community , data = siber_data [,c("iso1","iso2","community")])) 

# This is considering all sites and depths together. It does not answer our question. 

# Separating by depths but keeping all sites together.
Six_m_siber_data <- subset (siber_data, group == 6)
print (hotelling.test(.~community, data = Six_m_siber_data [,c(1,2,4)]))

Twenty_m_siber_data <- subset (siber_data, group == 20)
print (hotelling.test(.~community, data = Twenty_m_siber_data [,c(1,2,4)]))

Forty_m_siber_data <- subset (siber_data, group == 40)
print (hotelling.test(.~community, data = Forty_m_siber_data [,c(1,2,4)]))



# What answers our question is how it is changing across islands and according to depth
# We need to consider separate islands, depth and community 
# We are also adding the perm = T option to increase statistical power!

#### Bora_2 ####
# Extract 6 
Six_m_Bora_siber_data <- subset (Bora_siber_data, group == 6)
# with permutations
fit_Bora_sixm_perm <- hotelling.test(iso1+iso2 ~community, data  = Six_m_Bora_siber_data, perm =  TRUE)
fit_Bora_sixm_perm

# Extract 20 
Twenty_m_Bora_siber_data <- subset (Bora_siber_data, group == 20)
# with permutations
fit_Bora_twentym_perm <- hotelling.test(iso1+iso2 ~community, data  = Twenty_m_Bora_siber_data, perm =  TRUE)
fit_Bora_twentym_perm

Forty_m_Bora_siber_data <- subset (Bora_siber_data, group == 40)
# with permutations
fit_Bora_fortym_perm <- hotelling.test(iso1+iso2 ~community, data  = Forty_m_Bora_siber_data, perm =  TRUE)
fit_Bora_fortym_perm

Sixty_m_Bora_siber_data <- subset (Bora_siber_data, group == 60)
# with permutations
fit_Bora_sixtym_perm <- hotelling.test(iso1+iso2 ~community, data  = Sixty_m_Bora_siber_data, perm =  TRUE)
fit_Bora_sixtym_perm
#### Bora_2 ####


#### Moorea_2 ####
# Extract 6 
Six_m_Moorea_siber_data <- subset (Moorea_siber_data, group == 6)
# with permutations
fit_Moorea_sixm_perm <- hotelling.test(iso1+iso2 ~community, data  = Six_m_Moorea_siber_data, perm =  TRUE)
fit_Moorea_sixm_perm

# Extract 20 
Twenty_m_Moorea_siber_data <- subset (Moorea_siber_data, group == 20)
# with permutations
fit_Moorea_twentym_perm <- hotelling.test(iso1+iso2 ~community, data  = Twenty_m_Moorea_siber_data, perm =  TRUE)
fit_Moorea_twentym_perm

Forty_m_Moorea_siber_data <- subset (Moorea_siber_data, group == 40)
# with permutations
fit_Moorea_fortym_perm <- hotelling.test(iso1+iso2 ~community, data  = Forty_m_Moorea_siber_data, perm =  TRUE)
fit_Moorea_fortym_perm

Sixty_m_Moorea_siber_data <- subset (Moorea_siber_data, group == 60)
# with permutations
fit_Moorea_sixtym_perm <- hotelling.test(iso1+iso2 ~community, data  = Sixty_m_Moorea_siber_data, perm =  TRUE)
fit_Moorea_sixtym_perm
#### Moorea_2 ####


#### Tikehau_2 ####
# Extract 6 
Six_m_Tikehau_siber_data <- subset (Tikehau_siber_data, group == 6)
# with permutations
fit_Tikehau_sixm_perm <- hotelling.test(iso1+iso2 ~community, data  = Six_m_Tikehau_siber_data, perm =  TRUE)
fit_Tikehau_sixm_perm

# Extract 20 
Twenty_m_Tikehau_siber_data <- subset (Tikehau_siber_data, group == 20)
# with permutations
fit_Tikehau_twentym_perm <- hotelling.test(iso1+iso2 ~community, data  = Twenty_m_Tikehau_siber_data, perm =  TRUE)
fit_Tikehau_twentym_perm

Forty_m_Tikehau_siber_data <- subset (Tikehau_siber_data, group == 40)
# with permutations
fit_Tikehau_fortym_perm <- hotelling.test(iso1+iso2 ~community, data  = Forty_m_Tikehau_siber_data, perm =  TRUE)
fit_Tikehau_fortym_perm

# Tikehau does not have enough replicates of 60 m

#### Tikehau_2 ####

#### Rangiroa_3 ####
# Extract 6 
Six_m_Rangiroa_siber_data <- subset (Rangiroa_siber_data, group == 6)
# with permutations
fit_Rangiroa_sixm_perm <- hotelling.test(iso1+iso2 ~community, data  = Six_m_Rangiroa_siber_data, perm =  TRUE)
fit_Rangiroa_sixm_perm

# Extract 20 
Twenty_m_Rangiroa_siber_data <- subset (Rangiroa_siber_data, group == 20)
# with permutations
fit_Rangiroa_twentym_perm <- hotelling.test(iso1+iso2 ~community, data  = Twenty_m_Rangiroa_siber_data, perm =  TRUE)
fit_Rangiroa_twentym_perm

Forty_m_Rangiroa_siber_data <- subset (Rangiroa_siber_data, group == 40)
# with permutations
fit_Rangiroa_fortym_perm <- hotelling.test(iso1+iso2 ~community, data  = Forty_m_Rangiroa_siber_data, perm =  TRUE)
fit_Rangiroa_fortym_perm

# Rangiroa does not have enough replicates of 60 m
#### Rangiroa_3 ####

### End of Hotelling tests





### Finally, Relative Degree of Heterotrophy (RDH) Williams et al 2018, quick visualization. 
ggplot(Isotopes_Pocillopora, aes(x=factor(Depth), y=Delta13C_Pol_Zoox)) +
  geom_boxplot() + facet_grid(cols = vars(Island_Site)) + stat_summary (fun=mean, geom="point", shape=18, color="red", size=2) + theme_bw()  + 
  ylab ("Delta δ13C Polyps - Zoox") + xlab ("Depth (m)") + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                                                                 axis.text = element_text(size=10, colour="black"),
                                                                 axis.title = element_text(size=11, face="bold", colour="black")) 

# RDH = "Delta13C_Pol_Zoox" with depth. Additional tests on this RDH can be found in Script: "Pocillopora_All_Multivar_Bayesian.R"  



# Many other tests have been done, but for ease use of the data and interpretation, these were excluded from the manuscript and I have not included them here. 
# These extra tests are: 
# Additional tests from SIBER package: https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html and Jackson et al 2011 and Conti-Jerpe et al 2020
# Additional tests as suggested by Radice et al 2019
# Bayesian modelling of the isotope variables with depth. Also available in Script: "Pocillopora_All_Multivar_Bayesian.R" 

# The whole script with all different tests can be shared upon request. 


