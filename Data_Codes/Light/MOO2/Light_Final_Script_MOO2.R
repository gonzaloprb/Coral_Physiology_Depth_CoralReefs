
###### DEEPHOPE   -  "Light Reading and Graphs"

# Read straight from CSV format from logger. Important the name is "___.csv" format ("IslandSite-Prof-light")

# Only need to introduce, working directory [17], names of files [20], time of start and end reading [34]; Id_Site, archipelago, island, site, names to save plot and database [67, 70]. All things are " " 

# Information loggers: 
# // DEFI Series // CSV File // Firmware Version 1.02 // Software Version 1.03 [Head] SondeName=DEFI2-L SondeNo=0B2P015 SensorType=Q2B0 SensorType2=0102 SensorType3=60 Channel=2 Interval=900 SampleCnt=684 StartTime=2018/08/21 10:00:00 EndTime=2018/08/28 12:45:00 StopTime=2018/08/28 12:56:57 Status=01000 DepAdjRho=1.0250 ImmersionEN=1 CoefDate=2018/05/23 Immersion_Effect=1.39 Ch1=+8.768397e+03,-1.351438e-01,+0.000000e+00,+0.000000e+00,+0.000000e+00,+0.000000e+00,+0.000000e+00,+0.000000e+00, Ch2=+1.169766e-02,+8.109043e-04,+0.000000e+00,+0.000000e+00,+0.000000e+00,+0.000000e+00,+0.000000e+00,+0.000000e+00, [Item]

rm (list = ls())

#Necessary packages (in any case execute below)
library(ggplot2)
library(RColorBrewer)
library (hms)

Id_Site <- c("MOO-2")
Archipelago <- c("Society")
Island <- c("Moorea")
Site <- c("2")



# Introduce file names with the directory of the folder between " " where there are the CSV files
file.names = c("Data_Codes/Light/MOO2/MO2-20-light.csv","Data_Codes/Light/MOO2/MO2-40-light.csv","Data_Codes/Light/MOO2/MO2-60-light.csv","Data_Codes/Light/MOO2/MO2-90-light.csv","Data_Codes/Light/MOO2/MO2-120-light.csv")
prof = c(20,40,60,90,120)
df_all=data.frame()

# Introduce limits of time, to get to the depth and to take it out. 
# Date_Time initial >> "2018-09-04 12:00:00"
# Date_Time final >>  "2018-09-07 10:40:00"

for(x in 1:length(file.names)){
  prof_file = read.csv(file = file.names[x], header = T, dec = ".", sep = ",", skip = 24)
  df = data.frame(prof_file)
  colnames(df) <- c("Date_Time", "Quantum","Battery")
  df <- df[,-4]
  df$Date_Time <- as.POSIXct(df$Date_Time, format="%Y/%m/%d %H:%M:%S")
  df <- subset(df, Date_Time > "2018-09-04 12:00:00" & Date_Time < "2018-09-07 10:40:00")
  df$Depth <- prof[x]
  df_all = rbind(df_all,df)
  # str (df)
}

# Plot the graph
df_all$Depth <- as.factor (df_all$Depth)

myColors<- c('20' = 'darkslategray1', '40' = 'deepskyblue', '60' = 'blue', '90' = 'darkblue','120' = 'black')
names(myColors) <- levels(df$Depth)
colScale <- scale_colour_manual(name = "Depth",values = myColors)


ggplot(data = df_all, aes (y=Quantum, x = Date_Time, group = Depth, color = Depth)) +
  geom_point(size = 0.5) + geom_line () + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +    colScale +
  scale_x_datetime( date_breaks = "5 hour") +
  ylim (0,500)  + guides(colour = guide_legend(reverse=F)) + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Delete night values
df_all$Time <- format(df_all$Date_Time,"%H:%M:%S")
df_all_day <- subset(df_all, Time < "17:30:00" & Time > "07:00:00")

# Mean of light 
tapply (df_all_day$Quantum, INDEX = df_all_day$Depth, FUN = mean)


# Normalize respect to the light at 20 m 
df_mean_time <- aggregate(Quantum ~  Depth + Time,  df_all_day, mean)

# For unique island and time, divide Quantum by Quantum of 20 
normalize_light = data.frame()
for(i in unique(df_mean_time$Time)){
  new <-subset (df_mean_time, Time == i)
  # print(new)
  new$Quantum <- 1 - (new$Quantum / new$Quantum [new$Depth == 20])
  normalize_light = rbind (normalize_light, new)
  # print (new_recrue)
}


# Change format of time 
normalize_light$Time <- as.hms (normalize_light$Time,quiet = F, roll = F)

# Relative plot
ggplot(data = normalize_light, aes (y=Quantum, x = Time, group = Depth, color = Depth)) +
  geom_point(size = 0.5) + geom_line()  + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  colScale + scale_y_reverse() + 
  guides(colour = guide_legend(reverse=F)) + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Plot of the daytime only 

df_all_day <- aggregate(Quantum ~  Depth + Time,  df_all_day, mean)


df_all_day$Time <- as.POSIXct(df_all_day$Time, format="%H:%M:%S")

ggplot(data = df_all_day, aes (y=Quantum, x = Time, group = Depth, color = Depth)) +
  geom_point(size = 0.5) + geom_line (size = 0.05) +  geom_smooth(method = "loess", span = 0.1, se = F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +    colScale +
  scale_x_datetime( date_breaks = "2 hour") +
  ylim (0,250)  + 
  guides(colour = guide_legend(reverse=F)) + 
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Mean of the loss of light 
## Values of the loss of light per depth
tapply (normalize_light$Quantum, INDEX = list (normalize_light$Depth), FUN = mean)
Light_val <- 1- tapply (normalize_light$Quantum, INDEX = list (normalize_light$Depth), FUN = mean)


# Introduce the name of Id_Site, Archipelago, Island and Site
df_all$Id_Site <- "MO-2"
df_all$Archipelago <- "Society"
df_all$Island <- "Moorea"
df_all$Site <- "2"

# Re order them 
df_all <- df_all[,c(5,6,7,8,4,1,2,3)]




# Save the Light database
write.csv(df_all, "Data_Codes/Light/MOO2/MO2_Light_data.csv")




myColors<- c( '20' = 'aquamarine2', '40' = 'deepskyblue', '60' = 'blue', '90' = 'navyblue','120' = 'black') # '6' = 'cornsilk'
names(myColors) <- levels(df$Depth)
colScale <- scale_colour_manual(name = "Depth",values = myColors)

ggplot(data = df_all, aes (y=Quantum, x = Date_Time, group = Depth, color = Depth)) +
  geom_point(size = 0.1) + geom_line (size = 0.05) +  geom_smooth(method = "loess", span = 0.15, se = F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +    colScale +
  scale_x_datetime( date_breaks = "12 hour", date_labels = "%d %b %R") + xlab ("") + ylab ("Light (μmol m-2 s-1)") +
  ylim (0,250)  + 
  guides(colour = guide_legend(reverse=F)) + 
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Beer Lambert conversion
# Starting from df_all_day
df_all_day$Time <- strftime(df_all_day$Time, format="%H:%M:%S")


# For unique island and time, divide Quantum by Quantum of 20 
lambert_light = data.frame()
ref_prof = data.frame()
for(i in unique(df_all_day$Time)){
  new <-subset (df_all_day, Time == i)
  new20 <- subset (new, Depth == 20) # Take out reference value
  ref_prof = rbind (ref_prof,new20)
  new <- subset (new, Depth != 20)
  new$Depth <- as.numeric (as.character( new$Depth))
  for (j in unique (new$Depth)){
    new2 <- subset (new,  Depth == j)
    # print(new)
    new2$k <- -1 / (j*log(new20$Quantum [new20$Depth == 20]/new2$Quantum [new2$Depth == j]))
    lambert_light = rbind (lambert_light, new2)
    # print (new_recrue)
  }}


ref_prof$k <- NA

# Lambert light
lambert_light <- rbind (ref_prof,lambert_light)


# Complete database adding quantum and k NA

# Necessary to make --> Par(20m) * exp (K(40m) * 20)
Six_Meters <- data.frame(
  Depth=6,
  Time= unique (lambert_light$Time))

# Apply equation
Six_Meters$Quantum <- lambert_light$Quantum [lambert_light$Depth == 20] / exp (lambert_light$k [lambert_light$Depth == 40] * 20)  

Six_Meters$k <- NA

lambert_light_Final <- rbind (Six_Meters, lambert_light)


# Measure the index
# For unique island and time, divide Quantum by Quantum of 6
Relative_Index_Light = data.frame()
for(i in unique(lambert_light_Final$Time)){
  new <-subset (lambert_light_Final, Time == i)
  # print(new)
  new$IndexLoss <- (new$Quantum / new$Quantum [new$Depth == 6])
  Relative_Index_Light = rbind (Relative_Index_Light, new)
  # print (new_recrue)
}


# Transform to Index_light 
Index_light <- Relative_Index_Light [,c("Depth", "Time","Quantum", "k","IndexLoss")]

# To be in percentatge and as below!
Index_light$IndexLoss <- Index_light$IndexLoss*100

# Index loss using Beer Lambert Equation

Index_light <- subset(Index_light, Time < "16:00:00" & Time > "09:00:00")


#table2$Depth2<-as.numeric(table2$Depth2)
Index_light$Depth<-as.numeric(Index_light$Depth)

Median <- aggregate (Quantum ~  Depth,  Index_light, median)
Mean <- aggregate (Quantum ~  Depth,  Index_light, mean) 
Max <- aggregate (Quantum ~  Depth,  Index_light, max)
Min <- aggregate (Quantum ~  Depth,  Index_light, min)
colnames (Median) <- c ("Depth","Median")
colnames (Mean) <- c ("Depth","Mean")
colnames (Max) <- c ("Depth","Max")
colnames (Min) <- c ("Depth","Min")



Index_light <- merge (Index_light, Median, by = "Depth")
Index_light <- merge (Index_light, Mean, by = "Depth")
Index_light <- merge (Index_light, Max, by = "Depth")
Index_light <- merge (Index_light, Min, by = "Depth")

Index_light$Id_Site <- "MOO2"
Index_light$Archipelago <- "Society"
Index_light$Island <- "Moorea"
Index_light$Site <- "2"


myColors<- c('20' = 'darkslategray1', '40' = 'deepskyblue', '60' = 'blue', '90' = 'darkblue','120' = 'black')
names(myColors) <- levels(df$Depth)
colScale <- scale_colour_manual(name = "Depth",values = myColors)

# scale_colour_gradient(myColors)+
ggplot(data=Index_light, aes(y= Depth , x =Quantum, color =IndexLoss))+
  geom_point(size = 0.5)+
  geom_smooth( color="grey", fill="grey", size=0.2, alpha=0.2)+
  scale_y_reverse(name ="Depth (m)", limits=c(120,0),breaks = c(120,90, 60, 40, 20, 6,0))+
  scale_x_continuous (position = "top") + 
  scale_colour_gradient2(low ="blue4", mid = "deepskyblue",  high = "yellow", midpoint = 50)+
  labs(x = "Irradiance (μmol m-2 s-1)", y="Depth (m)")+
  geom_point (aes(y= Depth , x =Median, size = 1)) + 
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())




ggplot(data=Index_light, aes(y= Depth , x =Quantum, xmin = Min, xmax = Max,color =IndexLoss))+
  geom_errorbar (size = 1,width = 2) +  geom_point (aes(y= Depth , x =Quantum, size = 0.5)) + 
  #  geom_smooth(method = "loess", color="grey", fill="grey", size=0.2, alpha=0.2)+
  geom_path (aes(y= Depth,x =Mean,color=IndexLoss),size = 2,width = 2,alpha = 0.3) +
  geom_point (aes(y= Depth , x =Mean, size = 0.8)) +
  scale_y_reverse(name ="Depth (m)", limits=c(120,0),breaks = c(120,90, 60, 40, 20, 6,0))+
  scale_x_continuous (position = "top") + 
  scale_colour_gradient2(low ="blue4", mid = "deepskyblue",  high = "skyblue", midpoint = 50)+
  labs(x = "Irradiance (μmol m-2 s-1)", y="Depth (m)")+
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())






# geom_boxplot
Index_light2 <- Index_light
Index_light2$Depth <- Index_light2$Depth * (-1)
# Create five colours
colours <- c("black","blue4", "blue","skyblue", "cyan","yellow")
ggplot(Index_light2, aes(x = Quantum, y = factor (Depth), color = colours)) + 
  geom_boxplot(color= colours) + geom_point(size = 0.5, color = "grey") +
  scale_x_continuous(position = "top",limits=c(0,300),breaks = c(0,100,200,300)) + 
  scale_y_discrete(name ="Depth (m)", limits = c("-120","-90","-60","-40","-20","-6"), breaks = c(0,-6,-20,-40,-60,-90, -120)) + 
  labs(x = "Irradiance (μmol m-2 s-1)", y="Depth (m)")+
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(Index_light2, aes(x = IndexLoss, y = factor (Depth), color = colours)) + 
  geom_boxplot(color= colours) + geom_point(size = 0.5, color = "grey") +
  scale_x_continuous(position = "top",limits=c(0,100),breaks = c(0,25,50,75,100)) + 
  scale_y_discrete(name ="Depth (m)", limits = c("-120","-90","-60","-40","-20","-6"), breaks = c(0,-6,-20,-40,-60,-90, -120)) + 
  labs(x = "Irradiance Index Loss (μmol m-2 s-1)", y="Depth (m)")+
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Absolute value graph
ggplot(data = Index_light, aes(y= Depth , x =Quantum, color =Quantum))+
  geom_point (size = 0.5, colour = "yellow")  + geom_smooth(colour = "yellow", span = 1, method = "loess") +
  scale_x_continuous(position = "top",name ="Light", limits=c(0,400), breaks = c(0,100,200,300,400)) +
  scale_y_reverse(name ="Depth (m)", limits=c(120,0),breaks = c(120,90, 60, 40, 20,6, 0))+
  guides(colour = guide_legend(reverse=F)) + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Indexx loss
ggplot(data = Index_light, aes(y= Depth , x =IndexLoss, color =IndexLoss))+
  geom_point (size = 2, colour = "yellow") + geom_path(colour = "yellow") + geom_smooth(colour = "yellow", span = 4, method = "gam") +
  scale_x_continuous(position = "top",name ="Quantum_Index loss (%)", limits=c(0.,101), breaks = c(0,25,50,75,100)) +
  scale_y_reverse(name ="Depth (m)", limits=c(120,0),breaks = c(120,80, 60, 40, 20, 0))+
  guides(colour = guide_legend(reverse=F)) + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(data = Index_light, aes (y=IndexLoss, x = Depth)) +
  geom_point (size = 2, colour = "yellow") +
  geom_smooth(colour = "grey", span = 1, method = "loess") + 
  scale_x_continuous(name ="Depth (m)", limits=c(0,100), breaks = c(6,20,40,60,90)) +
  scale_y_continuous(name ="Quantum_Index loss (%)", limits=c(0,100), breaks = c(0,25,50,75,100)) +
  labs(x = "Quantum_Index loss (%)",  y = "Temperature",  title = "Light with depth") +
  theme_classic() + theme(strip.background = element_blank())


Index_light_Loss <- aggregate(IndexLoss ~  Depth,  Index_light, mean)
write.csv (Index_light_Loss, paste("Data_Codes/Light/MOO2/",Id_Site,"_", "Index_Light_loss.csv", sep = ""))
# Database used for figure in Physiology manuscript


