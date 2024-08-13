# Script to open data and save it all together. PhotoAutotropy and Morphology
# In the end, the Isotopes files are separated.

rm(list =ls())

# Open Photoautotrophy
Pachyseris_Photoautotrophy <- read.csv ("Data_Codes/Photoautotrophy/Pachyseris_Photoautotrophy.csv", header = T, dec = ",", sep = ",")
Pocillopora_Photoautotrophy <- read.csv ("Data_Codes/Photoautotrophy/Pocillopora_Photoautotrophy.csv", header = T, dec = ",", sep = ",")

# Open Morphology
Pachyseris_Morphology <- read.csv ("Data_Codes/Morphology/Pachyseris_Morphology.csv", header = T, dec = ",", sep = ",")
Pocillopora_Morphology <- read.csv ("Data_Codes/Morphology/Pocillopora_Morphology.csv", header = T, dec = ",", sep = ",")

# Open Isotopes
Pachyseris_Isotopes <- read.csv ("Data_Codes/Isotopes/Pachyseris_Isotopes.csv", header = T, dec = ",", sep = ",")
Pocillopora_Isotopes <- read.csv ("Data_Codes/Isotopes/Pocillopora_Isotopes.csv", header = T, dec = ",", sep = ",")


# For Pocillopora
# For Morphology
intersect(unique (Pocillopora_Photoautotrophy$ID_LAB),unique (Pocillopora_Morphology$ID_LAB))

keep <- intersect(unique (Pocillopora_Photoautotrophy$ID_LAB),unique (Pocillopora_Morphology$ID_LAB))
Pocillopora_Morphology <- Pocillopora_Morphology [Pocillopora_Morphology$ID_LAB %in% keep, ]
# For Isotopes
intersect(unique (Pocillopora_Photoautotrophy$ID_LAB),unique (Pocillopora_Isotopes$ID_LAB))

keep <- intersect(unique (Pocillopora_Photoautotrophy$ID_LAB),unique (Pocillopora_Isotopes$ID_LAB))
Pocillopora_Isotopes <- Pocillopora_Isotopes [Pocillopora_Isotopes$ID_LAB %in% keep, ]

# For Pachyseris
# For Morphology
intersect(unique (Pachyseris_Photoautotrophy$ID_LAB),unique (Pachyseris_Morphology$ID_LAB))

keep <- intersect(unique (Pachyseris_Photoautotrophy$ID_LAB),unique (Pachyseris_Morphology$ID_LAB))
Pachyseris_Morphology <- Pachyseris_Morphology [Pachyseris_Morphology$ID_LAB %in% keep, ]
# For Isotopes
intersect(unique (Pachyseris_Photoautotrophy$ID_LAB),unique (Pachyseris_Isotopes$ID_LAB))

keep <- intersect(unique (Pachyseris_Photoautotrophy$ID_LAB),unique (Pachyseris_Isotopes$ID_LAB))
Pachyseris_Isotopes <- Pachyseris_Isotopes [Pachyseris_Isotopes$ID_LAB %in% keep, ]


# Merge 
Physio_Pocillopora_all <- merge (Pocillopora_Photoautotrophy, Pocillopora_Morphology, by = c("ID_LAB","Island_Site","Depth") , all = T)
# Physio_Pocillopora_all <- merge (Physio_Pocillopora_all, Pocillopora_Isotopes, by = c("ID_LAB","Island_Site","Depth") , all = T)

Physio_Pachyseris_all <- merge (Pachyseris_Photoautotrophy, Pachyseris_Morphology, by = c("ID_LAB","Island_Site","Depth") , all = T)
# Physio_Pachyseris_all <- merge (Physio_Pachyseris_all, Pachyseris_Isotopes, by = c("ID_LAB","Island_Site","Depth") , all = T)


# Change the name of Index Loss Light
colnames(Physio_Pocillopora_all)[8] <- "Relative_Index_Light"
colnames(Physio_Pachyseris_all)[8] <- "Relative_Index_Light"

# Write databases
write_csv(Physio_Pocillopora_all, "Data_Codes/All_Features/Physio_Pocillopora_all.csv")
write_csv(Physio_Pachyseris_all, "Data_Codes/All_Features/Physio_Pachyseris_all.csv")
