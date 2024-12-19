##Load packages -----
library(tidyverse)
install.packages("BBmisc")
library(BBmisc)
library(DT)
library(plotly)
library(gt)
library(readxl)
library(limma)
library(edgeR)

##Pre-processing -----

df <- read_excel("/Users/evanvallenari/Library/CloudStorage/OneDrive-UniversitetetiOslo/Proteomics/LCM-FFPE_Prot/Data/Proteomics_data_FFPE_CFCN-MFMN.xlsx")

#Filter out peptides only IDd by site, reverse, and potential contaminants
df.filtered <- filter(df, is.na(`Only identified by site`) &
                        is.na(Reverse) &
                        is.na(`Potential contaminant`)
                      )

#Not clear why filtering for != '+' didn't work.
df.filtered.select <- df.filtered %>%
  select(`Protein names`, `Gene names`, contains('LFQ'))

#Change column names to sample_group
colnames(df.filtered.select) <- c("Protein_names", "Gene_names", "1562_EF", "1563_EN", "1564_MF", "1565_MN", "1570_EF", "1571_EN", "1572_MF", "1573_MN", "1578_EF", "1579_EN", "1580_MF", "1581_MN", "1586_EF", "1587_EN", "1588_MF", "1589_MN", "1594_EF", "1595_EN", "1596_MF", "1597_MN", "1602_EF", "1603_EN", "1604_MF", "1605_MN", "1610_EF", "1611_EN", "1612_MF", "1613_MN", "1618_EF", "1619_EN", "1620_MF", "1621_MN", "1622_EF", "1623_EN", "1624_MF", "1625_MN")

df.epi <- select(df.filtered.select, contains("E")) #Better to select all "E" to keep everything in order

df.m <- select(df.filtered.select, contains("M"))

##Epi -----
#Transform data (rows as columns) and keep rows numeric. Apply prot names as column names and remove first row (was column names, not numerical, so we have to remove)
df.epi.t <- t(sapply(df.epi, as.numeric))

colnames(df.epi.t) <- paste0(df.epi$Protein_names, " (", df.epi$Gene_names, ")")
df.epi.t <- df.epi.t[-1,] #x2

#Log2 transformation - add constant of 1 to avoid -inf from 0s
df.epi.t.log2 <- log2(df.epi.t)

#Uses function scale to Z-score normalize within protein column, confirm by taking mean/sd
df.epi.t.scale <- scale(df.epi.t.log2) %>%
  as.data.frame()
mean(df.epi.t.scale$`Keratin, type I cytoskeletal 14`) #Confirm z-score normalization - mean = 0, st.dev = 1
sd(df.epi.t.scale$`Keratin, type I cytoskeletal 14`)

#Transform (Each column is a sample, each row is a protein)
df.epi.revert <- t(df.epi.t.scale)
df.epi.r.log2 <- t(df.epi.t.log2_plus1) #without z-norm
#Omit NAs
df.scale.epi.NAomit <- na.omit(df.epi.r.log2)

df.epi.NAomit <- na.omit(df.epi.revert)
##PCA.Epi -----
#Heirarchical clustering
distance.epi <- dist(t(df.epi.NAomit), method = "euclidean")
clusters.epi <- hclust(distance.epi, method = "complete")
plot(clusters.epi)

#PCA
pca.res.epi <- prcomp(t(df.epi.NAomit), scale. = F, retx=T)
summary(pca.res.epi) #Prints variance summary 
screeplot(pca.res.epi) #Screeplot is standard way to view eigenvalues for each PCA, not great for presenting
pc.var.epi <- pca.res.epi$sdev^2 #Captures eigenvalues from PCA result
pc.per.epi <- round(pc.var.epi/sum(pc.var.epi)*100, 1) #Uses eigenvalues to calculate percentage variance accounted for by each PC
pc.per.epi
pca.res.epi.df <- as_tibble(pca.res.epi$x)

#Read in target data for epi-- identify variables of interest
targets.epi <- read.delim("Doc/targets_epi.txt")
group.epi <- factor(targets.epi$group)
sampleLabels.epi <- as.character(targets.epi$sample)
#samples <- colnames(df.epi.scale.NAomit) #Used to get order of columns to edit targets.txt file 
#could have been avoided by selecting "E" instead of "EF" and then "EN"

#Plotting with ggplot2 -- 2D XY plot
ggplot(pca.res.epi.df) + #Identify dataframe
  aes(x=PC1, y=PC2, label=sampleLabels.epi, color=group.epi) + #Aesthetics, X, Y, labels, color by group
  geom_point(size=4) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per.epi[1],"%",")")) + #Create x axis label by pasting together without spaces
  ylab(paste0("PC2 (",pc.per.epi[2],"%",")")) +
  labs(title="PCA plot",
       subtitle="Epithelium") +
  coord_fixed() +
  theme_bw()

ggsave("Res/Scaled_PCA_epi-update.png")

##Muscle -----
#Transform data (rows as columns) and keep rows numeric. Apply prot names as column names and remove first row (was column names, not numerical, so we have to remove)
df.m.t <- t(sapply(df.m, as.numeric))

colnames(df.m.t) <- df.m$Protein_names
df.m.t <- df.m.t[-1,]

df.m.t.log2 <- log2(df.m.t)

#Uses function scale to Z-score normalize within protein column, confirm by taking mean/sd
df.m.t.scale <- scale(df.m.t.log2) %>%
  as.data.frame()
mean(df.m.t.scale$`Keratin, type I cytoskeletal 14`) #Confirm z-score normalization - mean = 0, st.dev = 1
sd(df.m.t.scale$`Keratin, type I cytoskeletal 14`)

#Transform (Each column is a sample, each row is a protein)
df.m.revert <- t(df.m.t.scale)

#Omit rows with NA 
#is.na(df.m.revert) <- sapply(df.m.revert, is.infinite)
#df.m.scale.NAomit <- na.omit(df.m.revert)

#Change introduced NAs back to 0 instead of omitting
df.m.revert[is.na(df.m.revert)] <- 0

##PCA.m -----
#Heirarchical clustering
distance.m <- dist(t(df.m.revert), method = "euclidean")
clusters.m <- hclust(distance.m, method = "complete")
plot(clusters.m)
#samples <- colnames(df.m.scale.NAomit)

#PCA
pca.res.m <- prcomp(t(df.m.revert), scale. = F, retx=T)
summary(pca.res.m) #Prints variance summary 
screeplot(pca.res.m) #Screeplot is standard way to view eigenvalues for each PCA, not great for presenting
pc.var.m <- pca.res.m$sdev^2 #Captures eigenvalues from PCA result
pc.per.m <- round(pc.var.m/sum(pc.var.m)*100, 1) #Uses eigenvalues to calculate percentage variance accounted for by each PC
pc.per.m
pca.res.m.df <- as_tibble(pca.res.m$x)

#Read in target data for m-- identify variables of interest
targets.m <- read.delim("Doc/targets_m.txt", sep = ",")
group.m <- factor(targets.m$group)
sampleLabels.m <- as.character(targets.m$sample)

#Plotting with ggplot2 -- 2D XY plot
ggplot(pca.res.m.df) + #Identify dataframe
  aes(x=PC1, y=PC2, label=sampleLabels.m, color=group.m) + #Aesthetics, X, Y, labels, color by group
  geom_point(size=4) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per.m[1],"%",")")) + #Create x axis label by pasting together without spaces
  ylab(paste0("PC2 (",pc.per.m[2],"%",")")) +
  labs(title="PCA plot",
       subtitle="Muscle") +
  coord_fixed() +
  theme_bw()

ggsave("Res/Scaled_PCA_m.png")
