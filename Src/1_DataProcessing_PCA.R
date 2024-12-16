## Load/Install packeages -----
install.packages("tidyverse")
install.packages("readxl")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library(tidyverse)
library(readxl)
library(edgeR)


## Read in excel file and tidy -----
df <- read_excel("/Users/evanvallenari/Library/CloudStorage/OneDrive-UniversitetetiOslo/Proteomics/LCM-FFPE_Prot/Data/Proteomics_data_FFPE_CFCN-MFMN.xlsx")

#df with LFQ intensities
df.select <- df %>% 
  select(`Protein names`, `Gene names`, `Fasta headers`, contains('LFQ'))


## PCA Analysis -----

library(tidyverse)
library(DT)
library(plotly)
library(gt)

df.select
df.select.matrix <- df.select %>% #Create data matrix from data frame
  select(contains('LFQ')) %>% #Select columns containing LFQ values
  data.matrix() #As data matrix
  row.names(df.select.matrix) <- df.select$`Protein names` #w/ Protein names as row names
#Rename columns
  colnames(df.select.matrix) <- c("1562 EF", "1563 EN", "1564 MF", "1565 MN", "1570 EF", "1571 EN", "1572 MF", "157 MN", "1578 EF", "1579 EN", "1580 MF", "1581 MN", "1586 EF", "1587 EN", "1588 MF", "1589 MN", "1594 EF", "1595 EN", "1596 MF", "1597 MN", "1602 EF", "1603 EN", "1604 MF", "1605 MN", "1610 EF", "1611 EF", "1612 MF", "1613 MN", "1618 EF", "1619 EN", "1620 MF", "1621 MN", "1622 EF", "1623 EN", "1624 MF", "1625 MN")
  
#Read in target data -- identify variables of interest
targets <- read.delim("Doc/targets.txt")
group <- factor(targets$group)

#Hierarchical clustering
distance <- dist(t(df.select.matrix), method = "euclidean")
clusters <- hclust(distance, method = "complete")
plot(clusters, labels = paste(group, clusters$labels))
?plot
?hclust
#PCA
df.select.NAomit.matrix <- na.omit(df.select.matrix) #prcomp cannot run with NA variables. na.omit omits all NAs
pca.res <- prcomp(t(df.select.NAomit.matrix), scale. = F, retx=T)
summary(pca.res) #Prints variance summary 
screeplot(pca.res) #Screeplot is standard way to view eigenvalues for each PCA, not great for presenting
pc.var <- pca.res$sdev^2 #Captures eigenvalues from PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) #Uses eigenvalues to calculate percentage variance accounted for by each PC
pc.per
pca.res.df <- as_tibble(pca.res$x)

sampleLabels <- as.character(targets$sample)
#Plotting with ggplot2 -- 2D XY plot
ggplot(pca.res.df) + #Identify dataframe
  aes(x=PC1, y=PC2, label=sampleLabels, color=group) + #Aesthetics, X, Y, labels, color by group
  geom_point(size=4) +
  #geom_label(nudge_x = 20, nudge_y = 20) + #Haven't figured this out
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + #Create x axis label by pasting together without spaces
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

#Small multiples PCA plot
pca.res.df4 <- pca.res$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group) 
  pca.res.df4.rename <- rename(pca.res.df4, "PC1 (56.8%)" = "PC1", "PC2 (20.4%)" = "PC2", "PC3 (7.9%)" = "PC3", "PC4 (4.4%)" = "PC4") # Rename columns 'new = old' format
?rename
pca.pivot <- pivot_longer(pca.res.df4.rename, #Identify datafram to be pivoted
                          cols = 'PC1 (56.8%)':'PC4 (4.4%)', #Column names to be stored as single variable
                          names_to = "PC", #Column name of that variable
                          values_to = "loadings") #Name of new column storing data values
#labelsSM.PCA <- as.character(c(paste0("PC1 (",pc.per[1],"%",")"), paste0("PC2 (",pc.per[2],"%",")"), paste0("PC3 (",pc.per[3],"%",")"), paste0("PC4 (",pc.per[4],"%",")")))


ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) +
  geom_bar(stat='identity') +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       ) +
  theme_bw() +
  coord_flip() +
  
  
ggsave("Res/PCA_small_multiples.png")

##Separate muscle from epithelium - epi
df_lfq <- as.data.frame(df.select.NAomit.matrix)
Epi_matrix <- select(df_lfq, contains("E")) %>%
  as.matrix()


pca.res.epi <- prcomp(t(Epi_matrix), scale. = F, retx=T)
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

#Plotting with ggplot2 -- 2D XY plot
ggplot(pca.res.epi.df) + #Identify dataframe
  aes(x=PC1, y=PC2, label=sampleLabels.epi, color=group.epi) + #Aesthetics, X, Y, labels, color by group
  geom_point(size=4) +
  #geom_label(nudge_x = 20, nudge_y = 20) + #Haven't figured this out
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per.epi[1],"%",")")) + #Create x axis label by pasting together without spaces
  ylab(paste0("PC2 (",pc.per.epi[2],"%",")")) +
  labs(title="PCA plot",
       subtitle = "Epithelium"
       ) +
  coord_fixed() +
  theme_bw()

ggsave("Res/PCA_epi.png")
##Separate muscle from epithelium - muscle
df_lfq <- as.data.frame(df.select.NAomit.matrix)
Muscle_matrix <- select(df_lfq, contains("M")) %>%
  as.matrix()


pca.res.m <- prcomp(t(Muscle_matrix), scale. = F, retx=T)
summary(pca.res.m) #Prints variance summary 
screeplot(pca.res.m) #Screeplot is standard way to view eigenvalues for each PCA, not great for presenting
pc.var.m <- pca.res.m$sdev^2 #Captures eigenvalues from PCA result
pc.per.m <- round(pc.var.m/sum(pc.var.m)*100, 1) #Uses eigenvalues to calculate percentage variance accounted for by each PC
pc.per.m
pca.res.m.df <- as_tibble(pca.res.m$x)

#Read in target data for epi-- identify variables of interest
targets.m <- read.delim("Doc/targets_muscle.txt")
group.m <- factor(targets.m$group)
sampleLabels.m <- as.character(targets.m$sample)

#Plotting with ggplot2 -- 2D XY plot
ggplot(pca.res.m.df) + #Identify dataframe
  aes(x=PC1, y=PC2, label=sampleLabels.m, color=group.m) + #Aesthetics, X, Y, labels, color by group
  geom_point(size=4) +
  #geom_label(nudge_x = 20, nudge_y = 20) + #Haven't figured this out
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per.m[1],"%",")")) + #Create x axis label by pasting together without spaces
  ylab(paste0("PC2 (",pc.per.m[2],"%",")")) +
  labs(title="PCA plot",
       subtitle = "Muscle"
       ) +
  coord_fixed() +
  theme_bw()

ggsave("Res/PCA_muscle.png")
