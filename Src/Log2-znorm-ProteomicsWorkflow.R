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
df.filtered.tv <- filter(df, is.na(`Only identified by site`)) #Not clear why filtering for != '+' didn't work. 
?filter
df.filtered.select <- df.filtered.tv %>%
  select(`Protein names`, contains('LFQ'))

colnames(df.filtered.select) <- c("1562 EF", "1563 EN", "1564 MF", "1565 MN", "1570 EF", "1571 EN", "1572 MF", "157 MN", "1578 EF", "1579 EN", "1580 MF", "1581 MN", "1586 EF", "1587 EN", "1588 MF", "1589 MN", "1594 EF", "1595 EN", "1596 MF", "1597 MN", "1602 EF", "1603 EN", "1604 MF", "1605 MN", "1610 EF", "1611 EF", "1612 MF", "1613 MN", "1618 EF", "1619 EN", "1620 MF", "1621 MN", "1622 EF", "1623 EN", "1624 MF", "1625 MN")

#Transform data (rows as columns) and keep rows numeric. Apply prot names as column names and remove first row (was column names, not numerical, so we have to remove)
df.t <- t(sapply(df.filtered.select, as.numeric))

colnames(df.t) <- df.filtered.select$`Protein names`
df.t <- df.t[-1,]

df.t2 <- log2(df.t)

#Uses function scale to Z-score normalize within protein column, confirm by taking mean/sd
df.t.scale <- scale(df.t2) %>%
  as.data.frame()
mean(df.t.scale$`Keratin, type I cytoskeletal 14`) #Confirm z-score normalization - mean = 0, st.dev = 1
sd(df.t.scale$`Keratin, type I cytoskeletal 14`)

#Transform (Each column is a sample, each row is a protein)
df.revert <- t(df.t.scale)

#Omit rows with NA 
is.na(df.t.scale) <- sapply(df.t.scale, is.infinite)
df.scale.NAomit <- na.omit(df.revert)

colnames(df.scale.NAomit) <- c("1562 EF", "1563 EN", "1564 MF", "1565 MN", "1570 EF", "1571 EN", "1572 MF", "157 MN", "1578 EF", "1579 EN", "1580 MF", "1581 MN", "1586 EF", "1587 EN", "1588 MF", "1589 MN", "1594 EF", "1595 EN", "1596 MF", "1597 MN", "1602 EF", "1603 EN", "1604 MF", "1605 MN", "1610 EF", "1611 EF", "1612 MF", "1613 MN", "1618 EF", "1619 EN", "1620 MF", "1621 MN", "1622 EF", "1623 EN", "1624 MF", "1625 MN")


##PCA -----
#Heirarchical clustering
distance <- dist(t(df.scale.NAomit), method = "euclidean")
clusters <- hclust(distance, method = "single")
plot(clusters)

#PCA
pca.res <- prcomp(t(df.scale.NAomit), scale. = F, retx=T)
summary(pca.res) #Prints variance summary 
screeplot(pca.res) #Screeplot is standard way to view eigenvalues for each PCA, not great for presenting
pc.var <- pca.res$sdev^2 #Captures eigenvalues from PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) #Uses eigenvalues to calculate percentage variance accounted for by each PC
pc.per
pca.res.df <- as_tibble(pca.res$x)


sampleLabels <- as.character(targets$sample) #Creates character string of sample labels
#Plotting with ggplot2 -- 2D XY plot
ggplot(pca.res.df) + #Identify dataframe
  aes(x=PC1, y=PC2, label=sampleLabels, color=group) + #Aesthetics, X, Y, labels, color by group
  geom_point(size=4) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + #Create x axis label by pasting together without spaces
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggsave("Res/Scaled_PCA_1v2.png")


#Small multiples PCA plot
pca.res.df2 <- pca.res$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group) 
pca.res.df2.rename <- rename(pca.res.df2, "PC1 (56.8%)" = "PC1", "PC2 (8.5%)" = "PC2", "PC3 (4.8%)" = "PC3", "PC4 (3.6%)" = "PC4") # Rename columns 'new = old' format
?rename
pca.pivot <- pivot_longer(pca.res.df2.rename, #Identify dataframe to be pivoted
                          cols = 'PC1 (56.8%)':'PC4 (3.6%)', #Column names to be stored as single variable
                          names_to = "PC", #Column name of that variable
                          values_to = "loadings") #Name of new column storing data values

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) +
  geom_bar(stat='identity') +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       ) +
  theme_bw() +
  coord_flip()

ggsave("Res/Scaled_PCA_SmallMultiples.png")


## Calculate fold change -----

mydata.scaled.df <- df.scale.NAomit %>% #Create mydata.df with new columns for averages and fold change
  as_tibble()
  
#These didn't work inside mutate()
epi.far.AVG <- rowMeans(select(mydata.scaled.df, ends_with("EF")))
epi.near.AVG = rowMeans(select(mydata.scaled.df, ends_with("EN")))
muscle.far.AVG = rowMeans(select(mydata.scaled.df, ends_with("MF")))
muscle.near.AVG = rowMeans(select(mydata.scaled.df, ends_with("MN")))

mydata.scaled.df <- 
  mutate(epi.far.AVG = epi.far.AVG,
         epi.near.AVG = epi.near.AVG,
         muscle.far.AVG = muscle.far.AVG,
         muscle.near.AVG = muscle.near.AVG,
         LogFC.epi = epi.far.AVG - epi.near.AVG,
         LogFC.muscle = muscle.far.AVG - muscle.near.AVG) %>%
  mutate_if(is.numeric, round, 2) 

#Separate out cancer and muscle
mydata.df.cancer <- mydata.scaled.df %>%
  dplyr::select(Protein.names, cancer.far.AVG, cancer.near.AVG, LogFC.cancer) %>% #Pull out two columns
  dplyr::arrange(desc(LogFC.cancer)) #Arrange based on descending order of fold change

mydata.df.muscle <- mydata.scaled.df %>%
  dplyr::select(Protein.names, muscle.far.AVG, muscle.near.AVG, LogFC.muscle) %>%
  dplyr::arrange(desc(LogFC.muscle))

#Static table
gt(mydata.df.cancer)

#Create interactive table
datatable(mydata.df.cancer,
          extensions = c('KeyTable', 'FixedHeader'),
          filter = 'top',
          options = list(keys = TRUE,
                         searchHighlight = TRUE,
                         pageLength = 10,
                         lengthMenu = c("10", "25", "50", "100")))

datatable(mydata.df.muscle,
          extensions = c('KeyTable', 'FixedHeader'),
          filter = 'top',
          options = list(keys = TRUE,
                         searchHighlight = TRUE,
                         pageLength = 10,
                         lengthMenu = c("10", "25", "50", "100")))

#Create scatter plot of average values
canplot <- ggplot(mydata.df.cancer) +
  aes(x=cancer.far.AVG, y=cancer.near.AVG,
      text = paste(mydata.df.cancer$Protein.names)) + #Adds protein names for mouseover tooltip in ggplotly
  geom_point(shape=16, size=1) +
  ggtitle("near vs. far") +
  theme_bw()

ggplotly(canplot)

muscleplot <- ggplot(mydata.df.muscle) +
aes(x=muscle.far.AVG , y=muscle.near.AVG,
    text = paste(mydata.df.muscle$Protein.names)) + #Adds protein names for mouseover tooltip in ggplotly
  geom_point(shape=16, size=1) +
  ggtitle("near vs. far") +
  theme_bw()

ggplotly(muscleplot)

##Volcano plot -----

group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.cancer <- voom(mydata.df.cancer, design, plot = TRUE)
