##### Takes Differenitally expressed proteins in epi near v far (DEprot) and 
##### raw data (df.raw) as input for heatmap and distance-based clustering 

##Packages -----
library(Hmisc)
library(tidyverse)
library(RColorBrewer)
library(gplots) 


##Data -----
df2 <- read.delim("Data/ENEF_DE-proteins_filt.txt") #Filtered list of DE proteins from perseus
#Removed patient 9 before running students t-test
#Patient 9 is still present in raw data file and needs to be removed #Line 61

DEprot2 <- df2 %>% 
  mutate(FC = as.numeric(Student.s.T.test.Difference.EN_EF)) %>% 
  mutate(p.value = Student.s.T.test.p.value.EN_EF) %>% 
  filter(Q.value < 0.01) %>% 
  filter(FC > 0.5 | FC < -0.5) %>% 
  select(Protein.names, Gene.names, Majority.protein.IDs, FC, p.value, Q.value)

head(DEprot)

#Removed 2 unnamed proteins 
DEprot2 <- DEprot2[-c(5,22),]

write.csv(DEprot, "Res/sigDEs_filt.txt")


filt = DEprot$Protein.names
nofilt = df$ProteinNames_nofilt
same <- filt %in% nofilt
samesum <- sum(same)
diff <- filt %nin% nofilt
diffsum <- sum(diff)

DEprot.diff <- DEprot %>% 
  mutate(diff = diff) %>% 
  filter(diff == "TRUE")

write.csv(DEprot.diff, "Res/sigDEs_filt_diff.txt")


##Heatmap-----
DEproteins <- as.vector(DEprot2$Protein.names)

## Read in excel file and tidy -----
df.raw <- read_excel("Data/Proteomics_data_FFPE_CFCN-MFMN.xlsx")

df.select.rawInt <- df.raw %>% 
  select(`Protein names`, contains('Intensity')) %>% 
  select(!contains("LFQ")) %>% 
  filter(`Protein names` %in% DEproteins) %>% 
  select(c(`Protein names`, contains('slot')))

#Set col names with ID and group
colnames(df.select.rawInt) <- c("ProteinNames", "1562_EF", "1563_EN", "1564_MF", "1565_MN", "1570_EF", "1571_EN", "1572_MF", "1573_MN", "1578_EF", "1579_EN", "1580_MF", "1581_MN", "1586_EF", "1587_EN", "1588_MF", "1589_MN", "1594_EF", "1595_EN", "1596_MF", "1597_MN", "1602_EF", "1603_EN", "1604_MF", "1605_MN", "1610_EF", "1611_EN", "1612_MF", "1613_MN", "1618_EF", "1619_EN", "1620_MF", "1621_MN", "1622_EF", "1623_EN", "1624_MF", "1625_MN")

#Remove patient 9
df.select.rawInt <- df.select.rawInt[,1:33]

matrix.rawInt.E <- df.select.rawInt %>% 
  select(contains("E")) %>% 
  data.matrix()

row.names(matrix.rawInt.E) <- df.select.rawInt$ProteinNames #w/ Protein names as row names
#Remove column ProteinNames
matrix.rawInt.E <- matrix.rawInt.E[,-1]

## Cluster DE proteins ----
#Narrowed non-filtered DEGs (FC 0.5, padj 0.01)
clust_rows <- hclust(as.dist(1-cor(t(matrix.rawInt.E), method="pearson")), method="complete") 
clust_columns <- hclust(as.dist(1-cor(matrix.rawInt.E, method="spearman")), method="complete") #cluster columns by spearman correlation
#module_assign <- cutree(clust_rows, k=5)
module_assign <- cutree(clust_columns, k=2)
module_color <- rainbow(length(unique(module_assign)), start=0.1, end=0.9) 
module_color <- module_color[as.vector(module_assign)] 

## Heatmap -----

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))
#Narrow, non-filtered...
options(bitmapType = 'cairo')
heatmap <- heatmap.2(matrix.rawInt.E, 
                        Rowv=as.dendrogram(clust_rows), 
                        Colv=as.dendrogram(clust_columns),
                        #ColSideColors = module_color, #Pretty self-explanatory division
                        col=myheatcolors, scale='row',
                        density.info="none", trace="none",  
                        cexRow=1, cexCol=1.5, margins = c(7, 5),
                        key = TRUE, keysize = 1, labRow = NA
)

plot(clust_columns)

?heatmap.2

heatmap2 <- heatmap.2(matrix.rawInt.E,
                      col = myheatcolors,
                      )

dev.off()


##With all proteins and samples
## Read in excel file and tidy -----
df.raw <- read_excel("Data/Proteomics_data_FFPE_CFCN-MFMN.xlsx")

df.t <- df.raw %>% 
  select(`Protein names`, contains('Intensity')) %>% 
  select(!contains("LFQ")) %>% 
  select(c(`Protein names`, contains('slot')))

#Set col names with ID and group
colnames(df.t) <- c("ProteinNames", "1562_EF", "1563_EN", "1564_MF", "1565_MN", "1570_EF", "1571_EN", "1572_MF", "1573_MN", "1578_EF", "1579_EN", "1580_MF", "1581_MN", "1586_EF", "1587_EN", "1588_MF", "1589_MN", "1594_EF", "1595_EN", "1596_MF", "1597_MN", "1602_EF", "1603_EN", "1604_MF", "1605_MN", "1610_EF", "1611_EN", "1612_MF", "1613_MN", "1618_EF", "1619_EN", "1620_MF", "1621_MN", "1622_EF", "1623_EN", "1624_MF", "1625_MN")

matrix.t2 <- df.t %>% 
  select(contains("E")) %>% 
  data.matrix()

matrix.t <- data.matrix(df.t)

row.names(matrix.t2) <- df.raw$`Protein names` #w/ Protein names as row names
#Remove column ProteinNames
matrix.t2 <- matrix.t2[,-1]

## Cluster DE proteins
#Narrowed non-filtered DEGs (FC 0.5, padj 0.01)
 
clust_columns.t <- hclust(as.dist(1-cor(matrix.t2, method="spearman")), method="complete") #cluster columns by spearman correlation

plot(clust_columns.t)



##DAVID ----
#KEGG pathway cytoskeleton in muscle cells
kegg <- read_tsv("Res/GOterm_protLists/KEGG_pathway_cytoskeletonInMuscleCells.txt")

df.kegg <- df.raw %>% 
  select(`Protein names`, `Majority protein IDs`, contains('Intensity')) %>% 
  select(!contains("LFQ")) %>% 
  filter(`Majority protein IDs` %in% kegg$ID) %>% 
  select(c(`Protein names`, `Majority protein IDs`, contains('slot')))

#Set col names with ID and group
colnames(df.kegg) <- c("ProteinNames", "1562_EF", "1563_EN", "1564_MF", "1565_MN", "1570_EF", "1571_EN", "1572_MF", "1573_MN", "1578_EF", "1579_EN", "1580_MF", "1581_MN", "1586_EF", "1587_EN", "1588_MF", "1589_MN", "1594_EF", "1595_EN", "1596_MF", "1597_MN", "1602_EF", "1603_EN", "1604_MF", "1605_MN", "1610_EF", "1611_EN", "1612_MF", "1613_MN", "1618_EF", "1619_EN", "1620_MF", "1621_MN", "1622_EF", "1623_EN", "1624_MF", "1625_MN")

#Remove patient 9
df.kegg <- df.kegg[,1:33]

matrix.E <- df.kegg %>% 
  select(contains("E")) %>% 
  data.matrix()

row.names(matrix.E) <- df.kegg$ProteinNames #w/ Protein names as row names
#Remove column ProteinNames
matrix.E <- matrix.E[,-1]

dimnames(matrix.E)
rownames(matrix.E)[2] <- "Fibronectin"
#Heatmap

clust_rows <- hclust(as.dist(1-cor(t(matrix.E), method="pearson")), method="complete") 
clust_columns <- hclust(as.dist(1-cor(matrix.E, method="spearman")), method="complete")
plot(clust_columns)

myheatcolors <- rev(brewer.pal(name="RdBu", n=11))

heatmap <- heatmap.2(matrix.E, 
                     Rowv=as.dendrogram(clust_rows), 
                     Colv=as.dendrogram(clust_columns),
                     #ColSideColors = module_color, #Pretty self-explanatory division
                     col=myheatcolors, scale='row',
                     density.info="none", trace="none",  
                     cexRow=1, cexCol=1.5, margins = c(7, 15),
                     key = TRUE, keysize = 1, 
)
