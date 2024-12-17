library(dplyr)
## Calculate fold change -----
#Epi
mydata.epi <- df.epi.revert %>% #Create mydata.df with new columns for averages and fold change
  as_tibble()

#These didn't work inside mutate() / Works fine now?
epi.far.AVG <- rowMeans(select(mydata.epi, ends_with("EF")))
epi.near.AVG = rowMeans(select(mydata.epi, ends_with("EN")))
#muscle.far.AVG = rowMeans(select(mydata.scaled.df, ends_with("MF")))
#muscle.near.AVG = rowMeans(select(mydata.scaled.df, ends_with("MN")))


mydata.epi <- mydata.epi %>%
  mutate(far_AVG = epi.far.AVG,
         near_AVG = epi.near.AVG,
         LogFC = epi.far.AVG - epi.near.AVG) %>%
  mutate(Protein_names = rownames(df.epi.revert)) %>%
  relocate(Protein_names) %>%
  mutate_if(is.numeric, round, 2) 

mydata.df.epi <- mydata.epi %>%
  dplyr::select(Protein_names, far_AVG, near_AVG, LogFC) %>% #Pull out two columns
  dplyr::arrange(desc(LogFC)) #Arrange based on descending order of fold change

mydata.df.epi.f <- mydata.df.epi %>%
  filter(LogFC != 0)

data.no0.epi <- mydata.epi %>%
  filter(LogFC != 0) %>%
  mutate(pvalue = t.test(select(ends_with("EN")), 
                         y = select(ends_with("EF"))))


#Differentially expressed proteins Log2 FC greater than 1/-1
DEPe <- mydata.df.epi %>%
  filter(LogFC >= 1.0 | LogFC <= -1.0)

#Muscle
mydata.m <- df.m.revert %>%
  as_tibble()
#m.far.AVG <- rowMeans(select(mydata.m, ends_with("MF")))
mydata.m <- mydata.m %>%
  mutate(far_AVG = rowMeans(select(mydata.m, ends_with("MF"))),
         near_AVG = rowMeans(select(mydata.m, ends_with("MN"))),
         LogFC = far_AVG - near_AVG) %>%
  mutate(Protein_names = rownames(df.m.revert)) %>%
  relocate(Protein_names) %>%
  mutate_if(is.numeric, round, 2) 

mydata.df.m <- mydata.m %>%
  dplyr::select(Protein_names, far_AVG, near_AVG, LogFC) %>%
  dplyr::arrange(desc(LogFC))

mydata.df.m.f <- mydata.df.m %>%
  filter(LogFC != 0)

#Differentially expressed proteins Log2 FC greater than 1/-1
DEPm <- mydata.df.m %>%
  filter(LogFC >= 1.0 | LogFC <= -1.0)
  
#Static table
gt(DEPe)
gt(DEPm) #Something weird is happening with averages. Exact opposite values for near and far. 

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