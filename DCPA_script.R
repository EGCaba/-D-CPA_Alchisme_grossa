#Checking my r-version / Revisando la versión de R
R.version.string
install.packages('installr')
library(installr)
updateR()

#Installing adegenet package with dependencies / Instalar adegenet con las dependencias.
install.packages('dplyr')
install.packages('poppr')
install.packages('adegenet',dep=TRUE)
install.packages('scales')
scaleinstall.packages(c('dplyr','hierfstat','reshape2','s','RcolorBrewer'))
update.packages(ask=FALSE) #actualizar a la versión más reciente
packageVersion('adegenet')
library(poppr)
library(adegenet)
library(ape)
library(pegas)
library(seqinr)
library(ggplot2)
library(dplyr)
library(hierfstat)
library(reshape2)
library(RColorBrewer)
packageDescription('poppr', fields= 'Version')
packageDescription("adegenet", fields ="Version" )
help.search('Hardy-weinberg')
?adegenet
class ? genlight
data("nancycats")
is.genind(nancycats)
nancycats
getwd()
setwd("C:/Users/edgar/OneDrive/Documentos/R/DPCA_Alchisme_grossa)


#Loading snps in the console / Cargar bases de datos de snps a la consola
obj1 <- read.genepop('populations.snps.gen')
obj1
is.genind(obj1)
tab(obj1)
obj1$tab[1:5,1:10]

#calculate missing data per locus /cálculo de datos en blanco por locus
miss_loci= propTyped(obj1, by = "loc")
print(miss_loci[which(miss_loci < 0.60)])
barplot(miss_loci, ylim = c(0,1),ylab='complete genotypes (proportion)',xlab='locus',las = 2, cex.names = 0.4,main="Proportion of data (loci) with more than 60%")

#Remove loci with more than 20% of data missing /Quitar loci con más del 20% de datos pérdidos
remv_gen = missingno(obj1, type = 'loci', cutoff =  0.20)

#calculate unique genotypes / Calcular genotipos únicos
mlg(remv_gen)
mlg(rem_in)
mlg(obj1)
#calculate missing data per individual /cálculo de datos en blanco por individuo.
miss_ind= propTyped(obj1, by = "ind")
print(miss_ind[which(miss_ind < 0.60)])
barplot(miss_ind, ylim = c(0,1),ylab='complete genotypes (proportion)',xlab='individuals',las = 2, cex.names = 0.4,main="Proportion of data for individualswith more than 60%")

#remove individuals with more than 20% of missing data / Remover individuos con más de 20% de datos en blanco
rem_in = missingno(obj1, type = 'geno', cutoff =  0.40)

#Identify duplicates sample / Identificar duplicados

duplicados = mlg.id(obj1)
for (i in  duplicados) {
  if (length(obj1[i]) > 1){
    print(i)
  }
}

## Resumen de la tabla de datos.
obj1
table(obj1$loc.fac) #print the number of alleles per locus / Desplegar el número de alelos por locus
summary(obj1$pop) #print the sample size per population / Desplegar el tamaño muestral por población
private_alleles(obj1) %>% apply(MARGIN = 1, FUN = sum) #print the number of private alleles per population across loci/ Desplegar el número de alelos privados por población por todos los loci

##calculate mean oberved heterozygosity per site / calcular la heterocigosidad observada promedio por sitio

basic = basic.stats(obj1,diploid = T)
ho_basic = apply(basic$Ho,MARGIN = 2,FUN = mean, na.rm = T ) %>% round(digits = 2)
ho_basic

##calculate mean expected heterozygosity per site /calcular el promedio de la heterocigosidad esperada por sitio

he_basic = apply(basic$Hs,MARGIN = 2,FUN = mean, na.rm = T ) %>% round(digits = 2)
he_basic


##Visualize heterozygosity per site / Visualizar la hterocigosidad por sitio.

He_plot_df = data.frame(Site = c(1:5),Ho = ho_basic, He= he_basic) %>% melt(id.vars="Site")
He_plot_df

#Heterozygosity barplot

custom_theme = theme(
  axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, face = "bold"),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12),
  axis.title.x = element_blank(),
  axis.line.y = element_line(size = 0.5),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  panel.grid = element_blank(),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, size = 15, face="bold")
)

# Italic label
hetlab.o = expression(italic("H")[o])
hetlab.e = expression(italic("H")[e])

ggplot(data = He_plot_df, aes(x = Site, y = value, fill = variable))+ 
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), colour = "black" )+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.50))+
  scale_fill_manual(values = c("forestgreen", "gray47"), labels = c(hetlab.o, hetlab.e))+
  ylab("Heterozygosity")+
  ggtitle("Alchisme grossa population")+
            custom_theme
rlang::last_error()


##Inbreeding coefficient / Coeficiente de endogamia

apply(basic$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)

## FST
obj1_fst = genet.dist(obj1,method="WC84") %>% round(digits = 3)
obj1_fst


#Desire order of the labels /Orden Deseado de las etiquetas
lab_order = as.matrix(c("pop 1","pop 2","pop 3","pop 4", "pop 5"))
lb = as.vector(lab_order)
class(lab_order)
#Change order of rows and columns/ Cambiar el orden de las columnas y las filas.
fst.mat = as.matrix(obj1_fst)
class(lab_order)
fst.mat1 = fst.mat[obj1_fst, ]
fst.mat2 = fst.mat1[ ,obj1_fst]

# Create a data.frame
ind = which(upper.tri(fst.mat2), arr.ind = TRUE)
fst.df = data.frame(Site1 = dimnames(fst.mat2)[[2]][ind[,2]],
                    Site2 = dimnames(fst.mat2)[[1]][ind[,1]],
                    Fst = fst.mat2[ ind ])

# Keep the order of the levels in the data.frame for plotting 
fst.df$Site1 = factor(fst.df$Site1, levels = unique(fst.df$Site1))
fst.df$Site2 = factor(fst.df$Site2, levels = unique(fst.df$Site2))

# Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

# Print data.frame summary
fst.df %>% str
## 'data.frame':    15 obs. of  3 variables:
##  $ Site1: Factor w/ 5 levels "Brd","Pad","Vig",..: 1 2 2 3 3 3 4 4 4 4 ...
##  $ Site2: Factor w/ 5 levels "Ber","Brd","Pad",..: 1 1 2 1 2 3 1 2 3 4 ...
##  $ Fst  : num  0.007 0.025 0.008 0.064 0.038 0.018 0.174 0.171 0.161 0.112 ...

# Fst italic label
fst.label = expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid = max(fst.df$Fst) / 2

# Plot heatmap
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.05, 0.10, 0.15))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)
  )




##Perform the principal components analysis / Realizar el análisis de componentes principales

#Replace missing data with mean allele frequencies / Reemplazar datos en blancos con el promedio de las frencuencias alélicas.

x = tab(obj1, NA.method = "mean")

#Perform PCA / Hacer PCA

pca1 = dudi.pca(x,scannf = F, scale = F, nf = 3 )
pca1
# Analyse how much percent of genetic variance is explained by each axis/ Analizar que tanto porcentaje de varianza genética es explicada por cada X

percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1), main="Eigenvectors")


##Visualiza PCA
# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(obj1)

# Add a column with the site IDs
ind_coords$Site = obj1$pop

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define colour palette
cols = brewer.pal(nPop(obj1), "Set1")

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) 
)

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("Lobster PCA")+
  # custom theme
  ggtheme

##DAPC


# Perform cross validation to find the optimal number of PCs to retain in DAPC
set.seed(123)
x = tab(obj1, NA.method = "mean")
crossval = xvalDapc(x, obj1$pop, result = "groupMean", xval.plot = TRUE)

#Number of PCs with best stats (lower score = better)
crossval$`Root Mean Squared Error by Number of PCs of PCA`

crossval$`Number of PCs Achieving Highest Mean Success`

crossval$`Number of PCs Achieving Lowest MSE`

numPCs = as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)


# Run a DAPC using site IDs as priors
dapc1 = dapc(obj1, obj1$pop, n.pca = numPCs, n.da = 3)
dapc1

# Analyse how much percent of genetic variance is explained by each axis
percent = dapc1$eig/sum(dapc1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,60),
        names.arg = round(percent, 1),main="Explain how much genetic variance is explained by each axe")


#DAPC visualization.
# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(dapc1$ind.coord)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(obj1)

# Add a column with the site IDs
ind_coords$Site = obj1$pop

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define colour palette
cols = brewer.pal(nPop(obj1), "Set2")

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("Alchisme grossa DAPC")+
  # custom theme
  ggtheme

#Perform DPCA

obj1
grp <- find.clusters(obj1, max.n.clust = 10, n.pca = NULL )
grp
names(grp)
is.genind(obj1)
?find.clusters.genind()
getwd()
setwd("C:/Users/edgar/OneDrive/Documentos/R/DPCA_Alchisme_grossa") 
