profiles_36p <- read.delim("D:/2013_01_XAS_SEDIMENTS/profiles_36p.txt", dec=",")
attach(profiles_36p)
#Multivariate analysis
#PCA
pca <- prcomp (~ Position + Prepeak + WL + Intensity, data = profiles_36p, center = TRUE , scale = TRUE)

summary(pca)

# scree plot and biplot
plot(pca)
biplot(pca)

# scores and loadings
scores <-  pca$x
loadings <- pca$rotation
plot(scores[,1:2], pch=paste(1:5))

library(ggplot2)
qplot(x=scores[,1], y= scores[,2]) + geom_text(label=1:36) 
qplot(x=scores[,1], y= scores[,2]) + geom_text(label=Sample) 





samples_17 <- read.delim("D:/2013_01_XAS_SEDIMENTS/samples_17.txt", dec=",")

attach(samples_17)
#Multivariate analysis
#PCA
pca1 <- prcomp (~ Position + Prepeak + WL + Intensity, data = samples_17, center = TRUE , scale = TRUE)

summary(pca1)

# scree plot and biplot
plot(pca1)
biplot(pca1)

# scores and loadings
scores1 <-  pca1$x
loadings1 <- pca1$rotation
plot(scores1[,1:2])

library(ggplot2)
qplot(x=scores1[,2], y= scores1[,3]) + geom_text(label=1:17) 
qplot(x=scores1[,2], y= scores1[,3]) + geom_text(label=Sample) 




s01 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_01.txt", quote="\"")
s02 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_02.txt", quote="\"")
s03 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_03.txt", quote="\"")
s04 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_04.txt", quote="\"")
s05 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_05.txt", quote="\"")
s06 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_06.txt", quote="\"")
s07 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_07.txt", quote="\"")
s08 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_08.txt", quote="\"")
s09 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_09.txt", quote="\"")
s10 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_10.txt", quote="\"")
s11 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_11.txt", quote="\"")
s12 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_12.txt", quote="\"")
s13 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_13.txt", quote="\"")
s14 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_14.txt", quote="\"")
s15 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_15.txt", quote="\"")
s16 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_16.txt", quote="\"")
s17 <- read.table("D:/2013_01_XAS_SEDIMENTS/TXT_samples_17/sample_17.txt", quote="\"")


library(hyperSpec)
spc01 <- new("hyperSpec", spc = s01[,2], wavelength = s01[,1])
spc02 <- new("hyperSpec", spc = s02[,2], wavelength = s02[,1])
spc03 <- new("hyperSpec", spc = s03[,2], wavelength = s03[,1])
spc04 <- new("hyperSpec", spc = s04[,2], wavelength = s04[,1])
spc05 <- new("hyperSpec", spc = s05[,2], wavelength = s05[,1])
spc06 <- new("hyperSpec", spc = s06[,2], wavelength = s06[,1])
spc07 <- new("hyperSpec", spc = s07[,2], wavelength = s07[,1])
spc08 <- new("hyperSpec", spc = s08[,2], wavelength = s08[,1])
spc09 <- new("hyperSpec", spc = s09[,2], wavelength = s09[,1])
spc10 <- new("hyperSpec", spc = s10[,2], wavelength = s10[,1])
spc11 <- new("hyperSpec", spc = s01[,2], wavelength = s01[,1])
spc12 <- new("hyperSpec", spc = s02[,2], wavelength = s02[,1])
spc13 <- new("hyperSpec", spc = s03[,2], wavelength = s03[,1])
spc14 <- new("hyperSpec", spc = s04[,2], wavelength = s04[,1])
spc15 <- new("hyperSpec", spc = s05[,2], wavelength = s05[,1])
spc16 <- new("hyperSpec", spc = s06[,2], wavelength = s06[,1])
spc17 <- new("hyperSpec", spc = s07[,2], wavelength = s07[,1])



all <- collapse (spc01,spc02,spc03,spc04,spc04,spc06,spc07,spc08,spc09,spc10,spc11,spc12, spc13,spc14,spc15, spc16, spc17)
all01 <- orderwl(all)

all02 <- spc.loess (all01 , seq (7.04, 7.8, 0.001))
all02 <- spc.loess (all01 , seq (7.05, 7.8,  length.out = 400))

all03 <- sweep (all02, 2, mean, "-")

#PCA
pca <- prcomp (~ spc, data = all03$., center = TRUE , scale = TRUE)
# results to hyperSpec object
scores <- decomposition (all03, pca$x, label.wavelength = "PC", label.spc = "score / a.u.")
loadings <- decomposition (all03, t(pca$rotation), scores = FALSE, label.spc = "loading I / a.u.")
# scree plot
plot(pca)
#plot of loadings
plot (loadings [1:3], stacked = TRUE)


#extract to vectors scores values of pc1 pc2 and pc3
pc1 = scores [[,,1]]
pc2 = scores [[,,2]]
pc3 = scores [[,,3]]

# variable coding various colours for every cluster
cl1 = test5$clusters

# plot pc1 versus pc3 points in cluster colours
plot(pc1,pc2)
qplot(x= pc1[,1], y= pc2[,1]) + geom_text(label=1:17) 





plot(pca)
# distribution of scores on the map
plotmap (scores [,,1])
plotmap (scores [,,2])
plotmap (scores [,,3])
#plot of loadings
plot (loadings [1:3], stacked = TRUE, col = cols)


#extract to vectors scores values of pc1 pc2 and pc3
pc1 = scores [[,,1]]
pc2 = scores [[,,2]]
pc3 = scores [[,,3]]

# variable coding various colours for every cluster
cl1 = test5$clusters

# plot pc1 versus pc3 points in cluster colours
plot(pc1,pc2, col=cl1)


#export every spectrum to single file (automaticaly named) in long format - two columns with many rows - first column wavenumber second column absorbance no headers


for(i in 1:10)
{ 
  name = paste(i,".","csv", sep = "")
  write.txt.long(test[i,,],file=name, cols = c(".wavelength","spc"),col.names =FALSE)
}


#earlier one could cut from hyperSpec object only spectra from cluster1 for example with split function
clusters <- split(test5, test5$clusters)
#one can call the object with command, splitting can be reversed by rbind
clusters$Cl_02

























#wl.range = c(730~1200, 1400~1680, 3500~3780), xoffset = c(150,1700)
#testowanie baseline correction
#1 spc.fit.poly fits a polynomial baseline of the given order
#2 spc.fit.poly.below tries to find appropriate support points for the baseline #iteratively
#3 package baseline many methods

bl1 <- spc.fit.poly.below(test,poly.order=5)
test1 = test - bl1
plot(test1,"spcmeansd", wl.reverse =TRUE, col="red")
title(main="FIT POLY BELOW 5", col.main="blue")


#smoothing interpolation 
#spc.bin (averaging every by points)
#spc.loess(applies R loess function)

#namawiamy funkcje plotmap do wspolpracy definiujac col.regions
plotmap(test4, col.regions = pal(20))

write.table (mcaa,file="caa.txt", sep = "\t", row.names=F,col.names=F)

# change PC values to factors with cut function
scores[[,,1]] = cut(scores[[,,1]],c(-20,10,90))

#mozemy tez sklonic do wspolpracy plotspc definiujac col
plotspc(test [550:624,,820~2870], col=pal(10))

#Przyklady wizualizacji
plot(sample(test,600),"voronoi")




#preprocessing
bl <- spc.fit.poly.below (test)
test1 <- test - bl
plotmap(test1 [,,1550], col.regions = pal(200))
plotspc(test1 [550:624,,820~2870], col=pal(10))

#tak wyglada srednia +/-SD przed obrobka
plot(test,"spcmeansd", wl.reverse =TRUE)
# obrobka
#smooth
test1 = spc.loess(test, seq (604,4000,4))
plot(test1,"spcmeansd", wl.reverse =TRUE)

#odjecie 5 percentyla (wspolna komponenta)
test2 = sweep(test1, 2, apply(test1, 2, quantile, 0.05), "-")
plot(test2,"spcmeansd", wl.reverse =TRUE)

#normalizacja, podzielenie przez srednia
#test6 <- sweep(test1, 1, apply(test1,1,mean), "/")
#plot(test6,"spcmeansd", wl.reverse =TRUE)

#wyciecie fragmentow widma z obszarami: polisacharydy, amidy i aniony
test4 <- test2 [ , , c(730~1200, 1400~1680, 3500~3780)]
plot(test4,"spcmeansd", wl.reverse =TRUE, wl.range = c(730~1200, 1400~1680, 3500~3780), xoffset = c(150,1700))


#mapki intesywnosci dla trzech przykladowych pasm

plotmap(test4 [,,1031], col.regions = pal(200))
plotmap(test4 [,,1542], col.regions = pal(200))
plotmap(test4 [,,3620], col.regions = pal(200))


#wybieramy pozadana liczbe kolorow z ulubionej palety w przykladzie 4 z palety BuGn (odcienie zielonego)
colors <- brewer.pal(8, "Reds")
#tworzymy palete kolorow
pal <- colorRampPalette(colors)


plotmap(test4 [,,820~1200], col.regions = pal(50))
plotmap(test4 [,,1480~1680], col.regions = pal(50))
plotmap(test4 [,,3550~3750], col.regions = pal(50))


#wybieramy pozadana liczbe kolorow z ulubionej palety w przykladzie 4 z palety BuGn (odcienie zielonego)
colors <- brewer.pal(3, "Set1")
#tworzymy palete kolorow
pal <- colorRampPalette(colors)

#baseline correction
#test6 <-test4
#library(baseline)
#bl <- baseline (test6 [[]], method = "rollingBall", wm=500, ws=500)
#test6 [[]] <- getCorrected(bl)
#plot(test6,"spcmeansd", wl.reverse =TRUE, wl.range = c(730~1200, 1400~1680, 3500~3780), xoffset = c(150,1700))

#example of visualisation with ggplot2

df < as.t.df (apply(test5,2,mean_pm_sd))
ggplot(df, aes (x = .wavelength))+geom_ribbon(aes (ymin = mean.minus.sd, ymax = mean.plus.sd), alpha =0.25)+geom_line(aes(y=mean))




#HCA
#calculate distance between the spectra
dist <- dist ()
#construct dendrogram based on distance and linkage method
dendrogram <- hclust (dist, method ="ward")

#cut dendrogram for the expected number of clusters 
#add cluster membership to hyperSpec object
test5$clusters <- as.factor (cutree (dendrogram, k = 7))
#set colours and names for clusters
cols =pal5(7)
levels (test5$clusters) <- c ("Cl_01", "Cl_02", "Cl_03", "Cl_04", "Cl_05", "Cl_06", "Cl_07")
#map of the cluster location
print (plotmap (test5, clusters ~ x * y, col.regions = cols))
#dendrogram with marked cluster membership 
par (xpd = TRUE)
plot (dendrogram, labels = FALSE, hang = -1)
mark.dendrogram (dendrogram, test5$clusters, col = cols)
#mean spectra for clusters
cluster.means <- aggregate (test5, test5$clusters, mean_pm_sd)

plot(cluster.means, stacked = ".aggregate", fill = ".aggregate", col = cols, wl.reverse =TRUE,wl.range = c(3000~2700, 1800~670), xoffset = -800)

#voronoi type of cluster membership
plotvoronoi(sample(test5,600), clusters ~x *y, col.regions = cols)

#export cluster means to txt



