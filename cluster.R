# Sample dataset - University Admission Statistics

par(ask=FALSE)

# load relevant libraries
library(TeachingDemos)
library(fpc)
library(cluster)

x <- read.csv("Universities.csv", as.is=T)
# Minkowski distance - either Euclidean or Manhattan -  seems most appropriate 
# for the data,  which includes non-redundant continuous variables. 
# If the data were in the form of counts we would look toward canberra distance 
# instead, for instance. 
# Yes, will need to standardize data to perform distance calculations. 
# All variables are on separate scales. 

round(sqrt(apply(x[,-1],2,var)),2)

x.scale <- scale(x[,-1])
rownames(x.scale) = x[,1]
faces(x.scale, labels=x$University, ncol=5)
stars(x.scale, labels=as.character(x$University), ncol=5)
stars(x.scale, labels=as.character(x$University), ncol=5, draw.segments=T)

source("http://www.reuningscherer.net/STAT660/R/CSQPlot.r.txt")
CSQPlot(x.scale[,1:6], label="University Admission Data")
boxplot(x.scale)

## 2 
# Appears to be four clusters in the first method, which is supported by 
# the hclus_eval as well.
# The clusters also make "sense" when thinking about the schools included 
# together in the dendograms and other visualizations.  E.g. Yale and Harvard 
# are paired together, as are Texas A&M and Penn State.  
# Four clusters also captures the "outliers" of Caltech and Johns Hopkins. 
# In the second method, however, three clusters makes the most sense. 
# Using four splits Caltech and Johns Hopkins but combines Columbia and 
# Georgetown, which seem to have less in common. We are have less confidence
# in this method compared to the first since the cluster divisions fall in 
# in between Universities on the "lowest branches." e.g. Purdue and Texas A&M.

# Method 1
dist1 <- dist(x.scale, method = "euclidean", p=2)
clust1<-hclust(dist1,method="ward.D")
plot(clust1,labels=x$University,cex=0.5,xlab="",ylab="Distance"
     ,main="Clustering for Universities")
rect.hclust(clust1,k=5)

cuts1 = cutree(clust1,k=5)
cuts1

clusplot(x.scale, cuts1, color=TRUE, shade=TRUE, labels=2, lines=0,
         main="University Three Cluster Plot, Ward's Method, First two PC")

plotcluster(x.scale, cuts1, main="Three Cluster Solution in DA Space",
            xlab="First Discriminant Function", ylab="Second Discriminant Function")

# Evaluate Number of Clusters
source("http://reuningscherer.net/stat660/R/HClusEval.R.txt")
hclus_eval(x.scale, dist_m = 'euclidean', clus_m = 'ward.D', plot_op = T)

# Method 2 
library(vegan)
dist2 <- vegdist(x.scale, method="manhattan", upper=T)
clust2 <- hclust(dist2, method="complete")
plot(clust2, labels= x[,1], cex=0.5, xlab="", ylab="Distance", 
     main="Clustering for Universities")
rect.hclust(clust1,k=3)

cuts2 = cutree(clust2,k=3)
cuts2

clusplot(x.scale, cuts2, color=TRUE, shade=TRUE, labels=2, lines=0,
         main="Universities Two Cluster Plot, Complete Method, First two PC")
plotcluster(x.scale, cuts2, main="Two Cluster Solution in DA Space",
            xlab="First Discriminant Function", ylab="Second Discriminant Function")

source("http://reuningscherer.net/stat660/R/HClusEval.R.txt")
hclus_eval(x.scale, dist_m = 'manhattan', clus_m = 'complete', plot_op = T)

## 3 
# As mentioned above there appears to be three to five clusters. Attempting
# to use more than five leads to unusual groupings.  
# We plan to retain three clusters.    

## 4 
# We run k-means plot and several sum of squares plot below. In all cases
# the results indicate 3 clusters. This agrees with Method 1 used in the 
# hierarchical analysis from Q2. 

km1 <- kmeans(x.scale,centers=3, nstart=100)
pairs(x.scale, col=km1$cluster)
plot(x.scale, col = km1$cluster, pch = 19, frame = FALSE,
     main = "K-means with k = 3")
points(km1$centers, col = 1:5, pch = 8, cex = 3)

for (i in 1:3){
  print(paste("Universities in Cluster ",i))
  print(x$University[km1$cluster==i])
}

# Keep adding groups when the "Hartigan rule of thumb" (hrot) is > 10:

for (nclust in 1:20) {
  km.this <- kmeans(x.scale, centers=nclust, nstart=100)
  if (nclust >= 2) {
    klast <- nclust - 1
    hrot <- ( sum(km.last$withinss) / sum(km.this$withinss) - 1 ) *
      (nrow(x.scale) - klast - 1)
    cat(nclust, "-cluster solution HROT: ", hrot, "\n", sep="")
  }
  km.last <- km.this
}
# makes the case for 3 clusters

km.3 <- kmeans(x.scale, centers=3, nstart=100)
pairs(x.scale, col=km.3$cluster)

# Elbow Method
k.max <- 15 # Maximal number of clusters
wss <- sapply(1:k.max, 
              function(k){kmeans(x.scale, k, nstart=100 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
abline(v = 3, lty =2)

#kdata is just normalized input dataset
kdata=x.scale
n.lev=15  #set max value for k

# Calculate the within groups sum of squared error (SSE) for the number of cluster solutions selected by the user
wss <- rnorm(10)
while (prod(wss==sort(wss,decreasing=T))==0) {
  wss <- (nrow(kdata)-1)*sum(apply(kdata,2,var))
  for (i in 2:n.lev) wss[i] <- sum(kmeans(kdata, centers=i)$withinss)}

# Elbow method indicates 3 clusters 

# Calculate the within groups SSE for 250 randomized data sets (based on the original input data)
k.rand <- function(x){
  km.rand <- matrix(sample(x),dim(x)[1],dim(x)[2])
  rand.wss <- as.matrix(dim(x)[1]-1)*sum(apply(km.rand,2,var))
  for (i in 2:n.lev) rand.wss[i] <- sum(kmeans(km.rand, centers=i)$withinss)
  rand.wss <- as.matrix(rand.wss)
  return(rand.wss)
}

rand.mat <- matrix(0,n.lev,250)

k.1 <- function(x) { 
  for (i in 1:250) {
    r.mat <- as.matrix(suppressWarnings(k.rand(kdata)))
    rand.mat[,i] <- r.mat}
  return(rand.mat)
}

# Same function as above for data with < 3 column variables
k.2.rand <- function(x){
  rand.mat <- matrix(0,n.lev,250)
  km.rand <- matrix(sample(x),dim(x)[1],dim(x)[2])
  rand.wss <- as.matrix(dim(x)[1]-1)*sum(apply(km.rand,2,var))
  for (i in 2:n.lev) rand.wss[i] <- sum(kmeans(km.rand, centers=i)$withinss)
  rand.wss <- as.matrix(rand.wss)
  return(rand.wss)
}

k.2 <- function(x){
  for (i in 1:250) {
    r.1 <- k.2.rand(kdata)
    rand.mat[,i] <- r.1}
  return(rand.mat)
}

# Determine if the data data table has > or < 3 variables and call appropriate function above
if (dim(kdata)[2] == 2) { rand.mat <- k.2(kdata) } else { rand.mat <- k.1(kdata) }

# Plot within groups SSE against all tested cluster solutions for actual and randomized data - 1st: Log scale, 2nd: Normal scale

xrange <- range(1:n.lev)
yrange <- range(log(rand.mat),log(wss))
plot(xrange,yrange, type='n', xlab='Cluster Solution', ylab='Log of Within Group SSE', main='Cluster Solutions against Log of SSE')
for (i in 1:250) lines(log(rand.mat[,i]),type='l',col='red')
lines(log(wss), type="b", col='blue')
legend('topright',c('Actual Data', '250 Random Runs'), col=c('blue', 'red'), lty=1)

# Again indicates 3

yrange <- range(rand.mat,wss)
plot(xrange,yrange, type='n', xlab="Cluster Solution", ylab="Within Groups SSE", main="Cluster Solutions against SSE")
for (i in 1:250) lines(rand.mat[,i],type='l',col='red')
lines(1:n.lev, wss, type="b", col='blue')
legend('topright',c('Actual Data', '250 Random Runs'), col=c('blue', 'red'), lty=1)

# Less clear but still appears to indicate between 3 and 4

# Calculate the mean and standard deviation of difference between SSE of actual data and SSE of 250 randomized datasets
r.sse <- matrix(0,dim(rand.mat)[1],dim(rand.mat)[2])
wss.1 <- as.matrix(wss)
for (i in 1:dim(r.sse)[2]) {
  r.temp <- abs(rand.mat[,i]-wss.1[,1])
  r.sse[,i] <- r.temp}
r.sse.m <- apply(r.sse,1,mean)
r.sse.sd <- apply(r.sse,1,sd)
r.sse.plus <- r.sse.m + r.sse.sd
r.sse.min <- r.sse.m - r.sse.sd

# Plot differeince between actual SSE mean SSE from 250 randomized datasets - 1st: Log scale, 2nd: Normal scale 

xrange <- range(1:n.lev)
yrange <- range(log(r.sse.plus),log(r.sse.min))
plot(xrange,yrange, type='n',xlab='Cluster Solution', ylab='Log of SSE - Random SSE', main='Cluster Solustions against (Log of SSE - Random SSE)')
lines(log(r.sse.m), type="b", col='blue')
lines(log(r.sse.plus), type='l', col='red')
lines(log(r.sse.min), type='l', col='red')
legend('topright',c('SSE - random SSE', 'SD of SSE-random SSE'), col=c('blue', 'red'), lty=1)

# Indicates 3

xrange <- range(1:n.lev)
yrange <- range(r.sse.plus,r.sse.min)
plot(xrange,yrange, type='n',xlab='Cluster Solution', ylab='SSE - Random SSE', main='Cluster Solutions against (SSE - Random SSE)')
lines(r.sse.m, type="b", col='blue')
lines(r.sse.plus, type='l', col='red')
lines(r.sse.min, type='l', col='red')
legend('topright',c('SSE - random SSE', 'SD of SSE-random SSE'), col=c('blue', 'red'), lty=1)

# Again indicates 3

## 5 
# We believe there are three clusters present based on the information provided above. The 
# visualizations and k-means analysis correspond to this number.  

