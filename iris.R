this.dir <- dirname("parent.frame(2)$ofile")
setwd(this.dir)
library(FactoMineR)
require(FactoClass)
iris=data.frame(iris)
tiris=scale(iris[,1:4])
source("Kohonenprep.R")
source("Kohonennet.R")
source("plotcarte1.R")
source("evaluation.R")
base=matrix(rep(tiris,500),ncol=dim(tiris)[2], byrow=T)
prof=kohonenqualigo(25476,4,4,0.04,0.01,2.99,0.65,base,dim(base)[1])
m=kohonenqualiclass(prof,tiris,dim(tiris)[1])
irclust=cbind(iris,m)[order(m),]
table(irclust[,5],irclust[,6])
nmi(table(irclust[,5],irclust[,6]))
purity(table(irclust[,5],irclust[,6]))
# L'image de la carte
grillecarte(6,11,3,irclust[,5],irclust[,6])
par(xpd=TRUE)
nb=3
ncol=seq(0,240,length.out=nb)
legend("topright", inset=c(-0.16,0.2), title="Iris Type", c("setosa","versicolo","virginica"), pch=15,col=hcl(ncol,120,85),cex=0.55)
#clustering xyf
disj=as.matrix(acm.disjonctif(iris[,4:5]))
disj=disj[,(ncol(disj)-2):ncol(disj)]
tiris=cbind(tiris=scale(iris[,1:4]),disj)
w=c(rep(1/(4*4),4),rep(3/(4*3),3))
prof=kohonenqualigo.weight(17,3,3,0.04,0.01,2.99,0.65,tiris,dim(tiris)[1],w)
m=kohonenqualiclass.weight(prof,tiris,dim(tiris)[1],w)
irclust=cbind(iris,m)[order(m),]
table(irclust[,5],irclust[,6])
nmi(table(irclust[,5],irclust[,6]))
purity(table(irclust[,5],irclust[,6]))
#pCA
prc <- prcomp(iris[,1:2], center=TRUE, scale=TRUE)
newData <- tiris[,1:2] %*% prc$rotation 
base=matrix(rep(newData,500),ncol=dim(tiris)[2], byrow=T)
prof=kohonenqualigo(25476,4,4,0.04,0.01,2.99,0.65,base,dim(base)[1])
m=kohonenqualiclass(prof,newData,dim(tiris)[1])
irclust=cbind(iris,m)[order(m),]
table(irclust[,5],irclust[,6])
nmi(table(irclust[,5],irclust[,6]))
purity(table(irclust[,5],irclust[,6]))