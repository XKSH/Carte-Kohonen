#continu de cleve.R ligne 17
#La métrique de Shih(2010)
source("kohonennet_weight.R")
source("dist.R")
#Variance de population
pop.var=function(x){sum((x - mean(x))^2) / length(x)}
tw.heart=heart
tw.heart[,1:5]=scale(tw.heart[,1:5])
#Utiliser la métrique de Ming-Yi Shih(2010) pour données mixtes
nl=data.frame(lapply(tw.heart[,6:13], nlevels))
#Determiner l'attribut qualitative de base
base.attr=tw.heart[,5+which(nl[1,]==max(nl[1,]))[1]]
ssw=colSums(apply(data.frame(tw.heart[,1:5]),2,function(x) tapply(x,base.attr,pop.var))) 
#Determiner l'attribut numérique de base
base.num=tw.heart[,which(ssw==min(ssw))]
v=matrix(tapply(base.num,base.attr,mean))
disj=as.matrix(acm.disjonctif(heart[,6:13]))
burt=t(disj)%*%disj
tw.burt=burt[c(1:2,7:nrow(burt)),3:6]
nom=matrix(rep(rowSums(tw.burt),ncol(tw.burt)),ncol=ncol(tw.burt))+t(matrix(rep(t(colSums(tw.burt)),nrow(tw.burt)),ncol=nrow(tw.burt)))
nom=nom-tw.burt
a=tw.burt/nom
w=a%*%v
w=c(w[1:2],v,w[3:length(w)])
vm=t(disj*w)
vm=matrix(vm[vm!=0],nrow=8)
vm=t(vm)
# La métrique de shin
tw.heart[,6:13]=vm
tw.heart1=tw.heart[,c(1:13)]
tw.heart1=do.call(cbind,tw.heart1)
matd=matrix(rep(c(tw.heart1),60),ncol=dim(tw.heart1)[2], byrow=T)
nr=4;nc=4
prof=kohonenqualigo(17,nr,nc,0.4,0.1,3,3,matd[sample(nrow(matd)),],dim(matd)[1])
m=kohonenqualiclass(prof,tw.heart1,dim(tw.heart1)[1])
htw.clust=cbind(tw.heart,m)[order(m),]
table(htw.clust[,ncol(htw.clust)-1],htw.clust[,ncol(htw.clust)])
nmi(table(htw.clust[,ncol(htw.clust)-1],htw.clust[,ncol(htw.clust)]))
purity(table(htw.clust[,ncol(htw.clust)-1],htw.clust[,ncol(htw.clust)]))
#Repondération
w=wt[909,]#c(1.5,1,1,2,2,1,2,0.1,0.1,2,2,3,3)#c(1.5,1,1,2,2,1.5,2,1,1,2,2,2,2)#,rep(5/10,8))
matd=matrix(rep(c(tw.heart1),5),ncol=dim(tw.heart1)[2], byrow=T)
set.seed(156)
matd=matd[sample(nrow(matd)),]
prof=kohonenqualigo.weight(17,nr,nc,0.04,0.01,2.99,0.5,matd,dim(matd)[1],w)
m=kohonenqualiclass.weight(prof,tw.heart1,dim(tw.heart1)[1],w)
htw.clust=cbind(tw.heart,m)[order(m),]
table(htw.clust[,ncol(htw.clust)-1],htw.clust[,ncol(htw.clust)])
nmi(table(htw.clust[,ncol(htw.clust)-1],htw.clust[,ncol(htw.clust)]))
purity(table(htw.clust[,ncol(htw.clust)-1],htw.clust[,ncol(htw.clust)]))
grillecarte(nr,nc,2,htw.clust[,ncol(htw.clust)-1],htw.clust[,ncol(htw.clust)])
par(xpd=TRUE)
legend("topright", inset=c(-0.16,0.2), title="Patient Type", c("No disease","Disease"), pch=15,col=hcl(seq(0,240,length.out=2),120,85),cex=0.55)

#???calibrage (tuning) de paramètres
set.seed(15)
wt<-matrix(runif(13000),ncol=13)
wt=wt/rowSums(wt)
matd=matrix(rep(c(tw.heart1),5),ncol=dim(tw.heart1)[2], byrow=T)
set.seed(156)
matd=matd[sample(nrow(matd)),]
gridsearch=function(nit){
  w=wt[nit,]
  prof=kohonenqualigo.weight(17,nr,nc,0.04,0.01,2.99,0.5,matd,dim(matd)[1],w)
  m=kohonenqualiclass.weight(prof,tw.heart1,dim(tw.heart1)[1],w)
  htw.clust=cbind(tw.heart,m)[order(m),]
  t=nmi(table(htw.clust[,ncol(htw.clust)-1],htw.clust[,ncol(htw.clust)]))
  return(t)
  }
library(doSNOW)
threads =4
cl =makeCluster(threads)
registerDoSNOW(cl)
result = foreach(i=1:nrow(wt)) %dopar% gridsearch(i)
stopCluster(cl)
result=do.call(rbind,result)
write.table(cbind(wt,result), file =paste(paste("tuning cleve1"),".txt",sep=""), sep = " ")
write.table(result, file =paste(paste("tuning result cleve1"),".txt",sep=""), sep = " ")


#GSOM
require(GrowingSOM)
sFr=0.95
tw.heart1=tw.heart[,c(1:13)]
gsom_cleve <- train.gsom(tw.heart1,spreadFactor=sFr, keepdata=TRUE, iterations=60, 
                         alpha=0.5, gridsize = FALSE, nhood= "rect")
m.g=kohonenqualiclass(gsom_cleve$nodes$codes,tw.heart1,dim(tw.heart1)[1])
htw.clust.g=cbind(heart,m.g)[order(m.g),]
table(htw.clust.g[,ncol(htw.clust.g)-1],htw.clust.g[,ncol(htw.clust.g)])
nmi(table(htw.clust.g[,ncol(htw.clust.g)-1],htw.clust.g[,ncol(htw.clust.g)]))
purity(table(htw.clust.g[,ncol(htw.clust.g)-1],htw.clust.g[,ncol(htw.clust.g)]))