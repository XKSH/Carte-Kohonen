theart=scale(heart[,1:5])
disj=as.matrix(acm.disjonctif(heart[,6:13]))
burt=t(disj)%*%disj
# Tableaux corrigés
tab.cor=function(tab){
  x=matrix(sqrt(rowSums(tab)),ncol=1)%*%sqrt(colSums(tab))
  tabcor=tab/x
  return(tabcor)
}
disj.cor=tab.cor(disj)
theart=cbind(theart,disj)
#w=c(rep(1,5),rep(1/10000/c(do.call(rbind,lapply(heart[,6:13],nlevels))),c(do.call(rbind,lapply(heart[,6:13],nlevels)))))

matd=matrix(rep(theart,10),ncol=dim(theart)[2], byrow=T)
nr=8;nc=8
prof=kohonenqualigo.weight(17,nr,nc,0.04,0.01,2.99,0.65,matd[sample(nrow(matd)),],dim(matd)[1],w)
m=kohonenqualiclass.weight(prof,theart,dim(theart)[1])
htclust=cbind(theart,m)[order(m),]
htclust=cbind(htclust,heart[14])
table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1])
nmi(table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1]))
purity(table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1]))
# Le calcul de chi carré statistique pour les variables
# Convert les variables quantitatives
wheart=heart
wheart[,1:5]=apply(wheart[,1:5],2,function(x){cut(x,3)})
wheart[,1:5]=convert.magic(wheart[,1:5],c(rep(c("factor"),5)))
vchis=rep(0,13)
for(i in 1:13){
  wburt=as.matrix(acm.burt(data.frame(wheart[,i]),data.frame(wheart[,14])))
  e=rowSums(wburt)%*%t(colSums(wburt))/sum(wburt)
  vchis[i]=sum((wburt-e)^2/e)
}
w=vchis/sum(vchis)
#pour KACM
w=w/c(rep(1,5),rep(2,8))
w=rep(w,c(1,1,1,1,1,do.call(rbind,lapply(heart[,6:13],nlevels))))
#Transformer variables numériques en variables catégoriques
w=rep(w,c(rep(3,5),do.call(rbind,lapply(heart[,6:13],nlevels))))
disj=as.matrix(acm.disjonctif(wheart[,1:13]))
matd=matrix(rep(disj,10),ncol=dim(disj)[2], byrow=T)
nr=8;nc=8
prof=kohonenqualigo.weight(17,nr,nc,0.04,0.01,2.99,0.65,matd[sample(nrow(matd)),],dim(matd)[1],w)
m=kohonenqualiclass.weight(prof,disj,dim(disj)[1],w)
htclust=cbind(theart,m)[order(m),]
htclust=cbind(htclust,heart[14])
table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1])
nmi(table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1]))
purity(table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1]))