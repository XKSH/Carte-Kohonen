base.test=theart
n=303;nr=nc=5;k=46
nquan=c(as.numeric(which(sapply(theart, is.numeric))))
nqual=c(as.numeric(which(sapply(theart, is.factor))))
nmod=unlist(lapply(1:(length(nqual)),function(i){s=length(levels(as.factor(base.test[,nqual[i]]))) ;return(s)}))
mod=sum(nmod)
lmod= unlist(lapply(1:(length(nqual)),function(i){s=levels(as.factor(base.test[,nqual[i]])) ;return(s)}))
index=c(0,cumsum(nmod))
beta=para.expo(3,3,n)
alpha=para.expo(0.4,0.1,n)
profils=matrix(0,nrow=nr*nc,ncol=dim(base.test)[2])
set.seed(seed)
profils[,nquan]=matrix(runif(nquan*nr*nc,0,1),ncol=nquan)
set.seed(seed)
if(length(nqual)>0){
  profils[,nqual]=do.call(rbind,base.test[sample(1:n,nr*nc),nqual])
}else
  {
  profils[,nqual]=base.test[sample(1:n,nr*nc),nqual]
}
profils=as.matrix(profils)
neighboo=vois.constr(nr, nc, beta[1]) 
b=beta[1]
# La matrice de vote pour les variables qualitatives
qual.vote=matrix(0,nr*nc,mod)
if (beta[k]%/%1!=b%/%1) {neighboo=vois.constr(nr, nc, beta[1])} 
dista=apply(profils, 1, function(x, y=base.test) dist.mix(y[1,],x,nquan,nqual,w))
m=which(dista==min(dista))[1]