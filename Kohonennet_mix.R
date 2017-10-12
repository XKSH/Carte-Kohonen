

kohonenqualiclass.mix<-function(profils,dat,n,w,nquan,nqual)
{
  base.test=as.matrix(dat[1:n,])
  m=rep(0,n)
  for (k in 1:n)
  {
    dista=apply(profils, 1, function(x, y=base.test) dist.mix(y[k,],x,nquan,nqual,w))
    m[k]<-which(dista==min(dista))[1]
  }
  return(m) # Vecteur des neurones vainqueurs
}

# Pas de variable de classe dans 'dat'
kohonenqualigo.mix<-function(seed,nr,nc,alpha0,alphaT,beta0,betaT,dat,n,w) 
  # seed, dim carte 1, dim carte 2, taux d'apprentissage initial, final, rayon d'apprentissage initial, final, table des données, nombre d'observations lues, type de pas DTW, contrainte de fenêtre DTW
{
  base.test=dat[1:n,] # x=colonne PCO
  nquan=c(as.numeric(which(sapply(base.test, is.numeric))))
  nqual=c(as.numeric(which(sapply(base.test, is.factor))))
  #mod=Reduce("+", lapply(1:(length(nqual)-1),function(i){s=length(levels(as.factor(base.test[,nqual[i]]))) ;return(s)})) 
    nmod=unlist(lapply(1:(length(nqual)),function(i){s=length(levels(as.factor(base.test[,nqual[i]]))) ;return(s)}))
    mod=sum(nmod)
    lmod= unlist(lapply(1:(length(nqual)),function(i){s=levels(as.factor(base.test[,nqual[i]])) ;return(s)}))
    index=c(0,cumsum(nmod))
  beta=para.expo(beta0,betaT,n)
  alpha=para.expo(alpha0,alphaT,n)
  # Initialisation aléatoire des vecteurs poids
  profils=matrix(0,nrow=nr*nc,ncol=dim(base.test)[2])
  set.seed(seed)
  profils[,nquan]=matrix(runif(length(nquan)*nr*nc,0,1),ncol=length(nquan))
  
  set.seed(seed)
  if(length(nqual)>1){
  profils[,nqual]=do.call(cbind,base.test[sample(1:n,nr*nc,replace = TRUE),nqual])
  }else
  {
    profils[,nqual]=base.test[sample(1:n,nr*nc),nqual]
  }
  profils=as.matrix(profils)
  # Construction de la matrice de voisinage
  neighboo=vois.constr(nr, nc, beta[1]) 
  b=beta[1]
  # La matrice de vote pour les variables qualitatives
  qual.vote=matrix(0,nr*nc,mod)
  for (k in 1:n)
  {
    if (beta[k]%/%1!=b%/%1) {neighboo=vois.constr(nr, nc, beta[k])} 
    # On regarde si la décroissance de béta impacte le nombre de voisins
    b=beta[k]
    dista=apply(profils, 1, function(x, y=base.test) dist.mix(y[k,],x,nquan,nqual,w))
    m=which(dista==min(dista))[1] # Recherche du min de ces distances;si pas unique, on choisit 1er
    nf=exp(-apply(profils, 1, function(x, y=profils) dist.mix(y[m,],x,nquan,nqual,w)/(beta[k]^2)))
    nna=which(!is.na(base.test[k,]))
    n1=intersect(nna,nquan)
    n2=intersect(nna,nqual)
    # Le renouvellement de la partie qualitatiives
      qual.vote[m,]=qual.vote[m,]+c(as.numeric(acm.disjonctif(base.test[k,nqual])))
      candidat=lapply(1:length(nmod),function(i){im=index[i]+which.max(qual.vote[m,(index[i]+1):index[i+1]]);im})
      candidat=unlist(candidat)
      profils[m,nqual]=as.numeric(lmod[candidat])
    
   # Le renouvellement de la partie qualitatiives  
    updt=t(as.vector(do.call(rbind,as.list(base.test[k,n1])))-t(profils[,n1]))*alpha[k]*nf*neighboo[,m]
    updt[is.na(updt)==TRUE]=0
    profils[,n1]=profils[,n1]+updt
  }
  return(profils)
}