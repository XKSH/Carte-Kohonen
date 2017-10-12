# Fonction de voisinnage pour une grille
# nr:nrow,nc:ncol
vois.constr<-function(nr, nc, radius)
{
  # La grille des coordonn�es,x: axe horizontale  y:axe verticale 
  y=rep(1:nr,each=nc)
  x=c(rep(c(1:nc),nr))
  adr=cbind(x,y)
  # La matrice binaire pour la voisinnage  
  neighboo=matrix(rep(0,(nr*nc)^2), nr*nc,nr*nc)
  # moins de loop utilisant la symm�trie
  for (i in 1:(nr*nc))
  {
    for (j in i:(nr*nc))
    {
      if (max(abs(x[i]-x[j]),abs(y[i]-y[j]))<=radius){neighboo[i,j]=1} 
    }
  }
  neighboo=neighboo+t(neighboo)-diag(diag(neighboo))
  return(neighboo)
}
# Enhanced learning
# Fonction de d�croissance exponentielle des param�tres
para.expo<-function(beta0,betaT,n){
  return(beta0*((betaT/beta0)^seq(0,1,len=n)))
  } # BetaT=0.65 ca a l'air cool
