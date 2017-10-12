# Comparer les résulmtat de GSOM et SOM en utilisant le package GrowingSOM
# Entraîner le Modèle GSOM
sFr=0.1
require(GrowingSOM)
gsom_iris <- train.gsom(iris[,1:4],spreadFactor=sFr, keepdata=TRUE, iterations=30, 
                        alpha=0.5, gridsize = FALSE, nhood= "rect")
m.g=kohonenqualiclass(gsom_iris$nodes$codes,iris[,1:4],dim(iris)[1])
irclust.g=cbind(iris,m.g)[order(m.g),]
table(irclust.g[,5],irclust.g[,6])
nmi(table(irclust.g[,5],irclust.g[,6]))
purity(table(irclust.g[,5],irclust.g[,6]))
#Plot
xlim=c(min(gsom_iris$nodes$position$x),max(gsom_iris$nodes$position$x))
ylim=c(min(gsom_iris$nodes$position$y),max(gsom_iris$nodes$position$y))
adrg=expand.grid(xlim[1]:xlim[2],ylim[1]:ylim[2])
par(bg = "white",mai=c(0.2,0.2,0.6,0.8))
plot(c(xlim[1]-1,xlim[2]), c(ylim[1]-1, ylim[2]), xlab = "", ylab = "",axes=FALSE,xaxs="i",yaxs="i",main=paste("Carte Kohonen GSOM",sFr))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="white")
# Construire la palette de trois type
ncol=seq(0,240,length.out=3);
# Faire varier les couleurs en utilisant hcl
vcol=ncol%*%(prop.table(table(irclust.g[,5],irclust.g[,6]), 2))
vnames=c(as.numeric(colnames(vcol)))
rvcol=hcl(vcol,120,85)
# Affecter les couleurs aux noeuds
co=cbind(gsom_iris$nodes$position$x,gsom_iris$nodes$position$y)

for(i in 1:nrow(co)){
    rect(co[i,1]-1,co[i,2]-1,co[i,1],co[i,2],col=hcl(vcol[i],120,85),lwd=0)
}
grid(nx = xlim[2]+1-xlim[1],ny=ylim[2]+1-ylim[1])
# Ajouter le legend dans la marge
par(xpd=TRUE)
legend("topright", inset=c(-0.16,0.2), title="Iris Type", c("setosa","versicolo","virginica"), pch=15,col=hcl(ncol,120,85),cex=0.55)
