require(GrowingSOM)
sFr=0.96
theart1=tw.heart1
matd=matrix(rep(theart1,2),ncol=dim(theart1)[2], byrow=T)
gsom_cleve <- train.gsom(matd,spreadFactor=sFr, keepdata=TRUE, iterations=30, 
                        alpha=0.5, gridsize = FALSE, nhood= "rect")
m.g=kohonenqualiclass(gsom_cleve$nodes$codes,theart1,dim(theart1)[1])
htclust.g=cbind(heart,m.g)[order(m.g),]
table(htclust.g[,ncol(htclust.g)-1],htclust.g[,ncol(htclust.g)])
nmi(table(htclust.g[,ncol(htclust.g)-1],htclust.g[,ncol(htclust.g)]))
purity(table(htclust.g[,ncol(htclust.g)-1],htclust.g[,ncol(htclust.g)]))
# Plot
xlim=c(min(gsom_cleve$nodes$position$x),max(gsom_cleve$nodes$position$x))
ylim=c(min(gsom_cleve$nodes$position$y),max(gsom_cleve$nodes$position$y))
adrg=expand.grid(xlim[1]:xlim[2],ylim[1]:ylim[2])

par(bg = "white",mai=c(0.2,0.2,0.6,0.8))
plot(c(xlim[1]-1,xlim[2]), c(ylim[1]-1, ylim[2]), xlab = "", ylab = "",axes=FALSE,xaxs="i",yaxs="i",main=paste("Carte Kohonen GSOM",sFr))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="white")
# Construire la palette de trois type
ncol=seq(0,240,length.out=2);
# Faire varier les couleurs en utilisant hcl
vcol=ncol%*%(prop.table(table(htclust.g[,14],htclust.g[,15]), 2))
vnames=c(as.numeric(colnames(vcol)))
rvcol=hcl(vcol,120,85)
# Affecter les couleurs aux noeuds
co=cbind(gsom_cleve$nodes$position$x,gsom_cleve$nodes$position$y)

for(i in 1:nrow(co)){
  rect(co[i,1]-1,co[i,2]-1,co[i,1],co[i,2],col=hcl(vcol[i],120,85),lwd=0)
}
grid(nx = xlim[2]+1-xlim[1],ny=ylim[2]+1-ylim[1])
# Ajouter le legend dans la marge
par(xpd=TRUE)
legend("topright", inset=c(-0.16,0.2), title="Patient Type", c("No disease","Disease"), pch=15,col=hcl(seq(0,240,length.out=2),120,85),cex=0.55)
prop.table(table(htclust.g[,14],htclust.g[,15]), 2)