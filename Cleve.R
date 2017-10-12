library(FactoClass)
# Introduire les donées:
heart.data <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data",header=FALSE,sep=",",na.strings = '?')
chclass <-c("numeric","factor","factor","numeric","numeric","factor","factor","numeric","factor","numeric","factor","factor","factor","factor")
heart.data[,14][heart.data[,14] > 0] <- 1
# La conversion de variable
convert.magic <- function(obj, types){ 
   out <- lapply(1:length(obj),FUN = function(i){FUN1 <- switch(types[i],character = as.character,numeric = as.numeric,factor = as.factor); FUN1(obj[,i])})  
   names(out) <- colnames(obj)
   data.frame(out,stringsAsFactors = FALSE)
 } 
heart.data <- convert.magic(heart.data,chclass)
names(heart.data) <- c( "age", "sex", "cp", "trestbps", "chol","fbs", "restecg",
                        "thalach","exang", "oldpeak","slope", "ca", "thal", "num")
heart=heart.data[c(1,4,5,8,10,2:3,6:7,9,11:14)]
names(heart)=names(heart.data)[c(1,4,5,8,10,2:3,6:7,9,11:14)]
heart=heart[complete.cases(heart),]
theart=scale(heart[,1:5])
#theart=t(abs(t(heart[,3:4])-apply(heart[,3:4],2,min)))
#theart=t(apply(theart,1,"/",c(apply(heart[,3:4],2,max)-apply(heart[,3:4],2,min))))
# Les tableaux
disj=as.matrix(acm.disjonctif(heart[6:13]))
burt=t(disj)%*%disj
# Tableaux corrigés
tab.cor=function(tab){
  x=matrix(sqrt(rowSums(tab)),ncol=1)%*%sqrt(colSums(tab))
  tabcor=tab/x
  return(tabcor)
}
#disj.cor=tab.cor(disj)
disj=sqrt(0.5)*disj
#dis.pkq=0.7*disj/c(sqrt(colSums(disj)/dim(disj)[1]))
theart=cbind(theart,disj)
source("Kohonenprep.R")
source("Kohonennet.R")
source("plotcarte1.R")
source("evaluation.R")
matd=matrix(rep(theart,10),ncol=dim(theart)[2], byrow=T)
nr=8;nc=8
prof=kohonenqualigo(17,nr,nc,0.4,0.1,3,3,matd[sample(nrow(matd)),],dim(matd)[1])
m=kohonenqualiclass(prof,theart,dim(theart)[1])
htclust=cbind(theart,m)[order(m),]
htclust=cbind(htclust,heart[14])
table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1])
nmi(table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1]))
purity(table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1]))

grillecarte(nr,nc,2,htclust[,30],htclust[,29])
par(xpd=TRUE)
legend("topright", inset=c(-0.16,0.2), title="Patient Type", c("No disease","Disease"), pch=15,col=hcl(seq(0,240,length.out=2),120,85),cex=0.55)

# Feature selection; pas d'effet
ft=data.frame(disj,heart[,14])
fi=cnmi(ft)
tdisj.cor=disj.cor[,which(fi>0.1)]