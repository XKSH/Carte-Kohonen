heart.data =data.frame(read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data",header=FALSE,sep=",",na.strings = '?'))
chclass =c("numeric","factor","factor","numeric","numeric","factor","factor","numeric","factor","numeric","factor","factor","factor","factor")
heart.data[,14][heart.data[,14] > 0] <- 1
# La conversion de variable
convert.magic = function(obj, types){ 
  out = lapply(1:dim(obj)[2],FUN = function(i){FUN1 <- switch(types[i],character = as.character,numeric = as.numeric,factor = as.factor); FUN1(obj[,i])})  
  names(out) = colnames(obj)
  data.frame(out,stringsAsFactors = FALSE)
} 
heart.data =convert.magic(heart.data,chclass)
names(heart.data) = c( "age", "sex", "cp", "trestbps", "chol","fbs", "restecg",
                        "thalach","exang", "oldpeak","slope", "ca", "thal", "num")

theart=heart.data[,c(1:8)]
nquan=c(as.numeric(which(sapply(theart, is.numeric))))
nqual=c(as.numeric(which(sapply(theart, is.factor))))

theart[,nquan]=scale(theart[,nquan])

theart=do.call("rbind", replicate(1, theart, simplify = FALSE))
nr=8;nc=8
w=rep(1,8)
#w=c(1.5,1,1,2,2,1,2,0.1)#,0.1,2,2,3,3)
# Il faut deux variables qualitatives à cause de la limite de acm.disjoint
prof=kohonenqualigo.mix(17,nr,nc,0.4,0.1,3,3,theart[sample(nrow(theart)),],dim(theart)[1],w)
m=kohonenqualiclass(prof,theart,dim(theart)[1])
htclust=cbind(theart,m)[order(m),]
htclust=cbind(htclust,heart.data[14])
table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1])
nmi(table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1]))
purity(table(htclust[,ncol(htclust)],htclust[,ncol(htclust)-1]))

