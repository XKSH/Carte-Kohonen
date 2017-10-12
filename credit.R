test=read.csv("D:/Profiles/XU/Desktop/quan&qual/test.csv",header=TRUE)
nr=nc=8
matd=test
prof=kohonenqualigo(17,nr,nc,0.4,0.1,3,3,matd[sample(nrow(matd)),],dim(matd)[1])