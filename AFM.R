require(FactoMineR)
th=heart
th[,1:4]=scale(th[,1:4])
res = MFA(th, group=c(5,8,1), type=c("s","n","n"), ncp=8, name.group=c("Quantitative", "Qualitative","Malade"),num.group.sup=c(3))
