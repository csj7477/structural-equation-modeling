#Confirmatory Factor Analysis

install.packages("lavaan")
library('lavaan')
f2structure <- read.csv(file="C:/R_DATA/F2structure.csv",head=TRUE)

#Covariance Estimation
mycfa1 <- '
F1 =~ X1+X2+X3; F2 =~ X6+X7+X8
F1 ~~ F1; F2 ~~ F2; F1 ~~ F2
X1 ~~ X1; X2 ~~ X2; X3 ~~ X3 
X6 ~~ X6; X7 ~~ X7; X8 ~~ X8
'
#Intercept Including
mycfa2 <- paste(mycfa1,'X1 ~ 1; X2 ~ 1; X3 ~ 1
                X6 ~ 1; X7 ~ 1; X8 ~ 1',sep='')

obj.mycfa1 <- sem(mycfa1,data=f2structure)
obj.mycfa2 <- sem(mycfa2,data=f2structure)

summary(obj.mycfa1,fit.measures=T,standardized=T)
summary(obj.mycfa2,fit.measures=T,standardized=T)

library('igraph')
library('semPlot')
install.packages('semPlot')
install.packages('igraph')

semPaths(obj.mycfa1,what="par",intercepts=F,
         edge.label.cex=1,edge.color='black',edge.width=0.1)
semPaths(obj.mycfa1,what="std",intercepts=F,
         edge.label.cex=1,edge.color='black',edge.width=0.1)
semPaths(obj.mycfa2,what="par",intercepts=T,
         edge.label.cex=1,edge.color='black',edge.width=0.1)

#CR, AVE
CR.AVE.function <- function(myloadings,mythetas) {
  SS.L2 <- sum(myloadings)^2
  SS.Lsq <- sum(myloadings^2)
  SS.thetas <- sum(mythetas)
  myCR <- SS.L2/(SS.L2+SS.thetas)
  myAVE <- SS.Lsq/(SS.Lsq+SS.thetas)
  myresult <- cbind(myCR,myAVE)
  colnames(myresult) <- c('Composite.Reliability',
                          'Ave.variance.extracted')
  myresult 
}

myest <- standardizedSolution(obj.mycfa1)
myest

F1.loadings <- myest[1:3,4]
F2.loadings <- myest[4:6,4]
F1.var.theta <- myest[10:12,4]
F2.var.theta <- myest[13:15,4]

CR.AVE.function(F1.loadings,F1.var.theta)
CR.AVE.function(F2.loadings,F2.var.theta)

#MTMM
mtmm <- read.csv(file="C:/R_DATA/mtmm.csv",head=TRUE)
colnames(mtmm)
round(cor(mtmm),2)

myefa.f3 <- factanal(mtmm,factors=3,rotation='promax')
myefa.f3

mymtmm <- '
t1 =~ t1m1+t1m2+t1m3
t2 =~ t2m1+t2m2+t2m3
t3 =~ t3m1+t3m2+t3m3
t1 ~~ t1; t2 ~~ t2; t3 ~~ t3 
t1 ~~ t2; t1 ~~ t3; t2 ~~ t3 
t1m1 ~~ t1m1;t1m2 ~~ t1m2;t1m3 ~~ t1m3
t2m1 ~~ t2m1;t2m2 ~~ t2m2;t2m3 ~~ t2m3
t3m1 ~~ t3m1;t3m2 ~~ t3m2;t3m3 ~~ t3m3
t1m1 ~~ t2m1;t1m1 ~~ t3m1;t2m1 ~~ t3m1
t1m2 ~~ t2m2;t1m2 ~~ t3m2;t2m2 ~~ t3m2
t1m3 ~~ t2m3;t1m3 ~~ t3m3;t2m3 ~~ t3m3
'
obj.mtmm <- sem(mymtmm,data=mtmm)
summary(obj.mtmm,fit.measures=T,standardized=T)

#semPaths(obj.mtmm,what="std",intercepts=F,
#edge.label.cex=0.7,edge.color='black',edge.width=0.1)

#Equality Constraints
mymtmm.ec <- '
t1 =~ lmbd_1*t1m1+lmbd_1*t1m2+lmbd_1*t1m3
t2 =~ t2m1+t2m2+t2m3
t3 =~ t3m1+t3m2+t3m3
t1 ~~ t1; t2 ~~ t2; t3 ~~ t3 
t1 ~~ t2; t1 ~~ t3; t2 ~~ t3 
t1m1 ~~ t1m1;t1m2 ~~ t1m2;t1m3 ~~ t1m3
t2m1 ~~ t2m1;t2m2 ~~ t2m2;t2m3 ~~ t2m3
t3m1 ~~ t3m1;t3m2 ~~ t3m2;t3m3 ~~ t3m3
t1m1 ~~ t2m1;t1m1 ~~ t3m1;t2m1 ~~ t3m1
t1m2 ~~ t2m2;t1m2 ~~ t3m2;t2m2 ~~ t3m2
t1m3 ~~ t2m3;t1m3 ~~ t3m3;t2m3 ~~ t3m3
'
obj.mtmm.ec <- sem(mymtmm.ec,data=mtmm)
summary(obj.mtmm.ec,fit.measures=T,standardized=T)

anova(obj.mtmm, obj.mtmm.ec)

mymtmm.ec2 <- '
t1 =~ t1m1+t1m2+t1m3
t2 =~ t2m1+t2m2+t2m3
t3 =~ t3m1+t3m2+t3m3
t1 ~~ t1; t2 ~~ t2; t3 ~~ t3 
t1 ~~ t2; t1 ~~ t3; t2 ~~ t3 
t1m1 ~~ t1m1;t1m2 ~~ t1m2;t1m3 ~~ t1m3
t2m1 ~~ t2m1;t2m2 ~~ t2m2;t2m3 ~~ t2m3
t3m1 ~~ t3m1;t3m2 ~~ t3m2;t3m3 ~~ t3m3
t1m1 ~~ d1*t2m1;t1m1 ~~ d1*t3m1;t2m1 ~~ d1*t3m1
t1m2 ~~ d2*t2m2;t1m2 ~~ d2*t3m2;t2m2 ~~ d2*t3m2
t1m3 ~~ d3*t2m3;t1m3 ~~ d3*t3m3;t2m3 ~~ d3*t3m3
'
obj.mtmm.ec2 <- sem(mymtmm.ec2,data=mtmm)
summary(obj.mtmm.ec2,fit.measures=T,standardized=T)

anova(obj.mtmm, obj.mtmm.ec2)

#Modification Index
myMIexample <- '
F1 =~ X1+X2+X3+X6; F2 =~ X7+X8
F1 ~~ F1; F2 ~~ F2; F1 ~~ F2
X1 ~~ X1; X2 ~~ X2; X3 ~~ X3 
X6 ~~ X6; X7 ~~ X7; X8 ~~ X8 
'
obj.myMIexample <- sem(myMIexample,data=f2structure)
summary(obj.myMIexample,fit.measures=T,standardized=T)
myMI <- modindices(obj.myMIexample)
myMI
myMI[myMI$op=='=~',]


