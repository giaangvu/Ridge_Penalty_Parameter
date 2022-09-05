install.packages("R.matlab")
library(R.matlab)
library(MASS)
library(car)
library(clipr)

setwd("/Users/giangvu/Documents/Pitt/FALL20/STAT 2270 - Data Mining/Project")

data(Boston)

datpr <- Boston
set.seed(22700)
datpr_sp <- datpr[sample(nrow(datpr), 100), ]  
boxplot(datpr_sp$medv)
write_clip(datpr_sp)

cooks.distance(lm(medv~.,data = datpr_sp[,-c(2,4,9,10)]))
plot(cooks.distance(lm(medv~.,data = datpr_sp[,-c(2,4,9,10)])))
which(cooks.distance(lm(medv~.,data = datpr_sp[,-c(2,4,9,10)]))>(4/(100-9-1))) #4/(n-p-1) 
#cook's distance identified  points 11  18  22  30  40  76 to be outliers  
#(see how many of them are inlfuential using authors' method)


# Read covariates and outcome
prX <- datpr_sp[,-c(2,4,9,10,14)] 
pry <- datpr_sp$medv 
pr_p <- dim(prX)[2]
pr_n <- dim(prX)[1]

prX <- scale(prX,center=F) # Scale data
set.seed(2270)
prX = cbind(rep(1,pr_n),prX) # Add intercept to design matrix
pr_p=pr_p+1 # New dimension after intercept is added

#Do PCA.
pr_pca <- prcomp(prX)
pr_svd <- svd(prX)

# Percentage of variance explained by the 1st principal component: 0.63175
summary(pr_pca)$importance[2,1]

#Define explicit expression for leave-one-out CV error
tuning.cv.svd <- function(lambda,w){
  H <- pr_svd$u  %*%  diag(pr_svd$d^2/(pr_svd$d^2 +lambda)) %*% t(pr_svd$u)
  e <- (diag(pr_n) - H) %*% pry
  mean(w*(e/(1-diag(H)))^2)
}

#Visualize the 1. principal component and outcome
postscript(file='InfluentialBostonHousingScatterplotPCA.eps', width = 8, height = 8,
           horizontal = FALSE, onefile = FALSE)
plot(pr_pca$x[,1],pry,xlab='PC1',ylab='y',cex.lab=1,mgp=c(2.8,1,0),cex.axis=1)
points(pr_pca$x[11,1],pry[11],col='red',pch=16,cex=1.2)
points(pr_pca$x[8,1],pry[8],col='purple',pch=16,cex=1.2)
points(pr_pca$x[22,1],pry[22],col='green',pch=16,cex=1.2)
points(pr_pca$x[18,1],pry[18],col='orange',pch=16,cex=1.2)
points(pr_pca$x[33,1],pry[33],col='blue',pch=16,cex=1.2)

text(pr_pca$x[11,1],pry[11]-1,col='red',labels='11',cex=1)
text(pr_pca$x[8,1],pry[8]+1,col='purple',labels='8',cex=1)
text(pr_pca$x[22,1],pry[22]+1,col='green',labels='22',cex=1)
text(pr_pca$x[18,1],pry[18]+1,col='orange',labels='18',cex=1)
text(pr_pca$x[33,1],pry[33]-1,col='blue',labels='33',cex=1)
linreg=lm(y~x, data.frame(y=pry,x=pr_pca$x[,1]))
abline(h=mean(pry),lty=3)
abline(coef(linreg)[1],coef(linreg)[2],lty=2)

dev.off()

#8 11 18 22 33
# 8 (purple) 11 (red) 18 (orange) 22 (green) 33 (blue)

#Define grid of weight values, for each observation
nw=100 # no of grid points
weights <- seq(0,4.1,length.out = nw) # weight grid

# Find optimal cv-lambda for each observation and weight value
pr_lambdaMatrix <- matrix(,pr_n,nw)
for(i in 1:nw){
  for(j in 1:pr_n){
    w <- rep(1,pr_n)
    w[j] <- weights[i] # Weight to give observation i
    # find optimal lambda
    pr_lambdaMatrix[j,i]<- optim(par = 0.01,tuning.cv.svd,lower=-Inf,upper = Inf,
                              w=w/sum(w),method = "L-BFGS-B",control = list(factr=1e-3))$par
  }
}

# 8 (purple) 11 (red) 18 (orange) 22 (green) 33 (blue)

# Plot cv-lambda against weight
pr_lambdaMatrix[which(pr_lambdaMatrix<0,arr.ind = T)]=0 # Ensure no lambdas are negative
postscript(file='InfluentialBostonHousingLambdaCurve.eps', width = 8, height = 8,
           horizontal = FALSE, onefile = FALSE)
matplot(weights,t(pr_lambdaMatrix),type='l',ylim = c(-5e-5,max(pr_lambdaMatrix)+0.001),xlim=c(0,max(weights)+0.4),
        xlab='Weight of observation',ylab='Tuning parameter',
        col=c(rep('grey',7),'purple',rep('grey',2),'red',rep('grey',6),'orange',rep('grey',3),'green',rep('grey',10),'blue',rep('grey',67)),
        xaxs="i",yaxs='i',
        lty=1,lwd=c(rep(1,7),2,rep(1,2),2,rep(1,6),2,rep(1,3),2,rep(1,10),2,rep(1,67)),cex.lab=1,mgp=c(2.8,1,0),cex.axis=1)
abline(v=1) # Mark out weight=1
text(4.3,y=pr_lambdaMatrix[11,nw]-0.1,labels= as.character(11),col='red',cex=1) # Mark out some notable observations
text(4.3,y=pr_lambdaMatrix[18,nw],labels= as.character(18),col='orange',cex=1)
text(4.3,y=pr_lambdaMatrix[22,nw]-0.04,labels= as.character(22),col='green',cex=1)
text(4.3,y=pr_lambdaMatrix[33,nw]-0.01,labels= as.character(33),col='blue',cex=1)
text(4.3,y=pr_lambdaMatrix[8,nw]+0.02,labels= as.character(8),col='purple',cex=1)

dev.off()

# 8 (purple) 11 (red) 18 (orange) 22 (green) 33 (blue)
# Plot the effective degrees of freedom
df.bf = apply(pr_lambdaMatrix, c(1,2), function(lam) sum(pr_svd$d^2/(pr_svd$d^2 +lam))) # Formula for the effective degrees of freedom
postscript(file='InfluentialBostonDFCurve.eps', width = 8, height = 8,
           horizontal = FALSE, onefile = FALSE)
matplot(weights,t(df.bf),type='l',ylim = c(min(df.bf),max(df.bf)),xlim=c(0,max(weights)+0.3),
        xlab='Weight of observation',ylab='Effective degrees of freedom',
        col=c(rep('grey',7),'purple',rep('grey',2),'red',rep('grey',6),'orange',rep('grey',3),'green',rep('grey',10),'blue',rep('grey',67)),
        lty=1,
        lwd=c(rep(1,7),2,rep(1,2),2,rep(1,6),2,rep(1,3),2,rep(1,10),2,rep(1,67)),xaxs='i',yaxs='i',
        cex.lab=1,mgp=c(2.8,1,0),cex.axis=1)
abline(v=1)
text(4.2,y=df.bf[11,nw]+0.05,labels= as.character(11),col='red',cex=1)
text(4.2,y=df.bf[18,nw],labels= as.character(18),col='orange',cex=1)
text(4.2,y=df.bf[22,nw]+0.05,labels= as.character(22),col='green',cex=1)
text(4.2,y=df.bf[33,nw]+0.01,labels= as.character(33),col='blue',cex=1)
text(4.2,y=df.bf[8,nw]-0.01,labels= as.character(8),col='purple',cex=1)

dev.off()








