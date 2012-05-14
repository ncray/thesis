library(kernlab)
library(multicore)
x <- c(rnorm(100,mean=-.8),rnorm(100))
x <- c(rnorm(100,mean=0),rnorm(100))
x <- c(rpois(100,1),rpois(100,3))
y <- c(rep(-1,100),rep(1,100))

y.perm <- sample(y)
t.test(x[y.perm==1], x[y.perm==-1], var.equal=TRUE)
t.test(x[y.perm==1], x[y.perm==-1])
sqrt(var(x[y.perm==1])/100 + var(x[y.perm==-1])/100)
##(mean(x[y.perm==1]) - mean(x[y.perm==-1]))/sqrt(var(x[y.perm==1])/100 + var(x[y.perm==-1])/100)

y.perm <- sample(y)
sample(length(y.perm), 2)

x <- 1:10
swap <- sample(length(x), 2)
x[swap] = x[rev(swap)]

len <- 100000
T <- 1:len
diff <- 1:len
for(i in 1:len){
  y.perm <- sample(y)
  T[i] <- t.test(x[y.perm==1], x[y.perm==-1], var.equal=TRUE)$statistic
  swap <- sample(length(y.perm),2)
  y.perm[swap] <- y.perm[rev(swap)]
  diff[i] <- t.test(x[y.perm==1], x[y.perm==-1], var.equal=TRUE)$statistic - T[i]
}
plot(T, diff)
summary(lm(diff~T))

samp.size <- 40
x <- c(rnorm(samp.size,mean=1),rnorm(samp.size))
y <- c(rep(-1,samp.size),rep(1,samp.size))
num.perm <- 50
T <- 1:num.perm
diff <- 1:num.perm
for(i in 1:num.perm){
  y.perm <- sample(y)
  T[i] <- t.test(x[y.perm==1], x[y.perm==-1], var.equal=TRUE)$statistic
  ##T[i] <- mean(x[y.perm==1]) - mean(x[y.perm==-1])
  minus <- which(y.perm==-1)
  plus <- which(y.perm==1)
  diff.each <- 1:(samp.size^2)
  for(j in 1:samp.size){
    for(k in 1:samp.size){
      swap <- c(minus[j], plus[k])
      y.perm[swap] <- y.perm[rev(swap)]
      diff.each[samp.size*(j-1)+k] <- t.test(x[y.perm==1], x[y.perm==-1], var.equal=TRUE)$statistic - T[i]
      ##diff.each[samp.size*(j-1)+k] <- mean(x[y.perm==1]) - mean(x[y.perm==-1]) - T[i]
      ##swap it back
      y.perm[swap] <- y.perm[rev(swap)]
    }
  }
  diff[i] <- mean(diff.each)
}
plot(T, diff)
summary(lm(diff~T))


sk <- vanilladot()
km <- kernelMatrix(sk,as.matrix(x))

ksvm1 <- ksvm(km,y)
t.stat <- t.test(fitted(ksvm1)[y==1],fitted(ksvm1)[y==-1])$statistic

FSVMP <- function(i, km, y, ret.str="t.perm"){
  y.perm <- sample(y)
  ksvm2 <- ksvm(km,y.perm)
  t.perm <- t.test(fitted(ksvm2)[y.perm==1],fitted(ksvm2)[y.perm==-1], var.equal=TRUE)$statistic
  ##t.perm <- t.test(fitted(ksvm2)[y.perm==1],fitted(ksvm2)[y.perm==-1])$statistic
  t.perm2 <- t.test(x[y.perm==1],x[y.perm==-1], var.equal=TRUE)$statistic
  fits <- fitted(ksvm2)
  list("t.perm" = t.perm, "t.perm2" = t.perm2, "fits" = fits)[ret.str]
}

t.perm <- unlist(mclapply(1:1000, FSVMP, km, y, "t.perm", mc.cores = 4))
plot(sort(t.perm), sort(abs(rnorm(1000))));abline(a=0,b=1)

library(VGAM)
plot(qfnorm(seq(1/1001, 1000/1001, length.out = 1000)),sort(t.perm));abline(a=0,b=1)

checkunif <- function(i){
  x <- c(rnorm(100,mean=0),rnorm(100))
  y <- c(rep(-1,100),rep(1,100))
  sk <- vanilladot()
  km <- kernelMatrix(sk,as.matrix(x))
  t.perm <- unlist(mclapply(1:1000, FSVMP, km, y, "t.perm", mc.cores = 4))
  ks.test(t.perm, abs(rnorm(1000)))$p.value
}

p.vals <- unlist(lapply(1:1000, checkunif))








library(kernlab)
x <- c(rnorm(100,mean=-.8),rnorm(100))
x <- c(rnorm(100,mean=0),rnorm(100))
y <- c(rep(-1,100),rep(1,100))
t.test(x[y==-1],x[y==1])
##this t-stat is the same as the one on the fitted values

sk <- vanilladot()
km <- kernelMatrix(sk,as.matrix(x))

##ksvm1 <- ksvm(km,y,type="C-svc")
ksvm1 <- ksvm(km,y)
t.stat <- t.test(fitted(ksvm1)[y==1],fitted(ksvm1)[y==-1])$statistic
##because the data have been linearly transformed
plot(x,fitted(ksvm1))
t.perm <- rep(0,1000)
t.perm2 <- rep(0,1000)
fits <- matrix(0,nrow=1000,ncol=200)
for(i in 1:1000){
  y.perm <- sample(y)
  ksvm2 <- ksvm(km,y.perm)
  t.perm[i] <- t.test(fitted(ksvm2)[y.perm==1],fitted(ksvm2)[y.perm==-1], var.equal=TRUE)$statistic
  ##t.perm[i] <- t.test(fitted(ksvm2)[y.perm==1],fitted(ksvm2)[y.perm==-1])$statistic
  ##t.perm[i] <- sample(c(-1,1),1)*t.test(fitted(ksvm2)[y.perm==1],fitted(ksvm2)[y.perm==-1])$statistic
  t.perm2[i] <- t.test(x[y.perm==1],x[y.perm==-1], var.equal=TRUE)$statistic
  fits[i,] <- fitted(ksvm2)
}
summary(t.perm)

FSVMP <- function(i, km, y, ret.str="t.perm"){
  y.perm <- sample(y)
  ksvm2 <- ksvm(km,y.perm)
  t.perm <- t.test(fitted(ksvm2)[y.perm==1],fitted(ksvm2)[y.perm==-1], var.equal=TRUE)$statistic
  ##t.perm <- t.test(fitted(ksvm2)[y.perm==1],fitted(ksvm2)[y.perm==-1])$statistic
  ##t.perm <- sample(c(-1,1),1)*t.test(fitted(ksvm2)[y.perm==1],fitted(ksvm2)[y.perm==-1])$statistic
  t.perm2 <- t.test(x[y.perm==1],x[y.perm==-1], var.equal=TRUE)$statistic
  fits <- fitted(ksvm2)
  list("t.perm" = t.perm, "t.perm2" = t.perm2, "fits" = fits)[ret.str]
}

lapply(1:2, FSVMP, km, y)
library(multicore)
t.perm <- unlist(mclapply(1:1000, FSVMP, km, y, "t.perm", mc.cores = 4))

plot(fitted(ksvm1),ylim=c(-3,5))
points(x,col="red")
qqnorm(t.perm);qqline(t.perm)
qqnorm(t.perm2);qqline(t.perm2)
plot(t.perm2,t.perm)
qqnorm(abs(t.perm));qqline(abs(t.perm))

plot(sort(t.perm), qt(seq(.5,.999,length=1000),df=198))

plot(sort(t.perm), sort(abs(rnorm(1000))));abline(a=0,b=1)


library(VGAM)
detach("package:VGAM")
qqplot(qfnorm(seq(1/1000,1-1/1000,length.out=1000)),-t.perm)
fit=vglm(-t.perm~1,fam=fnormal1,trace=TRUE)
coefs <- Coef(fit)
#       mu        sd 
#0.4142370 0.9621249
#0.5320427 0.8391376
#0.5003790 0.9128801 
plot(seq(0,5,.01),dfnorm(seq(0,5,.01),mean=coefs[1],sd=coefs[2]),type="l")
lines(density(-t.perm,bw=.2),col="red")
lines(seq(0,5,.01),dfnorm(seq(0,5,.01)),type="l",col="blue")
pdf("p1.pdf")
hist(-t.perm,prob=TRUE,breaks=10,main="black:fitted;blue:0/1;red:density")
lines(seq(0,5,.01),dfnorm(seq(0,5,.01),mean=coefs[1],sd=coefs[2]),type="l")
lines(seq(0,5,.01),dfnorm(seq(0,5,.01)),type="l",col="blue")
lines(density(-t.perm,bw=.2),col="red")
dev.off()
plot(seq(0,5,.01),dfnorm(seq(0,5,.01),mean=sqrt(2/pi),sd=sqrt(1-2/pi)),type="l")
lines(density(-t.perm,bw=.1),col="red")
hist(-t.perm,prob=TRUE,breaks=10)
lines(seq(0,5,.01),dfnorm(seq(0,5,.01),mean=sqrt(2/pi),sd=sqrt(1-2/pi)),type="l")




##try classification
x <- c(rnorm(50,mean=-.8),rnorm(50))
y <- as.factor(c(rep(-1,50),rep(1,50)))
t.test(x[y==-1],x[y==1])
sk <- vanilladot()
km <- kernelMatrix(sk,as.matrix(x))
ksvm1 <- ksvm(km,y,prob.model=TRUE)
preds <- predict(ksvm1,km,type="probabilities")[,1]
t.stat <- t.test(preds[y==-1],preds[y==1])$statistic
t.perm <- rep(0,1000)
t.perm2 <- rep(0,1000)
fits <- matrix(0,nrow=1000,ncol=100)
for(i in 1:1000){
  y.perm <- sample(y)
  ksvm2 <- ksvm(km,y.perm,prob.model=TRUE)
  preds <- predict(ksvm2,km,type="probabilities")[,1]
  t.perm[i] <- t.test(preds[y.perm==-1],preds[y.perm==1])$statistic
  t.perm2[i] <- t.test(x[y.perm==-1],x[y.perm==1])$statistic
  fits[i,] <- preds
  print(i)
}



library(VGAM)
detach("package:VGAM")
qqplot(qfnorm(seq(1/1000,1-1/1000,length.out=1000)),-t.perm)
fit=vglm(-t.perm~1,fam=fnormal1,trace=TRUE)
coefs <- Coef(fit)
#       mu        sd 
#0.4142370 0.9621249
#0.5320427 0.8391376
#0.5003790 0.9128801 
plot(seq(0,5,.01),dfnorm(seq(0,5,.01),mean=coefs[1],sd=coefs[2]),type="l")
lines(density(-t.perm,bw=.2),col="red")
lines(seq(0,5,.01),dfnorm(seq(0,5,.01)),type="l",col="blue")
pdf("p1.pdf")
hist(-t.perm,prob=TRUE,breaks=10,main="black:fitted;blue:0/1;red:density")
lines(seq(0,5,.01),dfnorm(seq(0,5,.01),mean=coefs[1],sd=coefs[2]),type="l")
lines(seq(0,5,.01),dfnorm(seq(0,5,.01)),type="l",col="blue")
lines(density(-t.perm,bw=.2),col="red")
dev.off()
plot(seq(0,5,.01),dfnorm(seq(0,5,.01),mean=sqrt(2/pi),sd=sqrt(1-2/pi)),type="l")
lines(density(-t.perm,bw=.1),col="red")
hist(-t.perm,prob=TRUE,breaks=10)
lines(seq(0,5,.01),dfnorm(seq(0,5,.01),mean=sqrt(2/pi),sd=sqrt(1-2/pi)),type="l")

##mean with normal mean 0, sd 1
sqrt(2/pi)
##sd with normal mean 0, sd 1
sqrt(1-2/pi)
mean(rfnorm(1000))
sd(rfnorm(1000))

qqplot(-t.perm,qfnorm(seq(1/1000,1-1/1000,length.out=1000)));abline(a=0,b=1)
pdf("p2.pdf")
qqplot(-t.perm,qfnorm(seq(1/1000,1-1/1000,length.out=1000),coefs[1],coefs[2]));abline(a=0,b=1)
dev.off()
qqplot(abs(qnorm(seq(1/1000,1-1/1000,length.out=1000))),qfnorm(seq(1/1000,1-1/1000,length.out=1000)));abline(a=0,b=1)
qqplot(abs(qnorm(seq(1/1000,1-1/1000,length.out=1000))),qfnorm(seq(1/1000,1-1/1000,length.out=1000),coefs[1],coefs[2]));abline(a=0,b=1)

ks.test(-t.perm,"pfnorm")

t.test(x[y==-1],x[y==1])
1-pfnorm(-t.stat,mean=coefs[1],sd=coefs[2])
1-pfnorm(-t.stat)
1-pfnorm(-t.stat,mean(t.perm2),sd(t.perm2))
1-pfnorm(-t.stat,sqrt(2/pi),sqrt(1-2/pi))



plot(fits[,1],fits[,3])
cors1 <- sapply(2:200,function(x) cor(fits[,1],fits[,x]))
plot(cors1)
cors2 <- sapply(c(1,3:200),function(x) cor(fits[,2],fits[,x]))
plot(cors2)
plot(cors1[-2],cors2[-1])
cors3 <- sapply(c(1:2,4:200),function(x) cor(fits[,3],fits[,x]))
plot(cors3)
plot(cors1[-3],cors3[-1])
plot(cors2[-3],cors3[-2])
cors50 <- sapply(c(1:49,51:200),function(x) cor(fits[,50],fits[,x]))
pdf("p3.pdf")
par(mfrow=c(2,3))
plot(cors1[-1],cors2[-1])
plot(cors1[-2],cors3[-1])
plot(cors2[-2],cors3[-2])
plot(cors1[-49],cors50[-1])
plot(cors2[-49],cors50[-2])
plot(cors3[-49],cors50[-3])
dev.off()
cors1a <- sapply(2:200,function(x) cor(fits[1,],fits[x,]))
cors1a
plot(fits[1,],fits[2,])

corsx <- sapply(1:1000,function(y) cor(fits[y,],x))
plot(x,fits[1,])
plot(x,fits[6,])
sum(abs(t.perm*corsx-t.perm2))


pval <- rep(0,100)
for(j in 1:100){
x <- c(rnorm(10),rnorm(10))
y <- c(rep(0,10),rep(1,10))
sk <- vanilladot()
km <- kernelMatrix(sk,as.matrix(x))
ksvm1 <- ksvm(km,y)
t.stat <- t.test(fitted(ksvm1)[y==0],fitted(ksvm1)[y==1])$statistic
t.perm <- rep(0,100)
for(i in 1:100){
  y.perm <- sample(y)
  ksvm2 <- ksvm(km,y.perm)
  t.perm[i] <- t.test(fitted(ksvm2)[y.perm==0],fitted(ksvm2)[y.perm==1])$statistic
}
pval[j] <- mean(t.stat>=t.perm)
print(j)
}
hist(pval)
ks.test(pval,"punif")



##T2
x <- cbind(c(rnorm(100,mean=-.8),rnorm(100)),c(rnorm(100,mean=-.8),rnorm(100)))
y <- c(rep(-1,100),rep(1,100))
x <- matrix(rnorm(20*200,sd=10),nrow=200,ncol=20)



#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  when using below 
#-132.70 -107.90  -69.76  -79.20  -57.38  -35.66 
sigma=matrix(rep(.5,100^2),nrow=100)
diag(sigma)=rep(1,100)
sigma <- as.matrix(nearPD(sigma,corr=TRUE)$mat)
x <- rbind(t(rmvnorm(200,mean=rep(0,100),sigma=sigma)),t(rmvnorm(200,mean=rep(0,100),sigma=sigma)))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-140.30  -67.13  -48.01  -60.43  -41.66  -31.31 
sigma=matrix(rep(-.5,100^2),nrow=100)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -58.42  -40.63  -33.11  -35.63  -30.51  -21.78 
sigma=matrix(rep(.9,100^2),nrow=100)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -25.42  -22.41  -21.79  -21.21  -19.77  -17.92 
sigma=matrix(rep(.99,100^2),nrow=100)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-14.020 -11.920 -11.010 -11.040 -10.080  -7.741 
sigma=matrix(rep(.9,100^2),nrow=100)
diag(sigma)=rep(1,100)
sigma <- as.matrix(nearPD(sigma,corr=TRUE)$mat)
x <- rmvnorm(200,mean=rep(0,100),sigma=sigma)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -199.9  -196.6  -192.0  -191.4  -184.9  -182.2 
sigma=matrix(rep(.5,100^2),nrow=100)
diag(sigma)=rep(1,100)
sigma <- as.matrix(nearPD(sigma,corr=TRUE)$mat)
x <- rbind(t(rmvnorm(400,mean=rep(0,100),sigma=sigma)),t(rmvnorm(400,mean=rep(0,100),sigma=sigma)))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-1779.0  -893.6  -726.2  -899.1  -662.8  -496.3 
sigma=matrix(rep(.5,100^2),nrow=100)
diag(sigma)=rep(1,100)
sigma <- as.matrix(nearPD(sigma,corr=TRUE)$mat)
x <- rbind(t(rmvnorm(1000,mean=rep(0,100),sigma=sigma)),t(rmvnorm(1000,mean=rep(0,100),sigma=sigma)))
##even when independent, -1655.0  -727.1  -572.6  -747.9  -505.2  -473.8 

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -3.940  -3.443  -3.204  -3.188  -2.917  -2.467
x <- matrix(rnorm(200*10),nrow=200)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -12.99  -12.16  -11.58  -11.40  -10.60   -9.81 
x <- matrix(rnorm(200*100),nrow=200)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-142.70 -127.70  -70.76  -86.21  -59.74  -32.38
x <- matrix(rnorm(200*200),nrow=200)

library(ICSNP)
HotellingsT2(x[y==-1,],x[y==1,])

x <- matrix(rnorm(200*150),nrow=200)
x <- matrix(rnorm(200*20),nrow=200)
x <- matrix(rnorm(200*2),nrow=200)
x <- matrix(rnorm(200*1),nrow=200)
y <- c(rep(-1,100),rep(1,100))

sk <- vanilladot()
km <- kernelMatrix(sk,as.matrix(x))

ksvm1 <- ksvm(km,y)
t.stat <- t.test(fitted(ksvm1)[y==-1],fitted(ksvm1)[y==1])$statistic

nrep <- 1000
t.perm <- rep(0,nrep)
t.perm2 <- rep(0,nrep)
t.perm3 <- rep(0,nrep)
p.perm <- rep(0,nrep)
p.perm2 <- rep(0,nrep)
fits <- matrix(0,nrow=nrep,ncol=200)
for(i in 1:nrep){
  y.perm <- sample(y)
  ksvm2 <- ksvm(km,y.perm)
  t.perm[i] <- t.test(fitted(ksvm2)[y.perm==-1],fitted(ksvm2)[y.perm==1])$statistic
#  t.perm2[i] <- HotellingsT2(x[y.perm==-1,],x[y.perm==1,])$statistic[1]
  t.perm2[i] <- HotellingsT2(as.matrix(x[y.perm==-1,]),as.matrix(x[y.perm==1,]))$statistic[1]
  t.perm3[i] <- t.test(as.matrix(x[y.perm==-1,]),as.matrix(x[y.perm==1,]))$statistic
  p.perm[i] <- t.test(fitted(ksvm2)[y.perm==-1],fitted(ksvm2)[y.perm==1])$p.value
#  p.perm2[i] <- HotellingsT2(x[y.perm==-1,],x[y.perm==1,])$p.value
  fits[i,] <- fitted(ksvm2)
}
summary(t.perm)
plot(t.perm,sqrt(t.perm2));abline(a=0,b=-1)
lm1 <- lm(sqrt(t.perm2)~t.perm)
ad.test(lm1$residuals)
##are the residuals normal? p=.007
plot(sqrt(t.perm2),t.perm3)
plot(t.perm,t.perm3)
plot(p.perm,p.perm2)

library(MASS)
t.perm <- -t.perm2*(200-20-1)/((200-2)*20)
coefs <- c(20,179,0)

t.perm <- -t.perm2
coefs <- c(2,197,0)
#t.perm <- -1*rf(1000,10,10,10)

t.perm <- -t.perm^2

coefs <- fitdistr(-t.perm,"f",start=list("df1"=10,"df2"=10,"ncp"=5))$estimate
coefs <- fitdistr(-t.perm,"f",start=list("df1"=1,"df2"=1))$estimate
myseq <- seq(min(-t.perm),max(-t.perm),.01)
hist(-t.perm,prob=TRUE,main="black:fitted;red:density")
lines(myseq,df(myseq,coefs[1],coefs[2],coefs[3]),type="l")
lines(myseq,df(myseq,coefs[1],coefs[2]),type="l")
#lines(myseq,df(myseq,10,10,10),type="l",col="blue")
#lines(density(-t.perm,bw=2),col="red")

ks.test(-t.perm,"pf",coefs[1],coefs[2])
ks.test(-t.perm,"pf",coefs[1],coefs[2],coefs[3])

library(VGAM)
detach("package:VGAM")
qqplot(qfnorm(seq(1/nrep,1-1/nrep,length.out=nrep)),-t.perm)
fit=vglm(-t.perm~1,fam=fnormal1,trace=TRUE)
coefs <- Coef(fit)
myseq <- seq(floor(min(abs(t.perm))),ceiling(max(abs(t.perm))),.01)
plot(myseq,dfnorm(myseq,mean=coefs[1],sd=coefs[2]),type="l")
lines(density(-t.perm),col="red")
lines(myseq,dfnorm(myseq),type="l",col="blue")
1-pfnorm(-t.stat,mean=coefs[1],sd=coefs[2])




library(mvtnorm)
library(Matrix)
x <- rnorm(100)
y <- rnorm(100)
l <- c(rep(0,100),rep(1,100))
t.stat <- t.test(x,y)$statistic

t.stats <- 1:1000
for(i in 1:1000){
  l.perm <- sample(l)
  t.stats[i] <- t.test(c(x,y)*l.perm,c(x,y)*(1-l.perm))$statistic
}
sum(abs(t.stat)>abs(t.stats))
qqnorm(t.stats);qqline(t.stats)

sigma=matrix(rep(.5,100^2),nrow=100)
diag(sigma)=rep(1,100)
sigma <- as.matrix(nearPD(sigma,corr=TRUE)$mat)
sigma[1:5,1:5]
x <- rmvnorm(1,mean=rep(1,100),sigma=sigma)
y <- rmvnorm(1,mean=rep(0,100),sigma=sigma)
l <- c(rep(0,100),rep(1,100))
t.stat <- t.test(x,y)$statistic

t.stats <- 1:1000
for(i in 1:1000){
  l.perm <- sample(l)
  t.stats[i] <- t.test(c(x,y)*l.perm,c(x,y)*(1-l.perm))$statistic
}
sum(abs(t.stat)>abs(t.stats))
qqnorm(t.stats);qqline(t.stats)


sigma=matrix(rep(.5,5^2),nrow=5)
diag(sigma)=rep(1,5)
sigma <- as.matrix(nearPD(sigma,corr=TRUE)$mat)
sigma[1:5,1:5]
x <- as.vector(rmvnorm(20,mean=rep(0,5),sigma=sigma))
y <- as.vector(rmvnorm(20,mean=rep(0,5),sigma=sigma))
l <- c(rep(0,100),rep(1,100))
t.stat <- t.test(x,y)$statistic

t.stats <- 1:1000
for(i in 1:1000){
  l.perm <- sample(l)
  t.stats[i] <- t.test(c(x,y)*l.perm,c(x,y)*(1-l.perm))$statistic
}
sum(abs(t.stat)>abs(t.stats))
qqnorm(t.stats);qqline(t.stats)












