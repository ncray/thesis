library(kernlab)
library(stringkernels)

dat2A <- read.csv("P2A-2",header=FALSE,as.is=TRUE)
##unlist(lapply(dat2A[,2],nchar))
dat1B <- read.csv("P1B-2",header=FALSE,as.is=TRUE)
##unlist(lapply(dat1B[,2],nchar))
dat2 <- dat2A
dat1 <- dat1B

for(i in 1:nrow(dat2A)){
  dat2[i,2] <- tolower(gsub("\\.| ","",dat2A[i,2]))
}
for(i in 1:nrow(dat1B)){
  dat1[i,2] <- tolower(gsub("\\.| ","",dat1B[i,2]))
}

sk <- stringdot()
sk(dat2[1,2],dat2[2,2])
sk(dat2A[1,2],dat2A[2,2])
sk(".................... ee ....lVPGDlVllaaGDkvPAD",".......mdfptlssy...... ee .............kvPAD")
sk("eelVPGDlVllaaGDkvPAD","mdfptlssyeekvPAD")
sk("eelVPGDlVllaaGDkvPAD","mdfptlssyeekvPAd")
sk("eelvpgdlvllaagdkvpad","mdfptlssyeekvpad")

sk(dat1[1,2],dat1[2,2])
sk(dat2[1,2],dat1[2,2])

X <- c(dat1[,2],dat2[,2])
y <- c(rep(0,nrow(dat1)),rep(1,nrow(dat2)))
##y <- as.factor(c(rep(0,nrow(dat1)),rep(1,nrow(dat2))))

ksvm1 <- ksvm(as.list(X),y,kernel="stringdot")
y.perm <- sample(y)
ksvm2 <- ksvm(as.list(X),y.perm,kernel="stringdot")
predict(kvsm1,as.list(X))
fitted(ksvm1)
fitted(ksvm2)

t.test(fitted(ksvm1)[y==0],fitted(ksvm1)[y==1])$statistic
t.test(fitted(ksvm2)[y.perm==0],fitted(ksvm2)[y.perm==1])$statistic

t.perm <- rep(0,100)
for(i in 1:100){
  y.perm <- sample(y)
  ksvm2 <- ksvm(as.list(X),y.perm,kernel="stringdot")
  t.perm[i] <- t.test(fitted(ksvm2)[y.perm==0],fitted(ksvm2)[y.perm==1])$statistic
}

##kernelmatrix approach, much faster
sk <- stringdot()
km <- kernelMatrix(sk,as.list(X))
ksvm1 <- ksvm(km,y)
t.stat <- t.test(fitted(ksvm1)[y==0],fitted(ksvm1)[y==1])$statistic
t.perm <- rep(0,1000)
for(i in 1:1000){
  y.perm <- sample(y)
  ksvm2 <- ksvm(km,y.perm)
  t.perm[i] <- t.test(fitted(ksvm2)[y.perm==0],fitted(ksvm2)[y.perm==1])$statistic
  #t.perm[i] <- qnorm(pt(t.test(fitted(ksvm2)[y.perm==0],fitted(ksvm2)[y.perm==1],var.equal=TRUE)$statistic,length(y)-2))
}
library(ggplot2)
qplot(t.perm,geom="bar",main="Histogram of Permutation Null Distribution")
ggsave("pres1-1.pdf")
qplot(t.perm,geom="bar",binwidth=3,main="Permutation Null Distribution with Observed T-statistic (Red)") + scale_x_continuous(limits=c(-2400,0)) + geom_vline(xintercept=t.stat,col="red")
ggsave("pres2-1.pdf")

library(VGAM)
detach("package:VGAM")
qqplot(qfnorm(seq(1/1000,1-1/1000,length.out=1000)),-t.perm)
fit=vglm(-t.perm~1,fam=fnormal1,trace=TRUE)
coefs <- Coef(fit)
myseq <- seq(22,55,.01)
plot(myseq,dfnorm(myseq,mean=coefs[1],sd=coefs[2]),type="l")
lines(density(-t.perm,bw=1),col="red")
lines(myseq,dfnorm(myseq),type="l",col="blue")
1-pfnorm(-t.stat,mean=coefs[1],sd=coefs[2])
pdf("p7.pdf")
hist(-t.perm,prob=TRUE,breaks=10,main="black:fitted;red:density")
lines(myseq,dfnorm(myseq,mean=coefs[1],sd=coefs[2]),type="l")
lines(density(-t.perm,bw=10),col="red")
ks.test(-t.perm,"pfnorm")
ks.test(-t.perm,"pfnorm",coefs[1],coefs[2])
dev.off()
pdf("results-3.pdf")
plot(fitted(ksvm1))
dev.off()

pdf("results-4.pdf")
qqnorm(t.perm);qqline(t.perm)
dev.off()

plot(eigen(km)$vectors[,1])


library(twitteR)
obama <- userTimeline("BarackObama",n=500)
palin <- userTimeline("SarahPalinUSA",n=500)
obama2 <- obama
palin2 <- palin

for(i in 1:length(obama)){
  obama2[[i]] <- gsub("http.*","",obama[[i]]$getText())
  obama2[[i]] <- tolower(gsub("[^ a-zA-Z]","",obama2[[i]]))
}

for(i in 1:length(palin)){
  palin2[[i]] <- gsub("http.*","",palin[[i]]$getText())
  palin2[[i]] <- tolower(gsub("[^ a-zA-Z]","",palin2[[i]]))
}

which(lapply(obama2,nchar) < 5)
which(lapply(palin2,nchar) < 5)
obama2 <- unique(obama2)
palin2 <- unique(palin2)

X2 <- c(obama2,palin2)
y2 <- c(rep(0,length(obama2)),rep(1,length(palin2)))
y2 <- y2[-which(lapply(X2,nchar)<5)]
X2 <- X2[-which(lapply(X2,nchar)<5)]

sk <- stringdot()
km <- kernelMatrix(sk,as.list(X2))
ksvm1 <- ksvm(km,y2)
t.perm <- rep(0,1000)
t.stat <- t.test(fitted(ksvm1)[y2==0],fitted(ksvm1)[y2==1])$statistic
for(i in 1:1000){
  y.perm <- sample(y2)
  ksvm2 <- ksvm(km,y.perm)
  t.perm[i] <- t.test(fitted(ksvm2)[y.perm==0],fitted(ksvm2)[y.perm==1])$statistic
}
pdf("results-1.pdf")
plot(fitted(ksvm1))
dev.off()

t.test(fitted(ksvm1)[y2==0],fitted(ksvm1)[y2==1])$statistic
mean(t.test(fitted(ksvm1)[y2==0],fitted(ksvm1)[y2==1])$statistic > t.perm)

pdf("results-2.pdf")
qqnorm(t.perm);qqline(t.perm)
dev.off()
hist(t.perm)


qplot(t.perm,geom="bar",main="Histogram of Permutation Null Distribution")
ggsave("pres5.pdf")
qplot(t.perm,geom="bar",binwidth=3,main="Permutation Null Distribution with Observed T-statistic (Red)") + scale_x_continuous(limits=c(-200,-80)) + geom_vline(xintercept=t.stat,col="red")
ggsave("pres6.pdf")
