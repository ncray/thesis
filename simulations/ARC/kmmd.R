library(kernlab)
simNum2 <- function(num,mean){
  x <- rnorm(num)
  y <- rnorm(num,mean=mean)
  #c(kmmd(matrix(x),matrix(y),asymptotic=TRUE,kernel="vanilladot")@AsympH0,kmmd(matrix(x),matrix(y),kernel="vanilladot")@H0,t.test(x,y)$p.value<.05,simOne(x,y),twot.perm(x,y)<.05)
  c(kmmd(matrix(x),matrix(y),asymptotic=TRUE,kernel="vanilladot")@AsympH0,kmmd(matrix(x),matrix(y),kernel="vanilladot")@H0,t.test(x,y)$p.value<.05,twot.perm(x,y)<.05)
}

#equal number of samples
num <- 100
x <- rnorm(num)
y <- rnorm(num,mean=.6)
ker <- vanilladot()
km <- kernelMatrix(ker,as.matrix(c(x,y)))

nr <- nrow(km)
H <- diag(rep(1,nr))-1/nr*matrix(1,nr,nr)
km.c <- H%*%km%*%H
evs <- 1/(.5*nr)*as.real(eigen(km.c)$val)

kmmd1 <- kmmd(matrix(x), matrix(y), asymptotic=TRUE, kernel="vanilladot")

t.test(x,y)

m <- num
Kxx <- kernelMatrix(ker,as.matrix(x))
Kyy <- kernelMatrix(ker,as.matrix(y))
Kxy <- kernelMatrix(ker,as.matrix(x),as.matrix(y))
mmdu.sq <- 1/(m * (m-1)) * (sum(Kxx-diag(diag(Kxx))) + sum(Kyy-diag(diag(Kyy))) - 2*sum(Kxy-diag(diag(Kxy))))

kmmd1
sqrt(mmdu.sq)

n.sim <- 10000
evs.filt <- evs[evs > 1e-10]
dist.samps <- sapply(evs.filt, function(x) rnorm(n.sim, sd = sqrt(2 * x)) ^ 2 - 2 * x)
dist.ecdf <- ecdf(rowSums(dist.samps))
1 - dist.ecdf(m*mmdu.sq)


##try matrix
num <- 400
ndim <- 10
x <- matrix(rnorm(num * ndim), nrow = num)
y <- matrix(rnorm(num * ndim, mean=.3), nrow = num)
ker <- vanilladot()
km <- kernelMatrix(ker,rbind(x,y))

nr <- nrow(km)
H <- diag(rep(1,nr))-1/nr*matrix(1,nr,nr)
km.c <- H%*%km%*%H
evs <- 1/(.5*nr)*as.real(eigen(km.c)$val)

kmmd1 <- kmmd(x, y, asymptotic=TRUE, kernel="vanilladot")

library(ICSNP)
HotellingsT2(x,y)

m <- num
Kxx <- kernelMatrix(ker,x)
Kyy <- kernelMatrix(ker,y)
Kxy <- kernelMatrix(ker,x,y)
mmdu.sq <- 1/(m * (m-1)) * (sum(Kxx-diag(diag(Kxx))) + sum(Kyy-diag(diag(Kyy))) - 2*sum(Kxy-diag(diag(Kxy))))

kmmd1
sqrt(mmdu.sq)

n.sim <- 10000
evs.filt <- evs[evs > 1e-10]
dist.samps <- sapply(evs.filt, function(x) rnorm(n.sim, sd = sqrt(2 * x)) ^ 2 - 2 * x)
dist.ecdf <- ecdf(rowSums(dist.samps))
1 - dist.ecdf(m*mmdu.sq)


kmmd2 <- function(x, y, kernel, n.sim = 10000){
  num <- nrow(x)
  ker <- match.fun(kernel)()
  km <- kernelMatrix(ker,rbind(x,y))

  nr <- nrow(km)
  H <- diag(rep(1,nr))-1/nr*matrix(1,nr,nr)
  km.c <- H%*%km%*%H
  evs <- 1/(.5*nr)*as.real(eigen(km.c)$val)

  m <- num
  Kxx <- kernelMatrix(ker,x)
  Kyy <- kernelMatrix(ker,y)
  Kxy <- kernelMatrix(ker,x,y)
  mmdu.sq <- 1/(m * (m-1)) * (sum(Kxx-diag(diag(Kxx))) + sum(Kyy-diag(diag(Kyy))) - 2*sum(Kxy-diag(diag(Kxy))))

  evs.filt <- evs[evs > 1e-10]
  dist.samps <- sapply(evs.filt, function(x) rnorm(n.sim, sd = sqrt(2 * x)) ^ 2 - 2 * x)
  dist.ecdf <- ecdf(rowSums(dist.samps))
  1 - dist.ecdf(m*mmdu.sq)
}


simNum2 <- function(num,mean){
  x <- rnorm(num)
  y <- rnorm(num,mean=mean)
  #c(kmmd(matrix(x),matrix(y),asymptotic=TRUE,kernel="vanilladot")@AsympH0,kmmd(matrix(x),matrix(y),kernel="vanilladot")@H0,t.test(x,y)$p.value<.05,simOne(x,y),twot.perm(x,y)<.05)
  c(kmmd(matrix(x),matrix(y),asymptotic=TRUE,kernel="vanilladot")@AsympH0,kmmd(matrix(x),matrix(y),kernel="vanilladot")@H0,kmmd2(matrix(x),matrix(y),"vanilladot") < .05,t.test(x,y)$p.value<.05,twot.perm(x,y)<.05)
}

res.df <- data.frame()
deltas <- seq(0,2,.25)
for(del in deltas){
  res <- mclapply(rep(20,1000),simNum2,mean=del)
  means <- apply(do.call(rbind,res),2,mean)
  lowers <- means - sqrt(means*(1-means)/1000)
  uppers <- means + sqrt(means*(1-means)/1000)
  Tests <- c("MMDasymp","MMDradem","MMDgram","t","FSVMP")
  dels <- rep(del,5)
  res.df <- rbind(res.df,data.frame(dels,means,lowers,uppers,Tests))
}

qplot(dels,means,data=res.df,geom="line",col=Tests,linetype=Tests,ylab="Power",xlab=expression(Delta),main="Power for Normal Location Alternative") + geom_errorbar(aes(ymin=lowers,ymax=uppers))
ggsave("nips1-1.pdf",width=8,height=6)

res1 <- 1:1000
for(i in 1:1000){
 res1[i] <- kmmd2(matrix(rnorm(100)),matrix(rnorm(100)),"vanilladot")
 res1[i] <- kmmd2(matrix(rnorm(1000),nrow=10),matrix(rnorm(1000),nrow=10),"vanilladot") 
}
mean(res1 < .05)




res.df <- data.frame()
numobs <- seq(50,300,50)
for(num in numobs){
  res <- mclapply(rep(num,100),simNum)
  means <- apply(do.call(rbind,res),2,mean)
  lowers <- means - 2*sqrt(means*(1-means)/100)
  uppers <- means + 2*sqrt(means*(1-means)/100)
  Tests <- c("MMDasymp","MMDradem","MMDgram","FSVMP")
  num.obs <- rep(num,4)
  res.df <- rbind(res.df,data.frame(num.obs,means,lowers,uppers,Tests))
}
simNum <- function(num){
  #s1 <- obama2[sample(1:length(obama2),num)]
  #s2 <- palin2[sample(1:length(palin2),num)]
  #samp <- sample(1:length(obama2),2*num)
  #s1 <- obama2[samp[1:num]]
  #s2 <- obama2[samp[(num+1):(2*num)]]
  s1 <- obama2[sample(1:length(obama2),num,replace=TRUE)]
  s2 <- obama2[sample(1:length(obama2),num,replace=TRUE)]
  c(kmmd(s1,s2,asymptotic=TRUE)@AsympH0,kmmd(s1,s2)@H0,kmmd3(s1,s2,"stringdot")<.05,simOne(s1,s2))
}
simOne <- function(obama2,palin2){
  X2 <- c(obama2,palin2)
  y2 <- c(rep(0,length(obama2)),rep(1,length(palin2)))
  sk <- stringdot()
  km <- kernelMatrix(sk,as.list(X2))
  ksvm1 <- ksvm(km,y2)
  t.perm <- rep(0,19)
  for(i in 1:19){
    y.perm <- sample(y2)
    ksvm2 <- ksvm(km,y.perm)
    t.perm[i] <- t.test(fitted(ksvm2)[y.perm==0],fitted(ksvm2)[y.perm==1])$statistic
  }
  sum(t.test(fitted(ksvm1)[y2==0],fitted(ksvm1)[y2==1])$statistic > t.perm) <= 0
}

qplot(num.obs,means,data=res.df,geom="line",col=Tests,linetype=Tests,ylab="Power",xlab="Number of Observations",main="Power Versus Number of Observations (Twitter Data)") + geom_errorbar(aes(ymin=lowers,ymax=uppers))
ggsave("nips2-1.pdf",width=8,height=6)

qplot(num.obs,means,data=res.df,geom="line",col=Tests,linetype=Tests,ylab="Power",xlab="Number of Observations",main="Level Versus Number of Observations (Twitter Data)") + geom_errorbar(aes(ymin=lowers,ymax=uppers)) + geom_hline(aes(yintercept=.05))
ggsave("nips3-1.pdf",width=8,height=6)


kmmd3 <- function(x, y, kernel, n.sim = 1000){
  num <- length(x)
  ker <- match.fun(kernel)()
  km <- kernelMatrix(ker,c(x,y))

  nr <- nrow(km)
  H <- diag(rep(1,nr))-1/nr*matrix(1,nr,nr)
  km.c <- H%*%km%*%H
  evs <- 1/(.5*nr)*as.real(eigen(km.c)$val)

  m <- num
  Kxx <- kernelMatrix(ker,x)
  Kyy <- kernelMatrix(ker,y)
  Kxy <- kernelMatrix(ker,x,y)
  mmdu.sq <- 1/(m * (m-1)) * (sum(Kxx-diag(diag(Kxx))) + sum(Kyy-diag(diag(Kyy))) - 2*sum(Kxy-diag(diag(Kxy))))

  evs.filt <- evs[evs > 1e-10]
  dist.samps <- sapply(evs.filt, function(x) rnorm(n.sim, sd = sqrt(2 * x)) ^ 2 - 2 * x)
  dist.ecdf <- ecdf(rowSums(dist.samps))
  1 - dist.ecdf(m*mmdu.sq)
}

##example showing mmdu.sq can be negative
x <- matrix(c(1,-1))
y <- matrix(c(-1,1))
ker <- vanilladot()
m <- nrow(x)
Kxx <- kernelMatrix(ker,x)
Kyy <- kernelMatrix(ker,y)
Kxy <- kernelMatrix(ker,x,y)
mmdu.sq <- 1/(m * (m-1)) * (sum(Kxx-diag(diag(Kxx))) + sum(Kyy-diag(diag(Kyy))) - 2*sum(Kxy-diag(diag(Kxy))))

