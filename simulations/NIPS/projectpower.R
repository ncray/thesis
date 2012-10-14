library(kernlab)
##library(stringkernels) removed from CRAN
library(twitteR)
load(file = "obama")
load(file = "palin")
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
obama2 <- obama2[which(lapply(obama2,nchar) > 5)]
palin2 <- palin2[which(lapply(palin2,nchar) > 5)]

X2 <- c(obama2,palin2)
y2 <- c(rep(0,length(obama2)),rep(1,length(palin2)))
##y2 <- y2[-which(lapply(X2,nchar)<5)]
##X2 <- X2[-which(lapply(X2,nchar)<5)]

sk <- stringdot()
km <- kernelMatrix(sk,as.list(X2))
ksvm1 <- ksvm(km,y2)
t.perm <- rep(0,1000)
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

simOne(obama2[1:200],palin2[1:200])

res <- rep(0,10)
for(i in 1:10){
  res[i] <- simOne(obama2[sample(1:length(obama2),250)],palin2[sample(1:length(palin2),250)])
  print(i)
}

#.05 at 50, .08 at 100, .09 at 150, .35 at 200, .7 at 250, 1 at 300

simNum <- function(num){
  simOne(obama2[sample(1:length(obama2),num)],palin2[sample(1:length(palin2),num)])
  #samp <- sample(1:length(obama2),2*num)
  #simOne(obama2[samp[1:num]],obama2[samp[(num+1):(2*num)]])
}

res <- unlist(lapply(rep(150,50),simNum))

library(multicore)
res <- unlist(mclapply(rep(20,1000),simNum))

system.time(res <- unlist(mclapply(rep(300,20),simNum))) ##takes 77s

kmmd(obama2[sample(1:length(obama2),500)],palin2[sample(1:length(palin2),500)])
kmmd(obama2[sample(1:length(obama2),500)],palin2[sample(1:length(palin2),500)])
kmmd(as.matrix(rnorm(150)),as.matrix(rnorm(150,mean=1)),asymptotic=TRUE)

numobs <- seq(50,300,50)
res.df <- data.frame("numobs"=numobs,"power"=1:6,"lower"=1:6,"upper"=1:6)
for(i in 1:6){
  num <- 200
  res <- unlist(mclapply(rep(numobs[i],num),simNum))
  res.df[i,2] <- mean(res)
  res.df[i,3] <- mean(res) - sqrt(mean(res)*(1-mean(res))/num)
  res.df[i,4] <- mean(res) + sqrt(mean(res)*(1-mean(res))/num)
}

qplot(numobs,power,data=res.df,geom="line") + geom_errorbar(aes(ymin=lower,ymax=upper))
ggsave("pres7.pdf")

simNum2 <- function(num){
  #kmmd(as.matrix(rnorm(num)),as.matrix(rnorm(num)),asymptotic=TRUE)@AsympH0
  kmmd(matrix(rnorm(num*1),nrow=num),matrix(rnorm(num*1),nrow=num),asymptotic=TRUE)@AsympH0
}
res <- unlist(mclapply(rep(20,100),simNum2))
mean(res)

simNum2 <- function(num){
  #kmmd(obama2[sample(1:length(obama2),num)],palin2[sample(1:length(palin2),num)])@H0
  #kmmd(obama2[sample(1:length(obama2),num)],palin2[sample(1:length(palin2),num)],asymptotic=TRUE)@AsympH0
  samp <- sample(1:length(obama2),2*num)
  #kmmd(obama2[samp[1:num]],obama2[samp[(num+1):(2*num)]],asymptotic=TRUE)@AsympH0
  kmmd(obama2[samp[1:num]],obama2[samp[(num+1):(2*num)]],kernel="stringdot", kpar=list(type="spectrum",length=1),asymptotic=TRUE)@AsympH0
  #kmmd(obama2[samp[1:num]],obama2[samp[(num+1):(2*num)]])@H0
  #kmmd(obama2[sample(1:length(obama2),num)],obama2[sample(1:length(obama2),num)],asymptotic=TRUE)@AsympH0
}


numobs <- seq(50,300,50)
res.df2 <- data.frame("numobs"=numobs,"power"=1:6,"lower"=1:6,"upper"=1:6)
for(i in 1:6){
  num <- 200
  res <- unlist(mclapply(rep(numobs[i],num),simNum2))
  res.df2[i,2] <- mean(res)
  res.df2[i,3] <- mean(res) - sqrt(mean(res)*(1-mean(res))/num)
  res.df2[i,4] <- mean(res) + sqrt(mean(res)*(1-mean(res))/num)
}

samp <- sample(1:1196,1196/2)
kmmd(obama2[samp],obama2[-samp],asymptotic=TRUE)

simNum2 <- function(num,mean){
  x <- rnorm(num)
  y <- rnorm(num,mean=mean)
  #c(kmmd(matrix(x),matrix(y),asymptotic=TRUE,kernel="vanilladot")@AsympH0,kmmd(matrix(x),matrix(y),kernel="vanilladot")@H0,t.test(x,y)$p.value<.05,simOne(x,y),twot.perm(x,y)<.05)
  c(kmmd(matrix(x),matrix(y),asymptotic=TRUE,kernel="vanilladot")@AsympH0,kmmd(matrix(x),matrix(y),kernel="vanilladot")@H0,t.test(x,y)$p.value<.05,twot.perm(x,y)<.05)
}

do.call(rbind,res)

res <- mclapply(rep(20,50),simNum2,mean=1)
apply(do.call(rbind,res),2,mean)
#0.055 0.000 0.055 0.055 0.056
#[1] 0.864 0.000 0.869 0.798 0.870 19
#[1] 0.826 0.000 0.841 0.843 0.847 999

means <- apply(do.call(rbind,res),2,mean)
lowers <- means - sqrt(means*(1-means))/1000
uppers <- means + sqrt(means*(1-means))/1000

res.df <- data.frame()
deltas <- seq(0,2,.25)
for(del in deltas){
  res <- mclapply(rep(20,1000),simNum2,mean=del)
  means <- apply(do.call(rbind,res),2,mean)
  lowers <- means - sqrt(means*(1-means)/1000)
  uppers <- means + sqrt(means*(1-means)/1000)
  Tests <- c("MMDasymp","MMDradem","t","FSVMP")
  dels <- rep(del,4)
  res.df <- rbind(res.df,data.frame(dels,means,lowers,uppers,Tests))
}

qplot(dels,means,data=res.df,geom="line",col=Tests,linetype=Tests,ylab="Power",xlab=expression(Delta),main="Power for Normal Location Alternative") + geom_errorbar(aes(ymin=lowers,ymax=uppers))
ggsave("nips1.pdf",width=8,height=6)

simOne <- function(obama2,palin2){
  X2 <- c(obama2,palin2)
  y2 <- c(rep(0,length(obama2)),rep(1,length(palin2)))
  sk <- vanilladot()
  km <- kernelMatrix(sk,as.matrix(X2))
  ksvm1 <- ksvm(km,y2)
  t.perm <- rep(0,19)
  for(i in 1:19){
    y.perm <- sample(y2)
    ksvm2 <- ksvm(km,y.perm)
    t.perm[i] <- t.test(fitted(ksvm2)[y.perm==0],fitted(ksvm2)[y.perm==1])$statistic
  }
  mean(t.test(fitted(ksvm1)[y2==0],fitted(ksvm1)[y2==1])$statistic > t.perm) <= .05
}

twot.perm <- function (x1 = two65$ambient, x2 = two65$heated, nsim = 2000, 
    plotit = FALSE) 
{
    n1 <- length(x1)
    n2 <- length(x2)
    n <- n1 + n2
    x <- c(x1, x2)
    dbar <- mean(x2) - mean(x1)
    z <- array(, nsim)
    for (i in 1:nsim) {
        mn <- sample(n, n2, replace = FALSE)
        dbardash <- mean(x[mn]) - mean(x[-mn])
        z[i] <- dbardash
    }
    pval <- (sum(z >= abs(dbar)) + sum(z <= -abs(dbar)))/nsim
    pval
}



res.df <- data.frame()
numobs <- seq(50,300,50)
for(num in numobs){
  res <- mclapply(rep(num,100),simNum)
  means <- apply(do.call(rbind,res),2,mean)
  lowers <- means - 2*sqrt(means*(1-means)/100)
  uppers <- means + 2*sqrt(means*(1-means)/100)
  Tests <- c("MMDasymp","MMDradem","FSVMP")
  num.obs <- rep(num,3)
  res.df <- rbind(res.df,data.frame(num.obs,means,lowers,uppers,Tests))
}
simNum <- function(num){
  s1 <- obama2[sample(1:length(obama2),num)]
  s2 <- palin2[sample(1:length(palin2),num)]
  #samp <- sample(1:length(obama2),2*num)
  #s1 <- obama2[samp[1:num]]
  #s2 <- obama2[samp[(num+1):(2*num)]]
  #s1 <- obama2[sample(1:length(obama2),num,replace=TRUE)]
  #s2 <- obama2[sample(1:length(obama2),num,replace=TRUE)]
  c(kmmd(s1,s2,asymptotic=TRUE)@AsympH0,kmmd(s1,s2)@H0,simOne(s1,s2))
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
ggsave("nips2.pdf",width=8,height=6)

qplot(num.obs,means,data=res.df,geom="line",col=Tests,linetype=Tests,ylab="Power",xlab="Number of Observations",main="Level Versus Number of Observations (Twitter Data)") + geom_errorbar(aes(ymin=lowers,ymax=uppers)) + geom_hline(aes(yintercept=.05))
ggsave("nips3.pdf",width=8,height=6)




paste(letters[1:10],collapse="")
dist1 <- runif(26)
dist1 <- dist1/sum(dist1)
dist2 <- runif(26)
dist2 <- dist2/sum(dist2)

obama2 <- list()
palin2 <- list()
for(i in 1:500){
  obama2[[i]] <- paste(letters[sample(1:26,rpois(1,300),prob=dist1,replace=TRUE)],collapse="")
  palin2[[i]] <- paste(letters[sample(1:26,rpois(1,300),prob=dist2,replace=TRUE)],collapse="")
}
