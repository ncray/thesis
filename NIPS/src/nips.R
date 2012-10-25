library(kernlab)
library(plyr)
library(ggplot2)
library(ICSNP)
library(reshape)
library(boot)
library(doMC)
registerDoMC(8)

Npts <- 1000
Npwr <- 100

computeT <- function(u, l, ...) t.test(u[l == 1], u[l == -1], var.equal = TRUE)$statistic
computeT2 <- function(u, l, ...) as.numeric(HotellingsT2(X = data.frame(u[l == 1, ]), Y = data.frame(u[l == -1, ]))$statistic)
computeKMMD <- function(u, l, ...){
  if(class(u) == "list"){
    kmmd(x = u[l == 1], y = u[l == -1], ...)@mmdstats[1]
  } else {
    kmmd(x = as.matrix(u[l == 1, ]), y = as.matrix(u[l == -1, ]), ...)@mmdstats[1]
  }
}
computeFS <- function(u, l, ...){
  ksvm.fit <- ksvm(x = u, y = l, C = .1, ...)
  ##print(ksvm.fit)
  as.numeric(computeT(fitted(ksvm.fit), l))
}

getData <- function(N, D = 1, delta = 1){
  u <- rbind(matrix(rnorm(N * D, 0), ncol = D), matrix(rnorm(N * D, delta), ncol = D))
  l <- c(rep(-1, N), rep(1, N))
  list("u" = u, "l" = l)
}

dat <- getData(10, 2)
u <- dat$u
l <- dat$l
## computeT(u, l)^2
## computeT2(u, l)
## computeKMMD(u, l)
## computeFS(u, l)


nullDist <- function(D = 1, N = 100, ...){
  dat <- getData(200, D)
  u <- dat$u
  l <- dat$l
  print(paste("D:", D, "--degree:", ""))  
  ldply(1:N, function(x){
    print(x)
    l.perm <- sample(l)
    data.frame("D" = D, "N" = N,
               "T2" = computeT2(u, l.perm),
               "sqrtT2" = sqrt(computeT2(u, l.perm)),
               "KMMD-l" = computeKMMD(u, l.perm, kernel = "vanilladot"),
               "FS-l" = computeFS(u, l.perm, kernel = "vanilladot"),
               "KMMD-rbf" = computeKMMD(u, l.perm, kernel = "rbfdot"),
               "FS-rbf" = computeFS(u, l.perm, kernel = "rbfdot"))
  }, .parallel = TRUE)
}
res <- ldply(c(1, 5, 10), function(D) nullDist(D = D, N = Npts))
res$D <- factor(res$D)
res.m <- melt(res, id.vars = c("D", "N"))

## ggplot(data = res.m, aes(x = value, fill = variable)) +
##   geom_density(alpha = .4) +
##   facet_wrap(~D)

p1 <- ggplot(data = res.m, aes(x = value, fill = D)) +
  geom_density(alpha = .4) +
  facet_wrap(~variable, scales = "free") +
  opts(title = "Null Distributions (Faceted by Statistic)")
png("null_dist.png", width = 800, height = 600)
print(p1)
dev.off()
##ggsave(file = "null_dist.png", plot = p1)
rejectT2 <- function(u, l, ...) as.numeric(HotellingsT2(X = data.frame(u[l == 1, ]), Y = data.frame(u[l == -1, ]))$p.value < .05)
rejectKMMD <- function(u, l, ...) as.numeric(computeKMMD(u, l, ...) > max(laply(1:19, function(i) computeKMMD(u, sample(l), ...))))
rejectFS <- function(u, l, ...) as.numeric(computeFS(u, l, ...) > max(laply(1:19, function(i) computeFS(u, sample(l), ...))))

## km <- kernelMatrix(vanilladot(), x = u)
## system.time(rejectFS(u, l, kernel = "vanilladot"))
## system.time(rejectFS(km, l))
## computeFS(km, l)
## ksvm(x = km, y = l)
## ksvm(x = u, y = l, kernel = "vanilladot")
## cor(fitted(ksvm(x = km, y = l, kernel = "vanilladot")), fitted(ksvm(x = u, y = l, kernel = "vanilladot")))

## rejectT2(u, l)
## rejectKMMD(u, l)
## rejectFS(u, l)
## l <- c(rep(-1, 10), rep(1, 10))
## res <- laply(1:1000, function(x) rejectFS(rnorm(20), l))

power <- function(D = 1, delta = 1, ...){
  print(paste("D:", D, "--delta:", delta))
  ldply(1:Npwr, function(x){  
    dat <- getData(20, D, delta)
    u <- dat$u
    l <- dat$l
    print(x)
    data.frame("D" = D, "delta" = delta,
               "T2" = rejectT2(u, l, ...),
               "KMMD" = rejectKMMD(u, l, ...),
               "FS" = rejectFS(u, l, ...))
  }, .parallel = TRUE)
}

res <- mdply(expand.grid("delta" = seq(0, 1.5, .25), "D" = c(1, 5, 10, 20)), function(delta, D) power(delta = delta, D = D, kernel = "vanilladot"))

res2 <- ddply(res, .(delta, D), function(df){
  means <- apply(df[, -(1:2)], 2, mean)
  lims <- c(.025, .975)
  bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
  bounds <- rbind(bootM(df$T2),
                  bootM(df$KMMD),
                  bootM(df$FS))
  data.frame("delta" = df$delta[1], "D" = df$D[1], "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = names(df[-(1:2)]))
})

p2 <- ggplot(res2, aes(x = delta, y = value, color = group, linetype = group)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .07)) +
  xlab(expression(Delta)) +
  facet_wrap(~D) +
  opts(title = "Power (Faceted by Dimension)")
##ggsave(p2, "power_normal.png")
png("power_normal.png", width = 800, height = 600)
print(p2)
dev.off()

rejectFS <- function(km, l) as.numeric(computeFS(km, l) > max(laply(1:19, function(i) computeFS(km, sample(l)))))
load(file = "obama")
load(file = "palin")

for(i in 1:length(obama)){
  obama[[i]] <- gsub("http.*","",obama[[i]]$getText())
  obama[[i]] <- tolower(gsub("[^ a-zA-Z]","",obama[[i]]))
}

for(i in 1:length(palin)){
  palin[[i]] <- gsub("http.*","",palin[[i]]$getText())
  palin[[i]] <- tolower(gsub("[^ a-zA-Z]","",palin[[i]]))
}

obama <- unique(obama)
palin <- unique(palin)
obama <- obama[which(lapply(obama, nchar) > 5)]
palin <- palin[which(lapply(palin, nchar) > 5)]

getTextSamp <- function(obama, palin, samp.size){
  list(u = c(obama[sample(1:length(obama), samp.size)], palin[sample(1:length(palin), samp.size)]),
       l = c(rep(-1, samp.size), rep(1, samp.size)))
}

## system.time(kmmd(x = u[l == 1], y = u[l == -1], kernel = "stringdot", kpar = list(length = 10),
##                  asymptotic = FALSE, ntimes = 10000))

## N <- 50
## dat <- getTextSamp(obama, palin, N)
## u <- dat$u
## l <- dat$l
## fitted(ksvm(x = u, y = l, kernel = "stringdot"))
## computeFS(u, l, kernel = "stringdot")
## system.time(rejectFS(u, l, kernel = "stringdot", kpar = list(length = 3)))
## system.time(rejectKMMD(u, l, kernel = "stringdot", kpar = list(type = "spectrum", length = 3)))

power <- function(N = 10, length, ...){
  print(paste("N:", N, "--length:", length))
  ldply(1:Npwr, function(x){  
    dat <- getTextSamp(obama, palin, N)
    u <- dat$u
    l <- dat$l
    ker <- stringdot(list(length, lambda = .5))
    ## print(ker)
    km <- kernelMatrix(ker, x = u)
    print(x)
    data.frame("N" = N, "length" = length,
                     "KMMD" = rejectKMMD(u, l, kernel = "stringdot", kpar = list(length, lambda = .5)),
                     ## "FS" = rejectFS(u, l, kernel = "stringdot", kpar = list(length = length, lambda = .5)))
                     "FS" = rejectFS(km, l))
  }, .parallel = TRUE)
}

## for(i in 1:100){
##   l.perm <- sample(l)
##   ##kmmd(x = u[l.perm == 1], y = u[l.perm == -1], asymptotic = TRUE, ntimes = 0)@mmdstats[1]
##   kmmd(x = u[l.perm == 1], y = u[l.perm == -1])@mmdstats[1]
## }

Npwr <- 8
## 4 minutes for Npwr = 8
system.time(res <- mdply(expand.grid("N" = seq(10, 50, 10), "length" = 1:4), function(N, length) power(N = N, length)))
res2 <- ddply(res, .(N, length), function(df){
  means <- apply(df[, -(1:2)], 2, mean)
  lims <- c(.025, .975)
  bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
  bounds <- rbind(bootM(df$KMMD),
                  bootM(df$FS))
  data.frame("N" = df$N[1], "length" = df$length[1], "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = names(df[-(1:2)]))
})

p3 <- ggplot(res2, aes(x = N, y = value, color = group, linetype = group)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 1)) +
  xlab("N") +
  facet_wrap(~length) +
  opts(title = "Power on Text Data (Faceted by K)")
##ggsave(p3, "power_string.png")
png("power_string.png", width = 800, height = 600)
print(p3)
dev.off()

library(jpeg)
rooster.path <- list.files("~/101_ObjectCategories/rooster_resized", full.names = TRUE)[1:45]
pigeon.path <- list.files("~/101_ObjectCategories/pigeon_resized", full.names = TRUE)[1:45]

rooster <- laply(rooster.path, function(path) as.vector(readJPEG(path)))
pigeon <- laply(pigeon.path, function(path) as.vector(readJPEG(path)))

getImageSamp <- function(rooster, pigeon, samp.size){
  u <- rbind(rooster[sample(1:nrow(rooster), samp.size), ], pigeon[sample(1:nrow(pigeon), samp.size), ])
  u <- sweep(u, 1, apply(u, 1, mean), "-")
  u <- sweep(u, 1, apply(u, 1, function(vec) sqrt(sum(vec^2))), "/")
  l <- c(rep(-1, samp.size), rep(1, samp.size))
  list(u = u, l = l)
}

## dat <- getImageSamp(rooster, pigeon, 10)
## u <- dat$u
## l <- dat$l
## system.time(computeFS(u, l))
## system.time(computeKMMD(u, l))
## system.time(rejectFS(u, l))
## system.time(rejectFS(u, l))
## system.time(rejectKMMD(u, l))
## system.time(rejectFS(kernelMatrix(rbfdot(), x = u), l))

power <- function(N = 10){
  print(paste("N:", N))
  ldply(1:Npwr, function(x){  
    dat <- getImageSamp(rooster, pigeon, N)
    u <- dat$u
    l <- dat$l
    print(x)
    data.frame("N" = N, "length" = NA,
               "KMMDp1" = rejectKMMD(u, l, kernel = "polydot", kpar = list(degree = 1, offset = 1)),
               "FSp1" = rejectFS(kernelMatrix(polydot(degree = 1, offset = 1), x = u), l),
               "KMMDp2" = rejectKMMD(u, l, kernel = "polydot", kpar = list(degree = 2, offset = 1)),
               "FSp2" = rejectFS(kernelMatrix(polydot(degree = 2, offset = 1), x = u), l),
##               "KMMDp3" = rejectKMMD(u, l, kernel = "polydot", kpar = list(degree = 3, offset = 1)),
##               "FSp3" = rejectFS(kernelMatrix(polydot(degree = 3, offset = 1), x = u), l),
               "KMMDrbf" = rejectKMMD(u, l, kernel = "rbfdot"),
               "FSrbf" = rejectFS(kernelMatrix(rbfdot(sigest(u)[2]), x = u), l))
  }, .parallel = TRUE)
}

fitted(ksvm(x = km, y = l))
fitted(ksvm(x = km, y = l, C = .01))
fitted(ksvm(u, l, C = .01, kernel = "polydot", kpar = list(degree = 3, offset = 1)))
## fitted(ksvm(x = km, y = c(l, 1)))
## computeFS(km, c(l, 1))
computeFS(km, l, C = .01)
computeFS(u, l, kernel = "polydot", kpar = list(degree = 3, offset = 1))
rejectFS(u, l, kernel = "polydot", kpar = list(degree = 2, offset = 1))

km <- kernelMatrix(polydot(degree = 2, offset = 1), x = u)
km <- kernelMatrix(rbfdot(), x = u)
km <- kernelMatrix(rbfdot(sigest(u)[2]), x = u)
km <- kernelMatrix(polydot(degree = 3, offset = 1), x = scale(u))
computeFS(km, l)
## computeFS(km, l, C = .01)
sort(laply(1:19, function(i) computeFS(km, sample(l))))
## sort(laply(1:19, function(i) computeFS(km, sample(l), C = .01)))
## sort(laply(1:19, function(i) computeFS(u, sample(l)), .parallel = TRUE))
quantile(laply(1:1000, function(i) computeFS(km, sample(l))), .95)

##Npwr = 8 takes 9 minutes
Npwr <- 400
system.time(res <- ldply(seq(10, 40, 10), power))

res2 <- ddply(res, .(N), function(df){
  means <- apply(df[, -(1:2)], 2, mean)
  lims <- c(.025, .975)
  bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
  bounds <- rbind(bootM(df$KMMDp1),
                  bootM(df$FSp1),
                  bootM(df$KMMDp2),
                  bootM(df$FSp2),
##                  bootM(df$KMMDp3),
##                  bootM(df$FSp3),
                  bootM(df$KMMDrbf),
                  bootM(df$FSrbf))
  data.frame("N" = df$N[1], "length" = df$length[1], "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = names(df[-(1:2)]))
})
library(stringr)
res2 <- transform(res2, test = str_detect(group, "KMMD"))
res2$test[res2$test == TRUE] <- "KMMD"
res2$test[res2$test == FALSE] <- "FS"
ker.list <- str_split(res2$group, "KMMD|FS")
res2$ker <- laply(ker.list, `[[`, 2)

## p4 <- ggplot(res2, aes(x = N, y = value, color = group, linetype = group)) +
##   geom_line() + 
##   geom_errorbar(aes(ymin = lower, ymax = upper, width = 1)) +
##   xlab("N") +
##   opts(title = "Power on Image Data")

p4 <- ggplot(res2, aes(x = N, y = value, color = ker, linetype = ker)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 1)) +
  xlab("N") +
  opts(title = "Power on Image Data") +
  facet_grid(~test)
##ggsave(p4, "power_image.png")
png("power_image.png", width = 800, height = 600)
print(p4)
dev.off()



computeFS <- function(u, l, C, eps){
  ksvm.fit <- ksvm(x = u, y = l, C = C, epsilon = eps)
  ##print(ksvm.fit)
  as.numeric(computeT(fitted(ksvm.fit), l))
}
rejectFS <- function(km, l, C, eps) as.numeric(computeFS(km, l, C, eps) > max(laply(1:19, function(i) computeFS(km, sample(l), C, eps))))

power <- function(N, C, eps){
  print(paste("N:", N, "C:", C, "eps:", eps))
  ldply(1:Npwr, function(x){  
    dat <- getTextSamp(obama, palin, N)
    u <- dat$u
    l <- dat$l
    ker <- stringdot(list(length = 3))
    km <- kernelMatrix(ker, x = u)
    print(x)
    data.frame("N" = N, "C" = C, "eps" = eps,
                     "FS" = rejectFS(km, l, C, eps))
  }, .parallel = TRUE)
}

Npwr <- 1000
## 1 minutes for Npwr = 8
system.time(res <- mdply(expand.grid(N = seq(10, 50, 10), C = c(.1, 1, 10), eps = c(.01, .05, .1, .5)), power))
res2 <- ddply(res, .(N, C, eps), function(df){
  means <- mean(df$FS)
  lims <- c(.025, .975)
  bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
  bounds <- bootM(df$FS)
  data.frame("N" = df$N[1], "C" = df$C[1], "eps" = df$eps[1], "value" = means, "lower" = bounds[1], "upper" = bounds[2])
})

res2$C <- factor(res2$C)
res2$eps <- factor(res2$eps)
p3 <- ggplot(res2, aes(x = N, y = value, color = eps)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 1)) +
  xlab("N") +
  facet_wrap(~C) +  
  opts(title = "Power on Text Data")
##ggsave(p3, "power_string.png")
png("power_kpar.png", width = 800, height = 600)
print(p3)
dev.off()
