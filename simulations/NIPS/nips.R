library(kernlab)
library(plyr)
library(ggplot2)
library(ICSNP)
library(reshape)
library(boot)
library(doMC)
registerDoMC(2)

Npts <- 10
Npwr <- 2

computeT <- function(u, l, ...) t.test(u[l == 1], u[l == -1], var.equal = TRUE)$statistic
computeT2 <- function(u, l, ...) as.numeric(HotellingsT2(X = data.frame(u[l == 1, ]), Y = data.frame(u[l == -1, ]))$statistic)
computeKMMD <- function(u, l, ...){
  if(class(u) == "list"){
    kmmd(x = u[l == 1], y = u[l == -1], ...)@mmdstats[1]
  } else {
    kmmd(x = matrix(u[l == 1, ]), y = matrix(u[l == -1, ]), ...)@mmdstats[1]
  }
}
computeFS <- function(u, l, ...) as.numeric(computeT(fitted(ksvm(x = u, y = l, ...)), l))

getData <- function(N, D = 1, delta = 1){
  u <- rbind(matrix(rnorm(N * D, 0), ncol = D), matrix(rnorm(N * D, delta), ncol = D))
  l <- c(rep(-1, N), rep(1, N))
  list("u" = u, "l" = l)
}

## dat <- getData(10, 2)
## u <- dat$u
## l <- dat$l
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
res.m <- melt(res, id.vars = c("D", "N"))

## ggplot(data = res.m, aes(x = value, fill = variable)) +
##   geom_density(alpha = .4) +
##   facet_wrap(~D)

ggplot(data = res.m, aes(x = value, fill = factor(D))) +
  geom_density(alpha = .4) +
  facet_wrap(~variable, scales = "free") +
  opts(title = "Null Distributions (Faceted by Statistic)")
ggsave("null_dist.png")
rejectT2 <- function(u, l, ...) as.numeric(HotellingsT2(X = data.frame(u[l == 1, ]), Y = data.frame(u[l == -1, ]))$p.value < .05)
rejectKMMD <- function(u, l, ...) as.numeric(computeKMMD(u, l, ...) > max(laply(1:19, function(i) computeKMMD(u, sample(l), ...))))
rejectFS <- function(u, l, ...) as.numeric(computeFS(u, l, ...) > max(laply(1:19, function(i) computeFS(u, sample(l), ...))))
## rejectT2(u, l)
## rejectKMMD(u, l)
## rejectFS(u, l)
l <- c(rep(-1, 10), rep(1, 10))
res <- laply(1:1000, function(x) rejectFS(rnorm(20), l))

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

ggplot(res2, aes(x = delta, y = value, color = group, linetype = group)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .07)) +
  xlab(expression(Delta)) +
  facet_wrap(~D) +
  opts(title = "Power (Faceted by Dimension)")
ggsave("power_normal.png")

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

## dat <- getTextSamp(obama, palin, N)
## u <- dat$u
## l <- dat$l
## fitted(ksvm(x = u, y = l, kernel = "stringdot"))
## computeFS(u, l, kernel = "stringdot")
## rejectFS(u, l, kernel = "stringdot", kpar = list(length = 3))
## rejectKMMD(u, l, kernel = "stringdot", kpar = list(type = "spectrum", length = 3))

power <- function(N = 10, kpar = list(length = 4), ...){
  print(paste("N:", N, "--length:", kpar))
  ldply(1:Npwr, function(x){  
    dat <- getTextSamp(obama, palin, N)
    u <- dat$u
    l <- dat$l
    print(x)
    data.frame("N" = N, "length" = kpar,
               "KMMD" = rejectKMMD(u, l, kernel = "stringdot", kpar = kpar),
               "FS" = rejectFS(u, l, kernel = "stringdot", kpar = kpar))
  }, .parallel = TRUE)
}

res <- mdply(expand.grid("N" = seq(10, 20, 10), "length" = 1:4), function(N, length) power(N = N, kpar = list(length = length)))
res2 <- ddply(res, .(N, length), function(df){
  means <- apply(df[, -(1:2)], 2, mean)
  lims <- c(.025, .975)
  bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
  bounds <- rbind(bootM(df$KMMD),
                  bootM(df$FS))
  data.frame("N" = df$N[1], "length" = df$length[1], "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = names(df[-(1:2)]))
})

ggplot(res2, aes(x = N, y = value, color = group, linetype = group)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 1)) +
  xlab("N") +
  facet_wrap(~length) +
  opts(title = "Power on Text Data (Faceted by K)")
ggsave("power_string.png")


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

dat <- getImageSamp(rooster, pigeon, 10)
u <- dat$u
l <- dat$l
computeFS(u, l)
rejectFS(u, l)

power <- function(N = 10){
  print(paste("N:", N))
  ldply(1:Npwr, function(x){  
    dat <- getImageSamp(rooster, pigeon, N)
    u <- dat$u
    l <- dat$l
    print(x)
    data.frame("N" = N, "length" = NA,
               "KMMDp1" = rejectKMMD(u, l, kernel = "polydot", kpar = list(degree = 1, offset = 1)),
               "FSp1" = rejectFS(u, l, kernel = "polydot", kpar = list(degree = 1, offset = 1)),
               "KMMDp2" = rejectKMMD(u, l, kernel = "polydot", kpar = list(degree = 2, offset = 1)),
               "FSp2" = rejectFS(u, l, kernel = "polydot", kpar = list(degree = 2, offset = 1)),
               "KMMDp3" = rejectKMMD(u, l, kernel = "polydot", kpar = list(degree = 3, offset = 1)),
               "FSp3" = rejectFS(u, l, kernel = "polydot", kpar = list(degree = 3, offset = 1)),
               "KMMDrbf" = rejectKMMD(u, l, kernel = "rbfdot"),
               "FSrbf" = rejectFS(u, l, kernel = "rbfdot"))
  }, .parallel = TRUE)
}

res <- ldply(seq(10, 40, 10), power)

res2 <- ddply(res, .(N), function(df){
  means <- apply(df[, -(1:2)], 2, mean)
  lims <- c(.025, .975)
  bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
  bounds <- rbind(bootM(df$KMMDp1),
                  bootM(df$FSp1),
                  bootM(df$KMMDp2),
                  bootM(df$FSp2),
                  bootM(df$KMMDp3),
                  bootM(df$FSp3),
                  bootM(df$KMMDrbf),
                  bootM(df$FSrbf))
  data.frame("N" = df$N[1], "length" = df$length[1], "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = names(df[-(1:2)]))
})

ggplot(res2, aes(x = N, y = value, color = group, linetype = group)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 1)) +
  xlab("N") +
  opts(title = "Power on Image Data")
ggsave("power_image.png")
