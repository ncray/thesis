source("./twosample.R")
library(ggplot2)
library(reshape)
library(boot)
library(doMC)
registerDoMC(4)
parallel <- TRUE
Npwr <- 50

getData <- function(N, D = 1, delta = 1){
  u <- rbind(matrix(rnorm(N * D, 0), ncol = D), matrix(rnorm(N * D, delta), ncol = D))
  l <- factor(c(rep(-1, N), rep(1, N)))
  list("u" = u, "l" = l)
}

checkComputeReject <- function(){
  dat <- getData(10, 2)
  u <- dat$u
  km <- kernelMatrix(vanilladot(), x = u)
  l <- dat$l

  computeT(u, km, l)^2
  computeT2(u, km, l)
  computeKMMD(u, km, l)
  computeFS(u, km, l)

  reject(computeT2, TRUE)(u, km, l)
  reject(computeKMMD, TRUE)(u, km, l)
  reject(computeFS, TRUE)(u, km, l)
}

nullDist <- function(D = 1, N = 100, C = 1){
  print(unlist(as.list(environment())))
  dat <- getData(200, D)
  l <- dat$l
  u <- dat$u
  kmLinear <- kernelMatrix(vanilladot(), x = u)
  kmRBF <- kernelMatrix(rbfdot(sigma = 1), x = u)
  ldply(1:N, function(x){
    print(x)
    l.perm <- sample(l)
    data.frame("D" = D, "N" = N,
               "T2" = computeT2(u, kmLinear, l.perm, C),
               "sqrtT2" = sqrt(computeT2(u, kmLinear, l.perm, C)),
               "KMMD-l" = computeKMMD(u, kmLinear, l.perm, C),
               "FS-l" = computeFS(u, kmLinear, l.perm, C),
               "KMMD-rbf" = computeKMMD(u, kmRBF, l.perm, C),
               "FS-rbf" = computeFS(u, kmRBF, l.perm, C))
  }, .parallel = parallel)
}

nullDistSim <- function(){
  res <- mdply(expand.grid(D = c(1, 5, 10)), nullDist)
  res$D <- factor(res$D)
  res.m <- melt(res, id.vars = c("D", "N"))

  p1 <- ggplot(data = res.m, aes(x = value, fill = D)) +
    geom_density(alpha = .4) +
      facet_wrap(~variable, scales = "free") +
        opts(title = "Null Distributions (Faceted by Statistic)")
  png("null_dist.png", width = 800, height = 600)
  print(p1)
  dev.off()
}

powerMultivariate <- function(D = 1, delta = 1, C = 1){
  print(unlist(as.list(environment())))
  ldply(1:Npwr, function(x){
    print(x)
    dat <- getData(20, D, delta)
    l <- dat$l
    u <- dat$u
    kmLinear <- kernelMatrix(vanilladot(), x = u)
    kmRBF <- kernelMatrix(rbfdot(sigma = 1), x = u)
    data.frame("D" = D, "delta" = delta,
               "T2" = reject(computeT2)(u, kmLinear, l, C),
               "KMMD-l" = reject(computeKMMD)(u, kmLinear, l, C),
               "FS-l" = reject(computeFS)(u, kmLinear, l, C),
               "KMMD-rbf" = reject(computeKMMD)(u, kmRBF, l, C),
               "FS-rbf" = reject(computeFS)(u, kmRBF, l, C))
  }, .parallel = parallel)
}

powerSim <- function(){
  system.time(res <- mdply(expand.grid("delta" = seq(0, 1.5, .25), "D" = c(1, 5, 10, 20)),
                           powerMultivariate))

  res2 <- ddply(res, .(delta, D), function(df){
    ldply(names(df)[-(1:2)], function(name){
      dat <- df[, name]
      lims <- c(.025, .975)
      bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
      boot <- bootM(dat)
      data.frame("value" = mean(dat), "lower" = boot[1], "upper" = boot[2], "group" = name)
    })
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
}

twitter <- function(){
  source("./twitter.R")
  obama <- getObama()
  palin <- getPalin()
  getTextSamp <- function(obama, palin, n){
    list(u = c(obama[sample(1:length(obama), n)], palin[sample(1:length(palin), n)]),
         l = factor(c(rep(-1, n), rep(1, n))))
  }

  powerTwitter <- function(N = 10, len = 3, C = 1){
    print(unlist(as.list(environment())))
    ldply(1:Npwr, function(x){  
      dat <- getTextSamp(obama, palin, N)
      l <- dat$l
      u <- dat$u
      kmString <- kernelMatrix(stringdot(length = len), x = u)
      print(x)
      data.frame("N" = N, "length" = len,
                 "KMMD" = reject(computeKMMD)(u, kmString, l, C),
                 "FS" = reject(computeFS)(u, kmString, l, C))
    }, .parallel = parallel)
  }

  system.time(res <- mdply(expand.grid("N" = c(10, 15, 20, 30, 40), "len" = 3), powerTwitter))
  ddply(res, .(N, length), colMeans)
}

birds <- function(){
  library(jpeg)
  library(plyr)
  rooster.path <- list.files("../img/rooster_resized", full.names = TRUE)[1:45]
  pigeon.path <- list.files("../img/pigeon_resized", full.names = TRUE)[1:45]

  rooster <- laply(rooster.path, function(path) as.vector(readJPEG(path)))
  pigeon <- laply(pigeon.path, function(path) as.vector(readJPEG(path)))

  getImageSamp <- function(rooster, pigeon, N){
    u <- rbind(rooster[sample(1:nrow(rooster), N), ], pigeon[sample(1:nrow(pigeon), N), ])
    u <- sweep(u, 1, apply(u, 1, mean), "-")
    u <- sweep(u, 1, apply(u, 1, function(vec) sqrt(sum(vec^2))), "/")
    l <- factor(c(rep(-1, N), rep(1, N)))
    list(u = u, l = l)
  }

  powerBirds <- function(N = 10, C = 1){
    print(unlist(as.list(environment())))
    ldply(1:Npwr, function(x){
      print(x)
      dat <- getImageSamp(rooster, pigeon, N)
      l <- dat$l
      u <- dat$u
      kmp1 <- kernelMatrix(polydot(degree = 1, offset = 1), x = u)
      kmp2 <- kernelMatrix(polydot(degree = 2, offset = 1), x = u)
      kmp3 <- kernelMatrix(polydot(degree = 3, offset = 1), x = u)
      data.frame("N" = N, "length" = NA,
                 "KMMDp1" = reject(computeFS)(u, kmp1, l, C),
                 "FSp1" = reject(computeFS)(u, kmp1, l, C),
                 "KMMDp2" = reject(computeFS)(u, kmp2, l, C),
                 "FSp2" = reject(computeFS)(u, kmp2, l, C),
                 "KMMDp3" = reject(computeFS)(u, kmp3, l, C),
                 "FSp3" = reject(computeFS)(u, kmp3, l, C))
    }, .parallel = parallel)
  }
}
