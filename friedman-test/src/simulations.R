setwd("~/Dropbox/VMshare/thesis/friedman-test/src")
imgDir <- "../img/"
source("./twosample.R")
library(tikzDevice)
library(ggplot2)
library(reshape)
library(boot)
library(doMC)
registerDoMC(4)
parallel <- TRUE
Npwr <- 50

checkComputeReject <- function(){
  dat <- getDataNormal(N = 10, D = 10)
  ##dat <- getDataNormal(10, 1)
  u <- dat$u
  km <- kernelMatrix(vanilladot(), x = u)
  kmRBF <- kernelMatrix(rbfdot(sigma = 1), x = u)
  l <- dat$l

  computeT(u, km, l)^2
  computeT2(u, km, l)
  computeKMMD(u, km, l)
  computeFS(u, km, l)
  computeFS(u, km, l)^2

  laply(10^seq(-3, 3, 1), function(C) computeFS(u, km, l, C)) ##in 1 dimension, FS is independent of C, but not in > 1 dim

  reject(computeT2, TRUE)(u, km, l)
  reject(computeKMMD, TRUE)(u, km, l)
  reject(computeFS, TRUE)(u, km, l)

  set.seed(1)
  dat <- ldply(1:500, function(i){
    dat <- getDataNormal(10, 1)
    u <- dat$u
    km <- kernelMatrix(vanilladot(), x = u)
    l <- dat$l
    data.frame("t" = computeT(u, km, l), "fs" = computeFS(u, km, l))
  })
  which(abs(dat$t) - abs(dat$fs) > 1e-5)
  plot(dat$t, dat$fs)
}

regularizationDimension <- function(D){
  dat <- getDataNormal(N = 50, D = D)
  u <- dat$u
  km <- kernelMatrix(vanilladot(), x = u)
  kmRBF <- kernelMatrix(rbfdot(sigma = 1), x = u)
  l <- sample(dat$l)

  ldply(10^(seq(-3, 3, .25)), function(C) c("D" = D,
                       "C" = C,
                       "FS.l" = computeFS(u, km, l, C),
                       "FS.rbf" = computeFS(u, kmRBF, l, C)))
}

regularizationDimensionPlot <- function(){
  system.time(dat <- ldply(c(1, 2, 5, 10), regularizationDimension, .parallel = TRUE))
  dat.m <- melt(dat, id.vars = c("D", "C"))
  p1 <- ggplot(data = dat.m, aes(x = C, y = value, color = variable)) +
    geom_point() +
      geom_line() + 
        facet_grid(D~., scales = "free") +
          scale_x_log10() +
            ggtitle("Effect of C on t-statistic") +
              ylab("t-statistic")
  myTikz("c_param.tex", p1)
}

nullDist <- function(D = 1, N = 100, C = 1){
  print(unlist(as.list(environment())))
  dat <- getDataNormal(200, D)
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
  system.time(res <- mdply(expand.grid(D = c(1, 5, 10), N = 5000), nullDist))
  res$D <- factor(res$D)
  res.m <- melt(res, id.vars = c("D", "N"))

  library(nortest)
  pvals <- ddply(res.m, .(D, variable), function(df) ad.test(df$value)$p.value)
  subset(pvals, V1 > .001)
  ##ks.test(subset(res.m, D == 1 & variable == "FS.l")$value, "pt", 1)
  ##qplot(x = sort(subset(res.m, D == 1 & variable == "FS.l")$value), y = qt(p = seq(1/2001, 1-1/2001, 1 / 2001), df = 1))

  p1 <- ggplot(data = res.m, aes(x = value, fill = D)) +
    geom_density(alpha = .4) +
      facet_wrap(~variable, scales = "free") +
        ggtitle("Null Distributions (Faceted by Statistic)")
  myTikz("null_dist.tex", p1)
  myplot(p1, "null_dist.png")
}

##system.time(nullDistSim())

powerMultivariate <- function(D = 1, delta = 1, C = 1){
  print(unlist(as.list(environment())))
  ldply(1:Npwr, function(x){
    print(x)
    dat <- getDataNormal(20, D, delta)
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
  system.time(res <- mdply(expand.grid("delta" = seq(0, 1.5, .5), "D" = c(1, 5, 10, 20), "C" = c(.1, 1, 10)),
                           powerMultivariate)) ##700s for 4 cores, Npwr = 50

  res2 <- ddply(res, .(delta, D, C), function(df){
    ldply(names(df)[-(1:3)], function(name){
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
          facet_grid(C~D) +
            ggtitle("Power (Faceted by Dimension and C)")
  myTikz("power_normal.tex", p2)
  myplot(p2, "power_normal.png")
}

##powerSim()

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

  system.time(res <- mdply(expand.grid("N" = c(10, 15, 20, 30, 40), "len" = 1:3, "C" = c(.1, 1, 10)), powerTwitter))
  ddply(res, .(N, length), colMeans)

  res2 <- ddply(res, .(N, len, C), function(df){
    ldply(names(df)[-(1:4)], function(name){
      dat <- df[, name]
      lims <- c(.025, .975)
      bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
      boot <- bootM(dat)
      data.frame("value" = mean(dat), "lower" = boot[1], "upper" = boot[2], "group" = name)
    })
  })

  p2 <- ggplot(res2, aes(x = N, y = value, color = group, linetype = group)) +
    geom_line(size = 1.5) + 
      geom_errorbar(aes(ymin = lower, ymax = upper, width = 2), size = 1.5) +
        xlab("$n$") +
          facet_grid(C~len) +
            ggtitle("Power (Faceted by K and C)")
  myTikz("power_string.tex", p2)
  myplot(p2, "power_string.png")
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

  powerBirds <- function(N = 10, C = 1, deg = 3, offset = 1){
    print(unlist(as.list(environment())))
    ldply(1:Npwr, function(x){
      print(x)
      dat <- getImageSamp(rooster, pigeon, N)
      l <- dat$l
      u <- dat$u
      kmp <- kernelMatrix(polydot(degree = deg, offset = offset), x = u)
      data.frame("N" = N, "C" = C, "deg" = deg,
                 "KMMD" = reject(computeKMMD)(u, kmp, l, C),
                 "FS" = reject(computeFS)(u, kmp, l, C))
    }, .parallel = parallel)
  }

  system.time(res <- mdply(expand.grid("N" = seq(5, 30, 5), "C" = c(.1, 1, 10), "deg" = 1:4), powerBirds))
  res2 <- ddply(res, .(N, C, deg), function(df){
    ldply(names(df)[-(1:3)], function(name){
      dat <- df[, name]
      lims <- c(.025, .975)
      bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
      boot <- bootM(dat)
      data.frame("value" = mean(dat), "lower" = boot[1], "upper" = boot[2], "group" = name)
    })
  })

  p3 <- ggplot(res2, aes(x = N, y = value, color = group)) +
    geom_line(size = 1.5) + 
      geom_errorbar(aes(ymin = lower, ymax = upper, width = 2), size = 1.5) +
        ylab("power") +
          xlab("n") +
            facet_grid(C~deg) +
              ggtitle("Power (Faceted by C and Degree)")
  myTikz("power_birds.tex", p3)
  myplot(p3, "power_birds.png")
}

twitter()
birds()
