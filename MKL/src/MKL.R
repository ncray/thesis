setwd("~/Dropbox/VMshare/thesis/MKL/src")
imgDir <- "../img/"
source("./generate_data.R")
source("./MKL_code.R")
source("../../friedman-test/src/twosample.R")
library(ggplot2)
library(doMC)
registerDoMC(4)
parallel <- TRUE
##http://www.shogun-toolbox.org/doc/en/2.0.1/staticr.html
Npwr <- 200
Nwts <- 100

testShogunKernlab <- function(){
  dat <- ldply(1:20, function(D){
    dat <- getDataNormal(200, D)
    l <- dat$l
    u <- dat$u
    c(D, computeFS(u, kernelMatrix(vanilladot(), x = u), l, 1), compute(trainLinear)(u = t(u), l = l, r = 1, C = 1))
  })
  names(dat) <- c("D", "kernlab", "shogun")
  qplot(x = kernlab, y = shogun, data = dat) + geom_abline(slope = 1)

  dat <- getDataNormal(1000, 100)
  l <- dat$l
  u <- dat$u
  system.time(computeFS(u, kernelMatrix(vanilladot(), x = u), l, 1))
  system.time(compute(trainLinear)(u = t(u), l = l, r = 1, C = 1))

  sig <- 2
  x <- as.numeric(1:5)
  sg('clean_kernel')
  sg('set_features', 'TRAIN', matrix(x, nrow = 1))
  sg('set_kernel', 'GAUSSIAN', 'REAL', 10, sig)
  sg('get_kernel_matrix', 'TRAIN')
  kernelMatrix(rbfdot(sigma = 1 / sig), x)

  dat <- getDataNormal(20, 1)
  l <- dat$l
  u <- dat$u  
  sg('clean_kernel')
  sg('clean_features', 'TRAIN')
  sg('set_features', 'TRAIN', t(u)) ##takes numeric, not integer
  sg('set_kernel', 'GAUSSIAN', 'REAL', 10, 1)
  sg('set_labels', 'TRAIN', as.numeric(as.character(l)))
  sg('new_classifier', 'LIBSVM')
  sg('c', 1)
  sg('svm_use_bias', TRUE) ##default is TRUE
  ##sg('get_kernel_matrix', 'TRAIN')
  sg('train_classifier')
  svmparams <- sg('get_svm')
  b <- as.numeric(svmparams[[1]])
  aw <- svmparams[[2]][, 1]
  inds <- svmparams[[2]][, 2] + 1 ## 0 indexing
  km <- sg('get_kernel_matrix')
  km.sub <- km[, inds]
  ##mar <- as.numeric(km.sub %*% aw + b) ##old
  mar <- as.numeric(km.sub %*% (aw * as.numeric(as.character(l))[inds]) + b)
  ##print(all((mar < 0) == (sg('classify') == -1)))
  ##plot(mar, sg('classify'))
  ##sg('classify')
  mar

  kmRBF <- kernelMatrix(rbfdot(sigma = 1), x = u)
  ksvm.fit <- ksvm(x = kmRBF, y = l, C = 1, shrinking = FALSE, tol = .01)
  y <- as.numeric(as.character(l))
  alpha <- unlist(ksvm.fit@alpha)
  sv.ind <- ksvm.fit@SVindex
  km.sub <- km[, sv.ind]
  b <- ksvm.fit@b
  margins <- km.sub %*% alpha + b

  cor(mar, margins)
  qplot(x = mar, y = as.vector(margins)) + geom_abline(slope = 1)  

  dat <- getDataNormal(20, 1)
  l <- dat$l
  u <- dat$u
  compute(trainRBF)(u = t(u), l = l, r = 1, C = 1)
  kmRBF <- kernelMatrix(rbfdot(sigma = 1), x = u)
  computeFS(u, kmRBF, l, 1)
  
  dat <- mdply(expand.grid(D = 1:5, sig = 1:5), function(D, sig){
    dat <- getDataNormal(200, D)
    l <- dat$l
    u <- dat$u
    c(computeFS(u, kernelMatrix(rbfdot(sigma = 1 / sig), x = u), l, 1), compute(trainRBF)(u = t(u), l = l, r = as.numeric(sig), C = 1))
  })
  names(dat)[3:4] <- c("kernlab", "shogun")
  qplot(x = kernlab, y = shogun, data = dat) + geom_abline(slope = 1)

  dat <- getDataNormal(1000, 10)
  l <- dat$l
  u <- dat$u
  system.time(computeFS(u, kernelMatrix(rbfdot(sigma = 1), x = u), l, 1))
  system.time(compute(trainRBF)(u = t(u), l = l, r = 1, C = 1))
}

nullDistNormal <- function(D = 1, N = 100, C = 1){
  RBF.v <- c(5, 10, 20, 40) ##second run
  ##RBF.v <- 10^(-1:2)
  print(unlist(as.list(environment())))
  dat <- getDataNormal(100, D)
  l <- dat$l
  u <- dat$u

  ldply(1:N, function(x){
    print(x)
    l.perm <- sample(l)
    res <- data.frame("D" = D, "N" = N, "C" = C,
                      "FS-l" = compute(trainLinear)(u = t(u), l = l.perm, r = 1, C = C),
                      "FS-MKL1" = compute(trainMKL)(u1 = t(u), l = l.perm, RBF.v = c(RBF.v), mkl_norm = 1, C = C),
                      "FS-MKL2" = compute(trainMKL)(u1 = t(u), l = l.perm, RBF.v = c(RBF.v), mkl_norm = 2, C = C)
                      )
    res2 <- laply(RBF.v, function(r) compute(trainRBF)(u = t(u), l = l.perm, r = r, C = C))
    res2 <- data.frame(do.call(cbind, llply(RBF.v, function(r) compute(trainRBF)(u = t(u), l = l.perm, r = r, C = C))))
    names(res2) <- paste("FS-rbf", RBF.v, sep = "")
    cbind(res, res2)
  }, .parallel = parallel)
}

getNullDistPlot <- function(){
  system.time(res <- mdply(expand.grid(D = c(1, 5, 10), N = 100, C = c(.1, 1, 10)), nullDistNormal))
  library(nortest)
  res$D <- factor(res$D)
  res.m <- melt(res, id.vars = c("D", "N", "C"))
  pvals <- ddply(res.m, .(C, variable, D), function(df) c("pval" = ad.test(df$value)$p.value))
  arrange(pvals, desc(pval))

  p1 <- ggplot(data = res.m, aes(x = value, fill = D)) +
    geom_density(alpha = .4) +
      facet_grid(C~variable, scales = "free") +
        ggtitle("Null Distributions (Faceted by Statistic)") +
          geom_text(data = pvals, aes(x = 0, y = as.numeric(as.character(D)) / 20, label = round(pval, 5)), color = "red")
  p1
  ##myplot(p1, "null_dist.png")
  myplot(p1, "null_dist2.png")
}

testReject <- function(){
  dat <- getDataNormal(100, D = 10, delta = 1)
  l <- dat$l
  u <- dat$u
  RBF.v <- 10^(-2:2)
  ##RBF.v <- seq(.5, 5, 1)
  trainRBF(u = t(u), l = l, r = .01, C = 1) ## margin is always 1
  mar <- getMargins(l)
  computeT(u = mar, l = l)

  compute(trainLinear)(u = t(u), l = l, r = 1, C = 1)
  compute(trainMKL)(u1 = t(u), l = l, RBF.v = RBF.v, mkl_norm = 1, C = 1)
  compute(trainMKL)(u1 = t(u), l = l, RBF.v = RBF.v, mkl_norm = 2, C = 1)
  laply(RBF.v, function(r) compute(trainRBF)(u = t(u), km = NULL, l = l, r = r, C = 1))
  laply(c(.01, .1, 1, 10, 100), function(C) compute(trainMKL)(u1 = t(u), l = l, RBF.v = RBF.v, mkl_norm = 1, C = C))
  laply(c(.01, .1, 1, 10, 100), function(C) compute(trainMKL)(u1 = t(u), l = l, RBF.v = RBF.v, mkl_norm = 2, C = C))


  reject(compute(trainLinear), verbose = TRUE)(u = t(u), km = NULL, l = l, r = 1, C = 1)
  reject(compute(trainRBF), verbose = TRUE)(u = t(u), km = NULL, l = l, r = 1, C = 1)
  reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = c(RBF.v), mkl_norm = 1, C = 1)
  reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = c(RBF.v), mkl_norm = 2, C = 1)
  reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = c(RBF.v), mkl_norm = 2, C = 10)
  reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = 1, mkl_norm = 1, C = 1)
  reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = 1, mkl_norm = 2, C = 1)
  reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = c(RBF.v[-1]), mkl_norm = 2, C = 10)
}

powerNormal <- function(D = 1, delta = 1, C = 1){
  print(unlist(as.list(environment())))
  ldply(1:Npwr, function(x){
    RBF.v <- c(.5, 1, 5, 10)
    print(x)
    dat <- getDataNormal(20, D, delta)
    l <- dat$l
    u <- dat$u
    kmLinear <- kernelMatrix(vanilladot(), x = u)
    kmRBF <- kernelMatrix(rbfdot(sigma = 1), x = u)

    data.frame("D" = D, "delta" = delta,
               "T2" = reject(computeT2)(u, kmLinear, l, C),
               "KMMD-l" = reject(computeKMMD)(u, kmLinear, l, C),
               "FS-l" = reject(computeFS, parametric = TRUE)(u, kmLinear, l, C),
               "KMMD-rbf" = reject(computeKMMD)(u, kmRBF, l, C),
               "FS-rbf" = reject(computeFS, parametric = TRUE)(u, kmRBF, l, C),
               "FS-MKL1" = reject(compute(trainMKL), parametric = TRUE)(u1 = t(u), l = l, RBF.v = RBF.v, mkl_norm = 1, C = C),
               "FS-MKL2" = reject(compute(trainMKL), parametric = TRUE)(u1 = t(u), l = l, RBF.v = RBF.v, mkl_norm = 2, C = C))
  }, .parallel = parallel)
}
###
power23 <- function(dx1, dx2, C, p){
  print(paste("dx1:", dx1, "dx2:", dx2, "C:", C, "p:", p))
  ##RBF.v <- round(10^((0:6) / 2), 4)
  ##RBF.v <- 1 / round(10^((0:6) / 2), 4)
  RBF.v <- round(10^((-3:3) / 2), 4)
  ldply(1:Npwr, function(x){
    dat <- getData(50, dx1, dx2, p)
    u <- dat$u1
    u2 <- NULL
    l <- dat$l
    print(x)

    df1 <- data.frame("dx1" = dx1, "dx2" = dx2, "p" = p, "C" = C,
                      "T2" = reject(computeT2)(u, kmLinear, l, C),
                      "FS-l" = reject(compute(trainLinear), parametric = TRUE)(u = t(u), l = l, r = 1, C = C),
                      "FS-MKL1" = reject(compute(trainMKL), parametric = TRUE)(u1 = t(u), l = l, RBF.v = RBF.v, mkl_norm = 1, C = C),
                      "FS-MKL2" = reject(compute(trainMKL), parametric = TRUE)(u1 = t(u), l = l, RBF.v = RBF.v, mkl_norm = 2, C = C))
    df2 <- as.data.frame(matrix(laply(RBF.v, function(r) reject(compute(trainRBF), parametric = TRUE)
                 (u = t(u), l = l, r = r, C = C)), nrow = 1))
    cbind(df1, df2)
  }, .parallel = parallel)
}

test <- function(){
Npwr <- 50
##system.time(res <- mdply(expand.grid(dx1 = 1, dx2 = 10, C = .01, p = c(0, .5, 1)), power23))
system.time(res <- mdply(expand.grid(dx1 = 100, dx2 = 1000, C = .1, p = c(0, .5, 1)), power23))
res2 <- ddply(res, .(dx1, dx2, p, C), function(df){
  ldply(names(df)[-(1:4)], function(name){
    dat <- df[, name]
    lims <- c(.025, .975)
    bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
    boot <- bootM(dat)
    data.frame("value" = mean(dat), "lower" = boot[1], "upper" = boot[2], "group" = name)
  })
})

p4 <- ggplot(res2, aes(x = group, y = value)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .05)) +
  facet_grid(p~.)
p4
}
###
getPowerNormal <- function(){
  Npwr <- 1000
  system.time(res <- mdply(expand.grid("delta" = seq(0, 1.5, .5), "D" = c(1, 5, 10, 20), "C" = c(.1, 1, 10)),
                           powerNormal))

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
  p2
  myplot(p2, "normal_power.png")
}


plotStar <- function(){
  star <- generateStar(r1 = 10, r2 = 4, n = 1000)
  star.df <- data.frame(t(do.call(rbind, star)))
  names(star.df) <- c("x", "y", "label")
  star.df$label <- factor(star.df$label)
  qplot(x, y, data = star.df, geom = "point", color = label)
}

starTest <- function(){
  dat <- generateStar(r1 = 6, r2 = 4, n = 50)
  compute(trainMKL)(u1 = dat$u, l = dat$l, RBF.v = c(.1, 1, 10), mkl_norm = 2, C = 1)

  system.time(res <- laply(1:1000, function(i) compute(trainMKL)(u1 = dat$u, l = sample(dat$l), RBF.v = c(.1, 1, 10), mkl_norm = 2, C = .01), .parallel = TRUE))
  ad.test(res)

  RBF.v <- 10^(-2:2)
  foo <- function(C) compute(trainMKL)(u1 = dat$u, l = sample(dat$l), RBF.v = RBF.v, mkl_norm = 2, C = C)
  res <- laply(1:1000, function(i) laply(c(.01, .1, 1, 10), foo), .parallel = TRUE)
  res2 <- aaply(res, 1, max)
  ad.test(res2)
  qqnorm(res2); qqline(res2)

  trainMKL(u1 = dat$u, l = dat$l, RBF.v = c(.1, 1, 10, 100, 200))
  getMKLWeights()
}

MKLwtsStar <- function(r1, n, C){
  print(unlist(as.list(environment())))
  dat <- generateStar(r1 = r1, r2 = 4, n = n)
  trainMKL(u1 = dat$u, l = dat$l, RBF.v = RBF.v, C, mkl_norm = 1, linear = FALSE)
  wts <- getMKLWeights()
  df1 <- cbind(data.frame("r1" = r1, "n" = n, "C" = C, "perm" = 0, "mkl_norm" = 1), matrix(wts, nrow = 1))
  trainMKL(u1 = dat$u, l = dat$l, RBF.v = RBF.v, C, mkl_norm = 2, linear = FALSE)
  wts <- getMKLWeights()
  wts <- wts / sum(wts)
  df2 <- cbind(data.frame("r1" = r1, "n" = n, "C" = C, "perm" = 0, "mkl_norm" = 2), matrix(wts, nrow = 1))
  dfperm <- ldply(1:Nwts, function(x){
    trainMKL(u1 = dat$u, l = sample(dat$l), RBF.v = RBF.v, C, mkl_norm = 1, linear = FALSE)
    wts <- getMKLWeights()
    df1p <- cbind(data.frame("r1" = r1, "n" = n, "C" = C, "perm" = 1, "mkl_norm" = 1), matrix(wts, nrow = 1))
    trainMKL(u1 = dat$u, l = sample(dat$l), RBF.v = RBF.v, C, mkl_norm = 2, linear = FALSE)
    wts <- getMKLWeights()
    wts <- wts / sum(wts)
    df2p <- cbind(data.frame("r1" = r1, "n" = n, "C" = C, "perm" = 1, "mkl_norm" = 2), matrix(wts, nrow = 1))
    rbind(df1p, df2p)
  })
  df1 <- rbind(df1, df2, dfperm)
  len <- ncol(df1)
  names(df1)[(len - length(RBF.v) + 1):len] <- paste("rbf: ", RBF.v, sep = "")
  df1
}

MKLwtShiftStar <- function(){
  Nwts <- 100
  RBF.v <- round(10^(seq(.5, 2, .5)), 2)
  ##system.time(res <- mdply(expand.grid(r1 = round(10^(seq(-.5, 2, .5)) / 3, 2), n = 100, C = 1), MKLwtsStar, .parallel = TRUE))
  system.time(res <- mdply(expand.grid(r1 = c(1, 3, 7, 9, 11, 13, 15, 20), n = 100, C = 1), MKLwtsStar, .parallel = TRUE))

  res.m <- melt(res, id.vars = c(1:5))
  p1 <- qplot(variable, value, data = subset(res.m, perm == 1), geom = "boxplot") +
    facet_grid(mkl_norm~r1) +
      geom_point(data = subset(res.m, perm == 0), color = "red", size = 3) +
        xlab("Kernels") +
          ylab("Kernel Weights") +
            ggtitle("Boxplot of Null Distribution with Observed in Red Faceted by MKL Norm and Outer Radius")
  p1
  myplot(p1, "mkl_weights_star.png")
}

powerStar <- function(r1, n, C){
  print(unlist(as.list(environment())))
  RBF.v <- round(10^(seq(.5, 2, .5)), 2)
  ldply(1:Npwr, function(x){
    print(x)
    dat <- generateStar(r1 = r1, r2 = 4, n = n)
    l <- dat$l
    u <- dat$u
    kmLinear <- kernelMatrix(vanilladot(), x = t(u))
    kmRBF <- kernelMatrix(rbfdot(sigma = 1), x = t(u))
    data.frame("r1" = r1, "n" = n, "C" = C, 
               "T2" = reject(computeT2)(t(u), kmLinear, l, C),
               "KMMD-l" = reject(computeKMMD)(u, kmLinear, l, C),
               "FS-l" = reject(computeFS, parametric = TRUE)(u, kmLinear, l, C),
               "KMMD-rbf" = reject(computeKMMD)(u, kmRBF, l, C),
               "FS-rbf" = reject(computeFS, parametric = TRUE)(u, kmRBF, l, C),
               "FS-MKL1" = reject(compute(trainMKL), parametric = TRUE)(u1 = u, l = l, RBF.v = RBF.v, mkl_norm = 1, C = C, linear = FALSE),
               "FS-MKL2" = reject(compute(trainMKL), parametric = TRUE)(u1 = u, l = l, RBF.v = RBF.v, mkl_norm = 2, C = C, linear = FALSE))
  }, .parallel = parallel)
}

powerStarPlot <- function(){
  Npwr <- 1000
  system.time(res <- mdply(expand.grid(r1 = seq(4, 7, .5), n = 50, C = 1), powerStar))

  res2 <- ddply(res, .(r1, n, C), function(df){
    ldply(names(df)[-(1:3)], function(name){
      dat <- df[, name]
      lims <- c(.025, .975)
      bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
      boot <- bootM(dat)
      data.frame("value" = mean(dat), "lower" = boot[1], "upper" = boot[2], "group" = name)
    })
  })

  p2 <- ggplot(res2, aes(x = r1, y = value, color = group, linetype = group)) +
    geom_line() + 
      geom_errorbar(aes(ymin = lower, ymax = upper, width = .05)) +
        ##facet_grid(self~.) +
        ggtitle("Power (Christmas Star Example)") +
          xlab("Radius of Outer Star (Inner is 4)") +
            ylab("Power")
  p2
  myplot(p2, "star_power.png")
}


MKLwtsDNAStar <- function(r1, self, n, C){
  print(unlist(as.list(environment())))
  RBF.v <- round(10^(seq(.5, 2, .5)), 2)
  string.v <- 1:3
  dat <- getDataDNAStar(r1 = r1, self = self, n = n)
  trainMKL(u1 = dat$u1, u2 = dat$u2, l = dat$l, RBF.v = RBF.v, string.v = string.v, C, mkl_norm = 1, linear = FALSE)
  wts <- getMKLWeights()
  df1 <- cbind(data.frame("r1" = r1, "n" = n, "C" = C, "perm" = 0, "mkl_norm" = 1), matrix(wts, nrow = 1))
  trainMKL(u1 = dat$u1, u2 = dat$u2, l = dat$l, RBF.v = RBF.v, string.v = string.v, C, mkl_norm = 2, linear = FALSE)
  wts <- getMKLWeights()
  wts <- wts / sum(wts)
  df2 <- cbind(data.frame("r1" = r1, "n" = n, "C" = C, "perm" = 0, "mkl_norm" = 2), matrix(wts, nrow = 1))
  dfperm <- ldply(1:Nwts, function(x){
    trainMKL(u1 = dat$u1, u2 = dat$u2, l = sample(dat$l), RBF.v = RBF.v, string.v = string.v, C, mkl_norm = 1, linear = FALSE)
    wts <- getMKLWeights()
    df1p <- cbind(data.frame("r1" = r1, "n" = n, "C" = C, "perm" = 1, "mkl_norm" = 1), matrix(wts, nrow = 1))
    trainMKL(u1 = dat$u1, u2 = dat$u2, l = sample(dat$l), RBF.v = RBF.v, string.v = string.v, C, mkl_norm = 2, linear = FALSE)
    wts <- getMKLWeights()
    wts <- wts / sum(wts)
    df2p <- cbind(data.frame("r1" = r1, "n" = n, "C" = C, "perm" = 1, "mkl_norm" = 2), matrix(wts, nrow = 1))
    rbind(df1p, df2p)
  })
  df1 <- rbind(df1, df2, dfperm)
  len <- ncol(df1)
  names(df1)[6:(6 + length(RBF.v) - 1)] <- paste("rbf: ", RBF.v, sep = "")
  names(df1)[(6 + length(RBF.v)):len] <- paste("sk: ", string.v, sep = "")
  df1
}

MKLwtsDNAStarPlot <- function(){
  Nwts <- 100
  system.time(res <- mdply(expand.grid(r1 = 4.5, self = seq(.25, .4, .05), n = 200, C = .1), MKLwtsDNAStar))

  res.m <- melt(res, id.vars = c(1:6))
  p1 <- qplot(variable, value, data = subset(res.m, perm == 1), geom = "boxplot") +
    facet_grid(mkl_norm~self) +
      geom_point(data = subset(res.m, perm == 0), color = "red", size = 3) +
        xlab("Kernels") +
          ylab("Kernel Weights") +
            ggtitle("Boxplot of Null Distribution with Observed in Red Faceted by Self Transition Probability and MKL Norm")
  p1
  myplot(p1, "mkl_weights_star_dna.png")
}

powerDNAStar <- function(r1, self, n, C){
  print(unlist(as.list(environment())))
  RBF.v <- round(10^(seq(.5, 2, .5)), 2)
  string.v <- 1:3
  ldply(1:Npwr, function(x){
    print(x)
    dat <- getDataDNAStar(r1 = r1, self = self, n = n)
    l <- dat$l
    u1 <- dat$u1
    u2 <- dat$u2

    dfMKL <- data.frame("r1" = r1, "self" = self, "n" = n, "C" = C,
                        "FSMKL: 1" = reject(compute(trainMKL), parametric = TRUE)
                        (u1 = u1, u2 = u2, l = l, RBF.v = RBF.v, string.v = string.v, mkl_norm = 1, C = C, linear = FALSE),
                        "FSMKL: 2" = reject(compute(trainMKL), parametric = TRUE)
                        (u1 = u1, u2 = u2, l = l, RBF.v = RBF.v, string.v = string.v, mkl_norm = 2, C = C, linear = FALSE))
    dfRBF <- as.data.frame(matrix(laply(RBF.v, function(r) reject(compute(trainRBF), parametric = TRUE)
                                        (u = u1, l = l, r = r, C = C)), nrow = 1))
    names(dfRBF) <- paste("RBF: ", RBF.v, sep = "")
    dfSK <- as.data.frame(matrix(laply(string.v, function(order) reject(compute(trainString), parametric = TRUE)
                                       (u = u2, l = l, order = order, C = C)), nrow = 1))
    names(dfSK) <- paste("SK: ", string.v, sep = "")
    cbind(dfMKL, dfRBF, dfSK)
  }, .parallel = parallel)
}

test <- function(){
  C <- .1
  r1 <- 4.3
  self <- .35
  n <- 75
  dat <- getDataDNAStar(r1 = r1, self = self, n = n)
  l <- dat$l
  u1 <- dat$u1
  u2 <- dat$u2
  RBF.v <- 100
  string.v <- 3

  ##change kernel normalization from SQRTDIAG to IDENTITY
  reject(compute(trainMKL), parametric = TRUE)(u1 = u1, u2 = u2, l = l, RBF.v = RBF.v, string.v = string.v, mkl_norm = 1, C = C, linear = FALSE)
  getMKLWeights()
  reject(compute(trainMKL), parametric = TRUE)(u1 = u1, u2 = u2, l = l, RBF.v = RBF.v, string.v = string.v, mkl_norm = 2, C = C, linear = FALSE)
  getMKLWeights()
  
  reject(compute(trainMKL), parametric = FALSE, verbose = TRUE)(u1 = u1, u2 = u2, l = l, RBF.v = RBF.v, string.v = string.v, mkl_norm = 1, C = C, linear = FALSE)
  reject(compute(trainMKL), parametric = FALSE, verbose = TRUE)(u1 = u1, u2 = u2, l = l, RBF.v = RBF.v, string.v = string.v, mkl_norm = 2, C = C, linear = FALSE)

  reject(compute(trainRBF), parametric = TRUE)(u = u1, l = l, r = RBF.v, C = C)
  reject(compute(trainRBF), parametric = FALSE, verbose = TRUE)(u = u1, l = l, r = RBF.v, C = C)

  reject(compute(trainString), parametric = TRUE)(u = u2, l = l, order = string.v, C = C)
  reject(compute(trainString), parametric = FALSE, verbose = TRUE)(u = u2, l = l, order = string.v, C = C)
}

powerDNAStarPlot <- function(){
  Npwr <- 1000
  ##system.time(res <- mdply(expand.grid(r1 = c(4, 4.3), self = c(.25, .35, .45), n = 50, C = .1), powerDNAStar))
  system.time(res <- mdply(expand.grid(r1 = c(4, 4.3, 4.6), self = c(.25, .35, .45), n = 50, C = .1), powerDNAStar))

  ##system.time(res <- mdply(expand.grid(r1 = c(4.3), self = c(.335), n = 50, C = .1), powerDNAStar))
  ##colMeans(res)
  
  res2 <- ddply(res, .(r1, self, n, C), function(df){
    ldply(names(df)[-(1:4)], function(name){
      dat <- df[, name]
      lims <- c(.025, .975)
      bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
      boot <- bootM(dat)
      data.frame("value" = mean(dat), "lower" = boot[1], "upper" = boot[2], "group" = name)
    })
  })

  p2 <- ggplot(res2, aes(x = self, y = value, color = group, linetype = group)) +
    geom_line() + 
      geom_errorbar(aes(ymin = lower, ymax = upper, width = .005)) +
        facet_grid(r1~.) +
          ggtitle("Power (Christmas Star + DNA Example), Faceted on Outer Radius") +
            xlab("Self Transition Probability") +
              ylab("Power")
  p2
  myplot(p2, "dna_star_power.png")
}




getNullDistPlot()
getPowerNormal()
RBF.v <- round(10^(seq(.5, 2, .5)), 2)
MKLwtShiftStar()
powerStarPlot()
MKLwtsDNAStarPlot()
powerDNAStarPlot()
