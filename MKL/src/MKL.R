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
  RBF.v <- 10^(-1:2)
  print(unlist(as.list(environment())))
  dat <- getDataNormal(200, D)
  l <- dat$l
  u <- dat$u

  ldply(1:N, function(x){
    print(x)
    l.perm <- sample(l)
    res <- data.frame("D" = D, "N" = N, "C" = C,
                      "T2" = computeT2(u, km = NULL, l.perm, C),
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

getNullDist <- function(){
  system.time(res <- mdply(expand.grid(D = c(1, 5, 10), N = 100, C = c(.1, 1, 10)), nullDistNormal))
  library(nortest)
  apply(res[, -(1:3)], 2, ad.test)
  res$D <- factor(res$D)
  res.m <- melt(res, id.vars = c("D", "N", "C"))

  p1 <- ggplot(data = res.m, aes(x = value, fill = D)) +
    geom_density(alpha = .4) +
      facet_grid(C~variable, scales = "free") +
        ggtitle("Null Distributions (Faceted by Statistic)")
  p1
}

testReject <- function(){
  dat <- getDataNormal(50, 10, delta = 10)
  l <- dat$l
  u <- dat$u
  RBF.v <- 10^(-2:2)
  RBF.v <- seq(.5, 5, 1)

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
    RBF.v <- c(.5, 1, 5)
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

getPowerNormal <- function(){
  Npwr <- 20
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
  Nwts <- 50
  RBF.v <- round(10^(seq(.5, 2, .5)), 2)
  system.time(res <- mdply(expand.grid(r1 = round(10^(seq(-.5, 2, .5)) / 3, 2), n = 100, C = 1), MKLwtsStar, .parallel = TRUE))
  system.time(res <- mdply(expand.grid(r1 = c(1, 3, 7, 9, 11, 13, 15, 20), n = 100, C = 1), MKLwtsStar, .parallel = TRUE))

  res.m <- melt(res, id.vars = c(1:5))
  p1 <- qplot(variable, value, data = subset(res.m, perm == 1), geom = "boxplot") +
    facet_grid(mkl_norm~r1) +
      geom_point(data = subset(res.m, perm == 0), color = "red", size = 3) +
        xlab("Kernels") +
          ylab("Kernel Weights") +
            ggtitle("Boxplot of Null Distribution with Observed in Red Faceted by Self Transition Probability and Outer Radius")
  p1
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
  Npwr <- 200
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
}


