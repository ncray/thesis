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
    c(D, computeFS(u, kernelMatrix(vanilladot(), x = u), l, 1), computeFSLinear(u = t(u), l = l, r = 1, C = 1))
  })
  names(dat) <- c("D", "kernlab", "shogun")
  qplot(x = kernlab, y = shogun, data = dat) + geom_abline(slope = 1)

  dat <- getDataNormal(1000, 100)
  l <- dat$l
  u <- dat$u
  system.time(computeFS(u, kernelMatrix(vanilladot(), x = u), l, 1))
  system.time(computeFSLinear(u = t(u), l = l, r = 1, C = 1))

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
  mar <- as.numeric(km.sub %*% (aw * y[inds]) + b)
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

  dat <- getDataNormal(20, 1)
  l <- dat$l
  u <- dat$u
  computeFSRBF(u = t(u), l = l, r = 1, C = 1)
  kmRBF <- kernelMatrix(rbfdot(sigma = 1), x = u)
  computeFS(u, kmRBF, l, 1)
  
  dat <- mdply(expand.grid(D = 1:5, sig = 1:5), function(D, sig){
    dat <- getDataNormal(200, D)
    l <- dat$l
    u <- dat$u
    c(computeFS(u, kernelMatrix(rbfdot(sigma = 1 / sig), x = u), l, 1), computeFSRBF(u = t(u), l = l, r = as.numeric(sig), C = 1))
  })
  names(dat)[3:4] <- c("kernlab", "shogun")
  qplot(x = kernlab, y = shogun, data = dat) + geom_abline(slope = 1)

  dat <- getDataNormal(1000, 10)
  l <- dat$l
  u <- dat$u
  system.time(computeFS(u, kernelMatrix(rbfdot(sigma = 1), x = u), l, 1))
  system.time(computeFSRBF(u = t(u), l = l, r = 1, C = 1))
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
               "T2" = computeT2(u, kmLinear = NULL, l.perm, C),
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

system.time(res <- mdply(expand.grid(D = c(1, 5, 10), N = 10, C = c(.1, 1, 10)), nullDistNormal))
res$D <- factor(res$D)
res.m <- melt(res, id.vars = c("D", "N", "C"))

p1 <- ggplot(data = res.m, aes(x = value, fill = D)) +
  geom_density(alpha = .4) +
  facet_grid(C~variable, scales = "free") +
  opts(title = "Null Distributions (Faceted by Statistic)")
p1


compute(trainLinear)(u = t(u), l = l, r = 1, C = 1)
compute(trainMKL)(u1 = t(u), l = l, RBF.v = c(RBF.v), mkl_norm = 1, C = 1)

reject(compute(trainLinear), verbose = TRUE)(u = t(u), km = NULL, l = l, r = 1, C = 1)
reject(compute(trainRBF), verbose = TRUE)(u = t(u), km = NULL, l = l, r = 1, C = 1)
reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = c(RBF.v), mkl_norm = 1, C = 1)
reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = c(RBF.v), mkl_norm = 2, C = 1)
reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = c(RBF.v), mkl_norm = 2, C = 10)
reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = 1, mkl_norm = 1, C = 1)
reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = 1, mkl_norm = 2, C = 1)
reject(compute(trainMKL), verbose = TRUE)(u1 = t(u), l = l, RBF.v = c(RBF.v[-1]), mkl_norm = 2, C = 10)


powerNormal <- function(D = 1, delta = 1, C = 1){
  print(unlist(as.list(environment())))
  ldply(1:Npwr, function(x){
    print(x)
    dat <- getDataNormal(20, D, delta)
    l <- dat$l
    u <- dat$u

    data.frame("D" = D, "delta" = delta,
               "T2" = reject(computeT2)(u, kmLinear, l, C),
               "KMMD-l" = reject(computeKMMD)(u, kmLinear, l, C),
               "FS-l" = reject(computeFS)(u, kmLinear, l, C),
               "KMMD-rbf" = reject(computeKMMD)(u, kmRBF, l, C),
               "FS-rbf" = reject(computeFS)(u, kmRBF, l, C))
  }, .parallel = parallel)
}



star <- generateStar(r1 = 10, r2 = 4, n = 1000)
star.df <- data.frame(t(do.call(rbind, star)))
names(star.df) <- c("x", "y", "label")
star.df$label <- factor(star.df$label)
qplot(x, y, data = star.df, geom = "point", color = label)

dat <- generateStar(r1 = 10, r2 = 4, n = 20)
trainMKL(u1 = dat$u, l = dat$l, RBF.v = c(.1, 1, 10, 100, 200))
getMKLWeights()
