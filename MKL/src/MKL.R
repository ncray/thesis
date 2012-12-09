setwd("~/Dropbox/VMshare/thesis/MKL/src")
imgDir <- "../img/"
source("./generate_data.R")
source("./MKL_code.R")
source("../../friedman-test/src/twosample.R")
library(ggplot2)
##http://www.shogun-toolbox.org/doc/en/2.0.1/staticr.html

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








star <- generateStar(r1 = 10, r2 = 4, n = 1000)
star.df <- data.frame(t(do.call(rbind, star)))
names(star.df) <- c("x", "y", "label")
star.df$label <- factor(star.df$label)
qplot(x, y, data = star.df, geom = "point", color = label)

dat <- generateStar(r1 = 10, r2 = 4, n = 20)
trainMKL(u1 = dat$u, l = dat$l, RBF.v = c(.1, 1, 10, 100, 200))
getMKLWeights()
