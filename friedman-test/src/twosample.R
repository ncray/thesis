## Functions directly related to two-sample tests
library(kernlab)
library(ICSNP)
library(plyr)

computeT <- function(u, km, l, ...) t.test(u[l == 1], u[l == -1], var.equal = TRUE)$statistic

computeT2 <- function(u, km, l, ...) as.numeric(HotellingsT2(X = data.frame(u[l == 1, ]), Y = data.frame(u[l == -1, ]))$statistic)

computeKMMD <- function(u, km, l, ...){
  y <- as.numeric(as.character(l))
  x.ind <- which(y == -1)
  y.ind <- which(y == 1)
  makeKM <- function(mat){
    class(mat) <- "kernelMatrix"
    mat
  }
  kmmd(x = makeKM(km[x.ind, x.ind]), y = makeKM(km[y.ind, y.ind]), Kxy = makeKM(km[x.ind, y.ind]))@mmdstats[1]
}

computeFS <- function(u, km, l, C = 1, ...){
  ksvm.fit <- ksvm(x = km, y = l, C = C, shrinking = FALSE, tol = .01, ...)
  ##ksvm.fit <- ksvm(x = km, y = l, C = C)
  y <- as.numeric(as.character(l))
  alpha <- unlist(ksvm.fit@alpha)
  sv.ind <- ksvm.fit@SVindex
  km.sub <- km[, sv.ind]
  b <- ksvm.fit@b
  ## print(ksvm.fit)
  eps <- 1e-3
  if(!all(alpha >= -eps)) warning("alpha: ", alpha)
  if((abs(sum(alpha * y[sv.ind])) > eps)) warning("complementary slackness: ", sum(alpha * y[sv.ind]))
  ## support vector 0 doesn't enter in to the margin
  ##margins <- km.sub %*% alpha * y + b
  ##margins <- fitted(ksvm(x = km, y = y, C = C, shrinking = FALSE, tol = .01, ...))
  margins <- km.sub %*% alpha + b
  ##browser()
  as.numeric(computeT(margins, km, y))
}

## rejectT2 <- function(u, l, ...) as.numeric(HotellingsT2(X = data.frame(u[l == 1, ]), Y = data.frame(u[l == -1, ]))$p.value < .05)
## rejectKMMD <- function(u, l, ...) as.numeric(computeKMMD(u, l, ...) > max(laply(1:19, function(i) computeKMMD(u, sample(l), ...))))
## rejectFS <- function(u, l, ...) as.numeric(computeFS(u, l, ...) > max(laply(1:19, function(i) computeFS(u, sample(l), ...))))

reject <- function(compute, verbose = FALSE){
  function(u, km, l, ...){
    val <- compute(u, km, l, ...)
    ##perms <- laply(1:19, function(i) compute(u, km, sample(l), ...))
    perms <- laply(1:39, function(i) compute(u, km, sample(l), ...))
    if(verbose){
      print(paste("value:", val))
      print("permuted: ")
      print(round(sort(perms), 3))
    }
    ##ifelse((val > max(perms)), 1, 0)
    ifelse((val > max(perms) | val < min(perms)), 1, 0)
  }
}

myplot <- function(plot, name){
  png(paste(imgDir, name, sep = ""), width = 800, height = 600)
  print(plot)
  dev.off()
}

getDataNormal <- function(N, D = 1, delta = 1){
  u <- rbind(matrix(rnorm(N * D, 0), ncol = D), matrix(rnorm(N * D, delta), ncol = D))
  ##u <- rbind(matrix(rnorm(N * D, -delta / 2), ncol = D), matrix(rnorm(N * D, delta / 2), ncol = D))
  l <- factor(c(rep(-1, N), rep(1, N)))
  list("u" = u, "l" = l)
}
