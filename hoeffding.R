getMat <- function(n, normal = TRUE){
  if (normal) {
    matrix(rnorm(n^2), nrow = n)
  } else {
    mat <- matrix(0, nrow = n, ncol = n)
    mat[1, 1] <- 1
    mat
  }
}

standardizeMat <- function(mat){
  ai <- rowMeans(mat)
  aj <- colMeans(mat)
  aa <- 1 / prod(dim(mat)) * sum(mat)
  ##(mat - outer(ai, rep(1, nrow(mat))) - t(outer(aj, rep(1, nrow(mat)))) + aa) ##* sqrt(prod(dim(mat)) / (nrow(mat) - 1))
  ret <- mat - outer(ai, rep(1, nrow(mat))) - t(outer(aj, rep(1, nrow(mat)))) + aa
  sqrt((nrow(ret) - 1) / sum(ret^2)) * ret
}

mat <- standardizeMat(getMat(10, normal = TRUE))
mat <- standardizeMat(getMat(10, normal = FALSE))
rowMeans(mat)
colMeans(mat)
1 / (nrow(mat) - 1) * sum(mat^2)

c(sqrt(sum(mat^4)), sqrt(sum(abs(mat)^3)))

library(plyr)
library(reshape)
library(ggplot2)

dat <- mdply(expand.grid("n" = c(10, 20, 50, 100, 200, 500, 1000, 1500, 2000), "normal" = c(TRUE, FALSE)), function(n, normal){
  mat <- standardizeMat(getMat(n, normal))
  c(sqrt(sum(mat^4)), sqrt(sum(abs(mat)^3)))
})

dat.m <- melt(dat, id.vars = c("n", "normal"))
qplot(x = n, y = value, data = dat.m, geom = "path", color = normal, linetype = variable) + facet_grid(normal~., scales = "free")
