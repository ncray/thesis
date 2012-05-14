library(kernlab)
library(plyr)
library(ggplot2)
library(doMC)
registerDoMC(8)

computeT <- function(u, l) t.test(u[l == 1], u[l == -1], var.equal = TRUE)$statistic
getData <- function(N, D = 1){
  u <- rbind(matrix(rnorm(N * D, -1), ncol = D), matrix(rnorm(N * D, 1), ncol = D))
  l <- c(rep(-1, N), rep(1, N))
  list("u" = u, "l" = l)
}
friedStat <- function(km, l, C = 1, epsilon = .1){
  ksvm.mod <- ksvm(km, l, C = C, epsilon = epsilon)
  #sample(c(-1, 1), 1) * as.numeric(computeT(fitted(ksvm.mod), l))
  as.numeric(computeT(fitted(ksvm.mod), l))
}

D <- 2
N <- 10
dat <- getData(N, D)
u <- dat$u
l <- dat$l
u <- scale(u)
degree <- 5
ker <- polydot(degree = degree, offset = 0)
km <- kernelMatrix(ker, u)
sort(eigen(km)$values)
range(eigen(km)$values)
system.time(ksvm(km, l))

## dat <- getData(20, 1)
## u <- dat$u
## l <- dat$l
## ker <- vanilladot()
## km <- kernelMatrix(ker, u)
## image(km)
## ksvm1 <- ksvm(km, l)
## computeT(fitted(ksvm1), l)
## ksvm2 <- ksvm(km, l, C = 2, epsilon = .05)
## computeT(fitted(ksvm2), l)
##d = 1, C and epsilon don't affect t-stat, not true for d > 1

res <- ldply(10^(0:3), function(D){
  dat <- getData(100, D)
  u <- dat$u
  l <- dat$l
  ker <- vanilladot()
  km <- kernelMatrix(ker, u)
  ldply(1:100, function(x){
    l.perm <- sample(l)
    data.frame("D" = D, "T" = friedStat(km, l.perm))
  }, .parallel = TRUE)
})
ggplot(res, aes(sample = T)) +
  geom_point(stat = "qq", distribution = qnorm) + 
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~D, scales = "free") +
  opts(title = "Normal QQ Plots Faceted by Dimension")
ggsave("multivar_qq.png", width = 6, height = 4.5)

res <- ldply(1:4, function(degree){
  print(paste("degree:", degree))
  dat <- getData(100, 1)
  u <- dat$u
  l <- dat$l
  ker <- polydot(degree = degree)
  km <- kernelMatrix(ker, u)
  ldply(1:100, function(x){
    print(x)
    l.perm <- sample(l)
    data.frame("degree" = degree, "T" = friedStat(km, l.perm))
  }, .parallel = TRUE)
})
ggplot(res, aes(sample = T)) +
  geom_point(stat = "qq", distribution = qnorm) + 
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~degree, scales = "free") +
  opts(title = "Normal QQ Plots Faceted by Degree")
ggsave("poly_ker_qq.png", width = 6, height = 4.5)

## dat <- getData(20, 2)
## u <- dat$u
## l <- dat$l
## kernelMatrix(polydot(degree = 1), u) - kernelMatrix(vanilladot(), u)
## ker <- vanilladot()
## km <- kernelMatrix(ker, u)
## friedStat(kernelMatrix(vanilladot(), u), l)
## friedStat(kernelMatrix(polydot(degree = 1, offset = 0), u), l)

res <- mdply(expand.grid("D" = 10^(0:3), "degree" = 1:4), function(D, degree){
  print(paste("D:", D, "--", "degree", degree))
  dat <- getData(100, D)
  u <- dat$u
  l <- dat$l
  ker <- polydot(degree = degree)
  km <- kernelMatrix(ker, u)
  ldply(1:100, function(x){
    print(x)
    l.perm <- sample(l)
    data.frame("degree" = degree, "T" = friedStat(km, l.perm))
  }, .parallel = TRUE)
})
ggplot(res, aes(sample = T)) +
  geom_point(stat = "qq", distribution = qnorm) + 
  geom_abline(intercept = 0, slope = 1) +
  #facet_grid(D~degree, scales = "free") +
  facet_wrap(D~degree, scales = "free") +
  opts(title = "Normal QQ Plots Faceted by Dimension and Degree")
ggsave("multivar_poly_ker_qq.png", width = 6, height = 4.5)

getTp <- function(u, l, km){
  T <- friedStat(km, l)
  laply(1:100, function(i){
    print(i)
    l2 <- swap(l)
    friedStat(km, l2)
  }, .parallel = TRUE)
}
swap <- function(l){
  minus <- which(l == -1)
  plus <- which(l == 1)
  l[sample(minus, 1)] = 1
  l[sample(plus, 1)] = -1
  l
}
sim <- function(N, D){
  print(paste("N:", N, " -- D:", D))
  dat <- getData(N, D)
  u <- dat$u
  l <- dat$l
  ker <- vanilladot()
  km <- kernelMatrix(ker, u)
  nperm <- 8
  res <- ldply(1:nperm, function(i){
    print(i)
    l <- sample(l)
    T <- friedStat(km, l)
    Tp <- getTp(u, l, km)
    data.frame("N" = N, "T" = T, "Tprime" = Tp)
  })
  res
}
dat <- mdply(expand.grid("N" = 10^(1:3), "D" = 10^(0:2)), sim)

ggplot(dat, aes(T, Tprime)) +
  geom_point(alpha = .1) +
  xlab("T") +
  ylab("T'") +
  opts(title = "T' on T Faceted by Dimension and Sample Size") + 
  facet_wrap(D~N, scales = "free")
ggsave("multivar_ARC.png", width = 6, height = 4.5)
##there are some bad cases of km which end up in much slower simulations

sim <- function(N, D = 1, degree = 1){
  print(paste("N:", N, " -- D:", D, " -- degree:", degree))
  dat <- getData(N, D)
  u <- dat$u
  l <- dat$l
  u <- scale(u)  
  ker <- polydot(degree = degree, offset = 1)
  km <- kernelMatrix(ker, u)
  nperm <- 8
  res <- ldply(1:nperm, function(i){
    print(i)
    l <- sample(l)
    T <- friedStat(km, l)
    Tp <- getTp(u, l, km)
    data.frame("N" = N, "T" = T, "Tprime" = Tp)
  })
  res
}
dat <- mdply(expand.grid("N" = c(10, 100, 200), "degree" = c(1, 2, 5)), sim)
##scale, 0 offset 1.4, 1.4, 1.4 (11.5, 14.4)
##no scale, 0 offset 8, 113
##scale, 1 offset 1.4, 2.9, 3.7 (17.6, 17.4)
ggplot(dat, aes(T, Tprime)) +
  geom_point(alpha = .1) +
  xlab("T") +
  ylab("T'") +
  opts(title = "T' on T Faceted by Degree and Sample Size") + 
  facet_grid(degree~N, scales = "fixed")
ggsave("poly_ker_ARC.png", width = 6, height = 4.5)
