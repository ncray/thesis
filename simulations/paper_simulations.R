library(ggplot2)
library(boot)
library(reshape)
library(doMC)
registerDoMC(4)
#melt.data.frame(data.frame("ind" = 1:10, "x" = rnorm(10), "y" = rnorm(10)), id.vars = "ind")

getData <- function(N){
  u <- c(rnorm(N, -1), rnorm(N, 1))
  l <- c(rep(-1, N), rep(1, N))
  u <- u - mean(u)
  u <- u * sqrt(1 / sum(u^2) * 2 * N)
  list("u" = u, "l" = l)
}

getDataSpike <- function(N){
  u <- c(rcauchy(N, -1), rcauchy(N, 1))
  l <- c(rep(-1, N), rep(1, N))
  u <- u - mean(u)
  u <- u * sqrt(1 / sum(u^2) * 2 * N)
  list("u" = u, "l" = l)
}

range(getData(100)$u)
range(getDataSpike(100)$u)

u2BarP <- function(u, l, N, p) mean(u[l == 1])^p / N^(-p / 2)
dMinusP <- function(u, l, N, p) (2 * N * (1 - mean(u[l == 1])^2))^(-p / 2) / N^(-p / 2)
##fix swapping
hP <- function(u, l, N, p){
  samps <- sample(1:(2 * N), 2)
  i1 <- samps[1]
  i2 <- samps[2]
  d <- (u[i1] - u[i2])
  (abs(4 * mean(u[l == 1]) * d + 2 / N * d^2))^p / N^(-p / 2)
}

dDiffP <- function(u, l, N, p){
  samps <- sample(1:(2 * N), 2)
  d <- sqrt(2 * N * (1 - mean(u[l == 1])^2))
  u[samps] <- u[rev(samps)]
  dp <- sqrt(2 * N * (1 - mean(u[l == 1])^2))
  (d - dp)^p / N^(-p)
}
qP <- function(u, l, N, p) (-2 * N * mean(u[l == 1]))^p / N^(p / 2)
qddP <- function(u, l, N, p){
  samps <- sample(1:(2 * N), 2)
  d <- sqrt(2 * N * (1 - mean(u[l == 1])^2))
  u[samps] <- u[rev(samps)]
  dp <- sqrt(2 * N * (1 - mean(u[l == 1])^2))
  qp <- -2 * N * mean(u[l == 1])
  (qp / (d * dp))^p / N^(-p / 2)
}

sim <- function(N, p, funcs, labels){
  dat <- getData(N)
  u <- dat$u
  l <- dat$l
  nperm <- 10000
  res <- laply(1:nperm, function(i){
    u <- sample(u)
    funcs(u, l, N, p)
  })
  means <- apply(res, 2, mean)
  bounds <- aaply(res, 2, function(col)
                  quantile(as.vector(boot(col, function(dat, ind) mean(dat[ind]), 1000)$t),
                           c(.025, .975)))
  dat <- data.frame(N, p, "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2])
  data.frame(dat, "group" = labels)
}

#simOne <- function(...)  sim(..., funcs = each(u2BarP, dMinusP, hP), labels = c("u2BarP", "dMinusP", "hP"))
simOne <- function(...)  sim(..., funcs = each(u2BarP, dMinusP, hP), labels =
                             c("$N^{p/2}\\mathbb{E}[\\bar{u}_{2, \\Pi}^p]\\;$",
                               "$N^{p/2}(d_{\\Pi})^{-p}\\;$",
                               "$N^{p/2}h_{\\Pi}\\;$"))
#simTwo <- function(...)  sim(..., funcs = each(dDiffP, qP, qddP), labels = c("dDiffP", "qP", "qddP"))
simTwo <- function(...)  sim(..., funcs = each(dDiffP, qP, qddP), labels =
                             c("$N^{p}\\mathbb{E}|d_{\\Pi}-d'_{\\Pi}|^p\\;$",
                               "$N^{-p/2}\\mathbb{E}[q_{\\Pi}^p]\\;$",
                               "$N^{p/2}\\mathbb{E}\\left [ \\left ( \\frac{q'_{\\Pi}}{d_{\\Pi}d'_{\\Pi}} \\right ) ^p \\right ]\\;$"))
xbreaks <- floor(10^seq(1, 4, .5))
params <- expand.grid(xbreaks, seq(2, 8, 2))
names(params) <- c("N", "p")
system.time(dat <- mdply(params, simOne, .parallel = TRUE, .progress = "text"))
save(dat, file = "dat")
#dat <- ddply(dat, .(N, p), transform, value = value / N^(-p/2), lower = lower / N^(-p/2), upper = upper / N^(-p/2))
system.time(dat2 <- mdply(params, simTwo, .parallel = TRUE, .progress = "text"))
save(dat2, file = "dat2")

p1 <- ggplot(dat, aes(x = N, y = value, color = group)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .2)) +
  facet_wrap(~ p) +
  scale_y_log10(breaks = 10^seq(-1, 10, 2)) +
  scale_x_log10(breaks = xbreaks) + 
  xlab("$N$") +
  opts(legend.position = "bottom")
p1

p2 <- p1 %+% dat2

#dat2 <- unique(dat[, 1:2])
#dat2 <- transform(dat2, B2 = N^(-1), B4 = N^(-2), B6 = N^(-3), B8 = N^(-4), sq = N^2)
#dat2 <- melt(dat2, id.vars = c("N", "p"))
#geom_line(data = dat2, aes(linetype = variable), color = "black")


library(tikzDevice)
options(tikzLatexPackages = c(getOption("tikzLatexPackages"),
          "\\usepackage{amsfonts}"))
getOption("tikzLatexPackages")
tikz('sim1.tex', width = 6, height = 4.5)
print(p1)
dev.off()
tikz('sim2.tex', width = 6, height = 4.5)
print(p2)
dev.off()

computeT <- function(u, l) t.test(u[l == 1], u[l == -1], var.equal = TRUE)$statistic

getTpminusT <- function(u, l, N, p){
  x <- u[l == -1]
  y <- u[l == 1]
  del <- rep(y, length(x)) - rep(x, each = length(y))
  xbar <- mean(x)
  ybar <- mean(y)
  Tprime <- -(xbar - ybar + 2/N*del) /
    (sqrt(2/N)*sqrt(sum(u^2)/(2*(N-1)) - 1/2*N/(N-1)*(xbar^2 + ybar^2 + 2*del/N*(xbar-ybar) + 2*del^2/N^2)))
  Tprime - computeT(u, l)
}

swap <- function(l){
  minus <- which(l == -1)
  plus <- which(l == 1)
  l[sample(minus, 1)] = 1
  l[sample(plus, 1)] = -1
  l
}

## getTpminusTMC <- function(u, l, N, p){
##   T <- computeT(u, l)
##   mean(laply(1:1000, function(i){
##     l2 <- swap(l)
##     computeT(u, l2) - T
##   }))
## }

getTpminusTMC <- function(u, l, N, p){
  T <- computeT(u, l)
  laply(1:N, function(i){
    l2 <- swap(l)
    computeT(u, l2) - T
  })
}
      

#first <- function(u, l, N, p) mean(getTpminusT(u, l, N, p)^3)
#second <- function(u, l, N, p) mean(getTpminusT(u, l, N, p)^2)

simVar <- function(N){
  dat <- getData(N)
  u <- dat$u
  l <- dat$l
  nperm <- 100000
  res <- laply(1:nperm, function(i){computeT(sample(u), l)})
  f3 <- function(x, ind = 1:length(x)) abs(var(x[ind]) - 1) * N
  bootf <- function(x, f) quantile(as.vector(boot(x, f, 1000)$t), lims)
  #boot(res, f3, 1000)$t
  bounds <- bootf(res, f3)
  data.frame(N, "value" = f3(res), "lower" = bounds[1], "upper" = bounds[2], group = 1)
}

sim <- function(N){
  dat <- getData(N)
  u <- dat$u
  l <- dat$l
  nperm <- 100
  res <- llply(1:nperm, function(i){
    u <- sample(u)
    T <- computeT(u, l)
#    diff <- getTpminusT(u, l, N, p)
    diff <- getTpminusTMC(u, l, N, p)    
    R <- N / 2 * diff + T
    c(mean(abs(diff))^3, mean(diff^2), T, mean(T * R), mean(R))
  }, .parallel = TRUE)
  res <- do.call(rbind, res)
  f1 <- function(x, ind = 1:length(x)) (2 * pi)^(-1 / 4) * sqrt(mean(x[ind]) * N / 2) * N^(1 / 4)
  f2 <- function(x, ind = 1:length(x)) N * sqrt(var(x[ind])) * N
  f3 <- function(x, ind = 1:length(x)) abs(var(x[ind]) - 1) * N
  f4 <- function(x, ind = 1:length(x)) mean(abs(x[ind])) * N^(1 / 2)
  f5 <- function(x, ind = 1:length(x)) mean(abs(x[ind])) * N^(1 / 2)
  
  means <- c(f1(res[, 1]), f2(res[, 2]), f3(res[, 3]), f4(res[, 4]), f5(res[, 5]))
  lims <- c(.025, .975)
  bootf <- function(x, f) quantile(as.vector(boot(x, f, 1000)$t), lims)
  bounds <- rbind(bootf(res[, 1], f1),
                  bootf(res[, 2], f2),
                  bootf(res[, 3], f3),
                  bootf(res[, 4], f4),
                  bootf(res[, 5], f5))
  labels <- c("$(2\\pi)^{-1/4}\\sqrt{\\frac{\\mathbb{E}|T\'_{\\Pi}-T_{\\Pi}|^3}{\\lambda}}N^{1/4}\\quad $",
              "$\\frac{1}{2\\lambda}\\sqrt{\\mathrm{Var}(\\mathbb{E}[(T\'_{\\Pi}-T_{\\Pi})^2|T_{\\Pi}])}N\\quad $",
              "$|\\mathbb{E}T_{\\Pi}^2-1|N\\quad $",
              "$\\mathbb{E}|T_{\\Pi}R_{\\Pi}|N^{1/2}\\quad $",
              "$\\mathbb{E}|R_{\\Pi}|N^{1/2}\\quad $")
  data.frame(N, "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = labels)
}

xbreaks <- floor(10^(seq(1, 3, by = .5)))
system.time(dat3 <- ldply(xbreaks, sim, .progress = "text")) ##2 mins for 1k perm, 3 ##2 mins for 10k perm, 2.5
#save(dat3, file = "dat3")
dat <- ldply(floor(10^(seq(1, 5, by = .5))), simVar, .parallel = TRUE, .progress = "text")

p3 <- ggplot(dat3, aes(x = N, y = value, color = group)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .07)) +
  xlab("$N$") +
  #ylab("$\\log_{10}(\\mathrm{value})$") + 
  opts(legend.position = "bottom") +
  opts(legend.direction = "vertical") +
  #coord_trans(ytrans = "log10")
  scale_y_log10(breaks = round(10^seq(-2.5, 2.5, .5), 3)) +
  scale_x_log10(breaks = xbreaks)
p3
tikz('sim3.tex', width = 6, height = 4.5)
tikz('sim6.tex', width = 6, height = 4.5)
print(p3)
dev.off()
rescale <- function(df){
  if(df$group[[1]] == levels(dat3$group)[1]) df[, 2:4] <- df[, 2:4] / df$N^(1/4)
  if(df$group[[1]] == levels(dat3$group)[2]) df[, 2:4] <- df[, 2:4] / df$N
  if(df$group[[1]] == levels(dat3$group)[3]) df[, 2:4] <- df[, 2:4] / df$N^(1/2)
  if(df$group[[1]] == levels(dat3$group)[4]) df[, 2:4] <- df[, 2:4] / df$N
  if(df$group[[1]] == levels(dat3$group)[5]) df[, 2:4] <- df[, 2:4] / df$N^(1/2)
  df
}
dat3.rs <- ddply(dat3, .(group), rescale)
p3 %+% dat3.rs
dat3.sum <- data.frame(ddply(dat3.rs[, -5], .(N), function(df) colSums(df[, -1])), group = "Sum of Bounds")

##Approximate regression condition
sim <- function(N){
  dat <- getData(N)
  u <- dat$u
  l <- dat$l
  nperm <- 20
  res <- ldply(1:nperm, function(i){
    u <- sample(u)
    T <- computeT(u, l)
#    diff <- getTpminusT(u, l, N, p)
    diff <- getTpminusTMC(u, l, 100, p)    
    R <- N / 2 * diff + T
    data.frame("N" = N, "T" = T, "Tprime" = T + diff)
  }, .parallel = TRUE)
  res
}
dat <- ldply(10^(1:4), sim)
tikz('sim7.tex', width = 6, height = 4.5)
png('sim7.png', width = 6, height = 4.5, units = "in", res = 300)
ggplot(dat, aes(T, Tprime)) +
  geom_point(alpha = .1) +
  geom_line(aes(y = (1 - 2 / N) * T)) +
  xlab("$T_{\\Pi}$") +
  ylab("$T'_{\\Pi}$") +
  opts(title = "Approximate Regression Condition with $(1-\\lambda)T_{\\Pi} Line$") +
  facet_wrap(~ N)
dev.off()


###KS DISTANCE###
getStat <- function(vec, ind = 1:length(vec)) ks.test(vec[ind], "pnorm")$statistic

permTVar <- function(u){
  u.perm <- sample(u)
  N <- length(u) / 2
  t.test(u.perm[1:N], u.perm[-(1:N)], var.equal = TRUE)$statistic
}

permT <- function(N){
  x <- rnorm(N, 0)
  y <- rnorm(N, 4)
  u <- c(x, y)
  getStat(ldply(1:N, .fun = function(l) permTVar(u)))
}

oneSim <- function(N){
  c(getStat(rnorm(N)), permT(N), N)
}

sim <- function(N){
  dat <- getDataSpike(N)
  u <- dat$u
  l <- dat$l
  resT <- laply(1:N, function(i){computeT(sample(u), l)})
  resN <- rnorm(N)
  lims <- c(.025, .975)
  bootf <- function(x, f) quantile(as.vector(boot(x, f, 1000)$t), lims)
  bootT <- bootf(resT, getStat)
  bootN <- bootf(resN, getStat)
  bounds <- rbind(bootT, bootN)
  data.frame(N,
             "value" = c(getStat(resT), getStat(resN)),
             "lower" = bounds[, 1],
             "upper" = bounds[, 2],
             group = c("perm.t", "normal"))
}

sim <- function(N){
  getOne <- function(N){
    dat <- getData(N)
    dat2 <- getDataSpike(N)
    unlist(Map(getStat,
               list(laply(1:N, function(i){computeT(sample(dat$u), dat$l)}),
                    laply(1:N, function(i){computeT(sample(dat2$u), dat2$l)}),
                    rnorm(N))
               ))
  }
  res <- laply(rep(N, 500), getOne, .parallel = TRUE)
  data.frame(N,
             "value" = apply(res, 2, median),
             "lower" = apply(res, 2, quantile, .025),
             "upper" = apply(res, 2, quantile, .975),
             group = c("Permutation-T (Normal)", "Permutation-T (Cauchy)", "Normal"))
}

#sim(10)
dat4 <- ldply(xbreaks, sim, .progress = "text")
dat4 <- rbind(dat4, dat3.sum)

p4 <- ggplot(dat4, aes(x = N, y = value, color = group)) +
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .07)) +
  xlab("$N$") + 
  scale_x_log10(breaks = xbreaks) +
  scale_y_log10(breaks = c(.03, .05, .1, .5, 1, 2.5)) +
  geom_line(aes(y = 1 / sqrt(N)), color = "black", linetype = 2) +
  geom_line(aes(y = 1 / N^(1/4)), color = "black") +
  opts(title = "Log/Log Plot of Rates")
p4

tikz('sim4.tex', width = 6, height = 4.5)
print(p4)
dev.off()

getRange <- function(N){
  ldply(1:100,
        function(i) data.frame(N,
                               "Normal" = range(getData(N)$u),
                               "Cauchy" = range(getDataSpike(N)$u)))
}
dat5 <- ldply(xbreaks, getRange, .progress = "text")
dat5 <- melt(dat5, id.vars = "N")
p5 <- ggplot(dat5, aes(x = value, fill = variable)) +
  geom_density(size = .7, alpha = .5) +
  facet_grid(N~.) +
  opts(title = "Densities of Min/Max Values Faceted on N")
p5
ggsave('sim5.png', width = 6, height = 4.5)

###KS DISTANCE###







##taylor expansion
N <- 100
tayOne <- function(N){
  x <- c(rnorm(N, -3), rnorm(N, 3))
  x <- sample(x)
  l <- c(rep(-1, N), rep(1, N))
  x <- x - mean(x)
  x <- x * sqrt(1 / sum(x^2) * 2 * N)
  u2bar <- mean(x[which(l==1)])
  d <- 2*N*(1-u2bar^2)
  z <- x
  z[1] <- x[N+1]
  z[N+1] <- x[1]
  u2barp <- mean(z[which(l==1)])
  dp <- 2*N*(1-u2barp^2)
  h <- d^2-dp^2
  c(abs(d-dp),
    abs(h/(2*d)),
    h^2/(8*(d^2-abs(h))^(3/2)),
    abs(h)/(2*sqrt(d^2-max(0,h))),
    abs(h)/(2*min(d, dp)),
    max(abs(h)/(2*d), abs(h)/(2*dp)),
    abs(h)/(2*d) + abs(h)/(2*dp))
}
res <- ldply(rep(100, 1000), tayOne)
summary(with(res, V2+V3-V1))
summary(with(res, V4-V1))
summary(with(res, V5-V1))
summary(with(res, V6-V1))
summary(with(res, V7-V1))
##taylor expansion

##scaling sim
N <- 1000
x <- c(rnorm(N, -3), rnorm(N, 3))
x <- c(rep(-1, N), rep(1, N))
x <- x - mean(x)
x <- x * sqrt(1 / sum(x^2) * 2 * N)
plot(x, rep(0, 2*N))
mean(x[(N+1):(2*N)])
##scaling sim


##third moment sim
N <- 4
N <- 10
x <- rnorm(2*N)
x <- x - mean(x)
x <- x * sqrt(1 / sum(x^2) * 2 * N)
var(unlist(lapply(1:100000, function(i) mean(sample(x)[1:N]))))
1/(2*N-1)
summary(unlist(lapply(1:100000, function(i) mean(sample(x)[1:N])^3)))
y <- unlist(lapply(1:100000, function(i) mean(sample(x)[1:N])^4))
mean(y)
3 / N^4 * N * (N-1) / (2*N) / (2*N-1) * (4*N^2 - sum(x^4)) + sum(x^4)/(2*N^4)

mean(unlist(lapply(1:20000, function(i) mean((sample(x)[1:3])^c(2,1,1)))))

library(combinat)
combs <- combinat::combn(1:(2*N), N, simplify = FALSE)
combs <- c(combs, lapply(combs, function(l) (1:(2*N))[-l]))
combs <- permn(1:(2*N))
y <- ldply(combs, .fun = function(l) mean(x[l][1:N])^4)
y <- ldply(combs, .fun = function(l) mean(x[l][1]^4))
y <- ldply(combs, .fun = function(l) prod(x[l][1:2]^c(1,3)))
y <- ldply(combs, .fun = function(l) prod(x[l][1:2]^c(2,2)))
y <- ldply(combs, .fun = function(l) prod(x[l][1:3]^c(2,1,1)))
y <- ldply(combs, .fun = function(l) prod(x[l][1:4]))
##2,1,1
1/(2*N*(2*N-1)*(2*N-2))*(-4*N^2+2*sum(x^4))
##1,1,1,1
1/(2*N*(2*N-1)*(2*N-2)*(2*N-3))*(12*N^2-6*sum(x^4))
mean(y$V1)

(N*sum(x^4)/(2*N) + 
4*N*(N-1)*mean(ldply(combs, .fun = function(l) prod(x[l][1:2]^c(1,3)))$V1) + 
3*N*(N-1)*mean(ldply(combs, .fun = function(l) prod(x[l][1:2]^c(2,2)))$V1) +
6*N*(N-1)*(N-2)*mean(ldply(combs, .fun = function(l) prod(x[l][1:3]^c(2,1,1)))$V1) +
N*(N-1)*(N-2)*(N-3)*mean(ldply(combs, .fun = function(l) prod(x[l][1:4]))$V1)
)/N^4 - mean(ldply(combs, .fun = function(l) mean(x[l][1:N])^4)$V1)

combs <- combinat::combn2(1:(2*N))
combs <- rbind(combs, combs[, 2:1])
combs <- combs[order(combs[,1], combs[,2]),]
y <- adply(combs, 1, .fun = function(r) prod(x[r]^c(1,3)))
y <- adply(combs, 1, .fun = function(r) prod(x[r]^c(2,2)))
mean(y$V1)
mean(x^4)

##1,3
-1/(2*N)*1/(2*N-1) * sum(x^4)
##2,2
2*N/(2*N-1)-sum(x^4)/(2*N*(2*N-1))

## count <- 1
## z <- rep(0, 2*N*(2*N-1))
## for(i in 1:(2*N)){
##   for(j in (1:(2*N))[-i]){
## ##    z[count] <- x[i]^2*x[j]^2
##     z[count] <- x[i]^1*x[j]^3
##     count <- count + 1
##   }
## }
## mean(z)
##

###KS DISTANCE###
library(mvtnorm)
library(Matrix)
library(plyr)
library(ggplot2)
getPD <- function(mat){
  as.matrix(nearPD(mat, corr = TRUE)$mat)
}
getStat <- function(vec){
  ks.test(vec, "pnorm")$statistic
}

permTVar <- function(u){
  u.perm <- sample(u)
  N <- length(u) / 2
  t.test(u.perm[1:N], u.perm[-(1:N)], var.equal = TRUE)$statistic
}

fried1 <- function(N){
  x <- rnorm(N, 0)
  y <- rnorm(N, 4)
  u <- c(x, y)
  getStat(ldply(1:N, .fun = function(l) permTVar(u)))
}

fried2 <- function(N){
  x <- 1:N
  y <- (N+1):(2*N)
  u <- c(x, y)
  getStat(ldply(1:N, .fun = function(l) permTVar(u)))
}

oneSim <- function(N){
  c(getStat(rnorm(N)),
    getStat(rmvnorm(n = 1, sigma = getPD(matrix(rep(.1, N^2), nrow = N)))),
    getStat(rmvnorm(n = 1, sigma = getPD(matrix(rep(.9, N^2), nrow = N)))),
    getStat(rmvnorm(n = 1, sigma = getPD(matrix(rep(-.1, N^2), nrow = N)))),
    fried1(N),
    fried2(N),
    N)
}
oneSim(10)


dat <- ldply(rep(seq(10, 100, 10), 40), oneSim, .parallel = TRUE, .progress = "text")
dat <- ldply(rep(floor(10^(seq(1, 2.5, by = .25))), 40), oneSim, .parallel = TRUE, .progress = "text")
names(dat) <- c("ind", "1cor", "9cor", "-1cor", "fried1", "fried2", "N")
dat2 <- ddply(dat, .(N), mean)

dat.m <- melt(dat2, id.vars = "N")
qplot(N, log10(value), data = dat.m, geom = "line", color = variable) +
  geom_line(aes(y = log10(N^(-1/2))), col = "black", linetype = 2) + 
  geom_line(aes(y = log10(N^(-1/4))), col = "black", linetype = 2)
ggsave("ksconv.png", width = 8, height = 6)
###KS DISTANCE###


library(ggplot2)
library(doMC)
registerDoMC(4)

###VARIANCE SIMULATIONS###
library(combinat)
getTVarExact <- function(N){
  x <- rnorm(N, 0)
  y <- rnorm(N, 4)
  u <- c(x, y)
  combs <- combinat::combn(1:(2*N), N, simplify = FALSE)
  t.stats <- ldply(combs, .fun = function(l) t.test(u[l], u[-l], var.equal = TRUE)$statistic)
  c(N, var(t.stats))
}
dat <- ldply(seq(3, 10, by = 1), getTVarExact, .parallel = TRUE, .progress = "text")
names(dat) <- c("N", "var")
qplot(x = N, y = var, data = dat, geom = "line", main = "Exact Variance of T")
ggsave("varexact.png", width = 8, height = 6)

permTVar <- function(u){
  u.perm <- sample(u)
  N <- length(u) / 2
  t.test(u.perm[1:N], u.perm[-(1:N)], var.equal = TRUE)$statistic
}
getTVarMC <- function(N){
  x <- rnorm(N, 0)
  y <- rnorm(N, 4)
  u <- c(x, y)
  t.stats <- ldply(1:1000, .fun = function(l) permTVar(u))
  c(N, var(t.stats))
}
dat <- ldply(rep(floor(10^(seq(1, 3, by = .25))), 200), getTVarMC, .parallel = TRUE, .progress = "text")
names(dat) <- c("N", "var")
dat2 <- ddply(dat, .(N), function(df) mean(df$var))
names(dat2) <- c("N", "var")

qplot(x = N, y = var, data = dat2, geom = "line", main = "MC Variance of T (solid) and 1+70/N^2 (dashed)") +
  geom_line(aes(y = 1 + 70 / N^2), linetype = 2) +
  scale_y_continuous(limits = c(.9, 1.1))
ggsave("varMC.png", width = 8, height = 6)
###VARIANCE SIMULATIONS###

###ERROR BOUND SIMULATIONS###
swap <- function(l){
  minus <- which(l == -1)
  plus <- which(l == 1)
  l[sample(minus, 1)] = 1
  l[sample(plus, 1)] = -1
  l
}

computeT <- function(u, l){
  t.test(u[l==1], u[l==-1], var.equal=TRUE)$statistic
}

computeAllCond2 <- function(T, N, u, l, x, y){
  minus <- which(l == -1)
  plus <- which(l == 1)
  Tprime <- 1:(N^2)
  for(j in 1:N){
    for(k in 1:N){
      swap <- c(minus[j], plus[k])
      l[swap] <- l[rev(swap)]
      Tprime[N*(j-1)+k] <- t.test(u[l==1], u[l==-1], var.equal=TRUE)$statistic
      l[swap] <- l[rev(swap)]
    }
  }
  data.frame("T" = T, "Tprime" = Tprime, "N" = N, "lambda" = 2 / N)
}

computeAllCond <- function(T, N, u, l, x, y){
  del <- rep(y, length(x)) - rep(x, each = length(y))
  xbar <- mean(x)
  ybar <- mean(y)
  Tprime <- -(xbar - ybar + 2/N*del) /
    (sqrt(2/N)*sqrt(sum(u^2)/(2*(N-1)) - 1/2*N/(N-1)*(xbar^2 + ybar^2 + 2*del/N*(xbar-ybar) + 2*del^2/N^2)))
  data.frame("T" = T, "Tprime" = Tprime, "N" = N, "lambda" = 2 / N)
}

######IMPORTANT, TRY
  x <- u[l]
  y <- u[-l]
  del <- rep(y, length(x)) - rep(x, each = length(y))
  xbar <- mean(x)
  ybar <- mean(y)
#######


system.time(computeAllCond(T, N, u, l, x, y))
system.time(computeAllCond2(T, N, u, l, x, y))
sum((sort(computeAllCond(T, N, u, l, x, y)$Tprime) - sort(computeAllCond2(T, N, u, l, x, y)$Tprime))^2)

##First, write the funtion for N=M to look at dependence on N.
simulateBounds <- function(N = 20){
  x <- rnorm(N, 1)
  y <- rnorm(N, 4)
  u <- sample(c(x, y))
  x <- u[1:N]
  y <- u[-(1:N)]
  l <- rep(c(-1, 1), each = N)
  T <- computeT(u, l)
##  dat <- ldply(1:1000, function(x) c(T, computeT(u, swap(l)), N, 2 / N))
##  names(dat) <- c("T", "Tprime", "N", "lambda")
  dat <- computeAllCond(T, N, u, l, x, y)
  dat
}


dat <- ldply(rep(100, 100), simulateBounds, .parallel = TRUE, .progress = "text")
ddply(dat, .(T), function(df) sd(df$Tprime))
cond.dat <- ddply(dat, .(T), function(df) mean(df$Tprime))
names(cond.dat)[2] <- "Tprime"
summary(lm(T - Tprime ~ T, data = cond.dat))
qplot(T, T - Tprime, data = cond.dat, geom = c("point", "smooth"), method = "lm")

dat <- ldply(rep(floor(10^(seq(1, 2, by=.5))), each = 8), simulateBounds, .parallel = TRUE, .progress = "text")
dat <- ldply(rep(seq(10, 2000, 100), each = 20), simulateBounds, .parallel = TRUE, .progress = "text")

dat2 <- ddply(dat, .(N), function(df) max(abs(df$T-df$Tprime))) ##maybe zero bias approach?
qplot(N, log10(.13*V1), geom = "line", data = dat2) + geom_line(aes(y = log10(N^(-1/2))), col = "black", linetype = 2) + geom_line(aes(y = log10(N^(-1/4))), col = "black", linetype = 2)

dat2 <- ddply(dat, .(T, N), function(df)
              with(df,
                   summary(lm(T - Tprime ~ T))$sigma
                   )
              )
ddply(dat2, .(N), function(df) mean(df$V1))

dat2 <- ddply(dat, .(T, N), function(df)
              with(df,
                   c(mean((Tprime - T)^2),
                     mean(abs(Tprime - T)^3),
                     mean(T - Tprime),
                     lambda[1])
                   ),
              .parallel = TRUE,
              .progress = "text"
              )
names(dat2) <- c("T", "N", "one.T", "two.T", "three.T", "lambda")

##dat2.2 <- dat2
##dat2.2[, 3:5] <- .5 * dat2.2[, 3:5]
##dat2.2$lambda <- .5 * dat2.2$lambda

dat3 <- ddply(dat2, .(N), function(df)
              with(df,
                   c(1/(2*lambda[1]) * sqrt(var(one.T)),
                     (2*pi)^(-1/4) * sqrt(1/lambda[1] * mean(two.T)),
                     mean(abs(-1 / lambda[1] * three.T + T)))
                   ),
              .parallel = TRUE,
              .progress = "text"
              )
names(dat3) <- c("N", "one", "two", "three")
dat3.m <- melt(dat3, id.vars = "N")
qplot(N, log10(value), geom = "line", data = dat3.m, color = variable, linetype = variable, main = "Simulated Bounds with N^{-1/4, -1/2, -1} (dashed)") +
  geom_line(aes(y = log10(N^(-1))), col = "black", linetype = 2) + 
  geom_line(aes(y = log10(N^(-1/2))), col = "black", linetype = 2) + 
  geom_line(aes(y = log10(N^(-1/4))), col = "black", linetype = 2)
ggsave("boundsexact.png", width = 8, height = 6)

qplot(N, log10(value/(N^(-1))), geom = "line", data = dat3.m, color = variable)
ggsave("boundsone.png", width = 8, height = 6)
###ERROR BOUND SIMULATIONS###


###FIRST TERM SIMULATIONS###
computeAllCondOne <- function(T, N, x, y){
  uJ <- rep(y, length(x))
  uI <- rep(x, each = length(y))
  del <- uJ - uI
  u <- c(x, y)
  xbar <- mean(x)
  ybar <- mean(y)
  q <- sum(x - y)
  qprime <- q + 2*del
  d = sqrt(sum(u^2) - N * (xbar^2 + ybar^2))
##  sum((-sqrt((N-1)/N)*q/d-T)^2)
  dprime = sqrt(d^2 + 2/N*del*(-q-del))
##   Tprime <- -(xbar - ybar + 2/N*del) /
##     (sqrt(2/N)*sqrt(sum(u^2)/(2*(N-1)) - 1/2*N/(N-1)*(xbar^2 + ybar^2 + 2*del/N*(xbar-ybar) + 2*del^2/N^2)))
##   sum((sort(-sqrt((N-1)/N)*qprime/dprime) - sort(Tprime))^2)
  first <- 2*uJ-2*uI
  second <- qprime*(d/dprime - 1)
  data.frame("T" = T,
             "first.sq" = first^2,
             "all.sq" = (first+second)^2,
             "remainder" = second,
             "qprime" = qprime,
             "dprimes" = d/dprime - 1,
             "N" = N,
             "lambda" = 2 / N)
}

simulateBounds <- function(N = 20){
  x <- rnorm(N, 1)
  y <- rnorm(N, 4)
  u <- sample(c(x, y))
  x <- u[1:N]
  y <- u[-(1:N)]
  l <- rep(c(-1, 1), each = N)
  T <- computeT(u, l)
  dat <- computeAllCondOne(T, N, x, y)
  dat
}

dat <- ldply(rep(floor(10^(seq(1, 3.5, by=.5))), each = 4), simulateBounds, .parallel = TRUE, .progress = "text")
ddply(dat, .(N), function(df) summary(df$remainder), .parallel = TRUE, .progress = "text")
ddply(dat, .(N), function(df) summary(df$qprime), .parallel = TRUE, .progress = "text")
dat2 <- ddply(dat, .(N), function(df) max(abs(df$qprime)), .parallel = TRUE, .progress = "text")
names(dat2)[2] <- "max.abs.qprime"
qplot(N, max.abs.qprime / N, geom = "line", data = dat2)
ggsave("qprimeNone.png", width = 8, height = 6)
qplot(N, max.abs.qprime / N^(1/2), geom = "line", data = dat2)
ggsave("qprimeNhalf.png", width = 8, height = 6)

ddply(dat, .(N), function(df) summary(df$dprimes), .parallel = TRUE, .progress = "text")
dat2 <- ddply(dat, .(N), function(df) max(abs(df$dprimes)), .parallel = TRUE, .progress = "text")
names(dat2)[2] <- "max.abs.dprimes"
qplot(N, max.abs.dprimes * N, geom = "line", data = dat2)
ggsave("dprimesNone.png", width = 8, height = 6)
qplot(N, max.abs.dprimes * N^(3/2), geom = "line", data = dat2)
ggsave("dprimesNthreehalves.png", width = 8, height = 6)

dat2 <- ddply(dat, .(N, T), mean, .parallel = TRUE, .progress = "text")
qplot(N, log(abs(remainder)), data = dat2)
qplot(N, abs(dprimes * N^2), data = dat2)
dat3 <- ddply(dat2, .(N), mean, .parallel = TRUE,.progress = "text")[, -1]
###FIRST TERM SIMULATIONS###


###more sims
dat <- ldply(rep(10, 1000), simulateBounds, .parallel = TRUE, .progress = "text")
mean(dat$qprime)

###
