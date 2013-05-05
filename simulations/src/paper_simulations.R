setwd("~/Dropbox/VMshare/thesis/simulations/src")
imgDir <- "../img/"
source("./sim_funcs.R")

library(ggplot2)
library(boot)
library(reshape)
library(doMC)
#install.packages("tikzDevice", repos=c("http://r-forge.r-project.org", getOption("repos")))
library(tikzDevice)
registerDoMC(4)
options(tikzLatexPackages = c(getOption("tikzLatexPackages"),
          "\\usepackage{amsfonts}"))
getOption("tikzLatexPackages")

xbreaks <- floor(10^(seq(1, 2.5, by = .25)))


params <- expand.grid(xbreaks, seq(2, 8, 2))
names(params) <- c("n", "p")
system.time(sideRatesDF <- mdply(params, simOne, .parallel = TRUE, .progress = "text"))
save(sideRatesDF, file = "siderates_1.dat")
system.time(sideRates2DF <- mdply(params, simTwo, .parallel = TRUE, .progress = "text"))
save(sideRates2DF, file = "siderates_2.dat")

sideRatesPlot <- ggplot(sideRatesDF, aes(x = n, y = value, color = group, linetype = group)) +
  geom_line(size = 1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .2)) +
  facet_wrap(~ p) +
  scale_y_log10(breaks = 10^seq(-1, 10, 2)) +
  scale_x_log10(breaks = xbreaks) + 
  xlab("$n$") +
  theme(legend.position = "bottom")
sideRatesPlot
sideRatesPlot2 <- sideRatesPlot %+% sideRates2DF

myTikz("siderates_1.tex", sideRatesPlot)
myTikz("siderates_2.tex", sideRatesPlot2)

##Approximate regression condition
ARCDF <- ldply(10^(1:4), simARC)
save(ARCDF, file = "ARC.dat")

ARCPlot <- ggplot(ARCDF, aes(T, Tprime)) +
  geom_point(alpha = .1) +
  geom_line(aes(y = (1 - 2 / n) * T)) +
  xlab("$T_{\\Pi}$") +
  ylab("$T'_{\\Pi}$") +
  ggtitle("Approximate Regression Condition: $(1-\\lambda)T_{\\Pi}$ Line") +
  facet_wrap(~ n)

myTikz("ARC.tex", ARCPlot)



##system.time(origRateDF <- ldply(c(100, 200), simOrig, u = function(n) 1:(2*n), name = "integer", .progress = "text"))
system.time(origRateDF <- ldply(xbreaks, simOrig, .progress = "text", u = function(n) 1:(2*n), name = "integer") )
system.time(origRateDF <- ldply(xbreaks, simOrig, .progress = "text")) ##2 mins for 1k perm, 3 ##2 mins for 10k perm, 2.5
system.time(origRateMCDF <- ldply(xbreaks, simOrig, exact = FALSE, .progress = "text"))
system.time(origRateMCDF <- ldply(xbreaks, simOrig, exact = FALSE, .progress = "text", u = function(n) 1:(2*n), name = "integer"))
system.time(betterRateDF <- ldply(xbreaks, simBetterBound, .progress = "text"))
system.time(betterRateDF <- ldply(xbreaks, simBetterBound, .progress = "text", u = function(n) 1:(2*n), name = "integer"))

## system.time(betterRateDF <- ldply(xbreaks, simBetterBound, .progress = "text", u = function(n) c(rnorm(2*n-1,sd=.01), 1), name = "bad"))
## system.time(betterRateDF <- ldply(xbreaks, simBetterBound, .progress = "text", u = function(n) c(rnorm(n, -1, 1/n), rnorm(n, 1, 1/n)), name = "bad"))
## system.time(betterRateDF <- ldply(xbreaks, simBetterBound, .progress = "text", u = function(n) rcauchy(2*n), name = "cauchy"))

ratesPlot <- ggplot(subset(origRateDF, group != "Sum of Bounds"), aes(x = n, y = value, color = group, linetype = group)) +
  geom_line(size = 1.5) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .07), size = 1.5) +
  xlab("$n$") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  scale_y_log10(breaks = round(10^seq(-2.5, 2.5, .5), 3)) +
  scale_x_log10(breaks = xbreaks)
ratesPlot

save(origRateDF, file = "orig_rate.dat")
save(origRateMCDF, file = "better_rate.dat")
save(betterRateDF, file = "orig_rate_mc.dat")
myTikz("orig_rate.tex", ratesPlot)
myTikz("better_rate.tex", ratesPlot %+% subset(betterRateDF, group != "Sum of Bounds"))
myTikz("orig_rate_mc.tex", ratesPlot %+% subset(origRateMCDF, group != "Sum of Bounds"))

## df <- ldply(xbreaks, tDist, u = function(n) 1:(2*n), name = "integer", .parallel = TRUE)
## df <- ldply(xbreaks, tDist, u = function(n) c(rnorm(2*n-1,sd=.01), 1), name = "bad", .parallel = TRUE)
## df <- ldply(xbreaks, tDist, u = function(n) rcauchy(2*n), name = "cauchy", .parallel = TRUE)


simDelta <- function(n){
  calc <- function(x) .41 / 2 * n^(3/2) * getDelta(x)^3
  c("integer" = calc(1:(2 * n)),
    "big" = calc(c(1:(2 * n - 1), 2 * n)),
    "n" = n
    )
}

xbreaks <- floor(10^seq(1, 5, .5))
system.time(dat <- ldply(rep(xbreaks, 1), simDelta, .parallel = TRUE))
dat.m <- melt(dat, id.vars = "n")
qplot(x = n, y = value, data = dat.m, geom = "point", color = variable) +
  scale_x_log10(breaks = xbreaks) +
  scale_y_log10()


dat <- ldply(xbreaks, function(n) c(n = n, y = .41 / 2 * n^(3/2) * getDelta(1:(2*n))^3))
deltaplot <- qplot(x = n, y = y, data = dat, geom = "line", size = I(1.5)) +
  scale_x_log10(breaks = xbreaks) +
  scale_y_log10() +
  xlab("$n$") +
  ggtitle("$\\frac{.41\\delta^3}{\\lambda}n^{1/2}\\quad $")
myTikz("delta_plot.tex", deltaplot)

plot(laply(xbreaks, function(n) .41 / 2 * n^(3/2) * getDelta(rnorm(n))^3))




###KS DISTANCE###
getKSStat <- function(vec, ind = 1:length(vec)) ks.test(vec[ind], "pnorm")$statistic

permTVar <- function(u){
  u.perm <- sample(u)
  n <- length(u) / 2
  t.test(u.perm[1:n], u.perm[-(1:n)], var.equal = TRUE)$statistic
}

permT <- function(n){
  x <- rnorm(n, 0)
  y <- rnorm(n, 4)
  u <- c(x, y)
  getStat(ldply(1:n, .fun = function(l) permTVar(u)))
}

oneSim <- function(n){
  c(getStat(rnorm(n)), permT(n), n)
}

sim <- function(n){
  getOne <- function(n){
    dat <- getData(n)
    dat2 <- getDataSpike(n)
    unlist(Map(getStat,
               list(laply(1:n, function(i){computeT(sample(dat$u), dat$l)}),
                    laply(1:n, function(i){computeT(sample(dat2$u), dat2$l)}),
                    rnorm(n))
               ))
  }
  res <- laply(rep(n, 500), getOne, .parallel = TRUE)
  data.frame(n,
             "value" = apply(res, 2, median),
             "lower" = apply(res, 2, quantile, .025),
             "upper" = apply(res, 2, quantile, .975),
             group = c("Permutation-T (Normal)", "Permutation-T (Cauchy)", "Normal"))
}

###
library(combinat)
getT <- function(u, l) as.vector(t.test(u[l], u[-l], var.equal = TRUE)$statistic)
n <- 5
foo <- function(n){
  ##u <- 1:(2*n)
  ##u <- c(rnorm(n, mean = -1, sd = 1 / (5 * n)), rnorm(n, mean = 1, sd = 1 / (5 * n)))
  u <- c(1:(2*n-1), 10^n)
  combs <- combinat::combn(1:(2*n), n, simplify = FALSE)
  c(n = n, stat = getKSStat(laply(combs, function(l) getT(u, l))))
}
system.time(truedat <- ldply(2:10, foo, .parallel = TRUE))
system.time(dat <- ldply(2:8, foo, .parallel = TRUE))
system.time(dat <- ldply(1:5, function(i) ldply(2:8, foo, .parallel = TRUE)))
qplot(x = n, y = stat.D, data = dat, geom = "point")
qplot(x = n, y = stat.D * sqrt(n), data = dat, geom = "point")
ggplot(dat, aes(x = n, y = stat.D)) +
  geom_point() + 
  geom_smooth(method = "lm", formula = y ~ 1 + I(1/x^(1/2)), color = "blue") +
  geom_smooth(method = "lm", formula = y ~ 1 + I(1/x^(1/4)), color = "black")
###

getData <- function(n){
  u <- c(rnorm(n, -1, 1/n), rnorm(n, 1, 1/n))
  ##u <- 1:(2*n)
  u <- c(1:(2*n-1), 10^n)
  l <- c(rep(-1, n), rep(1, n))
  u <- u - mean(u)
  u <- u * sqrt(1 / sum(u^2) * 2 * n)
  list("u" = u, "l" = l)
}

KSSim <- function(n, getData, name, n1 = 100, n2 = 100){
  KSStats <- laply(1:n1, function(i){
    dat <- getData(n)
    u <- dat$u
    #u <- 1:(2*n)
    l <- dat$l
    getKSStat(laply(1:n2, function(i) computeT(sample(u), l)))
  }, .parallel = TRUE)
  qs <- as.numeric(quantile(KSStats, probs = c(.025, .5, .975)))
  data.frame(group = name,
             "value" = qs[2],
             "lower" = qs[1],
             "upper" = qs[3],
             "n" = n)
}

system.time(dat <- ldply(2:10, function(n) KSSim(n, getData, "integer", 1, 1000), .parallel = TRUE))
system.time(dat <- ldply(2:10, function(n) KSSim(n, getData, "integer", 1, 1000 * n), .parallel = TRUE))
plot(truedat$stat.D, dat$value)
##MC always overestimates the KS statistic
(truedat$stat.D - dat$value) / truedat$stat.D


KSSim(10, getData, "normal")

system.time(dat <- ldply(xbreaks, function(n) KSSim(n, getData, "normal", 1, 10000)))
system.time(dat <- ldply(xbreaks, function(n) KSSim(n, getData, "normal", 100, 100)))
p4 <- ggplot(dat, aes(x = n, y = value, color = group)) +
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .07)) 
p4

#sim(10)
dat4 <- ldply(xbreaks, sim, .progress = "text")
dat4 <- rbind(dat4, dat3.sum)

p4 <- ggplot(dat4, aes(x = n, y = value, color = group)) +
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .07)) +
  xlab("$n$") + 
  scale_x_log10(breaks = xbreaks) +
  scale_y_log10(breaks = c(.03, .05, .1, .5, 1, 2.5)) +
  geom_line(aes(y = 1 / sqrt(n)), color = "black", linetype = 2) +
  geom_line(aes(y = 1 / n^(1/4)), color = "black") +
  ggitle("Log/Log Plot of Rates")
p4

tikz('sim4.tex', width = 6, height = 4.5)
print(p4)
dev.off()

getRange <- function(n){
  ldply(1:100,
        function(i) data.frame(n,
                               "Normal" = range(getData(n)$u),
                               "Cauchy" = range(getDataSpike(n)$u)))
}
dat5 <- ldply(xbreaks, getRange, .progress = "text")
dat5 <- melt(dat5, id.vars = "n")
p5 <- ggplot(dat5, aes(x = value, fill = variable)) +
  geom_density(size = .7, alpha = .5) +
  facet_grid(n~.) +
  ggtitle("Densities of Min/Max Values Faceted on n")
p5
ggsave('sim5.png', width = 6, height = 4.5)

###KS DISTANCE###







##taylor expansion
n <- 100
tayOne <- function(n){
  x <- c(rnorm(n, -3), rnorm(n, 3))
  x <- sample(x)
  l <- c(rep(-1, n), rep(1, n))
  x <- x - mean(x)
  x <- x * sqrt(1 / sum(x^2) * 2 * n)
  u2bar <- mean(x[which(l==1)])
  d <- 2*n*(1-u2bar^2)
  z <- x
  z[1] <- x[n+1]
  z[n+1] <- x[1]
  u2barp <- mean(z[which(l==1)])
  dp <- 2*n*(1-u2barp^2)
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
n <- 1000
x <- c(rnorm(n, -3), rnorm(n, 3))
x <- c(rep(-1, n), rep(1, n))
x <- x - mean(x)
x <- x * sqrt(1 / sum(x^2) * 2 * n)
plot(x, rep(0, 2*n))
mean(x[(n+1):(2*n)])
##scaling sim


##third moment sim
n <- 4
n <- 10
x <- rnorm(2*n)
x <- x - mean(x)
x <- x * sqrt(1 / sum(x^2) * 2 * n)
var(unlist(lapply(1:100000, function(i) mean(sample(x)[1:n]))))
1/(2*n-1)
summary(unlist(lapply(1:100000, function(i) mean(sample(x)[1:n])^3)))
y <- unlist(lapply(1:100000, function(i) mean(sample(x)[1:n])^4))
mean(y)
3 / n^4 * n * (n-1) / (2*n) / (2*n-1) * (4*n^2 - sum(x^4)) + sum(x^4)/(2*n^4)

mean(unlist(lapply(1:20000, function(i) mean((sample(x)[1:3])^c(2,1,1)))))

library(combinat)
combs <- combinat::combn(1:(2*n), n, simplify = FALSE)
combs <- c(combs, lapply(combs, function(l) (1:(2*n))[-l]))
combs <- permn(1:(2*n))
y <- ldply(combs, .fun = function(l) mean(x[l][1:n])^4)
y <- ldply(combs, .fun = function(l) mean(x[l][1]^4))
y <- ldply(combs, .fun = function(l) prod(x[l][1:2]^c(1,3)))
y <- ldply(combs, .fun = function(l) prod(x[l][1:2]^c(2,2)))
y <- ldply(combs, .fun = function(l) prod(x[l][1:3]^c(2,1,1)))
y <- ldply(combs, .fun = function(l) prod(x[l][1:4]))
##2,1,1
1/(2*n*(2*n-1)*(2*n-2))*(-4*n^2+2*sum(x^4))
##1,1,1,1
1/(2*n*(2*n-1)*(2*n-2)*(2*n-3))*(12*n^2-6*sum(x^4))
mean(y$V1)

(n*sum(x^4)/(2*n) + 
4*n*(n-1)*mean(ldply(combs, .fun = function(l) prod(x[l][1:2]^c(1,3)))$V1) + 
3*n*(n-1)*mean(ldply(combs, .fun = function(l) prod(x[l][1:2]^c(2,2)))$V1) +
6*n*(n-1)*(n-2)*mean(ldply(combs, .fun = function(l) prod(x[l][1:3]^c(2,1,1)))$V1) +
n*(n-1)*(n-2)*(n-3)*mean(ldply(combs, .fun = function(l) prod(x[l][1:4]))$V1)
)/n^4 - mean(ldply(combs, .fun = function(l) mean(x[l][1:n])^4)$V1)

combs <- combinat::combn2(1:(2*n))
combs <- rbind(combs, combs[, 2:1])
combs <- combs[order(combs[,1], combs[,2]),]
y <- adply(combs, 1, .fun = function(r) prod(x[r]^c(1,3)))
y <- adply(combs, 1, .fun = function(r) prod(x[r]^c(2,2)))
mean(y$V1)
mean(x^4)

##1,3
-1/(2*n)*1/(2*n-1) * sum(x^4)
##2,2
2*n/(2*n-1)-sum(x^4)/(2*n*(2*n-1))

## count <- 1
## z <- rep(0, 2*n*(2*n-1))
## for(i in 1:(2*n)){
##   for(j in (1:(2*n))[-i]){
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
  n <- length(u) / 2
  t.test(u.perm[1:n], u.perm[-(1:n)], var.equal = TRUE)$statistic
}

fried1 <- function(n){
  x <- rnorm(n, 0)
  y <- rnorm(n, 4)
  u <- c(x, y)
  getStat(ldply(1:n, .fun = function(l) permTVar(u)))
}

fried2 <- function(n){
  x <- 1:n
  y <- (n+1):(2*n)
  u <- c(x, y)
  getStat(ldply(1:n, .fun = function(l) permTVar(u)))
}

oneSim <- function(n){
  c(getStat(rnorm(n)),
    getStat(rmvnorm(n = 1, sigma = getPD(matrix(rep(.1, n^2), nrow = n)))),
    getStat(rmvnorm(n = 1, sigma = getPD(matrix(rep(.9, n^2), nrow = n)))),
    getStat(rmvnorm(n = 1, sigma = getPD(matrix(rep(-.1, n^2), nrow = n)))),
    fried1(n),
    fried2(n),
    n)
}
oneSim(10)


dat <- ldply(rep(seq(10, 100, 10), 40), oneSim, .parallel = TRUE, .progress = "text")
dat <- ldply(rep(floor(10^(seq(1, 2.5, by = .25))), 40), oneSim, .parallel = TRUE, .progress = "text")
names(dat) <- c("ind", "1cor", "9cor", "-1cor", "fried1", "fried2", "n")
dat2 <- ddply(dat, .(n), mean)

dat.m <- melt(dat2, id.vars = "n")
qplot(n, log10(value), data = dat.m, geom = "line", color = variable) +
  geom_line(aes(y = log10(n^(-1/2))), col = "black", linetype = 2) + 
  geom_line(aes(y = log10(n^(-1/4))), col = "black", linetype = 2)
ggsave("ksconv.png", width = 8, height = 6)
###KS DISTANCE###


library(ggplot2)
library(doMC)
registerDoMC(4)

###VARIANCE SIMULATIONS###
library(combinat)
getTVarExact <- function(n){
  x <- rnorm(n, 0)
  y <- rnorm(n, 4)
  u <- c(x, y)
  combs <- combinat::combn(1:(2*n), n, simplify = FALSE)
  t.stats <- ldply(combs, .fun = function(l) t.test(u[l], u[-l], var.equal = TRUE)$statistic)
  c(n, var(t.stats))
}
dat <- ldply(seq(3, 10, by = 1), getTVarExact, .parallel = TRUE, .progress = "text")
names(dat) <- c("n", "var")
qplot(x = n, y = var, data = dat, geom = "line", main = "Exact Variance of T")
ggsave("varexact.png", width = 8, height = 6)

permTVar <- function(u){
  u.perm <- sample(u)
  n <- length(u) / 2
  t.test(u.perm[1:n], u.perm[-(1:n)], var.equal = TRUE)$statistic
}
getTVarMC <- function(n){
  x <- rnorm(n, 0)
  y <- rnorm(n, 4)
  u <- c(x, y)
  t.stats <- ldply(1:1000, .fun = function(l) permTVar(u))
  c(n, var(t.stats))
}
dat <- ldply(rep(floor(10^(seq(1, 3, by = .25))), 200), getTVarMC, .parallel = TRUE, .progress = "text")
names(dat) <- c("n", "var")
dat2 <- ddply(dat, .(n), function(df) mean(df$var))
names(dat2) <- c("n", "var")

qplot(x = n, y = var, data = dat2, geom = "line", main = "MC Variance of T (solid) and 1+70/n^2 (dashed)") +
  geom_line(aes(y = 1 + 70 / n^2), linetype = 2) +
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

computeAllCond2 <- function(T, n, u, l, x, y){
  minus <- which(l == -1)
  plus <- which(l == 1)
  Tprime <- 1:(n^2)
  for(j in 1:n){
    for(k in 1:n){
      swap <- c(minus[j], plus[k])
      l[swap] <- l[rev(swap)]
      Tprime[n*(j-1)+k] <- t.test(u[l==1], u[l==-1], var.equal=TRUE)$statistic
      l[swap] <- l[rev(swap)]
    }
  }
  data.frame("T" = T, "Tprime" = Tprime, "n" = n, "lambda" = 2 / n)
}

computeAllCond <- function(T, n, u, l, x, y){
  del <- rep(y, length(x)) - rep(x, each = length(y))
  xbar <- mean(x)
  ybar <- mean(y)
  Tprime <- -(xbar - ybar + 2/n*del) /
    (sqrt(2/n)*sqrt(sum(u^2)/(2*(n-1)) - 1/2*n/(n-1)*(xbar^2 + ybar^2 + 2*del/n*(xbar-ybar) + 2*del^2/n^2)))
  data.frame("T" = T, "Tprime" = Tprime, "n" = n, "lambda" = 2 / n)
}

######IMPORTANT, TRY
  x <- u[l]
  y <- u[-l]
  del <- rep(y, length(x)) - rep(x, each = length(y))
  xbar <- mean(x)
  ybar <- mean(y)
#######


system.time(computeAllCond(T, n, u, l, x, y))
system.time(computeAllCond2(T, n, u, l, x, y))
sum((sort(computeAllCond(T, n, u, l, x, y)$Tprime) - sort(computeAllCond2(T, n, u, l, x, y)$Tprime))^2)

##First, write the funtion for n=M to look at dependence on n.
simulateBounds <- function(n = 20){
  x <- rnorm(n, 1)
  y <- rnorm(n, 4)
  u <- sample(c(x, y))
  x <- u[1:n]
  y <- u[-(1:n)]
  l <- rep(c(-1, 1), each = n)
  T <- computeT(u, l)
##  dat <- ldply(1:1000, function(x) c(T, computeT(u, swap(l)), n, 2 / n))
##  names(dat) <- c("T", "Tprime", "n", "lambda")
  dat <- computeAllCond(T, n, u, l, x, y)
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

dat2 <- ddply(dat, .(n), function(df) max(abs(df$T-df$Tprime))) ##maybe zero bias approach?
qplot(n, log10(.13*V1), geom = "line", data = dat2) + geom_line(aes(y = log10(n^(-1/2))), col = "black", linetype = 2) + geom_line(aes(y = log10(n^(-1/4))), col = "black", linetype = 2)

dat2 <- ddply(dat, .(T, n), function(df)
              with(df,
                   summary(lm(T - Tprime ~ T))$sigma
                   )
              )
ddply(dat2, .(n), function(df) mean(df$V1))

dat2 <- ddply(dat, .(T, n), function(df)
              with(df,
                   c(mean((Tprime - T)^2),
                     mean(abs(Tprime - T)^3),
                     mean(T - Tprime),
                     lambda[1])
                   ),
              .parallel = TRUE,
              .progress = "text"
              )
names(dat2) <- c("T", "n", "one.T", "two.T", "three.T", "lambda")

##dat2.2 <- dat2
##dat2.2[, 3:5] <- .5 * dat2.2[, 3:5]
##dat2.2$lambda <- .5 * dat2.2$lambda

dat3 <- ddply(dat2, .(n), function(df)
              with(df,
                   c(1/(2*lambda[1]) * sqrt(var(one.T)),
                     (2*pi)^(-1/4) * sqrt(1/lambda[1] * mean(two.T)),
                     mean(abs(-1 / lambda[1] * three.T + T)))
                   ),
              .parallel = TRUE,
              .progress = "text"
              )
names(dat3) <- c("n", "one", "two", "three")
dat3.m <- melt(dat3, id.vars = "n")
qplot(n, log10(value), geom = "line", data = dat3.m, color = variable, linetype = variable, main = "Simulated Bounds with n^{-1/4, -1/2, -1} (dashed)") +
  geom_line(aes(y = log10(n^(-1))), col = "black", linetype = 2) + 
  geom_line(aes(y = log10(n^(-1/2))), col = "black", linetype = 2) + 
  geom_line(aes(y = log10(n^(-1/4))), col = "black", linetype = 2)
ggsave("boundsexact.png", width = 8, height = 6)

qplot(n, log10(value/(n^(-1))), geom = "line", data = dat3.m, color = variable)
ggsave("boundsone.png", width = 8, height = 6)
###ERROR BOUND SIMULATIONS###


###FIRST TERM SIMULATIONS###
computeAllCondOne <- function(T, n, x, y){
  uJ <- rep(y, length(x))
  uI <- rep(x, each = length(y))
  del <- uJ - uI
  u <- c(x, y)
  xbar <- mean(x)
  ybar <- mean(y)
  q <- sum(x - y)
  qprime <- q + 2*del
  d = sqrt(sum(u^2) - n * (xbar^2 + ybar^2))
##  sum((-sqrt((n-1)/n)*q/d-T)^2)
  dprime = sqrt(d^2 + 2/n*del*(-q-del))
##   Tprime <- -(xbar - ybar + 2/n*del) /
##     (sqrt(2/n)*sqrt(sum(u^2)/(2*(n-1)) - 1/2*n/(n-1)*(xbar^2 + ybar^2 + 2*del/n*(xbar-ybar) + 2*del^2/n^2)))
##   sum((sort(-sqrt((n-1)/n)*qprime/dprime) - sort(Tprime))^2)
  first <- 2*uJ-2*uI
  second <- qprime*(d/dprime - 1)
  data.frame("T" = T,
             "first.sq" = first^2,
             "all.sq" = (first+second)^2,
             "remainder" = second,
             "qprime" = qprime,
             "dprimes" = d/dprime - 1,
             "n" = n,
             "lambda" = 2 / n)
}

simulateBounds <- function(n = 20){
  x <- rnorm(n, 1)
  y <- rnorm(n, 4)
  u <- sample(c(x, y))
  x <- u[1:n]
  y <- u[-(1:n)]
  l <- rep(c(-1, 1), each = n)
  T <- computeT(u, l)
  dat <- computeAllCondOne(T, n, x, y)
  dat
}

dat <- ldply(rep(floor(10^(seq(1, 3.5, by=.5))), each = 4), simulateBounds, .parallel = TRUE, .progress = "text")
ddply(dat, .(n), function(df) summary(df$remainder), .parallel = TRUE, .progress = "text")
ddply(dat, .(n), function(df) summary(df$qprime), .parallel = TRUE, .progress = "text")
dat2 <- ddply(dat, .(n), function(df) max(abs(df$qprime)), .parallel = TRUE, .progress = "text")
names(dat2)[2] <- "max.abs.qprime"
qplot(n, max.abs.qprime / n, geom = "line", data = dat2)
ggsave("qprimenone.png", width = 8, height = 6)
qplot(n, max.abs.qprime / n^(1/2), geom = "line", data = dat2)
ggsave("qprimenhalf.png", width = 8, height = 6)

ddply(dat, .(n), function(df) summary(df$dprimes), .parallel = TRUE, .progress = "text")
dat2 <- ddply(dat, .(n), function(df) max(abs(df$dprimes)), .parallel = TRUE, .progress = "text")
names(dat2)[2] <- "max.abs.dprimes"
qplot(n, max.abs.dprimes * n, geom = "line", data = dat2)
ggsave("dprimesnone.png", width = 8, height = 6)
qplot(n, max.abs.dprimes * n^(3/2), geom = "line", data = dat2)
ggsave("dprimesnthreehalves.png", width = 8, height = 6)

dat2 <- ddply(dat, .(n, T), mean, .parallel = TRUE, .progress = "text")
qplot(n, log(abs(remainder)), data = dat2)
qplot(n, abs(dprimes * n^2), data = dat2)
dat3 <- ddply(dat2, .(n), mean, .parallel = TRUE,.progress = "text")[, -1]
###FIRST TERM SIMULATIONS###


###more sims
dat <- ldply(rep(10, 1000), simulateBounds, .parallel = TRUE, .progress = "text")
mean(dat$qprime)

###
