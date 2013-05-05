myTikz <- function(filename, plot){
  tikz(paste(imgDir, filename, sep = ""), width = 6, height = 4.5)
  print(plot)
  dev.off()
}

getData <- function(n, u = function(n) c(rnorm(n, -1), rnorm(n, 1)), name = "Normal"){
  l <- c(rep(-1, n), rep(1, n))
  u <- u(n)
  u <- u - mean(u)
  u <- u * sqrt(1 / sum(u^2) * 2 * n)
  list("u" = u, "l" = l, "name" = name)
}

getDataSpike <- function(n){
  u <- c(rcauchy(n, -1), rcauchy(n, 1))
  l <- c(rep(-1, n), rep(1, n))
  u <- u - mean(u)
  u <- u * sqrt(1 / sum(u^2) * 2 * n)
  list("u" = u, "l" = l)
}

swap <- function(l){
  minus <- which(l == -1)
  plus <- which(l == 1)
  l[sample(minus, 1)] = 1
  l[sample(plus, 1)] = -1
  l
}

u2BarP <- function(u, l, n, p) mean(u[l == 1])^p / n^(-p / 2)

dMinusP <- function(u, l, n, p) (2 * n * (1 - mean(u[l == 1])^2))^(-p / 2) / n^(-p / 2)

hP <- function(u, l, n, p){
  i1 <- sample(1:n, 1)
  i2 <- sample((n + 1):(2 * n), 1)
  d <- (u[i1] - u[i2])
  (abs(4 * mean(u[l == 1]) * d + 2 / n * d^2))^p / n^(-p / 2)
}

dDiffP <- function(u, l, n, p){
  d <- sqrt(2 * n * (1 - mean(u[l == 1])^2))
  l <- swap(l)
  dp <- sqrt(2 * n * (1 - mean(u[l == 1])^2))
  (d - dp)^p / n^(-p)
}

qP <- function(u, l, n, p) (-2 * n * mean(u[l == 1]))^p / n^(p / 2)

qddP <- function(u, l, n, p){
  d <- sqrt(2 * n * (1 - mean(u[l == 1])^2))
  l <- swap(l)
  dp <- sqrt(2 * n * (1 - mean(u[l == 1])^2))
  qp <- -2 * n * mean(u[l == 1])
  (qp / (d * dp))^p / n^(-p / 2)
}

sim <- function(n, p, funcs, labels){
  dat <- getData(n)
  u <- dat$u
  l <- dat$l
  nperm <- 10000
  res <- laply(1:nperm, function(i){
    u <- sample(u)
    funcs(u, l, n, p)
  })
  means <- apply(res, 2, mean)
  bounds <- aaply(res, 2, function(col)
                  quantile(as.vector(boot(col, function(dat, ind) mean(dat[ind]), 1000)$t),
                           c(.025, .975)))
  dat <- data.frame(n, p, "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2])
  data.frame(dat, "group" = labels)
}

#simOne <- function(...)  sim(..., funcs = each(u2BarP, dMinusP, hP), labels = c("u2BarP", "dMinusP", "hP"))
simOne <- function(...)  sim(...,
                             funcs = each(u2BarP, dMinusP, hP),
                             labels = c("$n^{p/2}\\mathbb{E} \\bar{u}_{2, \\Pi}^p \\;$",
                               "$n^{p/2}(d_{\\Pi})^{-p}\\;$",
                               "$n^{p/2}h_{\\Pi}\\;$")
                             )
#simTwo <- function(...)  sim(..., funcs = each(dDiffP, qP, qddP), labels = c("dDiffP", "qP", "qddP"))
simTwo <- function(...)  sim(...,
                             funcs = each(dDiffP, qP, qddP),
                             labels = c("$n^{p}\\mathbb{E}|d_{\\Pi}-d'_{\\Pi}|^p\\;\\;$",
                               "$n^{-p/2}\\mathbb{E} q_{\\Pi}^p \\;\\;$",
                               "$n^{p/2}\\mathbb{E}\\left [ \\left ( \\frac{q'_{\\Pi}}{d_{\\Pi}d'_{\\Pi}} \\right ) ^p \\right ]\\;\\;$")
                             )

computeT <- function(u, l) t.test(u[l == 1], u[l == -1], var.equal = TRUE)$statistic

simARC <- function(n){
  dat <- getData(n)
  u <- dat$u
  l <- dat$l
  nperm <- 20
  res <- ldply(1:nperm, function(i){
    u <- sample(u)
    T <- computeT(u, l)
    ##diff <- getTpminusT(u, l, n, p)
    diff <- getTpminusTMC(u, l, 50, p)
    R <- n / 2 * diff + T
    data.frame("n" = n, "T" = T, "Tprime" = T + diff)
  }, .parallel = TRUE)
  res
}

getTpminusT <- function(u, l, n, p){
  x <- u[l == -1]
  y <- u[l == 1]
  del <- rep(y, length(x)) - rep(x, each = length(y))
  xbar <- mean(x)
  ybar <- mean(y)
  Tprime <- -(xbar - ybar + 2/n*del) /
    (sqrt(2/n)*sqrt(sum(u^2)/(2*(n-1)) - 1/2*n/(n-1)*(xbar^2 + ybar^2 + 2*del/n*(xbar-ybar) + 2*del^2/n^2)))
  Tprime - computeT(u, l)
}

getTpminusTMC <- function(u, l, n, p){
  T <- computeT(u, l)
  laply(1:n, function(i){
    l2 <- swap(l)
    computeT(u, l2) - T
  })
}

simVar <- function(n){
  dat <- getData(n)
  u <- dat$u
  l <- dat$l
  nperm <- 100000
  res <- laply(1:nperm, function(i){computeT(sample(u), l)})
  f3 <- function(x, ind = 1:length(x)) abs(var(x[ind]) - 1) * n
  lims <- c(.025, .975)  
  bootf <- function(x, f) quantile(as.vector(boot(x, f, 1000)$t), lims)
  #boot(res, f3, 1000)$t
  bounds <- bootf(res, f3)
  data.frame(n, "value" = f3(res), "lower" = bounds[1], "upper" = bounds[2], group = 1)
}

simOrig <- function(n, exact = TRUE, scaled = TRUE, ...){
  dat <- getData(n, ...)
  u <- dat$u
  ##u <- 1:(2 * n)
  l <- dat$l
  nperm <- 10 * n
  res <- llply(1:nperm, function(i){
    u <- sample(u)
    T <- computeT(u, l)
    if(exact){
      diff <- getTpminusT(u, l, n, p)
    } else {
      diff <- getTpminusTMC(u, l, n, p)
    }
    R <- n / 2 * diff + T
    c(mean(abs(diff))^3, mean(diff^2), T, mean(T * R), mean(R))
  }, .parallel = TRUE)
  res <- do.call(rbind, res)
  f1 <- function(x, ind = 1:length(x)) (2 * pi)^(-1 / 4) * sqrt(mean(x[ind]) * n / 2) * if(scaled) n^(1 / 4) else 1
  f2 <- function(x, ind = 1:length(x)) 4 * n * sqrt(var(x[ind])) * if(scaled) n else 1
  f3 <- function(x, ind = 1:length(x)) abs(var(x[ind]) - 1) * if(scaled) n else 1
  f4 <- function(x, ind = 1:length(x)) mean(abs(x[ind])) * if(scaled) n^(1 / 2) else 1
  f5 <- function(x, ind = 1:length(x)) mean(abs(x[ind])) * if(scaled) n^(1 / 2) else 1
  allF <- function(df, ind = 1:nrow(df)) f1(df[, 1], ind) + f2(df[, 2], ind) + f3(df[, 3], ind) +
    f4(df[, 4], ind) + f5(df[, 5], ind)

  means <- c(f1(res[, 1]), f2(res[, 2]), f3(res[, 3]), f4(res[, 4]), f5(res[, 5]), allF(res))
  lims <- c(.025, .975)
  bootf <- function(x, f) quantile(as.vector(boot(x, f, 1000)$t), lims)
  bounds <- rbind(bootf(res[, 1], f1),
                  bootf(res[, 2], f2),
                  bootf(res[, 3], f3),
                  bootf(res[, 4], f4),
                  bootf(res[, 5], f5),
                  bootf(res, allF))
  labels <- c("$(2\\pi)^{-1/4}\\sqrt{\\frac{\\mathbb{E}|T\'_{\\Pi}-T_{\\Pi}|^3}{\\lambda}}n^{1/4}\\quad $",
              "$\\frac{1}{2\\lambda}\\sqrt{\\mathrm{Var}(\\mathbb{E}[(T\'_{\\Pi}-T_{\\Pi})^2|T_{\\Pi}])}n\\quad $",
              "$|\\mathbb{E}T_{\\Pi}^2-1|n\\quad $",
              "$\\mathbb{E}|T_{\\Pi}R_{\\Pi}|n^{1/2}\\quad $",
              "$\\mathbb{E}|R_{\\Pi}|n^{1/2}\\quad $",
              "Sum of Bounds")
  data.frame(n,
             "value" = means,
             "lower" = bounds[, 1],
             "upper" = bounds[, 2],
             group = factor(labels, levels = labels),
             distribution = dat$name)
}

getDelta <- function(u){
  n <- length(u) / 2
  getT <- function(u, l) as.vector(t.test(u[l], u[-l], var.equal = TRUE)$statistic)
  sorted <- sort(u)
  T <- getT(sorted, 1:n)
  Tprime <- getT(sorted, c(2:n, 2*n))
  abs(T - Tprime)
}

simBetterBound <- function(n, exact = TRUE, scaled = TRUE, ...){
  dat <- getData(n, ...)
  u <- dat$u
  l <- dat$l
  delta <- getDelta(u)
  nperm <- 10 * n
  res <- llply(1:nperm, function(i){
    u <- sample(u)
    T <- computeT(u, l)
    if(exact){
      diff <- getTpminusT(u, l, n, p)
    } else {
      diff <- getTpminusTMC(u, l, n, p)
    }
    R <- n / 2 * diff + T
    c(mean(abs(diff))^3, mean(diff^2), T, mean(T * R), mean(R))
  }, .parallel = TRUE)
  res <- do.call(rbind, res)
  f0 <- function(x, ind = 1:length(x)) n / 2 * .41 * delta^3 * if(scaled) n^(1 / 2) else 1
  f1 <- function(x, ind = 1:nrow(x)) 3 * delta * (sqrt(var(x[, 1][ind])) + mean(abs(x[, 2][ind]))) * if(scaled) n else 1
  f2 <- function(x, ind = 1:length(x)) 4 * n * sqrt(var(x[ind])) * if(scaled) n else 1
  f3 <- function(x, ind = 1:length(x)) abs(var(x[ind]) - 1) * if(scaled) n else 1
  f4 <- function(x, ind = 1:length(x)) mean(abs(x[ind])) * if(scaled) n^(1 / 2) else 1
  f5 <- function(x, ind = 1:length(x)) mean(abs(x[ind])) * if(scaled) n^(1 / 2) else 1
  allF <- function(df, ind = 1:nrow(df)) f0(0, ind) + f1(df[, c(2, 4)], ind) + f2(df[, 2], ind) + f3(df[, 3], ind) +
    f4(df[, 4], ind) + f5(df[, 5], ind)
  
  means <- c(f0(0), f1(res[, c(2, 4)]), f2(res[, 2]), f3(res[, 3]), f4(res[, 4]), f5(res[, 5]), allF(res))
  lims <- c(.025, .975)
  bootf <- function(x, f) quantile(as.vector(boot(x, f, 1000)$t), lims)
  bounds <- rbind(bootf(0, f0),
                  bootf(res[, c(2, 4)], f1),
                  bootf(res[, 2], f2),
                  bootf(res[, 3], f3),
                  bootf(res[, 4], f4),
                  bootf(res[, 5], f5),
                  bootf(res, allF))
  labels <- c("$\\frac{.41\\delta^3}{\\lambda}n^{1/2}\\quad $",
              "$3\\delta(\\sqrt{\\mathbb{E}T^2}+\\mathbb{E}|R|)n\\quad $",
              "$\\frac{1}{2\\lambda}\\sqrt{\\mathrm{Var}(\\mathbb{E}[(T\'_{\\Pi}-T_{\\Pi})^2|T_{\\Pi}])}n\\quad $",
              "$|\\mathbb{E}T_{\\Pi}^2-1|n\\quad $",
              "$\\mathbb{E}|T_{\\Pi}R_{\\Pi}|n^{1/2}\\quad $",
              "$\\mathbb{E}|R_{\\Pi}|n^{1/2}\\quad $",
              "Sum of Bounds")
  data.frame(n, "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = factor(labels, levels = labels))
}

library(nortest)
tDist <- function(n, exact = TRUE, scaled = TRUE, ...){
  dat <- getData(n, ...)
  u <- dat$u
  l <- dat$l
  nperm <- 1000
  samps <- laply(1:nperm, function(i) computeT(sample(u), l))
  c("n" = n, "pval" = ad.test(samps)$p.value, "ks" = ks.test(samps, "pnorm")$statistic)
}

