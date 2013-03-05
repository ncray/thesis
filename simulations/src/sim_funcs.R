myTikz <- function(filename, plot){
  tikz(paste(imgDir, filename, sep = ""), width = 6, height = 4.5)
  print(plot)
  dev.off()
}

getData <- function(N){
  u <- c(rnorm(N, -1), rnorm(N, 1))
  ##u <- 1:(2*N)
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

swap <- function(l){
  minus <- which(l == -1)
  plus <- which(l == 1)
  l[sample(minus, 1)] = 1
  l[sample(plus, 1)] = -1
  l
}

u2BarP <- function(u, l, N, p) mean(u[l == 1])^p / N^(-p / 2)

dMinusP <- function(u, l, N, p) (2 * N * (1 - mean(u[l == 1])^2))^(-p / 2) / N^(-p / 2)

hP <- function(u, l, N, p){
  i1 <- sample(1:N, 1)
  i2 <- sample((N + 1):(2 * N), 1)
  d <- (u[i1] - u[i2])
  (abs(4 * mean(u[l == 1]) * d + 2 / N * d^2))^p / N^(-p / 2)
}

dDiffP <- function(u, l, N, p){
  d <- sqrt(2 * N * (1 - mean(u[l == 1])^2))
  l <- swap(l)
  dp <- sqrt(2 * N * (1 - mean(u[l == 1])^2))
  (d - dp)^p / N^(-p)
}

qP <- function(u, l, N, p) (-2 * N * mean(u[l == 1]))^p / N^(p / 2)

qddP <- function(u, l, N, p){
  d <- sqrt(2 * N * (1 - mean(u[l == 1])^2))
  l <- swap(l)
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
simOne <- function(...)  sim(...,
                             funcs = each(u2BarP, dMinusP, hP),
                             labels = c("$N^{p/2}\\mathbb{E} \\bar{u}_{2, \\Pi}^p \\;$",
                               "$N^{p/2}(d_{\\Pi})^{-p}\\;$",
                               "$N^{p/2}h_{\\Pi}\\;$")
                             )
#simTwo <- function(...)  sim(..., funcs = each(dDiffP, qP, qddP), labels = c("dDiffP", "qP", "qddP"))
simTwo <- function(...)  sim(...,
                             funcs = each(dDiffP, qP, qddP),
                             labels = c("$N^{p}\\mathbb{E}|d_{\\Pi}-d'_{\\Pi}|^p\\;\\;$",
                               "$N^{-p/2}\\mathbb{E} q_{\\Pi}^p \\;\\;$",
                               "$N^{p/2}\\mathbb{E}\\left [ \\left ( \\frac{q'_{\\Pi}}{d_{\\Pi}d'_{\\Pi}} \\right ) ^p \\right ]\\;\\;$")
                             )

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

getTpminusTMC <- function(u, l, N, p){
  T <- computeT(u, l)
  laply(1:N, function(i){
    l2 <- swap(l)
    computeT(u, l2) - T
  })
}

simVar <- function(N){
  dat <- getData(N)
  u <- dat$u
  l <- dat$l
  nperm <- 100000
  res <- laply(1:nperm, function(i){computeT(sample(u), l)})
  f3 <- function(x, ind = 1:length(x)) abs(var(x[ind]) - 1) * N
  lims <- c(.025, .975)  
  bootf <- function(x, f) quantile(as.vector(boot(x, f, 1000)$t), lims)
  #boot(res, f3, 1000)$t
  bounds <- bootf(res, f3)
  data.frame(N, "value" = f3(res), "lower" = bounds[1], "upper" = bounds[2], group = 1)
}

simOrig <- function(N, exact = TRUE, scaled = TRUE){
  dat <- getData(N)
  u <- dat$u
  u <- 1:(2 * N)
  l <- dat$l
  nperm <- 10 * N
  res <- llply(1:nperm, function(i){
    u <- sample(u)
    T <- computeT(u, l)
    if(exact){
      diff <- getTpminusT(u, l, N, p)
    } else {
      diff <- getTpminusTMC(u, l, N, p)
    }
    R <- N / 2 * diff + T
    c(mean(abs(diff))^3, mean(diff^2), T, mean(T * R), mean(R))
  }, .parallel = TRUE)
  res <- do.call(rbind, res)
  f1 <- function(x, ind = 1:length(x)) (2 * pi)^(-1 / 4) * sqrt(mean(x[ind]) * N / 2) * if(scaled) N^(1 / 4) else 1
  f2 <- function(x, ind = 1:length(x)) 4 * N * sqrt(var(x[ind])) * if(scaled) N else 1
  f3 <- function(x, ind = 1:length(x)) abs(var(x[ind]) - 1) * if(scaled) N else 1
  f4 <- function(x, ind = 1:length(x)) mean(abs(x[ind])) * if(scaled) N^(1 / 2) else 1
  f5 <- function(x, ind = 1:length(x)) mean(abs(x[ind])) * if(scaled) N^(1 / 2) else 1
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
  labels <- c("$(2\\pi)^{-1/4}\\sqrt{\\frac{\\mathbb{E}|T\'_{\\Pi}-T_{\\Pi}|^3}{\\lambda}}N^{1/4}\\quad $",
              "$\\frac{1}{2\\lambda}\\sqrt{\\mathrm{Var}(\\mathbb{E}[(T\'_{\\Pi}-T_{\\Pi})^2|T_{\\Pi}])}N\\quad $",
              "$|\\mathbb{E}T_{\\Pi}^2-1|N\\quad $",
              "$\\mathbb{E}|T_{\\Pi}R_{\\Pi}|N^{1/2}\\quad $",
              "$\\mathbb{E}|R_{\\Pi}|N^{1/2}\\quad $",
              "Sum of Bounds")
  data.frame(N, "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = factor(labels, levels = labels))
}

getDelta <- function(u){
  N <- length(u) / 2
  getT <- function(u, l) as.vector(t.test(u[l], u[-l], var.equal = TRUE)$statistic)
  sorted <- sort(u)
  T <- getT(sorted, 1:N)
  Tprime <- getT(sorted, c(2:N, 2*N))
  abs(T - Tprime)
}

simBetterBound <- function(N, exact = TRUE, scaled = TRUE){
  dat <- getData(N)
  u <- dat$u
  u <- 1:(2 * N)
  l <- dat$l
  delta <- getDelta(u)
  nperm <- 10 * N
  res <- llply(1:nperm, function(i){
    u <- sample(u)
    T <- computeT(u, l)
    if(exact){
      diff <- getTpminusT(u, l, N, p)
    } else {
      diff <- getTpminusTMC(u, l, N, p)
    }
    R <- N / 2 * diff + T
    c(mean(abs(diff))^3, mean(diff^2), T, mean(T * R), mean(R))
  }, .parallel = TRUE)
  res <- do.call(rbind, res)
  f0 <- function(x, ind = 1:length(x)) N / 2 * .41 * delta^3 * if(scaled) N^(1 / 2) else 1
  f1 <- function(x, ind = 1:nrow(x)) 3 * delta * (sqrt(var(x[, 1][ind])) + mean(abs(x[, 2][ind]))) * if(scaled) N else 1
  f2 <- function(x, ind = 1:length(x)) 4 * N * sqrt(var(x[ind])) * if(scaled) N else 1
  f3 <- function(x, ind = 1:length(x)) abs(var(x[ind]) - 1) * if(scaled) N else 1
  f4 <- function(x, ind = 1:length(x)) mean(abs(x[ind])) * if(scaled) N^(1 / 2) else 1
  f5 <- function(x, ind = 1:length(x)) mean(abs(x[ind])) * if(scaled) N^(1 / 2) else 1
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
  labels <- c("$\\frac{.41\\delta^3}{\\lambda}N^{1/2}\\quad $",
              "$3\\delta(\\sqrt{\\mathbb{E}T^2}+\\mathbb{E}|R|)N\\quad $",
              "$\\frac{1}{2\\lambda}\\sqrt{\\mathrm{Var}(\\mathbb{E}[(T\'_{\\Pi}-T_{\\Pi})^2|T_{\\Pi}])}N\\quad $",
              "$|\\mathbb{E}T_{\\Pi}^2-1|N\\quad $",
              "$\\mathbb{E}|T_{\\Pi}R_{\\Pi}|N^{1/2}\\quad $",
              "$\\mathbb{E}|R_{\\Pi}|N^{1/2}\\quad $",
              "Sum of Bounds")
  data.frame(N, "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = factor(labels, levels = labels))
}

