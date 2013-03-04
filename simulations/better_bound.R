library(combinat)
library(plyr)
library(doMC)
registerDoMC(4)
parallel <- TRUE
setwd("better_bound_condition")

getPairedTprime <- function(T, N, u, l, x, y){
  x <- u[l]
  y <- u[-l]
  del <- rep(y, length(x)) - rep(x, each = length(y))
  xbar <- mean(x)
  ybar <- mean(y)
  Tprime <- (xbar - ybar + 2/N*del) /
    (sqrt(2/N)*sqrt(sum(u^2)/(2*(N-1)) - 1/2*N/(N-1)*(xbar^2 + ybar^2 + 2*del/N*(xbar-ybar) + 2*del^2/N^2)))
  Tprime[which.max(abs(T-Tprime))]
}

getT <- function(u, l) as.vector(t.test(u[l], u[-l], var.equal = TRUE)$statistic)

getDifferences <- function(N = 5, mu = 1){
  x <- rnorm(N, mu)
  y <- rnorm(N, 0)
  u <- c(x, y)
  
  combs <- combinat::combn(1:(2*N), N, simplify = FALSE)
  ret <- ldply(combs, .fun = function(l){
    T <- getT(u, l)
    Tprime <- getPairedTprime(T, N, u, l, x, y)
    c("T" = T, "Tprime" = Tprime, "absDiff" = abs(T - Tprime))
  })

  sorted <- sort(u)
  maxT <- max(abs(ret$T))
  shortcutT <- abs(getT(sorted, 1:N))
  if(abs(maxT - shortcutT) > 1e-10){
    print("T")
    print(c(maxT, shortcutT))
  }
  maxDiff <- max(ret$absDiff)
  shortcutDiff <- abs(getT(sorted, 1:N) - getT(sorted, c(2:N, 2*N)))
  if(abs(maxDiff - shortcutDiff) > 1e-10){
    print("diff")
    print(c(maxDiff, shortcutDiff))
  }
  
  ret
}

getDifferencesFast <- function(N = 5){
  std <- sqrt(N)
  ##std <- 1
  ##std <- 1 / sqrt(N)
  x <- rnorm(N, mean = 2, sd = std)
  y <- rnorm(N, mean = 0, sd = std)
  u <- c(x, y)
  ##u <- 1:(2*N)
  sorted <- sort(u)
  print(getT(u, 1:N))
  T <- getT(sorted, 1:N)
  Tprime <- getT(sorted, c(2:N, 2*N))
  diff <- abs(T - Tprime)
  c(T = abs(T), Tprime = abs(Tprime), diff = diff, N = N)
}

dat <- ldply(rep(floor(10^(seq(1, 2.5, by = .5))), each = 10), getDifferencesFast)
dat <- ldply(rep(seq(10, 100, 10), each = 5), getDifferencesFast)
dat <- ldply(rep(seq(100, 1000, 100), each = 5), getDifferencesFast)
dat <- ldply(rep(seq(10000, 100000, 10000), each = 5), getDifferencesFast)
dat <- ldply(rep(c(10, 100, 1000, 2500, 5000, 7500, 10000), each = 10), getDifferencesFast)
dat <- ldply(rep(seq(100000, 300000, 25000), each = 3), getDifferencesFast, .parallel = TRUE)
dat <- ldply(rep(seq(100000, 1000000, 25000), each = 3), getDifferencesFast, .parallel = TRUE)
dat <- ldply(rep(floor(10^(seq(1, 6, by = .2))), each = 5), getDifferencesFast, .parallel = TRUE)
dat <- ldply(rep(seq(1000000, 10000000, 100000), each = 1), getDifferencesFast, .parallel = TRUE)
arrange(dat, N)

dat <- ldply(rep(seq(1000, 10000, 10), each = 1), getDifferencesFast)
dat <- ldply(rep(floor(10^(seq(1, 6, by = .2))), each = 1), getDifferencesFast, .parallel = TRUE)
ggplot(dat, aes(x = N, y = diff)) + geom_point() +
  geom_smooth(method = "lm", formula = y ~ 1 + I(log(log(x))), color = "red") +
  geom_smooth(method = "lm", formula = y ~ 1 + I(log(x)), color = "blue")

ggplot(dat, aes(x = N, y = log(log(diff)))) + geom_point()
summary(lm(diff ~ -1 + I(log(log(N))), data = dat))

ggplot(dat, aes(x = N, y = diff)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ -1 + I(1/x), color = "red") +
  geom_smooth(method = "lm", formula = y ~ -1 + I(1/x^(1/2)), color = "blue") +
  geom_smooth(method = "lm", formula = y ~ -1 + I(1/x^(1/4)), color = "black") +
  ggtitle("SD: 1 / sqrt(N), Red: Fitted 1/N, Blue: Fitted 1/N^(1/2), Black: Fitted 1/N^(1/4)")
ggsave("rate_plot_3.png")

##scale_y_log10(breaks = c(.5, 1, 2, 3, 4, 5))
##geom_line(aes(y = log10(15 * N^(-1))), col = "black", linetype = 2) + 
##geom_line(aes(y = log10(15 * N^(-1/2))), col = "black", linetype = 2)

## lm(diff ~ -1 + I(1 / N), data = dat)
## lm(diff ~ -1 + I(1 / N^(1/2)), data = dat)
## lm(diff ~ -1 + I(1 / N^(1/4)), data = dat)

system.time(dat <- getDifferences(5))
arrange(dat, T)

library(ggplot2)

system.time(dat <- getDifferences(6))
system.time(dat <- mdply(expand.grid(N = 5:7, mu = c(1, 5, 10)), getDifferences, .parallel = TRUE))
dat <- ddply(dat, .(N, mu), transform, isMax = abs(absDiff - max(absDiff)) < 1e-10)
ggplot(data = dat, aes(x = T, y = Tprime, color = factor(isMax))) +
  geom_point(alpha = 1) +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(N~mu) + 
  ggtitle("Tprime on T, Faceted by N and Mu")
ggsave("t_tprime_plot.png")



getData <- function(N){
  u <- c(rnorm(N, -1), rnorm(N, 1))
  ##u <- 1:(2 * N)
  l <- c(rep(-1, N), rep(1, N))
  u <- u - mean(u)
  u <- u * sqrt(1 / sum(u^2) * 2 * N)
  list("u" = u, "l" = l)
}

getDifferencesFast <- function(N = 5){
  u <- getData(N)$u
  ##u <- rnorm(2 * N)
  ##u <- 1:(2 * N)
  sorted <- sort(u)
  ##print(getT(u, 1:N))
  T <- getT(sorted, 1:N)
  Tprime <- getT(sorted, c(2:N, 2*N))
  diff <- abs(T - Tprime)
  diff <- N / 2 * .41 * diff^3 * N^(1 / 2)
  c(T = abs(T), Tprime = abs(Tprime), diff = diff, N = N)
}

N <- 10
u <- 1:(2 * N)
sorted <- sort(u)
getT(sorted, 1:N)
f1 <- function(N) sqrt((N - 1) / N) * (1/2 * N * (N + 1) - 1/2 * N * (3 * N + 1)) / sqrt(N / 6 * (N^2 - 1))
f1(N)
getT(sorted, c(2:N, 2*N))
f2 <- function(N) sqrt((N - 1) / N) * (-N^2 + 4 * N - 2) / sqrt(2 * (N^3 / 12 + 2 * N^2 - 61/12 * N - 1 / N + 4))
f2(N)

f3 <- function(N) abs(f1(N) - f2(N))
laply(10^(1:10), f3)
getDelta(u)
f3(N)

laply(seq(1000000, 10000000, 100000), function(N) getDelta(1:(2 * N)) - f3(N))

f4 <- function(N) .41 * N / 2 * abs(f1(N) - f2(N))^3 * N^(1/2)
laply(10^(seq(1, 6, by = .5)), f4)

##http://www.wolframalpha.com/input/?i=limit+sqrt%28n%29+*+sqrt%28%28n-1%29%2Fn%29+*+abs%28n%5E2%2Fsqrt%28n%5E3%2F6-n%2F6%29+%2B+%28-n%5E2%2B4n-2%29%2Fsqrt%282+*+%28n%5E3%2F12%2B2n%5E2-61%2F12n-1%2Fn%2B4%29%29%29+as+n-%3Einfinity
