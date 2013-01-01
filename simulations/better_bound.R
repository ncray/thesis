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
  ##std <- sqrt(N)
  ##std <- 1
  std <- 1 / sqrt(N)
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

dat <- ldply(rep(c(10, 100, 1000, 2500, 5000, 7500, 10000), each = 10), getDifferencesFast)
dat <- ldply(rep(seq(100000, 300000, 25000), each = 3), getDifferencesFast, .parallel = TRUE)
dat <- ldply(rep(seq(100000, 1000000, 25000), each = 3), getDifferencesFast, .parallel = TRUE)
arrange(dat, N)

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

