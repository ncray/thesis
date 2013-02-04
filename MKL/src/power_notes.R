source("MKL_code.R")
r1 <- 4
self <- .4
C <- .5
cache_size <- 50
lambda <- 100
n <- 50
r.v <- c(0.01, 0.1, 1, 10, 1000)
dat <- getData(r1, self)
u1 <- dat$u1
u2 <- dat$u2
l <- dat$l
computeFSMKL(u1, u2, l, r.v, C)
##16.54504
sort(laply(1:19, function(i) computeFSMKL(u1, u2, sample(l), r.v, C)))
## [1]  8.728773  9.701665  9.886105 10.101656 10.610667 11.134634 13.247181
## [8] 14.236651 16.302945 20.411420 20.562508 20.886986 21.614235 22.498868
## [15] 22.618628 24.719914 28.143820 29.171920 32.178224
##poor MKL performance here

computeFSString(u2, l, 1, C)
## [1] 0.555421
sort(laply(1:19, function(i) computeFSString(u2, sample(l), 1, C)))
##  [1] 0.3649519 0.5314595 0.9279623 0.9464526 1.4266679 1.4269882 1.5652383
##  [8] 1.6040141 1.6062570 1.6795123 1.7875579 1.7954662 1.8442131 1.9568429
## [15] 1.9577881 1.9734489 2.1523404 2.2382138 2.2933613
##1-spectrum kernel has no power here

computeFSString(u2, l, 2, C)
## [1] 16.43727
sort(laply(1:19, function(i) computeFSString(u2, sample(l), 2, C)))
##  [1] 1.828907 1.919179 2.032788 2.037571 2.318295 2.341751 2.422413 2.483185
##  [9] 2.492075 2.512780 2.617341 2.713932 2.787643 2.918653 2.982431 3.230805
## [17] 3.549307 3.579574 3.888526
##but 2-spectrum has a lot

trainMKL(u1, u2, l, r.v, C)
wts <- getMKLWeights()
wts
##             [,1]         [,2]         [,3]        [,4]         [,5]
## [1,] 0.006245465 0.0008765602 0.0001058235 2.53889e-05 9.180381e-06
##              [,6]      [,7]
## [1,] 7.743008e-06 0.9927298
##the weights for MKL look great in that they completely pick out the 2-specturm kernel
##this jives with the t-statistic of 16.5

trainMKL(u1, u2, sample(l), r.v, C)
wts <- getMKLWeights()
wts
##           [,1]       [,2]         [,3]      [,4]         [,5]         [,6]
## [1,] 0.7152297 0.06961423 0.0006009961 0.2144602 2.546353e-05 2.741316e-05
##              [,7]
## [1,] 4.201437e-05
##randomizing the labels puts a lot of weight on the the .01 RBF kernel

r.v <- c(0.1, 1, 10, 1000)
computeFSMKL(u1, u2, l, r.v, C)
## [1] 17.09348
sort(laply(1:19, function(i) computeFSMKL(u1, u2, sample(l), r.v, C)))
##  [1]  7.166700  8.209707  9.069345  9.807837 10.199144 10.442660 10.771110
##  [8] 11.619869 11.951194 12.137014 12.212503 12.309588 12.404658 12.549092
## [15] 12.678796 12.793843 12.802443 13.850175 13.986753
##if we get rid of the .01 RBF kernel, the randomization-t distribution is better (shifted left)

r.v <- c(0.01, 0.1, 1, 10, 1000)
C <- .1
computeFSMKL(u1, u2, l, r.v, C)
## [1] 16.35192
sort(laply(1:19, function(i) computeFSMKL(u1, u2, sample(l), r.v, C)))
##  [1]  4.042382  7.263712  7.759659  8.015074  9.033516  9.409674  9.819325
##  [8] 10.198787 10.672614 10.914358 11.016767 11.754198 12.050675 12.176706
## [15] 12.219449 12.316027 13.561877 14.176406 17.000964
## decreasing C from .5 to .1 works pretty well too
## small C = small penalty on misclassified points (lenient, so not super wiggly)
## big C = overfit to .01 RBF kernel

## let's take a look at a whole range of C values:
C.v <- c(.001, .01, .1, .5, 1)
##C.v <- seq(.01, 5, by = .1)
res.l <- llply(C.v, function(C)laply(1:19, function(i) computeFSMKL(u1, u2, sample(l), r.v, C)))
res <- do.call(cbind, res.l)
res <- data.frame(res)
names(res) <- C.v
library(reshape)
res.m <- melt(res)
library(ggplot2)
actual <- laply(C.v, function(C) computeFSMKL(u1, u2, l, r.v, C))
qplot(variable, value, geom = "boxplot", data = res.m, xlab = "C", ylab = "T") +
  geom_point(aes(x = factor(C.v), y = actual), color = "red")
ggsave("tstat_on_c2.png")
