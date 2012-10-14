##mogrify -resize 300x197! *.jpg
##convert -monochrome image_0001.jpg image_0001mono.jpg (don't use black and white...use grayscale instead below)
##for i in *.jpg; convert $i -colorspace Gray $i; done
##in eshell...for i in *.jpg { convert $i -colorspace Gray $i }
##for pigeons and roosters...mogrify -resize 297x300! *.jpg 
## all.path <- list.files("~/101_ObjectCategories/", full.names = TRUE, recursive = TRUE)
## all.dims <- llply(all.path, function(path) dim(readJPEG(path)))

library(jpeg)
im <- readJPEG("~/101_ObjectCategories/car_side/image_0001.jpg")
flip <- function(im) image(t(im[nrow(im):1,]), col = gray((0:32) / 32))
flip(im)

one.path <- list.files("~/101_ObjectCategories/car_side", full.names = TRUE)[1:123]
two.path <- list.files("~/101_ObjectCategories/airplanes_resized", full.names = TRUE)[1:123]

one.path <- list.files("~/101_ObjectCategories/rooster_resized", full.names = TRUE)[1:45]
two.path <- list.files("~/101_ObjectCategories/pigeon_resized", full.names = TRUE)[1:45]

##llply(list.files("~/101_ObjectCategories/rooster", full.names = TRUE)[1:45], function(path) dim(readJPEG(path)))
##llply(list.files("~/101_ObjectCategories/pigeon", full.names = TRUE)[1:45], function(path) dim(readJPEG(path)))

library(ggplot2)
ldply(one.path, function(path) dim(readJPEG(path)))
ldply(two.path, function(path) dim(readJPEG(path)))

one <- laply(one.path, function(path) as.vector(readJPEG(path)))
two <- laply(two.path, function(path) as.vector(readJPEG(path)))

flip(matrix(one[1, ], nrow = 197))

n <- nrow(one)
one.samp <- one[1:n, ]
two.samp <- two[1:n, ]

library(kernlab)
X <- rbind(one.samp, two.samp)
X <- sweep(X, 1, apply(X, 1, mean), "-")
X <- sweep(X, 1, apply(X, 1, function(vec) sqrt(sum(vec^2))), "/")
pd <- polydot(degree = 1, scale = 1, offset = 0)
pd <- polydot(degree = 3, scale = 1, offset = 0)
pd <- polydot(degree = 4, scale = 1, offset = 1)
km <- kernelMatrix(pd, X)
png(filename = "car-plane-unscaled-1.png", width = 800, height = 600)
png(filename = "rooster-pigeon-unscaled-1.png", width = 800, height = 600)
png(filename = "car-plane-scaled-1.png", width = 800, height = 600)
png(filename = "rooster-pigeon-scaled-1.png", width = 800, height = 600)
png(filename = "rooster-pigeon-scaled-2.png", width = 800, height = 600)
png(filename = "rooster-pigeon-scaled-3.png", width = 800, height = 600)
png(filename = "rooster-pigeon-scaled-4.png", width = 800, height = 600)
image(km)
dev.off()
##X <- t(scale(t(X)))
apply(X, 1, mean)
apply(X, 1, sd)
apply(X, 1, function(vec) sum(vec^2))

flip(matrix(one.samp[1, ]-X[1,], nrow = 197))
flip(matrix(X[124, ], nrow = 197))

##km <- ((km-1)^2)
image(km)
## scramb <- sample(1:nrow(km))
## image(km[scramb, scramb])
## sum((km - km[scramb, scramb])^2)
## res <- 1:100
## for(i in 1:100){
##   scramb2 <- sample(1:nrow(km))
##   res[i] <- sum((km[scramb2, scramb2] - km[scramb, scramb])^2)
## }
## plot(log(eigen(km)$values))
## plot(log(eigen(km[scramb, scramb])$values))

## plot(eigen(km)$values)
## plot(eigen(km)$vectors[,1])
## plot(eigen(km[scramb, scramb])$vectors[,1])
## plot(km[1,])
## plot(log(km[1,]))
## plot(km[101,])
## plot(log(km[101,]))
## mean(km[1,1:100])
## mean(km[1,101:200])
## mean(km[101,1:100])
## mean(km[101,101:200])

## res <- adply(km, 1, function(row) c(mean(row[1:nrow(two.samp)]), mean(row[-(1:nrow(two.samp))])))
## res$X1 <- 1:nrow(km)
## res2 <- melt(res, id.vars = "X1")
## qplot(X1, value, data = res2, geom = "point") + facet_wrap(~variable, scale = "free")

y <- c(rep(0, nrow(one.samp)), rep(1, nrow(two.samp)))
ksvm1 <- ksvm(km, y)
t.stat <- t.test(fitted(ksvm1)[y==0],fitted(ksvm1)[y==1], variance.equal = TRUE)$statistic
t.perm <- rep(0,1000)
for(i in 1:1000){
  y.perm <- sample(y)
  ksvm2 <- ksvm(km,y.perm)
  t.perm[i] <- t.test(fitted(ksvm2)[y.perm==0],fitted(ksvm2)[y.perm==1], variance.equal = TRUE)$statistic
}
mean(t.stat > t.perm)
(t.stat - mean(t.perm)) / sd(t.perm)
qplot(t.perm,geom="bar",main="Permutation Null Distribution with Observed T-statistic (Red)") + geom_vline(xintercept=t.stat,col="red")
ggsave("car-airplane-twosamp.png", width = 6.4, height = 4.8)
ggsave("rooster-pigeon-twosamp.png", width = 6.4, height = 4.8)
ggsave("rooster-pigeon-twosamp-4-1.png", width = 6.4, height = 4.8)
ggsave("rooster-pigeon-twosamp-3-0.png", width = 6.4, height = 4.8)

##rooster-pigeon check 3/1 or 4/1
library(doMC)
registerDoMC(4)
getPval <- function(deg, offset){
  pd <- polydot(degree = deg, scale = 1, offset = offset)
  km <- kernelMatrix(pd, X)
  y <- c(rep(0, nrow(one.samp)), rep(1, nrow(two.samp)))
  ksvm1 <- ksvm(km, y)
  t.stat <- t.test(fitted(ksvm1)[y==0],fitted(ksvm1)[y==1], variance.equal = TRUE)$statistic
  t.perm <- rep(0,1000)
  for(i in 1:1000){
    y.perm <- sample(y)
    ksvm2 <- ksvm(km,y.perm)
    t.perm[i] <- t.test(fitted(ksvm2)[y.perm==0],fitted(ksvm2)[y.perm==1], variance.equal = TRUE)$statistic
  }
  data.frame("deg" = deg,
             "offset" = offset,
             "p.val" = 2 * min(mean(t.stat < t.perm), mean(t.stat > t.perm)),
             "standardized" = (t.stat - mean(t.perm)) / sd(t.perm))
}

dat.in <- data.frame(deg = rep(1:5, each = 2), offset = rep(0:1, 5))
dat.in <- do.call(rbind, llply(1:5, function(x) dat.in))

dat <- mdply(dat.in, getPval, .parallel = TRUE)
qplot(deg, p.val, data = dat, geom = "point", main = "P on Deg., Faceted on Offset (5 runs)") + facet_grid(.~offset)
ggsave("rooster-pigeon-degree.png", width = 6.4, height = 4.8)
