setwd("~/Dropbox/VMshare/thesis/MKL/src")
imgDir <- "../img/"
library(RJSONIO)
library(plyr)
library(stringr)
library(tikzDevice)

load("wine.dat")

getPercOkay <- function(dat) round(sort(aaply(dat, 2, function(x) sum(!is.na(x))), decreasing = TRUE) / nrow(dat) * 100, 2)

getPercOkay(dat)
sort(table(dat$info.Varietal), decreasing = TRUE)
sort(table(dat$info.Specific.Appellation.), decreasing = TRUE)

var1 <- "Pinot Noir"
var2 <- "Chardonnay"

dat2 <- subset(dat, info.Varietal. == var1 | info.Varietal. == var2)
sum(!is.na(dat2$ratings.RP) | !is.na(dat2$ratings.ST) | !is.na(dat2$ratings.BH))

dat3 <- subset(dat2, !is.na(ratings.RP) & !is.na(info.Sub.Region.))
dat3 <- dat3[, names(which(aaply(dat3, 2, function(x) sum(is.na(x))) == 0))]

dat3$ratings.RP <- as.numeric(str_match(dat3$ratings.RP, "[0-9]*"))
dat3$price <- as.numeric(str_replace_all(dat3$price, ",", ""))
dat3$desc <- str_replace_all(tolower(iconv(dat3$desc, "WINDOWS-1252", "UTF-8")), "[^a-z ]*", "")
dat3$year <- as.numeric(str_extract(dat3$title, "[0-9]{4}"))
dat3$title <- str_replace_all(tolower(iconv(dat3$title, "WINDOWS-1252", "UTF-8")), "[^a-z ]*", "")

dat3 <- subset(dat3, !is.na(year))
dat3 <- dat3[which(laply(dat3$desc, str_length) > 100), ]
dat3 <- subset(dat3, price < 600)

## years <- as.numeric(str_extract(dat3$title, "[0-9]{4}"))
## t.test(years[which(dat3$info.Varietal. == var2)], years[which(dat3$info.Varietal. == var1)])

getWineData <- function(n, dat){
  w1.ind <- sample(which(dat3$info.Varietal. == var1), size = n)
  w2.ind <- sample(which(dat3$info.Varietal. == var2), size = n)
  inds <- c(w1.ind, w2.ind)
  l <- c(-matrix(1, 1, n), matrix(1, 1, n))
  desc <- dat$desc[inds]

  subregion <- dat$info.Sub.Region.[inds]
  price <- dat$price[inds]
  ratings <- dat$ratings.RP[inds]
  year <- dat$year[inds]
  title <- dat$title[inds]
  list(desc = desc, subregion = subregion, price = price, ratings = ratings, year = year, title = title, l = l)
}

source("./generate_data.R")
source("./MKL_code.R")
source("../../friedman-test/src/twosample.R")
library(ggplot2)
library(doMC)
registerDoMC(4)
parallel <- TRUE

library(proxy)
trainMKL <- function(desc = NULL, subregion = NULL, price = NULL, ratings = NULL, year = NULL, title = NULL, l = NULL,
                     desc.v = NULL, title.v = NULL, 
                     mkl_norm = 2, C = .1, ...){
  dump <- sg('clean_kernel')
  dump <- sg('clean_features', 'TRAIN')
  sg('clean_preproc')
  if(length(desc.v) > 0){
    for(i in desc.v){
      dump <- sg('add_features','TRAIN', desc, 'RAW')
      dump <- sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', i, i - 1, 0, 'n')
      dump <- sg('add_preproc', 'SORTULONGSTRING')
      dump <- sg('attach_preproc', 'TRAIN')
    }
  }

  if(length(subregion) > 0) {
    x <- factor(subregion)
    subr <- as.matrix(dist(as.numeric(x), function(x, y) as.numeric(x == y), upper = TRUE, diag = TRUE))
    diag(subr) <- 1
  }
  
  if(length(price) > 0) dump <- sg('add_features', 'TRAIN', t(price))
  if(length(ratings) > 0) dump <- sg('add_features', 'TRAIN', t(ratings))
  if(length(year) > 0) dump <- sg('add_features', 'TRAIN', t(year))

  if(length(title.v) > 0){
    for(i in title.v){
      dump <- sg('add_features','TRAIN', title, 'RAW')
      dump <- sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', i, i - 1, 0, 'n')
      dump <- sg('add_preproc', 'SORTULONGSTRING')
      dump <- sg('attach_preproc', 'TRAIN')
    }
  }

  dump <- sg('set_labels','TRAIN', as.numeric(as.character(l)))
  dump <- sg('new_classifier', 'MKL_CLASSIFICATION')
  dump <- sg('mkl_parameters', mkl_eps, mkl_C, mkl_norm)
  dump <- sg('svm_epsilon', svm_eps)
  dump <- sg('set_kernel', 'COMBINED', 100)
  if(length(desc.v) > 0){
    for(i in desc.v){
      dump <- sg('add_kernel', 1, 'COMMSTRING', 'ULONG', cache_size, FALSE, 'FULL') ###NO,SQRT,LEN,SQLEN,FULL
    }
  }
  
  if(length(subregion) > 0) dump <- sg('add_kernel', 1, 'CUSTOM', subr, "FULL")

  ## if(length(price) > 0) dump <- sg('add_kernel', 1, 'LINEAR', 'REAL', cache_size)
  ## if(length(ratings) > 0) dump <- sg('add_kernel', 1, 'LINEAR', 'REAL', cache_size)
  ## if(length(year) > 0) dump <- sg('add_kernel', 1, 'LINEAR', 'REAL', cache_size)

  if(length(price) > 0) dump <- sg('add_kernel', 1, 'GAUSSIAN', 'REAL', cache_size, 1)
  if(length(ratings) > 0) dump <- sg('add_kernel', 1, 'GAUSSIAN', 'REAL', cache_size, 1)
  if(length(year) > 0) dump <- sg('add_kernel', 1, 'GAUSSIAN', 'REAL', cache_size, 1)

  if(length(title.v) > 0){
    for(i in title.v){
      dump <- sg('add_kernel', 1, 'COMMSTRING', 'ULONG', cache_size, FALSE, 'FULL') ###NO,SQRT,LEN,SQLEN,FULL
    }
  }

  dump <- sg('c', C)
  dump <- sg('set_kernel_normalization', 'SQRTDIAG') ##IDENTITY|AVGDIAG|SQRTDIAG|FIRSTELEMENT|VARIANCE|ZEROMEANCENTER
  ##dump <- sg('set_kernel_normalization', 'IDENTITY') ##IDENTITY|AVGDIAG|SQRTDIAG|FIRSTELEMENT|VARIANCE|ZEROMEANCENTER
  dump <- sg('train_classifier')
}

nvec <- rep(c(seq(50, 250, 50)), each = 200)

winePower1 <- function(n, mkl_norm = 2, C = .01){
  dat <- getWineData(n, dat3)
  rej <- with(dat, reject(compute(trainMKL), parametric = TRUE, verbose = TRUE)(desc = desc, subregion = subregion, price = price, ratings = ratings, year = year, title = title, l = l, desc.v = 4, title.v = 4, mkl_norm = mkl_norm, C = C))
  wts <- data.frame(getMKLWeights())
  names(wts) <- c("desc.K4", "subr.C", "price.R1", "ratings.R1", "year.R1", "title.K4")
  wts$reject <- rej
  wts$n <- n
  wts$C <- C
  wts$group <- "Desc./SubR./Price/Ratings/Year/Title"
  wts
}

system.time(res1 <- ldply(nvec, winePower1, .parallel = TRUE))

##qplot(x = variable, y = value, geom = "boxplot", data = melt(res1[, ], id.vars = c("n", "group")))

## system.time(res1 <- mdply(expand.grid(n = rep(c(seq(50, 200, 50)), each = 100),
##                                       C = c(.001, .01, .1, 1)), winePower, .parallel = TRUE))
## res.b <- ddply(res1, .(n, group, C), function(df){
##   x <- df$reject
##   lims <- c(.025, .975)
##   bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
##   boot <- bootM(x)
##   data.frame("value" = mean(x), "lower" = boot[1], "upper" = boot[2], "group" = df$group[1])
##   })
## qplot(x = n, y = value, data = res.b, geom = "line", color = factor(C)) +
##   geom_errorbar(aes(ymin = lower, ymax = upper, width = 7)) 

winePower2 <- function(n, mkl_norm = 2, C = .01){
  dat <- getWineData(n, dat3)
  rej <- with(dat, reject(compute(trainMKL), parametric = TRUE, verbose = TRUE)(desc = desc, subregion = subregion, price = price, ratings = ratings, year = year, title = NULL, l = l, desc.v = 4, title.v = NULL, mkl_norm = mkl_norm, C = C))
  wts <- data.frame(getMKLWeights())
  names(wts) <- c("desc.K4", "subr.C", "price.R1", "ratings.R1", "year.R1")
  wts$reject <- rej
  wts$n <- n
  wts$C <- C
  wts$group <- "Desc./SubR./Price/Ratings/Year"
  wts
}

system.time(res2 <- ldply(nvec, winePower2, .parallel = TRUE))
##qplot(x = variable, y = value, geom = "boxplot", data = melt(res2[, ], id.vars = c("n", "group")))

winePower3 <- function(n, mkl_norm = 2, C = .01){
  dat <- getWineData(n, dat3)
  rej <- with(dat, reject(compute(trainMKL), parametric = TRUE, verbose = TRUE)(desc = NULL, subregion = subregion, price = price, ratings = ratings, year = year, title = NULL, l = l, desc.v = NULL, title.v = NULL, mkl_norm = mkl_norm, C = C))
  wts <- data.frame(getMKLWeights())
  names(wts) <- c("subr.C", "price.R1", "ratings.R1", "year.R1")
  wts$reject <- rej
  wts$n <- n
  wts$C <- C
  wts$group <- "SubR./Price/Ratings/Year"
  wts
}

system.time(res3 <- ldply(nvec, winePower3, .parallel = TRUE))

##qplot(x = variable, y = value, geom = "boxplot", data = melt(res3[, ], id.vars = c("n", "group")))

##res.m <- melt(res3[, ], id.vars = c("n", "group"))
##qplot(x = variable, y = value, geom = "boxplot", data = res.m)

winePower4 <- function(n, mkl_norm = 2, C = .01){
  print(n)
  dat <- getWineData(n, dat3)
  rej <- with(dat, reject(compute(trainMKL), parametric = TRUE, verbose = TRUE)(desc = NULL, subregion = NULL, price = price, ratings = ratings, year = year, title = NULL, l = l, desc.v = NULL, title.v = NULL, mkl_norm = mkl_norm, C = C))
  wts <- data.frame(getMKLWeights())
  names(wts) <- c("price.R1", "ratings.R1", "year.R1")
  wts$reject <- rej
  wts$n <- n
  wts$C <- C
  wts$group <- "Price/Ratings/Year"
  wts
}

system.time(res4 <- ldply(nvec, winePower4, .parallel = TRUE))
##qplot(x = variable, y = value, geom = "boxplot", data = melt(res4[, ], id.vars = c("n", "group")))
## system.time(res4.2 <- mdply(expand.grid(n = rep(c(seq(50, 200, 50)), each = 50),
##                                         C = c(.001, .01, .1, 1, 3)), winePower, .parallel = TRUE))


winePower5 <- function(n, mkl_norm = 2, C = .01){
  print(n)
  dat <- getWineData(n, dat3)
  rej <- with(dat, reject(compute(trainMKL), parametric = TRUE, verbose = TRUE)(desc = NULL, subregion = NULL, price = price, ratings = NULL, year = year, title = NULL, l = l, desc.v = NULL, title.v = NULL, mkl_norm = mkl_norm, C = C))
  wts <- data.frame(getMKLWeights())
  names(wts) <- c("price.R1", "year.R1")
  wts$reject <- rej
  wts$n <- n
  wts$C <- C
  wts$group <- "Price/Year"
  wts
}

system.time(res5 <- ldply(nvec, winePower5, .parallel = TRUE))
##qplot(x = variable, y = value, geom = "boxplot", data = melt(res5[, ], id.vars = c("n", "group")))

res <- rbind(res1[, c("reject", "n", "group")],
             res2[, c("reject", "n", "group")],
             res3[, c("reject", "n", "group")],
             res4[, c("reject", "n", "group")],
             res5[, c("reject", "n", "group")])

res.b <- ddply(res, .(n, group), function(df){
  x <- df$reject
  lims <- c(.025, .975)
  bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
  boot <- bootM(x)
  data.frame("value" = mean(x), "lower" = boot[1], "upper" = boot[2], "group" = df$group[1])
  })

lvls <- names(table(res.b$group))[order(laply(names(table(res.b$group)), str_length), decreasing = TRUE)]

res.b$group <- factor(res.b$group, levels = lvls)
p1 <- ggplot(res.b, aes(x = n, y = value, color = group, linetype = group)) +
  geom_line(size = I(1.5)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 7), size = I(1.5)) +
  theme(legend.position = "bottom", legend.direction = "vertical")
p1
myTikz(filename = "wine_power.tex", plot = p1)
ggsave("wine_power.png")
  ##facet_grid(p~.) +
  ##ggtitle("Power (Christmas Star Example)") +
  ##xlab("Radius of Outer Star (Inner is 4)") +
  ##ylab("Power")


res.m <- rbind(melt(res1[, -which(names(res1) %in% c("reject", "n", "C"))], id.vars = "group"),
               melt(res2[, -which(names(res2) %in% c("reject", "n", "C"))], id.vars = "group"),
               melt(res3[, -which(names(res3) %in% c("reject", "n", "C"))], id.vars = "group"),
               melt(res4[, -which(names(res4) %in% c("reject", "n", "C"))], id.vars = "group"),
               melt(res5[, -which(names(res5) %in% c("reject", "n", "C"))], id.vars = "group"))
res.m$group <- factor(res.m$group, levels = lvls)
res.m$variable <- factor(res.m$variable, levels = names(table(res.m$variable))[order(table(res.m$variable))])
p2 <- qplot(x = variable, y = value, data = res.m, geom = "boxplot", color = group) +
  theme(legend.position = "bottom", legend.direction = "vertical")
myTikz(filename = "wine_weights.tex", plot = p2)
ggsave("wine_weights.png")

qplot(x = n, y = V1, data = ddply(res, .(n, group), function(df) mean(df$reject.t)), geom = "line", color = group)
##qplot(x = n, y = V1, data = subset(ddply(res, .(n, group), function(df) mean(df$reject.t)), n < 350), geom = "line", color = group)
##ggsave("wine_power.png")

##system.time(res3.2 <- ldply(rep(as.numeric(3), each = 200), function(C) winePower(350, 2, C), .parallel = TRUE))
##mean(as.logical(res3.2$reject.t))
##[1] 0.21


qplot(x = n, y = V1, data = ddply(res, .(n), function(df) mean(df$reject.t)), geom = "line")
res.m <- melt(res3[, ], id.vars = c("n", "group"))
qplot(x = variable, y = value, geom = "boxplot", data = res.m)




