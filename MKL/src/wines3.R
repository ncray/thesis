library(RJSONIO)
library(plyr)
library(stringr)

## wines <- fromJSON("klwines.json")
## dat <- ldply(wines, function(l) data.frame(t(unlist(l)), stringsAsFactors = FALSE))
## save(dat, file = "wine.dat")
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

##Sys.setlocale('LC_ALL','C') 
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

winePower <- function(n, mkl_norm = 2, C = .01){
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

system.time(res1 <- ldply(rep(c(seq(50, 200, 50)), each = 100), winePower, .parallel = TRUE))

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

winePower <- function(n, mkl_norm = 2, C = .01){
  dat <- getWineData(n, dat3)
  rej <- with(dat, reject(compute(trainMKL), parametric = TRUE, verbose = TRUE)(desc = desc, subregion = subregion, price = price, ratings = ratings, year = year, title = NULL, l = l, desc.v = 4, title.v = NULL, mkl_norm = mkl_norm, C = C))
  wts <- data.frame(getMKLWeights())
  names(wts) <- c("desc.K3", "subr.C", "price.R1", "ratings.R1", "year.R1")
  wts$reject <- rej
  wts$n <- n
  wts$C <- C
  wts$group <- "Desc./SubR./Price/Ratings/Year"
  wts
}

system.time(res2 <- ldply(rep(c(seq(50, 200, 50)), each = 100), winePower, .parallel = TRUE))
##qplot(x = variable, y = value, geom = "boxplot", data = melt(res2[, ], id.vars = c("n", "group")))

## winePower <- function(n, mkl_norm = 2, C = .01){
##   dat <- getWineData(n, dat3)
##   rej <- with(dat, reject(compute(trainMKL), parametric = TRUE, verbose = TRUE)(desc = desc, subregion = NULL, price = price, ratings = ratings, year = year, title = NULL, l = l, desc.v = 2, title.v = NULL, mkl_norm = mkl_norm, C = C))
##   wts <- data.frame(getMKLWeights())
##   names(wts) <- c("desc.K3", "price.L", "ratings.L", "year.L")
##   wts$reject <- rej
##   wts$n <- n
##   wts$C <- C
##   wts$group <- "Desc./Price/Ratings/Year"
##   wts
## }

winePower <- function(n, mkl_norm = 2, C = .01){
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

#system.time(res3.2 <- ldply(rep(c(.5, 1, 5, 10, 100), each = 5), function(C) winePower(500, 2, C), .parallel = TRUE))
#system.time(res3.2 <- ldply(rep(as.numeric(1:8), each = 100), function(C) winePower(350, 2, C), .parallel = TRUE))
##     n C   V1
## 1 350 1 0.08
## 2 350 2 0.10
## 3 350 3 0.34
## 1 350 4 0.28
## 2 350 5 0.26
## 3 350 6 0.18
## 4 350 7 0.20
## 5 350 8 0.12

##new
##     n C   V1
## 1 350 1 0.04
## 2 350 2 0.14
## 3 350 3 0.34
## 4 350 4 0.28
## 5 350 5 0.34
## 6 350 6 0.28
## 7 350 7 0.36
## 8 350 8 0.44
#ddply(res3.2, .(n, C), function(df) mean(as.logical(df$reject)))

system.time(res3 <- ldply(rep(c(seq(50, 200, 50)), each = 100), winePower, .parallel = TRUE))

##qplot(x = variable, y = value, geom = "boxplot", data = melt(res3[, ], id.vars = c("n", "group")))

##res.m <- melt(res3[, ], id.vars = c("n", "group"))
##qplot(x = variable, y = value, geom = "boxplot", data = res.m)

winePower <- function(n, mkl_norm = 2, C = .01){
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

system.time(res4 <- ldply(rep(c(seq(50, 200, 50)), each = 100), winePower, .parallel = TRUE))
##qplot(x = variable, y = value, geom = "boxplot", data = melt(res4[, ], id.vars = c("n", "group")))
## system.time(res4.2 <- mdply(expand.grid(n = rep(c(seq(50, 200, 50)), each = 50),
##                                         C = c(.001, .01, .1, 1, 3)), winePower, .parallel = TRUE))


winePower <- function(n, mkl_norm = 2, C = .01){
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

system.time(res5 <- ldply(rep(c(seq(50, 200, 50)), each = 100), winePower, .parallel = TRUE))
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

p1 <- ggplot(res.b, aes(x = n, y = value, color = group, linetype = group)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 7))
p1
  ##facet_grid(p~.) +
  ##ggtitle("Power (Christmas Star Example)") +
  ##xlab("Radius of Outer Star (Inner is 4)") +
  ##ylab("Power")

qplot(x = n, y = V1, data = ddply(res, .(n, group), function(df) mean(df$reject.t)), geom = "line", color = group)
##qplot(x = n, y = V1, data = subset(ddply(res, .(n, group), function(df) mean(df$reject.t)), n < 350), geom = "line", color = group)
##ggsave("wine_power.png")

##system.time(res3.2 <- ldply(rep(as.numeric(3), each = 200), function(C) winePower(350, 2, C), .parallel = TRUE))
##mean(as.logical(res3.2$reject.t))
##[1] 0.21


qplot(x = n, y = V1, data = ddply(res, .(n), function(df) mean(df$reject.t)), geom = "line")
res.m <- melt(res3[, ], id.vars = c("n", "group"))
qplot(x = variable, y = value, geom = "boxplot", data = res.m)




MKLwtsWine <- function(n, C){
  print(unlist(as.list(environment())))
  RBF.v <- round(10^(seq(.5, 2, .5)), 2)
  string.v <- 1:3
  dat <- getWineData(n = n, dat3)
  l <- dat$l
  u1 <- t(dat$price)
  u2 <- dat$desc
  
  trainMKL(u1 = u1, u2 = u2, l = l, RBF.v = RBF.v, string.v = string.v, C, mkl_norm = 1, linear = FALSE)
  wts <- getMKLWeights()
  df1 <- cbind(data.frame("n" = n, "C" = C, "perm" = 0, "mkl_norm" = 1), matrix(wts, nrow = 1))
  trainMKL(u1 = u1, u2 = u2, l = l, RBF.v = RBF.v, string.v = string.v, C, mkl_norm = 2, linear = FALSE)
  wts <- getMKLWeights()
  wts <- wts / sum(wts)
  df2 <- cbind(data.frame("n" = n, "C" = C, "perm" = 0, "mkl_norm" = 2), matrix(wts, nrow = 1))
  dfperm <- ldply(1:Nwts, function(x){
    trainMKL(u1 = u1, u2 = u2, l = sample(l), RBF.v = RBF.v, string.v = string.v, C, mkl_norm = 1, linear = FALSE)
    wts <- getMKLWeights()
    df1p <- cbind(data.frame("n" = n, "C" = C, "perm" = 1, "mkl_norm" = 1), matrix(wts, nrow = 1))
    trainMKL(u1 = u1, u2 = u2, l = sample(l), RBF.v = RBF.v, string.v = string.v, C, mkl_norm = 2, linear = FALSE)
    wts <- getMKLWeights()
    wts <- wts / sum(wts)
    df2p <- cbind(data.frame("n" = n, "C" = C, "perm" = 1, "mkl_norm" = 2), matrix(wts, nrow = 1))
    rbind(df1p, df2p)
  })
  df1 <- rbind(df1, df2, dfperm)
  len <- ncol(df1)
  names(df1)[5:(5 + length(RBF.v) - 1)] <- paste("rbf: ", RBF.v, sep = "")
  names(df1)[(5 + length(RBF.v)):len] <- paste("sk: ", string.v, sep = "")
  df1
}

MKLwtsWinePlot <- function(){
  Nwts <- 100
  system.time(res <- mdply(expand.grid(n = c(10, 50, 200), C = .1), MKLwtsWine))

  res.m <- melt(res, id.vars = c(1:4))
  p1 <- qplot(variable, value, data = subset(res.m, perm == 1), geom = "boxplot") +
    facet_grid(mkl_norm~n) +
      geom_point(data = subset(res.m, perm == 0), color = "red", size = 3) +
        xlab("Kernels") +
          ylab("Kernel Weights") +
            ggtitle("Boxplot of Null Distribution with Observed in Red Faceted by Self Transition Probability and MKL Norm")
  p1
  myplot(p1, "mkl_weights_wine.png")
}

powerWine <- function(n, C){
  print(unlist(as.list(environment())))
  RBF.v <- round(10^(seq(.5, 2, .5)), 2)
  string.v <- 1:3
  ldply(1:Npwr, function(x){
    print(x)
    dat <- getWineData(n = n, dat3)
    l <- dat$l
    u1 <- t(dat$price)
    u2 <- dat$desc

    dfMKL <- data.frame("n" = n, "C" = C,
                        "FSMKL: 1" = reject(compute(trainMKL), parametric = TRUE)
                        (u1 = u1, u2 = u2, l = l, RBF.v = RBF.v, string.v = string.v, mkl_norm = 1, C = C, linear = FALSE),
                        "FSMKL: 2" = reject(compute(trainMKL), parametric = TRUE)
                        (u1 = u1, u2 = u2, l = l, RBF.v = RBF.v, string.v = string.v, mkl_norm = 2, C = C, linear = FALSE))
    dfRBF <- as.data.frame(matrix(laply(RBF.v, function(r) reject(compute(trainRBF), parametric = TRUE)
                                        (u = u1, l = l, r = r, C = C)), nrow = 1))
    names(dfRBF) <- paste("RBF: ", RBF.v, sep = "")
    dfSK <- as.data.frame(matrix(laply(string.v, function(order) reject(compute(trainString), parametric = TRUE)
                                       (u = u2, l = l, order = order, C = C)), nrow = 1))
    names(dfSK) <- paste("SK: ", string.v, sep = "")
    cbind(dfMKL, dfRBF, dfSK)
  }, .parallel = parallel)
}

powerWinePlot <- function(){
  Npwr <- 100
  ##system.time(res <- mdply(expand.grid(r1 = c(4, 4.3), self = c(.25, .35, .45), n = 50, C = .1), powerDNAStar))
  system.time(res <- mdply(expand.grid(n = c(10, 20, 30), C = .1), powerWine))

  ##system.time(res <- mdply(expand.grid(r1 = c(4.3), self = c(.335), n = 50, C = .1), powerDNAStar))
  ##colMeans(res)
  
  res2 <- ddply(res, .(n, C), function(df){
    ldply(names(df)[-(1:2)], function(name){
      dat <- df[, name]
      lims <- c(.025, .975)
      bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
      boot <- bootM(dat)
      data.frame("value" = mean(dat), "lower" = boot[1], "upper" = boot[2], "group" = name)
    })
  })

  p2 <- ggplot(res2, aes(x = n, y = value, color = group, linetype = group)) +
    geom_line() + 
      geom_errorbar(aes(ymin = lower, ymax = upper, width = .005)) +
        ##facet_grid(n~.) +
          ggtitle("Power (Christmas Star + DNA Example), Faceted on Outer Radius") +
            xlab("Self Transition Probability") +
              ylab("Power")
  p2
  myplot(p2, "wine_power.png")
}



trainCustom <- function(){
  l <- as.numeric(c(rep(-1, 5), rep(1, 5)))
  ##x <- c((-5:-1), (1:5))
  x <- 1:10
  sg('clean_kernel')
  sg('clean_features', 'TRAIN')
  sg('set_kernel', 'CUSTOM', outer(x, x), "FULL")
  sg('set_labels', 'TRAIN', l)
  sg('new_classifier', 'LIBSVM')
  sg('c', 1)
  sg('svm_use_bias', TRUE) ##default is TRUE
  ##sg('get_kernel_matrix', 'TRAIN')
  sg('train_classifier')
  mar <- getMargins(l)
  (mar - mean(mar)) / 11 + 5.5
}
