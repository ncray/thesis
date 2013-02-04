library(RJSONIO)
library(plyr)
library(stringr)

## wines <- fromJSON("klwines.json")
## dat <- ldply(wines, function(l) data.frame(t(unlist(l)), stringsAsFactors = FALSE))
## save(dat, file = "wine.dat")
load("wine.dat")

sort(aaply(dat, 2, function(x) sum(is.na(x)))) / nrow(dat)

sort(table(dat$info.Varietal), decreasing = TRUE)

dat2 <- subset(dat, info.Varietal. == "Zinfandel" | info.Varietal. == "Shiraz/Syrah")
sort(aaply(dat2, 2, function(x) sum(is.na(x)))) / nrow(dat2)

dat3 <- dat2[which(!is.na(dat2$ratings.RP) & !is.na(dat2$info.Sub.Region.)), ]
dat3 <- dat3[, names(which(aaply(dat3, 2, function(x) sum(is.na(x))) == 0))]

##Sys.setlocale('LC_ALL','C') 
dat3$ratings.RP <- as.numeric(str_match(dat3$ratings.RP, "[0-9]*"))
dat3$price <- as.numeric(gsub(",", "", dat3$price))


getWineData <- function(n, dat){
  w1.ind <- sample(which(dat3$info.Varietal. == "Zinfandel"), size = n)
  w2.ind <- sample(which(dat3$info.Varietal. == "Shiraz/Syrah"), size = n)
  inds <- c(w1.ind, w2.ind)
  l <- c(-matrix(1, 1, n), matrix(1, 1, n))
  price <- dat$price[inds]
  desc <- dat$info.Sub.Region.[inds]
  ##desc <- dat$desc[inds]
  list(price = price, desc = desc, l = l)
}

source("./generate_data.R")
source("./MKL_code.R")
source("../../friedman-test/src/twosample.R")
library(ggplot2)
library(doMC)
registerDoMC(4)
parallel <- TRUE

dat <- getWineData(100, dat3)
t.test(dat3[which(dat3$info.Varietal. == "Zinfandel"), ]$price, dat3[which(dat3$info.Varietal. == "Shiraz/Syrah"), ]$price)
t.test(dat[which(dat$l == -1)]$price, dat[which(dat$l == 1)]$price)
trainMKL(u1 = t(dat$price), u2 = dat$desc, l = dat$l, RBF.v = RBF.v, string.v = string.v, C, mkl_norm = 1, linear = TRUE)
  wts <- getMKLWeights()

reject(compute(trainMKL), parametric = FALSE, verbose = TRUE)(u1 = t(dat$price), u2 = dat$desc, l = dat$l, RBF.v = RBF.v, string.v = string.v, mkl_norm = 2, C = .01, linear = TRUE)



trainMKL <- function(u1, u2, km = NULL, l, RBF.v = NULL, string.v = NULL, mkl_norm = 2, C = .1, linear = TRUE, ...){
  dump <- sg('clean_kernel')
  dump <- sg('clean_features', 'TRAIN')
  if(linear) dump <- sg('add_features','TRAIN', u1)
  if(length(RBF.v) > 0){
    for(i in RBF.v){dump <- sg('add_features','TRAIN', u1)}
  }
  if(length(string.v) > 0){
    for(i in string.v){
      dump <- sg('add_features','TRAIN', u2, 'RAW')
      dump <- sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', i, i - 1, 0, 'n')
      dump <- sg('add_preproc', 'SORTULONGSTRING')
      dump <- sg('attach_preproc', 'TRAIN')
    }
  }
  dump <- sg('set_labels','TRAIN', as.numeric(as.character(l)))
  dump <- sg('new_classifier', 'MKL_CLASSIFICATION')
  dump <- sg('mkl_parameters', mkl_eps, mkl_C, mkl_norm)
  dump <- sg('svm_epsilon', svm_eps)
  dump <- sg('set_kernel', 'COMBINED', 0)
  if(linear) dump <- sg('add_kernel', 1, 'LINEAR', 'REAL', cache_size)
  if(length(RBF.v) > 0){
    for(width in RBF.v){dump <- sg('add_kernel', 1, 'GAUSSIAN', 'REAL', cache_size, width)}
  }
  if(length(string.v) > 0){
    for(i in string.v){
      dump <- sg('add_kernel', 1, 'COMMSTRING', 'ULONG', cache_size, FALSE, 'FULL') ###NO,SQRT,LEN,SQLEN,FULL
    }
  }
  dump <- sg('c', C)
  dump <- sg('set_kernel_normalization', 'SQRTDIAG') ##IDENTITY|AVGDIAG|SQRTDIAG|FIRSTELEMENT|VARIANCE|ZEROMEANCENTER
  ##dump <- sg('set_kernel_normalization', 'IDENTITY') ##IDENTITY|AVGDIAG|SQRTDIAG|FIRSTELEMENT|VARIANCE|ZEROMEANCENTER
  dump <- sg('train_classifier')
}

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
