##http://shogun-toolbox.org/doc/en/latest/staticcmdline.html
library("sg")
library(plyr)
library(ICSNP)
library(boot)
library(ggplot2)
library(plyr)
library(reshape)
cache_size <- 100 ##cache size per kernel in MB
svm_eps <- 1e-3
mkl_eps <- 1e-3
mkl_C <- 0

trainString <- function(u, km = NULL, l, order, C){
  gap <- 0
  reverse <- 'n'
  use_sign <- FALSE
  normalization <- 'NO' #NO,SQRT,LEN,SQLEN,FULL
  ##normalization <- 'NO' #NO,SQRT,LEN,SQLEN,FULL
  sg('clean_kernel')
  sg('clean_preproc')
  sg('clean_features', 'TRAIN')
  sg('set_features', 'TRAIN', u, 'RAW')
  sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', order, order-1, gap, reverse)
  sg('add_preproc', 'SORTULONGSTRING')
  sg('attach_preproc', 'TRAIN')
  sg('set_kernel', 'COMMSTRING', 'ULONG', cache_size, use_sign, normalization)
  sg('set_labels', 'TRAIN', as.numeric(as.character(l)))
  sg('new_classifier', 'LIBSVM')
  sg('c', C)
  sg('svm_use_bias', TRUE) ##default is TRUE
  ##sg('get_kernel_matrix', 'TRAIN')
  sg('train_classifier')
}

trainRBF <- function(u, km = NULL, l, r, C){
  sg('clean_kernel')
  sg('clean_features', 'TRAIN')
  sg('set_features', 'TRAIN', u) ##takes numeric, not integer
  ##http://www.shogun-toolbox.org/doc/en/current/classshogun_1_1CGaussianKernel.html
  sg('set_kernel', 'GAUSSIAN', 'REAL', cache_size, r)
  sg('set_labels', 'TRAIN', as.numeric(as.character(l)))
  sg('new_classifier', 'LIBSVM')
  sg('c', C)
  sg('svm_use_bias', TRUE) ##default is TRUE
  ##sg('get_kernel_matrix', 'TRAIN')
  sg('train_classifier')
}
trainLinear <- function(u, km = NULL, l, r, C){
  sg('clean_kernel')
  sg('clean_features', 'TRAIN')
  sg('set_features', 'TRAIN', u) ##takes numeric, not integer
  sg('set_kernel', 'LINEAR', 'REAL', cache_size)
  sg('set_labels', 'TRAIN', as.numeric(as.character(l)))
  sg('new_classifier', 'LIBSVM')
  sg('c', C)
  sg('svm_use_bias', TRUE) ##default is TRUE
  ##sg('get_kernel_matrix', 'TRAIN')
  sg('train_classifier')
}
trainMKL <- function(u1, u2, km = NULL, l, RBF.v = NULL, string.v = NULL, mkl_norm = 2, C = .1, linear = TRUE, ...){
  dump <- sg('clean_kernel')
  dump <- sg('clean_features', 'TRAIN')
  sg('clean_preproc')
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
  dump <- sg('set_kernel', 'COMBINED', 100)
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
getMKLWeights <- function() sg('get_subkernel_weights')
getMargins <- function(l){
  svmparams <- sg('get_svm')
  b <- as.numeric(svmparams[[1]])
  aw <- svmparams[[2]][, 1]
  inds <- svmparams[[2]][, 2] + 1 ## 0 indexing
  km <- sg('get_kernel_matrix')
  km.sub <- km[, inds]
  mar <- as.numeric(km.sub %*% (aw * as.numeric(as.character(l))[inds]) + b)  
  ##print(all((mar < 0) == (sg('classify') == -1)))
  ##plot(mar, sg('classify'))
  ##sg('classify')
  mar
}

compute <- function(train){
  function(u, km, l, ...){
    train(u = u, km = NULL, l = l, ...)
    mar <- getMargins(l)
    computeT(u = mar, l = l)
  }
}

trainMKLRBF <- function(u1, u2, km = NULL, l, RBF.v = NULL, string.v = NULL, mkl_norm = 2, C = .1, linear = TRUE, ...){
  dump <- sg('clean_kernel')
  dump <- sg('clean_features', 'TRAIN')
  sg('clean_preproc')
  if(length(RBF.v) > 0){
    for(i in 1:length(RBF.v)){dump <- sg('add_features','TRAIN', u1[i, , drop = FALSE])}
  }
  dump <- sg('set_labels','TRAIN', as.numeric(as.character(l)))
  dump <- sg('new_classifier', 'MKL_CLASSIFICATION')
  dump <- sg('mkl_parameters', mkl_eps, mkl_C, mkl_norm)
  dump <- sg('svm_epsilon', svm_eps)
  dump <- sg('set_kernel', 'COMBINED', 100)
  if(length(RBF.v) > 0){
    for(width in RBF.v){dump <- sg('add_kernel', 1, 'GAUSSIAN', 'REAL', cache_size, width)}
  }
  dump <- sg('c', C)
  dump <- sg('set_kernel_normalization', 'SQRTDIAG') ##IDENTITY|AVGDIAG|SQRTDIAG|FIRSTELEMENT|VARIANCE|ZEROMEANCENTER
  ##dump <- sg('set_kernel_normalization', 'IDENTITY') ##IDENTITY|AVGDIAG|SQRTDIAG|FIRSTELEMENT|VARIANCE|ZEROMEANCENTER
  dump <- sg('train_classifier')
}


## MKLwts <- function(r1, self, C){
##   print(paste("r1:", r1, "self:", self, "C:", C))
##   dat <- getData(r1, self)
##   u1 <- dat$u1
##   u2 <- dat$u2
##   l <- dat$l
##   trainMKL(u1, u2, l, r.v, C)
##   wts <- getMKLWeights()
##   df1 <- cbind(data.frame("r1" = r1, "self" = self, "C" = C, "perm" = 0), matrix(wts, nrow = 1))
##   df2 <- ldply(1:Nwts, function(x){  
##     trainMKL(u1, u2, sample(l), r.v, C)
##     wts <- getMKLWeights()
##     cbind(data.frame("r1" = r1, "self" = self, "C" = C, "perm" = 1), matrix(wts, nrow = 1))
##   })
##   rbind(df1, df2)
## }
## power <- function(r1, self, C){
##   print(paste("r1:", r1, "self:", self, "C:", C))
##   ldply(1:Npwr, function(x){
##     dat <- getData(r1, self)
##     u1 <- dat$u1
##     u2 <- dat$u2
##     l <- dat$l
##     print(x)
##     data.frame(cbind(data.frame("r1" = r1, "self" = self, "C" = C, "FSMKL" = rejectFSMKL(u1, u2, l, r.v, C)),
##                      matrix(laply(r.v, function(r) rejectFSRBF(u1, l, r, C)), nrow = 1),
##                      matrix(laply(1:2, function(order) rejectFSString(u2, l, order, C)), nrow = 1)))
##   })
## }

