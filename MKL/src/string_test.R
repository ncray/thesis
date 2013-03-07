#http://shogun-toolbox.org/doc/en/latest/staticcmdline.html
library("sg")


alphabet <- c("A", "G", "T", "C")
trans.mat <- matrix(c(1/2, 0  , 1/4, 1/4,
                      1/3, 1/3, 1/3, 0,
                      0  , 1  , 0  , 0,
                      1/4, 1/4, 1/4, 1/4),
                    nrow = 4, byrow = TRUE)
##apply(trans.mat, 1, sum)
stat.dist <- Re(eigen(t(trans.mat))$vectors[,1])
stat.dist <- stat.dist / sum(stat.dist)
##stat.dist %*% trans.mat
##library(expm)
##library(MCMCpack)
##as.vector(rdirichlet(1, rep(1, 4))) %*% (trans.mat %^% 10) - stat.dist
generateMC <- function(lambda, alphabet, trans.mat){
  N <- rpois(1, lambda)
  stat.dist <- Re(eigen(t(trans.mat))$vectors[,1])
  stat.dist <- stat.dist / sum(stat.dist)
  res <- rep(0, N)
  res[1] <- sample(x = 1:length(alphabet), size = 1, prob = stat.dist)
  for(i in 2:N){
    res[i] <- sample(x = 1:length(alphabet), size = 1, prob = trans.mat[res[i-1], ])
  }
  paste(alphabet[res], sep = "", collapse = "")
}

n <- 10
library(plyr)
d1 <- c(laply(rep(50, n), function(lambda) generateMC(lambda, alphabet, trans.mat)),
        laply(rep(50, n), function(lambda) generateMC(lambda, alphabet, matrix(rep(1/4, 16), nrow = 4))))
##d2 <- matrix(as.numeric(1:(2*n)), nrow = 1)
d2 <- matrix(c(rnorm(n, -10), rnorm(n, 10)), nrow = 1)

d1 <- c("abcd", "abcd", "aaaa", "aabb", "abcdabcdabcd", "abcdabcdabcd", "ababab", "aaabbb", "bababa", "ababab")
d1 <- c(laply(1:10, function(i) generateMC(.25)),
laply(1:10, function(i) generateMC(.5)))

cache_size <- 10
order <- 1
gap <- 0
reverse <- 'n'
use_sign <- FALSE
normalization <- 'NO' #NO,SQRT,LEN,SQLEN,FULL
d1 <- c("GAGAGA", "GGGAAA", "AAGGAG", "GGAAGA")

sg('clean_kernel')
sg('add_preproc', 'SORTULONGSTRING')
sg('set_kernel', 'COMMSTRING', 'ULONG', cache_size, use_sign, normalization)
sg('set_features', 'TRAIN', d1, 'RAW')
sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', order, order-1, gap, reverse)
km1 <- sg('get_kernel_matrix', 'TRAIN')
##km1
image(km1)
sum(km1[1:10, 1:10])
sum(km1[11:20, 11:20])

sum(diag(km1[1:10, 1:10]))
sum(diag(km1[11:20, 11:20]))


library(plyr)
d1 <- c(laply(1:10, function(i) generateMC(.25)),
        laply(1:10, function(i) generateMC(.45)))
order <- 1
  gap <- 0
  reverse <- 'n'
  use_sign <- FALSE
  normalization <- 'FULL' #NO,SQRT,LEN,SQLEN,FULL
  ##normalization <- 'NO' #NO,SQRT,LEN,SQLEN,FULL
  sg('clean_kernel')
  sg('clean_features', 'TRAIN')
  sg('set_features', 'TRAIN', d1, 'RAW')
  sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', order, order-1, gap, reverse)
  sg('add_preproc', 'SORTULONGSTRING')
  sg('attach_preproc', 'TRAIN')
  sg('set_kernel', 'COMMSTRING', 'ULONG', 10, use_sign, normalization)
km1 <- sg('get_kernel_matrix', 'TRAIN')
image(km1)





d1 <- c("ACTG", "ACTG", "ACACAC", "GAGAGA", "AAAA", "CCCC", "GGGAAA", "AAGGAG")
d1 <- c("GAGAGA", "GGGAAA", "AAGGAG", "GGAAGA")
sg('clean_kernel')
sg('clean_preproc')
sg('add_preproc', 'SORTWORDSTRING')
sg('attach_preproc', 'TRAIN')
sg('set_kernel', 'COMMSTRING', 'WORD', 10, FALSE, 'NO')
sg('set_features', 'TRAIN', d1, 'DNA')
sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'WORD', 1, 0, 0, 'n')
##sg('attach_preproc', 'TRAIN')
sg('get_kernel_matrix', 'TRAIN')


sg('clean_kernel')
sg('clean_features','TRAIN'); 
sg('clean_features','TEST')
sg('set_features', 'TRAIN', d1, 'DNA');
sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'WORD', 1, 0);
sg('add_preproc', 'SORTWORDSTRING')
sg('attach_preproc', 'TRAIN'); 
sg('set_kernel','COMMSTRING','WORD', 10, FALSE, 'NO');
sg('init_kernel', 'TRAIN');
sg('get_kernel_matrix')


sg('clean_features', 'TRAIN')
sg('clean_kernel')
sg('set_features', 'TRAIN', d1, 'DNA')
sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'WORD', order, order-1)
sg('add_preproc', 'SORTWORDSTRING')
sg('attach_preproc', 'TRAIN')
sg('set_labels', 'TRAIN', c(1, -1, 1, -1))
sg('new_classifier', 'SVMLIGHT')
sg('set_kernel', 'COMMSTRING', 'WORD', 10, TRUE, 'NO')
sg('c', C)
km=sg('get_kernel_matrix', 'TRAIN')


### THIS ONE WORKS YAY
sg('clean_kernel')
sg('set_features', 'TRAIN', d1, 'RAW')
sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', order, order-1, gap, reverse)
sg('add_preproc', 'SORTULONGSTRING')
sg('attach_preproc', 'TRAIN')
sg('set_kernel', 'COMMSTRING', 'ULONG', cache_size, use_sign, normalization)
sg('get_kernel_matrix', 'TRAIN')
###

sg('clean_kernel')
sg('set_features', 'TRAIN', d2) ##takes numeric, not integer
##sg('set_kernel', 'LINEAR', 'REAL', cache_size)
sg('set_kernel', 'GAUSSIAN', 'REAL', cache_size, 1)
##sg('help', 'set_kernel_normalization')
sg('set_kernel_normalization', 'IDENTITY')
##sg('set_kernel_normalization', 'ZEROMEANCENTER')
km2 <- sg('get_kernel_matrix', 'TRAIN')

sg('clean_kernel')
sg('set_kernel', 'COMBINED', 0)
sg('add_preproc', 'SORTULONGSTRING')
sg('add_features', 'TRAIN', d1, 'RAW')
sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', order, order-1, gap, reverse)
sg('attach_preproc', 'TRAIN')
#sg('help', 'add_kernel')
sg('add_kernel', .5, 'COMMSTRING', 'ULONG', cache_size, use_sign, normalization)
sg('add_features', 'TRAIN', d2) ##takes numeric, not integer
sg('set_kernel_normalization', 'IDENTITY')
sg('add_kernel', .5, 'LINEAR', 'REAL', cache_size)
##sg('send_command', 'init_kernel TRAIN')
km <- sg('get_kernel_matrix', 'TRAIN')
km1[1:4, 1:4]
km2[1:4, 1:4]
km[1:4, 1:4]


sg('clean_kernel')
sg('add_features', 'TRAIN', d2)
sg('add_features', 'TRAIN', d1, 'DNA')
# this creates a combined kernel which is
# 0.5 * linear_kernel_on realvalued features plus
# 0.5 * wd kernel on dna
sg('send_command', paste('set_kernel COMBINED', cache_size))
sg('send_command', paste('add_kernel 0.5 LINEAR REAL', cache_size))
﻿sg('send_command', paste('add_kernel 0.5 WEIGHTEDDEGREE CHAR', cache_size, degree))
sg('send_command', 'init_kernel TRAIN')
sg('get_kernel_matrix')


sg('clean_kernel')
sg('clean_features', 'TRAIN')
sg('add_features','TRAIN', d2)
sg('add_features','TRAIN', d2)
sg('add_features','TRAIN', d1, 'RAW')
sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', order, order-1, gap, reverse)
sg('attach_preproc', 'TRAIN')
sg('set_labels','TRAIN', c(rep(-1, n), rep(1, n)))         # set the labels
sg('new_classifier', 'MKL_CLASSIFICATION')
sg('mkl_parameters', 1e-3, 0)
sg('svm_epsilon', 1e-3)
sg('set_kernel', 'COMBINED', 0)
sg('add_kernel', 1, 'GAUSSIAN', 'REAL', cache_size, .1)
sg('add_kernel', 1, 'GAUSSIAN', 'REAL', cache_size, .2)
sg('add_kernel', 1, 'COMMSTRING', 'ULONG', cache_size, use_sign, normalization)
sg('c', 1)
sg('train_classifier')
sg('get_svm')
sg('get_subkernel_weights')


