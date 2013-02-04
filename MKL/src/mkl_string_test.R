source("MKL_code.R")
cache_size <- 50
lambda <- 100
n <- 50
C <- .5
Nwts <- 100
Npwr <- 10
r.v <- c(0.01, 0.1, 1, 10, 1000)

dat <- getData(4, .25)
x <- dat$u2[1]
library(tau)
x.c <- textcnt(x, n = 2, tolower = FALSE)
x.c <- textcnt(x, n = 1, tolower = FALSE)
x.c[names(x.c) %in% c("A", "C", "G", "T")]

num <- 200
dat <- ldply(c(rep(.25, num), rep(.5, num)), function(self){
  dat <- getData(4, self)
  x <- dat$u2[1]
  x.c <- textcnt(x, n = 1, tolower = FALSE)
  c(x.c[names(x.c) %in% c("A", "C", "G", "T")], "self" = self)
})
dat.m <- melt(dat, id.vars = "self")
qplot(variable, value, data = dat.m, geom = "boxplot", color = factor(self))

library(stringr)
str_length(x)

n <- 10
dat <- getData(10000, 1)
u2 <- dat$u2
##u2 <- c("ACTG", "AAAA", "CCCC", "TTTT", "ACTG", "AA")
##u2 <- c("ACTG", "ACGT", "GTCA", "CTGA", "ACGT")
##u2 <- c("AT", "TA", "AA", "TT")

##http://www.shogun-toolbox.org/doc/en/current/classshogun_1_1CCommUlongStringKernel.html
dump <- sg('clean_kernel')
dump <- sg('clean_features', 'TRAIN')
dump <- sg('set_features','TRAIN', u2, 'RAW')
#dump <- sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', 1, 0, 0, 'n')
dump <- sg('convert', 'TRAIN', 'STRING', 'CHAR', 'STRING', 'ULONG', 2, 1, 0, 'n')
sg('add_preproc', 'SORTULONGSTRING')
dump <- sg('attach_preproc', 'TRAIN')
##dump <- sg('set_kernel', 'COMMSTRING', 'ULONG', cache_size, FALSE, 'NO')
dump <- sg('set_kernel', 'COMMSTRING', 'ULONG', cache_size, FALSE, 'FULL')
#NO,SQRT,LEN,SQLEN,FULL
km <- sg('get_kernel_matrix', 'TRAIN')
km
image(km[, nrow(km):1])

sg('clean_kernel')
sg('clean_features', 'TRAIN')
##sg('set_features', 'TRAIN', matrix(c(1.1, 1.2, 2.0), nrow = 1)) ##takes numeric, not integer
sg('set_features', 'TRAIN', dat$u1) ##takes numeric, not integer
sg('set_kernel', 'GAUSSIAN', 'REAL', 10, 1000)
sg('set_labels', 'TRAIN', l)
sg('new_classifier', 'LIBSVM')
sg('c', C)
sg('svm_use_bias', TRUE) ##default is TRUE
km <- sg('get_kernel_matrix', 'TRAIN')
km
image(km[, nrow(km):1])


generateMC <- function(self){
  alphabet <- c("A", "G", "T", "C")
  x <- rep(alphabet, each = 10)
  if(self == .25) return(paste(sample(x), collapse = ""))
  return(paste(x, collapse = ""))
}

system.time(res <- mdply(expand.grid(r1 = c(4, 8), self = seq(.25, .35, .1), C = .5), MKLwts))




power <- function(lambda, self, C){
  print(paste("lambda:", lambda, "self:", self, "C:", C))
  ldply(1:Npwr, function(x){
    dat <- getData(4, self)
    u2 <- dat$u2
    l <- dat$l
    print(x)
    data.frame(cbind(data.frame("lambda" = lambda, "self" = self, "C" = C,
                                matrix(laply(1:3, function(order) rejectFSString(u2, l, order, C)), nrow = 1))))
  })
}

system.time(res <- mdply(expand.grid(lambda = c(50, 200, 500), self = seq(.25, .5, .05), C = .5), power))
names(res) <- c("lambda", "self", "C", paste("sk", 1:3, sep = ": "))
res2 <- ddply(res, .(lambda, self, C), function(df){
  means <- apply(df[, -(1:3)], 2, mean)
  lims <- c(.025, .975)
  bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
  bounds <- rbind(bootM(df[, "sk: 1"]),
                  bootM(df[, "sk: 2"]),
                  bootM(df[, "sk: 2"])                                    
                  )
  data.frame("lambda" = df$lambda[1], "self" = df$self[1], "C" = df$C[1], 
             "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = names(df[-(1:3)]))
})

p3 <- ggplot(res2, aes(x = group, y = value, fill = group)) +
  geom_bar() +
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .5)) +  
  facet_grid(self~lambda) +
  opts(title = "Power (Christmas Star + DNA Example, C = .5, n = 50)") +
  xlab("Kernel") +
  ylab("Power")
p3
