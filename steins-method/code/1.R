library(combinat)
n <- 6
perms <- permn(1:n)
A <- matrix(rnorm(n^2), nrow = n)
A <- sweep(A, 1, apply(A, 1, mean))
A <- sweep(A, 2, apply(A, 2, mean))
A <- A / sqrt(sum(A^2)) * sqrt(n - 1)

calcW <- function(perm) sum(diag(A[, perm]))
qqnorm(unlist(lapply(perms, calcW)));abline(0,1)
mean(unlist(lapply(perms, calcW))^4-1)
sum(A^4)

calcX <- function(perm) sum(diag(A[, perm]))^2-sum(diag(A[, perm])^2)
mean(unlist(lapply(perms, calcX)))
var(unlist(lapply(perms, calcX)))

calcY <- function(perm) sum(diag(A[, perm] %*% A[, perm])) - sum(diag(A[, perm])^2)
mean(unlist(lapply(perms, calcY)))
var(unlist(lapply(perms, calcY)))
plot(unlist(lapply(perms, calcX)), unlist(lapply(perms, calcY)))
plot(sort(unlist(lapply(perms, calcX))), sort(unlist(lapply(perms, calcY))))

calcY <- function(perm) sum(diag(A[, perm] %*% A[, perm]))
unlist(lapply(perms, calcY))

diag(A[, perms[[1]]] %*% t(A[, perms[[1]]]))
calcY <- function(perm) c(A[, perm][1, 2] * A[, perm][2, 1], A[, perm][3, 4] * A[, perm][4, 3])
calcY <- function(perm) c(A[, perm][1, 2] * A[, perm][2, 1], A[, perm][1, 4] * A[, perm][4, 1])
dat <- do.call(rbind, lapply(perms, calcY))
plot(dat[, 1], dat[, 2])
cov(dat[, 1], dat[, 2])

calcY <- function(perm) as.vector(A[, perm]) * as.vector(t(A[, perm]))
calcY <- function(perm) as.vector(A[, perm] - diag(A[, perm])) * as.vector(t(A[, perm] - diag(A[, perm])))
dat <- do.call(rbind, lapply(perms, calcY))
sum(apply(dat, 2, var))
