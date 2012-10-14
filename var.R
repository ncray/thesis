getOneN <- function(n){
  print(n)
  x <- rnorm(2 * n)
  x <- x - mean(x)
  x <- x / sqrt(sum(x^2)) * sqrt(2*n)
  var(sapply(1:10000, function(a) 1 / (1 - mean(x[sample(1:(2*n), n)]))))
}
x <- sapply(10^(seq(1, 4, .5)), getOneN)

