generateStar <- function(r1 = 4.1, r2 = 4, n){
  #n <- 50 ## number of observations / data points (sum for train and test and both classes)
  k_star <- 20 ## number of "leaves" of the stars
  alpha <- 0.15 ## noise level of the data

  dummy <- matrix(0, 2, 4 * n)
  dummy[1,] <- runif(4 * n)
  noise <- alpha * rnorm(4 * n)

  dummy[2, ] <- sin(k_star * pi * dummy[1, ]) + noise ## sine
  dummy[2, 1:(2 * n)] <- dummy[2, 1:(2 * n)] + r1 ## distanz shift: first class
  dummy[2, (2 * n + 1):dim(dummy)[2]] <- dummy[2, (2 * n + 1):dim(dummy)[2]] + r2 ## distanz shift: second class   

  dummy[1,] <- 2 * pi * dummy[1, ]

  x <- matrix(0, dim(dummy)[1], dim(dummy)[2])
  x[1, ] <-  dummy[2, ] * sin(dummy[1, ])
  x[2, ] <-  dummy[2, ] * cos(dummy[1, ])

  train_y <- c(-matrix(1, 1, n), matrix(1, 1, n))
  test_y <- c(-matrix(1, 1, n), matrix(1, 1, n))

  train_x <- matrix(0, 0, seq(1, dim(x)[2] / 2))
  train_x <- x[, seq(1, dim(x)[2], 2)]
  test_x  <- x[, seq(2, dim(x)[2], 2)]
  list("u" = train_x, "l" = train_y)
  ##returns 2 x (2n) matrix as u, with first n r1 and second n r2
}

getTransMat <- function(self = .25){
  m <- matrix((1 - self) / 3, nrow = 4, ncol = 4)
  diag(m) <- self
  m
}

generateMC <- function(self){
  alphabet <- c("A", "G", "T", "C")  
  N <- rpois(1, lambda)
  trans.mat <- getTransMat(self)
  stat.dist <- Re(eigen(t(trans.mat))$vectors[,1])
  stat.dist <- stat.dist / sum(stat.dist)
  res <- rep(0, N)
  res[1] <- sample(x = 1:length(alphabet), size = 1, prob = stat.dist)
  for(i in 2:N){
    res[i] <- sample(x = 1:length(alphabet), size = 1, prob = trans.mat[res[i-1], ])
  }
  paste(alphabet[res], sep = "", collapse = "")
}
