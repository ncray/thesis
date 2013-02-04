N <- 10
x <- rnorm(N, 0)
z <- rnorm(N, 2)
u <- c(x, z)
y <- c(rep(1, N), rep(-1, N))
t.test(x, z)
library(ggplot2)
qplot(u, 1, geom = "point", color = factor(y))
s <- fitted(lm(as.numeric(y)~u))
color <- letters[factor(y)]
means1 <- tapply(u, y, mean)
means2 <- tapply(s, y, mean)
dat <- data.frame(u, y, s, color)
ggplot(dat, aes(u, 1, color = color)) +
  geom_point(aes(x = u, y = 1), size = 3) +
  geom_point(aes(x = s, y = -1), size = 3) +
  geom_segment(aes(x = u, y = 1, xend = s, yend = -1)) +
  geom_point(aes(x = means1[1], y = 1, color = "a"), shape = 2, size = 5) + 
  geom_point(aes(x = means1[2], y = 1, color = "b"), shape = 2, size = 5) + 
  geom_point(aes(x = means2[1], y = -1, color = "a"), shape = 2, size = 5) + 
  geom_point(aes(x = means2[2], y = -1, color = "b"), shape = 2, size = 5) + 
  geom_segment(aes(x = means1[1], y = 1, xend = means2[1], yend = -1, color = "a")) +
  geom_segment(aes(x = means1[2], y = 1, xend = means2[2], yend = -1, color = "b"))
ggsave("1D.png")
t.test(s[which(y==1)], s[which(y==-1)])
