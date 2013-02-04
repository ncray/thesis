source("MKL_code.R")
cache_size <- 50
lambda <- 100
##n <- 50
n <- 200
##C <- .5
##C <- .1
C <- .01
Nwts <- 100
Npwr <- 50
##r.v <- c(0.01, 0.1, 1, 10, 1000)
##r.v <- c(0.01, 0.1, 1, 10)
##r.v <- c(0.1, 1, 10)
r.v <- c(1, 10)

system.time(res <- mdply(expand.grid(r1 = 4.5, self = seq(.25, .35, .03), C = .001), MKLwts))
names(res) <- c("r1", "self", "C", "perm", paste("rbf", r.v, sep = ": "), paste("sk", 1:2, sep = ": "))
res.m <- melt(res, id.vars = c(1:4))
p1 <- qplot(variable, value, data = subset(res.m, perm == 1), geom = "boxplot") +
  facet_grid(self~r1) +
  geom_point(data = subset(res.m, perm == 0), color = "red", size = 5) +
  xlab("Kernels") +
  ylab("Kernel Weights") +
  opts(title = "Boxplot of Null Distribution with Observed in Red Faceted by Self Transition Probability and Outer Radius")
p1
##ggsave("mkl_weight_boxplot_christmas_star_DNA.png", width = 14, height = 12)

computeFS(u, l, 1)
computeFS(u, sample(l), 1)
system.time(rejectFS(u, l, 1))
system.time(rejectFSMKL(u, l, r.v))



#c(.05, .3, .5, .8, 5)
system.time(res <- mdply(expand.grid(r1 = seq(4, 6, 1), self = seq(.25, .7, .15), C = .5), power))
Npwr <- 100
system.time(res <- mdply(expand.grid(r1 = 4.4, self = 0.28, C = .001), power))
system.time(res <- mdply(expand.grid(r1 = 4.25, self = 0.27, C = .01), power))
###system.time(res <- mdply(expand.grid(r1 = 4.35, self = 0.28, C = .01), power))  ggsave("powerMKL_different_samples.png", width = 18, height = 12)
names(res) <- c("r1", "self", "C", "FSMKL", paste("rbf", r.v, sep = ": "), paste("sk", 1:2, sep = ": "))
res2 <- ddply(res, .(r1, self, C), function(df){
  exclude.inds <- -(1:3)
  means <- apply(df[, exclude.inds], 2, mean)
  lims <- c(.025, .975)
  bootM <- function(x) quantile(as.vector(boot(x, function(x, i) mean(x[i]), 1000)$t), lims)
  bounds <- do.call(rbind, llply(names(df)[exclude.inds], function(x) bootM(df[, x])))
  data.frame("r1" = df$r1[1], "self" = df$self[1], "C" = df$C[1], 
             "value" = means, "lower" = bounds[, 1], "upper" = bounds[, 2], group = names(df)[exclude.inds])
})

##for(i in 1:10){power(4, .25, .5)}
p2 <- ggplot(res2, aes(x = group, y = value)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .05))
p2
ggsave("powerMKL_same_samples.png", width = 18, height = 12)

p2 <- ggplot(res2, aes(x = r1, y = value, color = group, linetype = group)) +
  geom_line() + 
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .05)) +
  #facet_grid(.~self) +
  facet_grid(self~.) +
  opts(title = "Power (Christmas Star Example)") +
  xlab("Radius of Outer Star (Inner is 4)") +
  ylab("Power")
p2

p3 <- ggplot(res2, aes(x = group, y = value, fill = group)) +
  geom_bar() +
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .5)) +  
  facet_grid(self~r1) +
  opts(title = "Power (Christmas Star + DNA Example, C = .1, n = 50)") +
  xlab("Kernel") +
  ylab("Power")
p3

ggsave("power_christmas_star_DNA_face4.png", width = 18, height = 12)
