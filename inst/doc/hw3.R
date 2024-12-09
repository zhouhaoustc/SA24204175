## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
pBETA <- function(x, a, b, m = 10000){
  u <- rgamma(m, a, 1)
  v <- rgamma(m, b, 1)
  y <- u/(u + v)
  return(mean(y<=x))
}
x <- seq(0.1, 0.9, 0.1)
k <- length(x)
p <- numeric(k)
for (i in 1:k) {
  p[i] <- pBETA(x[i], 3, 3)
}
phi <- pbeta(x, 3, 3)
round(rbind(x, phi, p), 3)

## -----------------------------------------------------------------------------
Ray1 <- function(n, sigma){
  u <- runif(n)
  return(sigma * sqrt(-2 * log(u)))
}
Ray2 <- function(n, sigma){
  u <- runif(n/2)
  x1 <- sigma * sqrt(-2 * log(u))
  x2 <- sigma * sqrt(-2 * log(1-u))
  return(c(x1,x2))
}
m <- 10000
sigma <- 1
r1 <- replicate(1000, mean(Ray1(2, sigma)))
r2 <- replicate(1000, mean(Ray2(2, sigma)))
var(r1)
var(r2)

## -----------------------------------------------------------------------------
100 * (var(r1)-var(r2))/var(r1)

## -----------------------------------------------------------------------------
x <- seq(1, 10, 0.01)
y <- x^2 * exp(-x^2/2)/sqrt(2*pi)
plot(x, y, type = "l", ylim = c(0, 1))
lines(x, 2*dnorm(x,1), lty = 2)
lines(x, dgamma(x-1, 3/2, 2), lty = 3)
legend("topright", legend = c("g(x)","f1","f2"), lty = 1:3)

## -----------------------------------------------------------------------------
plot(x, y/(2 * dnorm(x, 1)), type = "l", lty =2, ylim = c(0,0.8), ylab = "")
lines(x, y/(dgamma(x-1, 3/2, 2)), lty = 3)
legend("topright", legend = c("f1", "f2"), lty = 2:3)

## -----------------------------------------------------------------------------
fast_sort <- function(x){
  if(length(x)<=1){
    return(x)
  }
  index <- sample(1:length(x), 1)
  pivot <- x[index]
  left <- x[x < pivot]
  right <- x[x > pivot]
  equal <- x[x == pivot]
  c(fast_sort(left), equal, fast_sort(right))
}

## -----------------------------------------------------------------------------
mc1 <- function(n, m=100){
  times <- numeric(m)
  for (i in 1:m) {
    x <- sample(n)
    start_time <- proc.time()
    sorted_vec <- fast_sort(x)
    end_time <- proc.time()
    times[i] <- end_time["elapsed"] - start_time["elapsed"]
  }
  mean(times)
}

## -----------------------------------------------------------------------------
n <- c(10^4, 2*10^4, 4*10^4, 6*10^4, 8*10^4)
a_n <- numeric(length(n))
t_n <- numeric(length(n))
for (i in seq_along(n)) {
  a_n[i] <- mc1(n[i])
  t_n[i] <- n[i] * log(n[i])
}

## -----------------------------------------------------------------------------
lm1 <- lm(a_n ~ t_n)
plot(t_n, a_n, pch=19, col="blue",
  xlab="t_n = n*log(n)",
  ylab=" 平均计算时间 (a_n)",
  main=" 蒙特卡洛实验：快速排序性能分析 ")
abline(lm1, col="red")

