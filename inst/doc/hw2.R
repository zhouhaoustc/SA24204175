## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 7,   # 调整图形宽度
  fig.height = 5,  # 调整图形高度
  fig.align = 'center',  # 居中对齐
  out.width = '80%',     # 输出宽度
  mar = c(4, 4, 2, 1)    # 调整图形边距 
)

## -----------------------------------------------------------------------------
n <- 3000
sigma <- seq(0.5,2,length = 4)

## -----------------------------------------------------------------------------
par(mfrow = c(2, 2))
for (i in 1:4){
  u <- runif(n)
  x <- sigma[i] * sqrt(-2*log(u))
  hist(x, breaks = "Scott", prob = TRUE, main = sigma[i])
  abline(v = sigma[i], lwd =1.5)
}

## -----------------------------------------------------------------------------
n <- 1000
p <- 0.75
mu <- sample(c(0, 3), size = 1000, replace =TRUE, prob = c(p, 1-p))
x <- rnorm(n, mu, 1)

## -----------------------------------------------------------------------------
hist(x, breaks = "Scott", prob = TRUE)
y <- sort(x)
fy <- p * dnorm(y) + (1-p) * dnorm(y, 3, 1)
lines(y, fy)

## -----------------------------------------------------------------------------
par(mfrow = c(3, 3))
p <- seq(0.1, 0.9, length = 9)
for (i in 1:9) {
  mu <- sample(c(0, 3),size = 1000, replace = TRUE, prob = c(p[i], 1-p[i]))
  x <- rnorm(n, mu, 1)
  hist(x, breaks = "Scott", prob = TRUE, main = paste("p=", p[i]))
  y <-sort(x)
  fy <- p[i] * dnorm(y) + (1-p[i]) * dnorm(y, 3, 1)
  lines(y, fy)
}

## -----------------------------------------------------------------------------
alpha <- 4
beta <- 2
lambda <- seq(2, 4, length = 3)
t <- 10
Pp <- function(lambda, t){
  tn <- rexp(1000, lambda)
  sn <- cumsum(tn)
  n <- min(which(sn > t))
  return(n-1)
}
x <- numeric(1000)
for (i in 1:3) {
  for (j in 1:1000) {
    N <- Pp(lambda[i], t)
    x[j] <- sum(rgamma(N, shape = alpha, rate = beta))
  }
  print(paste("Empirical value of the mean for lambda = ", lambda[i], "is", mean(x)))
  print(paste("Theoretical value of the mean for lambda = ", lambda[i], "is", lambda[i] * t * alpha/beta))
  print(paste("Empirical value of the variance for lambda = ", lambda[i], "is", var(x)))
  print(paste("Theoretical value of the varience for lambda = ", lambda[i], "is", lambda[i] * t * (alpha/beta^2 + (alpha/beta)^2)))
}

