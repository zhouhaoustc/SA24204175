## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(Rcpp)

# R版本的Gibbs采样器
gibbsSamplerR <- function(a, b, n, N, x_init) {
  chain <- matrix(0, nrow = N, ncol = 2)
  x <- x_init
  for (i in 1:N) {
    y <- rbeta(1, x + a, n - x + b)
    x <- rbinom(1, n, y)
    chain[i, ] <- c(x, y)
  }
  return(chain)
}

# Rcpp实现的Gibbs采样器
cppFunction('
NumericMatrix gibbsSampler(int a, int b, int n, int N, int x_init) {
  NumericMatrix chain(N, 2);
  int x = x_init;
  double y;
  for (int i = 0; i < N; i++) {
    y = R::rbeta(x + a, n - x + b);
    x = R::rbinom(n, y);
    chain(i, 0) = x;
    chain(i, 1) = y;
  }
  return chain;
}
')

# 设置参数
a <- 3
b <- 2
n <- 10
N <- 1000
x_init <- 5

# 使用R和Rcpp实现采样
set.seed(1)  
samples_R <- gibbsSamplerR(a, b, n, N, x_init)
samples_Rcpp <- gibbsSampler(a, b, n, N, x_init)

# 提取x和y的采样结果
x_R <- samples_R[, 1]
y_R <- samples_R[, 2]
x_Rcpp <- samples_Rcpp[, 1]
y_Rcpp <- samples_Rcpp[, 2]

# 使用qqplot比较生成的随机数
par(mfrow = c(1, 2))  

# 比较x的采样
qqplot(x_R, x_Rcpp, main = "QQ Plot of x", xlab = "R Samples (x)", ylab = "Rcpp Samples (x)")
abline(0, 1, col = "red")

# 比较y的采样
qqplot(y_R, y_Rcpp, main = "QQ Plot of y", xlab = "R Samples (y)", ylab = "Rcpp Samples (y)")
abline(0, 1, col = "red")

## -----------------------------------------------------------------------------
library(microbenchmark)
set.seed(1)
benchmark_results <- microbenchmark(
  R = gibbsSamplerR(a, b, n, N, x_init),
  Rcpp = gibbsSampler(a, b, n, N, x_init),
  times = 10  # 重复测试 10 次
)
print(benchmark_results)

## -----------------------------------------------------------------------------
boxplot(benchmark_results, main = "Computation Time Comparison", ylab = "Time (ms)")

