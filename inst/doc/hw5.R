## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 7,   # 调整图形宽度
  fig.height = 5,  # 调整图形高度
  fig.align = 'center',  # 居中对齐
  out.width = '80%',     # 输出宽度
  mar = c(4, 4, 2, 1)    # 调整图形边距 
)

## -----------------------------------------------------------------------------
library(bootstrap)
x <- as.matrix(scor)
n <- nrow(x)
theta_jack <- numeric(n)
lambda <- eigen(cov(x))$values
theta_hat <- max(lambda/sum(lambda))
for (i in 1:n) {
  y <- x[-i, ]
  z <- cov(y)
  lambda <- eigen(z)$values
  theta_jack[i] <- max(lambda/sum(lambda))
}
bias_jack <- (n - 1) * (mean(theta_jack) - theta_hat)
se_jack <- sqrt((n-1)/n*sum((theta_jack-mean(theta_jack))^2))
list(est = theta_hat, bias = bias_jack, se = se_jack)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
a <- seq(10, 40, 0.1)
par(mfrow = c(2, 2))
L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main = "Linear", pch = 16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd = 2)
L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main = "Quadratic", pch = 16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd = 2)
L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main = "Exponential", pch = 16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd = 2)
c2 <- chemical^2
c3 <- chemical^3
L4 <- lm(magnetic ~ chemical + c2 + c3)
plot(chemical, magnetic, main = "Cubic", pch = 16)
yhat4 <- L4$coef[1] + L4$coef[2] * a + L4$coef[3] * a^2 +L4$coef[4] * a^3
lines(a, yhat4, lwd = 2)
models <- list(L1, L2, L3, L4)  
Rsq <- sapply(models, function(model) summary(model)$adj.r.squared)
Rsq

## -----------------------------------------------------------------------------
n <- length(magnetic)  
e1 <- e2 <- e3 <- e4 <- numeric(n)  
for (k in seq_len(n)) {  
  y <- magnetic[-k]  
  x <- chemical[-k]  
  J1 <- lm(y ~ x)  
  e1[k] <- magnetic[k] - predict(J1, newdata = data.frame(x = chemical[k]))  
  J2 <- lm(y ~ x + I(x^2))  
  e2[k] <- magnetic[k] - predict(J2, newdata = data.frame(x = chemical[k]))  
  J3 <- lm(log(y) ~ x)  
  e3[k] <- magnetic[k] - exp(predict(J3, newdata = data.frame(x = chemical[k])))  
  J4 <- lm(y ~ x + I(x^2) + I(x^3))  
  e4[k] <- magnetic[k] - predict(J4, newdata = data.frame(x = chemical[k]))  
}  
  
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
detach(ironslag)
detach(package:DAAG)

## -----------------------------------------------------------------------------
cvm.test <- function(x, y, R = 199) {  
  z <- c(x, y)  
  N <- length(z)  
  n <- length(x)  
  m <- length(y)  
    
  # 计算原始的CvM统计量  
  Fn <- sapply(z, function(zi) mean(x >= zi))  
  Gm <- sapply(z, function(zi) mean(y >= zi))  
  cvm0 <- ((n * m) / N) * sum((Fn - Gm)^2)  
    
  # 使用replicate进行重复抽样并计算CvM统计量  
  cvm <- replicate(R, {  
    k <- sample(N)  
    Z <- z[k]  
    X <- Z[1:n]  
    Y <- Z[(n + 1):N]  
      
    Fn_perm <- sapply(Z, function(zi) mean(X >= zi))  
    Gm_perm <- sapply(Z, function(zi) mean(Y >= zi))  
      
    ((n * m) / N) * sum((Fn_perm - Gm_perm)^2)  
  })  
    
  # 计算p值  
  p.value <- mean(c(cvm, cvm0) >= cvm0)  
    
  # 返回包含统计量和p值的列表  
  return(list(statistic = cvm0, p.value = p.value))  
}
attach(chickwts)
x1 <- as.vector(weight[feed == "soybean"])
x2 <- as.vector(weight[feed == "sunflower"])
x3 <- as.vector(weight[feed == "linseed"])
detach(chickwts)
cvm.test(x1, x3)
cvm.test(x2, x3)

## -----------------------------------------------------------------------------
# 生成示例数据  
set.seed(123)  
x <- rnorm(50)  
y <- rnorm(50) + 0.5 * x  
  
# 计算原始Spearman相关系数  
rho_original <- cor(x, y, method = "spearman")  
  
# 设置置换次数  
R <- 999  
  
# 初始化一个向量来存储置换相关系数  
rho_perm <- numeric(R)  
  
# 执行置换检验  
for (i in 1:R) {  
  y_perm <- sample(y) # 置换y变量  
  rho_perm[i] <- cor(x, y_perm, method = "spearman") # 计算置换相关系数  
}  
  
# 计算置换p值  
p_perm <- mean(abs(rho_perm) >= abs(rho_original))  
  
# 使用cor.test计算p值  
p_cor_test <- cor.test(x, y, method = "spearman")$p.value  
  
# 打印结果  
cat("置换检验的p值:", p_perm, "\n")  
cat("cor.test的p值:", p_cor_test, "\n")

