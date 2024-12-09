## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
# 参数设置
N <- 1000
num_null <- 950
alpha <- 0.1
m <- 10000

# 结果统计
fwer_bonf <- fwer_bh <- fdr_bonf <- fdr_bh <- tpr_bonf <- tpr_bh <- 0

# 模拟
for (i in 1:m) {
  # 生成p值
  p_null <- runif(num_null)
  p_alt <- rbeta(N - num_null, 0.1, 1)
  p_values <- c(p_null, p_alt)
  
  # Bonferroni校正
  bonf_sig <- p_values <= alpha / N
  fwer_bonf <- fwer_bonf + any(bonf_sig[1:num_null])
  fdr_bonf <- fdr_bonf + sum(bonf_sig[1:num_null]) / max(1, sum(bonf_sig))
  tpr_bonf <- tpr_bonf + sum(bonf_sig[(num_null + 1):N]) / (N - num_null)
  
  # B-H校正
  p_sorted <- sort(p_values)
  bh_threshold <- (1:N) / N * alpha
  bh_sig <- p_sorted <= bh_threshold
  bh_sig <- p_values <= p_sorted[max(which(bh_sig))]
  fwer_bh <- fwer_bh + any(bh_sig[1:num_null])
  fdr_bh <- fdr_bh + sum(bh_sig[1:num_null]) / max(1, sum(bh_sig))
  tpr_bh <- tpr_bh + sum(bh_sig[(num_null + 1):N]) / (N - num_null)
}

# 平均化
fwer_bonf <- fwer_bonf / m
fdr_bonf <- fdr_bonf / m
tpr_bonf <- tpr_bonf / m
fwer_bh <- fwer_bh / m
fdr_bh <- fdr_bh / m
tpr_bh <- tpr_bh / m

# 输出结果
result <- matrix(c(fwer_bonf, fdr_bonf, tpr_bonf, fwer_bh, fdr_bh, tpr_bh), 
                 nrow = 3, byrow = FALSE,
                 dimnames = list(c("FWER", "FDR", "TPR"), 
                                 c("Bonferroni", "B-H")))
print(result)


## -----------------------------------------------------------------------------
library(boot)
x <- aircondit[1]
rate <- function(x, i){
  return(1/mean(as.matrix(x[i,])))
}
boot(x, statistic = rate, R = 1000)

## -----------------------------------------------------------------------------
meant <- function(x, i){
  return(mean(as.matrix(x[i,])))
}
a <- boot(x, statistic = meant, R = 1000)
a
boot.ci(a, type = c("norm", "perc", "basic", "bca"))
hist(a$t, prob = TRUE, main = "", xlab = "")

