## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(boot)
a <- c(4, 2, 9)
A1 <- rbind(c(2, 1, 1), c(1, -1, 3))
b1 <- c(2, 3)
simplex(a = a, A1 = A1, b1 = b1, maxi = TRUE)

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

la1 <- lapply(formulas, lm, data = mtcars)

la2 <- lapply(formulas, function(x) lm(formula = x, data = mtcars))

lf1 <- vector("list", length(formulas))
for (i in seq_along(formulas)){
  lf1[[i]] <- lm(formulas[[i]], data = mtcars)
}

## -----------------------------------------------------------------------------
la1

## -----------------------------------------------------------------------------
la2

## -----------------------------------------------------------------------------
lf1

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

la <- lapply(bootstraps, lm, formula = mpg ~ disp)

lf <- vector("list", length(bootstraps))
for (i in seq_along(bootstraps)){
  lf[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}

## -----------------------------------------------------------------------------
la

## -----------------------------------------------------------------------------
lf

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
# Exercise3
sapply(la1, rsq)
sapply(la2, rsq)
sapply(lf1, rsq)

## -----------------------------------------------------------------------------
# Exercise4
sapply(la, rsq)
sapply(lf, rsq)

## -----------------------------------------------------------------------------
trials <- replicate(
  100, 
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
sapply(trials, function(x) x[["p.value"]])

## -----------------------------------------------------------------------------
sapply(trials, "[[", "p.value")

## -----------------------------------------------------------------------------
# For example
testlist <- list(iris, mtcars, cars)
lapply(testlist, function(x) vapply(x, mean, numeric(1)))

## -----------------------------------------------------------------------------
# A combination of Map() and vapply() 
lmapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE){return(simplify2array(out))}
  out
}
lmapply(testlist, mean, numeric(1))

## -----------------------------------------------------------------------------
fast_chisq_test <- function(x, y) {
  stopifnot(is.integer(x), is.integer(y), !any(is.na(x)), !any(is.na(y)))
  tab <- table(x, y)
  E <- outer(rowSums(tab), colSums(tab), "*") / sum(tab)
  chi_sq_stat <- sum((tab - E)^2 / E)
  return(chi_sq_stat)
}

set.seed(1)
x <- sample(1:3, 100, replace = TRUE)
y <- sample(1:3, 100, replace = TRUE)
fast_chisq_test(x, y)

## -----------------------------------------------------------------------------
fast_table <- function(x, y) {
  stopifnot(is.integer(x), is.integer(y))
  n_x <- max(x)
  n_y <- max(y)
  result <- matrix(0, nrow = n_x, ncol = n_y)
  for (i in seq_along(x)) {
    result[x[i], y[i]] <- result[x[i], y[i]] + 1
  }
  return(result)
}

set.seed(1)
x <- sample(1:10, 1e6, replace = TRUE)
y <- sample(1:10, 1e6, replace = TRUE)
res <- fast_table(x, y)
res

