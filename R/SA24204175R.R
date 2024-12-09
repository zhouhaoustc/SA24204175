#' Oracle Trans-Lasso
#'
#' @param X target design matrix
#' @param y response
#' @param A0 informative auxiliary samples
#' @param n.vec sample index
#' @param lam.const tuning parameter
#' @param l1 l1=T means using l1 method
#'
#' @import Matrix glmnet  flare mvtnorm microbenchmark boot coda bootstrap DAAG Rcpp
#' @importFrom stats coef dnorm median pnorm qnorm symnum predict
#' @useDynLib SA24204175
#'
#' @export
las.kA<-function(X, y, A0, n.vec, lam.const=NULL, l1=T){
  p<-ncol(X)
  size.A0<- length(A0)
  if(size.A0 > 0){
    ind.kA<- ind.set(n.vec, c(1, A0+1))
    ind.1<-1:n.vec[1]
    if(l1){
      y.A<-y[ind.kA]
    }else{ #the l0-method
      y.A<- y[ind.1]
      Sig.hat<-t(X)%*%X/nrow(X)
      for(k in 1:size.A0){
        ind.k<- ind.set(n.vec,k+1)
        lam.k <- sqrt(mean(y[ind.1]^2)/n.vec[1]+mean(y[ind.k]^2)/n.vec[k]) * sqrt(2*log(p))
        delta.hat.k<-lassoshooting(XtX=Sig.hat,
                                   Xty=t(X[ind.k,])%*%y[ind.k]/n.vec[k+1]-t(X[1:n.vec[1],])%*%y[1:n.vec[1]]/n.vec[1],
                                   lambda=lam.k)$coef
        y.A<-c(y.A, y[ind.k]-X[ind.k,]%*%delta.hat.k)
      }

    }
    if(is.null(lam.const)){
      cv.init<-cv.glmnet(X[ind.kA,], y.A, nfolds=8, lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/length(ind.kA)))
      lam.const <- cv.init$lambda.min/sqrt(2*log(p)/length(ind.kA))
    }
    w.kA <- as.numeric(glmnet(X[ind.kA,], y.A, lambda=lam.const*sqrt(2*log(p)/length(ind.kA)))$beta)
    w.kA<-w.kA*(abs(w.kA)>=lam.const*sqrt(2*log(p)/length(ind.kA)))
    # cv.delta<-cv.glmnet(x=X[ind.1,],y=y[ind.1]-X[ind.1,]%*%w.kA, lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/length(ind.1)))
    #delta.kA<-predict(cv.delta, s='lambda.min', type='coefficients')[-1]
    delta.kA <- as.numeric(glmnet(x=X[ind.1,],y=y[ind.1]-X[ind.1,]%*%w.kA, lambda=lam.const*sqrt(2*log(p)/length(ind.1)))$beta)
    delta.kA<-delta.kA*(abs(delta.kA)>=lam.const*sqrt(2*log(p)/length(ind.1)))
    beta.kA <- w.kA + delta.kA
    lam.const=NA
  }else{
    cv.init<-cv.glmnet(X[1:n.vec[1],], y[1:n.vec[1]], nfolds=8, lambda=seq(1,0.1,length.out=20)*sqrt(2*log(p)/n.vec[1]))
    lam.const<-cv.init$lambda.min/sqrt(2*log(p)/n.vec[1])
    beta.kA <- predict(cv.init, s='lambda.min', type='coefficients')[-1]
    w.kA<-NA
  }
  list(beta.kA=as.numeric(beta.kA),w.kA=w.kA, lam.const=lam.const)

}

lassoshooting <- function(XtX, Xty, lambda) {
  p <- length(Xty)
  beta <- rep(0, p)
  for (j in 1:p) {
    r_j <- Xty[j] - sum(XtX[j, -j] * beta[-j])
    beta[j] <- soft_threshold(r_j, lambda)
  }
  return(beta)
}

soft_threshold <- function(x, lambda) {
  if (x > lambda) {
    return(x - lambda)
  } else if (x < -lambda) {
    return(x + lambda)
  } else {
    return(0)
  }
}

ind.set<- function(n.vec, k.vec){
  ind.re <- NULL
  for(k in k.vec){
    if(k==1){
      ind.re<-c(ind.re,1: n.vec[1])
    }else{
      ind.re<- c(ind.re, (sum(n.vec[1:(k-1)])+1): sum(n.vec[1:k]))
    }
  }
  ind.re
}

#' computing the MSE
#' @param beta coefficient vector
#' @param est estimator
#' @param X.test test dataset
#'
#' @export
mse.fun<- function(beta,est, X.test=NULL){
  pred.err<-NA
  est.err<- sum((beta-est)^2)

  if(!is.null(X.test)){
    pred.err<-  mean((X.test%*%(beta-est))^2)
  }
  return(list(est.err=est.err, pred.err= pred.err))
}
