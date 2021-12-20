#' @title Cross-validation
#' @description Use the cross-validation to choose the optimal bandwidth of the kernel function in the exponential generalized linear model
#' @param X the sample vector X (numeric)
#' @param Y the sample vector Y (numeric)
#' @param h the bandwidth of the kernel function (numeric)
#' @param n the size of the sample
#' @importFrom stats runif rexp
#' @importFrom graphics abline lines
#' @return the the optimal bandwidth of the kernel function 
#' @examples
#' \dontrun{
#'     n <- 200
#'     truebeta <- c(2,2)
#'     lambda <- function(x) exp(truebeta[1] + truebeta[2] * x)
#'     X <- runif(n = n, -3, 3)
#'     Y <- rexp(n = n, rate = lambda(X))
#'     h <- seq(0.1,1, by = .1)
#'     LCV(X,Y,h,n)
#' }
#' @export
LCV<-function(X,Y,h,n){
  res0<-matrix(0,ncol=n,nrow=length(h))
  res<-numeric(length(h))
  for(l in 1:length(h)){
   for(i in 1:n){
    Kernel <- dnorm(x = X[i], mean = X[-i], sd = h[l])
    beta<-nlm(f = function(beta) {
      sum(Kernel * (Y[-i] * exp(beta[1] + beta[2] * (X[-i] - X[i]))
                    - ( beta[1] + beta[2] * (X[-i] - X[i]))))
    },p = c(0, 0))$estimate[1]
    res0[l,i]<-(beta - exp(beta) * Y[i])
    }
  }
  res<-rowSums(res0)
  opt_h<-h[which.max(res)]
  plot(h, res, type = "o")
  abline(v = opt_h, col = 2)
  return(opt_h)
}

#' @title Estimate of the conditional expectation
#' @description Estimate of the conditional expectation in the exponential generalized linear model, Compare the true value with the estimate
#' @param X the sample vector X (numeric)
#' @param x evaluation grid
#' @param Y the sample vector Y (numeric)
#' @param h the bandwidth of the kernel function (numeric)
#' @param n the size of the sample
#' @importFrom stats dnorm nlm
#' @return the value of the (i) item of the Log-likelihood of beta
#' @examples
#' \dontrun{
#'     n <- 200
#'     truebeta <- c(2,2)
#'     lambda <- function(x) exp(truebeta[1] + truebeta[2] * x)
#'     X <- runif(n = n, -3, 3)
#'     Y <- rexp(n = n, rate = lambda(X))
#'     h <- seq(0.1,1, by = .1)
#'     x <- seq(-3, 3, l = 501)
#'     likelifunc(x,X,Y,h,n)
#' }
#' @export
likelifunc<-function(x,X,Y,h,n){
  opt_h<-LCV(X,Y,h,n)
  
  fitNlm.beta0 <- sapply(x, function(x) {
    K <- dnorm(x = x, mean = X, sd =opt_h)
    nlm(f = function(beta) {
      sum(K * (Y * exp(beta[1] + beta[2] * (X - x))
               - (beta[1] + beta[2] * (X - x))))
    }, p = c(0, 0))$estimate[1]
  })
  fitNlm.beta1 <- sapply(x, function(x) {
    K <- dnorm(x = x, mean = X, sd =opt_h)
    nlm(f = function(beta) {
      sum(K * (Y * exp(beta[1] + beta[2] * (X - x))
               - (beta[1] + beta[2] * (X - x))))
    }, p = c(0, 0))$estimate[2]
  })
  y<-numeric(length(fitNlm.beta1))
  for(i in 1:length(fitNlm.beta1)){
    y[i]<-exp(fitNlm.beta0+fitNlm.beta1*x[i])
  }
  plot(x,y, col = 2, lwd = 2, lty = 2)
  }
