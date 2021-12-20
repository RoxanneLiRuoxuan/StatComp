#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param x the origin value of the random variable X
//' @param y the origin value of the random variable Y
//' @param N the number of samples
//' @importFrom stats rbeta rbinom
//' @useDynLib StatComp21020
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' rnC <- Mygibbs(1,0.1,1000)
//' par(mfrow=c(2,1));
//' plot(rnC[1,],type='l')
//' plot(rnC[2,],type='l')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix Mygibbs(double x,double y,int N) {
  NumericMatrix GX(2,N);
  int a=2,b=2,n=10;
  GX(0,0)=x;
  GX(1,0)=y;  
  for(int i=1;i<N;i++){
    double x0=GX(0,i-1);
    GX(1,i)=rbeta(1,x0+a,n-x0+b)[0];
    double y0=GX(1,i);
    GX(0,i)=rbinom(1,n,y0)[0];
  }
  return GX;
}
