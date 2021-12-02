#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

colvec GdensV(const colvec x){
  return exp(-0.5 * x % x) / sqrt(2 * M_PI);
}

double Gdens(const double x){
  return exp(-0.5 * x * x) / sqrt(2 * M_PI);
}
Col<double> Gdensbis(Col<double> x){
  return exp(-0.5 * x % x) / sqrt(2 * M_PI);
}

// [[Rcpp::export]]
NumericVector obj2Cpp(const NumericVector& uR, const double xij, const NumericVector& xj, const NumericVector& weightsR, const double band){
  Col<double> u = uR;
  Col<double> weights = weightsR;
  colvec out = u;
  for (uword i=0; i<u.n_rows; i++){
    out(i) = log(sum(GdensV((xj - u[i]) / band) % weights) / (band * sum(weights))) * Gdens((xij - u[i]) / band) / band;
  }
  return wrap(out);
}


// [[Rcpp::export]]
NumericVector obj3Cpp(const NumericVector& uR,  const NumericVector& xj, const NumericVector& weightsR, const double band){
  Col<double> u = uR;
  Col<double> weights = weightsR;
  Col<double> obs = xj;
  Col<double> out = zeros(obs.n_rows);
  double tmp=0;
  for (uword r=0; r<u.n_rows; r++){
    tmp = log(sum(GdensV((obs - u[r]) / band) % weights) / (band * sum(weights)));
    if (tmp != log(0)){
      out +=  tmp* GdensV((u[r] - obs)/band)/band;
    }
  }
  return wrap(out);
}
