// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat reSAR1(arma::mat W, double rho) {
   int n = W.n_rows;
   arma::mat I = arma::eye<arma::mat>(n, n);
   arma::mat S1 = I - rho * W;
   return inv(sympd(S1.t() * S1));
}

// [[Rcpp::export]]
arma::mat reAR1(int nTime, double rho) {
   arma::mat Ome2(nTime, nTime);
   Ome2.fill(0.0);
   for(int i = 0; i < nTime; ++i) {
     Ome2.diag(i) += pow(rho, i);
   }
   Ome2 += Ome2.t();
   Ome2.diag() *= 0.0;
   Ome2.diag() += 1;
   return 1/(1-pow(rho, 2)) * Ome2;
}

// [[Rcpp::export]]
arma::mat makeBlockDiagonalMat(arma::mat X, int n) {
  arma::mat XX(X.n_cols*n, X.n_cols*n);
  XX.fill(0.0);
  for(int r = 0; r < n; ++r) {
    for(int i = 0; i < X.n_rows; ++i) { 
      for(int j = 0; j < X.n_cols; ++j) { 
        XX(i + X.n_rows * r, j + X.n_cols * r) = X(i, j);
        }
      }    
  }
  return XX;
}

// [[Rcpp::export]]
arma::mat matA(double sigma2, arma::mat Ome2, int nDomains, arma::colvec sigmaSamplingError) {
  arma::mat Ad = sigma2 * Ome2;
  arma::mat A = makeBlockDiagonalMat(Ad, nDomains);
  A.diag() += sigmaSamplingError;
  return A;
}

// [[Rcpp::export]]
arma::mat matV(arma::mat W, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError) {
  arma::mat Ome2 = reAR1(Z1.n_rows / W.n_rows, rho2);
  arma::mat A = matA(sigma2, Ome2, W.n_cols, sigmaSamplingError);
  arma::mat Ome1 = reSAR1(W, rho1);
  return sigma1 * Z1 * Ome1 * Z1.t() + A;
}

// [[Rcpp::export]]
Rcpp::List matVinv(arma::mat W, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError) {
  arma::mat V = matV(W, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError);
  arma::mat Vtmp = inv(trimatu(chol(V)));
  arma::mat Vinv = Vtmp * Vtmp.t();
  return Rcpp::List::create(Rcpp::Named("V", V),
  Rcpp::Named("Vinv", Vinv));
}

// [[Rcpp::export]]
arma::mat matP(arma::mat solvedV, arma::mat X) {
  arma::mat tmp1 = inv(sympd(X.t() * solvedV * X));
  return solvedV - solvedV * X * tmp1 * X.t() * solvedV;
}

// [[Rcpp::export]]
arma::colvec blue(arma::colvec y, arma::mat X, arma::mat Vinv) {
  arma::mat t1 = X.t() * Vinv;
  return inv(t1 * X) * t1 * y;
}

// [[Rcpp::export]]
double llr(arma::colvec y, arma::mat X, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError, arma::mat W) {
  arma::mat V = matV(W, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError);
  arma::mat Vtmp = inv(trimatu(chol(V)));
  arma::mat Vinv = Vtmp * Vtmp.t();
  arma::colvec beta = blue(y, X, Vinv);
  
  arma::mat resid = y - X * beta;
  double val;
  double sign;
  log_det(val, sign, V);
  arma::mat tmp = resid.t() * Vinv * resid;
  
  double val1;
  double sign1;
  log_det(val1, sign1, X.t() * Vinv * X);
  
  double llp = -0.5 * (val + tmp(0,0)); 
  double llr = llp - 0.5 * val1;
  return llr;
}

