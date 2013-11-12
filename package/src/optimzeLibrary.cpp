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
  arma::mat A = sigma2 * makeBlockDiagonalMat(Ome2, nDomains);
  A.diag() = A.diag() + sigmaSamplingError;
  return A;
}

// [[Rcpp::export]]
arma::mat matV1(double sigma1, arma::mat Ome1, arma::mat A, arma::mat Z1) {
  return sigma1 * Z1 * Ome1 * Z1.t() + A;
}

// [[Rcpp::export]]
arma::mat matV2(arma::mat W, double rho1, double sigma1, double sigma2, arma::mat Ome2, arma::mat Z1, int nDomains, arma::colvec sigmaSamplingError) {
  arma::mat A = matA(sigma2, Ome2, nDomains, sigmaSamplingError);
  arma::mat Ome1 = sigma1 * reSAR1(W, rho1);
  return Z1 * Ome1 * Z1.t() + A;
}


// [[Rcpp::export]]
arma::mat matVSolved1(double rho1, double sigma1, double rho2, double sigma2, arma::mat A, arma::mat Ome1, arma::mat Z1) {
  arma::mat solvedA = inv(A);
  return solvedA - solvedA * Z1 * inv(pow(sigma1, -1) * inv(Ome1) + Z1.t() * solvedA * Z1) * Z1.t() * solvedA;
}

// [[Rcpp::export]]
arma::mat matVSolved2(arma::mat W, double rho1, double sigma1, double sigma2, arma::mat Ome2, arma::mat Z1, int nDomains, arma::colvec sigmaSamplingError) {
  arma::mat V = matV2(W, rho1, sigma1, sigma2, Ome2, Z1, nDomains, sigmaSamplingError);
  arma::mat Vinv = inv(trimatu(chol(V)));
  return Rcpp::List::create(wrap(sigma2));
}


// [[Rcpp::export]]
arma::mat matP(arma::mat solvedV, arma::mat X) {
  arma::mat tmp1 = inv(sympd(X.t() * solvedV * X));
  return solvedV - solvedV * X * tmp1 * X.t() * solvedV;
}




