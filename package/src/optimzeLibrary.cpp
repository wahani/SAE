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
    XX.submat(arma::span(r * X.n_cols, (r + 1) * X.n_cols - 1), arma::span(r * X.n_cols, (r + 1) * X.n_cols - 1)) =
    X;
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
Rcpp::List matVinv(arma::mat W, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError) {
  arma::mat Ome2 = reAR1(Z1.n_rows / W.n_rows, rho2);
  arma::mat Ome1 = reSAR1(W, rho1);
  arma::mat Ad = sigma2 * Ome2;
  arma::mat A = makeBlockDiagonalMat(Ad, W.n_rows);
  A.diag() += sigmaSamplingError;
  arma::mat Ainv = A;
  
  for(int r = 0; r < W.n_rows; ++r) {
    Ainv.submat(arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1), arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1)) = 
    arma::inv(Ainv.submat(arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1), arma::span(r * Ad.n_cols, (r + 1) * Ad.n_cols - 1)));    
  }
  
  arma::mat V = sigma1 * Z1 * Ome1 * Z1.t() + A;
  arma::mat AinvZ1 = Ainv * Z1;
  arma::mat Ome1inv = arma::inv(sigma1 * Ome1);
  arma::mat Vinv = Ainv - AinvZ1 * arma::inv(Ome1inv + Z1.t() * AinvZ1) * AinvZ1.t();
  return Rcpp::List::create(Rcpp::Named("V", V),
    Rcpp::Named("Vinv", Vinv));
}

// [[Rcpp::export]]
arma::mat matV(arma::mat W, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError) {
  arma::mat Ome2 = reAR1(Z1.n_rows / W.n_rows, rho2);
  arma::mat A = matA(sigma2, Ome2, W.n_cols, sigmaSamplingError);
  arma::mat Ome1 = reSAR1(W, rho1);
  return sigma1 * Z1 * Ome1 * Z1.t() + A;
}

// [[Rcpp::export]]
Rcpp::List matVinv1(arma::mat W, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError) {
  arma::mat V = matV(W, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError);
  arma::mat Vtmp = arma::inv(trimatu(chol(V)));
  arma::mat Vinv = Vtmp * Vtmp.t();
  return Rcpp::List::create(Rcpp::Named("V", V),
  Rcpp::Named("Vinv", Vinv));
}

// [[Rcpp::export]]
arma::mat matP(arma::mat solvedV, arma::mat X) {
  arma::mat tmp1 = arma::inv(sympd(X.t() * solvedV * X));
  return solvedV - solvedV * X * tmp1 * X.t() * solvedV;
}

// [[Rcpp::export]]
arma::mat matVderS1(arma::mat Ome1, arma::mat Z1) {
  return Z1 * Ome1 * Z1.t();
}

// [[Rcpp::export]]
arma::mat matVderS2(arma::mat Ome2, int nDomains) {
  return makeBlockDiagonalMat(Ome2, nDomains);
}

// [[Rcpp::export]]
arma::mat matVderR1(double rho1, double sigma1, arma::mat Z1, arma::mat Ome1, arma::mat W) {
  return -sigma1 * Z1 * Ome1 * (-W-W.t() + 2 * rho1 * W.t() * W) * Ome1 * Z1.t();
}

// [[Rcpp::export]]
arma::mat matVderR2(double rho2, double sigma2, arma::mat Ome2, int nDomains) {
  arma::mat ome2derR2(Ome2.n_cols, Ome2.n_cols);
  ome2derR2.fill(0.0);
  for(int i = 1; i < Ome2.n_cols; ++i) {
     ome2derR2.diag(i) += i * pow(rho2, i-1);
  }
  ome2derR2 += ome2derR2.t();
  ome2derR2 = 1/(1-pow(rho2, 2)) * (ome2derR2 + 2 * rho2 * Ome2);
  return sigma2 * makeBlockDiagonalMat(ome2derR2, nDomains);
}


// [[Rcpp::export]]
arma::colvec blue(arma::colvec y, arma::mat X, arma::mat Vinv) {
  arma::mat t1 = X.t() * Vinv;
  return inv(t1 * X) * t1 * y;
}

// [[Rcpp::export]]
Rcpp::List llr(arma::colvec y, arma::mat X, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError, arma::mat W) {
  Rcpp::List Vlist = matVinv(W, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError);
  arma::mat V = Vlist(0);
  arma::mat Vinv = Vlist(1);
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
  
  arma::colvec gradient(4);
  arma::mat P = matP(Vinv, X);
  
  arma::mat VderR1 = matVderR1(rho1, sigma1, Z1, reSAR1(W, rho1), W);
  arma::mat VderS1 = matVderS1(reSAR1(W, rho1), Z1);
  arma::mat VderR2 = matVderR2(rho2, sigma2, reAR1(Z1.n_rows / W.n_rows, rho2), W.n_cols);
  arma::mat VderS2 = matVderS2(reAR1(Z1.n_rows / W.n_rows, rho2), W.n_cols);
  
  arma::mat tmpR1 = P*VderR1;
  VderR1 = -0.5 * trace(tmpR1) + 0.5 * y.t() * tmpR1 * P * y;
  gradient(0) = VderR1(0, 0);
  
  arma::mat tmpS1 = P*VderS1;
  VderS1 = -0.5 * trace(tmpS1) + 0.5 * y.t() * tmpS1 * P * y;
  gradient(1) = VderS1(0, 0);
  
  arma::mat tmpR2 = P*VderR2;
  VderR2 = -0.5 * trace(tmpR2) + 0.5 * y.t() * tmpR2 * P * y;
  gradient(2) = VderR2(0, 0);
  
  arma::mat tmpS2 = P*VderS2;
  VderS2 = -0.5 * trace(tmpS2) + 0.5 * y.t() * tmpS2 * P * y;
  gradient(3) = VderS2(0, 0);
  
  return Rcpp::List::create(Rcpp::Named("objective", -llr),
  Rcpp::Named("gradient", gradient));
}


arma::colvec psiOne(arma::colvec u, double k = 1.345){
  //arma::colvec weights = rep(1, length(u));
  //sm<-median(abs(u/sqrt(var.weights)))/0.6745
  double sm = median(abs(u)) / 0.6745;
  Rcpp::NumericVector tmp = as<NumericVector>(wrap(k/abs(u/sm)));
  arma::colvec w = Rcpp::pmin(1.0, tmp);
  return w%u;
}

// [[Rcpp::export]]
double optimizerRho(arma::colvec rho, arma::colvec y, arma::mat X, double sigma1, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError, arma::mat W, arma::colvec beta, double K) {
  
  // Variance-Covarianze matrix
  Rcpp::List Vlist = matVinv(W, rho(0), sigma1, rho(1), sigma2, Z1, sigmaSamplingError);
  arma::mat V = Vlist(0);
  arma::mat Vinv = Vlist(1);
  
  // sqrt of U + inverse
  arma::mat sqrtU(V.n_rows, V.n_rows);
  sqrtU.fill(0.0);
  sqrtU.diag() = sqrt(V.diag());
  
  arma::mat sqrtUinv(V.n_rows, V.n_rows);
  sqrtUinv.fill(0.0);
  sqrtUinv.diag() = 1/sqrtU.diag();
  
  // residuals and huber
  arma::colvec resid = sqrtUinv * (y - X * beta);
  arma::colvec phiR = psiOne(resid);
  
  arma::mat Ome1 = reSAR1(W, rho(0));
  arma::mat Ome2 = reAR1(Z1.n_rows / W.n_rows, rho(1));
  arma::mat derVSarCorr = matVderR1(rho(0), sigma1, Z1, Ome1, W);
  arma::mat derVArCorr = matVderR2(rho(1), sigma2, Ome2, W.n_rows);
  
  arma::mat tmp1 = phiR.t() * sqrtU * Vinv;
        
  double tmp2 = trace(K * Vinv * derVSarCorr);
  double tmp3 = trace(K * Vinv * derVArCorr);
  
  arma::mat tmp4 = tmp1 * derVSarCorr * tmp1.t();
  arma::mat tmp5 = tmp1 * derVArCorr * tmp1.t();
  
  double optRho1 = tmp4(0, 0) - tmp2;
  double optRho2 = tmp5(0, 0) - tmp3;
      
  return pow(optRho1, 2) + pow(optRho2, 2);
//    return Rcpp::List::create(Rcpp::Named("V", V),
//    Rcpp::Named("Vinv", Vinv),
//    Rcpp::Named("sqrtU", sqrtU),
//    Rcpp::Named("sqrtUinv", sqrtUinv),
//    Rcpp::Named("resid", resid),
//    Rcpp::Named("phiR", phiR),
//    Rcpp::Named("tmp1", tmp1),
//    Rcpp::Named("tmp2", tmp2),
//    Rcpp::Named("tmp3", tmp3),
//    Rcpp::Named("tmp4", tmp4),
//    Rcpp::Named("tmp5", tmp5),
//    Rcpp::Named("optRho1", optRho1),
//    Rcpp::Named("optRho2", sqrtU));
}

// [[Rcpp::export]]
arma::colvec optimizerSigma(arma::colvec sigma, arma::colvec rho, arma::colvec y, arma::mat X, arma::mat Z1, arma::colvec sigmaSamplingError, arma::mat W, arma::colvec beta, double K, arma::mat Z) {
            
  // Variance-Covarianze matrix
  Rcpp::List Vlist = matVinv(W, rho(0), sigma(0), rho(1), sigma(1), Z1, sigmaSamplingError);
  arma::mat V = Vlist(0);
  arma::mat Vinv = Vlist(1);
      
  // sqrt of U + inverse
  arma::mat sqrtU(V.n_rows, V.n_rows);
  sqrtU.fill(0.0);
  sqrtU.diag() = sqrt(V.diag());
  
  arma::mat sqrtUinv(V.n_rows, V.n_rows);
  sqrtUinv.fill(0.0);
  sqrtUinv.diag() = 1/sqrtU.diag();
  
  
  // residuals and huber
  arma::colvec resid = sqrtUinv * (y - X * beta);
  arma::colvec phiR = psiOne(resid);
  
  // Derivatives
  arma::mat Ome1 = reSAR1(W, rho(0));
  arma::mat Ome2 = reAR1(Z1.n_rows / W.n_rows, rho(1));
  arma::mat derVSigma1 = matVderS1(Ome1, Z1);
  arma::mat derVSigma2 = matVderS2(Ome2, W.n_rows);
  
  // ZVuZt + inverse
  arma::mat ZVuZt = V;
  ZVuZt.diag() -= sigmaSamplingError;
  arma::mat ZVuZtinv = arma::inv(trimatu(chol(ZVuZt)));
  ZVuZtinv = ZVuZtinv * ZVuZtinv.t();
  
  // OmegaBar
  arma::mat Ome1Bar(Z.n_cols, Z.n_cols);
  arma::mat Ome2Bar(Z.n_cols, Z.n_cols);
  Ome1Bar.fill(0.0);
  Ome2Bar.fill(0.0);
  
  Ome1Bar(arma::span(0, W.n_rows-1), arma::span(0, W.n_rows-1)) = Ome1;
  Ome2Bar(arma::span(W.n_rows, Ome2Bar.n_rows-1), arma::span(W.n_rows, Ome2Bar.n_rows-1)) = makeBlockDiagonalMat(Ome2, W.n_rows);
  
  // Compute a(theta) and A(theta)
  arma::mat tmp1 = phiR.t() * sqrtU * Vinv;
  
  arma::colvec a(2);
  
  a(0) = as_scalar(tmp1 * derVSigma1 * tmp1.t());
  a(1) = as_scalar(tmp1 * derVSigma2 * tmp1.t());
      
  arma::mat tmp2 = K * Vinv * derVSigma1 * ZVuZtinv;
  arma::mat tmp3 = K * Vinv * derVSigma2 * ZVuZtinv;
      
  arma::mat tmp4 = Z * Ome1Bar * Z.t();
  arma::mat tmp5 = Z * Ome2Bar * Z.t();
  
  arma::mat A(2,2);
  A(0,0) = trace(tmp2 * tmp4);
  A(0,1) = trace(tmp2 * tmp5);
  A(1,0) = trace(tmp3 * tmp4);
  A(1,1) = trace(tmp3 * tmp5);
  
  return inv(A) * a;
  
//    return Rcpp::List::create(Rcpp::Named("V", V),
//    Rcpp::Named("Vinv", Vinv),
//    Rcpp::Named("sqrtU", sqrtU),
//    Rcpp::Named("sqrtUinv", sqrtUinv),
//    Rcpp::Named("resid", resid),
//    Rcpp::Named("phiR", phiR),
//    Rcpp::Named("ZVuZt", ZVuZt),
//    Rcpp::Named("ZVuZtinv", ZVuZtinv),
//    Rcpp::Named("tmp1", tmp1),
//    Rcpp::Named("tmp2", tmp2),
//    Rcpp::Named("tmp3", tmp3),
//    Rcpp::Named("tmp4", tmp4),
//    Rcpp::Named("tmp5", tmp5),
//    Rcpp::Named("a", a),
//    Rcpp::Named("A", A));
}


// [[Rcpp::export]]
arma::colvec optimizeRE(arma::colvec sigma, arma::colvec rho, arma::colvec y, arma::mat X, arma::mat Z1, arma::colvec sigmaSamplingError, arma::mat W, arma::colvec beta, double K, arma::mat Z, double tol, int maxit) {
  
  int n = Z.n_rows;
  
  // sqrt of R inverse
  arma::mat sqrtRinv(n, n); sqrtRinv.fill(0.0); sqrtRinv.diag() = 1/sqrt(sigmaSamplingError);
  
  // RE - Spatio-Temporal
  arma::mat Ome1 = reSAR1(W, rho(0));
  arma::mat Ome2 = reAR1(Z1.n_rows / W.n_rows, rho(1));
  
  arma::mat G(Z.n_cols, Z.n_cols); G.fill(0.0);
  G(arma::span(0, W.n_rows-1), arma::span(0, W.n_rows-1)) = sigma(0) * Ome1;
  G(arma::span(W.n_rows, G.n_rows-1), arma::span(W.n_rows, G.n_rows-1)) = sigma(1) * makeBlockDiagonalMat(Ome2, W.n_rows);
    
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, G);

  arma::mat sqrtGinv = eigvec * arma::diagmat(sqrt(1/eigval)) * eigvec.t();
  
  // Variance-Covarianze matrix
  Rcpp::List Vlist = matVinv(W, rho(0), sigma(0), rho(1), sigma(1), Z1, sigmaSamplingError);
  arma::mat V = Vlist(0);
  arma::mat Vinv = Vlist(1);
  
  arma::colvec resid = y - X * beta;
  arma::colvec vv = G * Z.t() * Vinv * resid;
  
  // Algorithm
  int nDomains = W.n_rows;
  double diff = 1.0;
  
  int i = 0;
  arma::mat tmp = Z.t() * sqrtRinv;
  
  while (diff > tol)
  {
    i++;
    arma::colvec v_robust = vv;
    arma::colvec res1 = sqrtRinv * (resid - Z * v_robust);
    arma::colvec res2 = sqrtGinv * v_robust;
    
    arma::mat w2 = arma::diagmat(psiOne(res1)/res1);
    arma::mat w3 = arma::diagmat(psiOne(res2)/res2);
    
    arma::mat tmp1 = tmp * w2 * sqrtRinv;
    arma::mat Atmp1 =  tmp1 * Z;
    arma::mat Atmp2 = sqrtGinv * w3 * sqrtGinv;
    arma::mat A = Atmp1 + Atmp2;
    arma::mat B = tmp1 * resid;
    
    vv = inv(A) * B;
    
    diff = sum(pow(vv-v_robust, 2));
    if (i > maxit) break;
  }
  
  return Z * vv;
  
//    return Rcpp::List::create(Rcpp::Named("G", G),
//    Rcpp::Named("sqrtGinv", sqrtGinv),
//    Rcpp::Named("sqrtRinv", sqrtRinv),
//    Rcpp::Named("resid", resid),
//    Rcpp::Named("vv", vv),
//    Rcpp::Named("v_robust", v_robust),
//    Rcpp::Named("res1", res2),
//    Rcpp::Named("A", A),
//    Rcpp::Named("B", B));
}