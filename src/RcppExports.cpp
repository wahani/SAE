// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cppArmFunc
arma::mat cppArmFunc(arma::mat X);
RcppExport SEXP saedevel_cppArmFunc(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        arma::mat __result = cppArmFunc(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppArmFuncOpt
arma::mat cppArmFuncOpt(arma::mat X);
RcppExport SEXP saedevel_cppArmFuncOpt(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        arma::mat __result = cppArmFuncOpt(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppArmFuncPseudo
arma::mat cppArmFuncPseudo(arma::mat X);
RcppExport SEXP saedevel_cppArmFuncPseudo(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        arma::mat __result = cppArmFuncPseudo(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppArmFuncChol
arma::mat cppArmFuncChol(arma::mat X);
RcppExport SEXP saedevel_cppArmFuncChol(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        arma::mat __result = cppArmFuncChol(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppChol
arma::mat cppChol(arma::mat X);
RcppExport SEXP saedevel_cppChol(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        arma::mat __result = cppChol(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// cppChol2Inv
arma::mat cppChol2Inv(arma::mat X);
RcppExport SEXP saedevel_cppChol2Inv(SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        arma::mat __result = cppChol2Inv(X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// reSAR1
arma::mat reSAR1(arma::mat W, double rho);
RcppExport SEXP saedevel_reSAR1(SEXP WSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP );
        Rcpp::traits::input_parameter< double >::type rho(rhoSEXP );
        arma::mat __result = reSAR1(W, rho);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// reAR1
arma::mat reAR1(int nTime, double rho);
RcppExport SEXP saedevel_reAR1(SEXP nTimeSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type nTime(nTimeSEXP );
        Rcpp::traits::input_parameter< double >::type rho(rhoSEXP );
        arma::mat __result = reAR1(nTime, rho);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// makeBlockDiagonalMat
arma::mat makeBlockDiagonalMat(arma::mat X, int n);
RcppExport SEXP saedevel_makeBlockDiagonalMat(SEXP XSEXP, SEXP nSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< int >::type n(nSEXP );
        arma::mat __result = makeBlockDiagonalMat(X, n);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// matA
arma::mat matA(double sigma2, arma::mat Ome2, int nDomains, arma::colvec sigmaSamplingError);
RcppExport SEXP saedevel_matA(SEXP sigma2SEXP, SEXP Ome2SEXP, SEXP nDomainsSEXP, SEXP sigmaSamplingErrorSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Ome2(Ome2SEXP );
        Rcpp::traits::input_parameter< int >::type nDomains(nDomainsSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type sigmaSamplingError(sigmaSamplingErrorSEXP );
        arma::mat __result = matA(sigma2, Ome2, nDomains, sigmaSamplingError);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// matVinv
Rcpp::List matVinv(arma::mat W, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError);
RcppExport SEXP saedevel_matVinv(SEXP WSEXP, SEXP rho1SEXP, SEXP sigma1SEXP, SEXP rho2SEXP, SEXP sigma2SEXP, SEXP Z1SEXP, SEXP sigmaSamplingErrorSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP );
        Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP );
        Rcpp::traits::input_parameter< double >::type sigma1(sigma1SEXP );
        Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP );
        Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type sigmaSamplingError(sigmaSamplingErrorSEXP );
        Rcpp::List __result = matVinv(W, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// matV
arma::mat matV(arma::mat W, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError);
RcppExport SEXP saedevel_matV(SEXP WSEXP, SEXP rho1SEXP, SEXP sigma1SEXP, SEXP rho2SEXP, SEXP sigma2SEXP, SEXP Z1SEXP, SEXP sigmaSamplingErrorSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP );
        Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP );
        Rcpp::traits::input_parameter< double >::type sigma1(sigma1SEXP );
        Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP );
        Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type sigmaSamplingError(sigmaSamplingErrorSEXP );
        arma::mat __result = matV(W, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// matVinv1
Rcpp::List matVinv1(arma::mat W, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError);
RcppExport SEXP saedevel_matVinv1(SEXP WSEXP, SEXP rho1SEXP, SEXP sigma1SEXP, SEXP rho2SEXP, SEXP sigma2SEXP, SEXP Z1SEXP, SEXP sigmaSamplingErrorSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP );
        Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP );
        Rcpp::traits::input_parameter< double >::type sigma1(sigma1SEXP );
        Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP );
        Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type sigmaSamplingError(sigmaSamplingErrorSEXP );
        Rcpp::List __result = matVinv1(W, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// matP
arma::mat matP(arma::mat solvedV, arma::mat X);
RcppExport SEXP saedevel_matP(SEXP solvedVSEXP, SEXP XSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type solvedV(solvedVSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        arma::mat __result = matP(solvedV, X);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// matVderS1
arma::mat matVderS1(arma::mat Ome1, arma::mat Z1);
RcppExport SEXP saedevel_matVderS1(SEXP Ome1SEXP, SEXP Z1SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type Ome1(Ome1SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP );
        arma::mat __result = matVderS1(Ome1, Z1);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// matVderS2
arma::mat matVderS2(arma::mat Ome2, int nDomains);
RcppExport SEXP saedevel_matVderS2(SEXP Ome2SEXP, SEXP nDomainsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat >::type Ome2(Ome2SEXP );
        Rcpp::traits::input_parameter< int >::type nDomains(nDomainsSEXP );
        arma::mat __result = matVderS2(Ome2, nDomains);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// matVderR1
arma::mat matVderR1(double rho1, double sigma1, arma::mat Z1, arma::mat Ome1, arma::mat W);
RcppExport SEXP saedevel_matVderR1(SEXP rho1SEXP, SEXP sigma1SEXP, SEXP Z1SEXP, SEXP Ome1SEXP, SEXP WSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP );
        Rcpp::traits::input_parameter< double >::type sigma1(sigma1SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Ome1(Ome1SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP );
        arma::mat __result = matVderR1(rho1, sigma1, Z1, Ome1, W);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// matVderR2
arma::mat matVderR2(double rho2, double sigma2, arma::mat Ome2, int nDomains);
RcppExport SEXP saedevel_matVderR2(SEXP rho2SEXP, SEXP sigma2SEXP, SEXP Ome2SEXP, SEXP nDomainsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP );
        Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Ome2(Ome2SEXP );
        Rcpp::traits::input_parameter< int >::type nDomains(nDomainsSEXP );
        arma::mat __result = matVderR2(rho2, sigma2, Ome2, nDomains);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// blue
arma::colvec blue(arma::colvec y, arma::mat X, arma::mat Vinv);
RcppExport SEXP saedevel_blue(SEXP ySEXP, SEXP XSEXP, SEXP VinvSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP );
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Vinv(VinvSEXP );
        arma::colvec __result = blue(y, X, Vinv);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// llr
Rcpp::List llr(arma::colvec y, arma::mat X, double rho1, double sigma1, double rho2, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError, arma::mat W);
RcppExport SEXP saedevel_llr(SEXP ySEXP, SEXP XSEXP, SEXP rho1SEXP, SEXP sigma1SEXP, SEXP rho2SEXP, SEXP sigma2SEXP, SEXP Z1SEXP, SEXP sigmaSamplingErrorSEXP, SEXP WSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP );
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< double >::type rho1(rho1SEXP );
        Rcpp::traits::input_parameter< double >::type sigma1(sigma1SEXP );
        Rcpp::traits::input_parameter< double >::type rho2(rho2SEXP );
        Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type sigmaSamplingError(sigmaSamplingErrorSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP );
        Rcpp::List __result = llr(y, X, rho1, sigma1, rho2, sigma2, Z1, sigmaSamplingError, W);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// optimizerRho
double optimizerRho(arma::colvec rho, arma::colvec y, arma::mat X, double sigma1, double sigma2, arma::mat Z1, arma::colvec sigmaSamplingError, arma::mat W, arma::colvec beta, double K);
RcppExport SEXP saedevel_optimizerRho(SEXP rhoSEXP, SEXP ySEXP, SEXP XSEXP, SEXP sigma1SEXP, SEXP sigma2SEXP, SEXP Z1SEXP, SEXP sigmaSamplingErrorSEXP, SEXP WSEXP, SEXP betaSEXP, SEXP KSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::colvec >::type rho(rhoSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP );
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< double >::type sigma1(sigma1SEXP );
        Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type sigmaSamplingError(sigmaSamplingErrorSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< double >::type K(KSEXP );
        double __result = optimizerRho(rho, y, X, sigma1, sigma2, Z1, sigmaSamplingError, W, beta, K);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// optimizerSigma
arma::colvec optimizerSigma(arma::colvec sigma, arma::colvec rho, arma::colvec y, arma::mat X, arma::mat Z1, arma::colvec sigmaSamplingError, arma::mat W, arma::colvec beta, double K, arma::mat Z);
RcppExport SEXP saedevel_optimizerSigma(SEXP sigmaSEXP, SEXP rhoSEXP, SEXP ySEXP, SEXP XSEXP, SEXP Z1SEXP, SEXP sigmaSamplingErrorSEXP, SEXP WSEXP, SEXP betaSEXP, SEXP KSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::colvec >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type rho(rhoSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP );
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type sigmaSamplingError(sigmaSamplingErrorSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< double >::type K(KSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP );
        arma::colvec __result = optimizerSigma(sigma, rho, y, X, Z1, sigmaSamplingError, W, beta, K, Z);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// optimizeRESTR
arma::colvec optimizeRESTR(arma::colvec sigma, arma::colvec rho, arma::colvec y, arma::mat X, arma::mat Z1, arma::colvec sigmaSamplingError, arma::mat W, arma::colvec beta, int nDomains, int nTime, double K, arma::mat Z, double tol, int maxit);
RcppExport SEXP saedevel_optimizeRESTR(SEXP sigmaSEXP, SEXP rhoSEXP, SEXP ySEXP, SEXP XSEXP, SEXP Z1SEXP, SEXP sigmaSamplingErrorSEXP, SEXP WSEXP, SEXP betaSEXP, SEXP nDomainsSEXP, SEXP nTimeSEXP, SEXP KSEXP, SEXP ZSEXP, SEXP tolSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::colvec >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type rho(rhoSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP );
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type sigmaSamplingError(sigmaSamplingErrorSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< int >::type nDomains(nDomainsSEXP );
        Rcpp::traits::input_parameter< int >::type nTime(nTimeSEXP );
        Rcpp::traits::input_parameter< double >::type K(KSEXP );
        Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP );
        Rcpp::traits::input_parameter< double >::type tol(tolSEXP );
        Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP );
        arma::colvec __result = optimizeRESTR(sigma, rho, y, X, Z1, sigmaSamplingError, W, beta, nDomains, nTime, K, Z, tol, maxit);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// optimizeRER
Rcpp::List optimizeRER(double reVar, arma::colvec vardir, arma::colvec y, arma::mat X, arma::colvec beta, double K, double tol, int maxit);
RcppExport SEXP saedevel_optimizeRER(SEXP reVarSEXP, SEXP vardirSEXP, SEXP ySEXP, SEXP XSEXP, SEXP betaSEXP, SEXP KSEXP, SEXP tolSEXP, SEXP maxitSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type reVar(reVarSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type vardir(vardirSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP );
        Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP );
        Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< double >::type K(KSEXP );
        Rcpp::traits::input_parameter< double >::type tol(tolSEXP );
        Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP );
        Rcpp::List __result = optimizeRER(reVar, vardir, y, X, beta, K, tol, maxit);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}