library(RcppArmadillo)
library(inline)

FHcpp <- '
 Rcpp::NumericVector yr(ybar);
 Rcpp::NumericMatrix Xr(Xmean);
 Rcpp::NumericVector Dr(Dvec);

 int  m=Xr.nrow(), p=Xr.ncol();
 arma::mat X(Xr.begin(),m,p,false);
 arma::colvec y(yr.begin(),yr.size(),false);
 arma::colvec D(Dr.begin(),Dr.size(),false);
 arma::mat Xt=arma::trans(X);
 
 double Aold=0;
 double diff=1;
 int k=0;
 int IT=500;
 arma::colvec Vi(m);
 arma::mat XtVi(p,m);
 arma::mat Vimat(m,m);
 arma::mat P(m,m);
 arma::mat Q(p,p);
 arma::colvec Py(m);
 double s;
 double F;
 double Anew=arma::median(D);

 while((diff>0.0001 | diff< -0.0001) & (k<IT)){
   Aold = Anew;
   k = k+1;
     
   for (int i=0; i<m; i++) {
          Vi(i) = 1/(Aold + D[i]);
   }

   Vimat = arma::diagmat(Vi);
   XtVi = arma::trans(Vimat*X);
   Q = arma::inv(XtVi*X);
   P = Vimat-arma::trans(XtVi)*Q*XtVi;
   Py = P*y;
   s = (-0.5) * arma::as_scalar(arma::sum(arma::diagvec(P))) + 0.5 * arma::as_scalar(arma::trans(Py)*Py);
   F = 0.5 * arma::as_scalar(sum(diagvec(P*P)));
   Anew = Aold + s/F;
   diff = (Anew-Aold)/Aold;
 }
 
 double A=0;
 if(Anew > 0) {
  A = Anew;
 }

 for (int i=0; i<m; i++) {
          Vi(i) = 1/(A + D[i]);
 }

 Vimat = arma::diagmat(Vi);
 XtVi = arma::trans(Vimat*X);
 Q = arma::inv(XtVi*X);
 arma::mat beta = Q*XtVi*y;

 arma::mat Xbeta = X * beta;
 arma::colvec resid = y - Xbeta;
 arma::colvec theta = Xbeta + A*(Vi % resid);


 return Rcpp::List::create(Rcpp::Named("Beta") = Rcpp::wrap(beta),
                           Rcpp::Named("Residuals") = Rcpp::wrap(resid),
                           Rcpp::Named("Score") = Rcpp::wrap(s),
                           Rcpp::Named("Iterations") = Rcpp::wrap(k),
                           Rcpp::Named("FisherInformation") = Rcpp::wrap(F),
                           Rcpp::Named("Theta")= Rcpp::wrap(theta),
                           Rcpp::Named("Final Value for A", Rcpp::wrap(A))
                             );
'

FayHeriott <- cxxfunction(signature(ybar="numeric", Xmean="numeric", Dvec="numeric") , body=FHcpp, plugin = "RcppArmadillo" )

# Vergleich mit Implemenation aus sae-Paket anhand von  ausgedachten Daten

D <- 1000
Xmean <- cbind(rnorm(D),(1:D)/D)
beta <- rbind(1,1)

TIME <- matrix(NA,nrow=100,ncol=2)
colnames(TIME) <- c("Molina","Rcpp")

for (r in 1:100){
u <- rnorm(D,0,0.5)
vd <- (1/runif(D,0.5,1))^2
e <- rnorm(D,0,sqrt(vd))


DirEst <- 5 + Xmean %*%beta + u + e

TIME[r,"Molina"] <- system.time(resSAE    <- sae::eblupFH(DirEst ~ Xmean, vd))[3]
TIME[r,"Rcpp"]     <- system.time(resRcpp <- FayHeriott(ybar=DirEst,Xmean=cbind(1,Xmean),Dvec=vd))[3]
}

apply(TIME,2,summary)

# Test ob gleich
all.equal(as.numeric(resSAE$eblup),as.numeric(resRcpp$Theta))

