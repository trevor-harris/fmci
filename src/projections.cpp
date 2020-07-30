// #include<Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

// //' @export
// //[[Rcpp::export]]
// NumericVector fcmd (const NumericMatrix& fmat) {
//
//   int outcol = fmat.ncol();
//   int outrow = fmat.nrow();
//   int fmed = outrow/2;
//
//   NumericVector out(outcol);
//   NumericVector v(outrow);
//   NumericVector v1(fmed);
//   NumericVector v2(fmed);
//
//   for (int i = 0 ; i < outcol; i++){
//     v = fmat(_, i);
//     v1 = v[Range(0, fmed)];
//     v2 = v[Range(fmed+1, outrow)];
//     out(i) = mean(v1) - mean(v2);
//   }
//   return out;
// }


//' @export
//[[Rcpp::export]]
NumericVector fmag (const NumericMatrix& fmat) {

  int outcol = fmat.ncol();
  NumericVector out(outcol);
  // NumericVector m;

  for (int i = 0 ; i < outcol; i++) {
    out(i) = std::sqrt(mean(pow(fmat(_, i), 2)));
  }
  return out;
}


//' @export
//[[Rcpp::export]]
NumericVector fmean (const NumericMatrix& fmat){
  return colMeans(fmat);
}

//' @export
//[[Rcpp::export]]
NumericVector ftv (const NumericMatrix& fmat) {

  int outcol = fmat.ncol();
  NumericVector out(outcol);
  // NumericVector m;

  for (int i = 0 ; i < outcol; i++) {
    out(i) = sum(pow(diff(fmat(_, i)), 2));
  }
  return log(out + 1e-8);
}

// //' @export
// //[[Rcpp::export]]
// arma::vec fline (const arma::mat& fmat, const arma::vec& rhs) {
//   double n = fmat.n_cols;
//   return (fmat * rhs) / n;
// }


