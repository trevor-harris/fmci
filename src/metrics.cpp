#include<Rcpp.h>
#include<math.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

//' @export
//[[Rcpp::export]]
double annotation_error (const NumericVector& x, const NumericVector& y){

  int nx = x.size();
  int ny = y.size();
  return abs(nx - ny);
}

//' @export
//[[Rcpp::export]]
int hausdorff_dist (const NumericVector& x, const NumericVector& y){

  int nx = x.size();
  int ny = y.size();

  IntegerVector h1(nx);
  IntegerVector h2(ny);

  IntegerVector hmax(2);

  // estimated to actual
  for (int i = 0 ; i < nx; i++) {
    h1(i) = Rcpp::min(abs(x(i) - y));
  }
  hmax(0) = Rcpp::max(h1);

  // actual to estimated
  for (int i = 0 ; i < ny; i++) {
    h2(i) = Rcpp::min(abs(y(i) - x));
  }
  hmax(1) = Rcpp::max(h2);

  // return(Rcpp::max(Rcpp::max(h1), Rcpp::max(h2));
  return(Rcpp::max(hmax));
}

//' @export
//[[Rcpp::export]]
double adjusted_rand (const NumericVector& x, const NumericVector& y, const int& n){

  int nx = x.size();
  int ny = y.size();

  NumericVector r1(nx+2);
  NumericVector r2(ny+2);
  NumericVector r_est(n);
  NumericVector r_act(n);

  // fill out estimated boundaries
  r1(0) = 0;
  r1[Range(1, nx)] = x;
  r1(nx+1) = n;

  // convert estimated cp into regimes
  for (int i = 1 ; i < nx+2; i++) {
    for(int j = r1(i-1); j < r1(i); j++) {
      r_est(j) = i;
    }
  }

  // fill out actual boundaries
  r2(0) = 0;
  r2[Range(1, ny)] = y;
  r2(ny+1) = n;

  // convert actual cp into regimes
  for (int i = 1 ; i < ny+2; i++) {
    for(int j = r2(i-1); j < r2(i); j++) {
      r_act(j) = i;
    }
  }

  // create crosstab
  NumericMatrix tab(nx+1, ny+1);

  for (int i = 0; i < n; i++) {
    tab(r_est(i)-1, r_act(i)-1) += 1;
  }

  // compute ARI
  double a = sum(choose(tab, 2));
  double b = sum(choose(rowSums(tab), 2)) - a;
  double c = sum(choose(colSums(tab), 2)) - a;
  double d = Rf_choose(sum(tab), 2) - a - b - c;
  double ARI = (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d));

  return(ARI);
}



//' @export
//[[Rcpp::export]]
double energydist (const NumericVector& x, const NumericVector& y){

  int nx = x.size();
  int ny = y.size();

  // cross expectation
  double exy = 0.0;
  for (int i = 0 ; i < nx; i++){
    for (int j = 0 ; j < ny; j++){
      exy += abs(x(i) - y(j));
    }
  }
  exy /= (nx * ny);

  // x expectation
  double ex = 0.0;
  for (int i = 0 ; i < nx; i++){
    for (int j = 0 ; j < nx; j++){
      ex += abs(x(i) - x(j));
    }
  }
  ex /= (nx * nx);

  // y expectation
  double ey = 0.0;
  for (int i = 0 ; i < ny; i++){
    for (int j = 0 ; j < ny; j++){
      ey += abs(y(i) - y(j));
    }
  }
  ey /= (ny * ny);

  return (2*exy - ex - ey);
}


// //' @export
// //[[Rcpp::export]]
// double energy_error (const NumericVector& x, const NumericVector& y){
//
//   int nx = x.size();
//   int ny = y.size();
//
//   // cross expectation
//   double exy = 0.0;
//   for (int i = 0 ; i < nx; i++){
//     for (int j = 0 ; j < ny; j++){
//       exy += abs(x(i) - y(j));
//     }
//   }
//   exy /= (nx * ny);
//
//   // x expectation
//   double ex = 0.0;
//   for (int i = 0 ; i < nx; i++){
//     for (int j = 0 ; j < nx; j++){
//       ex += abs(x(i) - x(j));
//     }
//   }
//   ex /= (nx * nx);
//
//   // y expectation
//   double ey = 0.0;
//   for (int i = 0 ; i < ny; i++){
//     for (int j = 0 ; j < ny; j++){
//       ey += abs(y(i) - y(j));
//     }
//   }
//   ey /= (ny * ny);
//
//   return (2*exy - ex - ey) + abs(nx - ny);
// }



