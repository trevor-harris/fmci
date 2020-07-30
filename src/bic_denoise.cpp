#include<RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

//' @export
//[[Rcpp::export]]
IntegerVector which0 (NumericVector& x) {
  IntegerVector ind = seq(0, x.size()-1);
  return(ind[x != 0]);
}

//' @export
//[[Rcpp::export]]
NumericVector sharpenC (IntegerVector& cp_ind, NumericVector& w, double cval = 0.1, int reach = 10) {

  int ncp = cp_ind.size();
  NumericVector cp_out(ncp);

  int i = 0;
  while(i < ncp) {
    // link all cp within reach to form the group
    int k = 0;
    for(int j = i+1; j < ncp; j++) {
      if((cp_ind(j) - cp_ind(j-1)) < reach) {
        k += 1;
      } else {
        break;
      }
    }
    // if group has a collectively large enough jump then denote the largest one as the cp
    if(std::abs(sum(w[Range(i, i+k)])) > cval) {
      cp_out(i) = cp_ind(which_max(Rcpp::abs(w[Range(i, i+k)])) + i) + 1;
    }
    i += k + 1;

  }
  return(cp_out[which0(cp_out)]);
}

//' @export
//[[Rcpp::export]]
List mergeC (IntegerVector& cp_ind, NumericVector& w, int reach = 10) {

  int ncp = cp_ind.size();
  NumericVector cp_out(ncp);
  NumericVector w_out(ncp);

  int i = 0;
  int ind;
  while(i < ncp) {
    // link all cp within reach to form the group
    int k = 0;
    for(int j = i+1; j < ncp; j++) {
      if((cp_ind(j) - cp_ind(j-1)) < reach) {
        k += 1;
      } else {
        break;
      }
    }

    ind = which_max(abs(w[Range(i, i+k)])) + i;
    cp_out(i) = cp_ind(ind);
    w_out(i) = sum(w[Range(i, i+k)]);

    i += k + 1;
  }

  return(List::create(cp_out[which0(cp_out)],
                      w_out[which0(w_out)]));
  // return(cp_out[which0(cp_out)]);
}

// //' @export
// //[[Rcpp::export]]
// double vmult (const arma::vec& lhs,
//               const arma::vec& rhs) {
//   arma::mat prod = lhs.t() * rhs;
//   return prod(0,0);
// }
//
// //' @export
// //[[Rcpp::export]]
// arma::vec mergeM (const arma::vec& cp, const arma::vec& w, int reach = 10) {
//
//   int ncp = cp.n_elem;
//   arma::vec cp_out(ncp);
//   arma::vec cg(ncp);
//   arma::vec wg(ncp);
//
//
//   int i = 0;
//   int k = 0;
//   while(i < ncp) {
//     // link all cp within reach to form the group
//     k = 0;
//     cg(k) = cp(i);
//     wg(k) = w(i);
//
//     for(int j = i+1; j < ncp; j++) {
//       if((cp(j) - cp(j-1)) < reach) {
//         k += 1;
//         cg(k) = cp(j);
//         wg(k) = w(j);
//
//       } else {
//         break;
//       }
//     }
//     cp_out(i) = vmult(cg, wg) / sum(wg);
//
//     i += k + 1;
//   }
//
//   arma::uvec ids = find(cp_out > 1e-10); // Find indices
//   return(cp_out.elem(ids));
//   // return(cp_out[which0(cp_out)]);
// }


//' @export
//[[Rcpp::export]]
NumericVector poolC (IntegerVector& cp_ind, int reach = 10) {

  int ncp = cp_ind.size();
  NumericVector cp_out(ncp);

  int i = 0;
  // int ind;
  while(i < ncp) {
    // link all cp within reach to form the group
    int k = 0;
    for(int j = i+1; j < ncp; j++) {
      if((cp_ind(j) - cp_ind(j-1)) < reach) {
        k += 1;
      } else {
        break;
      }
    }

    // ind = mean(cp_ind[Range(i, i+k)]) + i;
    cp_out(i) = mean(cp_ind[Range(i, i+k)]);

    i += k + 1;
  }

  return(cp_out[which0(cp_out)]);
  // return(cp_out[which0(cp_out)]);
}



// //' @export
// //[[Rcpp::export]]
// NumericVector l1fit (NumericVector& m, NumericVector& cp) {
//   int n = m.size();
//   int cpn = cp.size();
//   double mhat;
//   NumericVector fit(n);
//
//   // if no cp then set fit to the mean of the m
//   if (cpn == 0) {
//     mhat = mean(m);
//     for(int j = 0; j < n; j++) {
//       fit(j) = mhat;
//     }
//     return(fit);
//   }
//
//   // else set fit to the segment means
//   // fill out regime boundaries
//   NumericVector regb(cpn+2);
//   regb(0) = 0;
//   regb[Range(1, cpn)] = cp;
//   regb(cpn+1) = n;
//
//   // fit based on change points
//   NumericVector mi;
//
//   for (int i = 1; i < (cpn+2); i++) {
//     mi = m[Range(regb(i-1), (regb(i)-1))];
//     mhat = mean(mi);
//
//     for(int j = regb[i-1]; j < regb[i]; j++) {
//       fit(j) = mhat;
//     }
//   }
//   return(fit);
// }

// //' @export
// //[[Rcpp::export]]
// double l1bic (NumericVector& m, NumericVector& cp) {
//   int n = m.size();
//   NumericVector fit(n);
//   double out;
//
//   fit = l1fit(m, cp);
//   out = n*log(var(m - fit)) + sum(diff(fit) != 0)*log(n);
//   return(out);
// }

// //' @export
// //[[Rcpp::export]]
// NumericVector bic_denoise (NumericVector& m, NumericVector& tv, const int r_init) {
//
//   // get cp locatioons and weights
//   NumericVector dtv = diff(tv);
//   IntegerVector cp_ind = which0(dtv);
//   NumericVector w = dtv[cp_ind];
//
//   // kseq are all change point jump sizes in the initial fit
//   // basically just backwards eliminate small change points using the sharpening strat
//   NumericVector adtv = Rcpp::abs(dtv);
//   NumericVector kseq = unique(adtv.sort());
//   int kl = kseq.size();
//
//   IntegerVector rseq = seq(r_init, r_init + kl);
//   NumericVector bics(kl);
//
//   for(int k = 0; k < kl; k++) {
//     NumericVector cp = sharpenC(cp_ind, w, kseq(k), r_init);
//     bics(k) = l1bic(m, cp);
//   }
//   double k = kseq[which_min(bics)];
//
//   for(int r = 0; r < kl; r++) {
//     NumericVector cp = sharpenC(cp_ind, w, k, rseq(r));
//     bics(r) = l1bic(m, cp);
//   }
//   int r = rseq[which_min(bics)];
//
//   NumericVector tvfit = sharpenC(cp_ind, w, k, r);
//   return(l1fit(m, tvfit));
// }
//
//

// //' @export
// //[[Rcpp::export]]
// NumericVector bic_cp (NumericVector& m, NumericVector& tv, const int r_init) {
//
//   // get cp locatioons and weights
//   NumericVector dtv = diff(tv);
//   IntegerVector cp_ind = which0(dtv);
//   NumericVector w = dtv[cp_ind];
//
//   // kseq are all change point jump sizes in the initial fit
//   // basically just backwards eliminate small change points using the sharpening strat
//   NumericVector adtv = unique(Rcpp::abs(dtv));
//   NumericVector kseq = adtv.sort(true);
//   int kl = kseq.size();
//
//   IntegerVector rseq = seq(r_init, r_init + kl);
//   NumericVector bics(kl);
//
//   double k = 0;
//   for(int i = 1; i < kl; i++) {
//     NumericVector cp = sharpenC(cp_ind, w, kseq(i), r_init);
//     bics(i) = l1bic(m, cp);
//
//     if(i > 1 && bics(i) > bics(i-1)) {
//       k = kseq(i-1);
//       break;
//     }
//   }
//   // NumericVector bic_sub = bics[which0(bics)];
//   // double k = kseq[which_min(bic_sub)];
//
//   // for(int r = 0; r < kl; r++) {
//   //   NumericVector cp = sharpenC(cp_ind, w, k, rseq(r));
//   //   bics(r) = l1bic(m, cp);
//   // }
//   // int r = rseq[which_min(bics)];
//
//   // NumericVector cp = sharpenC(cp_ind, w, k, r_init);
//   // NumericVector weight = diff(l1fit(m, cp));
//   //
//   // return(List::create(cp, weight[cp]));
//
//   return(sharpenC(cp_ind, w, k, r_init));
// }
