// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// which0
IntegerVector which0(NumericVector& x);
RcppExport SEXP _fmci_which0(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(which0(x));
    return rcpp_result_gen;
END_RCPP
}
// sharpenC
NumericVector sharpenC(IntegerVector& cp_ind, NumericVector& w, double cval, int reach);
RcppExport SEXP _fmci_sharpenC(SEXP cp_indSEXP, SEXP wSEXP, SEXP cvalSEXP, SEXP reachSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type cp_ind(cp_indSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type cval(cvalSEXP);
    Rcpp::traits::input_parameter< int >::type reach(reachSEXP);
    rcpp_result_gen = Rcpp::wrap(sharpenC(cp_ind, w, cval, reach));
    return rcpp_result_gen;
END_RCPP
}
// mergeC
List mergeC(IntegerVector& cp_ind, NumericVector& w, int reach);
RcppExport SEXP _fmci_mergeC(SEXP cp_indSEXP, SEXP wSEXP, SEXP reachSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type cp_ind(cp_indSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type reach(reachSEXP);
    rcpp_result_gen = Rcpp::wrap(mergeC(cp_ind, w, reach));
    return rcpp_result_gen;
END_RCPP
}
// poolC
NumericVector poolC(IntegerVector& cp_ind, int reach);
RcppExport SEXP _fmci_poolC(SEXP cp_indSEXP, SEXP reachSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type cp_ind(cp_indSEXP);
    Rcpp::traits::input_parameter< int >::type reach(reachSEXP);
    rcpp_result_gen = Rcpp::wrap(poolC(cp_ind, reach));
    return rcpp_result_gen;
END_RCPP
}
// annotation_error
double annotation_error(const NumericVector& x, const NumericVector& y);
RcppExport SEXP _fmci_annotation_error(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(annotation_error(x, y));
    return rcpp_result_gen;
END_RCPP
}
// hausdorff_dist
int hausdorff_dist(const NumericVector& x, const NumericVector& y);
RcppExport SEXP _fmci_hausdorff_dist(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(hausdorff_dist(x, y));
    return rcpp_result_gen;
END_RCPP
}
// adjusted_rand
double adjusted_rand(const NumericVector& x, const NumericVector& y, const int& n);
RcppExport SEXP _fmci_adjusted_rand(SEXP xSEXP, SEXP ySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(adjusted_rand(x, y, n));
    return rcpp_result_gen;
END_RCPP
}
// energydist
double energydist(const NumericVector& x, const NumericVector& y);
RcppExport SEXP _fmci_energydist(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(energydist(x, y));
    return rcpp_result_gen;
END_RCPP
}
// fmag
NumericVector fmag(const NumericMatrix& fmat);
RcppExport SEXP _fmci_fmag(SEXP fmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type fmat(fmatSEXP);
    rcpp_result_gen = Rcpp::wrap(fmag(fmat));
    return rcpp_result_gen;
END_RCPP
}
// fmean
NumericVector fmean(const NumericMatrix& fmat);
RcppExport SEXP _fmci_fmean(SEXP fmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type fmat(fmatSEXP);
    rcpp_result_gen = Rcpp::wrap(fmean(fmat));
    return rcpp_result_gen;
END_RCPP
}
// ftv
NumericVector ftv(const NumericMatrix& fmat);
RcppExport SEXP _fmci_ftv(SEXP fmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type fmat(fmatSEXP);
    rcpp_result_gen = Rcpp::wrap(ftv(fmat));
    return rcpp_result_gen;
END_RCPP
}
// tv1d
NumericVector tv1d(NumericVector& input, const float lambda);
RcppExport SEXP _fmci_tv1d(SEXP inputSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type input(inputSEXP);
    Rcpp::traits::input_parameter< const float >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(tv1d(input, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fmci_which0", (DL_FUNC) &_fmci_which0, 1},
    {"_fmci_sharpenC", (DL_FUNC) &_fmci_sharpenC, 4},
    {"_fmci_mergeC", (DL_FUNC) &_fmci_mergeC, 3},
    {"_fmci_poolC", (DL_FUNC) &_fmci_poolC, 2},
    {"_fmci_annotation_error", (DL_FUNC) &_fmci_annotation_error, 2},
    {"_fmci_hausdorff_dist", (DL_FUNC) &_fmci_hausdorff_dist, 2},
    {"_fmci_adjusted_rand", (DL_FUNC) &_fmci_adjusted_rand, 3},
    {"_fmci_energydist", (DL_FUNC) &_fmci_energydist, 2},
    {"_fmci_fmag", (DL_FUNC) &_fmci_fmag, 1},
    {"_fmci_fmean", (DL_FUNC) &_fmci_fmean, 1},
    {"_fmci_ftv", (DL_FUNC) &_fmci_ftv, 1},
    {"_fmci_tv1d", (DL_FUNC) &_fmci_tv1d, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_fmci(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}