// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// clo
Eigen::MatrixXd clo(Eigen::Map<Eigen::MatrixXd> X);
RcppExport SEXP _RcppCoDA_clo(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(clo(X));
    return rcpp_result_gen;
END_RCPP
}
// center
Eigen::MatrixXd center(Eigen::Map<Eigen::MatrixXd> X);
RcppExport SEXP _RcppCoDA_center(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(center(X));
    return rcpp_result_gen;
END_RCPP
}
// glr
Eigen::MatrixXd glr(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& V);
RcppExport SEXP _RcppCoDA_glr(SEXP XSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(glr(X, V));
    return rcpp_result_gen;
END_RCPP
}
// glrInv
Eigen::MatrixXd glrInv(Eigen::Map<Eigen::MatrixXd>& X, Eigen::Map<Eigen::MatrixXd>& V);
RcppExport SEXP _RcppCoDA_glrInv(SEXP XSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(glrInv(X, V));
    return rcpp_result_gen;
END_RCPP
}
// alrContrast
Eigen::MatrixXd alrContrast(int d, int D, bool inv);
RcppExport SEXP _RcppCoDA_alrContrast(SEXP dSEXP, SEXP DSEXP, SEXP invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< bool >::type inv(invSEXP);
    rcpp_result_gen = Rcpp::wrap(alrContrast(d, D, inv));
    return rcpp_result_gen;
END_RCPP
}
// clrContrast
Eigen::MatrixXd clrContrast(int D, bool inv);
RcppExport SEXP _RcppCoDA_clrContrast(SEXP DSEXP, SEXP invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< bool >::type inv(invSEXP);
    rcpp_result_gen = Rcpp::wrap(clrContrast(D, inv));
    return rcpp_result_gen;
END_RCPP
}
// ilrContrastDefault_internal
Eigen::MatrixXd ilrContrastDefault_internal(int D);
RcppExport SEXP _RcppCoDA_ilrContrastDefault_internal(SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(ilrContrastDefault_internal(D));
    return rcpp_result_gen;
END_RCPP
}
// ilrContrastSign_internal
Eigen::MatrixXd ilrContrastSign_internal(Eigen::Map<Eigen::MatrixXi>& S);
RcppExport SEXP _RcppCoDA_ilrContrastSign_internal(SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXi>& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(ilrContrastSign_internal(S));
    return rcpp_result_gen;
END_RCPP
}
// alr
Eigen::MatrixXd alr(Eigen::Map<Eigen::MatrixXd>& X, Rcpp::Nullable<int> d);
RcppExport SEXP _RcppCoDA_alr(SEXP XSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(alr(X, d));
    return rcpp_result_gen;
END_RCPP
}
// alrInv
Eigen::MatrixXd alrInv(Eigen::Map<Eigen::MatrixXd>& X, Rcpp::Nullable<int> d);
RcppExport SEXP _RcppCoDA_alrInv(SEXP XSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(alrInv(X, d));
    return rcpp_result_gen;
END_RCPP
}
// clr
Eigen::MatrixXd clr(Eigen::Map<Eigen::MatrixXd>& X);
RcppExport SEXP _RcppCoDA_clr(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(clr(X));
    return rcpp_result_gen;
END_RCPP
}
// clrInv
Eigen::MatrixXd clrInv(Eigen::Map<Eigen::MatrixXd>& X);
RcppExport SEXP _RcppCoDA_clrInv(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(clrInv(X));
    return rcpp_result_gen;
END_RCPP
}
// ilr
Eigen::MatrixXd ilr(Eigen::Map<Eigen::MatrixXd>& X, Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd>&> V);
RcppExport SEXP _RcppCoDA_ilr(SEXP XSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd>&> >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(ilr(X, V));
    return rcpp_result_gen;
END_RCPP
}
// ilrInv
Eigen::MatrixXd ilrInv(Eigen::Map<Eigen::MatrixXd>& X, Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd>&> V);
RcppExport SEXP _RcppCoDA_ilrInv(SEXP XSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd>& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Eigen::Map<Eigen::MatrixXd>&> >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(ilrInv(X, V));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RcppCoDA_clo", (DL_FUNC) &_RcppCoDA_clo, 1},
    {"_RcppCoDA_center", (DL_FUNC) &_RcppCoDA_center, 1},
    {"_RcppCoDA_glr", (DL_FUNC) &_RcppCoDA_glr, 2},
    {"_RcppCoDA_glrInv", (DL_FUNC) &_RcppCoDA_glrInv, 2},
    {"_RcppCoDA_alrContrast", (DL_FUNC) &_RcppCoDA_alrContrast, 3},
    {"_RcppCoDA_clrContrast", (DL_FUNC) &_RcppCoDA_clrContrast, 2},
    {"_RcppCoDA_ilrContrastDefault_internal", (DL_FUNC) &_RcppCoDA_ilrContrastDefault_internal, 1},
    {"_RcppCoDA_ilrContrastSign_internal", (DL_FUNC) &_RcppCoDA_ilrContrastSign_internal, 1},
    {"_RcppCoDA_alr", (DL_FUNC) &_RcppCoDA_alr, 2},
    {"_RcppCoDA_alrInv", (DL_FUNC) &_RcppCoDA_alrInv, 2},
    {"_RcppCoDA_clr", (DL_FUNC) &_RcppCoDA_clr, 1},
    {"_RcppCoDA_clrInv", (DL_FUNC) &_RcppCoDA_clrInv, 1},
    {"_RcppCoDA_ilr", (DL_FUNC) &_RcppCoDA_ilr, 2},
    {"_RcppCoDA_ilrInv", (DL_FUNC) &_RcppCoDA_ilrInv, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_RcppCoDA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
