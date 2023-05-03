// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// threshold_adaptive
Rcpp::NumericMatrix threshold_adaptive(Rcpp::NumericMatrix mat, double k, int windowsize, double maxsd);
RcppExport SEXP _pliman_threshold_adaptive(SEXP matSEXP, SEXP kSEXP, SEXP windowsizeSEXP, SEXP maxsdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type windowsize(windowsizeSEXP);
    Rcpp::traits::input_parameter< double >::type maxsd(maxsdSEXP);
    rcpp_result_gen = Rcpp::wrap(threshold_adaptive(mat, k, windowsize, maxsd));
    return rcpp_result_gen;
END_RCPP
}
// sobel_help
NumericMatrix sobel_help(NumericMatrix A);
RcppExport SEXP _pliman_sobel_help(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(sobel_help(A));
    return rcpp_result_gen;
END_RCPP
}
// rgb_to_hsb_help
NumericMatrix rgb_to_hsb_help(NumericVector r, NumericVector g, NumericVector b);
RcppExport SEXP _pliman_rgb_to_hsb_help(SEXP rSEXP, SEXP gSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g(gSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(rgb_to_hsb_help(r, g, b));
    return rcpp_result_gen;
END_RCPP
}
// rgb_to_srgb_help
arma::mat rgb_to_srgb_help(const arma::mat& rgb);
RcppExport SEXP _pliman_rgb_to_srgb_help(SEXP rgbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type rgb(rgbSEXP);
    rcpp_result_gen = Rcpp::wrap(rgb_to_srgb_help(rgb));
    return rcpp_result_gen;
END_RCPP
}
// help_edge_thinning
NumericMatrix help_edge_thinning(NumericMatrix img);
RcppExport SEXP _pliman_help_edge_thinning(SEXP imgSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type img(imgSEXP);
    rcpp_result_gen = Rcpp::wrap(help_edge_thinning(img));
    return rcpp_result_gen;
END_RCPP
}
// help_dist_transform
NumericMatrix help_dist_transform(const LogicalMatrix& bin);
RcppExport SEXP _pliman_help_dist_transform(SEXP binSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const LogicalMatrix& >::type bin(binSEXP);
    rcpp_result_gen = Rcpp::wrap(help_dist_transform(bin));
    return rcpp_result_gen;
END_RCPP
}
// help_watershed
IntegerMatrix help_watershed(IntegerMatrix binary, IntegerMatrix labels, IntegerMatrix distances);
RcppExport SEXP _pliman_help_watershed(SEXP binarySEXP, SEXP labelsSEXP, SEXP distancesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type binary(binarySEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type distances(distancesSEXP);
    rcpp_result_gen = Rcpp::wrap(help_watershed(binary, labels, distances));
    return rcpp_result_gen;
END_RCPP
}
// help_get_rgb
std::vector<std::vector<double>> help_get_rgb(const NumericMatrix& R, const NumericMatrix& G, const NumericMatrix& B, const IntegerMatrix& labels);
RcppExport SEXP _pliman_help_get_rgb(SEXP RSEXP, SEXP GSEXP, SEXP BSEXP, SEXP labelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type labels(labelsSEXP);
    rcpp_result_gen = Rcpp::wrap(help_get_rgb(R, G, B, labels));
    return rcpp_result_gen;
END_RCPP
}
// bounding_box
IntegerVector bounding_box(LogicalMatrix image, int edge);
RcppExport SEXP _pliman_bounding_box(SEXP imageSEXP, SEXP edgeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< LogicalMatrix >::type image(imageSEXP);
    Rcpp::traits::input_parameter< int >::type edge(edgeSEXP);
    rcpp_result_gen = Rcpp::wrap(bounding_box(image, edge));
    return rcpp_result_gen;
END_RCPP
}
// isolate_objects5
List isolate_objects5(NumericMatrix image, IntegerMatrix labels);
RcppExport SEXP _pliman_isolate_objects5(SEXP imageSEXP, SEXP labelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type image(imageSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type labels(labelsSEXP);
    rcpp_result_gen = Rcpp::wrap(isolate_objects5(image, labels));
    return rcpp_result_gen;
END_RCPP
}
// help_isolate_object
List help_isolate_object(NumericMatrix R, NumericMatrix G, NumericMatrix B, IntegerMatrix labels, bool remove_bg, int edge);
RcppExport SEXP _pliman_help_isolate_object(SEXP RSEXP, SEXP GSEXP, SEXP BSEXP, SEXP labelsSEXP, SEXP remove_bgSEXP, SEXP edgeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type R(RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type B(BSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< bool >::type remove_bg(remove_bgSEXP);
    Rcpp::traits::input_parameter< int >::type edge(edgeSEXP);
    rcpp_result_gen = Rcpp::wrap(help_isolate_object(R, G, B, labels, remove_bg, edge));
    return rcpp_result_gen;
END_RCPP
}
// help_shp
NumericMatrix help_shp(int rows, int cols, NumericVector dims);
RcppExport SEXP _pliman_help_shp(SEXP rowsSEXP, SEXP colsSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< int >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dims(dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(help_shp(rows, cols, dims));
    return rcpp_result_gen;
END_RCPP
}
// help_area
Rcpp::RObject help_area(Rcpp::RObject coord);
RcppExport SEXP _pliman_help_area(SEXP coordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type coord(coordSEXP);
    rcpp_result_gen = Rcpp::wrap(help_area(coord));
    return rcpp_result_gen;
END_RCPP
}
// help_slide
NumericMatrix help_slide(NumericMatrix coord, int fp);
RcppExport SEXP _pliman_help_slide(SEXP coordSEXP, SEXP fpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< int >::type fp(fpSEXP);
    rcpp_result_gen = Rcpp::wrap(help_slide(coord, fp));
    return rcpp_result_gen;
END_RCPP
}
// help_distpts
NumericVector help_distpts(NumericMatrix data);
RcppExport SEXP _pliman_help_distpts(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(help_distpts(data));
    return rcpp_result_gen;
END_RCPP
}
// help_centdist
NumericVector help_centdist(NumericMatrix data);
RcppExport SEXP _pliman_help_centdist(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(help_centdist(data));
    return rcpp_result_gen;
END_RCPP
}
// help_rotate
NumericMatrix help_rotate(NumericMatrix polygon, double angle);
RcppExport SEXP _pliman_help_rotate(SEXP polygonSEXP, SEXP angleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type polygon(polygonSEXP);
    Rcpp::traits::input_parameter< double >::type angle(angleSEXP);
    rcpp_result_gen = Rcpp::wrap(help_rotate(polygon, angle));
    return rcpp_result_gen;
END_RCPP
}
// help_align
NumericMatrix help_align(NumericMatrix coord);
RcppExport SEXP _pliman_help_align(SEXP coordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coord(coordSEXP);
    rcpp_result_gen = Rcpp::wrap(help_align(coord));
    return rcpp_result_gen;
END_RCPP
}
// help_lw
arma::mat help_lw(SEXP coord);
RcppExport SEXP _pliman_help_lw(SEXP coordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type coord(coordSEXP);
    rcpp_result_gen = Rcpp::wrap(help_lw(coord));
    return rcpp_result_gen;
END_RCPP
}
// help_eigen_ratio
Rcpp::RObject help_eigen_ratio(Rcpp::RObject coord);
RcppExport SEXP _pliman_help_eigen_ratio(SEXP coordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type coord(coordSEXP);
    rcpp_result_gen = Rcpp::wrap(help_eigen_ratio(coord));
    return rcpp_result_gen;
END_RCPP
}
// help_calliper
Rcpp::RObject help_calliper(Rcpp::RObject coord);
RcppExport SEXP _pliman_help_calliper(SEXP coordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type coord(coordSEXP);
    rcpp_result_gen = Rcpp::wrap(help_calliper(coord));
    return rcpp_result_gen;
END_RCPP
}
// help_elongation
Rcpp::RObject help_elongation(Rcpp::RObject coord);
RcppExport SEXP _pliman_help_elongation(SEXP coordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type coord(coordSEXP);
    rcpp_result_gen = Rcpp::wrap(help_elongation(coord));
    return rcpp_result_gen;
END_RCPP
}
// help_flip_y
arma::mat help_flip_y(arma::mat shape);
RcppExport SEXP _pliman_help_flip_y(SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(help_flip_y(shape));
    return rcpp_result_gen;
END_RCPP
}
// help_flip_x
arma::mat help_flip_x(arma::mat shape);
RcppExport SEXP _pliman_help_flip_x(SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(help_flip_x(shape));
    return rcpp_result_gen;
END_RCPP
}
// help_mc
arma::vec help_mc(const arma::mat& coords);
RcppExport SEXP _pliman_help_mc(SEXP coordsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type coords(coordsSEXP);
    rcpp_result_gen = Rcpp::wrap(help_mc(coords));
    return rcpp_result_gen;
END_RCPP
}
// help_limits
NumericVector help_limits(NumericMatrix mat);
RcppExport SEXP _pliman_help_limits(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(help_limits(mat));
    return rcpp_result_gen;
END_RCPP
}
// help_moments
NumericVector help_moments(const NumericMatrix& data);
RcppExport SEXP _pliman_help_moments(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(help_moments(data));
    return rcpp_result_gen;
END_RCPP
}
// get_area_mask
NumericVector get_area_mask(IntegerVector mask);
RcppExport SEXP _pliman_get_area_mask(SEXP maskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type mask(maskSEXP);
    rcpp_result_gen = Rcpp::wrap(get_area_mask(mask));
    return rcpp_result_gen;
END_RCPP
}
// polygon_to_binary
LogicalMatrix polygon_to_binary(NumericMatrix polygon);
RcppExport SEXP _pliman_polygon_to_binary(SEXP polygonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type polygon(polygonSEXP);
    rcpp_result_gen = Rcpp::wrap(polygon_to_binary(polygon));
    return rcpp_result_gen;
END_RCPP
}
// sum_true_cols
IntegerVector sum_true_cols(NumericMatrix x);
RcppExport SEXP _pliman_sum_true_cols(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_true_cols(x));
    return rcpp_result_gen;
END_RCPP
}
// help_poly_angles
NumericVector help_poly_angles(NumericMatrix coords);
RcppExport SEXP _pliman_help_poly_angles(SEXP coordsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coords(coordsSEXP);
    rcpp_result_gen = Rcpp::wrap(help_poly_angles(coords));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pliman_threshold_adaptive", (DL_FUNC) &_pliman_threshold_adaptive, 4},
    {"_pliman_sobel_help", (DL_FUNC) &_pliman_sobel_help, 1},
    {"_pliman_rgb_to_hsb_help", (DL_FUNC) &_pliman_rgb_to_hsb_help, 3},
    {"_pliman_rgb_to_srgb_help", (DL_FUNC) &_pliman_rgb_to_srgb_help, 1},
    {"_pliman_help_edge_thinning", (DL_FUNC) &_pliman_help_edge_thinning, 1},
    {"_pliman_help_dist_transform", (DL_FUNC) &_pliman_help_dist_transform, 1},
    {"_pliman_help_watershed", (DL_FUNC) &_pliman_help_watershed, 3},
    {"_pliman_help_get_rgb", (DL_FUNC) &_pliman_help_get_rgb, 4},
    {"_pliman_bounding_box", (DL_FUNC) &_pliman_bounding_box, 2},
    {"_pliman_isolate_objects5", (DL_FUNC) &_pliman_isolate_objects5, 2},
    {"_pliman_help_isolate_object", (DL_FUNC) &_pliman_help_isolate_object, 6},
    {"_pliman_help_shp", (DL_FUNC) &_pliman_help_shp, 3},
    {"_pliman_help_area", (DL_FUNC) &_pliman_help_area, 1},
    {"_pliman_help_slide", (DL_FUNC) &_pliman_help_slide, 2},
    {"_pliman_help_distpts", (DL_FUNC) &_pliman_help_distpts, 1},
    {"_pliman_help_centdist", (DL_FUNC) &_pliman_help_centdist, 1},
    {"_pliman_help_rotate", (DL_FUNC) &_pliman_help_rotate, 2},
    {"_pliman_help_align", (DL_FUNC) &_pliman_help_align, 1},
    {"_pliman_help_lw", (DL_FUNC) &_pliman_help_lw, 1},
    {"_pliman_help_eigen_ratio", (DL_FUNC) &_pliman_help_eigen_ratio, 1},
    {"_pliman_help_calliper", (DL_FUNC) &_pliman_help_calliper, 1},
    {"_pliman_help_elongation", (DL_FUNC) &_pliman_help_elongation, 1},
    {"_pliman_help_flip_y", (DL_FUNC) &_pliman_help_flip_y, 1},
    {"_pliman_help_flip_x", (DL_FUNC) &_pliman_help_flip_x, 1},
    {"_pliman_help_mc", (DL_FUNC) &_pliman_help_mc, 1},
    {"_pliman_help_limits", (DL_FUNC) &_pliman_help_limits, 1},
    {"_pliman_help_moments", (DL_FUNC) &_pliman_help_moments, 1},
    {"_pliman_get_area_mask", (DL_FUNC) &_pliman_get_area_mask, 1},
    {"_pliman_polygon_to_binary", (DL_FUNC) &_pliman_polygon_to_binary, 1},
    {"_pliman_sum_true_cols", (DL_FUNC) &_pliman_sum_true_cols, 1},
    {"_pliman_help_poly_angles", (DL_FUNC) &_pliman_help_poly_angles, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_pliman(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
