#include <RcppArmadillo.h>
#include <queue>
#include <cmath>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// adapted from imagerExtra  https://bit.ly/3HtxumB
Rcpp::NumericMatrix int_sum(Rcpp::NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  Rcpp::NumericMatrix res(nrow, ncol);

  res(0,0) = mat(0,0);
  for (int i = 1; i < nrow; ++i) {
    res(i,0) = mat(i,0) + res(i-1,0);
  }
  for (int j = 1; j < ncol; ++j) {
    res(0,j) = mat(0,j) + res(0,j-1);
  }
  for (int i = 1; i < nrow; ++i) {
    for (int j = 1; j < ncol; ++j) {
      res(i,j) = mat(i,j) + res(i-1,j) + res(i,j-1) - res(i-1,j-1);
    }
  }
  return res;
}

Rcpp::NumericMatrix int_sum_squared(Rcpp::NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  Rcpp::NumericMatrix mat_squared(nrow, ncol);
  Rcpp::NumericMatrix res(nrow, ncol);

  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      mat_squared(i,j) = mat(i,j) * mat(i,j);
    }
  }

  res(0,0) = mat_squared(0,0);
  for (int i = 1; i < nrow; ++i) {
    res(i,0) = mat_squared(i,0) + res(i-1,0);
  }
  for (int j = 1; j < ncol; ++j) {
    res(0,j) = mat_squared(0,j) + res(0,j-1);
  }
  for (int i = 1; i < nrow; ++i) {
    for (int j = 1; j < ncol; ++j) {
      res(i,j) = mat_squared(i,j) + res(i-1,j) + res(i,j-1) - res(i-1,j-1);
    }
  }
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix threshold_adaptive(Rcpp::NumericMatrix mat, double k, int windowsize, double maxsd) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  Rcpp::NumericMatrix res(nrow, ncol);
  Rcpp::NumericMatrix integ_sum = int_sum(mat);
  Rcpp::NumericMatrix int_sum_sqr = int_sum_squared(mat);
  int winhalf = windowsize / 2;
  int winsize_squared = windowsize * windowsize;
  int nrow_center = nrow - windowsize;
  int ncol_center = ncol - windowsize;


  for (int i = 0; i < winhalf; ++i) {
    for (int j = 0; j < winhalf; ++j) {
      int temp_winsize = (winhalf + i + 1) * (winhalf + j + 1);
      double mean_local = integ_sum(i+winhalf,j+winhalf) / temp_winsize;
      double sd_local = sqrt(int_sum_sqr(i+winhalf,j+winhalf) / temp_winsize - mean_local * mean_local);
      double threshold_local  = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = winhalf; i < nrow_center; ++i) {
    for (int j =0; j < winhalf; ++j) {
      int temp_winsize = windowsize * (winhalf + j + 1);
      double mean_local = (integ_sum(i+winhalf,j+winhalf) - integ_sum(i-winhalf,j+winhalf)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(i+winhalf,j+winhalf) - int_sum_sqr(i-winhalf,j+winhalf)) / temp_winsize - mean_local * mean_local);
      double threshold_local  = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = nrow_center; i < nrow; ++i) {
    for (int j = 0; j < winhalf; ++j) {
      int temp_winsize = (winhalf + nrow - i) * (winhalf + j + 1);
      double mean_local = (integ_sum(nrow-1,j+winhalf) - integ_sum(i-winhalf,j+winhalf)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(nrow-1,j+winhalf) - int_sum_sqr(i-winhalf,j+winhalf)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = 0; i < winhalf; ++i) {
    for (int j = winhalf; j < ncol_center; ++j) {
      int temp_winsize = (winhalf + i + 1) * windowsize;
      double mean_local = (integ_sum(i+winhalf,j+winhalf) - integ_sum(i+winhalf,j-winhalf)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(i+winhalf,j+winhalf) - int_sum_sqr(i+winhalf,j-winhalf)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = winhalf; i < nrow_center; ++i) {
    for (int j = winhalf; j < ncol_center; ++j) {
      double mean_local = (integ_sum(i+winhalf,j+winhalf) + integ_sum(i-winhalf,j-winhalf) - integ_sum(i+winhalf,j-winhalf) - integ_sum(i-winhalf,j+winhalf)) / winsize_squared;
      double sd_local = sqrt((int_sum_sqr(i+winhalf,j+winhalf) + int_sum_sqr(i-winhalf,j-winhalf) - int_sum_sqr(i+winhalf,j-winhalf) - int_sum_sqr(i-winhalf,j+winhalf)) / winsize_squared - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = nrow_center; i < nrow; ++i) {
    for (int j = winhalf; j < ncol_center; ++j) {
      int temp_winsize = (winhalf + nrow - i) * windowsize;
      double mean_local = (integ_sum(nrow-1,j+winhalf) + integ_sum(i-winhalf,j-winhalf) - integ_sum(nrow-1,j-winhalf) - integ_sum(i-winhalf,j+winhalf)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(nrow-1,j+winhalf) + int_sum_sqr(i-winhalf,j-winhalf) - int_sum_sqr(nrow-1,j-winhalf) - int_sum_sqr(i-winhalf,j+winhalf)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = 0; i < winhalf; ++i) {
    for (int j = ncol_center; j < ncol; ++j) {
      int temp_winsize = (winhalf + i + 1) * (winhalf + ncol - j);
      double mean_local = (integ_sum(i+winhalf,ncol-1) - integ_sum(i+winhalf,j-winhalf)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(i+winhalf,ncol-1) - int_sum_sqr(i+winhalf,j-winhalf)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = winhalf; i < nrow_center; ++i) {
    for (int j = ncol_center; j < ncol; ++j) {
      int temp_winsize = windowsize * (winhalf + ncol - j);
      double mean_local = (integ_sum(i+winhalf,ncol-1) + integ_sum(i-winhalf,j-winhalf) - integ_sum(i+winhalf,j-winhalf) - integ_sum(i-winhalf,ncol-1)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(i+winhalf,ncol-1) + int_sum_sqr(i-winhalf,j-winhalf) - int_sum_sqr(i+winhalf,j-winhalf) - int_sum_sqr(i-winhalf,ncol-1)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }

  for (int i = nrow_center; i < nrow; ++i) {
    for (int j = ncol_center; j < ncol; ++j) {
      int temp_winsize = (winhalf + nrow - i) * (winhalf + ncol - j);
      double mean_local = (integ_sum(nrow-1,ncol-1) + integ_sum(i-winhalf,j-winhalf) - integ_sum(nrow-1,j-winhalf) - integ_sum(i-winhalf,ncol-1)) / temp_winsize;
      double sd_local = sqrt((int_sum_sqr(nrow-1,ncol-1) + int_sum_sqr(i-winhalf,j-winhalf) - int_sum_sqr(nrow-1,j-winhalf) - int_sum_sqr(i-winhalf,ncol-1)) / temp_winsize - mean_local * mean_local);
      double threshold_local = mean_local * (1 + k * (sd_local / maxsd - 1));
      if (mat(i,j) <= threshold_local) {
        res(i,j) = 1;
      } else {
        res(i,j) = 0;
      }
    }
  }
  return res;
}


// adapted from https://en.wikipedia.org/wiki/Sobel_operator#MATLAB_implementation
// [[Rcpp::export]]
NumericMatrix sobel_help(NumericMatrix A) {
  NumericMatrix Gx(3, 3);
  NumericMatrix Gy(3, 3);
  Gx(0, 0) = -1; Gx(0, 1) = 0; Gx(0, 2) = 1;
  Gx(1, 0) = -2; Gx(1, 1) = 0; Gx(1, 2) = 2;
  Gx(2, 0) = -1; Gx(2, 1) = 0; Gx(2, 2) = 1;
  Gy(0, 0) = -1; Gy(0, 1) = -2; Gy(0, 2) = -1;
  Gy(1, 0) = 0; Gy(1, 1) = 0; Gy(1, 2) = 0;
  Gy(2, 0) = 1; Gy(2, 1) = 2; Gy(2, 2) = 1;

  int rows = A.nrow();
  int columns = A.ncol();
  NumericMatrix mag(rows, columns);

  for (int i = 0; i < rows - 2; i++) {
    for (int j = 0; j < columns - 2; j++) {
      double S1 = 0;
      double S2 = 0;

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          S1 += Gx(k, l) * A(i + k, j + l);
          S2 += Gy(k, l) * A(i + k, j + l);
        }
      }
      mag(i + 1, j + 1) = sqrt(S1 * S1 + S2 * S2);
    }
  }
  return mag;
}



// [[Rcpp::export]]
NumericMatrix rgb_to_hsb_help(NumericVector r, NumericVector g, NumericVector b) {
  NumericMatrix hsb(r.size(), 3);
  for (int i = 0; i < r.size(); i++) {
    double max_val = std::max(std::max(r[i], g[i]), b[i]);
    double min_val = std::min(std::min(r[i], g[i]), b[i]);
    double diff = max_val - min_val;
    if (max_val == r[i]) {
      hsb(i, 0) = 60 * ((g[i] - b[i]) / diff);
    } else if (max_val == g[i]) {
      hsb(i, 0) = 60 * (2 + (b[i] - r[i]) / diff);
    } else {
      hsb(i, 0) = 60 * (4 + (r[i] - g[i]) / diff);
    }
    hsb(i, 1) = (max_val - min_val) / max_val * 100;
    hsb(i, 2) = max_val * 100;
  }
  return hsb;
}


// [[Rcpp::export]]
arma::mat rgb_to_srgb_help(const arma::mat& rgb) {
  double gamma = 2.2;
  arma::mat rgb_gamma = pow(rgb, gamma);

  arma::mat matrix = {
    { 3.2406, -1.5372, -0.4986 },
    { -0.9689, 1.8758, 0.0415 },
    { 0.0557, -0.2040, 1.0570 }
  };
  arma::mat rgb_srgb = rgb_gamma * matrix;

  rgb_srgb.elem(find(rgb_srgb < 0)).zeros();
  rgb_srgb.elem(find(rgb_srgb > 1)).ones();
  return rgb_srgb;
}



// [[Rcpp::export]]
NumericMatrix help_edge_thinning(NumericMatrix img) {
  int rows = img.nrow();
  int cols = img.ncol();
  NumericMatrix thinned(rows, cols);

  for (int i = 1; i < rows - 1; i++) {
    for (int j = 1; j < cols - 1; j++) {
      int p2 = img(i-1, j);
      int p3 = img(i-1, j+1);
      int p4 = img(i, j+1);
      int p5 = img(i+1, j+1);
      int p6 = img(i+1, j);
      int p7 = img(i+1, j-1);
      int p8 = img(i, j-1);
      int p9 = img(i-1, j-1);

      int A  = (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) +
        (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) +
        (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
        (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1);

      int B  = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
      int m1 = (p2 * p4 * p8);
      int m2 = (p4 * p6 * p8);

      if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0) {
        thinned(i,j) = 0;
      } else {
        thinned(i,j) = img(i,j);
      }
    }
  }
  return thinned;
}


// DISTANCE TRANSFORM
// [[Rcpp::export]]
NumericMatrix help_dist_transform(const LogicalMatrix &bin) {
  int nrow = bin.nrow(), ncol = bin.ncol();
  NumericMatrix dist(nrow, ncol);
  std::fill(dist.begin(), dist.end(), -1.0);
  std::queue<std::pair<int, int> > q;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (bin(i, j)) {
        dist(i, j) = 0.0;
        q.push(std::make_pair(i, j));
      }
    }
  }
  int dx[] = {-1, 0, 1, 0};
  int dy[] = {0, 1, 0, -1};
  while (!q.empty()) {
    int x = q.front().first, y = q.front().second;
    q.pop();
    for (int d = 0; d < 4; d++) {
      int nx = x + dx[d], ny = y + dy[d];
      if (nx >= 0 && nx < nrow && ny >= 0 && ny < ncol && dist(nx, ny) < 0) {
        dist(nx, ny) = dist(x, y) + 1.0;
        q.push(std::make_pair(nx, ny));
      }
    }
  }
  return dist;
}


// WATERSHED
// [[Rcpp::export]]
IntegerMatrix help_watershed(IntegerMatrix binary, IntegerMatrix labels, IntegerMatrix distances) {
  int rows = binary.nrow();
  int cols = binary.ncol();

  // create the output matrix with zeros
  IntegerMatrix output(rows, cols);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      output(i, j) = 0;
    }
  }
  // find the location of all markers and store them in a queue
  std::queue<std::pair<int, int>> markers;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (labels(i, j) != 0) {
        markers.push(std::make_pair(i, j));
      }
    }
  }
  // loop through the markers and perform the watershed segmentation
  while (!markers.empty()) {
    std::pair<int, int> current = markers.front();
    markers.pop();
    int x = current.first;
    int y = current.second;
    int currentLabel = labels(x, y);
    if (x > 0 && binary(x-1, y) == 1 && labels(x-1, y) == 0) {
      labels(x-1, y) = currentLabel;
      markers.push(std::make_pair(x-1, y));
    }
    if (x < rows - 1 && binary(x+1, y) == 1 && labels(x+1, y) == 0) {
      labels(x+1, y) = currentLabel;
      markers.push(std::make_pair(x+1, y));
    }
    if (y > 0 && binary(x, y-1) == 1 && labels(x, y-1) == 0) {
      labels(x, y-1) = currentLabel;markers.push(std::make_pair(x, y-1));
    }
    if (y < cols - 1 && binary(x, y+1) == 1 && labels(x, y+1) == 0) {
      labels(x, y+1) = currentLabel;
      markers.push(std::make_pair(x, y+1));
    }
  }
  // perform the transform distance calculation
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (binary(i, j) == 1) {
        int minDistance = INT_MAX;
        int minLabel = -1;
        if (i > 0 && labels(i-1, j) != 0 && distances(i-1, j) < minDistance) {
          minDistance = distances(i-1, j);
          minLabel = labels(i-1, j);
        }
        if (i < rows - 1 && labels(i+1, j) != 0 && distances(i+1, j) < minDistance) {
          minDistance = distances(i+1, j);
          minLabel = labels(i+1, j);
        }
        if (j > 0 && labels(i, j-1) != 0 && distances(i, j-1) < minDistance) {
          minDistance = distances(i, j-1);
          minLabel = labels(i, j-1);
        }
        if (j < cols - 1 && labels(i, j+1) != 0 && distances(i, j+1) < minDistance) {
          minDistance = distances(i, j+1);
          minLabel = labels(i, j+1);
        }
        if (minLabel != -1) output(i, j) = minLabel;
      }
    }
  }
  return output;
}


// EXTRACT PIXELS
// [[Rcpp::export]]
std::vector<std::vector<double>> help_get_rgb(const NumericMatrix &R, const NumericMatrix &G, const NumericMatrix &B, const IntegerMatrix &labels) {
  int labelsCount = 0;
  int nrow = R.nrow();
  int ncol = R.ncol();
  // get the number of labels
  for (int i = 0; i < nrow * ncol; i++) {
    labelsCount = std::max(labelsCount, labels[i]);
  }
  labelsCount++;

  // create a list to store the RGB values for each label
  std::vector<std::vector<double>> result(labelsCount);

  // loop through each label
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int label = labels(i, j);
      if (label > 0) {
        result[label].push_back(label);
        result[label].push_back(R(i, j));
        result[label].push_back(G(i, j));
        result[label].push_back(B(i, j));
      }
    }
  }
  return result;
}

// EXTRACT RE and NIR
// [[Rcpp::export]]
std::vector<std::vector<double>> help_get_renir(const NumericMatrix &RE, const NumericMatrix &NIR, const IntegerMatrix &labels) {
  int labelsCount = 0;
  int nrow = RE.nrow();
  int ncol = RE.ncol();
  // get the number of labels
  for (int i = 0; i < nrow * ncol; i++) {
    labelsCount = std::max(labelsCount, labels[i]);
  }
  labelsCount++;

  // create a list to store the RGB values for each label
  std::vector<std::vector<double>> result(labelsCount);

  // loop through each label
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      int label = labels(i, j);
      if (label > 0) {
        result[label].push_back(label);
        result[label].push_back(RE(i, j));
        result[label].push_back(NIR(i, j));
      }
    }
  }
  return result;
}

// GET THE COORDINATES OF A BOUNDING BOX OF A BINARY IMAGE
// [[Rcpp::export]]
IntegerVector bounding_box(LogicalMatrix img, int edge) {
  int nrow = img.nrow();
  int ncol = img.ncol();

  int min_row = nrow;
  int max_row = 0;
  int min_col = ncol;
  int max_col = 0;

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      if (img(i, j)) {
        min_row = std::min(min_row, i);
        max_row = std::max(max_row, i);
        min_col = std::min(min_col, j);
        max_col = std::max(max_col, j);
      }
    }
  }
  min_row = std::max(0, min_row - edge);
  max_row = std::min(nrow - 1, max_row + edge);
  min_col = std::max(0, min_col - edge);
  max_col = std::min(ncol - 1, max_col + edge);
  return IntegerVector::create(min_row, max_row, min_col, max_col);
}

// [[Rcpp::export]]
List isolate_objects5(NumericMatrix img, IntegerMatrix labels) {

  int nrows = labels.nrow(), ncols = labels.ncol();

  // Get the unique labels in the mask, starting from 1
  IntegerVector unique_labels = sort_unique(labels);
  unique_labels = unique_labels[unique_labels != 0];

  // Create a list to store the isolated objects
  List isolated_objects(unique_labels.length());

  // Loop over the unique labels in the mask
  for (int i = 0; i < unique_labels.length(); i++) {
    int id = unique_labels[i];

    int top = nrows, bottom = 0, left = ncols, right = 0;

    // Loop over the rows and columns of the labels matrix
    for (int j = 0; j < nrows; j++) {
      for (int k = 0; k < ncols; k++) {
        if (labels(j,k) == id) {
          // Update the bounding box
          top = std::min(top, j);
          bottom = std::max(bottom, j);
          left = std::min(left, k);
          right = std::max(right, k);
        }
      }
    }

    // Crop the object
    int crop_nrows = bottom - top + 1;
    int crop_ncols = right - left + 1;
    NumericMatrix cropped(crop_nrows, crop_ncols);
    for (int j = 0; j < crop_nrows; j++) {
      for (int k = 0; k < crop_ncols; k++) {
        cropped(j,k) = img(top + j, left + k);
      }
    }

    // Add the isolated object to the list
    isolated_objects[i] = cropped;
  }

  return isolated_objects;
}




// HELPER FUNCTION TO ISOLATE OBJECTS BASED ON R-G-B and labels
// [[Rcpp::export]]
List help_isolate_object(NumericMatrix R, NumericMatrix G, NumericMatrix B, IntegerMatrix labels, bool remove_bg, int edge) {

  int nrows = labels.nrow(), ncols = labels.ncol();

  // Get the unique labels in the mask, starting from 1
  IntegerVector unique_labels = sort_unique(labels);
  unique_labels = unique_labels[unique_labels != 0];

  // Create a list to store the isolated objects
  List isolated_objects(unique_labels.length());

  // Loop over the unique labels in the mask
  for (int i = 0; i < unique_labels.length(); i++) {
    int id = unique_labels[i];

    int top = nrows, bottom = 0, left = ncols, right = 0;

    // Loop over the rows and columns of the labels matrix
    for (int j = 0; j < nrows; j++) {
      for (int k = 0; k < ncols; k++) {
        if (labels(j,k) == id) {
          // Update the bounding box
          top = std::min(top, j);
          bottom = std::max(bottom, j);
          left = std::min(left, k);
          right = std::max(right, k);

        }
      }
    }
    // Expand the bounding box by edge pixels
    top = std::max(0, top - edge);
    bottom = std::min(nrows - 1, bottom + edge);
    left = std::max(0, left - edge);
    right = std::min(ncols - 1, right + edge);

    // Crop the objects
    int crop_nrows = bottom - top + 1;
    int crop_ncols = right - left + 1;
    NumericMatrix croppedR(crop_nrows, crop_ncols);
    NumericMatrix croppedG(crop_nrows, crop_ncols);
    NumericMatrix croppedB(crop_nrows, crop_ncols);

    for (int j = 0; j < crop_nrows; j++) {
      for (int k = 0; k < crop_ncols; k++) {
        croppedR(j,k) = R(top + j, left + k);
        croppedG(j,k) = G(top + j, left + k);
        croppedB(j,k) = B(top + j, left + k);
      }
    }

    if(remove_bg){
      // Fill the pixels that are not part of the object with white
      for (int j = 0; j < crop_nrows; j++) {
        for (int k = 0; k < crop_ncols; k++) {
          if (labels(top + j, left + k) != id) {
            croppedR(j,k) = 1;
            croppedG(j,k) = 1;
            croppedB(j,k) = 1;
          }
        }
      }
    }

    // Store the isolated object in the list
    isolated_objects[i] = List::create(croppedR, croppedG, croppedB);
  }

  return isolated_objects;
}

// [[Rcpp::export]]
NumericMatrix help_shp(int rows, int cols, NumericVector dims, double buffer_x, double buffer_y) {
  double xmin = dims[0];
  double xmax = dims[1];
  double ymin = dims[2];
  double ymax = dims[3];
  double xr = xmax - xmin;
  double yr = ymax - ymin;

  double intx = xr / cols;
  double inty = yr / rows;

  NumericMatrix coords(rows * cols * 5, 2);
  int con = 0;

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      con++;

      double x_start = xmin + j * intx;
      double x_end = x_start + intx;
      double y_start = ymin + i * inty;
      double y_end = y_start + inty;

      double buffered_x_start = x_start + buffer_x * intx;
      double buffered_x_end = x_end - buffer_x * intx;
      double buffered_y_start = y_start + buffer_y * inty;
      double buffered_y_end = y_end - buffer_y * inty;

      coords((con - 1) * 5, 0) = buffered_x_start;
      coords((con - 1) * 5, 1) = buffered_y_start;
      coords((con - 1) * 5 + 1, 0) = buffered_x_end;
      coords((con - 1) * 5 + 1, 1) = buffered_y_start;
      coords((con - 1) * 5 + 2, 0) = buffered_x_end;
      coords((con - 1) * 5 + 2, 1) = buffered_y_end;
      coords((con - 1) * 5 + 3, 0) = buffered_x_start;
      coords((con - 1) * 5 + 3, 1) = buffered_y_end;
      coords((con - 1) * 5 + 4, 0) = buffered_x_start;
      coords((con - 1) * 5 + 4, 1) = buffered_y_start;
    }
  }
  return coords;
}





// Function to compute Otsu's threshold
// [[Rcpp::export]]
double help_otsu(const NumericVector& img) {
  int n = img.size();

  double x_max = max(img);
  double x_min = min(img);

  // Compute histogram
  std::vector<int> histogram(256, 0);
  for (int i = 0; i < n; i++) {
    int intensity = (int)(((1 - 0) / (x_max - x_min) * (img[i] - x_max) + 1) * 255);
    histogram[intensity]++;
  }

  // Compute total number of pixels
  int totalPixels = n;

  // Compute sum of intensities
  double sum = 0;
  for (int i = 0; i < 256; i++) {
    sum += i * histogram[i];
  }

  // Compute sum of background intensities
  double sumBackground = 0;
  int backgroundPixels = 0;

  // Initialize variables for storing optimal threshold and maximum between-class variance
  double maxVariance = 0;
  double threshold = 0;


  // Iterate through all possible thresholds
  for (int i = 0; i < 256; i++) {
    // Update background sum and number of background pixels
    backgroundPixels += histogram[i];
    sumBackground += i * histogram[i];

    // Calculate foreground and background weights
    double weightBackground = (double)backgroundPixels / totalPixels;
    double weightForeground = 1 - weightBackground;

    // Calculate mean intensities
    double meanBackground = sumBackground / backgroundPixels;
    double meanForeground = (sum - sumBackground) / (totalPixels - backgroundPixels);

    // Calculate between-class variance
    double variance = weightBackground * weightForeground * pow((meanBackground - meanForeground), 2);

    // Update maximum variance and threshold
    if (variance > maxVariance) {
      maxVariance = variance;
      threshold = i;
    }
  }

  // Scale the threshold value back to the range of 0-1

  return threshold * (x_max - x_min) / 255 + x_min;
}





// Function to apply Guo-Hall thinning algorithm to a binary image
// Adapted from https://observablehq.com/@esperanc/thinning#guoHall
// [[Rcpp::export]]
IntegerMatrix helper_guo_hall(IntegerMatrix image) {
  int wid = image.ncol();
  int hgt = image.nrow();
  IntegerMatrix data2 = Rcpp::clone(image);

  auto get = [&](int col, int row) { return image(row, col) != 0; };
  auto set = [&](int col, int row) { data2(row, col) = 255; };
  auto clear = [&](int col, int row) { data2(row, col) = 0; };

  IntegerMatrix stepCounter(wid, hgt);

  // Performs the conditional removal of one pixel. Even is true
  // if this is an even iteration.
  // Returns 1 if pixel was removed and 0 if not
  auto removePixel = [&](int col, int row, bool even) {
    if (!get(col, row)) return 0; // Not a 1-pixel
    int p2 = get(col - 1, row);
    int p3 = get(col - 1, row + 1);
    int p4 = get(col, row + 1);
    int p5 = get(col + 1, row + 1);
    int p6 = get(col + 1, row);
    int p7 = get(col + 1, row - 1);
    int p8 = get(col, row - 1);
    int p9 = get(col - 1, row - 1);
    int C = ((!p2) & (p3 | p4)) + ((!p4) & (p5 | p6)) + ((!p6) & (p7 | p8)) + ((!p8) & (p9 | p2));
    if (C != 1) return 0;
    int N1 = (p9 | p2) + (p3 | p4) + (p5 | p6) + (p7 | p8);
    int N2 = (p2 | p3) + (p4 | p5) + (p6 | p7) + (p8 | p9);
    int N = (N1 < N2) ? N1 : N2;
    if (N < 2 || N > 3) return 0;
    int m = even ? ((p6 | p7 | (!p9)) & p8) : ((p2 | p3 | (!p5)) & p4);
    if (m == 0) {
      clear(col, row);
      stepCounter(row, col) = 1;
      return 1;
    }
    return 0;
  };

  bool even = true;

  // Performs one thinning step.
  // Returns the number of removed pixels
  auto thinStep = [&]() {
    int result = 0;
    for (int row = 1; row < hgt - 1; row++) {
      for (int col = 1; col < wid - 1; col++) {
        result += removePixel(col, row, even);
      }
    }
    even = !even;
    image = clone(data2); // Copy data2 back to image
    return result;
  };

  // Performs the thinning algorithm
  int n = 0;
  do {
    stepCounter.fill(0);
    n = thinStep();
  } while (n > 0);

  return image;
}

// IDW interpolation function using Rcpp
// [[Rcpp::export]]
NumericVector idw_interpolation_cpp(NumericVector x, NumericVector y, NumericVector values,
                                    NumericVector new_x, NumericVector new_y, double power = 2) {
  // Calculate distances between new points and existing points
  NumericMatrix distances(new_x.size(), x.size());
  for (int i = 0; i < new_x.size(); ++i) {
    distances(i, _) = sqrt(pow(x - new_x[i], 2) + pow(y - new_y[i], 2));
  }

  // Initialize a NumericVector to store results
  NumericVector results(new_x.size(), NA_REAL);

  for (int i = 0; i < new_x.size(); ++i) {
    // Inverse distance weighting formula
    NumericVector weights = 1.0 / pow(distances(i, _), power);
    double weighted_sum = sum(weights * values);
    double total_weight = sum(weights);
    results[i] = total_weight > 0 ? weighted_sum / total_weight : NA_REAL;

  }
  return results;
}
// Function to adjust the bounding box around the centroid
NumericMatrix adjust_bbox(NumericMatrix coords, double width, double height) {
  NumericVector cent = colMeans(coords(Range(0, 3), _));
  double xmin = cent[0] - width / 2;
  double xmax = cent[0] + width / 2;
  double ymin = cent[1] - height / 2;
  double ymax = cent[1] + height / 2;

  NumericMatrix new_bbox(5, 2);
  new_bbox(0, 0) = xmin; new_bbox(0, 1) = ymin;
  new_bbox(1, 0) = xmin; new_bbox(1, 1) = ymax;
  new_bbox(2, 0) = xmax; new_bbox(2, 1) = ymax;
  new_bbox(3, 0) = xmax; new_bbox(3, 1) = ymin;
  new_bbox(4, 0) = xmin; new_bbox(4, 1) = ymin; // Closing the polygon

  return new_bbox;
}

// Function to rotate a polygon around its centroid
NumericMatrix rotate_polygon(NumericMatrix coords, double angle, NumericVector centroid) {
  NumericMatrix rot_mat(2, 2);
  rot_mat(0, 0) = cos(angle);
  rot_mat(0, 1) = -sin(angle);
  rot_mat(1, 0) = sin(angle);
  rot_mat(1, 1) = cos(angle);

  NumericMatrix rotated_coords(coords.nrow(), coords.ncol());

  for (int i = 0; i < coords.nrow(); ++i) {
    NumericVector point = coords(i, _) - centroid;
    NumericVector rotated_point = NumericVector::create(
      rot_mat(0, 0) * point[0] + rot_mat(0, 1) * point[1],
                                                      rot_mat(1, 0) * point[0] + rot_mat(1, 1) * point[1]
    ) + centroid;
    rotated_coords(i, _) = rotated_point;
  }

  return rotated_coords;
}

// [[Rcpp::export]]
List add_width_height_cpp(List grid, double width, double height, NumericVector points_align) {
  int n = grid.size();
  List adjusted_polygons(n);

  // Calculate the rotation angle
  double x1 = points_align[0];
  double y1 = points_align[2];
  double x2 = points_align[1];
  double y2 = points_align[3];
  double rotation_angle = atan2(y2 - y1, x2 - x1);

  for (int i = 0; i < n; ++i) {
    NumericMatrix coords = as<NumericMatrix>(grid[i]);
    NumericMatrix new_bbox = adjust_bbox(coords, width, height);

    // Calculate the centroid of the new bounding box
    NumericVector bbox_centroid = colMeans(new_bbox(Range(0, 3), _));

    // Rotate the bounding box polygon
    NumericMatrix rotated_coords = rotate_polygon(new_bbox, rotation_angle, bbox_centroid);

    adjusted_polygons[i] = rotated_coords;
  }

  return adjusted_polygons;
}
