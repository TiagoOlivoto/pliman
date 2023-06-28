#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::RObject help_area(Rcpp::RObject coord) {
  if(!Rcpp::is<Rcpp::List>(coord)) {
    Rcpp::NumericMatrix coord_matrix(coord);
    double area = 0;
    int n = coord_matrix.nrow();
    for (int i = 0; i < n; ++i) {
      int j = i + 1;
      if (j == n) {
        j = 0;
      }
      area += coord_matrix(i, 0) * coord_matrix(j, 1);
      area -= coord_matrix(j, 0) * coord_matrix(i, 1);
    }
    return Rcpp::wrap(std::abs(area / 2));
  } else if (Rcpp::is<Rcpp::List>(coord)) {
    Rcpp::List coord_list(coord);
    std::vector<double> result;
    for (int i = 0; i < coord_list.length(); i++) {
      Rcpp::NumericMatrix coord_matrix = coord_list[i];
      double area = 0;
      int n = coord_matrix.nrow();
      for (int i = 0; i < n; ++i) {
        int j = i + 1;
        if (j == n) {
          j = 0;
        }
        area += coord_matrix(i, 0) * coord_matrix(j, 1);
        area -= coord_matrix(j, 0) * coord_matrix(i, 1);
      }
      result.push_back(std::abs(area / 2));
    }
    return Rcpp::wrap(result);
  } else {
    Rcpp::stop("Invalid input. coord must be a matrix or a list of matrices");
  }
}




// [[Rcpp::export]]
NumericMatrix help_slide(NumericMatrix coord, int fp = 1) {
  int n = coord.nrow();
  NumericMatrix result(n, coord.ncol());
  for(int i = 0; i < n; i++) {
    int j = (i + fp - 1) % n;
    result(i, _) = coord(j, _);
  }
  return result;
}




// [[Rcpp::export]]
NumericVector help_distpts(NumericMatrix data) {
  int n = data.nrow();
  NumericVector distances(n - 1);
  double dx, dy;
  for (int i = 0; i < n - 1; i++) {
    dx = data(i + 1, 0) - data(i, 0);
    dy = data(i + 1, 1) - data(i, 1);
    distances[i] = sqrt(dx * dx + dy * dy);
  }
  return distances;
}





// [[Rcpp::export]]
NumericVector help_centdist(NumericMatrix data) {
  int n = data.nrow();
  int m = data.ncol();
  NumericVector centroid(m);
  NumericVector distances(n);
  // calculate centroid
  for (int j = 0; j < m; j++) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
      sum += data(i, j);
    }
    centroid[j] = sum / n;
  }
  // calculate Euclidean distances
  for (int i = 0; i < n; i++) {
    double dist = 0;
    for (int j = 0; j < m; j++) {
      dist += pow(data(i, j) - centroid[j], 2);
    }
    distances(i) = sqrt(dist);
  }
  return distances;
}



// [[Rcpp::export]]
NumericMatrix help_rotate(NumericMatrix polygon, double angle) {
  int n = polygon.nrow();
  NumericMatrix rotatedPolygon(n, 2);
  double theta = angle * datum::pi / 180;
  arma::mat R = { {cos(theta), -sin(theta)},
  {sin(theta),  cos(theta)} };

  arma::mat P = arma::mat(polygon.begin(), n, 2, false);
  arma::mat RP = P * R;
  rotatedPolygon = wrap(RP);

  return rotatedPolygon;
}



// [[Rcpp::export]]
NumericMatrix help_align(NumericMatrix coord) {
  arma::mat coord_arma(coord.begin(), coord.nrow(), coord.ncol(), false);
  arma::mat var_ = cov(coord_arma);
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, var_);
  arma::mat aligned = coord_arma * eigvec;
  return wrap(flipud(aligned));
}


// [[Rcpp::export]]
arma::mat help_lw(SEXP coord) {
  if(!Rcpp::is<Rcpp::List>(coord)) {
    NumericMatrix coord_matrix(Rcpp::as<Rcpp::NumericMatrix>(coord));
    arma::mat coord_arma(coord_matrix.begin(), coord_matrix.nrow(), coord_matrix.ncol(), false);
    arma::mat var_ = cov(coord_arma);
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, var_);
    arma::mat coordinates = flipud(coord_arma * eigvec);
    double y_min = arma::min(coordinates.col(0));
    double y_max = arma::max(coordinates.col(0));
    double x_min = arma::min(coordinates.col(1));
    double x_max = arma::max(coordinates.col(1));
    arma::vec lw = {x_max - x_min, y_max - y_min};
    return lw.t();
  } else if (Rcpp::is<Rcpp::List>(coord)) {
    List coord_list(coord);
    int n = coord_list.size();
    arma::mat lw_mat(n, 2);
    for (int i = 0; i < coord_list.length(); i++) {
      NumericMatrix coord_matrix = Rcpp::as<NumericMatrix>(coord_list[i]);
      arma::mat coord_arma(coord_matrix.begin(), coord_matrix.nrow(), coord_matrix.ncol(), false);
      arma::mat var_ = cov(coord_arma);
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, var_);
      arma::mat coordinates = flipud(coord_arma * eigvec);
      double y_min = arma::min(coordinates.col(0));
      double y_max = arma::max(coordinates.col(0));
      double x_min = arma::min(coordinates.col(1));
      double x_max = arma::max(coordinates.col(1));
      arma::vec lw = {x_max - x_min, y_max - y_min};
      lw_mat.row(i) = lw.t();
    }
    return lw_mat;
  } else {
    Rcpp::stop("Input must be either a matrix or a list of matrices");
  }
}




// [[Rcpp::export]]
Rcpp::RObject help_eigen_ratio(Rcpp::RObject coord) {
  if(!Rcpp::is<Rcpp::List>(coord)) {
    arma::mat mat_coord = Rcpp::as<arma::mat>(coord);
    arma::mat covmat = cov(mat_coord);
    arma::vec eigenval;
    arma::mat eigenvec;
    eig_sym(eigenval, eigenvec, covmat);
    return Rcpp::wrap(eigenval[0]/eigenval[1]);
  } else if (Rcpp::is<Rcpp::List>(coord)) {
    Rcpp::List list_coord = Rcpp::as<Rcpp::List>(coord);
    int n = list_coord.size();
    arma::vec result(n);
    for (int i = 0; i < n; i++) {
      arma::mat mat_coord = Rcpp::as<arma::mat>(list_coord[i]);
      arma::mat covmat = cov(mat_coord);
      arma::vec eigenval;
      arma::mat eigenvec;
      eig_sym(eigenval, eigenvec, covmat);
      result[i] = eigenval[0]/eigenval[1];
    }
    return Rcpp::wrap(result);
  } else {
    stop("Invalid input. coord must be a matrix or a list of matrices");
  }
}








// Function calculates distance
// between two points


long dist(pair<long, long> p1,
          pair<long, long> p2)
{
  long x0 = p1.first - p2.first;
  long y0 = p1.second - p2.second;
  return sqrt(x0 * x0 + y0 * y0);
}


// [[Rcpp::export]]
Rcpp::RObject help_calliper(Rcpp::RObject coord) {
  if(!Rcpp::is<Rcpp::List>(coord)) {
    NumericMatrix coord_matrix(Rcpp::as<Rcpp::NumericMatrix>(coord));
    int n = coord_matrix.nrow();
    double Max = 0;
    for(int i = 0; i < n; i++)
    {
      for(int j = i + 1; j < n; j++)
      {
        pair<long, long> point1 = make_pair(coord_matrix(i, 0), coord_matrix(i, 1));
        pair<long, long> point2 = make_pair(coord_matrix(j, 0), coord_matrix(j, 1));
        Max = max(Max, (double)dist(point1, point2));
      }
    }
    return Rcpp::wrap(sqrt(Max));
  } else if (Rcpp::is<Rcpp::List>(coord)) {
    List coord_list = as<List>(coord);
    std::vector<double> result;
    for (int i = 0; i < coord_list.length(); i++) {
      NumericMatrix coord_matrix = coord_list[i];
      int n = coord_matrix.nrow();
      double Max = 0;
      for(int i = 0; i < n; i++)
      {
        for(int j = i + 1; j < n; j++)
        {
          pair<long, long> point1 = make_pair(coord_matrix(i, 0), coord_matrix(i, 1));
          pair<long, long> point2 = make_pair(coord_matrix(j, 0), coord_matrix(j, 1));
          Max = max(Max, (double)dist(point1, point2));
        }
      }
      result.push_back(sqrt(Max));
    }
    return Rcpp::wrap(result);
  } else {
    stop("Invalid input. coord must be a matrix or a list of matrices");
  }
}





// [[Rcpp::export]]
Rcpp::RObject help_elongation(Rcpp::RObject coord) {
  if(!Rcpp::is<Rcpp::List>(coord)) {
    NumericMatrix coord_matrix(Rcpp::as<Rcpp::NumericMatrix>(coord));
    arma::mat coord_arma(coord_matrix.begin(), coord_matrix.nrow(), coord_matrix.ncol(), false);
    arma::mat var_ = cov(coord_arma);
    arma::vec eigval;
    arma::vec lw;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, var_);
    arma::mat coordinates = fliplr(coord_arma * eigvec);
    double x_min = arma::min(coordinates.col(0));
    double x_max = arma::max(coordinates.col(0));
    double y_min = arma::min(coordinates.col(1));
    double y_max = arma::max(coordinates.col(1));
    lw = 1 - (y_max - y_min) / (x_max - x_min);
    return Rcpp::wrap(lw);
  } else if (Rcpp::is<Rcpp::List>(coord)) {
    List coord_list(coord);
    int n = coord_list.size();
    arma::vec lw_vec(n);
    for (int i = 0; i < coord_list.length(); i++) {
      NumericMatrix coord_matrix = Rcpp::as<NumericMatrix>(coord_list[i]);
      arma::mat coord_arma(coord_matrix.begin(), coord_matrix.nrow(), coord_matrix.ncol(), false);
      arma::mat var_ = cov(coord_arma);
      arma::vec eigval;
      arma::mat eigvec;
      eig_sym(eigval, eigvec, var_);
      arma::mat coordinates = fliplr(coord_arma * eigvec);
      double x_min = arma::min(coordinates.col(0));
      double x_max = arma::max(coordinates.col(0));
      double y_min = arma::min(coordinates.col(1));
      double y_max = arma::max(coordinates.col(1));
      lw_vec(i) = 1 - (y_max - y_min) / (x_max - x_min);
    }
    return Rcpp::wrap(lw_vec);
  } else {
    Rcpp::stop("Input must be either a matrix or a list of matrices");
  }
}


// [[Rcpp::export]]
arma::mat help_flip_y(arma::mat shape) {
  shape.col(1) = -shape.col(1);
  return shape;
}


// [[Rcpp::export]]
arma::mat help_flip_x(arma::mat shape) {
  shape.col(0) = -shape.col(0);
  return shape;
}


// [[Rcpp::export]]
arma::vec help_mc(const arma::mat& coords) {
  arma::vec center(2);
  center.fill(0);
  double area = 0;
  int n = coords.n_rows;
  for (int i = 0; i < n; i++) {
    double x1 = coords(i, 0);
    double y1 = coords(i, 1);
    double x2 = coords((i + 1) % n, 0);
    double y2 = coords((i + 1) % n, 1);
    double a = x1 * y2 - x2 * y1;
    area += a;
    center(0) += (x1 + x2) * a;
    center(1) += (y1 + y2) * a;
  }
  area /= 2;
  center /= 6 * area;
  return center;
}


// [[Rcpp::export]]
NumericVector help_limits(NumericMatrix mat) {
  int nrow = mat.nrow();
  double minx = mat(0, 0), miny = mat(0, 1), maxx = mat(0, 0), maxy = mat(0, 1);
  for (int i = 0; i < nrow; i++) {
    if (mat(i, 0) < minx) minx = mat(i, 0);
    if (mat(i, 0) > maxx) maxx = mat(i, 0);
    if (mat(i, 1) < miny) miny = mat(i, 1);
    if (mat(i, 1) > maxy) maxy = mat(i, 1);
  }
  NumericVector result(4);
  result[0] = minx;
  result[1] = maxx;
  result[2] = miny;
  result[3] = maxy;
  return result;
}



// [[Rcpp::export]]
NumericVector help_moments(const NumericMatrix & data) {
  int n = data.nrow();
  double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0, sum_xy = 0;

  for (int i = 0; i < n; i++) {
    double x = data(i, 0);
    double y = data(i, 1);
    sum_x += x;
    sum_y += y;
    sum_x2 += x * x;
    sum_y2 += y * y;
    sum_xy += x * y;
  }

  double cov_xy = sum_xy / n - sum_x * sum_y / n / n;
  double var_x = sum_x2 / n - sum_x * sum_x / n / n;
  double var_y = sum_y2 / n - sum_y * sum_y / n / n;

  double theta = 0.5 * atan2(2 * cov_xy, var_x - var_y);
  double a = sqrt(0.5 * (var_x + var_y + sqrt(pow(var_x - var_y, 2) + 4 * pow(cov_xy, 2))));
  double b = sqrt(0.5 * (var_x + var_y - sqrt(pow(var_x - var_y, 2) + 4 * pow(cov_xy, 2))));

  NumericVector result(4);
  result[0] = fmax(a, b);
  result[1] = fmin(a, b);
  result[2] = sqrt(1 - pow(result[1]/result[0], 2));
  result[3] = theta;
  return result;
}


// Rcpp function to count the number of pixels in each connected object
// [[Rcpp::export]]
NumericVector get_area_mask(IntegerVector mask) {
  int n = mask.length();
  int max_value = max(mask);
  NumericVector area(max_value);
  std::fill(area.begin(), area.end(), 0);
  for (int i = 0; i < n; i++) {
    int x = mask[i];
    if (x > 0 && x <= max_value) {
      area[x - 1]++;
    }
  }
  return area;
}





// Function to check if a point is inside a polygon
bool pointInPolygon(NumericMatrix polygon, double x, double y) {
  int n = polygon.nrow();
  bool inside = false;
  for (int i = 0, j = n-1; i < n; j = i++) {
    if (((polygon(i, 1) > y) != (polygon(j, 1) > y)) &&
        (x < (polygon(j, 0) - polygon(i, 0)) * (y - polygon(i, 1)) / (polygon(j, 1) - polygon(i, 1)) + polygon(i, 0))) {
      inside = !inside;
    }
  }
  return inside;
}

// Function to create a binary image from a polygon
// Input: a matrix containing the coordinates of the polygon vertices (each row is a vertex)
// Output: a binary matrix where the polygon is marked as 1 and the remaining area is marked as 0
// Note: the function assumes that the polygon is closed (i.e., the first and last vertices are the same)
// [[Rcpp::export]]
LogicalMatrix polygon_to_binary(NumericMatrix polygon) {
  int xmin = floor(min(polygon(_, 0)));  // Compute minimum x-coordinate
  int ymin = floor(min(polygon(_, 1)));  // Compute minimum y-coordinate
  int xmax = ceil(max(polygon(_, 0)));   // Compute maximum x-coordinate
  int ymax = ceil(max(polygon(_, 1)));   // Compute maximum y-coordinate
  int width = xmax - xmin + 1;           // Compute width
  int height = ymax - ymin + 1;          // Compute height

  LogicalMatrix binaryImage(width, height);

  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      if (pointInPolygon(polygon, i + xmin, j + ymin)) {
        binaryImage(i, j) = true;
      } else {
        binaryImage(i, j) = false;
      }
    }
  }
  return binaryImage;
}



// [[Rcpp::export]]
IntegerVector sum_true_cols(NumericMatrix x) {
  int nrow = x.nrow();
  int ncol = x.ncol();
  IntegerVector col_sums(ncol);

  for (int j = 0; j < ncol; j++) {
    int col_sum = 0;
    for (int i = 0; i < nrow; i++) {
      if (x(i, j) == 1) {
        col_sum++;
      }
    }
    col_sums[j] = col_sum;
  }
  return col_sums;
}


// [[Rcpp::export]]
NumericVector help_poly_angles(NumericMatrix coords) {
  int num_vertices = coords.nrow();
  NumericMatrix sides(num_vertices, num_vertices);
  NumericVector angles(num_vertices);

  // Calculate the sides of the polygon using the distance formula
  for (int i = 0; i < num_vertices; i++) {
    for (int j = 0; j < num_vertices; j++) {
      sides(i,j) = sqrt(pow(coords(i,0) - coords(j,0), 2) + pow(coords(i,1) - coords(j,1), 2));
    }
  }
  // Calculate the internal angles of the polygon using the law of cosines
  for (int i = 0; i < num_vertices; i++) {
    int prev_vertex = (i == 0) ? num_vertices - 1 : i - 1;
    int next_vertex = (i == num_vertices - 1) ? 0 : i + 1;
    angles(i) = acos((pow(sides(prev_vertex,i), 2) + pow(sides(next_vertex,i), 2) - pow(sides(prev_vertex,next_vertex), 2)) /
      (2 * sides(prev_vertex,i) * sides(next_vertex,i)));
  }

  // Convert the angles from radians to degrees
  angles = angles * 180 / M_PI;
  return angles;
}


// [[Rcpp::export]]
NumericMatrix help_smoth(NumericMatrix coords, int niter) {
  int p = coords.nrow();
  NumericMatrix smoothedCoords(p, 2);

  for (int a = 0; a < niter; a++) {
    for (int i = 0; i < p; i++) {
      int prevIndex = (i == 0) ? (p - 1) : (i - 1);
      int nextIndex = (i == p - 1) ? 0 : (i + 1);

      smoothedCoords(i, 0) = (coords(i, 0) + coords(prevIndex, 0) + coords(nextIndex, 0)) / 3.0;
      smoothedCoords(i, 1) = (coords(i, 1) + coords(prevIndex, 1) + coords(nextIndex, 1)) / 3.0;
    }

    coords = smoothedCoords;
  }

  return smoothedCoords;
}
