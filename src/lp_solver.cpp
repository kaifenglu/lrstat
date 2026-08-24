#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

class Simplex {
public:
  Simplex(const std::vector<std::vector<double>> &a,
          const std::vector<double> &b, const std::vector<double> &c)
      : m_(static_cast<int>(b.size())), n_(static_cast<int>(c.size())),
        d_(m_ + 2, std::vector<double>(n_ + 2)), basis_(m_), nonbasis_(n_ + 1) {
    for (int i = 0; i < m_; ++i) {
      for (int j = 0; j < n_; ++j)
        d_[i][j] = a[i][j];
      basis_[i] = n_ + i;
      d_[i][n_] = -1.0;
      d_[i][n_ + 1] = b[i];
    }
    for (int j = 0; j < n_; ++j) {
      nonbasis_[j] = j;
      d_[m_][j] = -c[j];
    }
    nonbasis_[n_] = -1;
    d_[m_ + 1][n_] = 1.0;
  }

  double solve(std::vector<double> &x) {
    int row = 0;
    for (int i = 1; i < m_; ++i) {
      if (d_[i][n_ + 1] < d_[row][n_ + 1])
        row = i;
    }
    if (d_[row][n_ + 1] < -1e-10) {
      pivot(row, n_);
      if (!simplex(1) || d_[m_ + 1][n_ + 1] < -1e-10) {
        throw std::runtime_error("native LP is infeasible");
      }
      if (std::fabs(d_[m_ + 1][n_ + 1]) > 1e-8) {
        throw std::runtime_error("native LP is infeasible");
      }
      if (std::find(basis_.begin(), basis_.end(), -1) != basis_.end()) {
        row = static_cast<int>(std::find(basis_.begin(), basis_.end(), -1) -
                               basis_.begin());
        int col = 0;
        for (int j = 1; j <= n_ + 1; ++j) {
          if (d_[row][j] < d_[row][col] ||
              (std::fabs(d_[row][j] - d_[row][col]) < 1e-12 &&
               nonbasis_[j] < nonbasis_[col])) {
            col = j;
          }
        }
        pivot(row, col);
      }
    }
    if (!simplex(2))
      throw std::runtime_error("native LP is unbounded");
    x.assign(n_, 0.0);
    for (int i = 0; i < m_; ++i) {
      if (basis_[i] < n_)
        x[basis_[i]] = d_[i][n_ + 1];
    }
    return d_[m_][n_ + 1];
  }

private:
  static constexpr double eps_ = 1e-10;
  int m_, n_;
  std::vector<std::vector<double>> d_;
  std::vector<int> basis_, nonbasis_;

  void pivot(int row, int col) {
    double inv = 1.0 / d_[row][col];
    for (int i = 0; i < m_ + 2; ++i) {
      if (i == row)
        continue;
      for (int j = 0; j < n_ + 2; ++j) {
        if (j == col)
          continue;
        d_[i][j] -= d_[row][j] * d_[i][col] * inv;
      }
    }
    for (int j = 0; j < n_ + 2; ++j) {
      if (j != col)
        d_[row][j] *= inv;
    }
    for (int i = 0; i < m_ + 2; ++i) {
      if (i != row)
        d_[i][col] *= -inv;
    }
    d_[row][col] = inv;
    std::swap(basis_[row], nonbasis_[col]);
  }

  bool simplex(int phase) {
    int objective = phase == 1 ? m_ + 1 : m_;
    while (true) {
      int col = -1;
      for (int j = 0; j <= n_; ++j) {
        if (phase == 2 && nonbasis_[j] == -1)
          continue;
        if (col == -1 || d_[objective][j] < d_[objective][col] - eps_ ||
            (std::fabs(d_[objective][j] - d_[objective][col]) <= eps_ &&
             nonbasis_[j] < nonbasis_[col]))
          col = j;
      }
      if (d_[objective][col] >= -eps_)
        return true;
      int row = -1;
      for (int i = 0; i < m_; ++i) {
        if (d_[i][col] <= eps_)
          continue;
        if (row == -1) {
          row = i;
        } else {
          double lhs = d_[i][n_ + 1] / d_[i][col];
          double rhs = d_[row][n_ + 1] / d_[row][col];
          if (lhs < rhs - eps_ ||
              (std::fabs(lhs - rhs) <= eps_ && basis_[i] < basis_[row]))
            row = i;
        }
      }
      if (row == -1)
        return false;
      pivot(row, col);
    }
  }
};

} // namespace

// [[Rcpp::export]]
Rcpp::NumericVector lpMaxEqRcpp(const Rcpp::NumericVector &objective,
                                const Rcpp::NumericMatrix &equality,
                                const Rcpp::NumericVector &rhs) {

  const int m = equality.nrow();
  const int n = equality.ncol();
  if (m != rhs.size() || n != objective.size() || m == 0 || n == 0) {
    throw std::invalid_argument("incompatible native LP dimensions");
  }

  std::vector<std::vector<double>> a(2 * m, std::vector<double>(n));
  std::vector<double> b(2 * m);
  std::vector<double> c(n);
  for (int i = 0; i < m; ++i) {
    if (!std::isfinite(rhs[i]))
      throw std::invalid_argument("non-finite native LP input");
    b[i] = rhs[i];
    b[i + m] = -rhs[i];
    for (int j = 0; j < n; ++j) {
      if (!std::isfinite(equality(i, j)) || !std::isfinite(objective[j])) {
        throw std::invalid_argument("non-finite native LP input");
      }
      a[i][j] = equality(i, j);
      a[i + m][j] = -equality(i, j);
      c[j] = objective[j];
    }
  }

  std::vector<double> solution;
  Simplex solver(a, b, c);
  solver.solve(solution);
  return Rcpp::NumericVector(solution.begin(), solution.end());
}
