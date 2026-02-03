#ifndef __DATAFRAME_LIST__
#define __DATAFRAME_LIST__

// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>

#include "ska/flat_hash_map.hpp"

#include <algorithm>   // fill, min
#include <cmath>       // isnan
#include <cstddef>     // size_t
#include <cstring>     // memcpy
#include <iomanip>     // fixed, setprecision
#include <ios>         // ios::floatfield, ios::fmtflags
#include <iostream>    // cout, ostream
#include <memory>      // make_shared, shared_ptr
#include <stdexcept>   // invalid_argument, out_of_range, runtime_error
#include <string>      // string, to_string
#include <type_traits> // is_same_v
#include <variant>     // get, variant
#include <vector>      // vector
#include <utility>     // move


//
// FlatMatrix: contiguous column-major matrix representation (double)
//
struct FlatMatrix {
  std::vector<double> data; // column-major: element (r,c) => data[c*nrow + r]
  int nrow = 0;
  int ncol = 0;
  
  FlatMatrix() = default;
  FlatMatrix(int nr, int nc) : data(nr * nc), nrow(nr), ncol(nc) {}
  
  FlatMatrix(std::vector<double>&& d, int nr, int nc)
    : data(std::move(d)), nrow(nr), ncol(nc) {
    if (nr * nc != static_cast<int>(data.size()))
      throw std::runtime_error("FlatMatrix: data size mismatch with dimensions");
  }
  
  inline void resize(int nr, int nc) {
    nrow = nr; ncol = nc;
    data.resize(nr * nc);
  }
  
  inline void fill(double v) { std::fill(data.begin(), data.end(), v); }
  inline bool empty() const noexcept { 
    return data.empty() || nrow == 0 || ncol == 0; }
  inline std::size_t size() const noexcept { return data.size(); }
  
  // column-major index helper
  inline static int idx_col(int row, int col, int nrows) noexcept {
    return col * nrows + row;
  }
  
  inline double& operator()(int r, int c) {
    return data[idx_col(r, c, nrow)];
  }
  inline double operator()(int r, int c) const {
    return data[idx_col(r, c, nrow)];
  }
  
  // raw pointer accessors for parallel-friendly use
  inline const double* data_ptr() const noexcept { 
    return data.empty() ? nullptr : data.data(); }
  inline double* data_ptr() noexcept { return data.empty() ? nullptr : data.data(); }
  
  
  // Print helper: pretty-print a small view of the matrix to an ostream 
  // (default std::cout)
  inline void print(std::ostream& os = std::cout, int max_rows = 10, 
                    int max_cols = 10) const {
    os << "FlatMatrix: " << nrow << " x " << ncol << "\n";
    if (nrow == 0 || ncol == 0) return;
    const int rows = std::min(nrow, max_rows);
    const int cols = std::min(ncol, max_cols);
    os.setf(std::ios::fmtflags(0), std::ios::floatfield);
    os << std::fixed << std::setprecision(6);
    for (int r = 0; r < rows; ++r) {
      for (int c = 0; c < cols; ++c) {
        os << (*this)(r, c);
        if (c + 1 < cols) os << "\t";
      }
      if (cols < ncol) os << "\t...";
      os << "\n";
    }
    if (rows < nrow) os << "...\n";
    // restore formatting
    os.unsetf(std::ios::floatfield);
  }
  
  // ostream operator for convenience
  friend inline std::ostream& operator<<(std::ostream& os, const FlatMatrix& fm) {
    fm.print(os);
    return os;
  }
};

//
// IntMatrix: contiguous column-major integer matrix representation (int)
//
struct IntMatrix {
  std::vector<int> data; // column-major: element (r,c) => data[c * nrow + r]
  int nrow = 0;
  int ncol = 0;
  
  IntMatrix() = default;
  IntMatrix(int nr, int nc) : data(nr * nc), nrow(nr), ncol(nc) {}
  
  IntMatrix(std::vector<int>&& d, int nr, int nc)
    : data(std::move(d)), nrow(nr), ncol(nc) {
    if (nr * nc != static_cast<int>(data.size()))
      throw std::runtime_error("IntMatrix: data size mismatch with dimensions");
  }
  
  inline void resize(int nr, int nc) {
    nrow = nr; ncol = nc;
    data.resize(nr * nc);
  }
  
  inline void fill(int v) { std::fill(data.begin(), data.end(), v); }
  inline bool empty() const noexcept { 
    return data.empty() || nrow == 0 || ncol == 0; }
  inline std::size_t size() const noexcept { return data.size(); }
  
  // column-major index helper
  inline static int idx_col(int row, int col, int nrows) noexcept {
    return col * nrows + row;
  }
  
  inline int& operator()(int r, int c) {
    return data[idx_col(r, c, nrow)];
  }
  inline int operator()(int r, int c) const {
    return data[idx_col(r, c, nrow)];
  }
  
  // raw pointer accessors for parallel-friendly use
  inline const int* data_ptr() const noexcept { 
    return data.empty() ? nullptr : data.data(); }
  inline int* data_ptr() noexcept { return data.empty() ? nullptr : data.data(); }
};


//
// FlatArray: contiguous 3-way tensor representation (double)
// Layout: column-major within each 2-D slice, and slices stacked.
// Element (r,c,s) => data[ s * (nrow*ncol) + c * nrow + r ]
//
struct FlatArray {
  std::vector<double> data;
  int nrow = 0;
  int ncol = 0;
  int nslice = 0;
  
  FlatArray() = default;
  FlatArray(int nr, int nc, int ns) : data(nr * nc * ns), nrow(nr), ncol(nc), 
  nslice(ns) {}
  
  FlatArray(std::vector<double>&& d, int nr, int nc, int ns)
    : data(std::move(d)), nrow(nr), ncol(nc), nslice(ns) {
    if (nr * nc * ns != static_cast<int>(data.size()))
      throw std::runtime_error("FlatArray: data size mismatch with dimensions");
  }
  
  inline void resize(int nr, int nc, int ns) {
    nrow = nr; ncol = nc; nslice = ns;
    data.resize(nr * nc * ns);
  }
  
  inline void fill(double v) { std::fill(data.begin(), data.end(), v); }
  inline bool empty() const noexcept { return data.empty() || nrow == 0 || 
    ncol == 0 || nslice == 0; }
  inline std::size_t size() const noexcept { return data.size(); }
  
  // index helper: (row, col, slice)
  inline static int idx(int row, int col, int slice, int nrows, int ncols) noexcept {
    // slice outermost, then column, then row (rows contiguous)
    return slice * (nrows * ncols) + col * nrows + row;
  }
  
  inline double& operator()(int r, int c, int s) {
    return data[idx(r, c, s, nrow, ncol)];
  }
  
  inline double operator()(int r, int c, int s) const {
    return data[idx(r, c, s, nrow, ncol)];
  }
  
  // raw pointer access for parallel-friendly use
  inline const double* data_ptr() const noexcept { 
    return data.empty() ? nullptr : data.data(); }
  inline double* data_ptr() noexcept { 
    return data.empty() ? nullptr : data.data(); }
  
  // pointer to the start of slice s (returns pointer to element (0,0,s))
  inline double* slice_ptr(int s) noexcept {
    if (s < 0 || s >= nslice) return nullptr;
    return data.empty() ? nullptr : data.data() + s * (nrow * ncol);
  }
  inline const double* slice_ptr(int s) const noexcept {
    if (s < 0 || s >= nslice) return nullptr;
    return data.empty() ? nullptr : data.data() + s * (nrow * ncol);
  }
};

//
// DataFrameCpp using ska::flat_hash_map for columns
//
struct DataFrameCpp {
  std::vector<std::string> names_; // insertion order
  
  // use ska::flat_hash_map for better performance than unordered_map
  ska::flat_hash_map<std::string, std::vector<double>> numeric_cols;
  ska::flat_hash_map<std::string, std::vector<int>> int_cols;
  ska::flat_hash_map<std::string, std::vector<unsigned char>> bool_cols;
  ska::flat_hash_map<std::string, std::vector<std::string>> string_cols;
  
  DataFrameCpp() = default;
  
  inline static int idx_col(int row, int col, int nrows) noexcept {
    return FlatMatrix::idx_col(row, col, nrows);
  }
  
  std::size_t nrows() const {
    if (!names_.empty()) {
      const std::string& nm = names_.front();
      if (numeric_cols.count(nm)) return numeric_cols.at(nm).size();
      if (int_cols.count(nm)) return int_cols.at(nm).size();
      if (bool_cols.count(nm)) return bool_cols.at(nm).size();
      if (string_cols.count(nm)) return string_cols.at(nm).size();
    }
    return 0;
  }
  
  std::size_t size() const { return names_.size(); }
  const std::vector<std::string>& names() const { return names_; }
  
  bool containElementNamed(const std::string& name) const {
    return numeric_cols.count(name) || int_cols.count(name) ||
      bool_cols.count(name) || string_cols.count(name);
  }
  
  void check_row_size(std::size_t size, const std::string& name) const {
    std::size_t cur = nrows();
    if (cur > 0 && size != cur) 
      throw std::runtime_error("Column '" + name + 
                               "' has inconsistent number of rows");
  }
  
  // Push single column overloads (double)
  void push_back(const std::vector<double>& col, const std::string& name);
  void push_back(std::vector<double>&& col, const std::string& name);
  
  // int column overloads (const& + &&)
  void push_back(const std::vector<int>& col, const std::string& name);
  void push_back(std::vector<int>&& col, const std::string& name);
  
  // bool column overloads (const& + &&)
  void push_back(const std::vector<unsigned char>& col, const std::string& name);
  void push_back(std::vector<unsigned char>&& col, const std::string& name);
  
  // string column overloads (const& + &&)
  void push_back(const std::vector<std::string>& col, const std::string& name);
  void push_back(std::vector<std::string>&& col, const std::string& name);
  
  // Scalar expansions
  void push_back(double value, const std::string& name);
  void push_back(int value, const std::string& name);
  void push_back(bool value, const std::string& name);
  void push_back(const std::string& value, const std::string& name);
  
  // push_back_flat accepts a column-major flattened buffer with nrows*p values
  // create p new columns named base_name, base_name.1, ..., base_name.p (if p>1)
  void push_back_flat(const std::vector<double>& flat_col_major, int nrows, 
                      const std::string& base_name);
  
  void push_back(const FlatMatrix& fm, const std::string& base_name);
  void push_back(FlatMatrix&& fm, const std::string& base_name);
  
  // Push front variants (const& + && for all types)
  void push_front(const std::vector<double>& col, const std::string& name);
  void push_front(std::vector<double>&& col, const std::string& name);
  
  void push_front(const std::vector<int>& col, const std::string& name);
  void push_front(std::vector<int>&& col, const std::string& name);
  
  void push_front(const std::vector<unsigned char>& col, const std::string& name);
  void push_front(std::vector<unsigned char>&& col, const std::string& name);
  
  void push_front(const std::vector<std::string>& col, const std::string& name);
  void push_front(std::vector<std::string>&& col, const std::string& name);
  
  // Scalar expansions for push_front (requested)
  void push_front(double value, const std::string& name);
  void push_front(int value, const std::string& name);
  void push_front(bool value, const std::string& name);
  void push_front(const std::string& value, const std::string& name);
  
  // Erase column
  void erase(const std::string& name);
  
  // Accessors with type checking
  template <typename T>
  std::vector<T>& get(const std::string& name) {
    if constexpr (std::is_same_v<T, double>) {
      if (numeric_cols.count(name)) return numeric_cols.at(name);
    } else if constexpr (std::is_same_v<T, int>) {
      if (int_cols.count(name)) return int_cols.at(name);
    } else if constexpr (std::is_same_v<T, unsigned char>) {
      if (bool_cols.count(name)) return bool_cols.at(name);
    } else if constexpr (std::is_same_v<T, std::string>) {
      if (string_cols.count(name)) return string_cols.at(name);
    }
    throw std::runtime_error("Column '" + name + "' not found or type mismatch.");
  }
  
  template <typename T>
  const std::vector<T>& get(const std::string& name) const {
    if constexpr (std::is_same_v<T, double>) {
      if (numeric_cols.count(name)) return numeric_cols.at(name);
    } else if constexpr (std::is_same_v<T, int>) {
      if (int_cols.count(name)) return int_cols.at(name);
    } else if constexpr (std::is_same_v<T, unsigned char>) {
      if (bool_cols.count(name)) return bool_cols.at(name);
    } else if constexpr (std::is_same_v<T, std::string>) {
      if (string_cols.count(name)) return string_cols.at(name);
    }
    throw std::runtime_error("Column '" + name + "' not found or type mismatch.");
  }
  
  // Print helper for DataFrameCpp: prints basic table view to ostream
  // (default std::cout)
  inline void print(std::ostream& os = std::cout, int max_rows = 10, 
                    bool show_col_types = false) const {
    const std::size_t rows = nrows();
    const std::size_t cols = size();
    os << "DataFrameCpp: " << rows << " rows x " << cols << " cols\n";
    if (cols == 0) return;
    
    // Print header (column names and optional types)
    for (std::size_t c = 0; c < names_.size(); ++c) {
      const std::string& nm = names_[c];
      os << nm;
      if (show_col_types) {
        os << " (";
        if (numeric_cols.count(nm)) os << "double";
        else if (int_cols.count(nm)) os << "int";
        else if (bool_cols.count(nm)) os << "bool";
        else if (string_cols.count(nm)) os << "string";
        else os << "unknown";
        os << ")";
      }
      if (c + 1 < names_.size()) os << "\t";
    }
    os << "\n";
    
    if (rows == 0) return;
    
    const int rmax = static_cast<int>(std::min<std::size_t>(
      rows, static_cast<std::size_t>(max_rows)));
    
    for (int r = 0; r < rmax; ++r) {
      for (std::size_t c = 0; c < names_.size(); ++c) {
        const std::string& nm = names_[c];
        if (numeric_cols.count(nm)) {
          const auto& col = numeric_cols.at(nm);
          os << col[r];
        } else if (int_cols.count(nm)) {
          const auto& col = int_cols.at(nm);
          os << col[r];
        } else if (bool_cols.count(nm)) {
          const auto& col = bool_cols.at(nm);
          unsigned char v = col[r];
          if (v == 255) os << "NA";
          else os << (v ? "TRUE" : "FALSE");
        } else if (string_cols.count(nm)) {
          const auto& col = string_cols.at(nm);
          os << col[r];
        } else {
          os << "";
        }
        if (c + 1 < names_.size()) os << "\t";
      }
      os << "\n";
    }
    if (rmax < static_cast<int>(rows)) os << "...\n";
  }
  
  // ostream operator for convenience
  friend inline std::ostream& operator<<(std::ostream& os, const DataFrameCpp& df) {
    df.print(os);
    return os;
  }
};

//
// Forward declaration for recursive ListCpp usage
//
struct ListCpp;
using ListPtr = std::shared_ptr<ListCpp>;

//
// ListCpp: heterogeneous container using ska::flat_hash_map for storage
//
struct ListCpp {
  std::vector<std::string> names_; // insertion order
  ska::flat_hash_map<std::string, std::variant<
    double, int, bool, std::string,
    std::vector<double>,
    std::vector<int>,
    std::vector<unsigned char>,
    std::vector<std::string>,
    FlatMatrix,
    IntMatrix,
    FlatArray,
    DataFrameCpp,
    ListPtr,                            // pointer to nested ListCpp
    std::vector<DataFrameCpp>,
    std::vector<ListPtr>                // vector of pointers to nested ListCpp
    >> data;
    
    ListCpp() = default;
    
    std::size_t size() const { return data.size(); }
    std::vector<std::string> names() const { return names_; }
    bool containsElementNamed(const std::string& name) const { 
      return data.count(name) > 0; }
    
    // push_back for value-like types that match the variant alternatives
    template<typename T>
    void push_back(T value, const std::string& name) {
      if (containsElementNamed(name)) 
        throw std::runtime_error("Element '" + name + "' already exists.");
      data.emplace(name, std::move(value));
      names_.push_back(name);
    }
    
    // Convenience overloads for adding ListCpp by value or by shared_ptr
    void push_back(const ListCpp& l, const std::string& name);
    void push_back(ListCpp&& l, const std::string& name);
    void push_back(const ListPtr& p, const std::string& name);
    void push_back(ListPtr&& p, const std::string& name);
    
    // Convenience overloads for FlatMatrix / IntMatrix / FlatArray / DataFrameCpp
    void push_back(const FlatMatrix& fm, const std::string& name);
    void push_back(FlatMatrix&& fm, const std::string& name);
     
    void push_back(const IntMatrix& im, const std::string& name);
    void push_back(IntMatrix&& im, const std::string& name);

    void push_back(const FlatArray& fa, const std::string& name);
    void push_back(FlatArray&& fa, const std::string& name);
    
    void push_back(const DataFrameCpp& df, const std::string& name);
    void push_back(DataFrameCpp&& df, const std::string& name);
    
    // push_front variants
    template<typename T>
    void push_front(T value, const std::string& name) {
      if (containsElementNamed(name)) 
        throw std::runtime_error("Element '" + name + "' already exists.");
      data.emplace(name, std::move(value));
      names_.insert(names_.begin(), name);
    }
    
    void erase(const std::string& name);
    
    template <typename T>
    T& get(const std::string& name) {
      if (!containsElementNamed(name)) 
        throw std::runtime_error("Element with name '" + name + "' not found.");
      return std::get<T>(data.at(name));
    }
    
    template <typename T>
    const T& get(const std::string& name) const {
      if (!containsElementNamed(name)) 
        throw std::runtime_error("Element with name '" + name + "' not found.");
      return std::get<T>(data.at(name));
    }
    
    // Convenience: get nested ListCpp by reference (throws if not present)
    ListCpp& get_list(const std::string& name);
    const ListCpp& get_list(const std::string& name) const;
};

FlatMatrix subset_flatmatrix(const FlatMatrix& fm, const std::vector<int>& row_idx);
void subset_in_place_flatmatrix(FlatMatrix& fm, const std::vector<int>& row_idx);
FlatMatrix subset_flatmatrix(const FlatMatrix& fm, int start, int end);
void subset_in_place_flatmatrix(FlatMatrix& fm, int start, int end);
FlatArray subset_flatarray(const FlatArray& fa, const std::vector<int>& row_idx);
void subset_in_place_flatarray(FlatArray& fa, const std::vector<int>& row_idx);

FlatMatrix concat_flatmatrix(const FlatMatrix& fm1, const FlatMatrix& fm2);
void append_flatmatrix(FlatMatrix& fm1, const FlatMatrix& fm2);
void append_flatmatrix(std::vector<std::vector<double>>& fm1, const FlatMatrix& fm2);
void append_flatmatrix(std::vector<std::vector<double>>& fm1, 
                       const std::vector<std::vector<double>>& fm2);
FlatMatrix cols_to_flatmatrix(const std::vector<std::vector<double>>& cols);
std::vector<double> flatmatrix_get_column(const FlatMatrix& M, int col);
std::vector<int> intmatrix_get_column(const IntMatrix& M, int col);
void flatmatrix_set_column(FlatMatrix& M, int col, const std::vector<double>& src);
void intmatrix_set_column(IntMatrix& M, int col, const std::vector<int>& src);

// ------------------- Converters between R and C++ types (declarations) ----
DataFrameCpp convertRDataFrameToCpp(const Rcpp::DataFrame& r_df);
Rcpp::DataFrame convertDataFrameCppToR(const DataFrameCpp& df);
Rcpp::List convertListCppToR(const ListCpp& L);

// Rcpp wrap specializations
namespace Rcpp {
template <> inline SEXP wrap(const DataFrameCpp& df) {
  return Rcpp::wrap(convertDataFrameCppToR(df));
}
template <> inline SEXP wrap(const ListCpp& l) {
  return Rcpp::wrap(convertListCppToR(l));
}
template <> inline SEXP wrap(const FlatMatrix& fm) {
  if (fm.nrow == 0 || fm.ncol == 0) return R_NilValue;
  Rcpp::NumericMatrix M(fm.nrow, fm.ncol);
  std::memcpy(REAL(M), fm.data.data(), fm.data.size() * sizeof(double));
  return Rcpp::wrap(M);
}
template <> inline SEXP wrap(const IntMatrix& im) {
  if (im.nrow == 0 || im.ncol == 0) return R_NilValue;
  Rcpp::IntegerMatrix M(im.nrow, im.ncol);
  std::memcpy(INTEGER(M), im.data.data(), im.data.size() * sizeof(int));
  return Rcpp::wrap(M);
}
template <> inline SEXP wrap(const FlatArray& fa) {
  if (fa.nrow == 0 || fa.ncol == 0 || fa.nslice == 0) return R_NilValue;
  Rcpp::NumericVector vec(static_cast<R_xlen_t>(fa.data.size()));
  std::memcpy(REAL(vec), fa.data.data(), fa.data.size() * sizeof(double));
  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(fa.nrow, fa.ncol, fa.nslice);
  vec.attr("dim") = dims;
  return Rcpp::wrap(vec);
}
} // namespace Rcpp

// Declarations for helpers implemented in dataframe_list.cpp
FlatMatrix flatmatrix_from_Rmatrix(const Rcpp::NumericMatrix& M);
IntMatrix intmatrix_from_Rmatrix(const Rcpp::IntegerMatrix& M);
std::shared_ptr<ListCpp> listcpp_from_rlist(const Rcpp::List& rlist);

#endif // __DATAFRAME_LIST__