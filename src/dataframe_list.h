#pragma once

#include <algorithm>   // fill, min
#include <cmath>       // isnan
#include <cstddef>     // size_t
#include <cstring>     // memcpy
#include <iomanip>     // fixed, setprecision
#include <ios>         // ios::floatfield, ios::fmtflags
#include <iostream>    // cout, ostream
#include <limits>      // numeric_limits
#include <memory>      // make_shared, shared_ptr
#include <stdexcept>   // invalid_argument, out_of_range, runtime_error
#include <string>      // string, to_string
#include <type_traits> // is_same_v
#include <variant>     // get, variant
#include <vector>      // vector
#include <utility>     // move

#include <Rcpp.h>
#include "ska/flat_hash_map.hpp"
#include "utilities.h"


//
// FlatMatrix: contiguous column-major matrix representation (double)
//
struct FlatMatrix {
  std::vector<double> data; // column-major: element (r,c) => data[c*nrow + r]
  std::size_t nrow = 0;
  std::size_t ncol = 0;

  FlatMatrix() = default;
  FlatMatrix(std::size_t nr, std::size_t nc) : data(nr * nc), nrow(nr), ncol(nc) {}

  FlatMatrix(std::vector<double>&& d, std::size_t nr, std::size_t nc)
    : data(std::move(d)), nrow(nr), ncol(nc) {
    if (nr * nc != data.size())
      throw std::runtime_error("FlatMatrix: data size mismatch with dimensions");
  }

  inline void resize(std::size_t nr, std::size_t nc) {
    nrow = nr; ncol = nc;
    data.resize(nr * nc);
  }

  inline void fill(double v) { std::fill(data.begin(), data.end(), v); }
  inline bool empty() const noexcept {
    return data.empty() || nrow == 0 || ncol == 0; }
  inline std::size_t size() const noexcept { return data.size(); }

  // column-major index helper
  inline static std::size_t idx_col(std::size_t row, std::size_t col,
                                    std::size_t nrows) noexcept {
    return col * nrows + row;
  }

  inline double& operator()(std::size_t r, std::size_t c) {
    return data[idx_col(r, c, nrow)];
  }
  inline double operator()(std::size_t r, std::size_t c) const {
    return data[idx_col(r, c, nrow)];
  }

  // raw pointer accessors for parallel-friendly use
  inline const double* data_ptr() const noexcept {
    return data.empty() ? nullptr : data.data(); }
  inline double* data_ptr() noexcept { return data.empty() ? nullptr : data.data(); }


  // Print helper: pretty-print a small view of the matrix to an ostream
  // (default std::cout)
  inline void print(std::ostream& os = std::cout, std::size_t max_rows = 10,
                    std::size_t max_cols = 10) const {
    os << "FlatMatrix: " << nrow << " x " << ncol << "\n";
    if (nrow == 0 || ncol == 0) return;
    const std::size_t rows = std::min(nrow, max_rows);
    const std::size_t cols = std::min(ncol, max_cols);
    os.setf(std::ios::fmtflags(0), std::ios::floatfield);
    os << std::fixed << std::setprecision(6);
    for (std::size_t r = 0; r < rows; ++r) {
      for (std::size_t c = 0; c < cols; ++c) {
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
  std::size_t nrow = 0;
  std::size_t ncol = 0;

  IntMatrix() = default;
  IntMatrix(std::size_t nr, std::size_t nc) : data(nr * nc), nrow(nr), ncol(nc) {}

  IntMatrix(std::vector<int>&& d, std::size_t nr, std::size_t nc)
    : data(std::move(d)), nrow(nr), ncol(nc) {
    if (nr * nc != data.size())
      throw std::runtime_error("IntMatrix: data size mismatch with dimensions");
  }

  inline void resize(std::size_t nr, std::size_t nc) {
    nrow = nr; ncol = nc;
    data.resize(nr * nc);
  }

  inline void fill(int v) { std::fill(data.begin(), data.end(), v); }
  inline bool empty() const noexcept {
    return data.empty() || nrow == 0 || ncol == 0; }
  inline std::size_t size() const noexcept { return data.size(); }

  // column-major index helper
  inline static std::size_t idx_col(std::size_t row, std::size_t col,
                                    std::size_t nrows) noexcept {
    return col * nrows + row;
  }

  inline int& operator()(std::size_t r, std::size_t c) {
    return data[idx_col(r, c, nrow)];
  }
  inline int operator()(std::size_t r, std::size_t c) const {
    return data[idx_col(r, c, nrow)];
  }

  // raw pointer accessors for parallel-friendly use
  inline const int* data_ptr() const noexcept {
    return data.empty() ? nullptr : data.data(); }
  inline int* data_ptr() noexcept { return data.empty() ? nullptr : data.data(); }


  // Print helper for IntMatrix: pretty-print a small view of the matrix to an ostream
  // (default std::cout)
  inline void print(std::ostream& os = std::cout, std::size_t max_rows = 10,
                    std::size_t max_cols = 10) const {
    os << "IntMatrix: " << nrow << " x " << ncol << "\n";
    if (nrow == 0 || ncol == 0) return;
    const std::size_t rows = std::min(nrow, max_rows);
    const std::size_t cols = std::min(ncol, max_cols);
    for (std::size_t r = 0; r < rows; ++r) {
      for (std::size_t c = 0; c < cols; ++c) {
        os << (*this)(r, c);
        if (c + 1 < cols) os << "\t";
      }
      if (cols < ncol) os << "\t...";
      os << "\n";
    }
    if (rows < nrow) os << "...\n";
  }

  // ostream operator for convenience
  friend inline std::ostream& operator<<(std::ostream& os, const IntMatrix& im) {
    im.print(os);
    return os;
  }
};

//
// SztMatrix: contiguous column-major std::size_t matrix representation (std::size_t)
//
struct SztMatrix {
  std::vector<std::size_t> data; // column-major: element (r,c) => data[c * nrow + r]
  std::size_t nrow = 0;
  std::size_t ncol = 0;

  SztMatrix() = default;
  SztMatrix(std::size_t nr, std::size_t nc) : data(nr * nc), nrow(nr), ncol(nc) {}

  SztMatrix(std::vector<std::size_t>&& d, std::size_t nr, std::size_t nc)
    : data(std::move(d)), nrow(nr), ncol(nc) {
    if (nr * nc != data.size())
      throw std::runtime_error("SztMatrix: data size mismatch with dimensions");
  }

  inline void resize(std::size_t nr, std::size_t nc) {
    nrow = nr; ncol = nc;
    data.resize(nr * nc);
  }

  inline void fill(std::size_t v) { std::fill(data.begin(), data.end(), v); }
  inline bool empty() const noexcept {
    return data.empty() || nrow == 0 || ncol == 0; }
  inline std::size_t size() const noexcept { return data.size(); }

  // column-major index helper
  inline static std::size_t idx_col(std::size_t row, std::size_t col,
                                    std::size_t nrows) noexcept {
    return col * nrows + row;
  }

  inline std::size_t& operator()(std::size_t r, std::size_t c) {
    return data[idx_col(r, c, nrow)];
  }
  inline std::size_t operator()(std::size_t r, std::size_t c) const {
    return data[idx_col(r, c, nrow)];
  }

  // raw pointer accessors for parallel-friendly use
  inline const std::size_t* data_ptr() const noexcept {
    return data.empty() ? nullptr : data.data(); }
  inline std::size_t* data_ptr() noexcept { return data.empty() ? nullptr :
    data.data(); }


  // Print helper for SztMatrix: pretty-print a small view of the matrix to an ostream
  // (default std::cout)
  inline void print(std::ostream& os = std::cout, std::size_t max_rows = 10,
                    std::size_t max_cols = 10) const {
    os << "SztMatrix: " << nrow << " x " << ncol << "\n";
    if (nrow == 0 || ncol == 0) return;
    const std::size_t rows = std::min(nrow, max_rows);
    const std::size_t cols = std::min(ncol, max_cols);
    for (std::size_t r = 0; r < rows; ++r) {
      for (std::size_t c = 0; c < cols; ++c) {
        os << (*this)(r, c);
        if (c + 1 < cols) os << "\t";
      }
      if (cols < ncol) os << "\t...";
      os << "\n";
    }
    if (rows < nrow) os << "...\n";
  }

  // ostream operator for convenience
  friend inline std::ostream& operator<<(std::ostream& os, const SztMatrix& im) {
    im.print(os);
    return os;
  }
};

//
// BoolMatrix: contiguous column-major logical matrix representation (bool)
//
struct BoolMatrix {
  std::vector<unsigned char> data; // column-major: element(r,c) => data[c * nrow + r]
  std::size_t nrow = 0;
  std::size_t ncol = 0;

  BoolMatrix() = default;
  BoolMatrix(std::size_t nr, std::size_t nc) : data(nr * nc), nrow(nr), ncol(nc) {}

  BoolMatrix(std::vector<unsigned char>&& d, std::size_t nr, std::size_t nc)
    : data(std::move(d)), nrow(nr), ncol(nc) {
    if (nr * nc != data.size())
      throw std::runtime_error("BoolMatrix: data size mismatch with dimensions");
  }

  inline void resize(std::size_t nr, std::size_t nc) {
    nrow = nr; ncol = nc;
    data.resize(nr * nc);
  }

  inline void fill(unsigned char v) { std::fill(data.begin(), data.end(), v); }
  inline bool empty() const noexcept {
    return data.empty() || nrow == 0 || ncol == 0; }
  inline std::size_t size() const noexcept { return data.size(); }

  // column-major index helper
  inline static std::size_t idx_col(std::size_t row, std::size_t col,
                                    std::size_t nrows) noexcept {
    return col * nrows + row;
  }

  inline unsigned char& operator()(std::size_t r, std::size_t c) {
    return data[idx_col(r, c, nrow)];
  }
  inline unsigned char operator()(std::size_t r, std::size_t c) const {
    return data[idx_col(r, c, nrow)];
  }

  // raw pointer accessors for parallel-friendly use
  inline const unsigned char* data_ptr() const noexcept {
    return data.empty() ? nullptr : data.data(); }
  inline unsigned char* data_ptr() noexcept {
    return data.empty() ? nullptr : data.data(); }


  // Print helper for BoolMatrix: pretty-print a small view of the matrix to an ostream
  // Values printed as 0/1. (default std::cout)
  inline void print(std::ostream& os = std::cout, std::size_t max_rows = 10,
                    std::size_t max_cols = 10) const {
    os << "BoolMatrix: " << nrow << " x " << ncol << "\n";
    if (nrow == 0 || ncol == 0) return;
    const std::size_t rows = std::min(nrow, max_rows);
    const std::size_t cols = std::min(ncol, max_cols);
    for (std::size_t r = 0; r < rows; ++r) {
      for (std::size_t c = 0; c < cols; ++c) {
        os << static_cast<int>((*this)(r, c)); // print 0/1
        if (c + 1 < cols) os << "\t";
      }
      if (cols < ncol) os << "\t...";
      os << "\n";
    }
    if (rows < nrow) os << "...\n";
  }

  // ostream operator for convenience
  friend inline std::ostream& operator<<(std::ostream& os, const BoolMatrix& bm) {
    bm.print(os);
    return os;
  }
};


//
// FlatArray: contiguous 3-way tensor representation (double)
// Layout: column-major within each 2-D slice, and slices stacked.
// Element (r,c,s) => data[ s * (nrow*ncol) + c * nrow + r ]
//
struct FlatArray {
  std::vector<double> data;
  std::size_t nrow = 0;
  std::size_t ncol = 0;
  std::size_t nslice = 0;

  FlatArray() = default;
  FlatArray(std::size_t nr, std::size_t nc, std::size_t ns) :
    data(nr * nc * ns), nrow(nr), ncol(nc), nslice(ns) {}

  FlatArray(std::vector<double>&& d, std::size_t nr, std::size_t nc, std::size_t ns)
    : data(std::move(d)), nrow(nr), ncol(nc), nslice(ns) {
    if (nr * nc * ns != data.size())
      throw std::runtime_error("FlatArray: data size mismatch with dimensions");
  }

  inline void resize(std::size_t nr, std::size_t nc, std::size_t ns) {
    nrow = nr; ncol = nc; nslice = ns;
    data.resize(nr * nc * ns);
  }

  inline void fill(double v) { std::fill(data.begin(), data.end(), v); }
  inline bool empty() const noexcept { return data.empty() || nrow == 0 ||
    ncol == 0 || nslice == 0; }
  inline std::size_t size() const noexcept { return data.size(); }

  // index helper: (row, col, slice)
  inline static std::size_t idx(std::size_t row, std::size_t col, std::size_t slice,
                           std::size_t nrows, std::size_t ncols) noexcept {
    // slice outermost, then column, then row (rows contiguous)
    return slice * (nrows * ncols) + col * nrows + row;
  }

  inline double& operator()(std::size_t r, std::size_t c, std::size_t s) {
    return data[idx(r, c, s, nrow, ncol)];
  }

  inline double operator()(std::size_t r, std::size_t c, std::size_t s) const {
    return data[idx(r, c, s, nrow, ncol)];
  }

  // raw pointer access for parallel-friendly use
  inline const double* data_ptr() const noexcept {
    return data.empty() ? nullptr : data.data(); }
  inline double* data_ptr() noexcept {
    return data.empty() ? nullptr : data.data(); }

  // pointer to the start of slice s (returns pointer to element (0,0,s))
  inline double* slice_ptr(std::size_t s) noexcept {
    if (s >= nslice) return nullptr;
    return data.empty() ? nullptr : data.data() + s * (nrow * ncol);
  }
  inline const double* slice_ptr(std::size_t s) const noexcept {
    if (s >= nslice) return nullptr;
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
  ska::flat_hash_map<std::string, std::vector<std::size_t>> size_t_cols;

  DataFrameCpp() = default;

  inline static std::size_t idx_col(std::size_t row, std::size_t col,
                                    std::size_t nrows) noexcept {
    return FlatMatrix::idx_col(row, col, nrows);
  }

  std::size_t nrows() const {
    if (!names_.empty()) {
      const std::string& nm = names_.front();
      if (numeric_cols.count(nm)) return numeric_cols.at(nm).size();
      if (int_cols.count(nm)) return int_cols.at(nm).size();
      if (bool_cols.count(nm)) return bool_cols.at(nm).size();
      if (string_cols.count(nm)) return string_cols.at(nm).size();
      if (size_t_cols.count(nm)) return size_t_cols.at(nm).size();
    }
    return 0;
  }

  std::size_t size() const { return names_.size(); }
  const std::vector<std::string>& names() const { return names_; }

  bool containElementNamed(const std::string& name) const {
    return numeric_cols.count(name) || int_cols.count(name) ||
      bool_cols.count(name) || string_cols.count(name) || size_t_cols.count(name);
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

  // size_t column overloads (const& + &&)
  void push_back(const std::vector<std::size_t>& col, const std::string& name);
  void push_back(std::vector<std::size_t>&& col, const std::string& name);

  // Scalar expansions
  void push_back(double value, const std::string& name);
  void push_back(int value, const std::string& name);
  void push_back(bool value, const std::string& name);
  void push_back(const std::string& value, const std::string& name);
  void push_back(std::size_t value, const std::string& name);

  // push_back_flat accepts a column-major flattened buffer with nrows*p values
  // create p new columns named base_name, base_name.1, ..., base_name.p (if p>1)
  void push_back_flat(const std::vector<double>& flat_col_major, std::size_t nrows,
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

  void push_front(const std::vector<std::size_t>& col, const std::string& name);
  void push_front(std::vector<std::size_t>&& col, const std::string& name);

  // Scalar expansions for push_front (requested)
  void push_front(double value, const std::string& name);
  void push_front(int value, const std::string& name);
  void push_front(bool value, const std::string& name);
  void push_front(const std::string& value, const std::string& name);
  void push_front(std::size_t value, const std::string& name);

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
    } else if constexpr (std::is_same_v<T, std::size_t>) {
      if (size_t_cols.count(name)) return size_t_cols.at(name);
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
    } else if constexpr (std::is_same_v<T, std::size_t>) {
      if (size_t_cols.count(name)) return size_t_cols.at(name);
    }
    throw std::runtime_error("Column '" + name + "' not found or type mismatch.");
  }

  // Print helper for DataFrameCpp: prints basic table view to ostream
  // (default std::cout)
  inline void print(std::ostream& os = std::cout, std::size_t max_rows = 10,
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
        else if (size_t_cols.count(nm)) os << "size_t";
        else os << "unknown";
        os << ")";
      }
      if (c + 1 < names_.size()) os << "\t";
    }
    os << "\n";

    if (rows == 0) return;

    const std::size_t rmax = std::min(rows, max_rows);

    for (std::size_t r = 0; r < rmax; ++r) {
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
        } else if (size_t_cols.count(nm)) {
          const auto& col = size_t_cols.at(nm);
          os << col[r];
        } else {
          os << "";
        }
        if (c + 1 < names_.size()) os << "\t";
      }
      os << "\n";
    }
    if (rmax < rows) os << "...\n";
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
    double, int, bool, std::string, std::size_t,
    std::vector<double>,
    std::vector<int>,
    std::vector<unsigned char>,
    std::vector<std::string>,
    std::vector<std::size_t>,
    FlatMatrix,
    IntMatrix,
    SztMatrix,
    BoolMatrix,
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

    // Overloads for FlatMatrix / IntMatrix / BoolMatrix / FlatArray / DataFrameCpp
    void push_back(const FlatMatrix& fm, const std::string& name);
    void push_back(FlatMatrix&& fm, const std::string& name);

    void push_back(const IntMatrix& im, const std::string& name);
    void push_back(IntMatrix&& im, const std::string& name);

    void push_back(const SztMatrix& sm, const std::string& name);
    void push_back(SztMatrix&& sm, const std::string& name);

    void push_back(const BoolMatrix& im, const std::string& name);
    void push_back(BoolMatrix&& im, const std::string& name);

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

FlatMatrix subset_flatmatrix(const FlatMatrix& fm,
                             const std::vector<std::size_t>& row_idx);
void subset_in_place_flatmatrix(FlatMatrix& fm,
                                const std::vector<std::size_t>& row_idx);
FlatMatrix subset_flatmatrix(const FlatMatrix& fm, std::size_t start, std::size_t end);
void subset_in_place_flatmatrix(FlatMatrix& fm, std::size_t start, std::size_t end);
FlatMatrix subset_flatmatrix(const FlatMatrix& fm,
                             std::size_t row_start, std::size_t row_end,
                             std::size_t col_start, std::size_t col_end);
FlatArray subset_flatarray(const FlatArray& fa,
                           const std::vector<std::size_t>& row_idx);
void subset_in_place_flatarray(FlatArray& fa, const std::vector<std::size_t>& row_idx);

FlatMatrix concat_flatmatrix(const FlatMatrix& fm1, const FlatMatrix& fm2);
void append_flatmatrix(FlatMatrix& fm1, const FlatMatrix& fm2);
void append_flatmatrix(std::vector<std::vector<double>>& fm1, const FlatMatrix& fm2);
void append_flatmatrix(std::vector<std::vector<double>>& fm1,
                       const std::vector<std::vector<double>>& fm2);
FlatMatrix cols_to_flatmatrix(const std::vector<std::vector<double>>& cols);
std::vector<double> flatmatrix_get_column(const FlatMatrix& M, std::size_t col);
DoubleView flatmatrix_get_column_view(const FlatMatrix& M, std::size_t col);
std::vector<int> intmatrix_get_column(const IntMatrix& M, std::size_t col);
std::vector<std::size_t> sztmatrix_get_column(const SztMatrix& M, std::size_t col);
std::vector<unsigned char> boolmatrix_get_column(const BoolMatrix& M, std::size_t col);
void flatmatrix_set_column(FlatMatrix& M, std::size_t col,
                           const std::vector<double>& src);
void intmatrix_set_column(IntMatrix& M, std::size_t col, const std::vector<int>& src);
void sztmatrix_set_column(SztMatrix& M, std::size_t col,
                          const std::vector<std::size_t>& src);
void boolmatrix_set_column(BoolMatrix& M, std::size_t col,
                           const std::vector<unsigned char>& src);

// ------------------- Converters between R and C++ types (declarations) ----
DataFrameCpp convertRDataFrameToCpp(const Rcpp::DataFrame& r_df);
Rcpp::DataFrame convertDataFrameCppToR(const DataFrameCpp& df);
Rcpp::List convertListCppToR(const ListCpp& L);

// Rcpp wrap specializations
namespace Rcpp {
template <> inline SEXP wrap(const DataFrameCpp& df) {
  return convertDataFrameCppToR(df);
}
template <> inline SEXP wrap(const ListCpp& l) {
  return convertListCppToR(l);
}
template <> inline SEXP wrap(const FlatMatrix& fm) {
  if (fm.nrow == 0 || fm.ncol == 0) return R_NilValue;
  Rcpp::NumericMatrix M(fm.nrow, fm.ncol);
  std::memcpy(REAL(M), fm.data.data(), fm.data.size() * sizeof(double));
  return M;
}
template <> inline SEXP wrap(const IntMatrix& im) {
  if (im.nrow == 0 || im.ncol == 0) return R_NilValue;
  Rcpp::IntegerMatrix M(im.nrow, im.ncol);
  std::memcpy(INTEGER(M), im.data.data(), im.data.size() * sizeof(int));
  return M;
}
template <> inline SEXP wrap(const SztMatrix& sm) {
  if (sm.nrow == 0 || sm.ncol == 0) return R_NilValue;
  // Create an R integer matrix with same dims
  Rcpp::IntegerMatrix M(sm.nrow, sm.ncol);
  int* dst = INTEGER(M);
  const std::size_t* src = sm.data.data();
  const std::size_t N = sm.data.size();
  static bool warned_size_t_overflow = false;
  const std::size_t INTMAX = static_cast<std::size_t>(
    std::numeric_limits<int>::max());
  for (std::size_t i = 0; i < N; ++i) {
    std::size_t v = src[i];
    if (v > INTMAX) {
      dst[i] = NA_INTEGER;
      if (!warned_size_t_overflow) {
        Rcpp::warning("SztMatrix -> IntegerMatrix: some std::size_t values "
                        "exceed INT_MAX and were converted to NA_integer_");
        warned_size_t_overflow = true;
      }
    } else {
      dst[i] = static_cast<int>(v);
    }
  }
  return M;
}
template <> inline SEXP wrap(const BoolMatrix& im) {
  if (im.nrow == 0 || im.ncol == 0) return R_NilValue;
  Rcpp::LogicalMatrix M(im.nrow, im.ncol);
  int* dst = LOGICAL(M);
  const unsigned char* src = im.data.data();
  // Convert: 0->FALSE, 1->TRUE, 255->NA_LOGICAL
  for (std::size_t i = 0; i < im.data.size(); ++i) {
    if (src[i] == 255) {
      dst[i] = NA_LOGICAL;
    } else {
      dst[i] = static_cast<int>(src[i]);
    }
  }
  return M;
}
template <> inline SEXP wrap(const FlatArray& fa) {
  if (fa.nrow == 0 || fa.ncol == 0 || fa.nslice == 0) return R_NilValue;
  Rcpp::NumericVector vec(static_cast<R_xlen_t>(fa.data.size()));
  std::memcpy(REAL(vec), fa.data.data(), fa.data.size() * sizeof(double));
  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(fa.nrow, fa.ncol, fa.nslice);
  vec.attr("dim") = dims;
  return vec;
}
} // namespace Rcpp

// Declarations for helpers implemented in dataframe_list.cpp
FlatMatrix flatmatrix_from_Rmatrix(const Rcpp::NumericMatrix& M);
IntMatrix intmatrix_from_Rmatrix(const Rcpp::IntegerMatrix& M);
BoolMatrix boolmatrix_from_Rmatrix(const Rcpp::LogicalMatrix& M);
std::shared_ptr<ListCpp> listcpp_from_rlist(const Rcpp::List& rlist);
