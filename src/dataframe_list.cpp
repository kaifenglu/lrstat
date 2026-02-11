#include "dataframe_list.h"

#include <Rcpp.h>

#include <algorithm> // copy_n, fill, max_element, remove
#include <cmath>     // isnan
#include <cstddef>   // size_t
#include <cstring>   // memcpy
#include <memory>    // make_shared, shared_ptr
#include <stdexcept> // invalid_argument, out_of_range, runtime_error
#include <string>    // string
#include <utility>   // move
#include <variant>   // get_if, variant
#include <vector>    // vector

// ------------------------- DataFrameCpp members (small) -------------------

void DataFrameCpp::push_back(const std::vector<double>& col,
                             const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  numeric_cols[name] = col;
  names_.push_back(name);
}
void DataFrameCpp::push_back(std::vector<double>&& col, const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  numeric_cols.emplace(name, std::move(col));
  names_.push_back(name);
}

void DataFrameCpp::push_back(const std::vector<int>& col, const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  int_cols[name] = col;
  names_.push_back(name);
}
void DataFrameCpp::push_back(std::vector<int>&& col, const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  int_cols.emplace(name, std::move(col));
  names_.push_back(name);
}

void DataFrameCpp::push_back(const std::vector<unsigned char>& col,
                             const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  bool_cols[name] = col;
  names_.push_back(name);
}
void DataFrameCpp::push_back(std::vector<unsigned char>&& col,
                             const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  bool_cols.emplace(name, std::move(col));
  names_.push_back(name);
}

void DataFrameCpp::push_back(const std::vector<std::string>& col,
                             const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  string_cols[name] = col;
  names_.push_back(name);
}
void DataFrameCpp::push_back(std::vector<std::string>&& col,
                             const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  string_cols.emplace(name, std::move(col));
  names_.push_back(name);
}

// Scalar expansions (push_back)
void DataFrameCpp::push_back(double value, const std::string& name) {
  std::size_t cur = nrows();
  if (cur == 0 && !names_.empty())
    throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<double> col(cur, value);
  push_back(std::move(col), name);
}
void DataFrameCpp::push_back(int value, const std::string& name) {
  std::size_t cur = nrows();
  if (cur == 0 && !names_.empty())
    throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<int> col(cur, value);
  push_back(std::move(col), name);
}
void DataFrameCpp::push_back(bool value, const std::string& name) {
  std::size_t cur = nrows();
  if (cur == 0 && !names_.empty())
    throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<unsigned char> col(cur, value);
  push_back(std::move(col), name);
}
void DataFrameCpp::push_back(const std::string& value, const std::string& name) {
  std::size_t cur = nrows();
  if (cur == 0 && !names_.empty())
    throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<std::string> col(cur, value);
  push_back(std::move(col), name);
}

// push_back_flat accepts a column-major flattened buffer containing nrows * p values
void DataFrameCpp::push_back_flat(const std::vector<double>& flat_col_major,
                                  int nrows, const std::string& base_name) {
  if (flat_col_major.empty()) return;
  if (containElementNamed(base_name))
    throw std::runtime_error("Column '" + base_name + "' already exists.");
  if (nrows == 0) throw std::runtime_error("nrows must be > 0");
  if (flat_col_major.size() % nrows != 0)
    throw std::runtime_error("flattened data size is not divisible by nrows");
  int p = flat_col_major.size() / nrows;
  check_row_size(nrows, base_name);

  if (p == 1) {
    std::vector<double> col(nrows);
    std::copy_n(flat_col_major.begin(), nrows, col.begin());
    numeric_cols[base_name] = std::move(col);
    names_.push_back(base_name);
  } else {
    for (int c = 0; c < p; ++c) {
      std::string col_name = base_name + "." + std::to_string(c + 1);
      if (containElementNamed(col_name))
        throw std::runtime_error("Column '" + col_name + "' already exists.");
      const double* fcol = flat_col_major.data() + c * nrows;
      std::vector<double> col(fcol, fcol + nrows);
      numeric_cols.emplace(col_name, std::move(col));
      names_.push_back(col_name);
    }
  }
}

void DataFrameCpp::push_back(const FlatMatrix& fm, const std::string& base_name) {
  push_back_flat(fm.data, fm.nrow, base_name);
}

void DataFrameCpp::push_back(FlatMatrix&& fm, const std::string& base_name) {
  push_back_flat(fm.data, fm.nrow, base_name);
}

// Push front variants (vectors)
void DataFrameCpp::push_front(const std::vector<double>& col,
                              const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  numeric_cols[name] = col;
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(std::vector<double>&& col, const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  numeric_cols.emplace(name, std::move(col));
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(const std::vector<int>& col,
                              const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  int_cols[name] = col;
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(std::vector<int>&& col, const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  int_cols.emplace(name, std::move(col));
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(const std::vector<unsigned char>& col,
                              const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  bool_cols[name] = col;
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(std::vector<unsigned char>&& col,
                              const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  bool_cols.emplace(name, std::move(col));
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(const std::vector<std::string>& col,
                              const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  string_cols[name] = col;
  names_.insert(names_.begin(), name);
}
void DataFrameCpp::push_front(std::vector<std::string>&& col,
                              const std::string& name) {
  if (containElementNamed(name))
    throw std::runtime_error("Column '" + name + "' already exists.");
  check_row_size(col.size(), name);
  string_cols.emplace(name, std::move(col));
  names_.insert(names_.begin(), name);
}

// Scalar expansions for push_front
void DataFrameCpp::push_front(double value, const std::string& name) {
  std::size_t cur = nrows();
  if (cur == 0 && !names_.empty())
    throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<double> col(cur, value);
  push_front(std::move(col), name);
}
void DataFrameCpp::push_front(int value, const std::string& name) {
  std::size_t cur = nrows();
  if (cur == 0 && !names_.empty())
    throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<int> col(cur, value);
  push_front(std::move(col), name);
}
void DataFrameCpp::push_front(bool value, const std::string& name) {
  std::size_t cur = nrows();
  if (cur == 0 && !names_.empty())
    throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<unsigned char> col(cur, value);
  push_front(std::move(col), name);
}
void DataFrameCpp::push_front(const std::string& value, const std::string& name) {
  std::size_t cur = nrows();
  if (cur == 0 && !names_.empty())
    throw std::runtime_error("Cannot push scalar when DataFrame has 0 rows");
  if (cur == 0) cur = 1;
  std::vector<std::string> col(cur, value);
  push_front(std::move(col), name);
}

void DataFrameCpp::erase(const std::string& name) {
  if (numeric_cols.erase(name)) { }
  else if (int_cols.erase(name)) { }
  else if (bool_cols.erase(name)) { }
  else string_cols.erase(name);
  names_.erase(std::remove(names_.begin(), names_.end(), name), names_.end());
}


// -------------------------- ListCpp small members --------------------------

void ListCpp::push_back(const ListCpp& l, const std::string& name) {
  push_back(std::make_shared<ListCpp>(l), name);
}
void ListCpp::push_back(ListCpp&& l, const std::string& name) {
  push_back(std::make_shared<ListCpp>(std::move(l)), name);
}
void ListCpp::push_back(const ListPtr& p, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, p);
  names_.push_back(name);
}
void ListCpp::push_back(ListPtr&& p, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(p));
  names_.push_back(name);
}

// FlatMatrix overloads
void ListCpp::push_back(const FlatMatrix& fm, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, fm);
  names_.push_back(name);
}
void ListCpp::push_back(FlatMatrix&& fm, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(fm));
  names_.push_back(name);
}
// IntMatrix overloads
void ListCpp::push_back(const IntMatrix& im, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, im);
  names_.push_back(name);
}
void ListCpp::push_back(IntMatrix&& im, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(im));
  names_.push_back(name);
}
// BoolMatrix overloads
void ListCpp::push_back(const BoolMatrix& im, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, im);
  names_.push_back(name);
}
void ListCpp::push_back(BoolMatrix&& im, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(im));
  names_.push_back(name);
}

// FlatArray overloads
void ListCpp::push_back(const FlatArray& fa, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, fa);
  names_.push_back(name);
}
void ListCpp::push_back(FlatArray&& fa, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(fa));
  names_.push_back(name);
}

// DataFrameCpp overloads
void ListCpp::push_back(const DataFrameCpp& df, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, df);
  names_.push_back(name);
}
void ListCpp::push_back(DataFrameCpp&& df, const std::string& name) {
  if (containsElementNamed(name))
    throw std::runtime_error("Element '" + name + "' already exists.");
  data.emplace(name, std::move(df));
  names_.push_back(name);
}

void ListCpp::erase(const std::string& name) {
  data.erase(name);
  names_.erase(std::remove(names_.begin(), names_.end(), name), names_.end());
}

ListCpp& ListCpp::get_list(const std::string& name) {
  if (!containsElementNamed(name))
    throw std::runtime_error("Element with name '" + name + "' not found.");
  auto& var = data.at(name);
  if (auto p = std::get_if<ListPtr>(&var)) {
    if (!*p) throw std::runtime_error(
        "List pointer is null for element '" + name + "'");
    return **p;
  }
  throw std::runtime_error("Element '" + name + "' is not a ListCpp");
}
const ListCpp& ListCpp::get_list(const std::string& name) const {
  if (!containsElementNamed(name))
    throw std::runtime_error("Element with name '" + name + "' not found.");
  const auto& var = data.at(name);
  if (auto p = std::get_if<ListPtr>(&var)) {
    if (!*p) throw std::runtime_error(
        "List pointer is null for element '" + name + "'");
    return **p;
  }
  throw std::runtime_error("Element '" + name + "' is not a ListCpp");
}

// Validate row indices (throws if any out-of-range)
static inline void validate_row_idx(int nrows, const std::vector<int>& row_idx) {
  int n = row_idx.size();
  for (int k = 0; k < n; ++k) {
    if (row_idx[k] >= nrows)
      throw std::out_of_range("row_idx contains an out-of-range index");
  }
}

// Non-mutating subset: return a new FlatMatrix with only the rows in row_idx
FlatMatrix subset_flatmatrix(const FlatMatrix& fm, const std::vector<int>& row_idx) {
  int old_nrow = fm.nrow;
  int ncol = fm.ncol;
  int new_nrow = row_idx.size();
  if (old_nrow == 0 || ncol == 0 || new_nrow == 0) {
    return FlatMatrix(new_nrow, ncol);
  }
  validate_row_idx(old_nrow, row_idx);

  FlatMatrix out(new_nrow, ncol);
  for (int c = 0; c < ncol; ++c) {
    const double* fmcol = fm.data.data() + c * old_nrow;
    double* outcol = out.data.data() + c * new_nrow;
    for (int r = 0; r < new_nrow; ++r) {
      outcol[r] = fmcol[row_idx[r]];
    }
  }
  return out;
}

// In-place subset: keep only rows specified by row_idx (in that order).
void subset_in_place_flatmatrix(FlatMatrix& fm, const std::vector<int>& row_idx) {
  int old_nrow = fm.nrow;
  int ncol = fm.ncol;
  int new_nrow = row_idx.size();
  if (old_nrow == 0 || ncol == 0) {
    // nothing to do: ensure matrix becomes empty if row_idx is empty
    fm.data.clear();
    fm.nrow = new_nrow;
    fm.ncol = ncol;
    return;
  }
  validate_row_idx(old_nrow, row_idx);

  if (new_nrow == old_nrow) {
    // check whether row_idx is identity; if yes, nothing to do
    bool identity = true;
    for (int i = 0; i < new_nrow; ++i) {
      if (row_idx[i] != i) { identity = false; break; } }
    if (identity) return;
  }

  std::vector<double> newdata;
  newdata.resize(new_nrow * ncol);

  for (int c = 0; c < ncol; ++c) {
    const double* fmcol = fm.data.data() + c * old_nrow;
    double* newcol = newdata.data() + c * new_nrow;
    for (int r = 0; r < new_nrow; ++r) {
      newcol[r] = fmcol[row_idx[r]];
    }
  }

  fm.data = std::move(newdata);
  fm.nrow = new_nrow;
  // fm.ncol remains unchanged
}

// Non-mutating subset: return a new FlatMatrix with only rows [start, end)
FlatMatrix subset_flatmatrix(const FlatMatrix& fm, int start, int end) {
  int old_nrow = fm.nrow;
  int ncol = fm.ncol;
  if (start < 0 || end <= start || end > old_nrow)
    throw std::invalid_argument("subset_flatmatrix: invalid start/end");

  int new_nrow = end - start;
  if (old_nrow == 0 || ncol == 0 || new_nrow == 0) {
    return FlatMatrix(new_nrow, ncol);
  }

  FlatMatrix out(new_nrow, ncol);
  for (int c = 0; c < ncol; ++c) {
    const double* src = fm.data.data() + c * old_nrow  + start;
    double* dst = out.data.data() + c * new_nrow;
    std::memcpy(dst, src, new_nrow * sizeof(double));
  }
  return out;
}

// In-place subset: keep only rows [start, end) (in that order).
void subset_in_place_flatmatrix(FlatMatrix& fm, int start, int end) {
  FlatMatrix tmp = subset_flatmatrix(fm, start, end);
  fm = std::move(tmp);
}

// Non-mutating subset: return a new FlatMatrix with rows/columns in specified ranges
FlatMatrix subset_flatmatrix(const FlatMatrix& fm, int row_start, int row_end,
                             int col_start, int col_end) {
  // Calculate dimensions of result
  int nrows = row_end - row_start;
  int ncols = col_end - col_start;

  // Validate ranges
  if (row_start < 0 || row_end > fm.nrow || row_start > row_end) {
    throw std::out_of_range("Invalid row range.");
  }
  if (col_start < 0 || col_end > fm.ncol || col_start > col_end) {
    throw std::out_of_range("Invalid column range.");
  }

  FlatMatrix out(nrows, ncols);
  // Copy column by column (efficient for column-major storage)
  for (int c = 0; c < ncols; ++c) {
    int src_col = col_start + c;
    // Source: column src_col, rows [row_start, row_end]
    const double* src = fm.data.data() + src_col * fm.nrow + row_start;
    // Destination: column c, all rows
    double* dst = out.data.data() + c * nrows;
    // Copy contiguous column segment
    std::memcpy(dst, src, nrows * sizeof(double));
  }
  return out;
}


// Non-mutating subset for FlatArray: return a new FlatArray with rows in row_idx
FlatArray subset_flatarray(const FlatArray& fa, const std::vector<int>& row_idx) {
  int old_nrow = fa.nrow;
  int ncol = fa.ncol;
  int nslice = fa.nslice;
  int new_nrow = row_idx.size();

  if (old_nrow == 0 || ncol == 0 || nslice == 0 || new_nrow == 0) {
    return FlatArray(new_nrow, ncol, nslice);
  }

  validate_row_idx(old_nrow, row_idx);

  FlatArray out(new_nrow, ncol, nslice);
  for (int s = 0; s < nslice; ++s) {
    int off_old_slice = s * (old_nrow * ncol);
    int off_new_slice = s * (new_nrow * ncol);
    for (int c = 0; c < ncol; ++c) {
      int off_old = off_old_slice + c * old_nrow;
      int off_new = off_new_slice + c * new_nrow;
      const double* facol = fa.data_ptr() + off_old;
      double* outcol = out.data_ptr() + off_new;
      for (int r = 0; r < new_nrow; ++r) {
        outcol[r] = facol[row_idx[r]];
      }
    }
  }
  return out;
}

// In-place subset for FlatArray: keep only rows specified by row_idx.
void subset_in_place_flatarray(FlatArray& fa, const std::vector<int>& row_idx) {
  int old_nrow = fa.nrow;
  int ncol = fa.ncol;
  int nslice = fa.nslice;
  int new_nrow = row_idx.size();

  if (old_nrow == 0 || ncol == 0 || nslice == 0) {
    // nothing to do: ensure array becomes empty if row_idx is empty
    fa.data.clear();
    fa.nrow = new_nrow;
    fa.ncol = ncol;
    fa.nslice = nslice;
    return;
  }

  validate_row_idx(old_nrow, row_idx);

  if (new_nrow == old_nrow) {
    bool identity = true;
    for (int i = 0; i < new_nrow; ++i) {
      if (row_idx[i] != i) { identity = false; break; }
    }
    if (identity) return;
  }

  std::vector<double> newdata;
  newdata.resize(static_cast<std::size_t>(new_nrow) * ncol * nslice);

  // layout: slice outermost, then column, then row
  for (int s = 0; s < nslice; ++s) {
    int off_old_slice = s * (old_nrow * ncol);
    int off_new_slice = s * (new_nrow * ncol);
    for (int c = 0; c < ncol; ++c) {
      int off_old = off_old_slice + c * old_nrow;
      int off_new = off_new_slice + c * new_nrow;
      const double* facol = fa.data_ptr() + off_old;
      double* newcol = newdata.data() + off_new;
      for (int r = 0; r < new_nrow; ++r) {
        newcol[r] = facol[row_idx[r]];
      }
    }
  }

  fa.data = std::move(newdata);
  fa.nrow = new_nrow;
  // fa.ncol and fa.nslice remain unchanged
}

// Row-wise concatenation: both matrices must have same number of columns.
FlatMatrix concat_flatmatrix(const FlatMatrix& fm1, const FlatMatrix& fm2) {
  // trivial cases
  if (fm1.ncol != fm2.ncol)
    throw std::invalid_argument("concat_flatmatrix: column counts do not match");
  if (fm1.nrow == 0) return fm2;
  if (fm2.nrow == 0) return fm1;

  int ncol = fm1.ncol;
  int nrow1 = fm1.nrow;
  int nrow2 = fm2.nrow;
  int nrow = nrow1 + nrow2;

  FlatMatrix out(nrow, ncol);
  // For each column, copy rows from fm1 then fm2 into the destination column buffer
  for (int c = 0; c < ncol; ++c) {
    const double* src1 = fm1.data.empty() ? nullptr : fm1.data.data() + c * nrow1;
    const double* src2 = fm2.data.empty() ? nullptr : fm2.data.data() + c * nrow2;
    double* dst = out.data.data() + c * nrow;

    if (src1 && nrow1 > 0) {
      std::memcpy(dst, src1, nrow1 * sizeof(double));
    } else {
      if (nrow1 > 0) std::fill(dst, dst + nrow1, 0.0);
    }
    if (src2 && nrow2 > 0) {
      std::memcpy(dst + nrow1, src2, nrow2 * sizeof(double));
    } else {
      if (nrow2 > 0) std::fill(dst + nrow1, dst + nrow, 0.0);
    }
  }
  return out;
}

// Row-wise append: both matrices must have same number of columns.
void append_flatmatrix(FlatMatrix& fm1, const FlatMatrix& fm2) {
  FlatMatrix tmp = concat_flatmatrix(fm1, fm2);
  fm1 = std::move(tmp);
}


void append_flatmatrix(std::vector<std::vector<double>>& fm1,
                       const FlatMatrix& fm2) {
  int ncol = fm2.ncol;
  int nrow2 = fm2.nrow;
  if (fm1.size() != static_cast<std::size_t>(ncol))
    throw std::invalid_argument(
        "append_flatmatrix: fm1 column count does not match fm2");

  // Check that all columns in fm1 have the same number of rows & capture that count
  std::size_t nrow1 = 0;
  if (!fm1.empty()) {
    nrow1 = fm1[0].size();
    for (std::size_t c = 1; c < fm1.size(); ++c) {
      if (fm1[c].size() != nrow1) {
        throw std::invalid_argument(
            "append_flatmatrix: fm1 columns have inconsistent row counts");
      }
    }
  }

  if (nrow2 == 0) return; // nothing to append

  for (std::size_t c = 0; c < fm1.size(); ++c) {
    const double* src = fm2.data_ptr() + c * nrow2;
    auto &dest_col = fm1[c];
    dest_col.resize(nrow1 + nrow2);
    std::memcpy(dest_col.data() + nrow1, src, nrow2 * sizeof(double));
  }
}


void append_flatmatrix(std::vector<std::vector<double>>& fm1,
                       const std::vector<std::vector<double>>& fm2) {
  if (fm2.empty()) return;

  // If fm1 is empty, just copy fm2 (fast path)
  if (fm1.empty()) {
    fm1 = fm2;
    return;
  }

  // Require same number of columns
  if (fm1.size() != fm2.size()) {
    throw std::invalid_argument("append_flatmatrix: number of columns must match");
  }

  std::size_t ncol = fm1.size();
  std::size_t nrow1 = fm1[0].size();
  std::size_t nrow2 = fm2[0].size();

  for (std::size_t c = 0; c < ncol; ++c) {
    // Validate column sizes
    if (fm1[c].size() != nrow1) {
      throw std::invalid_argument(
          "append_flatmatrix: inconsistent row counts in fm1");
    }
    if (fm2[c].size() != nrow2) {
      throw std::invalid_argument(
          "append_flatmatrix: inconsistent row counts in fm2");
    }
    // Append column c
    auto &dest_col = fm1[c];
    std::size_t old_size = dest_col.size();
    dest_col.resize(old_size + nrow2);
    std::memcpy(dest_col.data() + old_size, fm2[c].data(), nrow2 * sizeof(double));
  }
}


/*
 Convert a vector-of-columns into a FlatMatrix (column-major).
 - 'cols' outer vector: columns
 - each inner vector: contiguous column values (length = nrow)
 Throws std::invalid_argument if columns have inconsistent lengths.
 */
FlatMatrix cols_to_flatmatrix(const std::vector<std::vector<double>>& cols) {
  const std::size_t ncol = cols.size();
  if (ncol == 0) {
    return FlatMatrix(); // empty
  }

  const std::size_t nrow = cols[0].size();
  // validate sizes
  for (std::size_t c = 1; c < ncol; ++c) {
    if (cols[c].size() != nrow) {
      throw std::invalid_argument("cols_to_flatmatrix: inconsistent column lengths");
    }
  }

  // allocate contiguous column-major buffer: nrow * ncol
  std::vector<double> data;
  data.resize(nrow * ncol);

  // copy each column into its place: column c starts at data + c * nrow
  const std::size_t bytes_per_col = nrow * sizeof(double);
  for (std::size_t c = 0; c < ncol; ++c) {
    if (nrow > 0) {
      std::memcpy(data.data() + c * nrow, cols[c].data(), bytes_per_col);
    }
  }

  // construct FlatMatrix by moving the data buffer (avoids another copy)
  return FlatMatrix(std::move(data), static_cast<int>(nrow), static_cast<int>(ncol));
}


std::vector<double> flatmatrix_get_column(const FlatMatrix& M, int col) {
  if (M.nrow == 0 || M.ncol == 0) return {};
  if (col < 0 || col >= M.ncol) throw std::out_of_range("column index out of range");
  const double* src = M.data_ptr() + FlatMatrix::idx_col(0, col, M.nrow);
  std::vector<double> out(static_cast<std::size_t>(M.nrow));
  std::memcpy(out.data(), src, static_cast<std::size_t>(M.nrow) * sizeof(double));
  return out;
}

std::vector<int> intmatrix_get_column(const IntMatrix& M, int col) {
  if (M.nrow == 0 || M.ncol == 0) return {};
  if (col < 0 || col >= M.ncol) throw std::out_of_range("column index out of range");
  const int* src = M.data_ptr() + IntMatrix::idx_col(0, col, M.nrow);
  std::vector<int> out(static_cast<std::size_t>(M.nrow));
  std::memcpy(out.data(), src, static_cast<std::size_t>(M.nrow) * sizeof(int));
  return out;
}

std::vector<unsigned char> boolmatrix_get_column(const BoolMatrix& M, int col) {
  if (M.nrow == 0 || M.ncol == 0) return {};
  if (col < 0 || col >= M.ncol) throw std::out_of_range("column index out of range");
  const unsigned char* src = M.data_ptr() + BoolMatrix::idx_col(0, col, M.nrow);
  std::vector<unsigned char> out(static_cast<std::size_t>(M.nrow));
  std::memcpy(out.data(), src, 
              static_cast<std::size_t>(M.nrow) * sizeof(unsigned char));
  return out;
}

void flatmatrix_set_column(FlatMatrix& M, int col, const std::vector<double>& src) {
  if (col < 0 || col >= M.ncol) throw std::out_of_range("col out of range");
  if (static_cast<int>(src.size()) != M.nrow)
    throw std::invalid_argument("src size != M.nrow");
  const double* src_ptr = src.data();
  double* dst_ptr = M.data_ptr() + FlatMatrix::idx_col(0, col, M.nrow);
  std::memcpy(dst_ptr, src_ptr, static_cast<std::size_t>(M.nrow) * sizeof(double));
}

void intmatrix_set_column(IntMatrix& M, int col, const std::vector<int>& src) {
  if (col < 0 || col >= M.ncol) throw std::out_of_range("col out of range");
  if (static_cast<int>(src.size()) != M.nrow)
    throw std::invalid_argument("src size != M.nrow");
  const int* src_ptr = src.data();
  int* dst_ptr = M.data_ptr() + IntMatrix::idx_col(0, col, M.nrow);
  std::memcpy(dst_ptr, src_ptr, static_cast<std::size_t>(M.nrow) * sizeof(int));
}

void boolmatrix_set_column(BoolMatrix& M, int col, 
                           const std::vector<unsigned char>& src) {
  if (col < 0 || col >= M.ncol) throw std::out_of_range("col out of range");
  if (static_cast<int>(src.size()) != M.nrow)
    throw std::invalid_argument("src size != M.nrow");
  const unsigned char* src_ptr = src.data();
  unsigned char* dst_ptr = M.data_ptr() + BoolMatrix::idx_col(0, col, M.nrow);
  std::memcpy(dst_ptr, src_ptr, 
              static_cast<std::size_t>(M.nrow) * sizeof(unsigned char));
}

// -------------------------- Converters implementations ---------------------
// Convert DataFrameCpp to Rcpp::DataFrame
Rcpp::DataFrame convertDataFrameCppToR(const DataFrameCpp& df) {
  Rcpp::List cols;
  Rcpp::CharacterVector names;
  for (const auto& nm : df.names_) {
    names.push_back(nm);
    if (df.numeric_cols.count(nm)) {
      cols.push_back(Rcpp::wrap(df.numeric_cols.at(nm)));
    } else if (df.int_cols.count(nm)) {
      cols.push_back(Rcpp::wrap(df.int_cols.at(nm)));
    } else if (df.bool_cols.count(nm)) {
      std::vector<unsigned char> v = df.bool_cols.at(nm);
      R_xlen_t n = v.size();
      Rcpp::LogicalVector lv(n);
      for (R_xlen_t i = 0; i < n; ++i) {
        unsigned char bv = v[i];
        if (bv == 255) lv[i] = NA_LOGICAL;
        else lv[i] = (bv != 0);
      }
      cols.push_back(std::move(lv));
    } else if (df.string_cols.count(nm)) {
      cols.push_back(Rcpp::wrap(df.string_cols.at(nm)));
    } else {
      cols.push_back(R_NilValue);
    }
  }
  cols.names() = names;
  return Rcpp::as<Rcpp::DataFrame>(cols);
}

// Visitor
struct RcppVisitor {
  template <typename T>
  Rcpp::RObject operator()(const T& x) const {
    return Rcpp::wrap(x);
  }

  Rcpp::RObject operator()(const FlatMatrix& fm) const {
    if (fm.nrow == 0 || fm.ncol == 0) return R_NilValue;
    Rcpp::NumericMatrix M(fm.nrow, fm.ncol);
    std::memcpy(REAL(M), fm.data.data(), fm.data.size() * sizeof(double));
    return M;
  }

  Rcpp::RObject operator()(const IntMatrix& im) const {
    if (im.nrow == 0 || im.ncol == 0) return R_NilValue;
    Rcpp::IntegerMatrix M(im.nrow, im.ncol);
    std::memcpy(INTEGER(M), im.data.data(), im.data.size() * sizeof(int));
    return M;
  }

  Rcpp::RObject operator()(const FlatArray& fa) const {
    if (fa.nrow == 0 || fa.ncol == 0 || fa.nslice == 0) return R_NilValue;
    Rcpp::NumericVector vec(static_cast<R_xlen_t>(fa.data.size()));
    std::memcpy(REAL(vec), fa.data.data(), fa.data.size() * sizeof(double));
    Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(fa.nrow, fa.ncol,
                                                           fa.nslice);
    vec.attr("dim") = dims;
    return vec;
  }

  Rcpp::RObject operator()(const DataFrameCpp& df) const {
    return convertDataFrameCppToR(df);
  }

  Rcpp::RObject operator()(const ListPtr& p) const {
    if (!p) return R_NilValue;
    return convertListCppToR(*p);
  }

  Rcpp::RObject operator()(const std::vector<DataFrameCpp>& dfs) const {
    Rcpp::List out;
    for (const auto& df : dfs) out.push_back(convertDataFrameCppToR(df));
    return out;
  }

  Rcpp::RObject operator()(const std::vector<ListPtr>& lists) const {
    Rcpp::List out;
    for (const auto& p : lists) {
      if (!p) out.push_back(R_NilValue);
      else out.push_back(convertListCppToR(*p));
    }
    return out;
  }
};

// convert ListCpp to Rcpp::List
Rcpp::List convertListCppToR(const ListCpp& L) {
  int p = L.names_.size();
  Rcpp::List out(p);
  Rcpp::CharacterVector names(p);
  RcppVisitor visitor;
  for (int i = 0; i < p; ++i) {
    const auto& nm = L.names_[i];
    const auto& var = L.data.at(nm);
    Rcpp::RObject obj = std::visit(visitor, var);
    out[i] = obj;
    names[i] = nm;
  }
  out.names() = names;
  return out;
}

// Convert Rcpp::NumericMatrix to FlatMatrix
FlatMatrix flatmatrix_from_Rmatrix(const Rcpp::NumericMatrix& M) {
  int nr = M.nrow(), nc = M.ncol();
  if (nr == 0 || nc == 0) return FlatMatrix();
  FlatMatrix fm(nr, nc);
  std::memcpy(fm.data.data(), REAL(M), nr * nc * sizeof(double));
  return fm;
}

// Convert Rcpp::IntegerMatrix to IntMatrix
IntMatrix intmatrix_from_Rmatrix(const Rcpp::IntegerMatrix& M) {
  int nr = M.nrow(), nc = M.ncol();
  if (nr == 0 || nc == 0) return IntMatrix();
  IntMatrix im(nr, nc);
  std::memcpy(im.data.data(), INTEGER(M), nr * nc * sizeof(int));
  return im;
}

// Convert Rcpp::LogicalMatrix to BoolMatrix
BoolMatrix boolmatrix_from_Rmatrix(const Rcpp::LogicalMatrix& M) {
  int nr = M.nrow(), nc = M.ncol();
  if (nr == 0 || nc == 0) return BoolMatrix();
  BoolMatrix im(nr, nc);
  const int* src = LOGICAL(M);
  unsigned char* dst = im.data.data();
  for (int i = 0, n = nr * nc; i < n; ++i) {
    if (src[i] == NA_LOGICAL) {
      dst[i] = 255;
    } else {
      dst[i] = static_cast<unsigned char>(src[i] != 0);
    }
  }
  return im;
}

// Convert Rcpp::DataFrame to DataFrameCpp
DataFrameCpp convertRDataFrameToCpp(const Rcpp::DataFrame& r_df) {
  DataFrameCpp df;
  Rcpp::CharacterVector cn = r_df.names();
  R_xlen_t nc = r_df.size();
  if (nc == 0) return df;

  for (R_xlen_t j = 0; j < nc; ++j) {
    std::string name;
    if (cn.size() > 0 && static_cast<R_xlen_t>(cn.size()) > j &&
        !Rcpp::StringVector::is_na(cn[j]) && std::string(cn[j]) != "") {
        name = Rcpp::as<std::string>(cn[j]);
    } else {
      name = "V" + std::to_string(j + 1);
    }

    Rcpp::RObject col = r_df[j];
    if (Rcpp::is<Rcpp::NumericVector>(col)) {
      Rcpp::NumericVector nv = Rcpp::as<Rcpp::NumericVector>(col); // wrapper
      R_xlen_t n = nv.size();
      std::vector<double> out;
      out.resize(static_cast<std::size_t>(n));
      if (n > 0) {
        // memcpy from R contiguous memory
        std::memcpy(out.data(), REAL(nv),
                    static_cast<std::size_t>(n) * sizeof(double));
      }
      df.push_back(std::move(out), name);

    } else if (Rcpp::is<Rcpp::IntegerVector>(col)) {
      Rcpp::IntegerVector iv = Rcpp::as<Rcpp::IntegerVector>(col);
      R_xlen_t n = iv.size();
      std::vector<int> out;
      out.resize(static_cast<std::size_t>(n));
      if (n > 0) {
        std::memcpy(out.data(), INTEGER(iv),
                    static_cast<std::size_t>(n) * sizeof(int));
      }
      df.push_back(std::move(out), name);
    } else if (Rcpp::is<Rcpp::LogicalVector>(col)) {
      Rcpp::LogicalVector lv = Rcpp::as<Rcpp::LogicalVector>(col);
      R_xlen_t n = lv.size();
      std::vector<unsigned char> out;
      out.resize(static_cast<std::size_t>(n));
      for (R_xlen_t i = 0; i < n; ++i) {
        int v = lv[i]; // 0,1 or NA_LOGICAL
        if (v == NA_LOGICAL) out[static_cast<std::size_t>(i)] = 255;
        else out[static_cast<std::size_t>(i)] = v ? 1 : 0;
      }
      df.push_back(std::move(out), name);
    } else if (Rcpp::is<Rcpp::CharacterVector>(col)) {
      std::vector<std::string> out = Rcpp::as<std::vector<std::string>>(col);
      df.push_back(std::move(out), name);
    } else {
      Rcpp::warning("Unsupported column type in DataFrame conversion: " + name);
    }
  }
  return df;
}

static bool convert_r_object_to_listcpp_variant(
    ListCpp& target, const std::string& name, const Rcpp::RObject& obj) {
  if (Rcpp::is<Rcpp::NumericMatrix>(obj)) {
    Rcpp::NumericMatrix M = Rcpp::as<Rcpp::NumericMatrix>(obj);
    FlatMatrix fm = flatmatrix_from_Rmatrix(M);
    target.push_back(std::move(fm), name);
    return true;
  }
  if (Rcpp::is<Rcpp::IntegerMatrix>(obj)) {
    Rcpp::IntegerMatrix M = Rcpp::as<Rcpp::IntegerMatrix>(obj);
    IntMatrix im = intmatrix_from_Rmatrix(M);
    target.push_back(std::move(im), name);
    return true;
  }
  if (Rcpp::is<Rcpp::DataFrame>(obj)) {
    Rcpp::DataFrame rdf = Rcpp::as<Rcpp::DataFrame>(obj);
    DataFrameCpp cppdf = convertRDataFrameToCpp(rdf);
    target.push_back(std::move(cppdf), name);
    return true;
  }
  if (Rcpp::is<Rcpp::List>(obj)) {
    Rcpp::List rl = Rcpp::as<Rcpp::List>(obj);
    std::shared_ptr<ListCpp> nested = listcpp_from_rlist(rl);
    target.push_back(std::move(nested), name);
    return true;
  }
  if (Rcpp::is<Rcpp::NumericVector>(obj)) {
    Rcpp::NumericVector nv = Rcpp::as<Rcpp::NumericVector>(obj);
    // Check for 'dim' attribute
    Rcpp::IntegerVector dim = nv.attr("dim");
    if (dim.size() == 3) {
      int nr = dim[0];
      int nc = dim[1];
      int ns = dim[2];
      std::vector<double> v = Rcpp::as<std::vector<double>>(nv);
      FlatArray fa(std::move(v), nr, nc, ns);
      target.push_back(std::move(fa), name);
      return true;
    }
    // fall through to treat as plain numeric vector (1D)
    {
      std::vector<double> v = Rcpp::as<std::vector<double>>(nv);
      target.push_back(std::move(v), name);
      return true;
    }
  }
  if (Rcpp::is<Rcpp::IntegerVector>(obj)) {
    Rcpp::IntegerVector iv = Rcpp::as<Rcpp::IntegerVector>(obj);
    std::vector<int> v = Rcpp::as<std::vector<int>>(iv);
    target.push_back(std::move(v), name);
    return true;
  }
  if (Rcpp::is<Rcpp::LogicalVector>(obj)) {
    Rcpp::LogicalVector lv = Rcpp::as<Rcpp::LogicalVector>(obj);
    R_xlen_t n = lv.size();
    std::vector<unsigned char> v;
    v.reserve(n);
    for (R_xlen_t i = 0; i < n; ++i) {
      int xi = lv[i];                       // returns int: 0,1 or NA_LOGICAL
      if (xi == NA_LOGICAL) v.push_back(255);
      else v.push_back(xi ? 1 : 0);
    }
    target.push_back(std::move(v), name);
    return true;
  }
  if (Rcpp::is<Rcpp::CharacterVector>(obj)) {
    Rcpp::CharacterVector cv = Rcpp::as<Rcpp::CharacterVector>(obj);
    std::vector<std::string> v = Rcpp::as<std::vector<std::string>>(cv);
    target.push_back(std::move(v), name);
    return true;
  }
  return false;
}

// Convert Rcpp::List to ListCpp
std::shared_ptr<ListCpp> listcpp_from_rlist(const Rcpp::List& rlist) {
  auto out = std::make_shared<ListCpp>();
  Rcpp::CharacterVector rn = rlist.names();
  R_xlen_t n = rlist.size();
  for (R_xlen_t i = 0; i < n; ++i) {
    std::string name;
    if (rn.size() > 0 && static_cast<R_xlen_t>(rn.size()) > i &&
        !Rcpp::StringVector::is_na(rn[i]) && std::string(rn[i]) != "")
        name = Rcpp::as<std::string>(rn[i]);
    else
      name = "V" + std::to_string(i + 1);
    Rcpp::RObject obj = rlist[i];
    bool ok = convert_r_object_to_listcpp_variant(*out, name, obj);
    if (!ok) Rcpp::warning(
        "listcpp_from_rlist: skipping unsupported element type for name '" +
          name + "'");
  }
  return out;
}
