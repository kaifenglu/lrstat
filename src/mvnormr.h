#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "utilities.h"
#include "dataframe_list.h"


// Prime bases. This list length sets the max supported dimension.
inline constexpr unsigned int primes[] = {
  2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
  103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
  199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
  313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,
  433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,
  563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,
  673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
  811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,
  941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,
  1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,
  1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249,1259,1277,
  1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,
  1399,1409,1423,1427,1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,
  1489,1493,1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,
  1601,1607,1609,1613,1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,
  1709,1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,
  1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,1949,
  1951,1973,1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,2063,
  2069,2081,2083,2087,2089,2099,2111,2113,2129,2131,2137,2141,2143,2153,2161,
  2179,2203,2207,2213,2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,2293,
  2297,2309,2311,2333,2339,2341,2347,2351,2357,2371,2377,2381,2383,2389,2393,
  2399,2411,2417,2423
};

inline constexpr unsigned int permTN2[] = {
  1,1,3,3,4,9,7,5,9,18,18,8,13,31,9,19,36,33,21,44,43,61,60,56,26,71,32,77,26,95,
  92,47,29,61,57,69,115,63,92,31,104,126,50,80,55,152,114,80,83,97,95,150,148,55,
  80,192,71,76,82,109,105,173,58,143,56,177,203,239,196,143,278,227,87,274,264,84,
  226,163,231,177,95,116,165,131,156,105,188,142,105,125,269,292,215,182,294,152,
  148,144,382,194,346,323,220,174,133,324,215,246,159,337,254,423,484,239,440,362,
  464,376,398,174,149,418,306,282,434,196,458,313,512,450,161,315,441,549,555,431,
  295,557,172,343,472,604,297,524,251,514,385,531,663,674,255,519,324,391,394,533,
  253,717,651,399,596,676,425,261,404,691,604,274,627,777,269,217,599,447,581,640,
  666,595,669,686,305,460,599,335,258,649,771,619,666,669,707,737,854,925,818,424,
  493,463,535,782,476,451,520,886,340,793,390,381,274,500,581,345,363,1024,514,
  773,932,556,954,793,294,863,393,827,527,1007,622,549,613,799,408,856,601,1072,
  938,322,1142,873,629,1071,1063,1205,596,973,984,875,918,1133,1223,933,1110,1228,
  1017,701,480,678,1172,689,1138,1022,682,613,635,984,526,1311,459,1348,477,716,
  1075,682,1245,401,774,1026,499,1314,743,693,1282,1003,1181,1079,765,815,1350,
  1144,1449,718,805,1203,1173,737,562,579,701,1104,1105,1379,827,1256,759,540,
  1284,1188,776,853,1140,445,1265,802,932,632,1504,856,1229,1619,774,1229,1300,
  1563,1551,1265,905,1333,493,913,1397,1250,612,1251,1765,1303,595,981,671,1403,
  820,1404,1661,973,1340,1015,1649,855,1834,1621,1704,893,1033,721,1737,1507,1851,
  1006,994,923,872,1860
};

constexpr std::size_t ghaltonMaxDim = 360;

inline std::uint64_t seed_from_clock_u64() {
  using namespace std::chrono;
  return static_cast<std::uint64_t>(
    high_resolution_clock::now().time_since_epoch().count());
}

inline std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// ---- Generalized Halton scramble (contiguous shift + inv_base) ----
struct GHaltonScramble {
  std::size_t J = 0;
  std::vector<unsigned int> base;   // length J-1 is used, store J for simplicity
  std::vector<unsigned int> f;      // length J
  std::vector<double> inv_base;     // length J
  std::vector<unsigned int> shift;  // length J*32, shift[j*32 + k]
};


inline unsigned int scrambled_digit(const GHaltonScramble& S,
                                    std::size_t j,
                                    std::size_t k,
                                    unsigned int digit) {
  return (S.f[j] * digit + S.shift[j * 32 + k]) % S.base[j];
}

// ---- Incremental-digit generalized Halton generator ----
// Only generates the first dim = J-1 coordinates
struct GHaltonIncremental {
  std::size_t dim = 0; // == J-1
  std::size_t i = 0;   // current index i (digits/u correspond to this i)

  // digits[j][k] in base base[j], k=0 least significant
  std::vector<std::array<unsigned int, 32>> digits; // length dim

  // current radical inverse value for each dimension
  std::vector<double> u; // length dim

  // invpow[j*32 + k] = 1 / base[j]^(k+1)
  std::vector<double> invpow; // length dim*32
};


struct ReplicateAccumulator {
  GHaltonScramble S;

  // incremental generator state for first J-1 dims
  GHaltonIncremental G;
  bool gen_initialized = false;

  std::size_t n_done = 0;
  double sum = 0.0;
};


struct PMVNResult {
  double prob;
  std::string method;
  double error;
  std::size_t nsamples;
};

// ---- Cached Cholesky factor for a fixed covariance matrix ----------------//
//                                                                            //
// Call precompute_chol() once per distinct sigma (O(J^3) work).  Then call  //
// pmvnorm_with_chol() for each new (lower, upper, mean) triple, replacing   //
// the O(J^3) factorisation inside pmvnormcpp with O(J) bound               //
// standardisation.  This is the key primitive for Rec 1 and Rec 2.         //
struct PMVNCholFactor {
  std::size_t J = 0;
  std::vector<double> sd; // sqrt of each diagonal entry of D in L D L^T
  std::vector<double> C;  // lower-triangular factor; packed, length J*(J-1)/2
};

// Returns true when sigma has compound symmetry with non-negative
// off-diagonal entries.  pmvnormcpp uses a fast 1-D integration path for
// such matrices; call pmvnormcpp rather than pmvnorm_with_chol in that case.
bool is_compound_symmetry(const FlatMatrix& sigma);

// Factorise sigma via Cholesky/LDL^T (O(J^3)) and return the result.
// Call this once per unique sigma; reuse the result across multiple bound
// evaluations that share the same sigma.
PMVNCholFactor precompute_chol(const FlatMatrix& sigma);

// Evaluate Pr(lower < X < upper) using a cached Cholesky factor.
// Only O(J) work per call (bound standardisation) rather than O(J^3).
// Caller is responsible for not passing a compound-symmetric sigma here;
// use pmvnormcpp for such matrices to benefit from its 1-D integration path.
PMVNResult pmvnorm_with_chol(const PMVNCholFactor& chol,
                              const std::vector<double>& lower,
                              const std::vector<double>& upper,
                              const std::vector<double>& mean,
                              std::size_t   n0      = 1024,
                              std::size_t   n_max   = 16384,
                              std::size_t   R       = 8,
                              double        abseps  = 1e-4,
                              double        releps  = 0.0,
                              std::uint64_t seed    = 314159,
                              bool          parallel = true);


PMVNResult pmvnormcpp(const std::vector<double>& lower,
                      const std::vector<double>& upper,
                      const std::vector<double>& mean,
                      const FlatMatrix& sigma,
                      std::size_t n0 = 1024,
                      std::size_t n_max = 16384,
                      std::size_t R = 8,
                      double abseps = 1e-4,
                      double releps = 0.0,
                      uint64_t seed = 0,
                      bool parallel = true);

double qmvnormcpp(const double p,
                  const std::vector<double>& mean,
                  const FlatMatrix& sigma,
                  std::size_t n0 = 1024,
                  std::size_t n_max = 16384,
                  std::size_t R = 8,
                  double abseps = 1e-4,
                  double releps = 0.0,
                  uint64_t seed = 0,
                  bool parallel = true);
