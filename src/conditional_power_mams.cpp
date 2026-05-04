#include "generic_design.h"
#include "mams_design.h"
#include "utilities.h"
#include "dataframe_list.h"

#include <algorithm>     // any_of, fill
#include <cctype>        // tolower
#include <cmath>         // isnan
#include <cstring>       // memcpy
#include <numeric>       // accumulate
#include <stdexcept>     // invalid_argument
#include <string>        // string
#include <vector>        // vector

#include <Rcpp.h>

using std::size_t;


// Helper function to compute conditional power
double getCP_mams_cpp(
    const double INew,
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t L,
    const std::vector<double>& zL,
    const std::vector<double>& theta,
    const double IMax,
    const size_t kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const bool MullerSchafer,
    const size_t MNew,
    const std::vector<int>& selected,
    const double rNew,
    const size_t kNew,
    const std::vector<double>& informationRatesNew,
    const std::vector<unsigned char>& efficacyStoppingNew,
    const std::string& typeAlphaSpendingNew,
    const double parameterAlphaSpendingNew,
    const std::vector<double>& spendingTimeNew) {

  // Basic validations
  if (std::isnan(INew)) throw std::invalid_argument("INew must be provided");
  if (INew <= 0.0) throw std::invalid_argument("INew must be positive");

  if (M < 1) throw std::invalid_argument("M must be at least 1");
  if (r <= 0.0) throw std::invalid_argument("r must be positive");
  if (L < 1) throw std::invalid_argument("L must be at least 1");
  if (!none_na(zL)) throw std::invalid_argument("zL must be provided");
  if (zL.size() != M) throw std::invalid_argument("Invalid length for zL");
  if (!none_na(theta)) throw std::invalid_argument("theta must be provided");
  if (theta.size() != M) throw std::invalid_argument("Invalid length for theta");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (kMax <= L) throw std::invalid_argument("kMax must be greater than L");

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(kMax);
  if (none_na(informationRates)) {
    if (informationRates.size() != kMax)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates.back() != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else {
    for (size_t i = 0; i < kMax; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(kMax);
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != kMax)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping.back() != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(kMax, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);

  if (!missingCriticalValues && criticalValues.size() != kMax) {
    throw std::invalid_argument("Invalid length for criticalValues");
  }
  if (missingCriticalValues && std::isnan(alpha)) {
    throw std::invalid_argument("alpha must be provided for missing criticalValues");
  }

  std::string asf = typeAlphaSpending;
  for (char &c : asf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingCriticalValues && !(asf == "of" || asf == "p" ||
      asf == "wt" || asf == "sfof" || asf == "sfp" ||
      asf == "sfkd" || asf == "sfhsd" || asf == "user" || asf == "none")) {
    throw std::invalid_argument("Invalid value for typeAlphaSpending");
  }
  if ((asf == "wt" || asf == "sfkd" || asf == "sfhsd") &&
      std::isnan(parameterAlphaSpending)) {
    throw std::invalid_argument("Missing value for parameterAlphaSpending");
  }
  if (asf == "sfkd" && parameterAlphaSpending <= 0.0) {
    throw std::invalid_argument ("parameterAlphaSpending must be positive for sfKD");
  }

  if (missingCriticalValues && asf == "user") {
    if (!none_na(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be specified");
    if (userAlphaSpending.size() != kMax)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending.back() != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != kMax)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime.back() != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // ----------- New Design Input Validation ----------- //
  std::vector<size_t> selectedNew;
  std::vector<double> infoRatesNew;
  std::vector<unsigned char> effStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (MNew < 1) throw std::invalid_argument("MNew must be at least 1");
  for (auto i : selected) {
    if (i < 1 || i > static_cast<int>(M)) {
      throw std::invalid_argument("Invalid value in selected");
    }
    selectedNew.push_back(static_cast<size_t>(i - 1));
  }
  selectedNew = unique_sorted(selectedNew);
  if (selectedNew.size() != MNew)
    throw std::invalid_argument("Length of selected does not match MNew");

  if (rNew <= 0.0) throw std::invalid_argument("rNew must be positive");

  if (MullerSchafer) {
    if (kNew < 1) throw std::invalid_argument("kNew must be at least 1");

    // informationRatesNew: default to (1:kNew)/kNew if missing
    infoRatesNew.resize(kNew);
    if (none_na(informationRatesNew)) {
      if (informationRatesNew.size() != kNew)
        throw std::invalid_argument("Invalid length for informationRatesNew");
      if (informationRatesNew[0] <= 0.0)
        throw std::invalid_argument("informationRatesNew must be positive");
      if (any_nonincreasing(informationRatesNew))
        throw std::invalid_argument("informationRatesNew must be increasing");
      if (informationRatesNew[kNew-1] != 1.0)
        throw std::invalid_argument("informationRatesNew must end with 1");
      infoRatesNew = informationRatesNew; // copy
    } else {
      for (size_t i = 0; i < kNew; ++i)
        infoRatesNew[i] = static_cast<double>(i+1) / static_cast<double>(kNew);
    }

    // effStoppingNew: default to all 1s if missing
    if (none_na(efficacyStoppingNew)) {
      if (efficacyStoppingNew.size() != kNew)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew[kNew-1] != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew; // copy
    } else {
      effStoppingNew.assign(kNew, 1);
    }

    if (!(asfNew == "of" || asfNew == "sfof" || asfNew == "sfp" ||
        asfNew == "sfkd" || asfNew == "sfhsd" || asfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpendingNew");
    }
    if ((asfNew == "sfkd" || asfNew == "sfhsd") &&
        std::isnan(parameterAlphaSpendingNew)) {
      throw std::invalid_argument("Missing value for parameterAlphaSpendingNew");
    }
    if (asfNew == "sfkd" && parameterAlphaSpendingNew <= 0.0) {
      throw std::invalid_argument (
          "parameterAlphaSpendingNew must be positive for sfKD");
    }

    // spendingTimeNew: default to informationRatesNew if missing
    if (none_na(spendingTimeNew)) {
      if (spendingTimeNew.size() != kNew)
        throw std::invalid_argument("Invalid length for spendingTimeNew");
      if (spendingTimeNew[0] <= 0.0)
        throw std::invalid_argument("spendingTimeNew must be positive");
      if (any_nonincreasing(spendingTimeNew))
        throw std::invalid_argument("spendingTimeNew must be increasing");
      if (spendingTimeNew[kNew-1] != 1.0)
        throw std::invalid_argument("spendingTimeNew must end with 1");
    } else {
      spendTimeNew = infoRatesNew;
    }
  }

  // ----------- End of Input Validation ----------- //

  // obtain critical values for the primary trial
  std::vector<double> zero(M, 0.0);
  std::vector<double> critValues = criticalValues;
  bool haybittle = false;
  if (missingCriticalValues) {
    if (kMax > 1 && criticalValues.size() == kMax) {
      bool hasNaN = false;
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[kMax-1])) haybittle = true;
    }

    if (haybittle) { // Haybittle & Peto
      FlatMatrix b(M, kMax);
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (!effStopping[i]) critValues[i] = 8.0;
        double* colptr = b.data_ptr() + i * M;
        std::fill_n(colptr, M, critValues[i]);
      }

      double* last_col = b.data_ptr() + (kMax - 1) * M;
      auto f = [&](double x)->double {
        std::fill_n(last_col, M, x);
        auto v = exitprob_mams_cpp(M, r, zero, corr_known, kMax, b, infoRates);
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, 0.0, 8.0, 1e-6);
    } else {
      critValues = getBound_mams_cpp(
        M, r, corr_known, kMax, infoRates, alpha, asf, parameterAlphaSpending,
        userAlphaSpending, spendTime, effStopping);
    }
  }

  // compute transformed quantities for adaptation
  size_t k1 = kMax - L;
  std::vector<double> s1(k1);
  FlatMatrix c1(M, k1);
  for (size_t i = 0; i < k1; ++i) {
    s1[i] = (infoRates[L + i] - infoRates[L - 1]) / (1.0 - infoRates[L - 1]);
    double* colptr = c1.data_ptr() + i * M; // start of column i
    if (effStopping[L + i]) {
      double cut = critValues[L + i];
      double r1 = infoRates[L - 1] / infoRates[L + i];
      double sqrt_r1 = std::sqrt(r1);
      double denom = std::sqrt(1.0 - r1);
      // write contiguous column
      for (size_t m = 0; m < M; ++m) {
        colptr[m] = (cut - zL[m] * sqrt_r1) / denom;
      }
    } else {
      std::fill_n(colptr, M, 8.0);
    }
  }

  // secondary trial information for original design
  double INew1 = IMax * (1.0 - infoRates[L - 1]);
  std::vector<double> I1(k1);
  for (size_t i = 0; i < k1; ++i) {
    I1[i] = INew1 * s1[i];
  }

  // conditional type I error
  auto v0 = exitprob_mams_cpp(M, r, zero, corr_known, k1, c1, I1);
  double c_alpha = std::accumulate(v0.begin(), v0.end(), 0.0);

  // secondary trial design and power calculations
  double IL = IMax * infoRates[L - 1];
  double sqrtIL = std::sqrt(IL);
  std::vector<double> zero2(MNew, 0.0);
  std::vector<double> theta2 = subset(theta, selectedNew);

  size_t k2 = MullerSchafer ? kNew : k1;
  std::vector<double> s2(k2);
  std::string asf2;
  std::vector<unsigned char> effStopping2(k2);
  std::vector<double> cpu0(k2);
  if (!MullerSchafer) {
    // number of secondary trial looks is the same as the original trial
    // spend the conditional type I error as in the original trial
    s2 = s1;
    asf2 = "user";
    std::memcpy(effStopping2.data(), effStopping.data() + L, k2);
    std::partial_sum(v0.begin(), v0.end(), cpu0.begin());
  } else {
    s2 = infoRatesNew;
    asf2 = asfNew;
    effStopping2 = effStoppingNew;
    if (asf2 != "none" && asf2 != "of") {
      for (size_t i = 0; i < k2; ++i) {
        cpu0[i] = errorSpentcpp(spendTimeNew[i], c_alpha, asfNew,
                                parameterAlphaSpendingNew);
      }
    }
  }

  std::vector<double> I2(k2); // information levels for secondary trial
  std::vector<double> Ic(k2); // information levels for integrated trial
  std::vector<double> sqrtI2(k2), sqrtIc(k2);
  FlatMatrix c2(MNew, k2);
  c2.fill(8.0);

  for (size_t i = 0; i < k2; ++i) {
    I2[i] = INew * s2[i];
    Ic[i] = I2[i] + IL;
    sqrtI2[i] = std::sqrt(I2[i]);
    sqrtIc[i] = std::sqrt(Ic[i]);
  }

  std::vector<double> zscaled(MNew);
  for (size_t j = 0; j < MNew; ++j) zscaled[j] = zL[selectedNew[j]] * sqrtIL;

  if (asf2 == "of") {
    auto g = [&c2, &I2, &sqrtI2, &sqrtIc, &zscaled, &zero2, &effStopping2,
              k2, c_alpha, MNew, rNew, corr_known]
    (double x)->double {
      double col_const = x * sqrtIc[k2 - 1];
      for (size_t i = 0; i < k2; ++i) {
        double* colptr = c2.data_ptr() + i * MNew;
        if (effStopping2[i]) {
          double denom = sqrtI2[i];
          for (size_t j = 0; j < MNew; ++j) {
            colptr[j] = (col_const - zscaled[j]) / denom;
          }
        } else {
          std::fill_n(colptr, MNew, 8.0);
        }
      }

      auto v = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, k2, c2, I2);
      double p0 = std::accumulate(v.begin(), v.end(), 0.0);
      return p0 - c_alpha;
    };

    double cof = brent(g, 0.0, 8.0, 1e-6);
    double col_const = cof * sqrtIc[k2 - 1];
    for (size_t i = 0; i < k2; ++i) {
      double* colptr = c2.data_ptr() + i * MNew;
      if (effStopping2[i]) {
        double denom = sqrtI2[i];
        for (size_t j = 0; j < MNew; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }
      } else {
        std::fill_n(colptr, MNew, 8.0);
      }
    }
  } else if (asf2 == "none") {
    double denom = sqrtI2[k2 - 1];

    auto g = [&c2, &I2, &sqrtIc, &zscaled, &zero2,
              denom, k2, c_alpha, MNew, rNew, corr_known]
    (double x)->double {
      double* colptr = c2.data_ptr() + (k2 - 1) * MNew;
      double col_const = x * sqrtIc[k2 - 1];
      for (size_t j = 0; j < MNew; ++j) {
        colptr[j] = (col_const - zscaled[j]) / denom;
      }

      auto v = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, k2, c2, I2);
      double p0 = std::accumulate(v.begin(), v.end(), 0.0);
      return p0 - c_alpha;
    };

    double cof = brent(g, 0.0, 8.0, 1e-6);
    double* colptr = c2.data_ptr() + (k2 - 1) * MNew;
    double col_const = cof * sqrtIc[k2 - 1];
    for (size_t j = 0; j < MNew; ++j) {
      colptr[j] = (col_const - zscaled[j]) / denom;
    }
  } else {
    for (size_t i = 0; i < k2; ++i) {
      if (!effStopping2[i]) continue;
      double denom = sqrtI2[i];

      auto g = [&c2, &I2, &sqrtIc, &zscaled, &cpu0, &zero2,
                denom, i, MNew, rNew, corr_known]
      (double x)->double {
        double col_const = x * sqrtIc[i];
        double* colptr = c2.data_ptr() + i * MNew;
        // update critical values of the secondary trial at current look
        for (size_t j = 0; j < MNew; ++j) {
          colptr[j] = (col_const - zscaled[j]) / denom;
        }

        auto v = exitprob_mams_cpp(MNew, rNew, zero2, corr_known, i + 1, c2, I2);
        double p0 = std::accumulate(v.begin(), v.end(), 0.0);
        return p0 - cpu0[i];
      };

      double cof = brent(g, 0.0, 8.0, 1e-6);
      double col_const = cof * sqrtIc[i];
      double* colptr = c2.data_ptr() + i * MNew;
      for (size_t j = 0; j < MNew; ++j) {
        colptr[j] = (col_const - zscaled[j]) / denom;
      }
    }
  }

  // compute conditional power of the secondary trial
  auto v2 = exitprob_mams_cpp(MNew, rNew, theta2, true, k2, c2, I2);
  return std::accumulate(v2.begin(), v2.end(), 0.0);
}


//' @title Conditional Power for a Multi-Arm Multi-Stage Design
//' @description Obtains the conditional power for specified incremental
//' information given the interim results, parameter values, and
//' data-dependent changes in the selected treatment(s),
//' the error spending function, as well as the
//' number and spacing of interim looks.
//'
//' @param INew The maximum information for any active arm versus the common
//'   control in the secondary trial.
//' @param M Number of active treatment arms in the primary trial.
//' @param r Randomization ratio of each active arm to the common control
//'   in the primary trial.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics is derived from the randomization ratio \eqn{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is assumed.
//' @param L The interim adaptation look of the primary trial.
//' @param zL The z-test statistics at the interim adaptation look of
//'   the primary trial.
//' @param theta A vector of length \eqn{M} representing the assumed treatment
//'   effects for each active arm versus the common control. The global null
//'   is \eqn{\theta_i = 0} for all \eqn{i}, and alternatives are one-sided:
//'   \eqn{\theta_i > 0} for at least one \eqn{i = 1, \ldots, M}.
//' @param IMax Maximum information for any active arm versus the common
//'   control for the primary trial. Must be provided.
//' @param kMax The maximum number of stages of the primary trial.
//' @param informationRates The information rates of the primary trial.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the max z-test statistic
//'   scale for efficacy stopping for the primary trial. If missing, boundaries
//'   will be computed based on the specified alpha spending function.
//' @param alpha The significance level of the primary trial.
//'   Defaults to 0.025.
//' @param typeAlphaSpending The type of alpha spending for the primary
//'   trial. One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function,
//'   \code{"user"} for user defined spending, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpending The parameter value of alpha spending
//'   for the primary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param userAlphaSpending The user-defined alpha spending for the
//'   primary trial. Represents the cumulative alpha spent up to each stage.
//' @param spendingTime The error spending time of the primary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRates}.
//' @param MullerSchafer Whether to use the Muller and Schafer (2001) method
//'   for trial adaptation.
//' @param MNew Number of active treatment arms in the secondary trial.
//' @param selected The indices of the selected active treatment arms for
//'   the secondary trial.
//' @param rNew Randomization ratio of each active arm to the common control
//'   in the secondary trial.
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param typeAlphaSpendingNew The type of alpha spending for the secondary
//'   trial. One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpendingNew The parameter value of alpha spending
//'   for the secondary trial. Corresponds to
//'   \eqn{\rho} for \code{"sfKD"}, and \eqn{\gamma} for \code{"sfHSD"}.
//' @param spendingTimeNew The error spending time of the secondary trial.
//'   Defaults to missing, in which case it is assumed to be the same as
//'   \code{informationRatesNew}.
//'
//' @return The conditional power given the interim results, parameter
//' values, and data-dependent design changes.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @references
//' Ping Gao, Yingqiu Li.
//' Adaptive multiple comparison sequential design (AMCSD) for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2024, 34(3), 424-440.
//'
//' @seealso \code{\link{adaptDesign_mams}}
//'
//' @examples
//'
//' getCP_mams(
//'   INew = 373 / 4, M = 2, r = 1, corr_known = FALSE,
//'   L = 1, zL = c(-log(0.91), -log(0.78)) * sqrt(324 / 4 / 2),
//'   theta = c(-log(0.91), -log(0.78)),
//'   IMax = 324 / 4, kMax = 2, informationRates = c(1/2, 1),
//'   alpha = 0.025, typeAlphaSpending = "OF",
//'   MNew = 1, selected = 2, rNew = 1)
//'
//' @export
// [[Rcpp::export]]
double getCP_mams(const double INew = NA_REAL,
                  const int M = NA_INTEGER,
                  const double r = 1,
                  const bool corr_known = true,
                  const int L = NA_INTEGER,
                  const Rcpp::NumericVector& zL = NA_REAL,
                  const Rcpp::NumericVector& theta = NA_REAL,
                  const double IMax = NA_REAL,
                  const int kMax = NA_INTEGER,
                  const Rcpp::NumericVector& informationRates = NA_REAL,
                  const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
                  const Rcpp::NumericVector& criticalValues = NA_REAL,
                  const double alpha = 0.025,
                  const std::string& typeAlphaSpending = "sfOF",
                  const double parameterAlphaSpending = NA_REAL,
                  const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
                  const Rcpp::NumericVector& spendingTime = NA_REAL,
                  const bool MullerSchafer = false,
                  const int MNew = NA_INTEGER,
                  const Rcpp::IntegerVector& selected = NA_INTEGER,
                  const double rNew = 1,
                  const int kNew = NA_INTEGER,
                  const Rcpp::NumericVector& informationRatesNew = NA_REAL,
                  const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
                  const std::string& typeAlphaSpendingNew = "sfOF",
                  const double parameterAlphaSpendingNew = NA_REAL,
                  const Rcpp::NumericVector& spendingTimeNew = NA_REAL) {

  auto zLVec = Rcpp::as<std::vector<double>>(zL);
  auto thetaVec = Rcpp::as<std::vector<double>>(theta);
  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto selectedNew = Rcpp::as<std::vector<int>>(selected);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);

  auto result = getCP_mams_cpp(
    INew,  static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(L), zLVec, thetaVec, IMax,
    static_cast<size_t>(kMax), infoRates, effStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, spendTime,
    MullerSchafer, static_cast<size_t>(MNew), selectedNew, rNew,
    static_cast<size_t>(kNew), infoRatesNew, effStoppingNew,
    typeAlphaSpendingNew, parameterAlphaSpendingNew, spendTimeNew
  );

  return result;
}
