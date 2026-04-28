#include "generic_design.h"
#include "seamless_design.h"
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
double getCP_seamless_cpp(
    const double INew,
    const size_t M,
    const double r,
    const bool corr_known,
    const size_t L,
    const double zL,
    const double theta,
    const double IMax,
    const size_t K,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& spendingTime,
    const bool MullerSchafer,
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
  if (std::isnan(zL)) throw std::invalid_argument("zL must be provided");
  if (std::isnan(theta)) throw std::invalid_argument("theta must be provided");
  if (std::isnan(IMax)) throw std::invalid_argument("IMax must be provided");
  if (IMax <= 0.0) throw std::invalid_argument("IMax must be positive");
  if (K <= L) throw std::invalid_argument("K must be greater than L");

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }

  size_t kMax = K + 1;

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
  std::vector<double> infoRatesNew;
  std::vector<unsigned char> effStoppingNew;
  std::vector<double> spendTimeNew = spendingTimeNew;

  std::string asfNew = typeAlphaSpendingNew;
  for (char &c : asfNew) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

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
      if (informationRatesNew.back() != 1.0)
        throw std::invalid_argument("informationRatesNew must end with 1");
    } else {
      for (size_t i = 0; i < kNew; ++i)
        infoRatesNew[i] = static_cast<double>(i+1) / static_cast<double>(kNew);
    }

    // effStoppingNew: default to all 1s if missing
    if (none_na(efficacyStoppingNew)) {
      if (efficacyStoppingNew.size() != kNew)
        throw std::invalid_argument("Invalid length for efficacyStoppingNew");
      if (efficacyStoppingNew.back() != 1)
        throw std::invalid_argument("efficacyStoppingNew must end with 1");
      effStoppingNew = efficacyStoppingNew; // copy
    } else {
      effStoppingNew.assign(kNew, 1);
    }

    if (!(asfNew == "of" || asfNew == "p" || asfNew == "wt" ||
        asfNew == "sfof" || asfNew == "sfp" ||
        asfNew == "sfkd" || asfNew == "sfhsd" || asfNew == "none")) {
      throw std::invalid_argument("Invalid value for typeAlphaSpendingNew");
    }
    if ((asfNew == "wt" || asfNew == "sfkd" || asfNew == "sfhsd") &&
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
      if (spendingTimeNew.back() != 1.0)
        throw std::invalid_argument("spendingTimeNew must end with 1");
    } else {
      spendTimeNew = infoRatesNew;
    }
  }

  // ----------- End of Input Validation ----------- //

  // set up efficacy bounds
  std::vector<double> zero(M, 0.0);
  std::vector<double> critValues = criticalValues;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (kMax > 1 && criticalValues.size() == kMax) {
      bool hasNaN = false;
      for (size_t i = 0; i < kMax - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[kMax-1])) haybittle = true;
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(kMax);
      for (size_t i = 0; i < kMax - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[kMax-1] = aval;
        ListCpp probs = exitprob_seamless_cpp(
          M, r, zero, corr_known, K, u, infoRates);
        auto v = probs.get<std::vector<double>>("exitProb");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[kMax-1] = brent(f, 0.0, 6.0, 1e-6);
    } else {
      critValues = getBound_seamless_cpp(
        M, r, corr_known, K, infoRates, alpha, asf, parameterAlphaSpending,
        userAlphaSpending, spendTime, effStopping);
    }
  }


  // compute transformed quantities for adaptation
  size_t k1 = K - L;
  std::vector<double> t1(k1), b1(k1), a1(k1, -6.0), zero1(k1, 0.0);
  for (size_t l = 0; l < k1; ++l) {
    t1[l] = (infoRates[l + L + 1] - infoRates[L]) / (1.0 - infoRates[L]);
    double r1 = infoRates[L] / infoRates[l + L + 1];
    b1[l] = (critValues[l + L + 1] - std::sqrt(r1) * zL) / std::sqrt(1.0 - r1);
    if (!effStopping[l + L + 1]) b1[l] = 6.0;
  }

  a1[k1-1] = b1[k1-1];

  double result = 0.0;
  if (!MullerSchafer) {
    std::vector<double> mu(k1, theta), I2(k1);
    for (size_t l = 0; l < k1; ++l) {
      I2[l] = INew * t1[l];
    }

    ListCpp probs = exitprobcpp(b1, a1, mu, I2);
    auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
    result = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);
  } else {
    // Muller-Schafer branch

    // conditional type I error
    ListCpp probs = exitprobcpp(b1, a1, zero1, t1);
    auto exitUpper = probs.get<std::vector<double>>("exitProbUpper");
    double alphaNew = std::accumulate(exitUpper.begin(), exitUpper.end(), 0.0);

    // efficacy boundaries for secondary trial
    auto b2 = getBoundcpp(kNew, informationRatesNew, alphaNew, asfNew,
                          parameterAlphaSpendingNew, std::vector<double>{},
                          spendTimeNew, effStoppingNew);

    // conditional power
    std::vector<double> mu(kNew, theta), I2(kNew);
    for (size_t l = 0; l < kNew; ++l) {
      I2[l] = INew * informationRatesNew[l];
    }

    std::string bsfNew = "none";
    double parameterBetaSpendingNew = NaN;
    std::vector<unsigned char> futStoppingNew(kNew, 1);
    std::vector<double> w1(kNew, 1.0);

    ListCpp out = getPower(alphaNew, kNew, b2, mu, I2, bsfNew,
                           parameterBetaSpendingNew, spendTimeNew,
                           futStoppingNew, w1);
    result = out.get<double>("power");
  }

  return result;
}


//' @title Conditional Power for Two-Stage Seamless Sequential Design
//' @description Obtains the conditional power for specified incremental
//' information given the interim results, parameter values, and
//' data-dependent changes in the error spending function, as well as the
//' number and spacing of interim looks.
//'
//' @param INew The maximum information for any active arm versus the common
//'   control in the secondary trial.
//' @param M Number of active treatment arms in Phase 2.
//' @param r Randomization ratio of each active arm to the common control
//'   in Phase 2.
//' @param corr_known Logical. If \code{TRUE}, the correlation between Wald
//'   statistics in Phase 2 is derived from the randomization ratio \code{r}
//'   as \eqn{r / (r + 1)}. If \code{FALSE}, a conservative correlation of
//'   0 is used.
//' @param L The interim adaptation look in Phase 3.
//' @param zL The z-test statistic at the interim adaptation look of
//'   Phase 3.
//' @param theta The treatment effect for the selected arm versus the
//'   common control.
//' @param IMax Maximum information for any active arm versus the common
//'   control for the original trial. Must be provided.
//' @param K Number of sequential looks in Phase 3.
//' @param informationRates A numeric vector of information rates fixed
//'   before the trial. If unspecified, defaults to \eqn{(1:(K+1)) / (K+1)}.
//' @param efficacyStopping Indicators of whether efficacy stopping is
//'   allowed at each stage of the primary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param criticalValues The upper boundaries on the max z-test statistic
//'   scale for Phase 2 and the z-test statistics for the selected arm
//'   in Phase 3 for the primary trial. If missing, boundaries
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
//' @param kNew The number of looks of the secondary trial.
//' @param informationRatesNew The spacing of looks of the secondary trial.
//' @param efficacyStoppingNew The indicators of whether efficacy stopping is
//'   allowed at each look of the secondary trial. Defaults to \code{TRUE}
//'   if left unspecified.
//' @param typeAlphaSpendingNew The type of alpha spending for the secondary
//'   trial. One of the following:
//'   \code{"OF"} for O'Brien-Fleming boundaries,
//'   \code{"P"} for Pocock boundaries,
//'   \code{"WT"} for Wang & Tsiatis boundaries,
//'   \code{"sfOF"} for O'Brien-Fleming type spending function,
//'   \code{"sfP"} for Pocock type spending function,
//'   \code{"sfKD"} for Kim & DeMets spending function,
//'   \code{"sfHSD"} for Hwang, Shi & DeCani spending function, and
//'   \code{"none"} for no early efficacy stopping.
//'   Defaults to \code{"sfOF"}.
//' @param parameterAlphaSpendingNew The parameter value of alpha spending
//'   for the secondary trial. Corresponds to \eqn{\Delta} for \code{"WT"},
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
//' Adaptive two-stage seamless sequential design for clinical trials.
//' Journal of Biopharmaceutical Statistics, 2025, 35(4), 565-587.
//'
//' @seealso \code{\link{adaptDesign_seamless}}
//'
//' @examples
//'
//' getCP_seamless(
//'   INew = 49.5, M = 2, r = 1, corr_known = FALSE,
//'   L = 1, zL = -log(0.67) * sqrt(80 / 4), theta = -log(0.691),
//'   IMax = 120 / 4, K = 2, informationRates = c(1/3, 2/3, 1),
//'   alpha = 0.025, typeAlphaSpending = "OF", kNew = 1)
//'
//' @export
// [[Rcpp::export]]
double getCP_seamless(const double INew = NA_REAL,
                      const int M = NA_INTEGER,
                      const double r = NA_REAL,
                      const bool corr_known = true,
                      const int L = NA_INTEGER,
                      const double zL = NA_REAL,
                      const double theta = NA_REAL,
                      const double IMax = NA_REAL,
                      const int K = NA_INTEGER,
                      const Rcpp::NumericVector& informationRates = NA_REAL,
                      const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
                      const Rcpp::NumericVector& criticalValues = NA_REAL,
                      const double alpha = 0.025,
                      const std::string& typeAlphaSpending = "sfOF",
                      const double parameterAlphaSpending = NA_REAL,
                      const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
                      const Rcpp::NumericVector& spendingTime = NA_REAL,
                      const bool MullerSchafer = false,
                      const int kNew = NA_INTEGER,
                      const Rcpp::NumericVector& informationRatesNew = NA_REAL,
                      const Rcpp::LogicalVector& efficacyStoppingNew = NA_LOGICAL,
                      const std::string& typeAlphaSpendingNew = "sfOF",
                      const double parameterAlphaSpendingNew = NA_REAL,
                      const Rcpp::NumericVector& spendingTimeNew = NA_REAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);
  auto infoRatesNew = Rcpp::as<std::vector<double>>(informationRatesNew);
  auto effStoppingNew = convertLogicalVector(efficacyStoppingNew);
  auto spendTimeNew = Rcpp::as<std::vector<double>>(spendingTimeNew);

  auto result = getCP_seamless_cpp(
    INew,  static_cast<size_t>(M), r, corr_known,
    static_cast<size_t>(L), zL, theta, IMax,
    static_cast<size_t>(K), infoRates, effStopping, critValues,
    alpha, typeAlphaSpending, parameterAlphaSpending, userAlpha, spendTime,
    MullerSchafer, static_cast<size_t>(kNew), infoRatesNew, effStoppingNew,
    typeAlphaSpendingNew, parameterAlphaSpendingNew, spendTimeNew
  );

  return result;
}
