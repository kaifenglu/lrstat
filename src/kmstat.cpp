#include "enrollment_event.h"
#include "generic_design.h"
#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Rcpp.h>

using std::size_t;

std::vector<double> make_breaks(
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& accrualTime,
    const double accrualDuration,
    const double maxFollowupTime,
    const double time,
    const double milestone) {

  double upper = std::min(std::min(milestone, time), maxFollowupTime);

  // Candidate set
  std::vector<double> pts;
  pts.reserve(piecewiseSurvivalTime.size() + accrualTime.size());

  // add piecewise survival cut points (exclude 0 and maxFollowupTime)
  for (double t : piecewiseSurvivalTime) {
    if (t <= 0.0 || t >= upper) continue;
    pts.push_back(t);
  }

  // add accrual-derived points: time - accrualTime[i]
  for (double at : accrualTime) {
    double s = time - at;
    if (s <= 0.0 || s >= upper) continue;
    pts.push_back(s);
  }

  double t_1 = time - accrualDuration;
  if (t_1 > 0.0 && t_1 < upper) pts.push_back(t_1);

  if (pts.empty()) {
    // no internal break points -> simple [0, max]
    return std::vector<double>{0.0, upper};
  }

  // sort and deduplicate with tolerance
  std::sort(pts.begin(), pts.end());
  double eps = std::max(1e-12, 1e-15 * std::max(1.0, upper));

  std::vector<double> uniques;
  uniques.reserve(pts.size());
  double last = pts.front();
  uniques.push_back(last);
  for (size_t i = 1; i < pts.size(); ++i) {
    double v = pts[i];
    if (v - last > eps) {
      uniques.push_back(v);
      last = v;
    }
  }

  // build final breaks with endpoints
  std::vector<double> breaks;
  breaks.reserve(uniques.size() + 2);
  breaks.push_back(0.0);
  for (double v : uniques) breaks.push_back(v);
  breaks.push_back(upper);

  // Optionally remove intervals of zero (or tiny) width:
  std::vector<double> final_breaks;
  final_breaks.reserve(breaks.size());
  final_breaks.push_back(breaks.front());
  for (size_t i = 1; i < breaks.size(); ++i) {
    if (breaks[i] - final_breaks.back() > eps) final_breaks.push_back(breaks[i]);
    // else skip near-duplicate boundary
  }
  return final_breaks;
}


double km_integrand(
  double x, // analysis time
  double time, // calendar time of analysis
  double phi,
  const std::vector<double>& accrualTime,
  const std::vector<double>& accrualIntensity,
  const std::vector<double>& piecewiseSurvivalTime,
  const std::vector<double>& lambda,
  const std::vector<double>& gamma,
  double accrualDuration) {

  // interval index j for x
  // findInterval1 returns an index in [0..nv], consistent with R's findInterval;
  int j = findInterval1(x, piecewiseSurvivalTime) - 1;
  if (j < 0) j = 0;
  if (j >= static_cast<int>(lambda.size())) j = static_cast<int>(lambda.size()) - 1;

  // p(x) = P(at risk at time x since randomization) under event+dropout hazards
  double p = patrisk1(x, piecewiseSurvivalTime, lambda, gamma);

  // N(time - x): number enrolled by calendar time (time-x)
  double tx = time - x;
  double N = accrual1(tx, accrualTime, accrualIntensity, accrualDuration);

  // lambda_j / (phi * N * p)
  // guard against division by 0
  double denom = phi * N * p;
  if (denom <= 0.0) return 0.0;

  return lambda[j] / denom;
}


DataFrameCpp kmstat1cpp(
    const double time,
    const double milestone,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<std::vector<double>>& lambda1,
    const std::vector<std::vector<double>>& lambda2,
    const std::vector<std::vector<double>>& gamma1,
    const std::vector<std::vector<double>>& gamma2,
    const double accrualDuration,
    const double followupTime,
    const bool fixedFollowup) {

  const std::size_t nstrata = stratumFraction.size();
  const std::size_t nintv = piecewiseSurvivalTime.size();
  std::vector<double> zerogam(nintv, 0.0);

  // phi = P(randomized to group 1)
  const double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // max follow-up time for first enrolled subject
  const double maxFollowupTime = fixedFollowup ? followupTime :
    (accrualDuration + followupTime);

  // total enrolled by calendar time 'time'
  const double a = accrual1(time, accrualTime, accrualIntensity, accrualDuration);

  // enrolled by (time - milestone) used for reaching milestone
  const double a2 = accrual1(time - milestone, accrualTime, accrualIntensity,
                             accrualDuration);

  // ---- outputs ----
  std::vector<int> stratum(nstrata);
  std::vector<double> calTime(nstrata, time);
  std::vector<double> mileTime(nstrata, milestone);

  std::vector<double> nsubjects(nstrata);
  std::vector<double> nevents(nstrata), nevents1(nstrata), nevents2(nstrata);
  std::vector<double> ndropouts(nstrata), ndropouts1(nstrata), ndropouts2(nstrata);

  std::vector<double> nmilestone(nstrata), nmilestone1(nstrata), nmilestone2(nstrata);

  std::vector<double> surv1(nstrata), surv2(nstrata), survDiff(nstrata);
  std::vector<double> vsurv1(nstrata), vsurv2(nstrata), vsurvDiff(nstrata);

  auto breaks = make_breaks(piecewiseSurvivalTime, accrualTime,
                            accrualDuration, maxFollowupTime,
                            time, milestone);
  double tol = 1e-6;

  for (std::size_t h = 0; h < nstrata; ++h) {
    stratum[h] = static_cast<int>(h + 1);
    const double frac = stratumFraction[h];

    // subset per-stratum hazards
    const std::vector<double>& lam1 = lambda1[h];
    const std::vector<double>& lam2 = lambda2[h];
    const std::vector<double>& gam1 = gamma1[h];
    const std::vector<double>& gam2 = gamma2[h];

    // scale accrualIntensity by stratum fraction
    std::vector<double> accrualIntensity_frac = accrualIntensity;
    for (double& v : accrualIntensity_frac) v *= frac;

    // events and dropouts at calendar time
    auto ne_row = nevent21cpp(
      time, allocationRatioPlanned, accrualTime, accrualIntensity_frac,
      piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
      accrualDuration, followupTime, maxFollowupTime);

    auto nd_row = nevent21cpp(
      time, allocationRatioPlanned, accrualTime, accrualIntensity_frac,
      piecewiseSurvivalTime, gam1, gam2, lam1, lam2,   // swapped
      accrualDuration, followupTime, maxFollowupTime);

    nsubjects[h]  = frac * a;
    nevents1[h]   = ne_row.first;
    nevents2[h]   = ne_row.second;
    nevents[h]    = nevents1[h] + nevents2[h];
    ndropouts1[h] = nd_row.first;
    ndropouts2[h] = nd_row.second;
    ndropouts[h]  = ndropouts1[h] + ndropouts2[h];

    // number reaching milestone (uses risk prob at milestone under event+dropout)
    double ncom = frac * a2;
    double p1 = patrisk1(milestone, piecewiseSurvivalTime, lam1, gam1);
    double p2 = patrisk1(milestone, piecewiseSurvivalTime, lam2, gam2);
    nmilestone1[h] = phi * ncom * p1;
    nmilestone2[h] = (1.0 - phi) * ncom * p2;
    nmilestone[h]  = nmilestone1[h] + nmilestone2[h];

    // milestone survival (event only: gamma = 0)
    surv1[h] = patrisk1(milestone, piecewiseSurvivalTime, lam1, zerogam);
    surv2[h] = patrisk1(milestone, piecewiseSurvivalTime, lam2, zerogam);
    survDiff[h] = surv1[h] - surv2[h];

    // variance integrals q1, q2
    auto f1 = [&](double x)->double {
      return km_integrand(
        x, time, phi, accrualTime, accrualIntensity_frac,
        piecewiseSurvivalTime, lam1, gam1, accrualDuration);
    };
    auto f2 = [&](double x)->double {
      return km_integrand(
        x, time, 1.0 - phi, accrualTime, accrualIntensity_frac,
        piecewiseSurvivalTime, lam2, gam2, accrualDuration);
    };

    double q1 = integrate3(f1, breaks, tol);
    double q2 = integrate3(f2, breaks, tol);

    vsurv1[h] = surv1[h] * surv1[h] * q1;
    vsurv2[h] = surv2[h] * surv2[h] * q2;
    vsurvDiff[h] = vsurv1[h] + vsurv2[h];
  }

  // build DataFrameCpp
  DataFrameCpp df;
  df.push_back(std::move(stratum), "stratum");
  df.push_back(std::move(calTime), "time");
  df.push_back(std::move(nsubjects), "subjects");
  df.push_back(std::move(nevents), "nevents");
  df.push_back(std::move(nevents1), "nevents1");
  df.push_back(std::move(nevents2), "nevents2");
  df.push_back(std::move(ndropouts), "ndropouts");
  df.push_back(std::move(ndropouts1), "ndropouts1");
  df.push_back(std::move(ndropouts2), "ndropouts2");
  df.push_back(std::move(mileTime), "milestone");
  df.push_back(std::move(nmilestone), "nmilestone");
  df.push_back(std::move(nmilestone1), "nmilestone1");
  df.push_back(std::move(nmilestone2), "nmilestone2");
  df.push_back(std::move(surv1), "surv1");
  df.push_back(std::move(surv2), "surv2");
  df.push_back(std::move(survDiff), "survDiff");
  df.push_back(std::move(vsurv1), "vsurv1");
  df.push_back(std::move(vsurv2), "vsurv2");
  df.push_back(std::move(vsurvDiff), "vsurvDiff");

  return df;
}

double extract_km(const DataFrameCpp& df, const char* name) {
  auto vec = df.get<double>(name);
  return std::accumulate(vec.begin(), vec.end(), 0.0);
};


DataFrameCpp kmstatcpp(
    const std::vector<double>& time,
    const double milestone,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double followupTime,
    const bool fixedFollowup) {

  if (!none_na(time)) throw std::invalid_argument("time must be provided");
  for (double v : time) {
    if (v <= 0.0) throw std::invalid_argument("time must be positive");
  }
  if (std::isnan(milestone)) throw std::invalid_argument("milestone must be provided");
  if (milestone <= 0.0) throw std::invalid_argument("milestone must be positive");
  if (allocationRatioPlanned <= 0.0)
    throw std::invalid_argument("allocationRatioPlanned must be positive");
  if (accrualTime[0] != 0.0)
    throw std::invalid_argument("accrualTime must start with 0");
  if (any_nonincreasing(accrualTime))
    throw std::invalid_argument("accrualTime should be increasing");
  if (!none_na(accrualIntensity))
    throw std::invalid_argument("accrualIntensity must be provided");
  if (accrualTime.size() != accrualIntensity.size())
    throw std::invalid_argument("Invalid length for accrualIntensity");
  for (double v : accrualIntensity) {
    if (v < 0.0) throw std::invalid_argument("accrualIntensity must be non-negative");
  }
  if (piecewiseSurvivalTime[0] != 0.0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");
  for (double v : stratumFraction) {
    if (v <= 0.0) throw std::invalid_argument("stratumFraction must be positive");
  }
  double sumf = std::accumulate(stratumFraction.begin(), stratumFraction.end(), 0.0);
  if (std::fabs(sumf - 1.0) > 1e-12)
    throw std::invalid_argument("stratumFraction must sum to 1");
  if (!none_na(lambda1)) throw std::invalid_argument("lambda1 must be provided");
  if (!none_na(lambda2)) throw std::invalid_argument("lambda2 must be provided");
  for (double v : lambda1) {
    if (v < 0.0) throw std::invalid_argument("lambda1 must be non-negative");
  }
  for (double v : lambda2) {
    if (v < 0.0) throw std::invalid_argument("lambda2 must be non-negative");
  }
  for (double v : gamma1) {
    if (v < 0.0) throw std::invalid_argument("gamma1 must be non-negative");
  }
  for (double v : gamma2) {
    if (v < 0.0) throw std::invalid_argument("gamma2 must be non-negative");
  }
  if (std::isnan(accrualDuration))
    throw std::invalid_argument("accrualDuration must be provided");
  if (accrualDuration <= 0.0)
    throw std::invalid_argument("accrualDuration must be positive");

  if (std::isnan(followupTime))
    throw std::invalid_argument("followupTime must be provided");
  if (fixedFollowup && followupTime <= 0.0)
    throw std::invalid_argument("followupTime must be positive for fixed follow-up");
  if (!fixedFollowup && followupTime < 0.0)
    throw std::invalid_argument(
        "followupTime must be non-negative for variable follow-up");
  for (double t : time) {
    if (t > accrualDuration + followupTime)
      throw std::invalid_argument("time cannot exceed accrualDuration + followupTime");
    if (t <= milestone)
      throw std::invalid_argument("time must be greater than milestone");
  }
  if (fixedFollowup && (milestone > followupTime))
    throw std::invalid_argument(
        "milestone must be <= followupTime for fixed follow-up");
  if (!fixedFollowup && (milestone >= accrualDuration + followupTime))
    throw std::invalid_argument(
        "milestone must be < accrualDuration + followupTime for variable follow-up");

  // expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambda1x = expand_stratified(lambda1, nstrata, nintv, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintv, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintv, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintv, "gamma2");

  // ---- outputs ----
  const size_t k = time.size();
  std::vector<double> calTime = time;
  std::vector<double> mileTime(k, milestone);
  std::vector<double> subjects(k), nevents(k), nevents1(k), nevents2(k);
  std::vector<double> ndropouts(k), ndropouts1(k), ndropouts2(k);
  std::vector<double> nmilestone(k), nmilestone1(k), nmilestone2(k);
  std::vector<double> surv1(k), surv2(k), vsurv1(k), vsurv2(k);
  std::vector<double> survDiff(k), vsurvDiff(k), information(k), survDiffZ(k);

  for (size_t j = 0; j < k; ++j) {
    DataFrameCpp df = kmstat1cpp(
      time[j], milestone, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup);

    subjects[j]     = extract_km(df, "subjects");
    nevents[j]      = extract_km(df, "nevents");
    nevents1[j]     = extract_km(df, "nevents1");
    nevents2[j]     = extract_km(df, "nevents2");
    ndropouts[j]    = extract_km(df, "ndropouts");
    ndropouts1[j]   = extract_km(df, "ndropouts1");
    ndropouts2[j]   = extract_km(df, "ndropouts2");
    nmilestone[j]   = extract_km(df, "nmilestone");
    nmilestone1[j]  = extract_km(df, "nmilestone1");
    nmilestone2[j]  = extract_km(df, "nmilestone2");

    const auto& s1  = df.get<double>("surv1");
    const auto& s2  = df.get<double>("surv2");
    const auto& vs1 = df.get<double>("vsurv1");
    const auto& vs2 = df.get<double>("vsurv2");
    for (size_t h = 0; h < nstrata; ++h) {
      surv1[j]  += stratumFraction[h] * s1[h];
      surv2[j]  += stratumFraction[h] * s2[h];
      vsurv1[j] += stratumFraction[h] * stratumFraction[h] * vs1[h];
      vsurv2[j] += stratumFraction[h] * stratumFraction[h] * vs2[h];
    }

    survDiff[j]    = surv1[j] - surv2[j];
    vsurvDiff[j]   = vsurv1[j] + vsurv2[j];
    information[j] = 1.0 / vsurvDiff[j];
    survDiffZ[j]   = survDiff[j] * std::sqrt(information[j]);
  }

  DataFrameCpp out;
  out.push_back(std::move(calTime), "time");
  out.push_back(std::move(subjects), "subjects");
  out.push_back(std::move(nevents), "nevents");
  out.push_back(std::move(nevents1), "nevents1");
  out.push_back(std::move(nevents2), "nevents2");
  out.push_back(std::move(ndropouts), "ndropouts");
  out.push_back(std::move(ndropouts1), "ndropouts1");
  out.push_back(std::move(ndropouts2), "ndropouts2");
  out.push_back(std::move(mileTime), "milestone");
  out.push_back(std::move(nmilestone), "nmilestone");
  out.push_back(std::move(nmilestone1), "nmilestone1");
  out.push_back(std::move(nmilestone2), "nmilestone2");
  out.push_back(std::move(surv1), "surv1");
  out.push_back(std::move(surv2), "surv2");
  out.push_back(std::move(survDiff), "survDiff");
  out.push_back(std::move(vsurv1), "vsurv1");
  out.push_back(std::move(vsurv2), "vsurv2");
  out.push_back(std::move(vsurvDiff), "vsurvDiff");
  out.push_back(std::move(information), "information");
  out.push_back(std::move(survDiffZ), "survDiffZ");

  return out;
}


//' @title Stratified Difference in Milestone Survival Probabilities
//' @description Obtains the stratified milestone survival probabilities
//' and difference in milestone survival probabilities at given
//' calendar times.
//'
//' @param time A vector of calendar times for data cut.
//' @param milestone The milestone time at which to calculate the
//'   survival probability.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//'
//' @return A data frame containing the following variables:
//'
//' * \code{time}: The calendar time since trial start.
//'
//' * \code{subjects}: The number of enrolled subjects.
//'
//' * \code{nevents}: The total number of events.
//'
//' * \code{nevents1}: The number of events in the active treatment group.
//'
//' * \code{nevents2}: The number of events in the control group.
//'
//' * \code{ndropouts}: The total number of dropouts.
//'
//' * \code{ndropouts1}: The number of dropouts in the active treatment
//'   group.
//'
//' * \code{ndropouts2}: The number of dropouts in the control group.
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{nmilestone}: The total number of subjects reaching milestone.
//'
//' * \code{nmilestone1}: The number of subjects reaching milestone
//'   in the active treatment group.
//'
//' * \code{nmiletone2}: The number of subjects reaching milestone
//'   in the control group.
//'
//' * \code{surv1}: The milestone survival probability for the treatment
//'   group.
//'
//' * \code{surv2}: The milestone survival probability for the control group.
//'
//' * \code{survDiff}: The difference in milestone survival probabilities,
//'   i.e., \code{surv1 - surv2}.
//'
//' * \code{vsurv1}: The variance for \code{surv1}.
//'
//' * \code{vsurv2}: The variance for \code{surv2}.
//'
//' * \code{vsurvDiff}: The variance for \code{survDiff}.
//'
//' * \code{information}: The information for \code{survDiff}, equal to
//'   \code{1/vsurvDiff}.
//'
//' * \code{survDiffZ}: The Z-statistic value, i.e.,
//'   \code{survDiff/sqrt(vsurvDiff)}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' kmstat(time = c(22, 40),
//'        milestone = 18,
//'        allocationRatioPlanned = 1,
//'        accrualTime = seq(0, 8),
//'        accrualIntensity = 26/9*seq(1, 9),
//'        piecewiseSurvivalTime = c(0, 6),
//'        stratumFraction = c(0.2, 0.8),
//'        lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'        lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'        gamma1 = -log(1-0.05)/12,
//'        gamma2 = -log(1-0.05)/12,
//'        accrualDuration = 22,
//'        followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame kmstat(
    const Rcpp::NumericVector& time = NA_REAL,
    const double milestone = NA_REAL,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double accrualDuration = NA_REAL,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = false) {

  auto time1 = Rcpp::as<std::vector<double>>(time);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  DataFrameCpp df = kmstatcpp(
    time1, milestone, allocationRatioPlanned,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup);

  return Rcpp::wrap(df);
}


ListCpp kmpowercpp(
    const int kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<unsigned char>& futilityStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& futilityBounds,
    const std::string& typeBetaSpending,
    const double parameterBetaSpending,
    const double milestone,
    const double survDiffH0,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double followupTime,
    const bool fixedFollowup,
    const std::vector<double>& spendingTime,
    const double studyDuration) {

  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }
  size_t K = static_cast<size_t>(kMax);

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(K);
  if (none_na(informationRates)) {
    if (informationRates.size() != K)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[K-1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else {
    for (size_t i = 0; i < K; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(K);
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != K)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[K-1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(K, 1);
  }

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != K)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping[K-1] != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(K, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds);

  if (!missingCriticalValues && criticalValues.size() != K) {
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
    if (userAlphaSpending.size() != K)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending[K-1] != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  if (!missingFutilityBounds) {
    if (!(futilityBounds.size() == K - 1 || futilityBounds.size() == K)) {
      throw std::invalid_argument("Invalid length for futilityBounds");
    }
  }
  if (!missingCriticalValues && !missingFutilityBounds) {
    for (size_t i = 0; i < K - 1; ++i) {
      if (futilityBounds[i] > criticalValues[i]) {
        throw std::invalid_argument("futilityBounds must lie below criticalValues");
      }
    }
    if (futilityBounds.size() == K && futilityBounds[K-1] != criticalValues[K-1]) {
      throw std::invalid_argument(
          "futilityBounds must meet criticalValues at the final look");
    }
  }

  std::string bsf = typeBetaSpending;
  for (char &c : bsf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
      bsf == "sfkd" || bsf == "sfhsd" || bsf == "none")) {
    throw std::invalid_argument("Invalid value for typeBetaSpending");
  }
  if ((bsf == "sfkd" || bsf == "sfhsd") && std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument("parameterBetaSpending must be positive for sfKD");
  }

  // milestone + survival diff null
  if (std::isnan(milestone)) throw std::invalid_argument("milestone must be provided");
  if (milestone <= 0.0) throw std::invalid_argument("milestone must be positive");
  if (survDiffH0 <= -1.0 || survDiffH0 >= 1.0)
    throw std::invalid_argument("survDiffH0 must lie between -1 and 1");

  if (allocationRatioPlanned <= 0.0)
    throw std::invalid_argument("allocationRatioPlanned must be positive");
  if (accrualTime[0] != 0.0)
    throw std::invalid_argument("accrualTime must start with 0");
  if (any_nonincreasing(accrualTime))
    throw std::invalid_argument("accrualTime should be increasing");

  if (!none_na(accrualIntensity))
    throw std::invalid_argument("accrualIntensity must be provided");
  if (accrualIntensity.size() != accrualTime.size())
    throw std::invalid_argument("Invalid length for accrualIntensity");
  for (double v : accrualIntensity) {
    if (v < 0.0) throw std::invalid_argument("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0.0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");

  for (double v : stratumFraction) {
    if (v <= 0.0) throw std::invalid_argument("stratumFraction must be positive");
  }
  double sumf = std::accumulate(stratumFraction.begin(), stratumFraction.end(), 0.0);
  if (std::fabs(sumf - 1.0) > 1e-12)
    throw std::invalid_argument("stratumFraction must sum to 1");

  if (!none_na(lambda1)) throw std::invalid_argument("lambda1 must be provided");
  if (!none_na(lambda2)) throw std::invalid_argument("lambda2 must be provided");

  for (double v : lambda1) {
    if (v < 0.0) throw std::invalid_argument("lambda1 must be non-negative");
  }
  for (double v : lambda2) {
    if (v < 0.0) throw std::invalid_argument("lambda2 must be non-negative");
  }
  for (double v : gamma1) {
    if (v < 0.0) throw std::invalid_argument("gamma1 must be non-negative");
  }
  for (double v : gamma2) {
    if (v < 0.0) throw std::invalid_argument("gamma2 must be non-negative");
  }

  // expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambda1x = expand_stratified(lambda1, nstrata, nintv, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintv, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintv, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintv, "gamma2");

  if (std::isnan(accrualDuration))
    throw std::invalid_argument("accrualDuration must be provided");
  if (accrualDuration <= 0.0)
    throw std::invalid_argument("accrualDuration must be positive");
  if (std::isnan(followupTime))
    throw std::invalid_argument("followupTime must be provided");
  if (fixedFollowup && followupTime <= 0.0)
    throw std::invalid_argument("followupTime must be positive for fixed follow-up");
  if (!fixedFollowup && followupTime < 0.0)
    throw std::invalid_argument(
        "followupTime must be non-negative for variable follow-up");

  // spendingTime default to informationRates
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != K)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[K-1] != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // fixed follow-up studyDuration constraints
  if (fixedFollowup && !std::isnan(studyDuration) &&
      (studyDuration < accrualDuration))
    throw std::invalid_argument(
        "studyDuration must be greater than or equal to accrualDuration");
  if (fixedFollowup && !std::isnan(studyDuration) &&
      (studyDuration > accrualDuration + followupTime))
    throw std::invalid_argument(
        "studyDuration must be less than or equal to accrualDuration + followupTime");

  if (milestone >= accrualDuration + followupTime)
    throw std::invalid_argument(
        "milestone must be less than accrualDuration + followupTime");
  if (fixedFollowup && milestone > followupTime)
    throw std::invalid_argument(
        "milestone cannot exceed followupTime for fixed follow-up");
  if (fixedFollowup && !std::isnan(studyDuration) && milestone >= studyDuration)
    throw std::invalid_argument(
        "milestone cannot exceed studyDuration for fixed follow-up");

  // --- Efficacy boundaries ---------------------------------------------------
  std::vector<double> l(K, -6.0), zero(K, 0.0);
  std::vector<double> critValues = criticalValues;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (K > 1 && criticalValues.size() == K) {
      bool hasNaN = false;
      for (size_t i = 0; i < K - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[K-1])) haybittle = true;
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(K);
      for (size_t i = 0; i < K - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[K-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[K-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  }

  ListCpp probs = exitprobcpp(critValues, l, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(K);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha :
    std::round(cumAlphaSpent.back() * 1e6) / 1e6;

  // --- Futility boundaries ---------------------------------------------------
  std::vector<double> futBounds = futilityBounds;
  if (K > 1) {
    if (missingFutilityBounds && bsf == "none") {
      futBounds = std::vector<double>(K, -6.0);
      futBounds[K-1] = critValues[K-1];
    } else if (!missingFutilityBounds && futBounds.size() == K-1) {
      futBounds.push_back(critValues[K-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  // Randomization probability for treatment group
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // ---- study duration ----
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || std::isnan(studyDuration)) {
    studyDuration1 = accrualDuration + followupTime;
  }

  std::vector<double> time(K), nsubjects(K), nsubjects1(K), nsubjects2(K);
  std::vector<double> nevents(K), nevents1(K), nevents2(K);
  std::vector<double> ndropouts(K), ndropouts1(K), ndropouts2(K);
  std::vector<double> nmilestone(K), nmilestone1(K), nmilestone2(K);
  std::vector<double> info(K);

  // ---- compute maxInformation and theta via kmstat at study end ----
  DataFrameCpp km_end = kmstat1cpp(
    studyDuration1, milestone, allocationRatioPlanned,
    accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction,
    lambda1x, lambda2x, gamma1x, gamma2x,
    accrualDuration, followupTime, fixedFollowup);

  const auto& s1    = km_end.get<double>("surv1");
  const auto& s2    = km_end.get<double>("surv2");
  const auto& vs1   = km_end.get<double>("vsurv1");
  const auto& vs2   = km_end.get<double>("vsurv2");
  double surv1 = 0.0, surv2 = 0.0, vsurv1 = 0.0, vsurv2 = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    surv1  += stratumFraction[h] * s1[h];
    surv2  += stratumFraction[h] * s2[h];
    vsurv1 += stratumFraction[h] * stratumFraction[h] * vs1[h];
    vsurv2 += stratumFraction[h] * stratumFraction[h] * vs2[h];
  }
  double survDiff = surv1 - surv2;
  double vsurvDiff = vsurv1 + vsurv2;
  double maxInformation = 1.0 / vsurvDiff;

  std::vector<double> theta(kMax, survDiff - survDiffH0);

  info[K - 1]        = maxInformation;
  time[K - 1]        = studyDuration1;
  nsubjects[K - 1]   = extract_km(km_end, "subjects");
  nsubjects1[K - 1]  = phi * nsubjects[K - 1];
  nsubjects2[K - 1]  = (1.0 - phi) * nsubjects[K - 1];
  nevents[K - 1]     = extract_km(km_end, "nevents");
  nevents1[K - 1]    = extract_km(km_end, "nevents1");
  nevents2[K - 1]    = extract_km(km_end, "nevents2");
  ndropouts[K - 1]   = extract_km(km_end, "ndropouts");
  ndropouts1[K - 1]  = extract_km(km_end, "ndropouts1");
  ndropouts2[K - 1]  = extract_km(km_end, "ndropouts2");
  nmilestone[K - 1]  = extract_km(km_end, "nmilestone");
  nmilestone1[K - 1] = extract_km(km_end, "nmilestone1");
  nmilestone2[K - 1] = extract_km(km_end, "nmilestone2");

  // ---- compute info, time, and other stats at interim analyses ----
  for (size_t i = 0; i < K - 1; ++i) {
    double information1 = maxInformation * infoRates[i];
    info[i] = information1;

    // solve for analysis time where total information equals information1
    auto g = [&](double t)->double {
      DataFrameCpp km1 = kmstat1cpp(
        t, milestone, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrualDuration, followupTime, fixedFollowup);

      const auto& vs1 = km1.get<double>("vsurv1");
      const auto& vs2 = km1.get<double>("vsurv2");
      double vsurv1 = 0.0, vsurv2 = 0.0;
      for (size_t h = 0; h < nstrata; ++h) {
        vsurv1 += stratumFraction[h] * stratumFraction[i] * vs1[h];
        vsurv2 += stratumFraction[h] * stratumFraction[i] * vs2[h];
      }
      double vsurvDiff = vsurv1 + vsurv2;
      double information = 1.0 / vsurvDiff;
      return information - information1;
    };

    time[i] = brent(g, milestone + 0.001, studyDuration1, 1e-6);

    DataFrameCpp km_i = kmstat1cpp(
      time[i], milestone, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup);

    nsubjects[i]   = extract_km(km_i, "subjects");
    nsubjects1[i]  = phi * nsubjects[i];
    nsubjects2[i]  = (1.0 - phi) * nsubjects[i];
    nevents[i]     = extract_km(km_i, "nevents");
    nevents1[i]    = extract_km(km_i, "nevents1");
    nevents2[i]    = extract_km(km_i, "nevents2");
    ndropouts[i]   = extract_km(km_i, "ndropouts");
    ndropouts1[i]  = extract_km(km_i, "ndropouts1");
    ndropouts2[i]  = extract_km(km_i, "ndropouts2");
    nmilestone[i]  = extract_km(km_i, "nmilestone");
    nmilestone1[i] = extract_km(km_i, "nmilestone");
    nmilestone2[i] = extract_km(km_i, "nmilestone");
  }

  // ---- compute exit probabilities ----
  ListCpp exit_probs;
  if (!missingFutilityBounds || bsf == "none" || K == 1) {
    exit_probs = exitprobcpp(critValues, futBounds, theta, info);
  } else {
    std::vector<double> w(K, 1.0);
    auto gp = getPower(alpha1, kMax, critValues, theta, info, bsf,
                       parameterBetaSpending, spendTime, futStopping, w);
    futBounds = gp.get<std::vector<double>>("futilityBounds");
    exit_probs = gp.get_list("probs");
  }

  std::vector<double> efficacyP(K), futilityP(K);
  for (size_t i = 0; i < K; ++i) {
    efficacyP[i] = 1 - boost_pnorm(critValues[i]);
    futilityP[i] = 1 - boost_pnorm(futBounds[i]);
  }

  auto pu = exit_probs.get<std::vector<double>>("exitProbUpper");
  auto pl = exit_probs.get<std::vector<double>>("exitProbLower");
  std::vector<double> ptotal(K);
  for (size_t i = 0; i < K; ++i) ptotal[i] = pu[i] + pl[i];

  // overall summaries
  double expectedNumberOfEvents = 0.0;
  double expectedNumberOfDropouts = 0.0;
  double expectedNumberOfSubjects = 0.0;
  double expectedNumberOfMiles = 0.0;
  double expectedNumberOfEvents1 = 0.0;
  double expectedNumberOfDropouts1 = 0.0;
  double expectedNumberOfSubjects1 = 0.0;
  double expectedNumberOfMiles1 = 0.0;
  double expectedNumberOfEvents2 = 0.0;
  double expectedNumberOfDropouts2 = 0.0;
  double expectedNumberOfSubjects2 = 0.0;
  double expectedNumberOfMiles2 = 0.0;
  double expectedStudyDuration = 0.0;
  double expectedInformation = 0.0;
  for (size_t i = 0; i < K; ++i) {
    expectedNumberOfEvents += ptotal[i] * nevents[i];
    expectedNumberOfDropouts += ptotal[i] * ndropouts[i];
    expectedNumberOfSubjects += ptotal[i] * nsubjects[i];
    expectedNumberOfMiles += ptotal[i] * nmilestone[i];
    expectedNumberOfEvents1 += ptotal[i] * nevents1[i];
    expectedNumberOfDropouts1 += ptotal[i] * ndropouts1[i];
    expectedNumberOfSubjects1 += ptotal[i] * nsubjects1[i];
    expectedNumberOfMiles1 += ptotal[i] * nmilestone1[i];
    expectedNumberOfEvents2 += ptotal[i] * nevents2[i];
    expectedNumberOfDropouts2 += ptotal[i] * ndropouts2[i];
    expectedNumberOfSubjects2 += ptotal[i] * nsubjects2[i];
    expectedNumberOfMiles2 += ptotal[i] * nmilestone2[i];
    expectedStudyDuration += ptotal[i] * time[i];
    expectedInformation += ptotal[i] * info[i];
  }

  std::vector<double> cpu(K), cpl(K);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());
  double overallReject = cpu[K-1];

  // efficacy/futility implied survDiff boundaries at each stage
  std::vector<double> sdu(K), sdl(K);
  for (size_t i = 0; i < K; ++i) {
    sdu[i] = survDiffH0 + critValues[i] / std::sqrt(info[i]);
    sdl[i] = survDiffH0 + futBounds[i] / std::sqrt(info[i]);
    if (critValues[i] == 6.0) { sdu[i] = NaN; effStopping[i] = 0; }
    if (futBounds[i] == -6.0) { sdl[i] = NaN; futStopping[i] = 0; }
  }

  // --- Build output DataFrames and Lists ------------------------------------
  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha1, "alpha");
  overallResults.push_back(nevents.back(), "numberOfEvents");
  overallResults.push_back(ndropouts.back(), "numberOfDropouts");
  overallResults.push_back(nsubjects.back(), "numberOfSubjects");
  overallResults.push_back(nmilestone.back(), "numberOfMilestone");
  overallResults.push_back(time.back(), "studyDuration");
  overallResults.push_back(maxInformation, "information");
  overallResults.push_back(expectedNumberOfEvents, "expectedNumberOfEvents");
  overallResults.push_back(expectedNumberOfDropouts, "expectedNumberOfDropouts");
  overallResults.push_back(expectedNumberOfSubjects, "expectedNumberOfSubjects");
  overallResults.push_back(expectedNumberOfMiles, "expectedNumberOfMilestone");
  overallResults.push_back(expectedStudyDuration, "expectedStudyDuration");
  overallResults.push_back(expectedInformation, "expectedInformation");
  overallResults.push_back(accrualDuration, "accrualDuration");
  overallResults.push_back(followupTime, "followupTime");
  overallResults.push_back(fixedFollowup, "fixedFollowup");
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(milestone, "milestone");
  overallResults.push_back(survDiffH0, "survDiffH0");
  overallResults.push_back(surv1, "surv1");
  overallResults.push_back(surv2, "surv2");
  overallResults.push_back(survDiff, "survDiff");

  DataFrameCpp byStageResults;
  byStageResults.push_back(std::move(infoRates), "informationRates");
  byStageResults.push_back(std::move(critValues), "efficacyBounds");
  byStageResults.push_back(std::move(futBounds), "futilityBounds");
  byStageResults.push_back(std::move(pu), "rejectPerStage");
  byStageResults.push_back(std::move(pl), "futilityPerStage");
  byStageResults.push_back(std::move(cpu), "cumulativeRejection");
  byStageResults.push_back(std::move(cpl), "cumulativeFutility");
  byStageResults.push_back(std::move(cumAlphaSpent), "cumulativeAlphaSpent");
  byStageResults.push_back(std::move(nevents), "numberOfEvents");
  byStageResults.push_back(std::move(ndropouts), "numberOfDropouts");
  byStageResults.push_back(std::move(nsubjects), "numberOfSubjects");
  byStageResults.push_back(std::move(nmilestone), "numberOfMilestone");
  byStageResults.push_back(std::move(time), "analysisTime");
  byStageResults.push_back(std::move(sdu), "efficacySurvDiff");
  byStageResults.push_back(std::move(sdl), "futilitySurvDiff");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(futilityP), "futilityP");
  byStageResults.push_back(std::move(info), "information");
  byStageResults.push_back(std::move(effStopping), "efficacyStopping");
  byStageResults.push_back(std::move(futStopping), "futilityStopping");

  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
  settings.push_back(typeBetaSpending, "typeBetaSpending");
  settings.push_back(parameterBetaSpending, "parameterBetaSpending");
  settings.push_back(allocationRatioPlanned, "allocationRatioPlanned");
  settings.push_back(accrualTime, "accrualTime");
  settings.push_back(accrualIntensity, "accrualIntensity");
  settings.push_back(piecewiseSurvivalTime, "piecewiseSurvivalTime");
  settings.push_back(stratumFraction, "stratumFraction");
  settings.push_back(lambda1, "lambda1");
  settings.push_back(lambda2, "lambda2");
  settings.push_back(gamma1, "gamma1");
  settings.push_back(gamma2, "gamma2");
  settings.push_back(spendingTime, "spendingTime");

  ListCpp byTreatmentCounts;
  byTreatmentCounts.push_back(std::move(nevents1), "numberOfEvents1");
  byTreatmentCounts.push_back(std::move(ndropouts1), "numberOfDropouts1");
  byTreatmentCounts.push_back(std::move(nsubjects1), "numberOfSubjects1");
  byTreatmentCounts.push_back(std::move(nmilestone1), "numberOfMilestone1");
  byTreatmentCounts.push_back(std::move(nevents2), "numberOfEvents2");
  byTreatmentCounts.push_back(std::move(ndropouts2), "numberOfDropouts2");
  byTreatmentCounts.push_back(std::move(nsubjects2), "numberOfSubjects2");
  byTreatmentCounts.push_back(std::move(nmilestone2), "numberOfMilestone2");
  byTreatmentCounts.push_back(expectedNumberOfEvents1, "expectedNumberOfEvents1");
  byTreatmentCounts.push_back(expectedNumberOfDropouts1, "expectedNumberOfDropouts1");
  byTreatmentCounts.push_back(expectedNumberOfSubjects1, "expectedNumberOfSubjects1");
  byTreatmentCounts.push_back(expectedNumberOfMiles1, "expectedNumberOfMilestone1");
  byTreatmentCounts.push_back(expectedNumberOfEvents2, "expectedNumberOfEvents2");
  byTreatmentCounts.push_back(expectedNumberOfDropouts2, "expectedNumberOfDropouts2");
  byTreatmentCounts.push_back(expectedNumberOfSubjects2, "expectedNumberOfSubjects2");
  byTreatmentCounts.push_back(expectedNumberOfMiles2, "expectedNumberOfMilestone2");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(settings), "settings");
  result.push_back(std::move(byTreatmentCounts), "byTreatmentCounts");

  return result;
}


//' @title Power for Difference in Milestone Survival Probabilities
//' @description Estimates the power for testing the difference in
//' milestone survival probabilities in a two-sample survival design.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @param typeBetaSpending The type of beta spending. One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock
//'   type spending function, "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and "none" for no
//'   early futility stopping. Defaults to "none".
//' @inheritParams param_parameterBetaSpending
//' @param milestone The milestone time at which to calculate the survival
//'   probability.
//' @param survDiffH0 The difference in milestone survival probabilities
//'   under the null hypothesis. Defaults to 0 for superiority test.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//'
//' @return An S3 class \code{kmpower} object with 4 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numberOfDropouts}: The total number of dropouts.
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{numberOfMilestone}: The total number of subjects reaching
//'       milestone.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
//'
//'     - \code{expectedNumberOfDropouts}: The expected number of dropouts.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedNumberOfMilestone}: The expected number of subjects
//'       reaching milestone.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up duration.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{milestone}: The milestone time relative to randomization.
//'
//'     - \code{survDiffH0}: The difference in milestone survival
//'       probabilities under the null hypothesis.
//'
//'     - \code{surv1}: The milestone survival probability for the
//'       treatment group.
//'
//'     - \code{surv2}: The milestone survival probability for the
//'       control group.
//'
//'     - \code{survDiff}: The difference in milestone survival
//'       probabilities, equal to \code{surv1 - surv2}.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
//'
//'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{futilityPerStage}: The probability for futility stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeFutility}: The cumulative probability for futility
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
//'
//'     - \code{numberOfEvents}: The number of events.
//'
//'     - \code{numberOfDropouts}: The number of dropouts.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{numberOfMilestone}: The number of subjects reaching
//'       milestone.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacySurvDiff}: The efficacy boundaries on the survival
//'       difference scale.
//'
//'     - \code{futilitySurvDiff}: The futility boundaries on the survival
//'       difference scale.
//'
//'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'
//'     - \code{futilityP}: The futility boundaries on the p-value scale.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'
//'     - \code{futilityStopping}: Whether to allow futility stopping.
//'
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{typeBetaSpending},
//'   \code{parameterBetaSpending}, \code{allocationRatioPlanned},
//'   \code{accrualTime}, \code{accuralIntensity},
//'   \code{piecewiseSurvivalTime}, \code{stratumFraction},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   and \code{spendingTime}.
//'
//' * \code{byTreatmentCounts}: A list containing the following counts by
//'   treatment group:
//'
//'     - \code{numberOfEvents1}: The number of events by stage for
//'       the treatment group.
//'
//'     - \code{numberOfDropouts1}: The number of dropouts by stage for
//'       the treatment group.
//'
//'     - \code{numberOfSubjects1}: The number of subjects by stage for
//'       the treatment group.
//'
//'     - \code{numberOfMilestone1}: The number of subjects reaching
//'       milestone by stage for the active treatment group.
//'
//'     - \code{numberOfEvents2}: The number of events by stage for
//'       the control group.
//'
//'     - \code{numberOfDropouts2}: The number of dropouts by stage for
//'       the control group.
//'
//'     - \code{numberOfSubjects2}: The number of subjects by stage for
//'       the control group.
//'
//'     - \code{numberOfMilestone2}: The number of subjects reaching
//'       milestone by stage for the control group.
//'
//'     - \code{expectedNumberOfEvents1}: The expected number of events for
//'       the treatment group.
//'
//'     - \code{expectedNumberOfDropouts1}: The expected number of dropouts
//'       for the active treatment group.
//'
//'     - \code{expectedNumberOfSubjects1}: The expected number of subjects
//'       for the active treatment group.
//'
//'     - \code{expectedNumberOfMilestone1}: The expected number of subjects
//'       reaching milestone for the active treatment group.
//'
//'     - \code{expectedNumberOfEvents2}: The expected number of events for
//'       control group.
//'
//'     - \code{expectedNumberOfDropouts2}: The expected number of dropouts
//'       for the control group.
//'
//'     - \code{expectedNumberOfSubjects2}: The expected number of subjects
//'       for the control group.
//'
//'     - \code{expectedNumberOfMilestone2}: The expected number of subjects
//'       reaching milestone for the control group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survival, and 5% dropout by
//' # the end of 1 year.
//'
//' kmpower(kMax = 2, informationRates = c(0.8, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         milestone = 18,
//'         allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'         accrualIntensity = 26/9*seq(1, 9),
//'         piecewiseSurvivalTime = c(0, 6),
//'         stratumFraction = c(0.2, 0.8),
//'         lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'         lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'         gamma1 = -log(1-0.05)/12,
//'         gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'         followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List kmpower(
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& futilityBounds = NA_REAL,
    const std::string& typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const double milestone = NA_REAL,
    const double survDiffH0 = 0,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double accrualDuration = NA_REAL,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = false,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const double studyDuration = NA_REAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto piecewiseTime = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  ListCpp out = kmpowercpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, survDiffH0, allocationRatioPlanned,
    accrualT, accrualInt, piecewiseTime, stratumFrac,
    lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, studyDuration);

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "kmpower";
  return result;
}


ListCpp kmsamplesizecpp(
    const double beta,
    const int kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<unsigned char>& futilityStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& futilityBounds,
    const std::string& typeBetaSpending,
    const double parameterBetaSpending,
    const std::vector<double>& userBetaSpending,
    const double milestone,
    const double survDiffH0,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    double accrualDuration, // may be solved
    double followupTime,    // may be solved
    const bool fixedFollowup,
    const std::vector<double>& spendingTime,
    const bool rounding) {

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1))
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  if (beta < 0.0001 || (!std::isnan(alpha) && beta >= 1.0 - alpha))
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");

  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");
  size_t K = static_cast<size_t>(kMax);

  // informationRates: default to (1:kMax)/kMax if missing
  std::vector<double> infoRates(K);
  if (none_na(informationRates)) {
    if (informationRates.size() != K)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[K-1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else {
    for (size_t i = 0; i < K; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(K);
  }

  // effStopping: default to all 1s if missing
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != K)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[K-1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(K, 1);
  }

  // futStopping: default to all 1s if missing
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != K)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping[K-1] != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(K, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds);

  if (!missingCriticalValues && criticalValues.size() != K) {
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
    if (userAlphaSpending.size() != K)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending[K-1] != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  if (!missingFutilityBounds) {
    if (!(futilityBounds.size() == K - 1 || futilityBounds.size() == K)) {
      throw std::invalid_argument("Invalid length for futilityBounds");
    }
  }
  if (!missingCriticalValues && !missingFutilityBounds) {
    for (size_t i = 0; i < K - 1; ++i) {
      if (futilityBounds[i] > criticalValues[i]) {
        throw std::invalid_argument("futilityBounds must lie below criticalValues");
      }
    }
    if (futilityBounds.size() == K && futilityBounds[K-1] != criticalValues[K-1]) {
      throw std::invalid_argument(
          "futilityBounds must meet criticalValues at the final look");
    }
  }

  std::string bsf = typeBetaSpending;
  for (char &c : bsf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
      bsf == "sfkd" || bsf == "sfhsd" || bsf == "user" || bsf == "none")) {
    throw std::invalid_argument("Invalid value for typeBetaSpending");
  }
  if ((bsf == "sfkd" || bsf == "sfhsd") && std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument("parameterBetaSpending must be positive for sfKD");
  }

  if (missingFutilityBounds && bsf=="user") {
    if (!none_na(userBetaSpending))
      throw std::invalid_argument("userBetaSpending must be specified");
    if (userBetaSpending.size() != K)
      throw std::invalid_argument("Invalid length of userBetaSpending");
    if (userBetaSpending[0] < 0.0)
      throw std::invalid_argument("userBetaSpending must be nonnegative");
    if (any_nonincreasing(userBetaSpending))
      throw std::invalid_argument("userBetaSpending must be nondecreasing");
    if (userBetaSpending[K-1] != beta)
      throw std::invalid_argument("userBetaSpending must end with specified beta");
  }

  // validate milestone and H0 diff
  if (std::isnan(milestone)) throw std::invalid_argument("milestone must be provided");
  if (milestone <= 0.0) throw std::invalid_argument("milestone must be positive");
  if (survDiffH0 <= -1.0 || survDiffH0 >= 1.0)
    throw std::invalid_argument("survDiffH0 must lie between -1 and 1");

  // validate accrual + survival inputs (same style as lrstat.cpp)
  if (allocationRatioPlanned <= 0.0)
    throw std::invalid_argument("allocationRatioPlanned must be positive");

  if (accrualTime[0] != 0.0)
    throw std::invalid_argument("accrualTime must start with 0");
  if (any_nonincreasing(accrualTime))
    throw std::invalid_argument("accrualTime should be increasing");

  if (!none_na(accrualIntensity))
    throw std::invalid_argument("accrualIntensity must be provided");
  if (accrualIntensity.size() != accrualTime.size())
    throw std::invalid_argument("Invalid length for accrualIntensity");
  for (double v : accrualIntensity) {
    if (v < 0.0) throw std::invalid_argument("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0.0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");

  for (double v : stratumFraction) {
    if (v <= 0.0) throw std::invalid_argument("stratumFraction must be positive");
  }
  double sumf = std::accumulate(stratumFraction.begin(), stratumFraction.end(), 0.0);
  if (std::fabs(sumf - 1.0) > 1e-12)
    throw std::invalid_argument("stratumFraction must sum to 1");

  if (!none_na(lambda1)) throw std::invalid_argument("lambda1 must be provided");
  if (!none_na(lambda2)) throw std::invalid_argument("lambda2 must be provided");

  for (double v : lambda1) {
    if (v < 0.0) throw std::invalid_argument("lambda1 must be non-negative");
  }
  for (double v : lambda2) {
    if (v < 0.0) throw std::invalid_argument("lambda2 must be non-negative");
  }
  for (double v : gamma1) {
    if (v < 0.0) throw std::invalid_argument("gamma1 must be non-negative");
  }
  for (double v : gamma2) {
    if (v < 0.0) throw std::invalid_argument("gamma2 must be non-negative");
  }

  // expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambda1x = expand_stratified(lambda1, nstrata, nintv, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintv, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintv, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintv, "gamma2");

  if (!std::isnan(accrualDuration) && accrualDuration <= 0.0)
    throw std::invalid_argument("accrualDuration must be positive");

  if (!std::isnan(followupTime) && fixedFollowup && followupTime <= 0.0)
    throw std::invalid_argument("followupTime must be positive for fixed follow-up");
  if (!std::isnan(followupTime) && !fixedFollowup && followupTime < 0.0)
    throw std::invalid_argument(
        "followupTime must be non-negative for variable follow-up");
  if (fixedFollowup && std::isnan(followupTime))
    throw std::invalid_argument("followupTime must be provided for fixed follow-up");

  if (!std::isnan(accrualDuration) && !std::isnan(followupTime) &&
      (milestone >= accrualDuration + followupTime)) {
    throw std::invalid_argument(
        "milestone must be less than accrualDuration + followupTime");
  }
  if (fixedFollowup && milestone > followupTime)
    throw std::invalid_argument(
        "milestone cannot exceed followupTime for fixed follow-up");

  // spendingTime default to informationRates
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != K)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[K-1] != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }


  // --- Efficacy boundaries ---------------------------------------------------
  std::vector<double> l(K, -6.0), zero(K, 0.0);
  std::vector<double> critValues = criticalValues;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (K > 1 && criticalValues.size() == K) {
      bool hasNaN = false;
      for (size_t i = 0; i < K - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[K-1])) haybittle = true;
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(K);
      for (size_t i = 0; i < K - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[K-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[K-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  }

  ListCpp probs = exitprobcpp(critValues, l, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(K);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha :
    std::round(cumAlphaSpent.back() * 1e6) / 1e6;


  // --- Futility boundaries ---------------------------------------------------
  std::vector<double> futBounds = futilityBounds;
  if (K > 1) {
    if (missingFutilityBounds && bsf == "none") {
      futBounds = std::vector<double>(K, -6.0);
      futBounds[K-1] = critValues[K-1];
    } else if (!missingFutilityBounds && futBounds.size() == K-1) {
      futBounds.push_back(critValues[K-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }


  // Which design parameter is unknown?
  enum Unknown { ACC_DUR, FUP_TIME, ACC_INT };
  Unknown unknown;
  bool missAccrual = std::isnan(accrualDuration);
  bool missFollow = std::isnan(followupTime);
  if (missAccrual && !missFollow) unknown = ACC_DUR;
  else if (!missAccrual && missFollow) unknown = FUP_TIME;
  else if (!missAccrual && !missFollow) unknown = ACC_INT;
  else throw std::invalid_argument(
      "accrualDuration and followupTime cannot be both missing");

  // local accrualIntensity copy (we might scale it)
  std::vector<double> accrualInt = accrualIntensity;

  // milestone survival under H1
  std::vector<double> zerogam(nintv, 0.0);
  double surv1 = 0.0, surv2 = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    double p1 = patrisk1(milestone, piecewiseSurvivalTime, lambda1x[h], zerogam);
    double p2 = patrisk1(milestone, piecewiseSurvivalTime, lambda2x[h], zerogam);
    surv1 += stratumFraction[h] * p1;
    surv2 += stratumFraction[h] * p2;
  }
  double theta1 = (surv1 - surv2 - survDiffH0);
  std::vector<double> theta(K, theta1);

  // --- determine target maxInformation from group sequential design ---
  ListCpp design = getDesigncpp(
    beta, NaN, theta1, kMax, infoRates,
    effStopping, futStopping,
    critValues, alpha1, asf,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, bsf, parameterBetaSpending,
    userBetaSpending, spendTime, 1.0);

  auto byStageResults = design.get<DataFrameCpp>("byStageResults");
  futBounds = byStageResults.get<double>("futilityBounds");

  auto overallResults = design.get<DataFrameCpp>("overallResults");
  double maxInformation = overallResults.get<double>("information")[0];

  // Helper: compute information under H1 given accrualDuration (accrDur),
  // followupTime (fu), and accrualIntensity (accrInt).
  auto info_minus_target_H1 = [&](double t, double accrDur, double fu,
                                  const std::vector<double>& accrInt) {
    DataFrameCpp km1 = kmstat1cpp(
      t, milestone, allocationRatioPlanned,
      accrualTime, accrInt,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrDur, fu, fixedFollowup);

    const auto& vs1 = km1.get<double>("vsurv1");
    const auto& vs2 = km1.get<double>("vsurv2");
    double vsurv1 = 0.0, vsurv2 = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      vsurv1 += stratumFraction[h] * stratumFraction[h] * vs1[h];
      vsurv2 += stratumFraction[h] * stratumFraction[h] * vs2[h];
    }
    double vsurvDiff = vsurv1 + vsurv2;
    return 1.0 / vsurvDiff - maxInformation;
  };


  // when info at zero follow-up (end of enrollment) already exceeds target,
  // set followupTime = 0 and find minimal accrualDuration achieving target
  bool curtailed = false;

  // when info at maximum follow-up (milestone) is still below target,
  // increase accrualDuration to achieve target
  bool expanded = false;

  // Root-finding function to solve for the unknown design parameter to achieve target
  std::function<double(double)> f_root;

  if (unknown == ACC_DUR) {
    f_root = [&](double x)->double {
      return info_minus_target_H1(x + followupTime, x, followupTime, accrualInt);
    };
  } else if (unknown == FUP_TIME) {
    if (!fixedFollowup && accrualDuration > milestone + 0.001 &&
        info_minus_target_H1(accrualDuration, accrualDuration, 0.0, accrualInt) > 0) {
      thread_utils::push_thread_warning(
        "Information at zero follow-up (end of enrollment) already exceeds target. "
        "Setting followupTime = 0 and finding minimal accrualDuration.");

      f_root = [&](double x)->double {
        return info_minus_target_H1(x, x, 0.0, accrualInt);
      };

      curtailed = true;
    } else if (info_minus_target_H1(accrualDuration + milestone, accrualDuration,
                                    milestone, accrualInt) < 0) {
      thread_utils::push_thread_warning(
        "The required information cannot be attained by increasing "
        "followupTime alone. accrualDuration is also increased to "
        "attain the required information.");

      f_root = [&](double x)->double {
        return info_minus_target_H1(x + milestone, x, milestone, accrualInt);
      };

      expanded = true;
    } else {
      f_root = [&](double x)->double {
        return info_minus_target_H1(accrualDuration + x, accrualDuration, x,
                                    accrualInt);
      };
    }
  } else {
    f_root = [&](double m)->double {
      std::vector<double> scaled = accrualInt;
      for (double &v : scaled) v *= m;
      return info_minus_target_H1(accrualDuration + followupTime,
                                  accrualDuration, followupTime, scaled);
    };
  }

  double lower, upper;
  if (unknown == ACC_DUR) {
    lower = std::max(milestone - followupTime, 0.0) + 0.001;
    upper = 2.0 * lower;
  } else if (unknown == FUP_TIME && curtailed) {
    lower = milestone + 0.001;
    upper = accrualDuration;
  } else if (unknown == FUP_TIME && expanded) {
    lower = accrualDuration;
    upper = 2.0 * lower;
  } else if (unknown == FUP_TIME && !curtailed && !expanded) {
    lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
    upper = milestone;
  } else { // unknown == ACC_INT
    lower = 0.001;
    upper = 120;
  }

  // expand upper if needed to ensure root is bracketed
  double fl_val = f_root(lower), fu_val = f_root(upper);
  if (unknown == ACC_DUR || (unknown == FUP_TIME && expanded) || unknown == ACC_INT) {
    int expand_iter = 0;
    while (fl_val * fu_val > 0.0 && expand_iter < 60) {
      lower = upper;
      fl_val = fu_val;
      upper *= 2.0;
      fu_val = f_root(upper);
      ++expand_iter;
    }
  }
  if (fl_val * fu_val > 0.0) throw std::runtime_error(
      "Unable to bracket root; check interval or inputs");

  // solve for root and apply solution for the unknown design parameter
  auto f_for_brent1 = [&](double x)->double {
    if (x == lower) return fl_val;
    if (x == upper) return fu_val;
    return f_root(x);
  };
  double solution = brent(f_for_brent1, lower, upper, 1e-6);

  if (unknown == ACC_DUR) {
    accrualDuration = solution;
  } else if (unknown == FUP_TIME && curtailed) {
    followupTime = 0.0;
    accrualDuration = solution;
  } else if (unknown == FUP_TIME && expanded) {
    followupTime = milestone;
    accrualDuration = solution;
  } else if (unknown == FUP_TIME && !curtailed && !expanded) {
    followupTime = solution;
  } else { // scaled multiplier for accrualIntensity
    for (double &v : accrualInt) v *= solution;
  }

  double studyDuration = accrualDuration + followupTime;

  // --- rounding to integer N (adjust intensity or accrualDuration) ---
  if (rounding) {
    double n0 = accrual1(studyDuration, accrualTime, accrualInt, accrualDuration);
    double n = std::ceil(n0 - 1.0e-12);

    if (n - n0 > 1e-6) {
      if (unknown == ACC_INT) { // scale intensity to hit integer n
        double mult = n / n0;
        for (double& v : accrualInt) v *= mult;
      } else { // adjust accrualDuration to hit integer n
        accrualDuration = getAccrualDurationFromN1(n, accrualTime, accrualInt);
      }

      if (!fixedFollowup) {
        // variable follow-up: adjust follow-up time to match maxInformation
        auto h = [&](double x)->double {
          return info_minus_target_H1(accrualDuration + x, accrualDuration, x,
                                      accrualInt);
        };
        double lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // fixed follow-up: adjust studyDuration by extending post-accrual time
        auto h = [&](double x)->double {
          return info_minus_target_H1(accrualDuration + x, accrualDuration,
                                      followupTime, accrualInt);
        };
        double lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
        double extra = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + extra;
      }
    }
  }

  // --- Results under H1 ---
  ListCpp resultsH1 = kmpowercpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha1, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, survDiffH0, allocationRatioPlanned,
    accrualTime, accrualInt,
    piecewiseSurvivalTime, stratumFraction,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  // refresh maxInformation from H1 power results    {
  DataFrameCpp overallResultsH1 = resultsH1.get<DataFrameCpp>("overallResults");
  maxInformation = overallResultsH1.get<double>("information")[0];


  // --- Construct lambda1 under H0 (scale lambda2 to match survDiffH0) ---
  // Solve for multiplier a s.t. weighted_surv1(a*lambda2) - surv2 - survDiffH0 == 0
  auto f_surv = [&](double aval)->double {
    double surv1 = 0.0;
    for (std::size_t h = 0; h < nstrata; ++h) {
      auto lam1H0 = lambda2x[h]; // build per-stratum scaled lambda
      for (double &v : lam1H0) v *= aval;
      double p = patrisk1(milestone, piecewiseSurvivalTime, lam1H0, zerogam);
      surv1 += stratumFraction[h] * p;
    }
    return surv1 - surv2 - survDiffH0;
  };

  double lo = 0.001, hi = 2.0;
  double flo = f_surv(lo);
  double fhi = f_surv(hi);
  int it = 0;
  while (flo * fhi > 0.0 && it < 60) {
    lo = hi;
    flo = fhi;
    hi *= 2.0;
    fhi = f_surv(hi);
  }
  if (flo * fhi > 0.0) throw std::runtime_error(
      "Unable to bracket root for lambda1H0 under H0; check inputs");

  auto f_for_brent2 = [&](double x)->double {
    if (x == lo) return flo;
    if (x == hi) return fhi;
    return f_surv(x);
  };

  double multH0 = brent(f_for_brent2, lo, hi, 1.0e-6);

  std::vector<double> lambda1H0 = lambda2;
  for (double &v : lambda1H0) v *= multH0;
  auto lambda1H0x = expand_stratified(lambda1H0, nstrata, nintv, "lambda1H0");

  // Helper: compute information under H0 given accrualDuration (accrDur),
  // followupTime (fu), and accrualIntensity (accrInt).
  auto info_minus_target_H0 = [&](double t, double accrDur, double fu,
                                  const std::vector<double>& accrInt) {
    DataFrameCpp km1 = kmstat1cpp(
      t, milestone, allocationRatioPlanned,
      accrualTime, accrInt,
      piecewiseSurvivalTime, stratumFraction,
      lambda1H0x, lambda2x, gamma1x, gamma2x,
      accrDur, fu, fixedFollowup);

    const auto& vs1 = km1.get<double>("vsurv1");
    const auto& vs2 = km1.get<double>("vsurv2");
    double vsurv1 = 0.0, vsurv2 = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      vsurv1 += stratumFraction[h] * stratumFraction[h] * vs1[h];
      vsurv2 += stratumFraction[h] * stratumFraction[h] * vs2[h];
    }
    double vsurvDiff = vsurv1 + vsurv2;
    return 1.0 / vsurvDiff - maxInformation;
  };

  if (!fixedFollowup) {
    auto h_follow = [&](double x)->double {
      return info_minus_target_H0(accrualDuration + x, accrualDuration, x, accrualInt);
    };

    if (accrualDuration > milestone + 0.001 && h_follow(0.0) > 0.0) {
      thread_utils::push_thread_warning(
        "Information at zero follow-up (end of enrollment) already exceeds target. "
        "Setting followupTime = 0 and finding minimal accrualDuration.");

      auto h_accr = [&](double x)->double {
        return info_minus_target_H0(x, x, 0.0, accrualInt);
      };

      double lo = milestone + 0.001;
      double hi = accrualDuration;
      accrualDuration = brent(h_accr, lo, hi, 1.0e-6);
      followupTime = 0.0;
      studyDuration = accrualDuration;
    } else if (h_follow(followupTime) < 0.0) {
      thread_utils::push_thread_warning(
        "The required information cannot be attained by increasing "
        "followupTime alone. accrualDuration is also increased to "
        "attain the required information.");

      auto h_accr = [&](double x)->double {
        return info_minus_target_H0(x + milestone, x, milestone, accrualInt);
      };

      // adjust accrualDuration upward to match target
      double lo = accrualDuration;
      double hi = 2.0 * accrualDuration;
      double flo = h_accr(lo);
      double fhi = h_accr(hi);
      int it = 0;
      while (flo * fhi > 0.0 && it < 60) {
        lo = hi;
        flo = fhi;
        hi *= 2.0;
        fhi = h_accr(hi);
      }
      if (flo * fhi > 0.0) throw std::runtime_error(
          "Unable to bracket accrualDuration under H0; check inputs");

      auto h_for_brent = [&](double x)->double {
        if (x == lo) return flo;
        if (x == hi) return fhi;
        return h_accr(x);
      };

      accrualDuration = brent(h_for_brent, lo, hi, 1.0e-6);
      followupTime = milestone;
      studyDuration = accrualDuration + followupTime;
    } else {
      // adjust follow-up time to match target
      double lo = std::max(milestone - accrualDuration, 0.0) + 0.001;
      double hi = milestone;

      followupTime = brent(h_follow, lo, hi, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    }
  } else {
    auto h_study = [&](double x)->double {
      return info_minus_target_H0(accrualDuration + x, accrualDuration,
                                  followupTime, accrualInt);
    };

    auto h_accr = [&](double x)->double {
      return info_minus_target_H0(x + followupTime, x, followupTime, accrualInt);
    };

    if (accrualDuration > milestone + 0.001 && h_study(0.0) > 0.0) {
      thread_utils::push_thread_warning(
        "Information at zero follow-up (end of enrollment) already exceeds target. "
        "Decrease accrual druation.");

      double lo = milestone + 0.001;
      double hi = accrualDuration;
      accrualDuration = brent(h_accr, lo, hi, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (h_study(followupTime) < 0.0) {
      thread_utils::push_thread_warning(
        "The required information cannot be attained by increasing "
        "followupTime alone. accrualDuration is also increased to "
        "attain the required information.");

      double lo = accrualDuration;
      double hi = 2.0 * accrualDuration;
      double flo = h_accr(lo);
      double fhi = h_accr(hi);
      int it = 0;
      while (flo * fhi > 0.0 && it < 60) {
        lo = hi;
        flo = fhi;
        hi *= 2.0;
        fhi = h_accr(hi);
      }
      if (flo * fhi > 0.0) throw std::runtime_error(
          "Unable to bracket accrualDuration under H0; check inputs");

      auto h_for_brent = [&](double x)->double {
        if (x == lo) return flo;
        if (x == hi) return fhi;
        return h_accr(x);
      };

      accrualDuration = brent(h_for_brent, lo, hi, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else {
      double lo = std::max(milestone - accrualDuration, 0.0) + 0.001;
      double hi = milestone;
      double extra = brent(h_study, lo, hi, 1.0e-6);
      studyDuration = accrualDuration + extra;
    }
  }

  // --- Results under H0 with same boundaries as H1 ---
  ListCpp resultsH0 = kmpowercpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha1, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, survDiffH0, allocationRatioPlanned,
    accrualTime, accrualInt,
    piecewiseSurvivalTime, stratumFraction,
    lambda1H0, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  // add userBetaSpending to settings
  if (bsf == "user") {
    ListCpp& settingsH1 = resultsH1.get_list("settings");
    ListCpp& settingsH0 = resultsH0.get_list("settings");
    settingsH1.push_back(userBetaSpending, "userBetaSpending");
    settingsH0.push_back(userBetaSpending, "userBetaSpending");
  }

  ListCpp result;
  result.push_back(std::move(resultsH1), "resultsUnderH1");
  result.push_back(std::move(resultsH0), "resultsUnderH0");

  return result;
}


//' @title Sample Size for Difference in Milestone Survival Probabilities
//' @description Obtains the needed accrual duration given power,
//' accrual intensity, and follow-up time, the needed follow-up time
//' given power, accrual intensity, and accrual duration, or the needed
//' absolute accrual intensity given power, relative accrual intensity,
//' accrual duration, and follow-up time in a two-group survival design.
//'
//' @param beta Type II error. Defaults to 0.2.
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @inheritParams param_typeBetaSpending
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_userBetaSpending
//' @param milestone The milestone time at which to calculate the survival
//'   probability.
//' @param survDiffH0 The difference in milestone survival probabilities
//'   under the null hypothesis. Defaults to 0 for superiority test.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupTime, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return A list of two components:
//'
//' * \code{resultsUnderH1}: An S3 class \code{kmpower} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{kmpower} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{kmpower}}
//'
//' @examples
//' # Example 1: Obtains follow-up time given power, accrual intensity,
//' # and accrual duration for variable follow-up. Of note, the power
//' # reaches the maximum when the follow-up time equals milestone.
//'
//' kmsamplesize(beta = 0.25, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'              followupTime = NA, fixedFollowup = FALSE)
//'
//' # Example 2: Obtains accrual intensity given power, accrual duration, and
//' # follow-up time for variable follow-up
//'
//' kmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//'
//' # Example 3: Obtains accrual duration given power, accrual intensity, and
//' # follow-up time for fixed follow-up
//'
//' kmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = NA,
//'              followupTime = 18, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List kmsamplesize(
    const double beta = 0.2,
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& futilityBounds = NA_REAL,
    const std::string& typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const Rcpp::NumericVector& userBetaSpending = NA_REAL,
    const double milestone = NA_REAL,
    const double survDiffH0 = 0,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    double accrualDuration = NA_REAL,
    double followupTime = NA_REAL,
    const bool fixedFollowup = false,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool rounding = true) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto userBeta = Rcpp::as<std::vector<double>>(userBetaSpending);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  auto out = kmsamplesizecpp(
    beta, kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    userBeta, milestone, survDiffH0, allocationRatioPlanned,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, rounding);

  thread_utils::drain_thread_warnings_to_R();

  ListCpp resultsUnderH1 = out.get_list("resultsUnderH1");
  ListCpp resultsUnderH0 = out.get_list("resultsUnderH0");

  Rcpp::List resultsH1 = Rcpp::wrap(resultsUnderH1);
  Rcpp::List resultsH0 = Rcpp::wrap(resultsUnderH0);

  resultsH1.attr("class") = "kmpower";
  resultsH0.attr("class") = "kmpower";

  return Rcpp::List::create(
    Rcpp::Named("resultsUnderH1") = resultsH1,
    Rcpp::Named("resultsUnderH0") = resultsH0
  );
}


ListCpp kmpower1scpp(
    const int kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<unsigned char>& futilityStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& futilityBounds,
    const std::string& typeBetaSpending,
    const double parameterBetaSpending,
    const double milestone,
    const double survH0,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda,
    const std::vector<double>& gamma,
    const double accrualDuration,
    const double followupTime,
    const bool fixedFollowup,
    const std::vector<double>& spendingTime,
    const double studyDuration) {

  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");
  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1.0))
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");

  size_t K = static_cast<size_t>(kMax);

  // informationRates default
  std::vector<double> infoRates(K);
  if (none_na(informationRates)) {
    if (informationRates.size() != K)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[K - 1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates;
  } else {
    for (size_t i = 0; i < K; ++i)
      infoRates[i] = static_cast<double>(i + 1) / static_cast<double>(K);
  }

  // efficacyStopping default
  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != K)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[K - 1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping;
  } else {
    effStopping.assign(K, 1);
  }

  // futilityStopping default
  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != K)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping[K - 1] != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping;
  } else {
    futStopping.assign(K, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds);

  if (!missingCriticalValues && criticalValues.size() != K)
    throw std::invalid_argument("Invalid length for criticalValues");
  if (missingCriticalValues && std::isnan(alpha))
    throw std::invalid_argument("alpha must be provided for missing criticalValues");

  // canonicalize alpha-spending name (lowercase)
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
    if (userAlphaSpending.size() != K)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending[K-1] != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  if (!missingFutilityBounds) {
    if (!(futilityBounds.size() == K - 1 || futilityBounds.size() == K)) {
      throw std::invalid_argument("Invalid length for futilityBounds");
    }
  }
  if (!missingCriticalValues && !missingFutilityBounds) {
    for (size_t i = 0; i < K - 1; ++i) {
      if (futilityBounds[i] > criticalValues[i]) {
        throw std::invalid_argument("futilityBounds must lie below criticalValues");
      }
    }
    if (futilityBounds.size() == K && futilityBounds[K-1] != criticalValues[K-1]) {
      throw std::invalid_argument(
          "futilityBounds must meet criticalValues at the final look");
    }
  }

  std::string bsf = typeBetaSpending;
  for (char &c : bsf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
      bsf == "sfkd" || bsf == "sfhsd" || bsf == "none")) {
    throw std::invalid_argument("Invalid value for typeBetaSpending");
  }
  if ((bsf == "sfkd" || bsf == "sfhsd") && std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument("parameterBetaSpending must be positive for sfKD");
  }

  // milestone + survival diff null
  if (std::isnan(milestone)) throw std::invalid_argument("milestone must be provided");
  if (milestone <= 0.0) throw std::invalid_argument("milestone must be positive");
  if (survH0 <= 0.0 || survH0 >= 1.0)
    throw std::invalid_argument("survDiffH0 must lie between 0 and 1");

  if (accrualTime[0] != 0.0)
    throw std::invalid_argument("accrualTime must start with 0");
  if (any_nonincreasing(accrualTime))
    throw std::invalid_argument("accrualTime should be increasing");

  if (!none_na(accrualIntensity))
    throw std::invalid_argument("accrualIntensity must be provided");
  if (accrualIntensity.size() != accrualTime.size())
    throw std::invalid_argument("Invalid length for accrualIntensity");
  for (double v : accrualIntensity) {
    if (v < 0.0) throw std::invalid_argument("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0.0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");

  for (double v : stratumFraction) {
    if (v <= 0.0) throw std::invalid_argument("stratumFraction must be positive");
  }
  double sumf = std::accumulate(stratumFraction.begin(), stratumFraction.end(), 0.0);
  if (std::fabs(sumf - 1.0) > 1e-12)
    throw std::invalid_argument("stratumFraction must sum to 1");

  if (none_na(lambda) == false) throw std::invalid_argument("lambda must be provided");
  for (double v : lambda) {
    if (v < 0.0) throw std::invalid_argument("lambda must be non-negative");
  }
  for (double v : gamma) {
    if (v < 0.0) throw std::invalid_argument("gamma must be non-negative");
  }

  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambdax = expand_stratified(lambda, nstrata, nintv, "lambda");
  auto gammax = expand_stratified(gamma, nstrata, nintv, "gamma");

  if (std::isnan(accrualDuration))
    throw std::invalid_argument("accrualDuration must be provided");
  if (accrualDuration <= 0.0)
    throw std::invalid_argument("accrualDuration must be positive");
  if (std::isnan(followupTime))
    throw std::invalid_argument("followupTime must be provided");
  if (fixedFollowup && followupTime <= 0.0)
    throw std::invalid_argument("followupTime must be positive for fixed follow-up");
  if (!fixedFollowup && followupTime < 0.0)
    throw std::invalid_argument(
        "followupTime must be non-negative for variable follow-up");

  // spendingTime default to informationRates
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != K)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[K-1] != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }

  // fixed follow-up studyDuration constraints
  if (fixedFollowup && !std::isnan(studyDuration) &&
      (studyDuration < accrualDuration))
    throw std::invalid_argument(
        "studyDuration must be greater than or equal to accrualDuration");
  if (fixedFollowup && !std::isnan(studyDuration) &&
      (studyDuration > accrualDuration + followupTime))
    throw std::invalid_argument(
        "studyDuration must be less than or equal to accrualDuration + followupTime");

  if (milestone >= accrualDuration + followupTime)
    throw std::invalid_argument(
        "milestone must be less than accrualDuration + followupTime");
  if (fixedFollowup && milestone > followupTime)
    throw std::invalid_argument(
        "milestone cannot exceed followupTime for fixed follow-up");
  if (fixedFollowup && !std::isnan(studyDuration) && milestone >= studyDuration)
    throw std::invalid_argument(
        "milestone cannot exceed studyDuration for fixed follow-up");

  // --- obtain criticalValues if missing ---
  std::vector<double> l(K, -6.0), zero(K, 0.0);
  std::vector<double> critValues = criticalValues;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (K > 1 && criticalValues.size() == K) {
      bool hasNaN = false;
      for (size_t i = 0; i < K - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[K-1])) haybittle = true;
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(K);
      for (size_t i = 0; i < K - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[K-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[K-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  }

  // compute cumulative alpha spent (and alpha1)
  ListCpp probs = exitprobcpp(critValues, l, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(K);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha :
    std::round(cumAlphaSpent.back() * 1e6) / 1e6;

  // futility bounds handling
  std::vector<double> futBounds = futilityBounds;
  if (K > 1) {
    if (missingFutilityBounds && bsf == "none") {
      futBounds = std::vector<double>(K, -6.0);
      futBounds[K-1] = critValues[K-1];
    } else if (!missingFutilityBounds && futBounds.size() == K-1) {
      futBounds.push_back(critValues[K-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }

  // studyDuration1
  double studyDuration1 = studyDuration;
  if (!fixedFollowup || std::isnan(studyDuration))
    studyDuration1 = accrualDuration + followupTime;

  // Prepare containers for stagewise results
  std::vector<double> time(K), info(K);
  std::vector<double> nsubjects(K), nevents(K), ndropouts(K), nmilestone(K);

  // compute information using twin group trick
  std::vector<double> accrualInt2 = accrualIntensity;
  for (double& v : accrualInt2) v *= 2.0;

  DataFrameCpp km_end = kmstat1cpp(
    studyDuration1, milestone, 1,
    accrualTime, accrualInt2,
    piecewiseSurvivalTime, stratumFraction,
    lambdax, lambdax, gammax, gammax,
    accrualDuration, followupTime, fixedFollowup);

  const auto& s  = km_end.get<double>("surv1");
  const auto& vs = km_end.get<double>("vsurv1");
  double surv = 0.0, vsurv = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    surv  += stratumFraction[h] * s[h];
    vsurv += stratumFraction[h] * stratumFraction[h] * vs[h];
  }
  double maxInformation = 1.0 / vsurv;
  std::vector<double> theta(K, surv - survH0);

  // set final-stage quantities
  info[K - 1] = maxInformation;
  time[K - 1] = studyDuration1;

  nsubjects[K - 1]  = extract_km(km_end, "subjects") * 0.5; // half because twin-group
  nevents[K - 1]    = extract_km(km_end, "nevents1");
  ndropouts[K - 1]  = extract_km(km_end, "ndropouts1");
  nmilestone[K - 1] = extract_km(km_end, "nmilestone1");

  // Compute interim analysis times and stagewise counts
  for (size_t i = 0; i < K - 1; ++i) {
    double information1 = maxInformation * infoRates[i];
    info[i] = information1;

    // function to match predicted information at candidate time
    auto g = [&](double t)->double {
      DataFrameCpp km1 = kmstat1cpp(
        t, milestone, 1,
        accrualTime, accrualInt2,
        piecewiseSurvivalTime, stratumFraction,
        lambdax, lambdax, gammax, gammax,
        accrualDuration, followupTime, fixedFollowup);

      const auto& vs = km1.get<double>("vsurv1");
      double vsurv = 0.0;
      for (size_t h = 0; h < nstrata; ++h) {
        vsurv += stratumFraction[h] * stratumFraction[h] * vs[h];
      }
      return 1.0 / vsurv - information1;
    };

    // bracket for brent: between slightly after milestone and studyDuration1
    time[i] = brent(g, milestone + 0.001, studyDuration1, 1e-6);

    // get counts at time[i]
    DataFrameCpp km_i = kmstat1cpp(
      time[i], milestone, 1,
      accrualTime, accrualInt2,
      piecewiseSurvivalTime, stratumFraction,
      lambdax, lambdax, gammax, gammax,
      accrualDuration, followupTime, fixedFollowup);

    nsubjects[i]  = 0.5 * extract_km(km_i, "subjects"); // half because twin-group
    nevents[i]    = extract_km(km_i, "nevents1");
    ndropouts[i]  = extract_km(km_i, "ndropouts1");
    nmilestone[i] = extract_km(km_i, "nmilestone1");
  }

  // compute exit probabilities
  ListCpp exit_probs;
  if (!missingFutilityBounds || bsf == "none" || K == 1) {
    exit_probs = exitprobcpp(critValues, futBounds, theta, info);
  } else {
    std::vector<double> w(K, 1.0);
    auto gp = getPower(alpha1, kMax, critValues, theta, info, bsf,
                       parameterBetaSpending, spendTime, futStopping, w);
    futBounds = gp.get<std::vector<double>>("futilityBounds");
    exit_probs = gp.get_list("probs");
  }

  // efficacy/futility p-values (implied)
  std::vector<double> efficacyP(K), futilityP(K);
  for (size_t i = 0; i < K; ++i) {
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
    futilityP[i] = 1.0 - boost_pnorm(futBounds[i]);
  }

  auto pu = exit_probs.get<std::vector<double>>("exitProbUpper");
  auto pl = exit_probs.get<std::vector<double>>("exitProbLower");
  std::vector<double> ptotal(K);
  for (size_t i = 0; i < K; ++i) ptotal[i] = pu[i] + pl[i];

  double expectedNumberOfEvents = 0.0;
  double expectedNumberOfDropouts = 0.0;
  double expectedNumberOfSubjects = 0.0;
  double expectedNumberOfMilestone = 0.0;
  double expectedStudyDuration = 0.0;
  double expectedInformation = 0.0;
  for (size_t i = 0; i < K; ++i) {
    expectedNumberOfEvents += ptotal[i] * nevents[i];
    expectedNumberOfDropouts += ptotal[i] * ndropouts[i];
    expectedNumberOfSubjects += ptotal[i] * nsubjects[i];
    expectedNumberOfMilestone += ptotal[i] * nmilestone[i];
    expectedStudyDuration += ptotal[i] * time[i];
    expectedInformation += ptotal[i] * info[i];
  }

  std::vector<double> cpu(K), cpl(K);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());
  double overallReject = cpu[K-1];

  // implied survival boundaries
  std::vector<double> survu(K), survl(K);
  for (size_t i = 0; i < K; ++i) {
    survu[i] = survH0 + critValues[i] / std::sqrt(info[i]);
    survl[i] = survH0 + futBounds[i] / std::sqrt(info[i]);
    if (critValues[i] == 6.0) { survu[i] = NaN; effStopping[i] = 0; }
    if (futBounds[i] == -6.0) { survl[i] = NaN; futStopping[i] = 0; }
  }

  // Build byStageResults DataFrameCpp
  // overallResults
  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha1, "alpha");
  overallResults.push_back(nevents.back(), "numberOfEvents");
  overallResults.push_back(ndropouts.back(), "numberOfDropouts");
  overallResults.push_back(nsubjects.back(), "numberOfSubjects");
  overallResults.push_back(nmilestone.back(), "numberOfMilestone");
  overallResults.push_back(time.back(), "studyDuration");
  overallResults.push_back(maxInformation, "information");
  overallResults.push_back(expectedNumberOfEvents, "expectedNumberOfEvents");
  overallResults.push_back(expectedNumberOfDropouts, "expectedNumberOfDropouts");
  overallResults.push_back(expectedNumberOfSubjects, "expectedNumberOfSubjects");
  overallResults.push_back(expectedNumberOfMilestone, "expectedNumberOfMilestone");
  overallResults.push_back(expectedStudyDuration, "expectedStudyDuration");
  overallResults.push_back(expectedInformation, "expectedInformation");
  overallResults.push_back(accrualDuration, "accrualDuration");
  overallResults.push_back(followupTime, "followupTime");
  overallResults.push_back(fixedFollowup, "fixedFollowup");
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(milestone, "milestone");
  overallResults.push_back(survH0, "survH0");
  overallResults.push_back(surv, "surv");

  DataFrameCpp byStageResults;
  byStageResults.push_back(std::move(infoRates), "informationRates");
  byStageResults.push_back(std::move(critValues), "efficacyBounds");
  byStageResults.push_back(std::move(futBounds), "futilityBounds");
  byStageResults.push_back(std::move(pu), "rejectPerStage");
  byStageResults.push_back(std::move(pl), "futilityPerStage");
  byStageResults.push_back(std::move(cpu), "cumulativeRejection");
  byStageResults.push_back(std::move(cpl), "cumulativeFutility");
  byStageResults.push_back(std::move(cumAlphaSpent), "cumulativeAlphaSpent");
  byStageResults.push_back(std::move(nevents), "numberOfEvents");
  byStageResults.push_back(std::move(ndropouts), "numberOfDropouts");
  byStageResults.push_back(std::move(nsubjects), "numberOfSubjects");
  byStageResults.push_back(std::move(nmilestone), "numberOfMilestone");
  byStageResults.push_back(std::move(time), "analysisTime");
  byStageResults.push_back(std::move(survu), "efficacySurv");
  byStageResults.push_back(std::move(survl), "futilitySurv");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(futilityP), "futilityP");
  byStageResults.push_back(std::move(info), "information");
  byStageResults.push_back(std::move(effStopping), "efficacyStopping");
  byStageResults.push_back(std::move(futStopping), "futilityStopping");

  // settings list
  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
  settings.push_back(typeBetaSpending, "typeBetaSpending");
  settings.push_back(parameterBetaSpending, "parameterBetaSpending");
  settings.push_back(accrualTime, "accrualTime");
  settings.push_back(accrualIntensity, "accrualIntensity");
  settings.push_back(piecewiseSurvivalTime, "piecewiseSurvivalTime");
  settings.push_back(stratumFraction, "stratumFraction");
  settings.push_back(lambda, "lambda");
  settings.push_back(gamma, "gamma");
  settings.push_back(spendingTime, "spendingTime");

  // assemble result
  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(settings), "settings");

  return result;
}


//' @title Power for One-Sample Milestone Survival Probability
//' @description Estimates the power, stopping probabilities, and expected
//' sample size in a one-group survival design.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @param typeBetaSpending The type of beta spending. One of the following:
//'   "sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock
//'   type spending function, "sfKD" for Kim & DeMets spending function,
//'   "sfHSD" for Hwang, Shi & DeCani spending function, and "none" for no
//'   early futility stopping. Defaults to "none".
//' @inheritParams param_parameterBetaSpending
//' @param milestone The milestone time at which to calculate the survival
//'   probability.
//' @param survH0 The milestone survival probability under the null
//'   hypothesis.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param lambda A vector of hazard rates for the event in each analysis
//'  time interval by stratum under the alternative hypothesis.
//' @param gamma The hazard rate for exponential dropout or a vector of
//'   hazard rates for piecewise exponential dropout. Defaults to 0 for
//'   no dropout.
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//'
//' @return An S3 class \code{kmpower1s} object with 3 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numbeOfSubjects}: The total number of subjects.
//'
//'     - \code{numberOfMilestone}: The total number of subjects reaching
//'       milestone.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedNumberOfMilestone}: The expected number of subjects
//'       reaching milestone.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up duration.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{milestone}: The milestone time to calculate the survival
//'       probability.
//'
//'     - \code{survH0}: The milestone survival probability under the null
//'       hypothesis.
//'
//'     - \code{surv}: The milestone survival probability under the
//'       alternative hypothesis.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale.
//'
//'     - \code{futilityBounds}: The futility boundaries on the Z-scale.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{futilityPerStage}: The probability for futility stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeFutility}: The cumulative probability for futility
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha spent.
//'
//'     - \code{numberOfEvents}: The number of events.
//'
//'     - \code{numberOfDropouts}: The number of dropouts.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//'     - \code{numberOfMilestone}: The number of subjects reaching
//'       milestone.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacySurv}: The efficacy boundaries on the milestone
//'       survival probability scale.
//'
//'     - \code{futilitySurv}: The futility boundaries on the milestone
//'       survival probability scale.
//'
//'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'
//'     - \code{futilityP}: The futility boundaries on the p-value scale.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{efficacyStopping}: Whether to allow efficacy stopping.
//'
//'     - \code{futilityStopping}: Whether to allow futility stopping.
//'
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{typeBetaSpending},
//'   \code{parameterBetaSpending}, \code{accrualTime},
//'   \code{accuralIntensity}, \code{piecewiseSurvivalTime},
//'   \code{stratumFraction}, \code{lambda}, \code{gamma},
//'   and \code{spendingTime}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{kmstat}}
//'
//' @examples
//'
//' kmpower1s(kMax = 2, informationRates = c(0.8, 1),
//'           alpha = 0.025, typeAlphaSpending = "sfOF",
//'           milestone = 18, survH0 = 0.30,
//'           accrualTime = seq(0, 8),
//'           accrualIntensity = 26/9*seq(1, 9),
//'           piecewiseSurvivalTime = c(0, 6),
//'           stratumFraction = c(0.2, 0.8),
//'           lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'           gamma = -log(1-0.05)/12, accrualDuration = 22,
//'           followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List kmpower1s(
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& futilityBounds = NA_REAL,
    const std::string& typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const double milestone = NA_REAL,
    const double survH0 = NA_REAL,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& lambda = NA_REAL,
    const Rcpp::NumericVector& gamma = 0,
    const double accrualDuration = NA_REAL,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = false,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const double studyDuration = NA_REAL) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam = Rcpp::as<std::vector<double>>(lambda);
  auto gam = Rcpp::as<std::vector<double>>(gamma);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  ListCpp out = kmpower1scpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, survH0,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    lam, gam,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, studyDuration);

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "kmpower1s";
  return result;
}


ListCpp kmsamplesize1scpp(
    const double beta,
    const int kMax,
    const std::vector<double>& informationRates,
    const std::vector<unsigned char>& efficacyStopping,
    const std::vector<unsigned char>& futilityStopping,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const std::vector<double>& futilityBounds,
    const std::string& typeBetaSpending,
    const double parameterBetaSpending,
    const std::vector<double>& userBetaSpending,
    const double milestone,
    const double survH0,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda,
    const std::vector<double>& gamma,
    double accrualDuration,
    double followupTime,
    const bool fixedFollowup,
    const std::vector<double>& spendingTime,
    const bool rounding) {

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1))
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  if (beta < 0.0001 || (!std::isnan(alpha) && beta >= 1.0 - alpha))
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");
  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");
  size_t K = static_cast<size_t>(kMax);

  std::vector<double> infoRates(K);
  if (none_na(informationRates)) {
    if (informationRates.size() != K)
      throw std::invalid_argument("Invalid length for informationRates");
    if (informationRates[0] <= 0.0)
      throw std::invalid_argument("informationRates must be positive");
    if (any_nonincreasing(informationRates))
      throw std::invalid_argument("informationRates must be increasing");
    if (informationRates[K-1] != 1.0)
      throw std::invalid_argument("informationRates must end with 1");
    infoRates = informationRates; // copy
  } else {
    for (size_t i = 0; i < K; ++i)
      infoRates[i] = static_cast<double>(i+1) / static_cast<double>(K);
  }

  std::vector<unsigned char> effStopping;
  if (none_na(efficacyStopping)) {
    if (efficacyStopping.size() != K)
      throw std::invalid_argument("Invalid length for efficacyStopping");
    if (efficacyStopping[K-1] != 1)
      throw std::invalid_argument("efficacyStopping must end with 1");
    effStopping = efficacyStopping; // copy
  } else {
    effStopping.assign(K, 1);
  }

  std::vector<unsigned char> futStopping;
  if (none_na(futilityStopping)) {
    if (futilityStopping.size() != K)
      throw std::invalid_argument("Invalid length for futilityStopping");
    if (futilityStopping[K-1] != 1)
      throw std::invalid_argument("futilityStopping must end with 1");
    futStopping = futilityStopping; // copy
  } else {
    futStopping.assign(K, 1);
  }

  bool missingCriticalValues = !none_na(criticalValues);
  bool missingFutilityBounds = !none_na(futilityBounds);

  if (!missingCriticalValues && criticalValues.size() != K) {
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
    if (userAlphaSpending.size() != K)
      throw std::invalid_argument("Invalid length of userAlphaSpending");
    if (userAlphaSpending[0] < 0.0)
      throw std::invalid_argument("userAlphaSpending must be nonnegative");
    if (any_nonincreasing(userAlphaSpending))
      throw std::invalid_argument("userAlphaSpending must be nondecreasing");
    if (userAlphaSpending[K-1] != alpha)
      throw std::invalid_argument("userAlphaSpending must end with specified alpha");
  }

  if (!missingFutilityBounds) {
    if (!(futilityBounds.size() == K - 1 || futilityBounds.size() == K)) {
      throw std::invalid_argument("Invalid length for futilityBounds");
    }
  }
  if (!missingCriticalValues && !missingFutilityBounds) {
    for (size_t i = 0; i < K - 1; ++i) {
      if (futilityBounds[i] > criticalValues[i]) {
        throw std::invalid_argument("futilityBounds must lie below criticalValues");
      }
    }
    if (futilityBounds.size() == K &&
        futilityBounds[K-1] != criticalValues[K-1]) {
      throw std::invalid_argument(
          "futilityBounds must meet criticalValues at the final look");
    }
  }

  std::string bsf = typeBetaSpending;
  for (char &c : bsf) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }

  if (missingFutilityBounds && !(bsf == "sfof" || bsf == "sfp" ||
      bsf == "sfkd" || bsf == "sfhsd" || bsf=="user" || bsf == "none")) {
    throw std::invalid_argument("Invalid value for typeBetaSpending");
  }
  if ((bsf == "sfkd" || bsf == "sfhsd") && std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument("parameterBetaSpending must be positive for sfKD");
  }

  if (missingFutilityBounds && bsf=="user") {
    if (!none_na(userBetaSpending))
      throw std::invalid_argument("userBetaSpending must be specified");
    if (userBetaSpending.size() != K)
      throw std::invalid_argument("Invalid length of userBetaSpending");
    if (userBetaSpending[0] < 0.0)
      throw std::invalid_argument("userBetaSpending must be nonnegative");
    if (any_nonincreasing(userBetaSpending))
      throw std::invalid_argument("userBetaSpending must be nondecreasing");
    if (userBetaSpending[K-1] != beta)
      throw std::invalid_argument("userBetaSpending must end with specified beta");
  }

  if (std::isnan(milestone)) throw std::invalid_argument("milestone must be provided");
  if (milestone <= 0.0) throw std::invalid_argument("milestone must be positive");
  if (survH0 <= 0.0 || survH0 >= 1.0)
    throw std::invalid_argument("survDiffH0 must lie between 0 and 1");

  if (accrualTime[0] != 0.0)
    throw std::invalid_argument("accrualTime must start with 0");
  if (any_nonincreasing(accrualTime))
    throw std::invalid_argument("accrualTime should be increasing");

  if (!none_na(accrualIntensity))
    throw std::invalid_argument("accrualIntensity must be provided");
  if (accrualIntensity.size() != accrualTime.size())
    throw std::invalid_argument("Invalid length for accrualIntensity");
  for (double v : accrualIntensity) {
    if (v < 0.0) throw std::invalid_argument("accrualIntensity must be non-negative");
  }

  if (piecewiseSurvivalTime[0] != 0.0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(piecewiseSurvivalTime))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");

  for (double v : stratumFraction) {
    if (v <= 0.0) throw std::invalid_argument("stratumFraction must be positive");
  }
  double sumf = std::accumulate(stratumFraction.begin(), stratumFraction.end(), 0.0);
  if (std::fabs(sumf - 1.0) > 1e-12)
    throw std::invalid_argument("stratumFraction must sum to 1");

  if (none_na(lambda) == false) throw std::invalid_argument("lambda must be provided");
  for (double v : lambda) {
    if (v < 0.0) throw std::invalid_argument("lambda must be non-negative");
  }
  for (double v : gamma) {
    if (v < 0.0) throw std::invalid_argument("gamma must be non-negative");
  }

  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambdax = expand_stratified(lambda, nstrata, nintv, "lambda");
  auto gammax = expand_stratified(gamma, nstrata, nintv, "gamma");

  if (!std::isnan(accrualDuration) && accrualDuration <= 0.0)
    throw std::invalid_argument("accrualDuration must be positive");
  if (!std::isnan(followupTime)) {
    if (fixedFollowup && followupTime <= 0.0)
      throw std::invalid_argument("followupTime must be positive for fixed follow-up");
    if (!fixedFollowup && followupTime < 0.0)
      throw std::invalid_argument(
          "followupTime must be non-negative for variable follow-up");
  }
  if (fixedFollowup && std::isnan(followupTime))
    throw std::invalid_argument("followupTime must be provided for fixed follow-up");

  if (!std::isnan(accrualDuration) && !std::isnan(followupTime) &&
      (milestone >= accrualDuration + followupTime))
    throw std::invalid_argument(
        "milestone must be less than accrualDuration + followupTime");
  if (fixedFollowup && milestone > followupTime)
    throw std::invalid_argument(
        "milestone cannot exceed followupTime for fixed follow-up");

  // spendingTime default to informationRates
  std::vector<double> spendTime;
  if (none_na(spendingTime)) {
    if (spendingTime.size() != K)
      throw std::invalid_argument("Invalid length for spendingTime");
    if (spendingTime[0] <= 0.0)
      throw std::invalid_argument("spendingTime must be positive");
    if (any_nonincreasing(spendingTime))
      throw std::invalid_argument("spendingTime must be increasing");
    if (spendingTime[K-1] != 1.0)
      throw std::invalid_argument("spendingTime must end with 1");
    spendTime = spendingTime; // copy
  } else {
    spendTime = infoRates;
  }


  // --- obtain criticalValues if missing ---
  std::vector<double> l(K, -6.0), zero(K, 0.0);
  std::vector<double> critValues = criticalValues;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (K > 1 && criticalValues.size() == K) {
      bool hasNaN = false;
      for (size_t i = 0; i < K - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[K-1])) haybittle = true;
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(K);
      for (size_t i = 0; i < K - 1; ++i) {
        u[i] = criticalValues[i];
        if (!effStopping[i]) u[i] = 6.0;
      }

      auto f = [&](double aval)->double {
        u[K-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[K-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  }

  // compute cumulative alpha spent (and alpha1)
  ListCpp probs = exitprobcpp(critValues, l, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(K);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());
  double alpha1 = missingCriticalValues ? alpha :
    std::round(cumAlphaSpent.back() * 1e6) / 1e6;


  // futility bounds handling
  std::vector<double> futBounds = futilityBounds;
  if (K > 1) {
    if (missingFutilityBounds && bsf == "none") {
      futBounds = std::vector<double>(K, -6.0);
      futBounds[K-1] = critValues[K-1];
    } else if (!missingFutilityBounds && futBounds.size() == K-1) {
      futBounds.push_back(critValues[K-1]);
    }
  } else {
    if (missingFutilityBounds) {
      futBounds = critValues;
    }
  }


  // Which design parameter is unknown?
  enum Unknown { ACC_DUR, FUP_TIME, ACC_INT };
  Unknown unknown;
  bool missAccrual = std::isnan(accrualDuration);
  bool missFollow = std::isnan(followupTime);
  if (missAccrual && !missFollow) unknown = ACC_DUR;
  else if (!missAccrual && missFollow) unknown = FUP_TIME;
  else if (!missAccrual && !missFollow) unknown = ACC_INT;
  else throw std::invalid_argument(
      "accrualDuration and followupTime cannot be both missing");

  // local accrualIntensity copy (we might scale it)
  std::vector<double> accrualInt = accrualIntensity;

  // milestone survival under H1
  std::vector<double> zerogam(nintv, 0.0);
  double surv = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    double p = patrisk1(milestone, piecewiseSurvivalTime, lambdax[h], zerogam);
    surv += stratumFraction[h] * p;
  }
  double theta1 = surv - survH0;
  std::vector<double> theta(K, theta1);

  // get design under H1 (delegates much work)
  ListCpp design = getDesigncpp(
    beta, NaN, theta1, kMax, infoRates,
    effStopping, futStopping,
    critValues, alpha1, asf,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, bsf, parameterBetaSpending,
    userBetaSpending, spendTime, 1.0);

  auto byStageResults = design.get<DataFrameCpp>("byStageResults");
  futBounds = byStageResults.get<double>("futilityBounds");

  auto overallResults = design.get<DataFrameCpp>("overallResults");
  double maxInformation = overallResults.get<double>("information")[0];

  // Helper: compute information under H1 given accrualDuration (accrDur),
  // followupTime (fu), and accrualIntensity (accrInt).
  auto info_minus_target_H1 = [&](double t, double accrDur, double fu,
                                  const std::vector<double>& accrInt) {
    // compute information using twin group trick
    std::vector<double> accrInt2 = accrInt;
    for (double& v : accrInt2) v *= 2.0;

    DataFrameCpp km1 = kmstat1cpp(
      t, milestone, 1,
      accrualTime, accrInt2,
      piecewiseSurvivalTime, stratumFraction,
      lambdax, lambdax, gammax, gammax,
      accrDur, fu, fixedFollowup);

    const auto& vs = km1.get<double>("vsurv1");
    double vsurv = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      vsurv += stratumFraction[h] * stratumFraction[h] * vs[h];
    }
    return 1.0 / vsurv - maxInformation;
  };

  // when info at zero follow-up (end of enrollment) already exceeds target,
  // set followupTime = 0 and find minimal accrualDuration achieving target
  bool curtailed = false;

  // when info at maximum follow-up (milestone) is still below target,
  // increase accrualDuration to achieve target
  bool expanded = false;

  // Root-finding function to solve for the unknown design parameter to achieve target
  std::function<double(double)> f_root;

  if (unknown == ACC_DUR) {
    f_root = [&](double x)->double {
      return info_minus_target_H1(x + followupTime, x, followupTime, accrualInt);
    };
  } else if (unknown == FUP_TIME) {
    if (!fixedFollowup && accrualDuration > milestone + 0.001 &&
        info_minus_target_H1(accrualDuration, accrualDuration, 0.0, accrualInt) > 0) {
      thread_utils::push_thread_warning(
        "Information at zero follow-up (end of enrollment) already exceeds target. "
        "Setting followupTime = 0 and finding minimal accrualDuration.");

      f_root = [&](double x)->double {
        return info_minus_target_H1(x, x, 0.0, accrualInt);
      };

      curtailed = true;
    } else if (info_minus_target_H1(accrualDuration + milestone, accrualDuration,
                                    milestone, accrualInt) < 0) {
      thread_utils::push_thread_warning(
        "The required information cannot be attained by increasing "
        "followupTime alone. accrualDuration is also increased to "
        "attain the required information.");

      f_root = [&](double x)->double {
        return info_minus_target_H1(x + milestone, x, milestone, accrualInt);
      };

      expanded = true;
    } else {
      f_root = [&](double x)->double {
        return info_minus_target_H1(accrualDuration + x, accrualDuration, x,
                                    accrualInt);
      };
    }
  } else {
    f_root = [&](double m)->double {
      std::vector<double> scaled = accrualInt;
      for (double &v : scaled) v *= m;
      return info_minus_target_H1(accrualDuration + followupTime,
                                  accrualDuration, followupTime, scaled);
    };
  }

  double lower, upper;
  if (unknown == ACC_DUR) {
    lower = std::max(milestone - followupTime, 0.0) + 0.001;
    upper = 2.0 * lower;
  } else if (unknown == FUP_TIME && curtailed) {
    lower = milestone + 0.001;
    upper = accrualDuration;
  } else if (unknown == FUP_TIME && expanded) {
    lower = accrualDuration;
    upper = 2.0 * lower;
  } else if (unknown == FUP_TIME && !curtailed && !expanded) {
    lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
    upper = milestone;
  } else { // unknown == ACC_INT
    lower = 0.001;
    upper = 120;
  }

  // expand upper if needed to ensure root is bracketed
  double fl_val = f_root(lower), fu_val = f_root(upper);
  if (unknown == ACC_DUR || (unknown == FUP_TIME && expanded) || unknown == ACC_INT) {
    int expand_iter = 0;
    while (fl_val * fu_val > 0.0 && expand_iter < 60) {
      lower = upper;
      fl_val = fu_val;
      upper *= 2.0;
      fu_val = f_root(upper);
      ++expand_iter;
    }
  }
  if (fl_val * fu_val > 0.0) throw std::runtime_error(
      "Unable to bracket root; check interval or inputs");

  // solve for root and apply solution for the unknown design parameter
  auto f_for_brent1 = [&](double x)->double {
    if (x == lower) return fl_val;
    if (x == upper) return fu_val;
    return f_root(x);
  };
  double solution = brent(f_for_brent1, lower, upper, 1e-6);

  if (unknown == ACC_DUR) {
    accrualDuration = solution;
  } else if (unknown == FUP_TIME && curtailed) {
    followupTime = 0.0;
    accrualDuration = solution;
  } else if (unknown == FUP_TIME && expanded) {
    followupTime = milestone;
    accrualDuration = solution;
  } else if (unknown == FUP_TIME && !curtailed && !expanded) {
    followupTime = solution;
  } else { // scaled multiplier for accrualIntensity
    for (double &v : accrualInt) v *= solution;
  }

  double studyDuration = accrualDuration + followupTime;

  // --- rounding to integer N (adjust intensity or accrualDuration) ---
  if (rounding) {
    double n0 = accrual1(studyDuration, accrualTime, accrualInt, accrualDuration);
    double n = std::ceil(n0 - 1.0e-12);

    if (n - n0 > 1e-6) {
      if (unknown == ACC_INT) { // scale intensity to hit integer n
        double mult = n / n0;
        for (double& v : accrualInt) v *= mult;
      } else { // adjust accrualDuration to hit integer n
        accrualDuration = getAccrualDurationFromN1(n, accrualTime, accrualInt);
      }

      if (!fixedFollowup) {
        // variable follow-up: adjust follow-up time to match maxInformation
        auto h = [&](double x)->double {
          return info_minus_target_H1(accrualDuration + x, accrualDuration, x,
                                      accrualInt);
        };
        double lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // fixed follow-up: adjust studyDuration by extending post-accrual time
        auto h = [&](double x)->double {
          return info_minus_target_H1(accrualDuration + x, accrualDuration,
                                      followupTime, accrualInt);
        };
        double lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
        double extra = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + extra;
      }
    }
  }

  // --- Results under H1 ---
  ListCpp resultsH1 = kmpower1scpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha1, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, survH0,
    accrualTime, accrualInt,
    piecewiseSurvivalTime, stratumFraction,
    lambda, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  // refresh maxInformation from H1 power results    {
  DataFrameCpp overallResultsH1 = resultsH1.get<DataFrameCpp>("overallResults");
  maxInformation = overallResultsH1.get<double>("information")[0];

  // --- Construct lambda under H0 (scale lambda to match survH0) ---
  // Solve for multiplier a s.t. weighted_surv(a*lambda) - survH0 == 0
  auto f_surv = [&](double aval)->double {
    double surv = 0.0;
    for (std::size_t h = 0; h < nstrata; ++h) {
      auto lamH0 = lambdax[h]; // build per-stratum scaled lambda
      for (double &v : lamH0) v *= aval;
      double p = patrisk1(milestone, piecewiseSurvivalTime, lamH0, zerogam);
      surv += stratumFraction[h] * p;
    }
    return surv - survH0;
  };

  double lo = 0.001, hi = 2.0;
  double flo = f_surv(lo);
  double fhi = f_surv(hi);
  int it = 0;
  while (flo * fhi > 0.0 && it < 60) {
    lo = hi;
    flo = fhi;
    hi *= 2.0;
    fhi = f_surv(hi);
  }
  if (flo * fhi > 0.0) throw std::runtime_error(
      "Unable to bracket root for lambda under H0; check inputs");

  auto f_for_brent2 = [&](double x)->double {
    if (x == lo) return flo;
    if (x == hi) return fhi;
    return f_surv(x);
  };

  double multH0 = brent(f_for_brent2, lo, hi, 1.0e-6);

  std::vector<double> lambdaH0 = lambda;
  for (double &v : lambdaH0) v *= multH0;
  auto lambdaH0x = expand_stratified(lambdaH0, nstrata, nintv, "lambdaH0");

  // Helper: compute information under H0 given accrualDuration (accrDur),
  // followupTime (fu), and accrualIntensity (accrInt).
  auto info_minus_target_H0 = [&](double t, double accrDur, double fu,
                                  const std::vector<double>& accrInt) {
    // compute information using twin group trick
    std::vector<double> accrInt2 = accrInt;
    for (double& v : accrInt2) v *= 2.0;

    DataFrameCpp km1 = kmstat1cpp(
      t, milestone, 1,
      accrualTime, accrInt2,
      piecewiseSurvivalTime, stratumFraction,
      lambdaH0x, lambdaH0x, gammax, gammax,
      accrDur, fu, fixedFollowup);

    const auto& vs = km1.get<double>("vsurv1");
    double vsurv = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      vsurv += stratumFraction[h] * stratumFraction[h] * vs[h];
    }
    return 1.0 / vsurv - maxInformation;
  };

  if (!fixedFollowup) {
    auto h_follow = [&](double x)->double {
      return info_minus_target_H0(accrualDuration + x, accrualDuration, x, accrualInt);
    };

    if (accrualDuration > milestone + 0.001 && h_follow(0.0) > 0.0) {
      thread_utils::push_thread_warning(
        "Information at zero follow-up (end of enrollment) already exceeds target. "
        "Setting followupTime = 0 and finding minimal accrualDuration.");

      auto h_accr = [&](double x)->double {
        return info_minus_target_H0(x, x, 0.0, accrualInt);
      };

      double lo = milestone + 0.001;
      double hi = accrualDuration;
      accrualDuration = brent(h_accr, lo, hi, 1.0e-6);
      followupTime = 0.0;
      studyDuration = accrualDuration;
    } else if (h_follow(followupTime) < 0.0) {
      thread_utils::push_thread_warning(
        "The required information cannot be attained by increasing "
        "followupTime alone. accrualDuration is also increased to "
        "attain the required information.");

      auto h_accr = [&](double x)->double {
        return info_minus_target_H0(x + milestone, x, milestone, accrualInt);
      };

      // adjust accrualDuration upward to match target
      double lo = accrualDuration;
      double hi = 2.0 * accrualDuration;
      double flo = h_accr(lo);
      double fhi = h_accr(hi);
      int it = 0;
      while (flo * fhi > 0.0 && it < 60) {
        lo = hi;
        flo = fhi;
        hi *= 2.0;
        fhi = h_accr(hi);
      }
      if (flo * fhi > 0.0) throw std::runtime_error(
          "Unable to bracket accrualDuration under H0; check inputs");

      auto h_for_brent = [&](double x)->double {
        if (x == lo) return flo;
        if (x == hi) return fhi;
        return h_accr(x);
      };

      accrualDuration = brent(h_for_brent, lo, hi, 1.0e-6);
      followupTime = milestone;
      studyDuration = accrualDuration + followupTime;
    } else {
      // adjust follow-up time to match target
      double lo = std::max(milestone - accrualDuration, 0.0) + 0.001;
      double hi = milestone;

      followupTime = brent(h_follow, lo, hi, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    }
  } else {
    auto h_study = [&](double x)->double {
      return info_minus_target_H0(accrualDuration + x, accrualDuration,
                                  followupTime, accrualInt);
    };

    auto h_accr = [&](double x)->double {
      return info_minus_target_H0(x + followupTime, x, followupTime, accrualInt);
    };

    if (accrualDuration > milestone + 0.001 && h_study(0.0) > 0.0) {
      thread_utils::push_thread_warning(
        "Information at zero follow-up (end of enrollment) already exceeds target. "
        "Decrease accrual druation.");

      double lo = milestone + 0.001;
      double hi = accrualDuration;
      accrualDuration = brent(h_accr, lo, hi, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (h_study(followupTime) < 0.0) {
      thread_utils::push_thread_warning(
        "The required information cannot be attained by increasing "
        "followupTime alone. accrualDuration is also increased to "
        "attain the required information.");

      double lo = accrualDuration;
      double hi = 2.0 * accrualDuration;
      double flo = h_accr(lo);
      double fhi = h_accr(hi);
      int it = 0;
      while (flo * fhi > 0.0 && it < 60) {
        lo = hi;
        flo = fhi;
        hi *= 2.0;
        fhi = h_accr(hi);
      }
      if (flo * fhi > 0.0) throw std::runtime_error(
          "Unable to bracket accrualDuration under H0; check inputs");

      auto h_for_brent = [&](double x)->double {
        if (x == lo) return flo;
        if (x == hi) return fhi;
        return h_accr(x);
      };

      accrualDuration = brent(h_for_brent, lo, hi, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else {
      double lo = std::max(milestone - accrualDuration, 0.0) + 0.001;
      double hi = followupTime;
      double extra = brent(h_study, lo, hi, 1.0e-6);
      studyDuration = accrualDuration + extra;
    }
  }

  // --- Results under H0 with same boundaries as H1 ---
  ListCpp resultsH0 = kmpower1scpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha1, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, survH0,
    accrualTime, accrualInt,
    piecewiseSurvivalTime, stratumFraction,
    lambdaH0, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  // add userBetaSpending to settings
  if (bsf == "user") {
    ListCpp& settingsH1 = resultsH1.get_list("settings");
    ListCpp& settingsH0 = resultsH0.get_list("settings");
    settingsH1.push_back(userBetaSpending, "userBetaSpending");
    settingsH0.push_back(userBetaSpending, "userBetaSpending");
  }

  ListCpp result;
  result.push_back(std::move(resultsH1), "resultsUnderH1");
  result.push_back(std::move(resultsH0), "resultsUnderH0");

  return result;
}


//' @title Sample Size for One-Sample Milestone Survival Probability
//' @description Obtains the needed accrual duration given power and
//' follow-up time, the needed follow-up time given power and
//' accrual duration, or the needed absolute accrual rates given
//' power, accrual duration, follow-up duration, and relative accrual
//' rates in a one-group survival design.
//'
//' @param beta Type II error. Defaults to 0.2.
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_efficacyStopping
//' @inheritParams param_futilityStopping
//' @inheritParams param_criticalValues
//' @inheritParams param_alpha
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @inheritParams param_futilityBounds
//' @inheritParams param_typeBetaSpending
//' @inheritParams param_parameterBetaSpending
//' @inheritParams param_userBetaSpending
//' @param milestone The milestone time at which to calculate the survival
//'   probability.
//' @param survH0 The milestone survival probability under the null
//'   hypothesis.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param lambda A vector of hazard rates for the event in each analysis
//'  time interval by stratum under the alternative hypothesis.
//' @param gamma The hazard rate for exponential dropout or a vector of
//'   hazard rates for piecewise exponential dropout. Defaults to 0 for
//'   no dropout.
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupDuration, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return A list of two components:
//'
//' * \code{resultsUnderH1}: An S3 class \code{kmpower1s} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{kmpower1s} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{kmpower1s}}
//'
//' @examples
//' # Example 1: Obtains follow-up duration given power, accrual intensity,
//' # and accrual duration for variable follow-up
//'
//' kmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, survH0 = 0.30,
//'                accrualTime = seq(0, 8),
//'                accrualIntensity = 26/9*seq(1, 9),
//'                piecewiseSurvivalTime = c(0, 6),
//'                stratumFraction = c(0.2, 0.8),
//'                lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'                gamma = -log(1-0.05)/12, accrualDuration = 22,
//'                followupTime = NA, fixedFollowup = FALSE)
//'
//' # Example 2: Obtains accrual intensity given power, accrual duration, and
//' # follow-up duration for variable follow-up
//'
//' kmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, survH0 = 0.30,
//'                accrualTime = seq(0, 8),
//'                accrualIntensity = 26/9*seq(1, 9),
//'                piecewiseSurvivalTime = c(0, 6),
//'                stratumFraction = c(0.2, 0.8),
//'                lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'                gamma = -log(1-0.05)/12, accrualDuration = 22,
//'                followupTime = 18, fixedFollowup = FALSE)
//'
//'
//' # Example 3: Obtains accrual duration given power, accrual intensity, and
//' # follow-up duration for fixed follow-up
//'
//' kmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, survH0 = 0.30,
//'                accrualTime = seq(0, 8),
//'                accrualIntensity = 26/9*seq(1, 9),
//'                piecewiseSurvivalTime = c(0, 6),
//'                stratumFraction = c(0.2, 0.8),
//'                lambda = c(0.0533, 0.0309, 1.5*0.0533, 1.5*0.0309),
//'                gamma = -log(1-0.05)/12, accrualDuration = NA,
//'                followupTime = 18, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List kmsamplesize1s(
    const double beta = 0.2,
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::LogicalVector& efficacyStopping = NA_LOGICAL,
    const Rcpp::LogicalVector& futilityStopping = NA_LOGICAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.025,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& futilityBounds = NA_REAL,
    const std::string& typeBetaSpending = "none",
    const double parameterBetaSpending = NA_REAL,
    const Rcpp::NumericVector& userBetaSpending = NA_REAL,
    const double milestone = NA_REAL,
    const double survH0 = NA_REAL,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& lambda = NA_REAL,
    const Rcpp::NumericVector& gamma = 0,
    double accrualDuration = NA_REAL,
    double followupTime = NA_REAL,
    const bool fixedFollowup = false,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool rounding = true) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam = Rcpp::as<std::vector<double>>(lambda);
  auto gam = Rcpp::as<std::vector<double>>(gamma);
  auto userBeta = Rcpp::as<std::vector<double>>(userBetaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  ListCpp out = kmsamplesize1scpp(
    beta, kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    userBeta, milestone, survH0,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    lam, gam,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, rounding);

  thread_utils::drain_thread_warnings_to_R();

  ListCpp resultsUnderH1 = out.get_list("resultsUnderH1");
  ListCpp resultsUnderH0 = out.get_list("resultsUnderH0");

  Rcpp::List resultsH1 = Rcpp::wrap(resultsUnderH1);
  Rcpp::List resultsH0 = Rcpp::wrap(resultsUnderH0);

  resultsH1.attr("class") = "kmpower1s";
  resultsH0.attr("class") = "kmpower1s";

  return Rcpp::List::create(
    Rcpp::Named("resultsUnderH1") = resultsH1,
    Rcpp::Named("resultsUnderH0") = resultsH0
  );
}

