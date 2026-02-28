#include "enrollment_event.h"
#include "generic_design.h"
#include "utilities.h"
#include "dataframe_list.h"

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

std::vector<double> rm_make_breaks(
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


double rmstcpp(const double t1,
               const double t2,
               const std::vector<double>& piecewiseSurvivalTime,
               const std::vector<double>& lambda) {

  // alias for convenience
  const std::vector<double>& t = piecewiseSurvivalTime;

  // identify the time interval containing the specified analysis time
  // t[m1] <= t1 < t[m1+1], t[m2] <= t2 < t[m2+1]
  int m1 = findInterval1(t1, t) - 1;
  int m2 = findInterval1(t2, t) - 1;

  double ch = 0.0;
  for (int j = 0; j < m1; ++j) {
    ch += lambda[j] * (t[j+1] - t[j]);
  }

  double result;
  if (m1 == m2) { // t1 and t2 lie in the same interval
    double lama = lambda[m1];
    double ta = t[m1];
    if (lama == 0.0) {
      result = std::exp(-ch) * (t2 - t1);
    } else {
      double s1 = std::exp(-lama * (t1 - ta));
      double s2 = std::exp(-lama * (t2 - ta));
      result = std::exp(-ch) * (s1 - s2) / lama;
    }
  } else {
    double lama = lambda[m1], lamb = lambda[m2];
    double ta = t[m1], tb = t[m2];

    // Partial interval from t1 to end of its interval
    if (lama == 0.0) {
      result = std::exp(-ch) * (t[m1+1] - t1);
    } else {
      double s1 = std::exp(-lama * (t1 - ta));
      double s2 = std::exp(-lama * (t[m1+1] - ta));
      result = std::exp(-ch) * (s1 - s2) / lama;
    }

    // Add full intervals between m1+1 and m2-1 (if any).
    // Before entering interval j we add the previous interval's full hazard to ch:
    for (int j = m1 + 1; j < m2; ++j) {
      ch += lambda[j-1] * (t[j] - t[j-1]);
      double lam = lambda[j];
      if (lam == 0.0) {
        result += std::exp(-ch) * (t[j+1] - t[j]);
      } else {
        double s2 = std::exp(-lam * (t[j+1] - t[j]));
        result += std::exp(-ch) * (1.0 - s2) / lam;
      }
    }

    // Add last partial interval from piecewise[m2] to t2.
    // Before doing that we must add the full hazard for interval m2-1
    ch += lambda[m2 - 1] * (tb - t[m2 - 1]);
    if (lamb == 0.0) {
      result += std::exp(-ch) * (t2 - tb);
    } else {
      double s2 = std::exp(-lamb * (t2 - tb));
      result += std::exp(-ch) * (1.0 - s2) / lamb;
    }
  }

  return result;
}


//' @title Restricted Mean Survival Time
//' @description Obtains the restricted mean survival time over an interval.
//'
//' @param t1 Lower bound of the analysis time interval.
//' @param t2 Upper bound of the analysis time interval.
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda
//'
//' @return The integral of the survival function from \code{t1} to \code{t2}
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' rmst(t1 = 0, t2 = 7, piecewiseSurvivalTime = c(0, 6),
//'      lambda = c(0.0533, 0.0309))
//'
//' @export
// [[Rcpp::export]]
double rmst(const double t1 = 0,
            const double t2 = NA_REAL,
            const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
            const Rcpp::NumericVector& lambda = NA_REAL) {

  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lam = Rcpp::as<std::vector<double>>(lambda);

  if (std::isnan(t2)) throw std::invalid_argument("t2 must be provided");
  if (t1 < 0.0) throw std::invalid_argument("t1 must be non-negative");
  if (t2 < t1) throw std::invalid_argument("t1 must be less than or equal to t2");
  if (pwSurvT[0] != 0.0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(pwSurvT))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");
  if (!none_na(lam)) {
    throw std::invalid_argument("lambda must be provided");
  }
  for (double v : lam) {
    if (v < 0.0) throw std::invalid_argument("lambda must be non-negative");
  }

  return rmstcpp(t1, t2, pwSurvT, lam);
}



double rm_integrand(
  const double x,    // analysis time (integration variable)
  const double t2,   // calendar time for look 2
  const double tau1, // milestone for look 1
  const double tau2, // milestone for look 2
  const double phi,
  const std::vector<double>& accrualTime,
  const std::vector<double>& accrualIntensity,
  const std::vector<double>& piecewiseSurvivalTime,
  const std::vector<double>& lambda,
  const std::vector<double>& gamma,
  const double accrualDuration) {

  double a1 = rmstcpp(x, tau1, piecewiseSurvivalTime, lambda);
  double a2 = rmstcpp(x, tau2, piecewiseSurvivalTime, lambda);

  // interval index for x in piecewiseSurvivalTime
  int j = findInterval1(x, piecewiseSurvivalTime) - 1;

  // p(x) = P(at risk at time x since randomization) under event+dropout hazards
  double p = patrisk1(x, piecewiseSurvivalTime, lambda, gamma);

  // N(t2 - x): number enrolled by calendar time (t2 - x)
  double N = accrual1(t2 - x, accrualTime, accrualIntensity, accrualDuration);

  // a1 * a2 * lambda_j / (phi * N * p)
  // guard against division by 0
  double denom = phi * N * p;
  if (denom <= 0.0) return 0.0;

  return a1 * a2 * lambda[j] / denom;
}


std::pair<double, double> covrmstcpp(
    double t2,
    double tau1,
    double tau2,
    double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    double accrualDuration,
    double maxFollowupTime) {

  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);
  auto breaks = rm_make_breaks(piecewiseSurvivalTime, accrualTime,
                               accrualDuration, maxFollowupTime,
                               t2, tau1);
  double tol = 1e-6;

  auto f1 = [&](double x)->double {
    return rm_integrand(x, t2, tau1, tau2, phi, accrualTime, accrualIntensity,
                       piecewiseSurvivalTime, lambda1, gamma1, accrualDuration);
  };
  auto f2 = [&](double x)->double {
    return rm_integrand(x, t2, tau1, tau2, 1.0 - phi, accrualTime, accrualIntensity,
                       piecewiseSurvivalTime, lambda2, gamma2, accrualDuration);
  };

  double q1 = integrate3(f1, breaks, tol);
  double q2 = integrate3(f2, breaks, tol);

  return std::make_pair(q1, q2);
}



//' @title Covariance Between Restricted Mean Survival Times
//' @description Obtains the covariance between restricted mean survival
//' times at two different time points.
//'
//' @param t2 The calendar time for analysis 2.
//' @param tau1 The milestone time for analysis 1.
//' @param tau2 The milestone time for analysis 2.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//' @inheritParams param_accrualDuration
//' @inheritParams param_maxFollowupTime
//'
//' @return The covariance between the restricted mean survival times
//' for each treatment group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' covrmst(t2 = 25, tau1 = 16, tau2 = 18, allocationRatioPlanned = 1,
//'         accrualTime = c(0, 3), accrualIntensity = c(10, 20),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 12, maxFollowupTime = 30)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector covrmst(
    const double t2 = NA_REAL,
    const double tau1 = NA_REAL,
    const double tau2 = NA_REAL,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double accrualDuration = NA_REAL,
    const double maxFollowupTime = NA_REAL) {

  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  // input validation
  if (std::isnan(t2)) throw std::invalid_argument("t2 must be provided");
  if (std::isnan(tau1)) throw std::invalid_argument("tau1 must be provided");
  if (tau1 <= 0.0) throw std::invalid_argument("tau1 must be positive");
  if (std::isnan(tau2)) throw std::invalid_argument("tau2 must be provided");
  if (tau2 < tau1)
    throw std::invalid_argument("tau2 must be greater than or equal to tau1");
  if (t2 <= tau2) throw std::invalid_argument("t2 must be greater than tau2");
  if (allocationRatioPlanned <= 0.0)
    throw std::invalid_argument("allocationRatioPlanned must be positive");
  if (accrualT[0] != 0.0)
    throw std::invalid_argument("accrualTime must start with 0");
  if (any_nonincreasing(accrualT))
    throw std::invalid_argument("accrualTime should be increasing");
  if (!none_na(accrualInt))
    throw std::invalid_argument("accrualIntensity must be provided");
  if (accrualT.size() != accrualInt.size())
    throw std::invalid_argument("Invalid length for accrualIntensity");
  for (double v : accrualInt) {
    if (v < 0.0) throw std::invalid_argument("accrualIntensity must be non-negative");
  }
  if (pwSurvT[0] != 0.0)
    throw std::invalid_argument("piecewiseSurvivalTime must start with 0");
  if (any_nonincreasing(pwSurvT))
    throw std::invalid_argument("piecewiseSurvivalTime should be increasing");
  if (!none_na(lam1)) throw std::invalid_argument("lambda1 must be provided");
  if (!none_na(lam2)) throw std::invalid_argument("lambda2 must be provided");
  for (double v : lam1) {
    if (v < 0.0) throw std::invalid_argument("lambda1 must be non-negative");
  }
  for (double v : lam2) {
    if (v < 0.0) throw std::invalid_argument("lambda2 must be non-negative");
  }
  for (double v : gam1) {
    if (v < 0.0) throw std::invalid_argument("gamma1 must be non-negative");
  }
  for (double v : gam2) {
    if (v < 0.0) throw std::invalid_argument("gamma2 must be non-negative");
  }
  if (std::isnan(accrualDuration))
    throw std::invalid_argument("accrualDuration must be provided");
  if (accrualDuration <= 0.0)
    throw std::invalid_argument("accrualDuration must be positive");
  if (std::isnan(maxFollowupTime))
    throw std::invalid_argument("maxFollowupTime must be provided");
  if (tau2 > maxFollowupTime)
    throw std::invalid_argument("tau2 must be less than or equal to maxFollowupTime");

  size_t nintv = pwSurvT.size();
  auto lambda1x = expand1(lam1, nintv, "lambda1");
  auto lambda2x = expand1(lam2, nintv, "lambda2");
  auto gamma1x = expand1(gam1, nintv, "gamma1");
  auto gamma2x = expand1(gam2, nintv, "gamma2");

  auto out = covrmstcpp(
    t2, tau1, tau2, allocationRatioPlanned,
    accrualT, accrualInt, pwSurvT,
    lambda1x, lambda2x, gamma1x, gamma2x,
    accrualDuration, maxFollowupTime);

  return Rcpp::NumericVector::create(out.first, out.second);
}


DataFrameCpp rmstat1cpp(
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

  std::vector<double> rmst1(nstrata), rmst2(nstrata), rmstDiff(nstrata);
  std::vector<double> vrmst1(nstrata), vrmst2(nstrata), vrmstDiff(nstrata);

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

    // rmst (event only: gamma = 0)
    rmst1[h] = rmstcpp(0, milestone, piecewiseSurvivalTime, lam1);
    rmst2[h] = rmstcpp(0, milestone, piecewiseSurvivalTime, lam2);
    rmstDiff[h] = rmst1[h] - rmst2[h];

    auto v = covrmstcpp(
      time, milestone, milestone, allocationRatioPlanned,
      accrualTime, accrualIntensity_frac, piecewiseSurvivalTime,
      lam1, lam2, gam1, gam2, accrualDuration, maxFollowupTime);
    vrmst1[h] = v.first;
    vrmst2[h] = v.second;
    vrmstDiff[h] = vrmst1[h] + vrmst2[h];
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
  df.push_back(std::move(rmst1), "rmst1");
  df.push_back(std::move(rmst2), "rmst2");
  df.push_back(std::move(rmstDiff), "rmstDiff");
  df.push_back(std::move(vrmst1), "vrmst1");
  df.push_back(std::move(vrmst2), "vrmst2");
  df.push_back(std::move(vrmstDiff), "vrmstDiff");

  return df;
}


DataFrameCpp rmstatcpp(
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
  if (fixedFollowup && milestone > followupTime)
    throw std::invalid_argument(
        "milestone cannot exceed followupTime for fixed follow-up");
  if (milestone >= accrualDuration + followupTime)
    throw std::invalid_argument(
        "milestone cannot exceed accrualDuration + followupTime");

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
  std::vector<double> rmst1(k), rmst2(k), vrmst1(k), vrmst2(k);
  std::vector<double> rmstDiff(k), vrmstDiff(k), information(k), rmstDiffZ(k);

  for (size_t j = 0; j < k; ++j) {
    DataFrameCpp df = rmstat1cpp(
      time[j], milestone, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup);

    subjects[j]     = extract_sum(df, "subjects");
    nevents[j]      = extract_sum(df, "nevents");
    nevents1[j]     = extract_sum(df, "nevents1");
    nevents2[j]     = extract_sum(df, "nevents2");
    ndropouts[j]    = extract_sum(df, "ndropouts");
    ndropouts1[j]   = extract_sum(df, "ndropouts1");
    ndropouts2[j]   = extract_sum(df, "ndropouts2");
    nmilestone[j]   = extract_sum(df, "nmilestone");
    nmilestone1[j]  = extract_sum(df, "nmilestone1");
    nmilestone2[j]  = extract_sum(df, "nmilestone2");

    const auto& s1  = df.get<double>("rmst1");
    const auto& s2  = df.get<double>("rmst2");
    const auto& vs1 = df.get<double>("vrmst1");
    const auto& vs2 = df.get<double>("vrmst2");
    for (size_t h = 0; h < nstrata; ++h) {
      rmst1[j]  += stratumFraction[h] * s1[h];
      rmst2[j]  += stratumFraction[h] * s2[h];
      vrmst1[j] += stratumFraction[h] * stratumFraction[h] * vs1[h];
      vrmst2[j] += stratumFraction[h] * stratumFraction[h] * vs2[h];
    }

    rmstDiff[j]    = rmst1[j] - rmst2[j];
    vrmstDiff[j]   = vrmst1[j] + vrmst2[j];
    information[j] = 1.0 / vrmstDiff[j];
    rmstDiffZ[j]   = rmstDiff[j] * std::sqrt(information[j]);
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
  out.push_back(std::move(rmst1), "rmst1");
  out.push_back(std::move(rmst2), "rmst2");
  out.push_back(std::move(rmstDiff), "rmstDiff");
  out.push_back(std::move(vrmst1), "vrmst1");
  out.push_back(std::move(vrmst2), "vrmst2");
  out.push_back(std::move(vrmstDiff), "vrmstDiff");
  out.push_back(std::move(information), "information");
  out.push_back(std::move(rmstDiffZ), "rmstDiffZ");

  return out;
}


//' @title Stratified Difference in Restricted Mean Survival Times
//' @description Obtains the stratified restricted mean survival times
//' and difference in restricted mean survival times at given calendar
//' times.
//'
//' @param time A vector of calendar times for data cut.
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
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
//' * \code{rmst1}: The restricted mean survival time for the treatment
//'   group.
//'
//' * \code{rmst2}: The restricted mean survival time for the control group.
//'
//' * \code{rmstDiff}: The difference in restricted mean survival times,
//'   i.e., \code{rmst1 - rmst2}.
//'
//' * \code{vrmst1}: The variance for \code{rmst1}.
//'
//' * \code{vrmst2}: The variance for \code{rmst2}.
//'
//' * \code{vrmstDiff}: The variance for \code{rmstDiff}.
//'
//' * \code{information}: The information for \code{rmstDiff}, equal to
//'   \code{1/vrmstDiff}.
//'
//' * \code{rmstDiffZ}: The Z-statistic value, i.e.,
//'   \code{rmstDiff/sqrt(vrmstDiff)}.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' rmstat(time = c(22, 40),
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
Rcpp::DataFrame rmstat(
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

  DataFrameCpp df = rmstatcpp(
    time1, milestone, allocationRatioPlanned,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup);

  return Rcpp::wrap(df);
}


ListCpp rmpowercpp(
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
    const double rmstDiffH0,
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

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1))
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
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
      bsf == "sfkd" || bsf == "sfhsd" || bsf == "none")) {
    throw std::invalid_argument("Invalid value for typeBetaSpending");
  }
  if ((bsf == "sfkd" || bsf == "sfhsd") && std::isnan(parameterBetaSpending)) {
    throw std::invalid_argument("Missing value for parameterBetaSpending");
  }
  if (bsf == "sfkd" && parameterBetaSpending <= 0.0) {
    throw std::invalid_argument("parameterBetaSpending must be positive for sfKD");
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
  if (!std::isnan(studyDuration) && studyDuration < accrualDuration)
    throw std::invalid_argument(
        "studyDuration must be greater than or equal to accrualDuration");
  if (!std::isnan(studyDuration) && studyDuration > accrualDuration + followupTime)
    throw std::invalid_argument(
        "studyDuration cannot exceed accrualDuration + followupTime");
  if (milestone >= accrualDuration + followupTime)
    throw std::invalid_argument(
        "milestone must be less than accrualDuration + followupTime");
  if (fixedFollowup && milestone > followupTime)
    throw std::invalid_argument(
        "milestone cannot exceed followupTime for fixed follow-up");
  if (!std::isnan(studyDuration) && milestone >= studyDuration)
    throw std::invalid_argument("milestone cannot exceed studyDuration");

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

  // expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambda1x = expand_stratified(lambda1, nstrata, nintv, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintv, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintv, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintv, "gamma2");


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
  std::vector<double> I(K);

  // ---- compute maxInformation and theta via rmstat at study end ----
  DataFrameCpp rm_end = rmstat1cpp(
    studyDuration1, milestone, allocationRatioPlanned,
    accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction,
    lambda1x, lambda2x, gamma1x, gamma2x,
    accrualDuration, followupTime, fixedFollowup);

  const auto& s1  = rm_end.get<double>("rmst1");
  const auto& s2  = rm_end.get<double>("rmst2");
  const auto& vs1 = rm_end.get<double>("vrmst1");
  const auto& vs2 = rm_end.get<double>("vrmst2");
  double rmst1 = 0.0, rmst2 = 0.0, vrmst1 = 0.0, vrmst2 = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    rmst1  += stratumFraction[h] * s1[h];
    rmst2  += stratumFraction[h] * s2[h];
    vrmst1 += stratumFraction[h] * stratumFraction[h] * vs1[h];
    vrmst2 += stratumFraction[h] * stratumFraction[h] * vs2[h];
  }
  double rmstDiff = rmst1 - rmst2;
  double vrmstDiff = vrmst1 + vrmst2;
  double maxInformation = 1.0 / vrmstDiff;

  std::vector<double> theta(K, rmstDiff - rmstDiffH0);

  I[K - 1]           = maxInformation;
  time[K - 1]        = studyDuration1;
  nsubjects[K - 1]   = extract_sum(rm_end, "subjects");
  nsubjects1[K - 1]  = phi * nsubjects[K - 1];
  nsubjects2[K - 1]  = (1.0 - phi) * nsubjects[K - 1];
  nevents[K - 1]     = extract_sum(rm_end, "nevents");
  nevents1[K - 1]    = extract_sum(rm_end, "nevents1");
  nevents2[K - 1]    = extract_sum(rm_end, "nevents2");
  ndropouts[K - 1]   = extract_sum(rm_end, "ndropouts");
  ndropouts1[K - 1]  = extract_sum(rm_end, "ndropouts1");
  ndropouts2[K - 1]  = extract_sum(rm_end, "ndropouts2");
  nmilestone[K - 1]  = extract_sum(rm_end, "nmilestone");
  nmilestone1[K - 1] = extract_sum(rm_end, "nmilestone1");
  nmilestone2[K - 1] = extract_sum(rm_end, "nmilestone2");

  // ---- compute information, time, and other stats at interim analyses ----
  for (size_t i = 0; i < K - 1; ++i) {
    double information1 = maxInformation * infoRates[i];
    I[i] = information1;

    // solve for analysis time where total information equals information1
    auto g = [&](double t)->double {
      DataFrameCpp rm1 = rmstat1cpp(
        t, milestone, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrualDuration, followupTime, fixedFollowup);

      const auto& vs1 = rm1.get<double>("vrmst1");
      const auto& vs2 = rm1.get<double>("vrmst2");
      double vrmst1 = 0.0, vrmst2 = 0.0;
      for (size_t h = 0; h < nstrata; ++h) {
        vrmst1 += stratumFraction[h] * stratumFraction[h] * vs1[h];
        vrmst2 += stratumFraction[h] * stratumFraction[h] * vs2[h];
      }
      double vrmstDiff = vrmst1 + vrmst2;
      return 1.0 / vrmstDiff - information1;
    };

    time[i] = brent(g, milestone + 0.001, studyDuration1, 1e-6);

    DataFrameCpp rm_i = rmstat1cpp(
      time[i], milestone, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup);

    nsubjects[i]   = extract_sum(rm_i, "subjects");
    nsubjects1[i]  = phi * nsubjects[i];
    nsubjects2[i]  = (1.0 - phi) * nsubjects[i];
    nevents[i]     = extract_sum(rm_i, "nevents");
    nevents1[i]    = extract_sum(rm_i, "nevents1");
    nevents2[i]    = extract_sum(rm_i, "nevents2");
    ndropouts[i]   = extract_sum(rm_i, "ndropouts");
    ndropouts1[i]  = extract_sum(rm_i, "ndropouts1");
    ndropouts2[i]  = extract_sum(rm_i, "ndropouts2");
    nmilestone[i]  = extract_sum(rm_i, "nmilestone");
    nmilestone1[i] = extract_sum(rm_i, "nmilestone");
    nmilestone2[i] = extract_sum(rm_i, "nmilestone");
  }

  // ---- compute exit probabilities ----
  ListCpp exit_probs;
  if (!missingFutilityBounds || bsf == "none" || K == 1) {
    exit_probs = exitprobcpp(critValues, futBounds, theta, I);
  } else {
    std::vector<double> w(K, 1.0);
    auto gp = getPower(alpha1, kMax, critValues, theta, I, bsf,
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
    expectedInformation += ptotal[i] * I[i];
  }

  std::vector<double> cpu(K), cpl(K);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());
  double overallReject = cpu[K-1];

  // efficacy/futility implied rmstDiff boundaries at each stage
  std::vector<double> rdu(K), rdl(K);
  for (size_t i = 0; i < K; ++i) {
    rdu[i] = rmstDiffH0 + critValues[i] / std::sqrt(I[i]);
    rdl[i] = rmstDiffH0 + futBounds[i] / std::sqrt(I[i]);
    if (critValues[i] == 6.0) { rdu[i] = NaN; effStopping[i] = 0; }
    if (futBounds[i] == -6.0) { rdl[i] = NaN; futStopping[i] = 0; }
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
  overallResults.push_back(rmstDiffH0, "rmstDiffH0");
  overallResults.push_back(rmst1, "rmst1");
  overallResults.push_back(rmst2, "rmst2");
  overallResults.push_back(rmstDiff, "rmstDiff");

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
  byStageResults.push_back(std::move(rdu), "efficacyRmstDiff");
  byStageResults.push_back(std::move(rdl), "futilityRmstDiff");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(futilityP), "futilityP");
  byStageResults.push_back(std::move(I), "information");
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


//' @title Power for Difference in Restricted Mean Survival Times
//' @description Estimates the power for testing the difference in
//' restricted mean survival times in a two-sample survival design.
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
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffH0 The difference in restricted mean survival times
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
//' @return An S3 class \code{rmpower} object with 4 components:
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
//'     - \code{rmstDiffH0}: The difference in restricted mean survival
//'       times under the null hypothesis.
//'
//'     - \code{rmst1}: The restricted mean survival time for the
//'       treatment group.
//'
//'     - \code{rmst2}: The restricted mean survival time for the
//'       control group.
//'
//'     - \code{rmstDiff}: The difference in restricted mean survival times,
//'       equal to \code{rmst1 - rmst2}.
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
//'     - \code{efficacyRmstDiff}: The efficacy boundaries on the restricted
//'       mean survival time difference scale.
//'
//'     - \code{futilityRmstDiff}: The futility boundaries on the restricted
//'       mean survival time difference scale.
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
//' rmpower(kMax = 2, informationRates = c(0.8, 1),
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
Rcpp::List rmpower(
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
    const double rmstDiffH0 = 0,
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

  ListCpp out = rmpowercpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, rmstDiffH0, allocationRatioPlanned,
    accrualT, accrualInt, piecewiseTime, stratumFrac,
    lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, studyDuration);

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "rmpower";
  return result;
}


ListCpp rmsamplesizecpp(
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
    const double rmstDiffH0,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    std::vector<double>& accrualIntensity,
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
      milestone >= accrualDuration + followupTime)
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

  // expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambda1x = expand_stratified(lambda1, nstrata, nintv, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintv, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintv, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintv, "gamma2");


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

  // RMST under H1
  double rmst1 = 0.0, rmst2 = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    double s1 = rmstcpp(0, milestone, piecewiseSurvivalTime, lambda1x[h]);
    double s2 = rmstcpp(0, milestone, piecewiseSurvivalTime, lambda2x[h]);
    rmst1 += stratumFraction[h] * s1;
    rmst2 += stratumFraction[h] * s2;
  }
  double theta1 = (rmst1 - rmst2 - rmstDiffH0);
  std::vector<double> theta(K, theta1);

  // --- determine target maxInformation from group sequential design ---
  ListCpp design = getDesigncpp(
    beta, NaN, theta1, kMax, infoRates,
    effStopping, futStopping,
    critValues, alpha, asf,
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
    DataFrameCpp rm1 = rmstat1cpp(
      t, milestone, allocationRatioPlanned,
      accrualTime, accrInt,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrDur, fu, fixedFollowup);

    const auto& vs1 = rm1.get<double>("vrmst1");
    const auto& vs2 = rm1.get<double>("vrmst2");
    double vrmst1 = 0.0, vrmst2 = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      vrmst1 += stratumFraction[h] * stratumFraction[h] * vs1[h];
      vrmst2 += stratumFraction[h] * stratumFraction[h] * vs2[h];
    }
    double vrmstDiff = vrmst1 + vrmst2;
    return 1.0 / vrmstDiff - maxInformation;
  };


  // when information at zero follow-up (end of enrollment) already exceeds target,
  // set followupTime = 0 and find minimal accrualDuration achieving target
  bool curtailed = false;

  // when information at maximum follow-up (milestone) is still below target,
  // increase accrualDuration to achieve target
  bool expanded = false;

  // Root-finding function to solve for the unknown design parameter to achieve target
  std::function<double(double)> f_root;

  if (unknown == ACC_DUR) {
    f_root = [&](double x)->double {
      return info_minus_target_H1(x + followupTime, x, followupTime, accrualIntensity);
    };
  } else if (unknown == FUP_TIME) {
    if (!fixedFollowup && accrualDuration > milestone + 0.001 &&
        info_minus_target_H1(accrualDuration, accrualDuration, 0.0,
                             accrualIntensity) > 0) {
      std::clog << "WARNING: Information at zero follow-up (end of enrollment) "
                   "already exceeds target. Setting followupTime = 0 and "
                   "finding minimal accrualDuration.\n";
      f_root = [&](double x)->double {
        return info_minus_target_H1(x, x, 0.0, accrualIntensity);
      };

      curtailed = true;
    } else if (info_minus_target_H1(accrualDuration + milestone, accrualDuration,
                                    milestone, accrualIntensity) < 0) {
      std::clog << "WARNING: The required information cannot be attained by "
                   "increasing followupTime alone. accrualDuration is also "
                   "increased to attain the required information.\n";
      f_root = [&](double x)->double {
        return info_minus_target_H1(x + milestone, x, milestone, accrualIntensity);
      };

      expanded = true;
    } else {
      f_root = [&](double x)->double {
        return info_minus_target_H1(accrualDuration + x, accrualDuration, x,
                                    accrualIntensity);
      };
    }
  } else {
    f_root = [&](double m)->double {
      std::vector<double> scaled = accrualIntensity;
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
    for (double &v : accrualIntensity) v *= solution;
  }

  double studyDuration = accrualDuration + followupTime;

  // --- rounding to integer N (adjust intensity or accrualDuration) ---
  if (rounding) {
    double n0 = accrual1(studyDuration, accrualTime, accrualIntensity,
                         accrualDuration);
    double n = std::ceil(n0 - 1.0e-12);

    if (n - n0 > 1e-6) {
      if (unknown == ACC_INT) { // scale intensity to hit integer n
        double mult = n / n0;
        for (double& v : accrualIntensity) v *= mult;
      } else { // adjust accrualDuration to hit integer n
        accrualDuration = getAccrualDurationFromN1(n, accrualTime, accrualIntensity);
      }

      if (!fixedFollowup) {
        // variable follow-up: adjust follow-up time to match maxInformation
        auto h = [&](double x)->double {
          return info_minus_target_H1(accrualDuration + x, accrualDuration, x,
                                      accrualIntensity);
        };
        double lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // fixed follow-up: adjust studyDuration by extending post-accrual time
        auto h = [&](double x)->double {
          return info_minus_target_H1(accrualDuration + x, accrualDuration,
                                      followupTime, accrualIntensity);
        };
        double lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
        double extra = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + extra;
      }
    }
  }

  // --- Results under H1 ---
  ListCpp resultsH1 = rmpowercpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, rmstDiffH0, allocationRatioPlanned,
    accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  // refresh maxInformation from H1 power results    {
  DataFrameCpp overallResultsH1 = resultsH1.get<DataFrameCpp>("overallResults");
  maxInformation = overallResultsH1.get<double>("information")[0];


  // --- Construct lambda1 under H0 (scale lambda2 to match rmstDiffH0) ---
  // Solve for multiplier a s.t. weighted_rmst1(a*lambda2) - rmst2 - rmstDiffH0 == 0
  auto f_rmst = [&](double aval)->double {
    double rmst1 = 0.0;
    for (std::size_t h = 0; h < nstrata; ++h) {
      auto lam1H0 = lambda2x[h]; // build per-stratum scaled lambda
      for (double &v : lam1H0) v *= aval;
      double p = rmstcpp(0, milestone, piecewiseSurvivalTime, lam1H0);
      rmst1 += stratumFraction[h] * p;
    }
    return rmst1 - rmst2 - rmstDiffH0;
  };

  double lo = 0.001, hi = 2.0;
  double flo = f_rmst(lo);
  double fhi = f_rmst(hi);
  int it = 0;
  while (flo * fhi > 0.0 && it < 60) {
    lo = hi;
    flo = fhi;
    hi *= 2.0;
    fhi = f_rmst(hi);
  }
  if (flo * fhi > 0.0) throw std::runtime_error(
      "Unable to bracket root for lambda1H0 under H0; check inputs");

  auto f_for_brent2 = [&](double x)->double {
    if (x == lo) return flo;
    if (x == hi) return fhi;
    return f_rmst(x);
  };

  double multH0 = brent(f_for_brent2, lo, hi, 1.0e-6);

  std::vector<double> lambda1H0 = lambda2;
  for (double &v : lambda1H0) v *= multH0;
  auto lambda1H0x = expand_stratified(lambda1H0, nstrata, nintv, "lambda1H0");

  // Helper: compute information under H0 given accrualDuration (accrDur),
  // followupTime (fu), and accrualIntensity (accrInt).
  auto info_minus_target_H0 = [&](double t, double accrDur, double fu,
                                  const std::vector<double>& accrInt) {
    DataFrameCpp rm1 = rmstat1cpp(
      t, milestone, allocationRatioPlanned,
      accrualTime, accrInt,
      piecewiseSurvivalTime, stratumFraction,
      lambda1H0x, lambda2x, gamma1x, gamma2x,
      accrDur, fu, fixedFollowup);

    const auto& vs1 = rm1.get<double>("vrmst1");
    const auto& vs2 = rm1.get<double>("vrmst2");
    double vrmst1 = 0.0, vrmst2 = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      vrmst1 += stratumFraction[h] * stratumFraction[h] * vs1[h];
      vrmst2 += stratumFraction[h] * stratumFraction[h] * vs2[h];
    }
    double vrmstDiff = vrmst1 + vrmst2;
    return 1.0 / vrmstDiff - maxInformation;
  };

  if (!fixedFollowup) {
    auto h_follow = [&](double x)->double {
      return info_minus_target_H0(accrualDuration + x, accrualDuration, x,
                                  accrualIntensity);
    };

    if (accrualDuration > milestone + 0.001 && h_follow(0.0) > 0.0) {
      std::clog << "WARNING: Information at zero follow-up (end of enrollment) "
                   "already exceeds target. Setting followupTime = 0 and "
                   "finding minimal accrualDuration.\n";
      auto h_accr = [&](double x)->double {
        return info_minus_target_H0(x, x, 0.0, accrualIntensity);
      };

      double lo = milestone + 0.001;
      double hi = accrualDuration;
      accrualDuration = brent(h_accr, lo, hi, 1.0e-6);
      followupTime = 0.0;
      studyDuration = accrualDuration;
    } else if (h_follow(followupTime) < 0.0) {
      std::clog << "WARNING: The required information cannot be attained by "
                   "increasing followupTime alone. accrualDuration is also "
                   "increased to attain the required information.\n";
      auto h_accr = [&](double x)->double {
        return info_minus_target_H0(x + milestone, x, milestone, accrualIntensity);
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
                                  followupTime, accrualIntensity);
    };

    auto h_accr = [&](double x)->double {
      return info_minus_target_H0(x + followupTime, x, followupTime, accrualIntensity);
    };

    if (accrualDuration > milestone + 0.001 && h_study(0.0) > 0.0) {
      std::clog << "WARNING: Information at zero follow-up (end of enrollment) "
                   "already exceeds target. Decrease accrual druation.\n";
      double lo = milestone + 0.001;
      double hi = accrualDuration;
      accrualDuration = brent(h_accr, lo, hi, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (h_study(followupTime) < 0.0) {
      std::clog << "WARNING: The required information cannot be attained by "
                   "increasing followupTime alone. accrualDuration is also "
                   "increased to attain the required information.\n";
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
  ListCpp resultsH0 = rmpowercpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, rmstDiffH0, allocationRatioPlanned,
    accrualTime, accrualIntensity,
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




//' @title Sample Size for Difference in Restricted Mean Survival Times
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
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffH0 The difference in restricted mean survival times
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
//' * \code{resultsUnderH1}: An S3 class \code{rmpower} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{rmpower} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmpower}}
//'
//' @examples
//' # Example 1: Obtains follow-up time given power, accrual intensity,
//' # and accrual duration for variable follow-up. Of note, the power
//' # reaches the maximum when the follow-up time equals milestone.
//'
//' rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
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
//' rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
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
//' rmsamplesize(beta = 0.2, kMax = 2, informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 100/9*seq(1, 9),
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
Rcpp::List rmsamplesize(
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
    const double rmstDiffH0 = 0,
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

  auto out = rmsamplesizecpp(
    beta, kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    userBeta, milestone, rmstDiffH0, allocationRatioPlanned,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, rounding);

  ListCpp resultsUnderH1 = out.get_list("resultsUnderH1");
  ListCpp resultsUnderH0 = out.get_list("resultsUnderH0");

  Rcpp::List resultsH1 = Rcpp::wrap(resultsUnderH1);
  Rcpp::List resultsH0 = Rcpp::wrap(resultsUnderH0);

  resultsH1.attr("class") = "rmpower";
  resultsH0.attr("class") = "rmpower";

  return Rcpp::List::create(
    Rcpp::Named("resultsUnderH1") = resultsH1,
    Rcpp::Named("resultsUnderH0") = resultsH0
  );
}


ListCpp rmpower1scpp(
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
    const double rmstH0,
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

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1.0))
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");
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
  if (std::isnan(milestone)) throw std::invalid_argument("milestone must be provided");
  if (milestone <= 0.0) throw std::invalid_argument("milestone must be positive");
  if (std::isnan(rmstH0)) throw std::invalid_argument("rmstH0 must be provided");
  if (rmstH0 <= 0.0) throw std::invalid_argument("rmstH0 must be positive");
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
  if (!std::isnan(studyDuration) && studyDuration < accrualDuration)
    throw std::invalid_argument(
        "studyDuration must be greater than or equal to accrualDuration");
  if (!std::isnan(studyDuration) && studyDuration > accrualDuration + followupTime)
    throw std::invalid_argument(
        "studyDuration cannot exceed accrualDuration + followupTime");
  if (milestone >= accrualDuration + followupTime)
    throw std::invalid_argument(
        "milestone must be less than accrualDuration + followupTime");
  if (fixedFollowup && milestone > followupTime)
    throw std::invalid_argument(
        "milestone cannot exceed followupTime for fixed follow-up");
  if (!std::isnan(studyDuration) && milestone >= studyDuration)
    throw std::invalid_argument("milestone cannot exceed studyDuration");

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

  // expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambdax = expand_stratified(lambda, nstrata, nintv, "lambda");
  auto gammax = expand_stratified(gamma, nstrata, nintv, "gamma");


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
  std::vector<double> time(K), I(K);
  std::vector<double> nsubjects(K), nevents(K), ndropouts(K), nmilestone(K);

  // compute information using twin group trick
  std::vector<double> accrualInt2 = accrualIntensity;
  for (double& v : accrualInt2) v *= 2.0;

  DataFrameCpp rm_end = rmstat1cpp(
    studyDuration1, milestone, 1,
    accrualTime, accrualInt2,
    piecewiseSurvivalTime, stratumFraction,
    lambdax, lambdax, gammax, gammax,
    accrualDuration, followupTime, fixedFollowup);

  const auto& s  = rm_end.get<double>("rmst1");
  const auto& vs = rm_end.get<double>("vrmst1");
  double rmst = 0.0, vrmst = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    rmst  += stratumFraction[h] * s[h];
    vrmst += stratumFraction[h] * stratumFraction[h] * vs[h];
  }
  double maxInformation = 1.0 / vrmst;
  std::vector<double> theta(K, rmst - rmstH0);

  // set final-stage quantities
  I[K - 1] = maxInformation;
  time[K - 1] = studyDuration1;

  nsubjects[K - 1]  = extract_sum(rm_end, "subjects") * 0.5; // half because twin-group
  nevents[K - 1]    = extract_sum(rm_end, "nevents1");
  ndropouts[K - 1]  = extract_sum(rm_end, "ndropouts1");
  nmilestone[K - 1] = extract_sum(rm_end, "nmilestone1");

  // Compute interim analysis times and stagewise counts
  for (size_t i = 0; i < K - 1; ++i) {
    double information1 = maxInformation * infoRates[i];
    I[i] = information1;

    // function to match predicted information at candidate time
    auto g = [&](double t)->double {
      DataFrameCpp rm1 = rmstat1cpp(
        t, milestone, 1,
        accrualTime, accrualInt2,
        piecewiseSurvivalTime, stratumFraction,
        lambdax, lambdax, gammax, gammax,
        accrualDuration, followupTime, fixedFollowup);

      const auto& vs = rm1.get<double>("vrmst1");
      double vrmst = 0.0;
      for (size_t h = 0; h < nstrata; ++h) {
        vrmst += stratumFraction[h] * stratumFraction[h] * vs[h];
      }
      return 1.0 / vrmst - information1;
    };

    // bracket for brent: between slightly after milestone and studyDuration1
    time[i] = brent(g, milestone + 0.001, studyDuration1, 1e-6);

    // get counts at time[i]
    DataFrameCpp rm_i = rmstat1cpp(
      time[i], milestone, 1,
      accrualTime, accrualInt2,
      piecewiseSurvivalTime, stratumFraction,
      lambdax, lambdax, gammax, gammax,
      accrualDuration, followupTime, fixedFollowup);

    nsubjects[i]  = 0.5 * extract_sum(rm_i, "subjects"); // half because twin-group
    nevents[i]    = extract_sum(rm_i, "nevents1");
    ndropouts[i]  = extract_sum(rm_i, "ndropouts1");
    nmilestone[i] = extract_sum(rm_i, "nmilestone1");
  }

  // compute exit probabilities
  ListCpp exit_probs;
  if (!missingFutilityBounds || bsf == "none" || K == 1) {
    exit_probs = exitprobcpp(critValues, futBounds, theta, I);
  } else {
    std::vector<double> w(K, 1.0);
    auto gp = getPower(alpha1, kMax, critValues, theta, I, bsf,
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
    expectedInformation += ptotal[i] * I[i];
  }

  std::vector<double> cpu(K), cpl(K);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());
  double overallReject = cpu[K-1];

  // implied rmst boundaries
  std::vector<double> rmstu(K), rmstl(K);
  for (size_t i = 0; i < K; ++i) {
    rmstu[i] = rmstH0 + critValues[i] / std::sqrt(I[i]);
    rmstl[i] = rmstH0 + futBounds[i] / std::sqrt(I[i]);
    if (critValues[i] == 6.0) { rmstu[i] = NaN; effStopping[i] = 0; }
    if (futBounds[i] == -6.0) { rmstl[i] = NaN; futStopping[i] = 0; }
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
  overallResults.push_back(rmstH0, "rmstH0");
  overallResults.push_back(rmst, "rmst");

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
  byStageResults.push_back(std::move(rmstu), "efficacyRmst");
  byStageResults.push_back(std::move(rmstl), "futilityRmst");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(futilityP), "futilityP");
  byStageResults.push_back(std::move(I), "information");
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


//' @title Power for One-Sample Restricted Mean Survival Time
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
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstH0 The restricted mean survival time under the null
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
//' @return An S3 class \code{rmpower1s} object with 3 components:
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
//'     - \code{milestone}: The milestone time to calculate the restricted
//'       mean survival time.
//'
//'     - \code{rmstH0}: The restricted mean survival time under the null
//'       hypothesis.
//'
//'     - \code{rmst}: The restricted mean survival time under the
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
//'     - \code{efficacyRmst}: The efficacy boundaries on the restricted
//'       mean survival time.
//'
//'     - \code{futilityRmst}: The futility boundaries on the restricted
//'       mean survival time.
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
//' @seealso \code{\link{rmstat}}
//'
//' @examples
//'
//' rmpower1s(kMax = 2, informationRates = c(0.8, 1),
//'           alpha = 0.025, typeAlphaSpending = "sfOF",
//'           milestone = 18, rmstH0 = 10,
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
Rcpp::List rmpower1s(
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
    const double rmstH0 = NA_REAL,
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

  ListCpp out = rmpower1scpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, rmstH0,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    lam, gam,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, studyDuration);

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "rmpower1s";
  return result;
}


ListCpp rmsamplesize1scpp(
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
    const double rmstH0,
    const std::vector<double>& accrualTime,
    std::vector<double>& accrualIntensity,
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
  if (std::isnan(rmstH0)) throw std::invalid_argument("rmstH0 must be provided");
  if (rmstH0 <= 0.0) throw std::invalid_argument("rmstH0 must be positive");
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

  // expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambdax = expand_stratified(lambda, nstrata, nintv, "lambda");
  auto gammax = expand_stratified(gamma, nstrata, nintv, "gamma");


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

  // RMST under H1
  double rmst = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    double p = rmstcpp(0, milestone, piecewiseSurvivalTime, lambdax[h]);
    rmst += stratumFraction[h] * p;
  }
  double theta1 = rmst - rmstH0;
  std::vector<double> theta(K, theta1);

  // get design under H1 (delegates much work)
  ListCpp design = getDesigncpp(
    beta, NaN, theta1, kMax, infoRates,
    effStopping, futStopping,
    critValues, alpha, asf,
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

    DataFrameCpp rm1 = rmstat1cpp(
      t, milestone, 1,
      accrualTime, accrInt2,
      piecewiseSurvivalTime, stratumFraction,
      lambdax, lambdax, gammax, gammax,
      accrDur, fu, fixedFollowup);

    const auto& vs = rm1.get<double>("vrmst1");
    double vrmst = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      vrmst += stratumFraction[h] * stratumFraction[h] * vs[h];
    }
    return 1.0 / vrmst - maxInformation;
  };

  // when information at zero follow-up (end of enrollment) already exceeds target,
  // set followupTime = 0 and find minimal accrualDuration achieving target
  bool curtailed = false;

  // when information at maximum follow-up (milestone) is still below target,
  // increase accrualDuration to achieve target
  bool expanded = false;

  // Root-finding function to solve for the unknown design parameter to achieve target
  std::function<double(double)> f_root;

  if (unknown == ACC_DUR) {
    f_root = [&](double x)->double {
      return info_minus_target_H1(x + followupTime, x, followupTime, accrualIntensity);
    };
  } else if (unknown == FUP_TIME) {
    if (!fixedFollowup && accrualDuration > milestone + 0.001 &&
        info_minus_target_H1(accrualDuration, accrualDuration, 0.0,
                             accrualIntensity) > 0) {
      std::clog << "WARNING: Information at zero follow-up (end of enrollment) "
                   "already exceeds target. Setting followupTime = 0 and "
                   "finding minimal accrualDuration.\n";
      f_root = [&](double x)->double {
        return info_minus_target_H1(x, x, 0.0, accrualIntensity);
      };

      curtailed = true;
    } else if (info_minus_target_H1(accrualDuration + milestone, accrualDuration,
                                    milestone, accrualIntensity) < 0) {
      std::clog << "WARNING: The required information cannot be attained by "
                   "increasing followupTime alone. accrualDuration is also "
                   "increased to attain the required information.\n";
      f_root = [&](double x)->double {
        return info_minus_target_H1(x + milestone, x, milestone, accrualIntensity);
      };

      expanded = true;
    } else {
      f_root = [&](double x)->double {
        return info_minus_target_H1(accrualDuration + x, accrualDuration, x,
                                    accrualIntensity);
      };
    }
  } else {
    f_root = [&](double m)->double {
      std::vector<double> scaled = accrualIntensity;
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
    for (double &v : accrualIntensity) v *= solution;
  }

  double studyDuration = accrualDuration + followupTime;

  // --- rounding to integer N (adjust intensity or accrualDuration) ---
  if (rounding) {
    double n0 = accrual1(studyDuration, accrualTime, accrualIntensity,
                         accrualDuration);
    double n = std::ceil(n0 - 1.0e-12);

    if (n - n0 > 1e-6) {
      if (unknown == ACC_INT) { // scale intensity to hit integer n
        double mult = n / n0;
        for (double& v : accrualIntensity) v *= mult;
      } else { // adjust accrualDuration to hit integer n
        accrualDuration = getAccrualDurationFromN1(n, accrualTime, accrualIntensity);
      }

      if (!fixedFollowup) {
        // variable follow-up: adjust follow-up time to match maxInformation
        auto h = [&](double x)->double {
          return info_minus_target_H1(accrualDuration + x, accrualDuration, x,
                                      accrualIntensity);
        };
        double lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // fixed follow-up: adjust studyDuration by extending post-accrual time
        auto h = [&](double x)->double {
          return info_minus_target_H1(accrualDuration + x, accrualDuration,
                                      followupTime, accrualIntensity);
        };
        double lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
        double extra = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + extra;
      }
    }
  }

  // --- Results under H1 ---
  ListCpp resultsH1 = rmpower1scpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, rmstH0,
    accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction,
    lambda, gamma,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  // refresh maxInformation from H1 power results    {
  DataFrameCpp overallResultsH1 = resultsH1.get<DataFrameCpp>("overallResults");
  maxInformation = overallResultsH1.get<double>("information")[0];

  // --- Construct lambda under H0 (scale lambda to match rmstH0) ---
  // Solve for multiplier a s.t. weighted_rmst(a*lambda) - rmstH0 == 0
  auto f_rmst = [&](double aval)->double {
    double rmst = 0.0;
    for (std::size_t h = 0; h < nstrata; ++h) {
      auto lamH0 = lambdax[h]; // build per-stratum scaled lambda
      for (double &v : lamH0) v *= aval;
      double p = rmstcpp(0, milestone, piecewiseSurvivalTime, lamH0);
      rmst += stratumFraction[h] * p;
    }
    return rmst - rmstH0;
  };

  double lo = 0.001, hi = 2.0;
  double flo = f_rmst(lo);
  double fhi = f_rmst(hi);
  int it = 0;
  while (flo * fhi > 0.0 && it < 60) {
    lo = hi;
    flo = fhi;
    hi *= 2.0;
    fhi = f_rmst(hi);
  }
  if (flo * fhi > 0.0) throw std::runtime_error(
      "Unable to bracket root for lambda under H0; check inputs");

  auto f_for_brent2 = [&](double x)->double {
    if (x == lo) return flo;
    if (x == hi) return fhi;
    return f_rmst(x);
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

    DataFrameCpp rm1 = rmstat1cpp(
      t, milestone, 1,
      accrualTime, accrInt2,
      piecewiseSurvivalTime, stratumFraction,
      lambdaH0x, lambdaH0x, gammax, gammax,
      accrDur, fu, fixedFollowup);

    const auto& vs = rm1.get<double>("vrmst1");
    double vrmst = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      vrmst += stratumFraction[h] * stratumFraction[h] * vs[h];
    }
    return 1.0 / vrmst - maxInformation;
  };

  if (!fixedFollowup) {
    auto h_follow = [&](double x)->double {
      return info_minus_target_H0(accrualDuration + x, accrualDuration, x,
                                  accrualIntensity);
    };

    if (accrualDuration > milestone + 0.001 && h_follow(0.0) > 0.0) {
      std::clog << "WARNING: Information at zero follow-up (end of enrollment) "
                   "already exceeds target. Setting followupTime = 0 and "
                   "finding minimal accrualDuration.\n";
      auto h_accr = [&](double x)->double {
        return info_minus_target_H0(x, x, 0.0, accrualIntensity);
      };

      double lo = milestone + 0.001;
      double hi = accrualDuration;
      accrualDuration = brent(h_accr, lo, hi, 1.0e-6);
      followupTime = 0.0;
      studyDuration = accrualDuration;
    } else if (h_follow(followupTime) < 0.0) {
      std::clog << "WARNING: The required information cannot be attained by "
                   "increasing followupTime alone. accrualDuration is also "
                   "increased to attain the required information.\n";
      auto h_accr = [&](double x)->double {
        return info_minus_target_H0(x + milestone, x, milestone, accrualIntensity);
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
                                  followupTime, accrualIntensity);
    };

    auto h_accr = [&](double x)->double {
      return info_minus_target_H0(x + followupTime, x, followupTime, accrualIntensity);
    };

    if (accrualDuration > milestone + 0.001 && h_study(0.0) > 0.0) {
      std::clog << "WARNING: Information at zero follow-up (end of enrollment) "
                   "already exceeds target. Decrease accrual druation.\n";
      double lo = milestone + 0.001;
      double hi = accrualDuration;
      accrualDuration = brent(h_accr, lo, hi, 1.0e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (h_study(followupTime) < 0.0) {
      std::clog << "WARNING: The required information cannot be attained by "
                   "increasing followupTime alone. accrualDuration is also "
                   "increased to attain the required information.\n";
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
  ListCpp resultsH0 = rmpower1scpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, typeBetaSpending, parameterBetaSpending,
    milestone, rmstH0,
    accrualTime, accrualIntensity,
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


//' @title Sample Size for One-Sample Restricted Mean Survival Time
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
//' @param milestone The milestone time at which to calculate the
//'   restricted survival time.
//' @param rmstH0 The restricted mean survival time under the null
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
//' * \code{resultsUnderH1}: An S3 class \code{rmpower1s} object under the
//'   alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{rmpower1s} object under the
//'   null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmpower1s}}
//'
//' @examples
//' # Example 1: Obtains follow-up duration given power, accrual intensity,
//' # and accrual duration for variable follow-up
//'
//' rmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, rmstH0 = 10,
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
//' rmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, rmstH0 = 10,
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
//' rmsamplesize1s(beta = 0.2, kMax = 2,
//'                informationRates = c(0.8, 1),
//'                alpha = 0.025, typeAlphaSpending = "sfOF",
//'                milestone = 18, rmstH0 = 10,
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
Rcpp::List rmsamplesize1s(
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
    const double rmstH0 = NA_REAL,
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

  ListCpp out = rmsamplesize1scpp(
    beta, kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    userBeta, milestone, rmstH0,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    lam, gam,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, rounding);

  ListCpp resultsUnderH1 = out.get_list("resultsUnderH1");
  ListCpp resultsUnderH0 = out.get_list("resultsUnderH0");

  Rcpp::List resultsH1 = Rcpp::wrap(resultsUnderH1);
  Rcpp::List resultsH0 = Rcpp::wrap(resultsUnderH0);

  resultsH1.attr("class") = "rmpower1s";
  resultsH0.attr("class") = "rmpower1s";

  return Rcpp::List::create(
    Rcpp::Named("resultsUnderH1") = resultsH1,
    Rcpp::Named("resultsUnderH0") = resultsH0
  );
}


ListCpp rmpowerequivcpp(
    const int kMax,
    const std::vector<double>& informationRates,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const double milestone,
    const double rmstDiffLower,
    const double rmstDiffUpper,
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

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1))
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
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

  bool missingCriticalValues = !none_na(criticalValues);
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
  if (std::isnan(milestone)) throw std::invalid_argument("milestone must be provided");
  if (milestone <= 0.0) throw std::invalid_argument("milestone must be positive");
  if (std::isnan(rmstDiffLower))
    throw std::invalid_argument("rmstDiffLower must be provided");
  if (std::isnan(rmstDiffUpper))
    throw std::invalid_argument("rmstDiffUpper must be provided");
  if (rmstDiffLower >= rmstDiffUpper)
    throw std::invalid_argument("rmstDiffLower must be less than rmstDiffUpper");
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
  if (!std::isnan(studyDuration) && studyDuration < accrualDuration)
    throw std::invalid_argument(
        "studyDuration must be greater than or equal to accrualDuration");
  if (!std::isnan(studyDuration) && studyDuration > accrualDuration + followupTime)
    throw std::invalid_argument(
        "studyDuration cannot exceed accrualDuration + followupTime");
  if (milestone >= accrualDuration + followupTime)
    throw std::invalid_argument(
        "milestone must be less than accrualDuration + followupTime");
  if (fixedFollowup && milestone > followupTime)
    throw std::invalid_argument(
        "milestone cannot exceed followupTime for fixed follow-up");
  if (!std::isnan(studyDuration) && milestone >= studyDuration)
    throw std::invalid_argument("milestone cannot exceed studyDuration");

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

  // expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambda1x = expand_stratified(lambda1, nstrata, nintv, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintv, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintv, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintv, "gamma2");


  // obtain criticalValues if missing
  std::vector<double> u(K), l(K, -6.0), zero(K, 0.0);
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
      for (size_t i = 0; i < K - 1; ++i) u[i] = criticalValues[i];

      auto f = [&](double aval)->double {
        u[K-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[K-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      std::vector<unsigned char> effStopping(K, 1);
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
    }
  }

  std::vector<double> li(K, -6.0), ui(K, 6.0);
  ListCpp probs = exitprobcpp(critValues, li, zero, infoRates);
  auto v = probs.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cumAlphaSpent(K);
  std::partial_sum(v.begin(), v.end(), cumAlphaSpent.begin());

  std::vector<double> efficacyP(K);
  for (size_t i = 0; i < K; ++i) {
    efficacyP[i] = 1.0 - boost_pnorm(critValues[i]);
  }

  // --- timing, events, information -----------
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  double studyDuration1 = studyDuration;
  if (!fixedFollowup || std::isnan(studyDuration))
    studyDuration1 = accrualDuration + followupTime;

  std::vector<double> time(K), nsubjects(K), nsubjects1(K), nsubjects2(K);
  std::vector<double> nevents(K), nevents1(K), nevents2(K);
  std::vector<double> ndropouts(K), ndropouts1(K), ndropouts2(K);
  std::vector<double> nmilestone(K), nmilestone1(K), nmilestone2(K);

  // ---- compute maxInformation and theta via rmstat at study end ----
  DataFrameCpp rm_end = rmstat1cpp(
    studyDuration1, milestone, allocationRatioPlanned,
    accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction,
    lambda1x, lambda2x, gamma1x, gamma2x,
    accrualDuration, followupTime, fixedFollowup);

  const auto& s1  = rm_end.get<double>("rmst1");
  const auto& s2  = rm_end.get<double>("rmst2");
  const auto& vs1 = rm_end.get<double>("vrmst1");
  const auto& vs2 = rm_end.get<double>("vrmst2");
  double rmst1 = 0.0, rmst2 = 0.0, vrmst1 = 0.0, vrmst2 = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    rmst1  += stratumFraction[h] * s1[h];
    rmst2  += stratumFraction[h] * s2[h];
    vrmst1 += stratumFraction[h] * stratumFraction[h] * vs1[h];
    vrmst2 += stratumFraction[h] * stratumFraction[h] * vs2[h];
  }
  double rmstDiff = rmst1 - rmst2;
  double vrmstDiff = vrmst1 + vrmst2;
  double maxInformation = 1.0 / vrmstDiff;

  // stage information and analysis times
  std::vector<double> I(K);
  for (size_t i = 0; i < K; ++i) I[i] = maxInformation * infoRates[i];

  time[K - 1]        = studyDuration1;
  nsubjects[K - 1]   = extract_sum(rm_end, "subjects");
  nsubjects1[K - 1]  = phi * nsubjects[K - 1];
  nsubjects2[K - 1]  = (1.0 - phi) * nsubjects[K - 1];
  nevents[K - 1]     = extract_sum(rm_end, "nevents");
  nevents1[K - 1]    = extract_sum(rm_end, "nevents1");
  nevents2[K - 1]    = extract_sum(rm_end, "nevents2");
  ndropouts[K - 1]   = extract_sum(rm_end, "ndropouts");
  ndropouts1[K - 1]  = extract_sum(rm_end, "ndropouts1");
  ndropouts2[K - 1]  = extract_sum(rm_end, "ndropouts2");
  nmilestone[K - 1]  = extract_sum(rm_end, "nmilestone");
  nmilestone1[K - 1] = extract_sum(rm_end, "nmilestone1");
  nmilestone2[K - 1] = extract_sum(rm_end, "nmilestone2");

  // ---- compute information, time, and other stats at interim analyses ----
  for (size_t i = 0; i < K - 1; ++i) {
    double information1 = maxInformation * infoRates[i];

    // solve for analysis time where total information equals information1
    auto g = [&](double t)->double {
      DataFrameCpp rm1 = rmstat1cpp(
        t, milestone, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrualDuration, followupTime, fixedFollowup);

      const auto& vs1 = rm1.get<double>("vrmst1");
      const auto& vs2 = rm1.get<double>("vrmst2");
      double vrmst1 = 0.0, vrmst2 = 0.0;
      for (size_t h = 0; h < nstrata; ++h) {
        vrmst1 += stratumFraction[h] * stratumFraction[h] * vs1[h];
        vrmst2 += stratumFraction[h] * stratumFraction[h] * vs2[h];
      }
      double vrmstDiff = vrmst1 + vrmst2;
      return 1.0 / vrmstDiff - information1;
    };

    time[i] = brent(g, milestone + 0.001, studyDuration1, 1e-6);

    DataFrameCpp rm_i = rmstat1cpp(
      time[i], milestone, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup);

    nsubjects[i]   = extract_sum(rm_i, "subjects");
    nsubjects1[i]  = phi * nsubjects[i];
    nsubjects2[i]  = (1.0 - phi) * nsubjects[i];
    nevents[i]     = extract_sum(rm_i, "nevents");
    nevents1[i]    = extract_sum(rm_i, "nevents1");
    nevents2[i]    = extract_sum(rm_i, "nevents2");
    ndropouts[i]   = extract_sum(rm_i, "ndropouts");
    ndropouts1[i]  = extract_sum(rm_i, "ndropouts1");
    ndropouts2[i]  = extract_sum(rm_i, "ndropouts2");
    nmilestone[i]  = extract_sum(rm_i, "nmilestone");
    nmilestone1[i] = extract_sum(rm_i, "nmilestone");
    nmilestone2[i] = extract_sum(rm_i, "nmilestone");
  }

  // compute cumulative rejection probabilities under H1
  std::vector<double> sqrtI(K), b(K), a(K);
  for (size_t i = 0; i < K; ++i) {
    sqrtI[i] = std::sqrt(I[i]);
    l[i] = critValues[i] + (rmstDiffLower - rmstDiff) * sqrtI[i];
    u[i] = -critValues[i] + (rmstDiffUpper - rmstDiff) * sqrtI[i];
    b[i] = std::max(l[i], li[i]);
    a[i] = std::min(u[i], ui[i]);
  }

  std::vector<double> cpl(K), cpu(K);
  ListCpp probs1 = exitprobcpp(b, li, zero, I);
  ListCpp probs2 = exitprobcpp(ui, a, zero, I);
  auto v1 = probs1.get<std::vector<double>>("exitProbUpper");
  auto v2 = probs2.get<std::vector<double>>("exitProbLower");
  std::partial_sum(v1.begin(), v1.end(), cpl.begin());
  std::partial_sum(v2.begin(), v2.end(), cpu.begin());

  // index for the first crossing look (0-based)
  size_t kk = K;
  for (size_t i = 0; i < K; ++i) {
    if (l[i] <= u[i]) { kk = i; break; }
  }
  int k = static_cast<int>(kk);

  std::vector<double> cp(K);
  if (k == 0) { // crossing at the first look
    for (size_t i = 0; i < K; ++i) {
      cp[i] = cpl[i] + cpu[i] - 1.0;
    }
  } else {
    std::vector<double> cplx(kk), cpux(kk);
    std::vector l1 = subset(l, 0, k);
    std::vector u1 = subset(u, 0, k);
    std::vector d1 = subset(zero, 0, k);
    std::vector I1 = subset(I, 0, k);
    ListCpp probs = exitprobcpp(l1, u1, d1, I1);
    auto v1x = probs.get<std::vector<double>>("exitProbUpper");
    auto v2x = probs.get<std::vector<double>>("exitProbLower");
    std::partial_sum(v1x.begin(), v1x.end(), cplx.begin());
    std::partial_sum(v2x.begin(), v2x.end(), cpux.begin());
    for (size_t i = 0; i < kk; ++i) {
      cp[i] = cpl[i] + cpu[i] - cplx[i] - cpux[i];
    }
    for (size_t i = kk; i < K; ++i) {
      cp[i] = cpl[i] + cpu[i] - 1.0;
    }
  }

  // incremental rejection probability at each stage
  std::vector<double> rejectPerStage(K);
  rejectPerStage[0] = cp[0];
  for (size_t i = 1; i < K; ++i) {
    rejectPerStage[i] = cp[i] - cp[i-1];
  }

  std::vector<double> q = rejectPerStage;
  if (K > 1) q[K - 1] = 1.0 - cp[K - 2];

  // efficacy rmst diff bounds
  std::vector<double> efficacyRmstDiffLower(K), efficacyRmstDiffUpper(K);
  for (size_t i=0;i<K;++i) {
    double rmstDiffBound = critValues[i] / sqrtI[i];
    efficacyRmstDiffLower[i] = rmstDiffLower + rmstDiffBound;
    efficacyRmstDiffUpper[i] = rmstDiffUpper - rmstDiffBound;
  }

  // cumulative attained alpha under H10 (at rmstDiffLower)
  for (size_t i = 0; i < K; ++i) {
    l[i] = critValues[i];
    u[i] = -critValues[i] + (rmstDiffUpper - rmstDiffLower) * sqrtI[i];
    a[i] = std::min(u[i], ui[i]);
  }
  ListCpp probsH10 = exitprobcpp(ui, a, zero, I);
  auto vH10 = probsH10.get<std::vector<double>>("exitProbLower");
  std::vector<double> cpuH10(K);
  std::partial_sum(vH10.begin(), vH10.end(), cpuH10.begin());
  std::vector<double> cplH10 = cumAlphaSpent;

  std::vector<double> cpH10(K);
  if (k == 0) {
    for (size_t i = 0; i < K; ++i) {
      cpH10[i] = cplH10[i] + cpuH10[i] - 1.0;
    }
  } else {
    std::vector<double> cplH10x(kk), cpuH10x(kk);
    std::vector l1 = subset(l, 0, k);
    std::vector u1 = subset(u, 0, k);
    std::vector d1 = subset(zero, 0, k);
    std::vector I1 = subset(I, 0, k);
    ListCpp probs = exitprobcpp(l1, u1, d1, I1);
    auto v1x = probs.get<std::vector<double>>("exitProbUpper");
    auto v2x = probs.get<std::vector<double>>("exitProbLower");
    std::partial_sum(v1x.begin(), v1x.end(), cplH10x.begin());
    std::partial_sum(v2x.begin(), v2x.end(), cpuH10x.begin());
    for (size_t i = 0; i < kk; ++i) {
      cpH10[i] = cplH10[i] + cpuH10[i] - cplH10x[i] - cpuH10x[i];
    }
    for (size_t i = kk; i < K; ++i) {
      cpH10[i] = cplH10[i] + cpuH10[i] - 1.0;
    }
  }

  // cumulative attained alpha under H20 (at rmstDiffUpper)
  for (size_t i = 0; i < K; ++i) {
    l[i] = critValues[i] + (rmstDiffLower - rmstDiffUpper) * sqrtI[i];
    u[i] = -critValues[i];
    b[i] = std::max(l[i], li[i]);
  }

  ListCpp probsH20 = exitprobcpp(b, li, zero, I);
  auto vH20 = probsH20.get<std::vector<double>>("exitProbUpper");
  std::vector<double> cplH20(K);
  std::partial_sum(vH20.begin(), vH20.end(), cplH20.begin());
  std::vector<double> cpuH20 = cumAlphaSpent;

  std::vector<double> cpH20(K);
  if (k == 0) {
    for (size_t i = 0; i < K; ++i) {
      cpH20[i] = cplH20[i] + cpuH20[i] - 1.0;
    }
  } else {
    std::vector<double> cplH20x(kk), cpuH20x(kk);
    std::vector l1 = subset(l, 0, k);
    std::vector u1 = subset(u, 0, k);
    std::vector d1 = subset(zero, 0, k);
    std::vector I1 = subset(I, 0, k);
    ListCpp probs = exitprobcpp(l1, u1, d1, I1);
    auto v1x = probs.get<std::vector<double>>("exitProbUpper");
    auto v2x = probs.get<std::vector<double>>("exitProbLower");
    std::partial_sum(v1x.begin(), v1x.end(), cplH20x.begin());
    std::partial_sum(v2x.begin(), v2x.end(), cpuH20x.begin());
    for (size_t i = 0; i < kk; ++i) {
      cpH20[i] = cplH20[i] + cpuH20[i] - cplH20x[i] - cpuH20x[i];
    }
    for (size_t i = kk; i < K; ++i) {
      cpH20[i] = cplH20[i] + cpuH20[i] - 1.0;
    }
  }

  double overallReject = cp[K-1];
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
    expectedNumberOfEvents += q[i] * nevents[i];
    expectedNumberOfDropouts += q[i] * ndropouts[i];
    expectedNumberOfSubjects += q[i] * nsubjects[i];
    expectedNumberOfMiles += q[i] * nmilestone[i];
    expectedNumberOfEvents1 += q[i] * nevents1[i];
    expectedNumberOfDropouts1 += q[i] * ndropouts1[i];
    expectedNumberOfSubjects1 += q[i] * nsubjects1[i];
    expectedNumberOfMiles1 += q[i] * nmilestone1[i];
    expectedNumberOfEvents2 += q[i] * nevents2[i];
    expectedNumberOfDropouts2 += q[i] * ndropouts2[i];
    expectedNumberOfSubjects2 += q[i] * nsubjects2[i];
    expectedNumberOfMiles2 += q[i] * nmilestone2[i];
    expectedStudyDuration += q[i] * time[i];
    expectedInformation += q[i] * I[i];
  }

  // Build output
  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha, "alpha");
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
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(milestone, "milestone");
  overallResults.push_back(rmstDiffLower, "rmstDiffLower");
  overallResults.push_back(rmstDiffUpper, "rmstDiffUpper");
  overallResults.push_back(rmst1, "rmst1");
  overallResults.push_back(rmst2, "rmst2");
  overallResults.push_back(rmstDiff, "rmstDiff");
  overallResults.push_back(accrualDuration, "accrualDuration");
  overallResults.push_back(followupTime, "followupTime");
  overallResults.push_back(fixedFollowup, "fixedFollowup");

  // Build byStageResults DataFrameCpp
  DataFrameCpp byStageResults;
  byStageResults.push_back(std::move(infoRates), "informationRates");
  byStageResults.push_back(std::move(critValues), "efficacyBounds");
  byStageResults.push_back(std::move(rejectPerStage), "rejectPerStage");
  byStageResults.push_back(std::move(cp), "cumulativeRejection");
  byStageResults.push_back(std::move(cumAlphaSpent), "cumulativeAlphaSpent");
  byStageResults.push_back(std::move(cpH10), "cumulativeAttainedAlphaH10");
  byStageResults.push_back(std::move(cpH20), "cumulativeAttainedAlphaH20");
  byStageResults.push_back(std::move(nevents), "numberOfEvents");
  byStageResults.push_back(std::move(ndropouts), "numberOfDropouts");
  byStageResults.push_back(std::move(nsubjects), "numberOfSubjects");
  byStageResults.push_back(std::move(nmilestone), "numberOfMilestone");
  byStageResults.push_back(std::move(time), "analysisTime");
  byStageResults.push_back(std::move(efficacyRmstDiffLower), "efficacyRmstDiffLower");
  byStageResults.push_back(std::move(efficacyRmstDiffUpper), "efficacyRmstDiffUpper");
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(I), "information");

  // settings
  ListCpp settings;
  settings.push_back(typeAlphaSpending, "typeAlphaSpending");
  settings.push_back(parameterAlphaSpending, "parameterAlphaSpending");
  settings.push_back(userAlphaSpending, "userAlphaSpending");
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

  // byTreatmentCounts
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

  // assemble final result
  ListCpp result;
  result.push_back(byStageResults, "byStageResults");
  result.push_back(overallResults, "overallResults");
  result.push_back(settings, "settings");
  result.push_back(byTreatmentCounts, "byTreatmentCounts");

  return result;
}


//' @title Power for Equivalence in Restricted Mean Survival Time Difference
//' @description Obtains the power for equivalence in restricted mean
//' survival time difference.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_criticalValues
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffLower The lower equivalence limit of restricted mean
//'   survival time difference.
//' @param rmstDiffUpper The upper equivalence limit of restricted mean
//'   survival time difference.
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
//' @return An S3 class \code{rmpowerequiv} object with 4 components:
//'
//' * \code{overallResults}: A data frame containing the following variables:
//'
//'     - \code{overallReject}: The overall rejection probability.
//'
//'     - \code{alpha}: The overall significance level.
//'
//'     - \code{numberOfEvents}: The total number of events.
//'
//'     - \code{numberOfSubjects}: The total number of subjects.
//'
//'     - \code{studyDuration}: The total study duration.
//'
//'     - \code{information}: The maximum information.
//'
//'     - \code{expectedNumberOfEvents}: The expected number of events.
//'
//'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
//'
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{milestone}: The milestone time relative to randomization.
//'
//'     - \code{rmstDiffLower}: The lower equivalence limit of restricted
//'       mean survival time difference.
//'
//'     - \code{rmstDiffUpper}: The upper equivalence limit of restricted
//'       mean survival time difference.
//'
//'     - \code{rmst1}: The restricted mean survival time for the
//'       treatment group.
//'
//'     - \code{rmst2}: The restricted mean survival time for the
//'       control group.
//'
//'     - \code{rmstDiff}: The restricted mean survival time difference.
//'
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up duration.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//' * \code{byStageResults}: A data frame containing the following variables:
//'
//'     - \code{informationRates}: The information rates.
//'
//'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
//'       each of the two one-sided tests.
//'
//'     - \code{rejectPerStage}: The probability for efficacy stopping.
//'
//'     - \code{cumulativeRejection}: The cumulative probability for efficacy
//'       stopping.
//'
//'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
//'       the two one-sided tests.
//'
//'     - \code{cumulativeAttainedAlphaH10}: The cumulative alpha attained
//'       under \code{H10}.
//'
//'     - \code{cumulativeAttainedAlphaH20}: The cumulative alpha attained
//'       under \code{H20}.
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
//'     - \code{efficacyRmstDiffLower}: The efficacy boundaries on the
//'       restricted mean survival time difference scale for the one-sided
//'       null hypothesis at the lower equivalence limit.
//'
//'     - \code{efficacyRmstDiffUpper}: The efficacy boundaries on the
//'       restricted mean survival time difference scale for the one-sided
//'       null hypothesis at the upper equivalence limit.
//'
//'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
//'       each of the two one-sided tests.
//'
//'     - \code{information}: The cumulative information.
//'
//' * \code{settings}: A list containing the following input parameters:
//'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
//'   \code{userAlphaSpending}, \code{allocationRatioPlanned},
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
//' @seealso \code{\link{rmstat}}
//'
//' @examples
//'
//' rmpowerequiv(kMax = 2, informationRates = c(0.5, 1),
//'              alpha = 0.05, typeAlphaSpending = "sfOF",
//'              milestone = 18,
//'              rmstDiffLower = -2, rmstDiffUpper = 2,
//'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'              accrualIntensity = 29/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              stratumFraction = c(0.2, 0.8),
//'              lambda1 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List rmpowerequiv(
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.05,
    const std::string& typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const double milestone = NA_REAL,
    const double rmstDiffLower = NA_REAL,
    const double rmstDiffUpper = NA_REAL,
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
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  ListCpp out = rmpowerequivcpp(
    kMax, infoRates,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    milestone, rmstDiffLower, rmstDiffUpper, allocationRatioPlanned,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, studyDuration);

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "rmpowerequiv";
  return result;
}


ListCpp rmsamplesizeequivcpp(
    const double beta,
    const int kMax,
    const std::vector<double>& informationRates,
    const std::vector<double>& criticalValues,
    const double alpha,
    const std::string& typeAlphaSpending,
    const double parameterAlphaSpending,
    const std::vector<double>& userAlphaSpending,
    const double milestone,
    const double rmstDiffLower,
    const double rmstDiffUpper,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
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

  // informationRates default
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

  bool missingCriticalValues = !none_na(criticalValues);
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
  if (std::isnan(milestone)) throw std::invalid_argument("milestone must be provided");
  if (milestone <= 0.0) throw std::invalid_argument("milestone must be positive");
  if (std::isnan(rmstDiffLower))
    throw std::invalid_argument("rmstDiffLower must be provided");
  if (std::isnan(rmstDiffUpper))
    throw std::invalid_argument("rmstDiffUpper must be provided");
  if (rmstDiffLower >= rmstDiffUpper)
    throw std::invalid_argument("rmstDiffLower must be less than rmstDiffUpper");
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

  // expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();
  auto lambda1x = expand_stratified(lambda1, nstrata, nintv, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintv, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintv, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintv, "gamma2");


  // obtain criticalValues if missing
  std::vector<double> u(K), l(K, -6.0), zero(K, 0.0);
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
      for (size_t i = 0; i < K - 1; ++i) u[i] = criticalValues[i];

      auto f = [&](double aval)->double {
        u[K-1] = aval;
        ListCpp probs = exitprobcpp(u, l, zero, infoRates);
        auto v = probs.get<std::vector<double>>("exitProbUpper");
        double cpu = std::accumulate(v.begin(), v.end(), 0.0);
        return cpu - alpha;
      };

      critValues[K-1] = brent(f, -5.0, 6.0, 1e-6);
    } else {
      std::vector<unsigned char> effStopping(K, 1);
      critValues = getBoundcpp(kMax, infoRates, alpha, asf,
                               parameterAlphaSpending, userAlphaSpending,
                               spendTime, effStopping);
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

  // RMST under H1
  double rmst1 = 0.0, rmst2 = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    double p1 = rmstcpp(0, milestone, piecewiseSurvivalTime, lambda1x[h]);
    double p2 = rmstcpp(0, milestone, piecewiseSurvivalTime, lambda2x[h]);
    rmst1 += stratumFraction[h] * p1;
    rmst2 += stratumFraction[h] * p2;
  }
  double rmstDiff = rmst1 - rmst2;

  // getDesignEquiv to compute maxInformation and other design items
  ListCpp design = getDesignEquivcpp(
    beta, NaN, rmstDiffLower, rmstDiffUpper, rmstDiff,
    kMax, infoRates, critValues,
    alpha, asf, parameterAlphaSpending,
    userAlphaSpending, spendTime);

  auto overallResults = design.get<DataFrameCpp>("overallResults");
  double maxInformation = overallResults.get<double>("information")[0];

  // Helper: compute information under H1 given accrualDuration (accrDur),
  // followupTime (fu), and accrualIntensity (accrInt).
  auto info_minus_target_H1 = [&](double t, double accrDur, double fu,
                                  const std::vector<double>& accrInt) {
    DataFrameCpp rm1 = rmstat1cpp(
      t, milestone, allocationRatioPlanned,
      accrualTime, accrInt,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrDur, fu, fixedFollowup);

    const auto& vs1 = rm1.get<double>("vrmst1");
    const auto& vs2 = rm1.get<double>("vrmst2");
    double vrmst1 = 0.0, vrmst2 = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      vrmst1 += stratumFraction[h] * stratumFraction[h] * vs1[h];
      vrmst2 += stratumFraction[h] * stratumFraction[h] * vs2[h];
    }
    double vrmstDiff = vrmst1 + vrmst2;
    return 1.0 / vrmstDiff - maxInformation;
  };

  // when information at zero follow-up (end of enrollment) already exceeds target,
  // set followupTime = 0 and find minimal accrualDuration achieving target
  bool curtailed = false;

  // when information at maximum follow-up (milestone) is still below target,
  // increase accrualDuration to achieve target
  bool expanded = false;

  // Root-finding function to solve for the unknown design parameter to achieve target
  std::function<double(double)> f_root;

  if (unknown == ACC_DUR) {
    f_root = [&](double x)->double {
      return info_minus_target_H1(x + followupTime, x, followupTime, accrualIntensity);
    };
  } else if (unknown == FUP_TIME) {
    if (!fixedFollowup && accrualDuration > milestone + 0.001 &&
        info_minus_target_H1(accrualDuration, accrualDuration, 0.0,
                             accrualIntensity) > 0) {
      std::clog << "WARNING: Information at zero follow-up (end of enrollment) "
                   "already exceeds target. Setting followupTime = 0 and "
                   "finding minimal accrualDuration.\n";
      f_root = [&](double x)->double {
        return info_minus_target_H1(x, x, 0.0, accrualIntensity);
      };

      curtailed = true;
    } else if (info_minus_target_H1(accrualDuration + milestone, accrualDuration,
                                    milestone, accrualIntensity) < 0) {
      std::clog << "WARNING: The required information cannot be attained by "
                   "increasing followupTime alone. accrualDuration is also "
                   "increased to attain the required information.\n";
      f_root = [&](double x)->double {
        return info_minus_target_H1(x + milestone, x, milestone, accrualIntensity);
      };

      expanded = true;
    } else {
      f_root = [&](double x)->double {
        return info_minus_target_H1(accrualDuration + x, accrualDuration, x,
                                    accrualIntensity);
      };
    }
  } else {
    f_root = [&](double m)->double {
      std::vector<double> scaled = accrualIntensity;
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
    for (double &v : accrualIntensity) v *= solution;
  }

  double studyDuration = accrualDuration + followupTime;

  // --- rounding to integer N (adjust intensity or accrualDuration) ---
  if (rounding) {
    double n0 = accrual1(studyDuration, accrualTime, accrualIntensity,
                         accrualDuration);
    double n = std::ceil(n0 - 1.0e-12);

    if (n - n0 > 1e-6) {
      if (unknown == ACC_INT) { // scale intensity to hit integer n
        double mult = n / n0;
        for (double& v : accrualIntensity) v *= mult;
      } else { // adjust accrualDuration to hit integer n
        accrualDuration = getAccrualDurationFromN1(n, accrualTime, accrualIntensity);
      }

      if (!fixedFollowup) {
        // variable follow-up: adjust follow-up time to match maxInformation
        auto h = [&](double x)->double {
          return info_minus_target_H1(accrualDuration + x, accrualDuration, x,
                                      accrualIntensity);
        };
        double lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
        followupTime = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + followupTime;
      } else {
        // fixed follow-up: adjust studyDuration by extending post-accrual time
        auto h = [&](double x)->double {
          return info_minus_target_H1(accrualDuration + x, accrualDuration,
                                      followupTime, accrualIntensity);
        };
        double lower = std::max(milestone - accrualDuration, 0.0) + 0.001;
        double extra = brent(h, lower, followupTime, 1.0e-6);
        studyDuration = accrualDuration + extra;
      }
    }
  }

  // call rmpowerequivcpp to compute final results
  ListCpp result = rmpowerequivcpp(
    kMax, infoRates,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    milestone, rmstDiffLower, rmstDiffUpper, allocationRatioPlanned,
    accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction,
    lambda1, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    spendingTime, studyDuration);

  return result;
}


//' @title Sample Size for Equivalence in Restricted Mean Survival Time
//' Difference
//' @description Obtains the sample size for equivalence in restricted
//' mean survival time difference.
//'
//' @param beta The type II error.
//' @inheritParams param_kMax
//' @param informationRates The information rates.
//'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
//' @inheritParams param_criticalValues
//' @param alpha The significance level for each of the two one-sided
//'   tests. Defaults to 0.05.
//' @inheritParams param_typeAlphaSpending
//' @inheritParams param_parameterAlphaSpending
//' @inheritParams param_userAlphaSpending
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffLower The lower equivalence limit of restricted mean
//'   survival time difference.
//' @param rmstDiffUpper The upper equivalence limit of restricted mean
//'   survival time difference.
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
//'   accrualDuration, followupDuration, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size.
//'   Defaults to 1 for sample size rounding.
//'
//' @return An S3 class \code{rmpowerequiv} object
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{rmpowerequiv}}
//'
//' @examples
//'
//' rmsamplesizeequiv(beta = 0.1, kMax = 2, informationRates = c(0.5, 1),
//'                   alpha = 0.05, typeAlphaSpending = "sfOF",
//'                   milestone = 18,
//'                   rmstDiffLower = -2, rmstDiffUpper = 2,
//'                   allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'                   accrualIntensity = 26/9*seq(1, 9),
//'                   piecewiseSurvivalTime = c(0, 6),
//'                   stratumFraction = c(0.2, 0.8),
//'                   lambda1 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'                   lambda2 = c(0.0533, 0.0533, 1.5*0.0533, 1.5*0.0533),
//'                   gamma1 = -log(1-0.05)/12,
//'                   gamma2 = -log(1-0.05)/12, accrualDuration = NA,
//'                   followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List rmsamplesizeequiv(
    const double beta = 0.2,
    const int kMax = 1,
    const Rcpp::NumericVector& informationRates = NA_REAL,
    const Rcpp::NumericVector& criticalValues = NA_REAL,
    const double alpha = 0.05,
    const std::string typeAlphaSpending = "sfOF",
    const double parameterAlphaSpending = NA_REAL,
    const Rcpp::NumericVector& userAlphaSpending = NA_REAL,
    const double milestone = NA_REAL,
    const double rmstDiffLower = NA_REAL,
    const double rmstDiffUpper = NA_REAL,
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
    const bool fixedFollowup = 0,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const bool rounding = 1) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  ListCpp out = rmsamplesizeequivcpp(
    beta, kMax, infoRates,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    milestone, rmstDiffLower, rmstDiffUpper, allocationRatioPlanned,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, rounding);

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "rmpowerequiv";
  return result;
}
