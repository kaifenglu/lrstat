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


double kmsurv1cpp(
    const double time,
    const double allocationRatioPlanned,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2) {

  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // lamgam = lambda + gamma for each group and interval
  std::size_t J = piecewiseSurvivalTime.size();
  std::vector<double> lamgam1(J), lamgam2(J);
  for (std::size_t j = 0; j < J; ++j) {
    lamgam1[j] = lambda1[j] + gamma1[j];
    lamgam2[j] = lambda2[j] + gamma2[j];
  }

  std::size_t m = findInterval1(time, piecewiseSurvivalTime);

  // compute v then return exp(-v)
  double v = 0.0, ch1 = 0.0, ch2 = 0.0;
  // iterate intervals j = 0..m-1
  for (std::size_t j = 0; j < m; ++j) {
    if (j > 0) {
      double dt = piecewiseSurvivalTime[j] - piecewiseSurvivalTime[j - 1];
      ch1 += lamgam1[j - 1] * dt;
      ch2 += lamgam2[j - 1] * dt;
    }

    double b1 = phi * std::exp(-ch1);
    double b2 = (1.0 - phi) * std::exp(-ch2);

    // u = length of interval considered
    double u;
    if (j < m - 1) {
      u = piecewiseSurvivalTime[j + 1] - piecewiseSurvivalTime[j];
    } else {
      u = time - piecewiseSurvivalTime[j];
    }

    double d = lamgam1[j] - lamgam2[j]; // can be zero
    double v1 = lambda2[j] * u;
    double v2 = lambda1[j] - lambda2[j];

    double v3;
    if (d == 0.0) {
      // limit when d->0: log((b2 + b1*exp(-d*u)) / (b1 + b2)) / d
      v3 = - (b1 * u) / (b1 + b2);
    } else {
      double numer = b2 + b1 * std::exp(-d * u);
      double denom = b1 + b2;
      v3 = std::log(numer / denom) / d;
    }

    v += v1 - v2 * v3;
  }

  // convert to survival: exp(-v)
  return std::exp(-v);
}


std::vector<double> kmsurvcpp(
    const std::vector<double>& time,
    const double allocationRatioPlanned,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2) {

  std::size_t J = static_cast<int>(piecewiseSurvivalTime.size());
  auto lambda1x = expand_stratified(lambda1, 1, J, "lambda1");
  auto lambda2x = expand_stratified(lambda2, 1, J, "lambda2");
  auto gamma1x = expand_stratified(gamma1, 1, J, "gamma1");
  auto gamma2x = expand_stratified(gamma2, 1, J, "gamma2");

  std::size_t k = time.size();
  std::vector<double> v(k);
  for (std::size_t i = 0; i < k; ++i) {
    v[i] = kmsurv1cpp(time[i], allocationRatioPlanned, piecewiseSurvivalTime,
                      lambda1x, lambda2x, gamma1x, gamma2x);
  }
  return v;
}


//' @title Kaplan-Meier Survival Probability Based on Pooled Sample
//' @description Obtains the limit of Kaplan-Meier estimate of the survival
//' probabilities based on the pooled sample.
//'
//' @param time A vector of analysis times at which to calculate the
//'   Kaplan-Meier Survival Probability.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_lambda1
//' @inheritParams param_lambda2
//' @inheritParams param_gamma1
//' @inheritParams param_gamma2
//'
//' @return A vector of Kaplan-Meier survival probabilities at the
//' specified analysis times for piecewise exponential survival and
//' dropout distributions.
//'
//' @keywords internal
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise exponential survivals, and 5% dropout by the end of
//' # 1 year.
//'
//' kmsurv(t = c(2, 8), allocationRatioPlanned = 1,
//'        piecewiseSurvivalTime = c(0, 6),
//'        lambda1 = c(0.0533, 0.0309), lambda2 = c(0.0533, 0.0533),
//'        gamma1 = -log(1-0.05)/12, gamma2 = -log(1-0.05)/12)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector kmsurv(
    const Rcpp::NumericVector& time = NA_REAL,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0) {

  auto time1 = Rcpp::as<std::vector<double>>(time);
  auto piecewiseTime = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  auto result = kmsurvcpp(
    time1, allocationRatioPlanned, piecewiseTime, lam1, lam2, gam1, gam2);

  return Rcpp::wrap(result);
}



std::vector<double> make_breaks(
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& accrualTime,
    double accrualDuration,
    double minFollowupTime,
    double maxFollowupTime) {

  if (!(maxFollowupTime > 0.0)) return std::vector<double>{0.0, maxFollowupTime};

  // Candidate set
  std::vector<double> pts;
  pts.reserve(piecewiseSurvivalTime.size() + accrualTime.size());

  // add piecewise survival cut points (exclude 0 and maxFollowupTime)
  for (double t : piecewiseSurvivalTime) {
    if (t <= 0.0 || t >= maxFollowupTime) continue;
    pts.push_back(t);
  }

  // add accrual-derived points: accrualDuration + minFollowupTime - accrualTime[i]
  for (double at : accrualTime) {
    double s = accrualDuration + minFollowupTime - at;
    if (s <= 0.0 || s >= maxFollowupTime) continue;
    pts.push_back(s);
  }

  // Also include the two special switch points where u(t) hits its clamp bounds:
  // s(t)=accrualDuration  => t = minFollowupTime
  // s(t)=0                => t = accrualDuration + minFollowupTime
  // They may not fall on the accrualTime grid; include them if inside (0,max)
  double t_1 = minFollowupTime;
  double t_2 = accrualDuration + minFollowupTime;
  if (t_1 > 0.0 && t_1 < maxFollowupTime) pts.push_back(t_1);
  if (t_2 > 0.0 && t_2 < maxFollowupTime) pts.push_back(t_2);

  if (pts.empty()) {
    // no internal break points -> simple [0, max]
    return std::vector<double>{0.0, maxFollowupTime};
  }

  // sort and deduplicate with tolerance
  std::sort(pts.begin(), pts.end());
  double eps = std::max(1e-12, 1e-15 * std::max(1.0, std::fabs(maxFollowupTime)));

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
  breaks.push_back(maxFollowupTime);

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


DataFrameCpp lrstat1cpp(
    double time,
    double hazardRatioH0,
    double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    double accrualDuration,
    double followupTime,
    bool fixedFollowup,
    double rho1,
    double rho2,
    bool predictEventOnly) {

  std::size_t nstrata = stratumFraction.size();
  std::size_t nintervals = piecewiseSurvivalTime.size();

  // determine maxFollowupTime like original
  double maxFollowupTime = fixedFollowup ? followupTime :
    (accrualDuration + followupTime);

  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // total enrolled by calendar time
  double a = accrual1(time, accrualTime, accrualIntensity, accrualDuration);

  // precompute a2 for nfmax
  double a2 = accrual1(time - maxFollowupTime, accrualTime, accrualIntensity,
                       accrualDuration);

  std::vector<int> stratum(nstrata);
  std::vector<double> times(nstrata), nsubjects(nstrata);
  std::vector<double> nevents1(nstrata), nevents2(nstrata), nevents(nstrata);
  std::vector<double> ndropouts1(nstrata), ndropouts2(nstrata), ndropouts(nstrata);
  std::vector<double> nfmax1(nstrata), nfmax2(nstrata), nfmax(nstrata);
  std::vector<double> uscore(nstrata), vscore(nstrata), iscore(nstrata);

  DataFrameCpp result;
  double tol = 1e-6;

  // loop strata
  for (std::size_t h = 0; h < nstrata; ++h) {
    stratum[h] = h + 1;
    times[h] = time;

    double frac = stratumFraction[h];

    // build per-stratum vectors for the intervals
    int start = h * nintervals;
    int end = start + nintervals;
    std::vector<double> lam1 = subset(lambda1, start, end);
    std::vector<double> lam2 = subset(lambda2, start, end);
    std::vector<double> gam1 = subset(gamma1, start, end);
    std::vector<double> gam2 = subset(gamma2, start, end);

    // number of events in the stratum at calendar time
    std::vector<double> accrualIntensity_frac = accrualIntensity;
    for (double &v : accrualIntensity_frac) v *= frac;

    auto ne_row = nevent21cpp(
      time, allocationRatioPlanned, accrualTime, accrualIntensity_frac,
      piecewiseSurvivalTime, lam1, lam2, gam1, gam2, accrualDuration,
      followupTime, maxFollowupTime);

    // obtain dropouts by swapping hazard roles (as in original)
    auto nd_row = nevent21cpp(
      time, allocationRatioPlanned, accrualTime, accrualIntensity_frac,
      piecewiseSurvivalTime, gam1, gam2, lam1, lam2, accrualDuration,
      followupTime, maxFollowupTime);

    // number of subjects, events, dropouts in the stratum at calendar time
    nsubjects[h] = frac * a;
    nevents1[h] = ne_row.first;
    nevents2[h] = ne_row.second;
    nevents[h] = nevents1[h] + nevents2[h];
    ndropouts1[h] = nd_row.first;
    ndropouts2[h] = nd_row.second;
    ndropouts[h] = ndropouts1[h] + ndropouts2[h];

    // nfmax: those reaching max follow-up (censoring at max follow-up)
    double p1 = patrisk1(maxFollowupTime, piecewiseSurvivalTime, lam1, gam1);
    double p2 = patrisk1(maxFollowupTime, piecewiseSurvivalTime, lam2, gam2);
    double ncom = frac * a2;
    nfmax1[h] = phi * ncom * p1;
    nfmax2[h] = (1.0 - phi) * ncom * p2;
    nfmax[h]  = nfmax1[h] + nfmax2[h];

    if (!predictEventOnly) {

      double accrualDuration0 = std::min(time, accrualDuration);
      double minFollowupTime0 = std::max(time - accrualDuration, 0.0);
      double maxFollowupTime0 = std::min(time, maxFollowupTime);

      auto fu = [&](double t)->double {
        // call natriskcpp for a single time point
        auto risk = natrisk1cpp(
          t, allocationRatioPlanned, accrualTime, accrualIntensity_frac,
          piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
          accrualDuration0, minFollowupTime0, maxFollowupTime0);

        // extract r1,r2 for row 0 (time t) and columns 0,1 (groups 1,2)
        double r1 = risk.first;
        double r2 = risk.second;

        // find interval j for time t
        int j = findInterval1(t, piecewiseSurvivalTime) - 1;

        // weight w
        double w = 1.0;
        if (rho1 != 0.0 || rho2 != 0.0) {
          double s = kmsurv1cpp(
            t, allocationRatioPlanned, piecewiseSurvivalTime,
            lam1, lam2, gam1, gam2);
          w = std::pow(s, rho1) * std::pow(1.0 - s, rho2);
        }

        double denom = r1 * hazardRatioH0 + r2;
        if (denom <= 0.0) return 0.0;

        double N = (r1 * hazardRatioH0) * r2 / denom;
        double d = lam1[j] / hazardRatioH0 - lam2[j];

        double val = w * N * d;
        return val;
      };

      auto fv = [&](double t)->double {
        auto risk = natrisk1cpp(
          t, allocationRatioPlanned, accrualTime, accrualIntensity_frac,
          piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
          accrualDuration0, minFollowupTime0, maxFollowupTime0);

        double r1 = risk.first;
        double r2 = risk.second;

        int j = findInterval1(t, piecewiseSurvivalTime) - 1;

        double w = 1.0;
        if (rho1 != 0.0 || rho2 != 0.0) {
          double s = kmsurv1cpp(
            t, allocationRatioPlanned, piecewiseSurvivalTime,
            lam1, lam2, gam1, gam2);
          w = std::pow(s, rho1) * std::pow(1.0 - s, rho2);
        }

        double denom = r1 * hazardRatioH0 + r2;
        if (denom <= 0.0) return 0.0;

        double N = (r1 * hazardRatioH0) * r2 / (denom * denom);
        double d = r1 * lam1[j] + r2 * lam2[j];

        double val = w * w * N * d;
        return val;
      };

      auto fi = [&](double t)->double {
        auto risk = natrisk1cpp(
          t, allocationRatioPlanned, accrualTime, accrualIntensity_frac,
          piecewiseSurvivalTime, lam1, lam2, gam1, gam2,
          accrualDuration0, minFollowupTime0, maxFollowupTime0);

        double r1 = risk.first;
        double r2 = risk.second;

        int j = findInterval1(t, piecewiseSurvivalTime) - 1;

        double w = 1.0;
        if (rho1 != 0.0 || rho2 != 0.0) {
          double s = kmsurv1cpp(
            t, allocationRatioPlanned, piecewiseSurvivalTime,
            lam1, lam2, gam1, gam2);
          w = std::pow(s, rho1) * std::pow(1.0 - s, rho2);
        }

        double denom = r1 * hazardRatioH0 + r2;
        if (denom <= 0.0) return 0.0;

        double N = (r1 * hazardRatioH0) * r2 / (denom * denom);
        double d = r1 * lam1[j] + r2 *lam2[j];

        double val = w * N * d;
        return val;
      };

      // integrate fu, fv, fi over [0, maxFollowupTime0] by intervals
      auto breaks = make_breaks(piecewiseSurvivalTime, accrualTime,
                                accrualDuration0, minFollowupTime0,
                                maxFollowupTime0);

      uscore[h] = integrate3(fu, breaks, tol);
      vscore[h] = integrate3(fv, breaks, tol);
      iscore[h] = integrate3(fi, breaks, tol);
    }
  }

  result.push_back(std::move(stratum), "stratum");
  result.push_back(std::move(times), "time");
  result.push_back(std::move(nsubjects), "subjects");
  result.push_back(std::move(nevents), "nevents");
  result.push_back(std::move(nevents1), "nevents1");
  result.push_back(std::move(nevents2), "nevents2");
  result.push_back(std::move(ndropouts), "ndropouts");
  result.push_back(std::move(ndropouts1), "ndropouts1");
  result.push_back(std::move(ndropouts2), "ndropouts2");
  result.push_back(std::move(nfmax), "nfmax");
  result.push_back(std::move(nfmax1), "nfmax1");
  result.push_back(std::move(nfmax2), "nfmax2");

  if (!predictEventOnly) {
    result.push_back(std::move(uscore), "uscore");
    result.push_back(std::move(vscore), "vscore");
    result.push_back(std::move(iscore), "iscore");
  }

  return result;
}


double extract_lr(const DataFrameCpp& df, const char* name) {
  auto vec = df.get<double>(name);
  return std::accumulate(vec.begin(), vec.end(), 0.0);
};


DataFrameCpp lrstatcpp(
    const std::vector<double>& time,
    double hazardRatioH0,
    double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    double accrualDuration,
    double followupTime,
    bool fixedFollowup,
    double rho1,
    double rho2,
    int predictTarget) {

  if (!none_na(time)) throw std::invalid_argument("time must be provided");
  for (double v : time) {
    if (v < 0.0) throw std::invalid_argument("time must be non-negative");
  }

  if (hazardRatioH0 <= 0.0)
    throw std::invalid_argument("hazardRatioH0 must be positive");

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
  if (!none_na(gamma1)) throw std::invalid_argument("gamma1 must be provided");
  if (!none_na(gamma2)) throw std::invalid_argument("gamma2 must be provided");

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

  if (rho1 < 0.0) throw std::invalid_argument("rho1 must be non-negative");
  if (rho2 < 0.0) throw std::invalid_argument("rho2 must be non-negative");

  if (predictTarget != 1 && predictTarget != 2 && predictTarget != 3)
    throw std::invalid_argument("predictTarget must be equal to 1, 2, or 3");
  const bool predictEventOnly = (predictTarget == 1);

  // expand stratified rates
  std::size_t nstrata = stratumFraction.size();
  std::size_t nintervals = piecewiseSurvivalTime.size();
  auto lambda1x = expand_stratified(lambda1, nstrata, nintervals, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintervals, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintervals, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintervals, "gamma2");

  // prepare output
  std::size_t k = time.size();
  std::vector<double> subjects(k), nevents(k), nevents1(k), nevents2(k);
  std::vector<double> ndropouts(k), ndropouts1(k), ndropouts2(k);
  std::vector<double> nfmax(k), nfmax1(k), nfmax2(k);
  std::vector<double> uscore(k), vscore(k), logRankZ(k);
  std::vector<double> HR(k), vlogHR(k), zlogHR(k);

  // For each requested time call lrstat1cpp and aggregate across strata
  for (std::size_t i = 0; i < k; ++i) {
    DataFrameCpp df = lrstat1cpp(
      time[i], hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup,
      rho1, rho2, predictEventOnly
    );

    // sums over strata
    subjects[i] = extract_lr(df, "subjects");
    nevents[i] = extract_lr(df, "nevents");
    nevents1[i] = extract_lr(df, "nevents1");
    nevents2[i] = extract_lr(df, "nevents2");
    ndropouts[i] = extract_lr(df, "ndropouts");
    ndropouts1[i] = extract_lr(df, "ndropouts1");
    ndropouts2[i] = extract_lr(df, "ndropouts2");
    nfmax[i] = extract_lr(df, "nfmax");
    nfmax1[i] = extract_lr(df, "nfmax1");
    nfmax2[i] = extract_lr(df, "nfmax2");

    if (!predictEventOnly) {
      double us = extract_lr(df, "uscore");
      double vs = extract_lr(df, "vscore");
      uscore[i] = us;
      vscore[i] = vs;
      logRankZ[i] = vs > 0.0 ? us / std::sqrt(vs) : 0.0;
    }
  }

  // If predictTarget == 3, solve for weighted Cox estimator (logHR) for each time
  if (predictTarget == 3) {
    double logHR0 = std::log(hazardRatioH0);

    // g(beta) = sum_uscore at hazardRatio = exp(beta)
    for (std::size_t j = 0; j < k; ++j) {
      auto g = [&](double beta)->double {
        double hazardRatio = std::exp(beta);
        DataFrameCpp rr = lrstat1cpp(
          time[j], hazardRatio, allocationRatioPlanned,
          accrualTime, accrualIntensity,
          piecewiseSurvivalTime, stratumFraction,
          lambda1x, lambda2x, gamma1x, gamma2x,
          accrualDuration, followupTime, fixedFollowup,
          rho1, rho2, false
        );
        return extract_lr(rr, "uscore");
      };

      double logHR = brent(g, -4.6, 4.6, 1.0e-6);
      HR[j] = std::exp(logHR);

      DataFrameCpp rr2 = lrstat1cpp(
        time[j], HR[j], allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrualDuration, followupTime, fixedFollowup,
        rho1, rho2, false
      );

      double vs = extract_lr(rr2, "vscore");
      double is = extract_lr(rr2, "iscore");
      double vlog = 0.0, zlog = 0.0;
      if (is != 0.0) {
        vlog = vs / (is * is);
        if (vlog > 0.0) zlog = (logHR - logHR0) / std::sqrt(vlog);
      }
      vlogHR[j] = vlog;
      zlogHR[j] = zlog;
    }
  }

  DataFrameCpp out;
  out.push_back(time, "time");
  out.push_back(std::move(subjects), "subjects");
  out.push_back(std::move(nevents), "nevents");
  out.push_back(std::move(nevents1), "nevents1");
  out.push_back(std::move(nevents2), "nevents2");
  out.push_back(std::move(ndropouts), "ndropouts");
  out.push_back(std::move(ndropouts1), "ndropouts1");
  out.push_back(std::move(ndropouts2), "ndropouts2");
  out.push_back(std::move(nfmax), "nfmax");
  out.push_back(std::move(nfmax1), "nfmax1");
  out.push_back(std::move(nfmax2), "nfmax2");

  if (!predictEventOnly) {
    out.push_back(std::move(uscore), "uscore");
    out.push_back(std::move(vscore), "vscore");
    out.push_back(std::move(logRankZ), "logRankZ");
    out.push_back(hazardRatioH0, "hazardRatioH0");
  }

  if (predictTarget == 3) {
    out.push_back(std::move(HR), "HR");
    out.push_back(std::move(vlogHR), "vlogHR");
    out.push_back(std::move(zlogHR), "zlogHR");
  }

  return out;
}


//' @title Number of Subjects Having an Event and Log-Rank Statistics
//' @description Obtains the number of subjects accrued, number of events,
//' number of dropouts, and number of subjects reaching the maximum
//' follow-up in each group, mean and variance of weighted log-rank
//' score statistic, estimated hazard ratio from weighted Cox regression
//' and variance of log hazard ratio estimate at given calendar times.
//'
//' @param time A vector of calendar times at which to calculate the number
//'   of events and the mean and variance of log-rank test score statistic.
//' @inheritParams param_hazardRatioH0
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
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @param predictTarget The target of prediction.
//'   Set \code{predictTarget = 1} to predict the number of events only.
//'   Set \code{predictTarget = 2} (default) to predict the number of events
//'   and log-rank score statistic mean and variance.
//'   Set \code{predictTarget = 3} to predict the number of events,
//'   log-rank score statistic mean and variance, and
//'   hazard ratio and variance of log hazard ratio.
//'
//' @return A data frame containing the following variables if
//' \code{predictTarget = 1}:
//'
//' * \code{time}: The analysis time since trial start.
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
//' * \code{nfmax}: The total number of subjects reaching maximum follow-up.
//'
//' * \code{nfmax1}: The number of subjects reaching maximum follow-up in
//'   the active treatment group.
//'
//' * \code{nfmax2}: The number of subjects reaching maximum follow-up in
//'   the control group.
//'
//' If \code{predictTarget = 2}, the following variables will also
//' be included:
//'
//' * \code{uscore}: The numerator of the log-rank test statistic.
//'
//' * \code{vscore}: The variance of the log-rank score test statistic.
//'
//' * \code{logRankZ}: The log-rank test statistic on the Z-scale.
//'
//' * \code{hazardRatioH0}: The hazard ratio under the null hypothesis.
//'
//' Furthermore, if \code{predictTarget = 3}, the following additional
//' variables will also be included:
//'
//' * \code{HR}: The average hazard ratio from weighted Cox regression.
//'
//' * \code{vlogHR}: The variance of log hazard ratio.
//'
//' * \code{zlogHR}: The Z-statistic for log hazard ratio.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' lrstat(time = c(22, 40), allocationRatioPlanned = 1,
//'        accrualTime = seq(0, 8),
//'        accrualIntensity = 26/9*seq(1, 9),
//'        piecewiseSurvivalTime = c(0, 6),
//'        lambda1 = c(0.0533, 0.0309),
//'        lambda2 = c(0.0533, 0.0533),
//'        gamma1 = -log(1-0.05)/12,
//'        gamma2 = -log(1-0.05)/12,
//'        accrualDuration = 22,
//'        followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame lrstat(
    const Rcpp::NumericVector& time = NA_REAL,
    const double hazardRatioH0 = 1,
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
    const double rho1 = 0,
    const double rho2 = 0,
    const int predictTarget = 2) {

  auto time1 = Rcpp::as<std::vector<double>>(time);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto piecewiseTime = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  DataFrameCpp df = lrstatcpp(
    time1, hazardRatioH0, allocationRatioPlanned,
    accrualT, accrualInt, piecewiseTime,
    stratumFrac, lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup,
    rho1, rho2, predictTarget
  );

  return Rcpp::wrap(df);
}


// caltime1cpp: returns calendar time matching target event count (scalar)
double caltime1cpp(
    const double nevents,
    double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    double accrualDuration,
    double followupTime,
    bool fixedFollowup) {

  // given candidate calendar time t returns predicted total events - nevents
  auto f = [&](double t)->double {
    DataFrameCpp df = lrstat1cpp(
      t, 1.0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1, lambda2, gamma1, gamma2,
      accrualDuration, followupTime, fixedFollowup,
      0.0, 0.0, true
    );
    return extract_lr(df, "nevents") - nevents;
  };

  // initial feasibility check: use maximum target event
  const double studyTime = accrualDuration + followupTime;
  if (f(studyTime) < 0.0) throw std::invalid_argument(
      "followupTime is too short to reach the target number of events");

  double out = brent(f, 0.0, studyTime, 1e-6);
  return out;
}


// caltimecpp: returns vector of calendar times matching target event counts
std::vector<double> caltimecpp(
    const std::vector<double>& nevents,
    double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    double accrualDuration,
    double followupTime,
    bool fixedFollowup) {

  if (!none_na(nevents)) throw std::invalid_argument("nevents must be provided");
  for (double v : nevents) {
    if (v <= 0.0) throw std::invalid_argument("nevents must be positive");
  }

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
  if (!none_na(gamma1)) throw std::invalid_argument("gamma1 must be provided");
  if (!none_na(gamma2)) throw std::invalid_argument("gamma2 must be provided");

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
  std::size_t nstrata = stratumFraction.size();
  std::size_t nintervals = piecewiseSurvivalTime.size();
  auto lambda1x = expand_stratified(lambda1, nstrata, nintervals, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintervals, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintervals, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintervals, "gamma2");

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

  // Prepare output
  std::size_t k = nevents.size();
  std::vector<double> out(k);
  for (std::size_t i = 0; i < k; ++i) {
    double eventTarget = std::max(nevents[i], 0.0);

    out[i] = caltime1cpp(
      eventTarget, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup);
  }

  return out;
}


//' @title Calendar Times for Target Number of Events
//' @description Obtains the calendar times needed to reach the target
//' number of subjects experiencing an event.
//'
//' @param nevents A vector of target number of events.
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
//' @return A vector of calendar times expected to yield the target
//' number of events.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' caltime(nevents = c(24, 80), allocationRatioPlanned = 1,
//'         accrualTime = seq(0, 8),
//'         accrualIntensity = 26/9*seq(1, 9),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309),
//'         lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12,
//'         gamma2 = -log(1-0.05)/12,
//'         accrualDuration = 22,
//'         followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector caltime(
    const Rcpp::NumericVector& nevents = NA_REAL,
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

  auto nevents1 = Rcpp::as<std::vector<double>>(nevents);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto piecewiseTime = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  auto out = caltimecpp(
    nevents1, allocationRatioPlanned, accrualT, accrualInt,
    piecewiseTime, stratumFrac, lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup
  );

  return Rcpp::wrap(out);
}


DataFrameCpp getDurationFromNeventscpp(
    double nevents,
    double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    double followupTime,
    bool fixedFollowup,
    const int npoints = 23) {

  // Input validation (mirror R checks)
  if (std::isnan(nevents)) throw std::invalid_argument("nevents must be provided");
  if (nevents <= 0.0) throw std::invalid_argument("nevents must be positive");

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
  if (!none_na(gamma1)) throw std::invalid_argument("gamma1 must be provided");
  if (!none_na(gamma2)) throw std::invalid_argument("gamma2 must be provided");

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
  std::size_t nstrata = stratumFraction.size();
  std::size_t nintervals = piecewiseSurvivalTime.size();
  auto lambda1x = expand_stratified(lambda1, nstrata, nintervals, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintervals, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintervals, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintervals, "gamma2");

  if (fixedFollowup && std::isnan(followupTime))
    throw std::invalid_argument("followupTime must be provided for fixed follow-up");
  if (fixedFollowup && followupTime <= 0.0)
    throw std::invalid_argument("followupTime must be positive for fixed follow-up");

  if (npoints < 2)
    throw std::invalid_argument("npoints must be greater than or equal to 2");

  // total predicted events at analysisTime, given accrualDuration and followupTime
  auto total_events_at = [&](double t, double accrdur, double futime)->double {
    DataFrameCpp df = lrstat1cpp(
      t, 1.0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrdur, futime, fixedFollowup,
      0.0, 0.0, true
    );
    return extract_lr(df, "nevents");
  };

  // 1) find minimal accrualDuration t[0] such that with Tf1 we reach nevents
  double Tf1 = fixedFollowup ? followupTime : 1000.0;
  // f1(accrD) = total_events_at(accrD + Tf1, accrD, Tf1) - nevents
  auto f1 = [&](double accrdur)->double {
    return total_events_at(accrdur + Tf1, accrdur, Tf1) - nevents;
  };

  // bracket and find lower bound root for f1
  double lower = 0.001;
  double upper = 240.0;
  // ensure bracket where f1(upper) >= 0
  int iter = 0;
  const int max_iter_doubling = 60;
  while (f1(upper) < 0.0) {
    lower = upper;
    upper *= 2.0;
    if (++iter > max_iter_doubling)
      throw std::runtime_error("Unable to bracket minimal accrualDuration");
  }
  double t0 = brent(f1, lower, upper, 1e-6);

  // 2) find maximal accrualDuration t[1] such that with Tf2 we reach nevents
  double Tf2 = fixedFollowup ? followupTime : 0.0;
  auto f2 = [&](double accrdur)->double {
    return total_events_at(accrdur, accrdur, Tf2) - nevents;
  };

  // bracket between t0 and previous upper (upper already >= necessary)
  double t1 = brent(f2, t0, upper, 1e-6);

  // 3) Prepare grid of accrual durations ta (npoints) between t0 and t1
  std::size_t npts = static_cast<std::size_t>(npoints);
  std::vector<double> ta(npts), ts(npts), tf(npts);
  double dt = (t1 - t0) / static_cast<double>(npts - 1);
  for (std::size_t i = 0; i < npts; ++i) ta[i] = t0 + i * dt;

  for (std::size_t i = 0; i < npts; ++i) {
    if (i == 0) {
      ts[i] = ta[i] + Tf1;
    } else if (i == npts - 1) {
      ts[i] = ta[i];
    } else {
      ts[i] = caltime1cpp(nevents, allocationRatioPlanned,
                          accrualTime, accrualIntensity,
                          piecewiseSurvivalTime, stratumFraction,
                          lambda1x, lambda2x, gamma1x, gamma2x,
                          ta[i], Tf1, fixedFollowup);
    }
    tf[i] = fixedFollowup ? followupTime : (ts[i] - ta[i]);
  }

  // use large accrualDurationLarge to avoid truncation
  auto subjects = accrual(ta, accrualTime, accrualIntensity, 1000.0);

  // Build DataFrameCpp result
  DataFrameCpp out;

  std::vector<double> nevents_vec(npts, nevents);
  out.push_back(nevents_vec, "nevents");
  out.push_back(fixedFollowup, "fixedFollowup");
  out.push_back(ta, "accrualDuration");
  out.push_back(subjects, "subjects");
  out.push_back(tf, "followupTime");
  out.push_back(ts, "studyDuration");

  return out;
}


//' @title Range of Accrual Duration for Target Number of Events
//' @description Obtains a range of accrual duration to reach the
//' target number of events.
//'
//' @param nevents The target number of events.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @inheritParams param_lambda1_stratified
//' @inheritParams param_lambda2_stratified
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @param followupTime Follow-up time for the last enrolled subjects.
//'   Must be provided for fixed follow-up design.
//' @inheritParams param_fixedFollowup
//' @param npoints The number of accrual duration time points.
//'   Defaults to 23.
//'
//' @return A data frame of the following variables:
//'
//' * \code{nevents}: The target number of events.
//'
//' * \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//' * \code{accrualDuration}: The accrual duration.
//'
//' * \code{subjects}: The total number of subjects.
//'
//' * \code{followupTime}: The follow-up time for the last enrolled subject.
//'
//' * \code{studyDuration}: The study duration.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survivals, and 5% dropout by
//' # the end of 1 year.
//'
//' getDurationFromNevents(
//'   nevents = 80, allocationRatioPlanned = 1,
//'   accrualTime = seq(0, 8),
//'   accrualIntensity = 26/9*seq(1, 9),
//'   piecewiseSurvivalTime = c(0, 6),
//'   lambda1 = c(0.0533, 0.0309),
//'   lambda2 = c(0.0533, 0.0533),
//'   gamma1 = -log(1-0.05)/12,
//'   gamma2 = -log(1-0.05)/12,
//'   fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame getDurationFromNevents(
    const double nevents = NA_REAL,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = false,
    const int npoints = 23) {

  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto piecewiseTime = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  auto df = getDurationFromNeventscpp(
    nevents, allocationRatioPlanned, accrualT, accrualInt,
    piecewiseTime, stratumFrac, lam1, lam2, gam1, gam2,
    followupTime, fixedFollowup, npoints
  );

  return Rcpp::wrap(df);
}


ListCpp lrpowercpp(
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
    const double hazardRatioH0,
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
    const double rho1,
    const double rho2,
    const bool estimateHazardRatio,
    const std::string& typeOfComputation,
    const std::vector<double>& spendingTime,
    const double studyDuration) {

  // --- Input validation & normalization -------------------
  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }

  std::size_t K = static_cast<std::size_t>(kMax);

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
    for (std::size_t i = 0; i < K; ++i)
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
    for (std::size_t i = 0; i < K - 1; ++i) {
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

  if (hazardRatioH0 <= 0.0)
    throw std::invalid_argument("hazardRatioH0 must be positive");

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
  if (!none_na(gamma1)) throw std::invalid_argument("gamma1 must be provided");
  if (!none_na(gamma2)) throw std::invalid_argument("gamma2 must be provided");

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
  std::size_t nstrata = stratumFraction.size();
  std::size_t nintervals = piecewiseSurvivalTime.size();
  std::size_t nsi = nstrata * nintervals;
  auto lambda1x = expand_stratified(lambda1, nstrata, nintervals, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintervals, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintervals, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintervals, "gamma2");

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

  if (rho1 < 0.0) throw std::invalid_argument("rho1 must be non-negative");
  if (rho2 < 0.0) throw std::invalid_argument("rho2 must be non-negative");

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

  if (fixedFollowup && !std::isnan(studyDuration) &&
      (studyDuration < accrualDuration))
    throw std::invalid_argument("studyDuration must be >= accrualDuration");
  if (fixedFollowup && !std::isnan(studyDuration) &&
      (studyDuration > accrualDuration + followupTime))
    throw std::invalid_argument(
        "studyDuration must be <= accrualDuration + followupTime");


  // --- Determine if Schoenfeld method is eligible-----------------------------
  double hazardRatio = 1.0;
  std::vector<double> hrx(nsi);
  for (std::size_t i = 0; i < nsi; ++i) hrx[i] = lambda1x[i] / lambda2x[i];
  bool proportional = true;
  double hrx0 = hrx[0];
  for (std::size_t i = 1; i < nsi; ++i) {
    if (std::fabs(hrx[i] - hrx0) > 1e-8) { proportional = false; break; }
  }
  if (proportional) hazardRatio = hrx0;

  bool schoenfeld_eligible = proportional && rho1 == 0.0 && rho2 == 0.0;
  std::string su = typeOfComputation;
  if (su.empty()) {
    su = schoenfeld_eligible ? "schoenfeld" : "direct";
  } else {
    for (char &c : su) {
      c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
  }
  char su1 = su.front();
  if (su1 != 'd' && su1 != 's')
    throw std::invalid_argument("typeOfComputation must be direct or schoenfeld");
  if (!schoenfeld_eligible && su1 == 's') {
    throw std::invalid_argument(
        "Schoenfeld method can only be used for ordinary log-rank test "
        "with proportional hazards");
  }

  // --- Efficacy boundaries ---------------------------------------------------
  std::vector<double> l(K, -6.0), zero(K, 0.0);
  std::vector<double> critValues = criticalValues;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (K > 1 && criticalValues.size() == K) {
      bool hasNaN = false;
      for (std::size_t i = 0; i < K - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[K-1])) {
        haybittle = true;
      }
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(K);
      for (std::size_t i = 0; i < K - 1; ++i) {
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
  double alpha1 = missingCriticalValues ? alpha : cumAlphaSpent[K-1];

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

  // --- Analysis timing, number of events, and information --------------------
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  double studyDuration1 = studyDuration;
  if (!fixedFollowup || std::isnan(studyDuration))
    studyDuration1 = accrualDuration + followupTime;

  std::vector<double> time(K);
  std::vector<double> nsubjects(K), nsubjects1(K), nsubjects2(K);
  std::vector<double> nevents(K), nevents1(K), nevents2(K);
  std::vector<double> ndropouts(K), ndropouts1(K), ndropouts2(K);
  std::vector<double> theta(K), info(K);

  if (rho1 == 0 && rho2 == 0) {
    // ordinary log-rank test: can directly use information = phi*(1-phi)*events
    double vtrt = phi * (1.0 - phi);
    double theta1 = -std::log(hazardRatio / hazardRatioH0);

    bool predictEventOnly = (su1 == 's');

    DataFrameCpp lr = lrstat1cpp(
      studyDuration1, hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup,
      rho1, rho2, predictEventOnly);

    double totalEvents = extract_lr(lr, "nevents");

    time[K - 1] = studyDuration1;
    nsubjects[K - 1] = extract_lr(lr, "subjects");
    nsubjects1[K - 1] = phi * nsubjects[K - 1];
    nsubjects2[K - 1] = (1.0 - phi) * nsubjects[K - 1];
    nevents[K - 1] = totalEvents;
    nevents1[K - 1] = extract_lr(lr, "nevents1");
    nevents2[K - 1] = extract_lr(lr, "nevents2");
    ndropouts[K - 1] = extract_lr(lr, "ndropouts");
    ndropouts1[K - 1] = extract_lr(lr, "ndropouts1");
    ndropouts2[K - 1] = extract_lr(lr, "ndropouts2");

    if (su1 == 's') {
      theta[K - 1] = theta1;
      info[K - 1] = vtrt * totalEvents;
    } else {
      double uscore = extract_lr(lr, "uscore");
      double vscore = extract_lr(lr, "vscore");
      theta[K - 1] = -uscore / vscore;
      info[K - 1] = vscore;
    }

    // for interim analyses
    for (std::size_t i = 0; i < K - 1; ++i) {
      double nevents_target = totalEvents * infoRates[i];

      time[i] = caltime1cpp(
        nevents_target, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrualDuration, followupTime, fixedFollowup);

      DataFrameCpp lr1 = lrstat1cpp(
        time[i], hazardRatioH0, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrualDuration, followupTime, fixedFollowup,
        rho1, rho2, predictEventOnly);

      nsubjects[i] = extract_lr(lr1, "subjects");
      nsubjects1[i] = phi * nsubjects[i];
      nsubjects2[i] = (1.0 - phi) * nsubjects[i];
      nevents[i] = nevents_target;
      nevents1[i] = extract_lr(lr1, "nevents1");
      nevents2[i] = extract_lr(lr1, "nevents2");
      ndropouts[i] = extract_lr(lr1, "ndropouts");
      ndropouts1[i] = extract_lr(lr1, "ndropouts1");
      ndropouts2[i] = extract_lr(lr1, "ndropouts2");

      if (su1 == 's') {
        theta[i] = theta1;
        info[i] = vtrt * nevents_target;
      } else {
        double uscore1 = extract_lr(lr1, "uscore");
        double vscore1 = extract_lr(lr1, "vscore");
        theta[i] = -uscore1 / vscore1;
        info[i] = vscore1;
      }
    }
  } else {
    // general case: need maxInformation at study duration for matching
    DataFrameCpp lr = lrstat1cpp(
      studyDuration1, hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup,
      rho1, rho2, false);

    double maxInformation = extract_lr(lr, "vscore");

    time[K - 1] = studyDuration1;
    nsubjects[K - 1] = extract_lr(lr, "subjects");
    nsubjects1[K - 1] = phi * nsubjects[K - 1];
    nsubjects2[K - 1] = (1.0 - phi) * nsubjects[K - 1];
    nevents[K - 1] = extract_lr(lr, "nevents");
    nevents1[K - 1] = extract_lr(lr, "nevents1");
    nevents2[K - 1] = extract_lr(lr, "nevents2");
    ndropouts[K - 1] = extract_lr(lr, "ndropouts");
    ndropouts1[K - 1] = extract_lr(lr, "ndropouts1");
    ndropouts2[K - 1] = extract_lr(lr, "ndropouts2");

    double uscore = extract_lr(lr, "uscore");
    theta[K - 1] = -uscore / maxInformation;
    info[K - 1] =  maxInformation;

    // compute times by matching informationRates1 * maxInformation
    for (std::size_t i = 0; i < K - 1; ++i) {
      double information1 = maxInformation * infoRates[i];

      // solve for analysis time where total information equals information1
      auto g = [&](double t)->double {
        DataFrameCpp lr1 = lrstat1cpp(
          t, hazardRatioH0, allocationRatioPlanned,
          accrualTime, accrualIntensity,
          piecewiseSurvivalTime, stratumFraction,
          lambda1x, lambda2x, gamma1x, gamma2x,
          accrualDuration, followupTime, fixedFollowup,
          rho1, rho2, false);

        return extract_lr(lr1, "vscore") - information1;
      };

      time[i] = brent(g, 1e-6, studyDuration1, 1e-6);

      DataFrameCpp lr1 = lrstat1cpp(
        time[i], hazardRatioH0, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrualDuration, followupTime, fixedFollowup,
        rho1, rho2, false);

      nsubjects[i] = extract_lr(lr1, "subjects");
      nsubjects1[i] = phi * nsubjects[i];
      nsubjects2[i] = (1.0 - phi) * nsubjects[i];
      nevents[i] = extract_lr(lr1, "nevents");
      nevents1[i] = extract_lr(lr1, "nevents1");
      nevents2[i] = extract_lr(lr1, "nevents2");
      ndropouts[i] = extract_lr(lr1, "ndropouts");
      ndropouts1[i] = extract_lr(lr1, "ndropouts1");
      ndropouts2[i] = extract_lr(lr1, "ndropouts2");

      double uscore1 = extract_lr(lr1, "uscore");
      theta[i] = -uscore1 / information1;
      info[i] = information1;
    }
  }

  // --- Compute hazard-ratio estimates ----------------------------------------
  std::vector<double> HR(K), vlogHR(K);
  if (estimateHazardRatio) {
    if (su1 == 's') {
      for (std::size_t i = 0; i < K; ++i) {
        HR[i] = hazardRatio;
        vlogHR[i] = 1.0 / info[i];
      }
    } else {
      DataFrameCpp df = lrstatcpp(
        time, hazardRatioH0, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrualDuration, followupTime, fixedFollowup,
        rho1, rho2, 3);

      HR = df.get<double>("HR");
      vlogHR = df.get<double>("vlogHR");
    }
  }

  // --- compute stagewise exit probabilities and related metrics --------------
  ListCpp exit_probs;
  if (!missingFutilityBounds || bsf == "none" || K == 1) {
    exit_probs = exitprobcpp(critValues, futBounds, theta, info);
  } else {
    std::vector<double> w(K, 1.0);
    auto gp = getPower(alpha1, kMax, critValues, theta, info, bsf,
                       parameterBetaSpending, spendTime, futStopping, w);

    // gp structure: [power, futilityBounds, probs]
    // extract futility bounds and probs appropriately
    futBounds = gp.get<std::vector<double>>("futilityBounds");
    exit_probs = gp.get_list("probs");
  }

  std::vector<double> efficacyP(K), futilityP(K);
  for (std::size_t i = 0; i < K; ++i) {
    efficacyP[i] = 1 - boost_pnorm(critValues[i]);
    futilityP[i] = 1 - boost_pnorm(futBounds[i]);
  }

  auto pu = exit_probs.get<std::vector<double>>("exitProbUpper");
  auto pl = exit_probs.get<std::vector<double>>("exitProbLower");
  std::vector<double> ptotal(K);
  for (std::size_t i = 0; i < K; ++i) ptotal[i] = pu[i] + pl[i];

  double overallReject = std::accumulate(pu.begin(), pu.end(), 0.0);
  double expectedNumberOfEvents = 0.0;
  double expectedNumberOfDropouts = 0.0;
  double expectedNumberOfSubjects = 0.0;
  double expectedNumberOfEvents1 = 0.0;
  double expectedNumberOfDropouts1 = 0.0;
  double expectedNumberOfSubjects1 = 0.0;
  double expectedNumberOfEvents2 = 0.0;
  double expectedNumberOfDropouts2 = 0.0;
  double expectedNumberOfSubjects2 = 0.0;
  double expectedStudyDuration = 0.0;
  double expectedInformation = 0.0;

  for (std::size_t i = 0; i < K; ++i) {
    expectedNumberOfEvents += ptotal[i] * nevents[i];
    expectedNumberOfDropouts += ptotal[i] * ndropouts[i];
    expectedNumberOfSubjects += ptotal[i] * nsubjects[i];
    expectedNumberOfEvents1 += ptotal[i] * nevents1[i];
    expectedNumberOfDropouts1 += ptotal[i] * ndropouts1[i];
    expectedNumberOfSubjects1 += ptotal[i] * nsubjects1[i];
    expectedNumberOfEvents2 += ptotal[i] * nevents2[i];
    expectedNumberOfDropouts2 += ptotal[i] * ndropouts2[i];
    expectedNumberOfSubjects2 += ptotal[i] * nsubjects2[i];
    expectedStudyDuration += ptotal[i] * time[i];
    expectedInformation += ptotal[i] * info[i];
  }

  // cummulative sums cpu and cpl
  std::vector<double> cpu(K), cpl(K);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());

  // hazard ratio transforms for display if estimateHazardRatio
  std::vector<double> hru(K), hrl(K);
  if (estimateHazardRatio) {
    for (std::size_t i = 0; i < K; ++i) {
      hru[i] = hazardRatioH0 * std::exp(-critValues[i] * std::sqrt(vlogHR[i]));
      hrl[i] = hazardRatioH0 * std::exp(-futBounds[i] * std::sqrt(vlogHR[i]));
      if (critValues[i] == 6.0) { hru[i] = NaN; effStopping[i] = 0; }
      if (futBounds[i] == -6.0) { hrl[i] = NaN; futStopping[i] = 0; }
    }
  }

  // --- Build output DataFrames and Lists ------------------------------------
  std::string typeOfComputationDisplay = su1 == 's' ? "Schoenfeld" : "Direct";
  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha1, "alpha");
  overallResults.push_back(nevents[K-1], "numberOfEvents");
  overallResults.push_back(ndropouts[K-1], "numberOfDropouts");
  overallResults.push_back(nsubjects[K-1], "numberOfSubjects");
  overallResults.push_back(time[K-1], "studyDuration");
  overallResults.push_back(info[K-1], "information");
  overallResults.push_back(expectedNumberOfEvents, "expectedNumberOfEvents");
  overallResults.push_back(expectedNumberOfDropouts, "expectedNumberOfDropouts");
  overallResults.push_back(expectedNumberOfSubjects, "expectedNumberOfSubjects");
  overallResults.push_back(expectedStudyDuration, "expectedStudyDuration");
  overallResults.push_back(expectedInformation, "expectedInformation");
  overallResults.push_back(accrualDuration, "accrualDuration");
  overallResults.push_back(followupTime, "followupTime");
  overallResults.push_back(fixedFollowup, "fixedFollowup");
  overallResults.push_back(rho1, "rho1");
  overallResults.push_back(rho2, "rho2");
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(hazardRatioH0, "hazardRatioH0");
  overallResults.push_back(typeOfComputationDisplay, "typeOfComputation");

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
  byStageResults.push_back(std::move(time), "analysisTime");
  if (estimateHazardRatio) {
    byStageResults.push_back(std::move(hru), "efficacyHR");
    byStageResults.push_back(std::move(hrl), "futilityHR");
  }
  byStageResults.push_back(std::move(efficacyP), "efficacyP");
  byStageResults.push_back(std::move(futilityP), "futilityP");
  byStageResults.push_back(std::move(info), "information");
  if (estimateHazardRatio) byStageResults.push_back(std::move(HR), "HR");
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
  settings.push_back(estimateHazardRatio, "estimateHazardRatio");
  settings.push_back(spendingTime, "spendingTime");

  ListCpp byTreatmentCounts;
  byTreatmentCounts.push_back(std::move(nevents1), "numberOfEvents1");
  byTreatmentCounts.push_back(std::move(ndropouts1), "numberOfDropouts1");
  byTreatmentCounts.push_back(std::move(nsubjects1), "numberOfSubjects1");
  byTreatmentCounts.push_back(std::move(nevents2), "numberOfEvents2");
  byTreatmentCounts.push_back(std::move(ndropouts2), "numberOfDropouts2");
  byTreatmentCounts.push_back(std::move(nsubjects2), "numberOfSubjects2");
  byTreatmentCounts.push_back(expectedNumberOfEvents1, "expectedNumberOfEvents1");
  byTreatmentCounts.push_back(expectedNumberOfDropouts1, "expectedNumberOfDropouts1");
  byTreatmentCounts.push_back(expectedNumberOfSubjects1, "expectedNumberOfSubjects1");
  byTreatmentCounts.push_back(expectedNumberOfEvents2, "expectedNumberOfEvents2");
  byTreatmentCounts.push_back(expectedNumberOfDropouts2, "expectedNumberOfDropouts2");
  byTreatmentCounts.push_back(expectedNumberOfSubjects2, "expectedNumberOfSubjects2");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(settings), "settings");
  result.push_back(std::move(byTreatmentCounts), "byTreatmentCounts");

  return result;
}


//' @title Log-Rank Test Power
//' @description Estimates the power, stopping probabilities, and expected
//' sample size in a two-group survival design.
//'
//' @inheritParams param_kMax
//' @param informationRates The information rates in terms of number
//'   of events for the conventional log-rank test and in terms of
//'   the actual information for weighted log-rank tests.
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
//' @inheritParams param_hazardRatioH0
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
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @inheritParams param_estimateHazardRatio
//' @inheritParams param_typeOfComputation
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param studyDuration Study duration for fixed follow-up design.
//'   Defaults to missing, which is to be replaced with the sum of
//'   \code{accrualDuration} and \code{followupTime}. If provided,
//'   the value is allowed to be less than the sum of \code{accrualDuration}
//'   and \code{followupTime}.
//'
//' @return An S3 class \code{lrpower} object with 4 components:
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
//'     - \code{expectedStudyDuration}: The expected study duration.
//'
//'     - \code{expectedInformation}: The expected information.
//'
//'     - \code{accrualDuration}: The accrual duration.
//'
//'     - \code{followupTime}: The follow-up time.
//'
//'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
//'
//'     - \code{rho1}: The first parameter of the Fleming-Harrington family
//'       of weighted log-rank test.
//'
//'     - \code{rho2}: The second parameter of the Fleming-Harrington family
//'       of weighted log-rank test.
//'
//'     - \code{kMax}: The number of stages.
//'
//'     - \code{hazardRatioH0}: The hazard ratio under the null hypothesis.
//'
//'     - \code{typeOfComputation}: The type of computation,
//'       either "direct" for the direct approximation method,
//'       or "schoenfeld" for the Schoenfeld method.
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
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyHR}: The efficacy boundaries on the hazard ratio
//'       scale if \code{estimateHazardRatio}.
//'
//'     - \code{futilityHR}: The futility boundaries on the hazard ratio
//'       scale if \code{estimateHazardRatio}.
//'
//'     - \code{efficacyP}: The efficacy boundaries on the p-value scale.
//'
//'     - \code{futilityP}: The futility boundaries on the p-value scale.
//'
//'     - \code{information}: The cumulative information.
//'
//'     - \code{HR}: The average hazard ratio.
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
//'   \code{estimateHazardRatio}, and \code{spendingTime}.
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
//'     - \code{numberOfEvents2}: The number of events by stage for
//'       the control group.
//'
//'     - \code{numberOfDropouts2}: The number of dropouts by stage for
//'       the control group.
//'
//'     - \code{numberOfSubjects2}: The number of subjects by stage for
//'       the control group.
//'
//'     - \code{expectedNumberOfEvents1}: The expected number of events for
//'       the treatment group.
//'
//'     - \code{expectedNumberOfDropouts1}: The expected number of dropouts
//'       for the treatment group.
//'
//'     - \code{expectedNumberOfSubjects1}: The expected number of subjects
//'       for the treatment group.
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
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survival, and 5% dropout by
//' # the end of 1 year.
//'
//' lrpower(kMax = 2, informationRates = c(0.8, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         allocationRatioPlanned = 1, accrualTime = seq(0, 8),
//'         accrualIntensity = 26/9*seq(1, 9),
//'         piecewiseSurvivalTime = c(0, 6),
//'         lambda1 = c(0.0533, 0.0309),
//'         lambda2 = c(0.0533, 0.0533),
//'         gamma1 = -log(1-0.05)/12,
//'         gamma2 = -log(1-0.05)/12, accrualDuration = 22,
//'         followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List lrpower(
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
    const double hazardRatioH0 = 1,
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
    const bool fixedFollowup = 0,
    const double rho1 = 0,
    const double rho2 = 0,
    const bool estimateHazardRatio = 1,
    const std::string& typeOfComputation = "",
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

  auto out = lrpowercpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    hazardRatioH0, allocationRatioPlanned,
    accrualT, accrualInt,
    piecewiseTime, stratumFrac,
    lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup,
    rho1, rho2, estimateHazardRatio,
    typeOfComputation, spendTime, studyDuration);

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "lrpower";
  return result;
}


double getNeventsFromHazardRatiocpp(
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
    const std::vector<double>& spendingTime,
    const double hazardRatioH0,
    const double hazardRatio,
    const double allocationRatioPlanned,
    const bool rounding) {

  if (beta >= 1.0 - alpha || beta < 0.0001) {
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");
  }

  if (hazardRatioH0 <= 0.0)
    throw std::invalid_argument("hazardRatioH0 must be positive");

  if (std::isnan(hazardRatio))
    throw std::invalid_argument("hazardRatio must be provided");
  if (hazardRatio <= 0.0)
    throw std::invalid_argument("hazardRatio must be positive");

  if (allocationRatioPlanned <= 0.0)
    throw std::invalid_argument("allocationRatioPlanned must be positive");

  // compute theta = | -log(hazardRatio / hazardRatioH0) |
  const double theta = std::fabs(-std::log(hazardRatio / hazardRatioH0));

  // call getDesigncpp with IMax = NaN (solve for IMax given beta)
  ListCpp design = getDesigncpp(
    beta, NaN, theta, kMax, informationRates, efficacyStopping,
    futilityStopping, criticalValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending, futilityBounds,
    typeBetaSpending, parameterBetaSpending, userBetaSpending,
    spendingTime, 1.0 // varianceRatio default 1
  );

  // extract information from overallResults
  DataFrameCpp overall = design.get<DataFrameCpp>("overallResults");
  double maxInformation = overall.get<double>("information")[0];

  // compute required number of events
  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);
  double D = maxInformation / (phi * (1.0 - phi));

  if (rounding) D = std::ceil(D - 1.0e-12);
  return D;
}


//' @title Required Number of Events Given Hazard Ratio
//' @description Obtains the required number of events given the hazard
//' ratios under the null and alternative hypotheses for a group
//' sequential design.
//'
//' @param beta Type II error. Defaults to 0.2.
//' @inheritParams param_kMax
//' @inheritParams param_informationRates
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
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @inheritParams param_hazardRatioH0
//' @param hazardRatio Hazard ratio under the alternative hypothesis
//'   for the active treatment versus control.
//' @inheritParams param_allocationRatioPlanned
//' @param rounding Whether to round up the number of events.
//'   Defaults to 1 for rounding.
//'
//' @return The required number of events.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' getNeventsFromHazardRatio(
//'   beta = 0.2, kMax = 2,
//'   informationRates = c(0.5,1),
//'   alpha = 0.025, typeAlphaSpending = "sfOF",
//'   typeBetaSpending = "sfP",
//'   hazardRatio = 0.673)
//'
//' @export
// [[Rcpp::export]]
double getNeventsFromHazardRatio(
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
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const double hazardRatioH0 = 1,
    const double hazardRatio = NA_REAL,
    const double allocationRatioPlanned = 1,
    const bool rounding = 1) {

  auto infoRates = Rcpp::as<std::vector<double>>(informationRates);
  auto effStopping = convertLogicalVector(efficacyStopping);
  auto futStopping = convertLogicalVector(futilityStopping);
  auto critValues = Rcpp::as<std::vector<double>>(criticalValues);
  auto userAlpha = Rcpp::as<std::vector<double>>(userAlphaSpending);
  auto futBounds = Rcpp::as<std::vector<double>>(futilityBounds);
  auto userBeta = Rcpp::as<std::vector<double>>(userBetaSpending);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  auto D = getNeventsFromHazardRatiocpp(
    beta, kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    userBeta, spendTime, hazardRatioH0, hazardRatio,
    allocationRatioPlanned, rounding);

  return D;
}


ListCpp lrsamplesizecpp(
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
    const double hazardRatioH0,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    double accrualDuration,
    double followupTime,
    const bool fixedFollowup,
    const double rho1,
    const double rho2,
    const bool estimateHazardRatio,
    const std::string& typeOfComputation,
    const std::vector<double>& interval,
    const std::vector<double>& spendingTime,
    const bool rounding = true) {

  if (!std::isnan(alpha) && (alpha < 0.00001 || alpha >= 1)) {
    throw std::invalid_argument("alpha must lie in [0.00001, 1)");
  }

  if (beta >= 1.0 - alpha || beta < 0.0001)
    throw std::invalid_argument("beta must lie in [0.0001, 1-alpha)");

  if (kMax < 1) throw std::invalid_argument("kMax must be a positive integer");

  std::size_t K = static_cast<std::size_t>(kMax);

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
    for (std::size_t i = 0; i < K; ++i)
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
    for (std::size_t i = 0; i < K - 1; ++i) {
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

  if (hazardRatioH0 <= 0.0)
    throw std::invalid_argument("hazardRatioH0 must be positive");

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
  if (!none_na(gamma1)) throw std::invalid_argument("gamma1 must be provided");
  if (!none_na(gamma2)) throw std::invalid_argument("gamma2 must be provided");

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
  std::size_t nstrata = stratumFraction.size();
  std::size_t nintervals = piecewiseSurvivalTime.size();
  std::size_t nsi = nstrata * nintervals;
  auto lambda1x = expand_stratified(lambda1, nstrata, nintervals, "lambda1");
  auto lambda2x = expand_stratified(lambda2, nstrata, nintervals, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintervals, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintervals, "gamma2");

  if (!std::isnan(accrualDuration) && accrualDuration <= 0.0)
    throw std::invalid_argument("accrualDuration must be positive");

  if (!std::isnan(followupTime) && fixedFollowup && followupTime <= 0.0)
    throw std::invalid_argument("followupTime must be positive for fixed follow-up");
  if (!std::isnan(followupTime) && !fixedFollowup && followupTime < 0.0)
    throw std::invalid_argument(
        "followupTime must be non-negative for variable follow-up");
  if (fixedFollowup && std::isnan(followupTime))
    throw std::invalid_argument("followupTime must be provided for fixed follow-up");

  if (rho1 < 0.0) throw std::invalid_argument("rho1 must be non-negative");
  if (rho2 < 0.0) throw std::invalid_argument("rho2 must be non-negative");

  if (interval.size() != 2)
    throw std::invalid_argument("interval must have 2 elements");
  if (interval[0] < 0)
    throw std::invalid_argument("lower limit of interval must be positive");
  if (interval[0] >= interval[1])
    throw std::invalid_argument("upper limit must exceed lower limit for interval");

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

  // --- Determine if Schoenfeld method is eligible-----------------------------
  double hazardRatio = 1.0;
  std::vector<double> hrx(nsi);
  for (std::size_t i = 0; i < nsi; ++i) hrx[i] = lambda1x[i] / lambda2x[i];
  bool proportional = true;
  double hrx0 = hrx[0];
  for (std::size_t i = 1; i < nsi; ++i) {
    if (std::fabs(hrx[i] - hrx0) > 1e-8) { proportional = false; break; }
  }
  if (proportional) hazardRatio = hrx0;

  bool schoenfeld_eligible = proportional && rho1 == 0.0 && rho2 == 0.0;
  std::string su = typeOfComputation;
  if (su.empty()) {
    su = schoenfeld_eligible ? "schoenfeld" : "direct";
  } else {
    for (char &c : su) {
      c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
  }
  char su1 = su.front();
  if (su1 != 'd' && su1 != 's')
    throw std::invalid_argument("typeOfComputation must be direct or schoenfeld");
  if (!schoenfeld_eligible && su1 == 's') {
    throw std::invalid_argument(
        "Schoenfeld method can only be used for ordinary log-rank test "
        "with proportional hazards");
  }

  // --- Efficacy boundaries ---------------------------------------------------
  std::vector<double> l(K, -6.0), zero(K, 0.0);
  std::vector<double> critValues = criticalValues;
  if (missingCriticalValues) {
    bool haybittle = false;
    if (K > 1 && criticalValues.size() == K) {
      bool hasNaN = false;
      for (std::size_t i = 0; i < K - 1; ++i) {
        if (std::isnan(criticalValues[i])) { hasNaN = true; break; }
      }
      if (!hasNaN && std::isnan(criticalValues[K-1])) {
        haybittle = true;
      }
    }

    if (haybittle) { // Haybittle & Peto
      std::vector<double> u(K);
      for (std::size_t i = 0; i < K - 1; ++i) {
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
  double alpha1 = missingCriticalValues ? alpha : cumAlphaSpent[K-1];

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

  double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

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


  // Helper: compute number of events under H1 given accrualDuration (accrDur),
  // followupTime (fu), and accrualIntensity (accrInt).
  auto events_under_H1 = [&](double t, double accrDur, double fu,
                             const std::vector<double>& accrInt) {
    DataFrameCpp lr = lrstat1cpp(
      t, hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrInt,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrDur, fu, fixedFollowup,
      rho1, rho2, true);
    return extract_lr(lr, "nevents");
  };

  // Helper: difference between the calculated cumulative type II error and
  // the pre-specified cumulative beta spending under H1 given
  // accrualDuration (accrDur), followupTime (fu), and accrualIntensity (accrI).
  // NOTE: we cannot directly call lrpowercpp to solve for unknown design
  // parameters since lrpowercpp does not support user-specified beta spending.
  auto betadiff_under_H1 = [&](double accrDur, double fu,
                               const std::vector<double>& accrInt) {

    double studyDuration1 = accrDur + fu;

    // --- obtain theta and information at each look under H1 as in lrpowercpp
    std::vector<double> theta(K), info(K);
    if (rho1 == 0 && rho2 == 0) {
      double vtrt = phi * (1.0 - phi);
      double theta1 = -std::log(hazardRatio / hazardRatioH0);
      bool predictEventOnly = (su1 == 's');

      DataFrameCpp lr = lrstat1cpp(
        studyDuration1, hazardRatioH0, allocationRatioPlanned,
        accrualTime, accrInt,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrDur, fu, fixedFollowup,
        rho1, rho2, predictEventOnly);

      double totalEvents = extract_lr(lr, "nevents");

      if (su1 == 's') {
        theta[K - 1] = theta1;
        info[K - 1] = vtrt * totalEvents;
      } else {
        double uscore = extract_lr(lr, "uscore");
        double vscore = extract_lr(lr, "vscore");
        theta[K - 1] = -uscore / vscore;
        info[K - 1] = vscore;
      }

      for (std::size_t i = 0; i < K - 1; ++i) {
        double nevents_target = totalEvents * infoRates[i];

        double time = caltime1cpp(
          nevents_target, allocationRatioPlanned,
          accrualTime, accrInt,
          piecewiseSurvivalTime, stratumFraction,
          lambda1x, lambda2x, gamma1x, gamma2x,
          accrDur, fu, fixedFollowup);

        DataFrameCpp lr1 = lrstat1cpp(
          time, hazardRatioH0, allocationRatioPlanned,
          accrualTime, accrInt,
          piecewiseSurvivalTime, stratumFraction,
          lambda1x, lambda2x, gamma1x, gamma2x,
          accrDur, fu, fixedFollowup,
          rho1, rho2, predictEventOnly);

        if (su1 == 's') {
          theta[i] = theta1;
          info[i] = vtrt * nevents_target;
        } else {
          double uscore1 = extract_lr(lr1, "uscore");
          double vscore1 = extract_lr(lr1, "vscore");
          theta[i] = -uscore1 / vscore1;
          info[i] = vscore1;
        }
      }
    } else {
      DataFrameCpp lr = lrstat1cpp(
        studyDuration1, hazardRatioH0, allocationRatioPlanned,
        accrualTime, accrInt,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrDur, fu, fixedFollowup,
        rho1, rho2, false);

      double maxInformation = extract_lr(lr, "vscore");

      double uscore = extract_lr(lr, "uscore");
      theta[K - 1] = -uscore / maxInformation;
      info[K - 1] =  maxInformation;

      for (std::size_t i = 0; i < K - 1; ++i) {
        double information1 = maxInformation * infoRates[i];

        // solve for analysis time where total information equals information1
        auto g = [&](double t)->double {
          DataFrameCpp lr1 = lrstat1cpp(
            t, hazardRatioH0, allocationRatioPlanned,
            accrualTime, accrInt,
            piecewiseSurvivalTime, stratumFraction,
            lambda1x, lambda2x, gamma1x, gamma2x,
            accrDur, fu, fixedFollowup,
            rho1, rho2, false);

          return extract_lr(lr1, "vscore") - information1;
        };

        double time = brent(g, 1e-6, studyDuration1, 1e-6);

        DataFrameCpp lr1 = lrstat1cpp(
          time, hazardRatioH0, allocationRatioPlanned,
          accrualTime, accrInt,
          piecewiseSurvivalTime, stratumFraction,
          lambda1x, lambda2x, gamma1x, gamma2x,
          accrDur, fu, fixedFollowup,
          rho1, rho2, false);

        double uscore1 = extract_lr(lr1, "uscore");
        theta[i] = -uscore1 / information1;
        info[i] = information1;
      }
    }

    // --- compute futility bounds and cumulative beta spending under H1
    if (!missingFutilityBounds || bsf == "none" || K == 1) {
      // compute overall beta using the specified futility bounds
      ListCpp probs = exitprobcpp(critValues, futBounds, theta, info);
      auto v = probs.get<std::vector<double>>("exitProbUpper");
      double overallReject = std::accumulate(v.begin(), v.end(), 0.0);
      return (1.0 - overallReject) - beta;
    } else {
      std::vector<double> u1; u1.reserve(K);
      std::vector<double> l1; l1.reserve(K);
      double thetaSqrtI0 = theta[0] * std::sqrt(info[0]);

      // compute cumulative beta spending under H1 similar to getPower
      std::fill(futBounds.begin(), futBounds.end(), -6.0);
      double eps = 0.0;

      // first stage
      double cb = (bsf == "user") ? userBetaSpending[0] :
        errorSpentcpp(spendTime[0], beta, bsf, parameterBetaSpending);

      if (futStopping[0]) {
        eps = boost_pnorm(critValues[0] - thetaSqrtI0) - cb;
        if (eps < 0.0) return -1.0;
        futBounds[0] = boost_qnorm(cb) + thetaSqrtI0;
      }

      // subsequent stages
      for (std::size_t k = 1; k < K; ++k) {
        if (futStopping[k]) {
          cb = (bsf == "user") ? userBetaSpending[k] :
          errorSpentcpp(spendTime[k], beta, bsf, parameterBetaSpending);

          u1.resize(k + 1);
          l1.resize(k + 1);
          std::memcpy(u1.data(), critValues.data(), k * sizeof(double));
          std::memcpy(l1.data(), futBounds.data(), k * sizeof(double));
          u1[k] = 6.0;

          // lambda expression for finding futility bound at stage k
          // g is an increasing function of the futility bound at stage k,
          // and we want to find the root
          auto g = [&](double aval) -> double {
            l1[k] = aval;
            ListCpp probs = exitprobcpp(u1, l1, theta, info);
            auto v = probs.get<std::vector<double>>("exitProbLower");
            double cpl = std::accumulate(v.begin(), v.end(), 0.0);
            return cpl - cb;
          };

          double bk = critValues[k];
          eps = g(bk);
          double g_minus6 = g(-6.0);
          if (g_minus6 > 0.0) { // no beta spent at current visit
            futBounds[k] = -6.0;
          } else if (eps > 0.0) {
            auto g_for_brent = [&](double x)->double {
              if (x == -6.0) return g_minus6;
              if (x == bk) return eps;
              return g(x);
            };
            futBounds[k] = brent(g_for_brent, -6.0, bk, 1e-6);
          } else if (k < K - 1) {
            return -1.0;
          }
        }
      }

      return eps;
    }
  };


  // For special handling of variable follow-up with unknown == FUP_TIME
  // when power at zero follow-up (end of enrollment) already exceeds target,
  // set followupTime = 0 and find minimal accrualDuration achieving target
  bool curtailed = false;

  // Root-finding function to solve for the unknown design parameter to achieve target
  std::function<double(double)> f_root;
  if (su1 == 's') {
    // --- Schoenfeld method if eligible and requested -------------------------
    double theta = -std::log(hazardRatio / hazardRatioH0);

    ListCpp design = getDesigncpp(
      beta, NaN, theta, kMax, infoRates,
      effStopping, futStopping,
      critValues, alpha1, asf,
      parameterAlphaSpending, userAlphaSpending,
      futBounds, bsf, parameterBetaSpending,
      userBetaSpending, spendTime, 1.0);

    auto byStageResults = design.get<DataFrameCpp>("byStageResults");
    futBounds = byStageResults.get<double>("futilityBounds");

    auto overallResults = design.get<DataFrameCpp>("overallResults");
    double maxInformation = overallResults.get<double>("information")[0];
    double D = maxInformation / (phi * (1.0 - phi));

    // Build root function f_root depending on unknown parameter
    if (unknown == ACC_DUR) {
      // find minimal accrualDuration achieving target events
      f_root = [&](double x)->double {
        double t = x + followupTime;
        return events_under_H1(t, x, followupTime, accrualInt) - D;
      };
    } else if (unknown == FUP_TIME) {
      if (!fixedFollowup &&
          events_under_H1(accrualDuration, accrualDuration, 0, accrualInt) > D) {
        thread_utils::push_thread_warning(
          "Events at zero follow-up (end of enrollment) already exceeds target. "
          "Setting followupTime = 0 and finding minimal accrualDuration.");

        f_root = [&](double x)->double {
          return events_under_H1(x, x, 0.0, accrualInt) - D;
        };

        curtailed = true;
      } else {
        // find minimal followupTime achieving target events
         f_root = [&](double x)->double {
          double t = accrualDuration + x;
          return events_under_H1(t, accrualDuration, x, accrualInt) - D;
        };
      }
    } else { // ACC_INT: find multiplier m s.t. scaled accrualInt*m attains power
      f_root = [&](double m)->double {
        std::vector<double> scaled = accrualInt;
        for (double &v : scaled) v *= m;
        double t = accrualDuration + followupTime;
        return events_under_H1(t, accrualDuration, followupTime, scaled) - D;
      };
    }
  } else {
    // --- Direct method: solve for design parameter to achieve target power ---
    if (unknown == ACC_DUR) {
      f_root = [&](double x)->double {
        return betadiff_under_H1(x, followupTime, accrualInt);
      };
    } else if (unknown == FUP_TIME) {
      if (!fixedFollowup &&
          betadiff_under_H1(accrualDuration, 0, accrualInt) < 0.0) {
        thread_utils::push_thread_warning(
          "Power at zero follow-up (end of enrollment) already exceeds target. "
          "Setting followupTime=0 and finding minimal accrualDuration.");
        f_root = [&](double x)->double {
          return betadiff_under_H1(x, 0.0, accrualInt);
        };
        curtailed = true;
      } else {
        f_root = [&](double x)->double {
          return betadiff_under_H1(accrualDuration, x, accrualInt);
        };
      }
    } else { // ACC_INT: find multiplier m s.t. scaled accrualInt*m attains power
      f_root = [&](double m)->double {
        std::vector<double> scaled = accrualInt;
        for (double &v : scaled) v *= m;
        return betadiff_under_H1(accrualDuration, followupTime, scaled);
      };
    }
  }

  // determine root-finding interval
  double lower, upper;
  if (curtailed) {
    // we will search for minimal accrualDuration with followupTime fixed at 0
    lower = 1e-6;
    upper = accrualDuration;
  } else {
    lower = interval[0];
    upper = interval[1];
  }

  // expand upper if needed to ensure root is bracketed
  double fl_val = f_root(lower), fu_val = f_root(upper);
  if (!curtailed) {
    int expand_iter = 0;
    while (fl_val * fu_val > 0.0 && expand_iter < 60) {
      upper *= 2.0;
      fu_val = f_root(upper);
      ++expand_iter;
    }
  }
  if (fl_val * fu_val > 0.0) throw std::runtime_error(
      "Unable to bracket root; check interval or inputs");

  // solve for root and apply solution for the unknown design parameter
  auto f_for_brent = [&](double x)->double {
    if (x == lower) return fl_val;
    if (x == upper) return fu_val;
    return f_root(x);
  };
  double solution = brent(f_for_brent, lower, upper, 1e-6);

  if (unknown == ACC_DUR) accrualDuration = solution;
  else if (unknown == FUP_TIME && !curtailed) followupTime = solution;
  else if (unknown == FUP_TIME && curtailed) {
    followupTime = 0.0;
    accrualDuration = solution;
  } else {
    // scaled multiplier for accrualIntensity
    for (double &v : accrualInt) v *= solution;
  }

  double studyDuration = accrualDuration + followupTime;

  // NOTE: futility bounds are produced as a by-product of the solver
  // for the unknown design parameter.
  // The futility bound must meet the efficacy bound at the final look.
  futBounds[K - 1] = critValues[K - 1];


  ListCpp resultsH1; // results under H1 with final design parameters

  // Rounding adjustments: make integer events/subjects
  if (rounding) {
    // Get number of events at study end under H1
    DataFrameCpp lr_end = lrstat1cpp(
      studyDuration, hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrualInt,
      piecewiseSurvivalTime, stratumFraction,
      lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup,
      rho1, rho2, true);

    // unrounded number of events at study end under H1
    double D0 = extract_lr(lr_end, "nevents");
    // round up to nearest integer (subtract small epsilon to
    // avoid rounding up when D0 is very close to an integer)
    double D = std::ceil(D0 - 1e-12);

    if (!fixedFollowup) { // variable follow-up
      // adjust accrualDuration / accrualIntensity to get integer subjects
      double n0 = extract_lr(lr_end, "subjects");
      double n = std::ceil(n0 - 1e-12);
      if (n - n0 > 1e-6) {
        if (unknown == ACC_INT) {
          // scale accrual intensity
          double scale = n / n0;
          for (double &v : accrualInt) v *= scale;
        } else {
          // solve for accrualDuration that achieves n subjects
          accrualDuration = getAccrualDurationFromN1(n, accrualTime, accrualInt);
        }
      }

      // adjust follow-up time to obtain integer number of events
      auto h_follow = [&](double x)->double {
        double t = accrualDuration + x;
        return events_under_H1(t, accrualDuration, x, accrualInt) - D;
      };

      // bracket and solve for follow-up time
      double lower_h = 0.0, upper_h = std::max(1.1*followupTime, 1.0);
      while (h_follow(upper_h) < 0.0) {
        lower_h = upper_h;
        upper_h *= 2.0;
      }
      followupTime = brent(h_follow, lower_h, upper_h, 1e-6);
      studyDuration = accrualDuration + followupTime;
    } else { // fixed follow-up
      // adjust accrualDuration / accrualIntensity to get integer events
      if (unknown == ACC_INT) {
        double scale = D / D0;
        for (double &v : accrualInt) v *= scale;
      } else {
        // solve for accrualDuration giving D events at fixed followup
        auto h_accr = [&](double x)->double {
          double t = x + followupTime;
          return events_under_H1(t, x, followupTime, accrualInt) - D;
        };

        double lower_h = std::max(1e-6, accrualDuration);
        double upper_h = accrualDuration * 1.1;
        while (h_accr(upper_h) < 0.0) {
          lower_h = upper_h;
          upper_h *= 2.0;
        }
        accrualDuration = brent(h_accr, lower_h, upper_h, 1e-6);
      }

      // recompute integer n
      double n0 = accrual1(accrualDuration, accrualTime, accrualInt, accrualDuration);
      double n = std::ceil(n0 - 1e-12);
      if (n - n0 > 1e-6) {
        if (unknown == ACC_INT) {
          double scale = n / n0;
          for (double &v : accrualInt) v *= scale;
        } else {
          accrualDuration = getAccrualDurationFromN1(n, accrualTime, accrualInt);
        }
      }

      // adjust study duration (no change to followupTime) to get integer events
      auto h_follow2 = [&](double x)->double {
        double t = accrualDuration + x;
        return events_under_H1(t, accrualDuration, followupTime, accrualInt) - D;
      };

      double aval = brent(h_follow2, 0.0, followupTime, 1e-6);
      studyDuration = accrualDuration + aval;
    }

    // update informationRates
    if (rho1 == 0.0 && rho2 == 0.0) {
      for (std::size_t i = 0; i < K - 1; ++i) {
        double nevents = std::floor(D * infoRates[i] + 0.5);
        infoRates[i] = nevents / D;
      }
    } else {
      // obtain maximum information
      DataFrameCpp lr = lrstat1cpp(
        studyDuration, hazardRatioH0, allocationRatioPlanned,
        accrualTime, accrualInt,
        piecewiseSurvivalTime, stratumFraction,
        lambda1x, lambda2x, gamma1x, gamma2x,
        accrualDuration, followupTime, fixedFollowup,
        rho1, rho2, false);

      double maxInformation = extract_lr(lr, "vscore");

      for (std::size_t i = 0; i < K - 1; ++i) {
        // recompute analysis time from information target and then events
        double information1 = maxInformation * infoRates[i];

        auto g_info = [&](double t)->double {
          DataFrameCpp lr1 = lrstat1cpp(
            t, hazardRatioH0, allocationRatioPlanned,
            accrualTime, accrualInt,
            piecewiseSurvivalTime, stratumFraction,
            lambda1x, lambda2x, gamma1x, gamma2x,
            accrualDuration, followupTime, fixedFollowup,
            rho1, rho2, false);

          return extract_lr(lr1, "vscore") - information1;
        };

        double time = brent(g_info, 1e-6, studyDuration, 1e-6);

        // get event count at the calendar time and round up
        DataFrameCpp lr_events = lrstat1cpp(
          time, hazardRatioH0, allocationRatioPlanned,
          accrualTime, accrualInt,
          piecewiseSurvivalTime, stratumFraction,
          lambda1x, lambda2x, gamma1x, gamma2x,
          accrualDuration, followupTime, fixedFollowup,
          rho1, rho2, true);

        double nevents = extract_lr(lr_events, "nevents");
        nevents = std::floor(nevents + 0.5);

        // recompute analysis time from rounded event count
        time = caltime1cpp(
          nevents, allocationRatioPlanned,
          accrualTime, accrualInt,
          piecewiseSurvivalTime, stratumFraction,
          lambda1x, lambda2x, gamma1x, gamma2x,
          accrualDuration, followupTime, fixedFollowup);

        // update informationRates based on new analysis time
        DataFrameCpp lr_info = lrstat1cpp(
          time, hazardRatioH0, allocationRatioPlanned,
          accrualTime, accrualInt,
          piecewiseSurvivalTime, stratumFraction,
          lambda1x, lambda2x, gamma1x, gamma2x,
          accrualDuration, followupTime, fixedFollowup,
          rho1, rho2, false);

        double info = extract_lr(lr_info, "vscore");
        infoRates[i] = info / maxInformation;
      }
    }


    // --- Recalculate boundaries --------------------------------------
    if (bsf != "user") {
      // criticalValues, futilityBounds, spendingTime are all recalculated
      // if not user-specified
      resultsH1 = lrpowercpp(
        kMax, infoRates, effStopping, futStopping,
        criticalValues, alpha1, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futilityBounds, typeBetaSpending, parameterBetaSpending,
        hazardRatioH0, allocationRatioPlanned,
        accrualTime, accrualInt,
        piecewiseSurvivalTime, stratumFraction,
        lambda1, lambda2, gamma1, gamma2,
        accrualDuration, followupTime, fixedFollowup,
        rho1, rho2, estimateHazardRatio,
        su, spendingTime, studyDuration);
    } else {
      // criticalValues and spendingTime are recalculated, but not futilityBounds
      resultsH1 = lrpowercpp(
        kMax, infoRates, effStopping, futStopping,
        criticalValues, alpha1, typeAlphaSpending,
        parameterAlphaSpending, userAlphaSpending,
        futBounds, typeBetaSpending, parameterBetaSpending,
        hazardRatioH0, allocationRatioPlanned,
        accrualTime, accrualInt,
        piecewiseSurvivalTime, stratumFraction,
        lambda1, lambda2, gamma1, gamma2,
        accrualDuration, followupTime, fixedFollowup,
        rho1, rho2, estimateHazardRatio,
        su, spendingTime, studyDuration);
    }
  } else { // no rounding adjustments
    // boundaries are not recalculated since analysis times are unchanged,
    // but we still call lrpowercpp to get results under H1 with final
    // design parameters
    resultsH1 = lrpowercpp(
      kMax, infoRates, effStopping, futStopping,
      critValues, alpha1, typeAlphaSpending,
      parameterAlphaSpending, userAlphaSpending,
      futBounds, typeBetaSpending, parameterBetaSpending,
      hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrualInt,
      piecewiseSurvivalTime, stratumFraction,
      lambda1, lambda2, gamma1, gamma2,
      accrualDuration, followupTime, fixedFollowup,
      rho1, rho2, estimateHazardRatio,
      su, spendingTime, studyDuration);
  }


  // obtain results under H0 by matching the total number of events
  // for conventional log-rank test and maximum information for
  // weighted log-rank tests
  DataFrameCpp byStageResults = resultsH1.get<DataFrameCpp>("byStageResults");
  critValues = byStageResults.get<double>("efficacyBounds");
  futBounds = byStageResults.get<double>("futilityBounds");

  DataFrameCpp overallResults = resultsH1.get<DataFrameCpp>("overallResults");
  double D = overallResults.get<double>("numberOfEvents")[0];
  double maxI = overallResults.get<double>("information")[0];

  // Under H0, we assume the same baseline hazard in the control arm as under H1,
  // and the hazard in the treatment arm is scaled by hazardRatioH0.
  std::vector<double> lambda1H0 = lambda2;
  std::vector<double> lambda1H0x = lambda2x;
  for (double &v : lambda1H0) v *= hazardRatioH0;
  for (double &v : lambda1H0x) v *= hazardRatioH0;

  // Helper: compute number of events under H0 given accrualDuration (accrDur),
  // followupTime (fu), and accrualIntensity (accrInt).
  auto events_under_H0 = [&](double t, double accrDur, double fu,
                             const std::vector<double>& accrInt) {
    DataFrameCpp lr = lrstat1cpp(
      t, hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrInt,
      piecewiseSurvivalTime, stratumFraction,
      lambda1H0x, lambda2x, gamma1x, gamma2x,
      accrDur, fu, fixedFollowup,
      rho1, rho2, true);
    return extract_lr(lr, "nevents");
  };

  // Helper: compute information under H0 given accrualDuration (accrDur),
  // followupTime (fu), and accrualIntensity (accrInt).
  auto info_under_H0 = [&](double t, double accrDur, double fu,
                           const std::vector<double>& accrInt) {
    DataFrameCpp lr = lrstat1cpp(
      t, hazardRatioH0, allocationRatioPlanned,
      accrualTime, accrInt,
      piecewiseSurvivalTime, stratumFraction,
      lambda1H0x, lambda2x, gamma1x, gamma2x,
      accrDur, fu, fixedFollowup,
      rho1, rho2, false);
    return extract_lr(lr, "vscore");
  };


  if (!fixedFollowup) { // variable follow-up
    auto h_follow = [&](double x)->double {
      double t = accrualDuration + x;
      if (rho1 == 0 && rho2 == 0) {
        return events_under_H0(t, accrualDuration, x, accrualInt) - D;
      } else {
        return info_under_H0(t, accrualDuration, x, accrualInt) - maxI;
      }
    };

    double h_0 = h_follow(0.0);
    if (h_0 < 0.0) { // normal case: adjust followupTime
      double lower_h = 0.0, upper_h = std::max(1.0, followupTime);
      double hl_val = h_0, hu_val = h_follow(upper_h);
      int expand_iter = 0;
      while (hu_val < 0.0 && expand_iter < 60) {
        lower_h = upper_h;
        hl_val = hu_val;
        upper_h *= 2.0;
        hu_val = h_follow(upper_h);
        ++expand_iter;
      }
      if (hl_val * hu_val > 0.0) throw std::runtime_error(
          "Unable to bracket root for H0; check interval or inputs");

      auto h_for_brent = [&](double x)->double {
        if (x == lower_h) return hl_val;
        if (x == upper_h) return hu_val;
        return h_follow(x);
      };
      followupTime = brent(h_for_brent, lower_h, upper_h, 1e-6);
      studyDuration = accrualDuration + followupTime;
    } else {
      thread_utils::push_thread_warning(
        "Events at zero follow-up (end of enrollment) already exceeds target. "
        "Setting followupTime = 0 and finding minimal accrualDuration.");
      auto h_accr = [&](double x)->double {
        if (rho1 == 0 && rho2 == 0) {
          return events_under_H0(x, x, 0.0, accrualInt) - D;
        } else {
          return info_under_H0(x, x, 0.0, accrualInt) - maxI;
        }
      };

      accrualDuration = brent(h_accr, 1e-6, accrualDuration, 1e-6);
      followupTime = 0.0;
      studyDuration = accrualDuration + followupTime;
    }
  } else { // fixed follow-up
    auto h_minfu = [&](double x)->double {
      double t = accrualDuration + x;
      if (rho1 == 0 && rho2 == 0) {
        return events_under_H0(t, accrualDuration, followupTime, accrualInt) - D;
      } else {
        return info_under_H0(t, accrualDuration, followupTime, accrualInt) - maxI;
      }
    };

    double h_0 = h_minfu(0.0);
    double h_fu = h_minfu(followupTime);

    if (h_fu < 0.0) { // need to increase accrualDuration
      auto h_accr = [&](double x)->double {
        double t = x + followupTime;
        if (rho1 == 0 && rho2 == 0) {
          return events_under_H0(t, x, followupTime, accrualInt) - D;
        } else {
          return info_under_H0(t, x, followupTime, accrualInt) - maxI;
        }
      };

      double lower_h = std::max(1e-6, accrualDuration);
      double upper_h = accrualDuration * 2.0;
      double hl_val = h_accr(lower_h), hu_val = h_accr(upper_h);
      int expand_iter = 0;
      while (hu_val < 0.0 && expand_iter < 60) {
        lower_h = upper_h;
        hl_val = hu_val;
        upper_h *= 2.0;
        hu_val = h_accr(upper_h);
        ++expand_iter;
      }
      if (hl_val * hu_val > 0.0) throw std::runtime_error(
          "Unable to bracket root for H0; check interval or inputs");

      auto h_for_brent = [&](double x)->double {
        if (x == lower_h) return hl_val;
        if (x == upper_h) return hu_val;
        return h_accr(x);
      };
      accrualDuration = brent(h_for_brent, lower_h, upper_h, 1e-6);
      studyDuration = accrualDuration + followupTime;
    } else if (h_0 < 0.0) {
      // study duration between accrualDuration and accrualDuration + followupTime
       auto h_for_brent = [&](double x)->double {
        if (x == 0.0) return h_0;
        if (x == followupTime) return h_fu;
        return h_minfu(x);
      };
      double aval = brent(h_for_brent, 0.0, followupTime, 1e-6);
      studyDuration = accrualDuration + aval;
    } else {
      thread_utils::push_thread_warning(
        "Events at zero follow-up (end of enrollment) already exceeds target. "
        "finding minimal accrualDuration as followupTime is fixed.");
       auto h_accr = [&](double x)->double {
        if (rho1 == 0 && rho2 == 0) {
          return events_under_H0(x, x, followupTime, accrualInt) - D;
        } else {
          return info_under_H0(x, x, followupTime, accrualInt) - maxI;
        }
      };

      accrualDuration = brent(h_accr, 1e-6, accrualDuration, 1e-6);
      studyDuration = accrualDuration + followupTime;
    }
  }

  ListCpp resultsH0 = lrpowercpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha1, typeAlphaSpending,
    parameterAlphaSpending, userAlphaSpending,
    futBounds, typeBetaSpending, parameterBetaSpending,
    hazardRatioH0, allocationRatioPlanned,
    accrualTime, accrualInt,
    piecewiseSurvivalTime, stratumFraction,
    lambda1H0, lambda2, gamma1, gamma2,
    accrualDuration, followupTime, fixedFollowup,
    rho1, rho2, estimateHazardRatio,
    su, spendingTime, studyDuration);

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


//' @title Log-Rank Test Sample Size
//' @description Obtains the needed accrual duration given power and
//' follow-up time, the needed follow-up time given power and
//' accrual duration, or the needed absolute accrual rates given
//' power, accrual duration, follow-up time, and relative accrual
//' rates in a two-group survival design.
//'
//' @param beta Type II error. Defaults to 0.2.
//' @inheritParams param_kMax
//' @param informationRates The information rates in terms of number
//'   of events for the conventional log-rank test and in terms of
//'   the actual information for weighted log-rank tests.
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
//' @inheritParams param_hazardRatioH0
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
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//' @inheritParams param_estimateHazardRatio
//' @inheritParams param_typeOfComputation
//' @param interval The interval to search for the solution of
//'   accrualDuration, followupTime, or the proportionality constant
//'   of accrualIntensity. Defaults to \code{c(0.001, 240)}. Adjustment
//'   may be needed for non-monotone relationship with study power.
//' @param spendingTime A vector of length \code{kMax} for the error spending
//'   time at each analysis. Defaults to missing, in which case, it is the
//'   same as \code{informationRates}.
//' @param rounding Whether to round up sample size and events.
//'   Defaults to 1 for sample size rounding.
//'
//' @return A list of two components:
//'
//' * \code{resultsUnderH1}: An S3 class \code{lrpower} object under the
//' alternative hypothesis.
//'
//' * \code{resultsUnderH0}: An S3 class \code{lrpower} object under the
//' null hypothesis.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{lrpower}}
//'
//' @examples
//' # Piecewise accrual, piecewise exponential survival, and 5% dropout by
//' # the end of 1 year.
//'
//' # Example 1: Obtains accrual duration given power and follow-up time
//'
//' lrsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              lambda1 = c(0.0533, 0.0309),
//'              lambda2 = c(0.0533, 0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12,
//'              accrualDuration = NA,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//'
//' # Example 2: Obtains follow-up time given power and accrual duration
//'
//' lrsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              lambda1 = c(0.0533, 0.0309),
//'              lambda2 = c(0.0533, 0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12,
//'              accrualDuration = 22,
//'              followupTime = NA, fixedFollowup = FALSE)
//'
//'
//' # Example 3: Obtains absolute accrual intensity given power,
//' # accrual duration, follow-up time, and relative accrual intensity
//'
//' lrsamplesize(beta = 0.2, kMax = 2,
//'              informationRates = c(0.8, 1),
//'              alpha = 0.025, typeAlphaSpending = "sfOF",
//'              accrualTime = seq(0, 8),
//'              accrualIntensity = 26/9*seq(1, 9),
//'              piecewiseSurvivalTime = c(0, 6),
//'              lambda1 = c(0.0533, 0.0309),
//'              lambda2 = c(0.0533, 0.0533),
//'              gamma1 = -log(1-0.05)/12,
//'              gamma2 = -log(1-0.05)/12,
//'              accrualDuration = 22,
//'              followupTime = 18, fixedFollowup = FALSE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List lrsamplesize(
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
    const double hazardRatioH0 = 1,
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
    const double rho1 = 0,
    const double rho2 = 0,
    const bool estimateHazardRatio = true,
    const std::string& typeOfComputation = "",
    const Rcpp::NumericVector& interval = Rcpp::NumericVector::create(0.001, 240),
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
  auto piecewiseTime = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);
  auto intv = Rcpp::as<std::vector<double>>(interval);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  auto out = lrsamplesizecpp(
    beta, kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    userBeta, hazardRatioH0, allocationRatioPlanned,
    accrualT, accrualInt,
    piecewiseTime, stratumFrac,
    lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup,
    rho1, rho2, estimateHazardRatio,
    typeOfComputation, intv, spendTime, rounding);

  ListCpp resultsUnderH1 = out.get_list("resultsUnderH1");
  ListCpp resultsUnderH0 = out.get_list("resultsUnderH0");

  Rcpp::List resultsH1 = Rcpp::wrap(resultsUnderH1);
  Rcpp::List resultsH0 = Rcpp::wrap(resultsUnderH0);

  resultsH1.attr("class") = "lrpower";
  resultsH0.attr("class") = "lrpower";

  return Rcpp::List::create(
    Rcpp::Named("resultsUnderH1") = resultsH1,
    Rcpp::Named("resultsUnderH0") = resultsH0
  );
}

//
//
// //' @title Power for Equivalence in Hazard Ratio
// //' @description Obtains the power for equivalence in hazard ratio.
// //'
// //' @inheritParams param_kMax
// //' @param informationRates The information rates.
// //'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
// //' @inheritParams param_criticalValues
// //' @param alpha The significance level for each of the two one-sided
// //'   tests. Defaults to 0.05.
// //' @inheritParams param_typeAlphaSpending
// //' @inheritParams param_parameterAlphaSpending
// //' @inheritParams param_userAlphaSpending
// //' @param hazardRatioLower The lower equivalence limit of hazard ratio.
// //' @param hazardRatioUpper The upper equivalence limit of hazard ratio.
// //' @inheritParams param_allocationRatioPlanned
// //' @inheritParams param_accrualTime
// //' @inheritParams param_accrualIntensity
// //' @inheritParams param_piecewiseSurvivalTime
// //' @inheritParams param_stratumFraction
// //' @inheritParams param_lambda1_stratified
// //' @inheritParams param_lambda2_stratified
// //' @inheritParams param_gamma1_stratified
// //' @inheritParams param_gamma2_stratified
// //' @inheritParams param_accrualDuration
// //' @inheritParams param_followupTime
// //' @inheritParams param_fixedFollowup
// //' @inheritParams param_typeOfComputation
// //' @param spendingTime A vector of length \code{kMax} for the error spending
// //'   time at each analysis. Defaults to missing, in which case, it is the
// //'   same as \code{informationRates}.
// //' @param studyDuration Study duration for fixed follow-up design.
// //'   Defaults to missing, which is to be replaced with the sum of
// //'   \code{accrualDuration} and \code{followupTime}. If provided,
// //'   the value is allowed to be less than the sum of \code{accrualDuration}
// //'   and \code{followupTime}.
// //'
// //' @return An S3 class \code{lrpowerequiv} object with 4 components:
// //'
// //' * \code{overallResults}: A data frame containing the following variables:
// //'
// //'     - \code{overallReject}: The overall rejection probability.
// //'
// //'     - \code{alpha}: The overall significance level.
// //'
// //'     - \code{numberOfEvents}: The total number of events.
// //'
// //'     - \code{numberOfDropouts}: The total number of dropouts.
// //'
// //'     - \code{numbeOfSubjects}: The total number of subjects.
// //'
// //'     - \code{studyDuration}: The total study duration.
// //'
// //'     - \code{information}: The maximum information.
// //'
// //'     - \code{expectedNumberOfEvents}: The expected number of events.
// //'
// //'     - \code{expectedNumberOfDropouts}: The expected number of dropouts.
// //'
// //'     - \code{expectedNumberOfSubjects}: The expected number of subjects.
// //'
// //'     - \code{expectedStudyDuration}: The expected study duration.
// //'
// //'     - \code{expectedInformation}: The expected information.
// //'
// //'     - \code{kMax}: The number of stages.
// //'
// //'     - \code{hazardRatioLower}: The lower equivalence limit of hazard
// //'       ratio.
// //'
// //'     - \code{hazardRatioUpper}: The upper equivalence limit of hazard
// //'       ratio.
// //'
// //'     - \code{accrualDuration}: The accrual duration.
// //'
// //'     - \code{followupTime}: The follow-up time.
// //'
// //'     - \code{fixedFollowup}: Whether a fixed follow-up design is used.
// //'
// //' * \code{byStageResults}: A data frame containing the following variables:
// //'
// //'     - \code{informationRates}: The information rates.
// //'
// //'     - \code{efficacyBounds}: The efficacy boundaries on the Z-scale for
// //'       each of the two one-sided tests.
// //'
// //'     - \code{rejectPerStage}: The probability for efficacy stopping.
// //'
// //'     - \code{cumulativeRejection}: The cumulative probability for efficacy
// //'       stopping.
// //'
// //'     - \code{cumulativeAlphaSpent}: The cumulative alpha for each of
// //'       the two one-sided tests.
// //'
// //'     - \code{cumulativeAttainedAlphaH10}: The cumulative alpha attained
// //'       under \code{H10}.
// //'
// //'     - \code{cumulativeAttainedAlphaH20}: The cumulative alpha attained
// //'       under \code{H20}.
// //'
// //'     - \code{numberOfEvents}: The number of events.
// //'
// //'     - \code{numberOfDropouts}: The number of dropouts.
// //'
// //'     - \code{numberOfSubjects}: The number of subjects.
// //'
// //'     - \code{analysisTime}: The average time since trial start.
// //'
// //'     - \code{efficacyHRLower}: The efficacy boundaries on the
// //'       hazard ratio scale for the one-sided null hypothesis
// //'       at the lower equivalence limit.
// //'
// //'     - \code{efficacyHRUpper}: The efficacy boundaries on the
// //'       hazard ratio scale for the one-sided null hypothesis
// //'       at the upper equivalence limit.
// //'
// //'     - \code{efficacyP}: The efficacy bounds on the p-value scale for
// //'       each of the two one-sided tests.
// //'
// //'     - \code{information}: The cumulative information.
// //'
// //'     - \code{HR}: The average hazard ratio.
// //'
// //' * \code{settings}: A list containing the following input parameters:
// //'   \code{typeAlphaSpending}, \code{parameterAlphaSpending},
// //'   \code{userAlphaSpending}, \code{allocationRatioPlanned},
// //'   \code{accrualTime}, \code{accuralIntensity},
// //'   \code{piecewiseSurvivalTime}, \code{stratumFraction},
// //'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
// //'   \code{typeOfComputation}, and \code{spendingTime}.
// //'
// //' * \code{byTreatmentCounts}: A list containing the following counts by
// //'   treatment group:
// //'
// //'     - \code{numberOfEvents1}: The number of events by stage for
// //'       the treatment group.
// //'
// //'     - \code{numberOfDropouts1}: The number of dropouts by stage for
// //'       the treatment group.
// //'
// //'     - \code{numberOfSubjects1}: The number of subjects by stage for
// //'       the treatment group.
// //'
// //'     - \code{numberOfEvents2}: The number of events by stage for
// //'       the control group.
// //'
// //'     - \code{numberOfDropouts2}: The number of dropouts by stage for
// //'       the control group.
// //'
// //'     - \code{numberOfSubjects2}: The number of subjects by stage for
// //'       the control group.
// //'
// //'     - \code{expectedNumberOfEvents1}: The expected number of events for
// //'       the treatment group.
// //'
// //'     - \code{expectedNumberOfDropouts1}: The expected number of dropouts
// //'       for the treatment group.
// //'
// //'     - \code{expectedNumberOfSubjects1}: The expected number of subjects
// //'       for the treatment group.
// //'
// //'     - \code{expectedNumberOfEvents2}: The expected number of events for
// //'       control group.
// //'
// //'     - \code{expectedNumberOfDropouts2}: The expected number of dropouts
// //'       for the control group.
// //'
// //'     - \code{expectedNumberOfSubjects2}: The expected number of subjects
// //'       for the control group.
// //'
// //' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
// //'
// //' @seealso \code{\link{rmstat}}
// //'
// //' @examples
// //'
// //' lrpowerequiv(kMax = 2, informationRates = c(0.5, 1),
// //'              alpha = 0.05, typeAlphaSpending = "sfOF",
// //'              hazardRatioLower = 0.71, hazardRatioUpper = 1.4,
// //'              allocationRatioPlanned = 1, accrualTime = seq(0, 8),
// //'              accrualIntensity = 100/9*seq(1, 9),
// //'              piecewiseSurvivalTime = c(0, 6),
// //'              lambda1 = c(0.0533, 0.0533),
// //'              lambda2 = c(0.0533, 0.0533),
// //'              gamma1 = -log(1-0.05)/12,
// //'              gamma2 = -log(1-0.05)/12, accrualDuration = 22,
// //'              followupTime = 18, fixedFollowup = FALSE)
// //'
// //' @export
// // [[Rcpp::export]]
// List lrpowerequiv(const int kMax = 1,
//                   const NumericVector& informationRates = NA_REAL,
//                   const NumericVector& criticalValues = NA_REAL,
//                   const double alpha = 0.05,
//                   const std::string typeAlphaSpending = "sfOF",
//                   const double parameterAlphaSpending = NA_REAL,
//                   const NumericVector& userAlphaSpending = NA_REAL,
//                   const double hazardRatioLower = NA_REAL,
//                   const double hazardRatioUpper = NA_REAL,
//                   const double allocationRatioPlanned = 1,
//                   const NumericVector& accrualTime = 0,
//                   const NumericVector& accrualIntensity = NA_REAL,
//                   const NumericVector& piecewiseSurvivalTime = 0,
//                   const NumericVector& stratumFraction = 1,
//                   const NumericVector& lambda1 = NA_REAL,
//                   const NumericVector& lambda2 = NA_REAL,
//                   const NumericVector& gamma1 = 0,
//                   const NumericVector& gamma2 = 0,
//                   const double accrualDuration = NA_REAL,
//                   const double followupTime = NA_REAL,
//                   const bool fixedFollowup = 0,
//                   const std::string typeOfComputation = "direct",
//                   const NumericVector& spendingTime = NA_REAL,
//                   const double studyDuration = NA_REAL) {
//
//   NumericVector informationRates1 = clone(informationRates);
//   NumericVector criticalValues1 = clone(criticalValues);
//   NumericVector spendingTime1 = clone(spendingTime);
//
//   std::string asf = typeAlphaSpending;
//   std::for_each(asf.begin(), asf.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//
//   double asfpar = parameterAlphaSpending;
//
//   int nstrata = static_cast<int>(stratumFraction.size());
//   int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
//   int nsi = nstrata*nintervals;
//
//   if (kMax < 1) {
//     stop("kMax must be a positive integer");
//   }
//
//   if (is_false(any(is_na(informationRates)))) {
//     if (informationRates.size() != kMax) {
//       stop("Invalid length for informationRates");
//     } else if (informationRates[0] <= 0) {
//       stop("Elements of informationRates must be positive");
//     } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
//       stop("Elements of informationRates must be increasing");
//     } else if (informationRates[kMax-1] != 1) {
//       stop("informationRates must end with 1");
//     }
//   } else {
//     IntegerVector tem = seq_len(kMax);
//     informationRates1 = NumericVector(tem)/(kMax+0.0);
//   }
//
//   if (is_false(any(is_na(criticalValues)))) {
//     if (criticalValues.size() != kMax) {
//       stop("Invalid length for criticalValues");
//     }
//   }
//
//   if (!std::isnan(alpha)) {
//     if (alpha < 0.00001 || alpha >= 1) {
//       stop("alpha must lie in [0.00001, 1)");
//     }
//   }
//
//   if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
//     stop("alpha must be provided when criticalValues is missing");
//   }
//
//   if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
//       asf=="wt" || asf=="sfof" || asf=="sfp" ||
//       asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
//     stop("Invalid value for typeAlphaSpending");
//   }
//
//   if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
//     stop("Missing value for parameterAlphaSpending");
//   }
//
//   if (asf=="sfkd" && asfpar <= 0) {
//     stop ("parameterAlphaSpending must be positive for sfKD");
//   }
//
//   if (is_true(any(is_na(criticalValues))) && asf=="user") {
//     if (is_true(any(is_na(userAlphaSpending)))) {
//       stop("userAlphaSpending must be specified");
//     } else if (userAlphaSpending.size() < kMax) {
//       stop("Insufficient length of userAlphaSpending");
//     } else if (userAlphaSpending[0] < 0) {
//       stop("Elements of userAlphaSpending must be nonnegative");
//     } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
//       stop("Elements of userAlphaSpending must be nondecreasing");
//     } else if (userAlphaSpending[kMax-1] != alpha) {
//       stop("userAlphaSpending must end with specified alpha");
//     }
//   }
//
//   if (std::isnan(hazardRatioLower)) {
//     stop("hazardRatioLower must be provided");
//   }
//
//   if (std::isnan(hazardRatioUpper)) {
//     stop("hazardRatioUpper must be provided");
//   }
//
//   if (hazardRatioLower <= 0) {
//     stop("hazardRatioLower must be positive");
//   }
//
//   if (hazardRatioLower >= hazardRatioUpper) {
//     stop("hazardRatioLower must be less than hazardRatioUpper");
//   }
//
//   if (allocationRatioPlanned <= 0) {
//     stop("allocationRatioPlanned must be positive");
//   }
//
//   if (accrualTime[0] != 0) {
//     stop("accrualTime must start with 0");
//   }
//
//   if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
//     stop("accrualTime should be increasing");
//   }
//
//   if (is_true(any(is_na(accrualIntensity)))) {
//     stop("accrualIntensity must be provided");
//   }
//
//   if (accrualTime.size() != accrualIntensity.size()) {
//     stop("accrualTime must have the same length as accrualIntensity");
//   }
//
//   if (is_true(any(accrualIntensity < 0))) {
//     stop("accrualIntensity must be non-negative");
//   }
//
//   if (piecewiseSurvivalTime[0] != 0) {
//     stop("piecewiseSurvivalTime must start with 0");
//   }
//
//   if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
//     stop("piecewiseSurvivalTime should be increasing");
//   }
//
//   if (is_true(any(stratumFraction <= 0))) {
//     stop("stratumFraction must be positive");
//   }
//
//   if (sum(stratumFraction) != 1) {
//     stop("stratumFraction must sum to 1");
//   }
//
//   if (is_true(any(is_na(lambda1)))) {
//     stop("lambda1 must be provided");
//   }
//
//   if (is_true(any(is_na(lambda2)))) {
//     stop("lambda2 must be provided");
//   }
//
//   if (is_true(any(lambda1 < 0))) {
//     stop("lambda1 must be non-negative");
//   }
//
//   if (is_true(any(lambda2 < 0))) {
//     stop("lambda2 must be non-negative");
//   }
//
//   if (is_true(any(gamma1 < 0))) {
//     stop("gamma1 must be non-negative");
//   }
//
//   if (is_true(any(gamma2 < 0))) {
//     stop("gamma2 must be non-negative");
//   }
//
//   if (lambda1.size() != 1 && lambda1.size() != nintervals &&
//       lambda1.size() != nsi) {
//     stop("Invalid length for lambda1");
//   }
//
//   if (lambda2.size() != 1 && lambda2.size() != nintervals &&
//       lambda2.size() != nsi) {
//     stop("Invalid length for lambda2");
//   }
//
//   if (gamma1.size() != 1 && gamma1.size() != nintervals &&
//       gamma1.size() != nsi) {
//     stop("Invalid length for gamma1");
//   }
//
//   if (gamma2.size() != 1 && gamma2.size() != nintervals &&
//       gamma2.size() != nsi) {
//     stop("Invalid length for gamma2");
//   }
//
//   if (std::isnan(accrualDuration)) {
//     stop("accrualDuration must be provided");
//   }
//
//   if (accrualDuration <= 0) {
//     stop("accrualDuration must be positive");
//   }
//
//   if (std::isnan(followupTime)) {
//     stop("followupTime must be provided");
//   }
//
//   if (fixedFollowup && followupTime <= 0) {
//     stop("followupTime must be positive for fixed follow-up");
//   }
//
//   if (!fixedFollowup && followupTime < 0) {
//     stop("followupTime must be non-negative for variable follow-up");
//   }
//
//   if (fixedFollowup && !std::isnan(studyDuration) &&
//       studyDuration < accrualDuration) {
//     stop("studyDuration must be greater than or equal to accrualDuration");
//   }
//
//   if (fixedFollowup && !std::isnan(studyDuration) &&
//       studyDuration > accrualDuration + followupTime) {
//     stop("studyDuration cannot exceed accrualDuration + followupTime");
//   }
//
//   std::string su = typeOfComputation;
//   std::for_each(su.begin(), su.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//
//   char su1 = su[0];
//
//   if (su1 != 'd' && su1 != 's') {
//     stop("typeOfComputation must be direct or Schoenfeld");
//   }
//
//   double hazardRatio = 1.0;
//   if (su1 == 's') {
//     NumericVector lambda1x = rep(lambda1, nsi/lambda1.size());
//     NumericVector lambda2x = rep(lambda2, nsi/lambda2.size());
//     NumericVector hrx = lambda1x / lambda2x;
//
//     bool proportionalHazards = 1;
//     for (int i=1; i<nsi; ++i) {
//       if (fabs(hrx[i] - hrx[0]) > 1e-8) {
//         proportionalHazards = 0;
//         break;
//       }
//     }
//
//     if (!proportionalHazards) {
//       stop("Schoenfeld method can only be used for proportional hazards");
//     } else {
//       hazardRatio = hrx[0];
//     }
//   }
//
//   if (is_false(any(is_na(spendingTime)))) {
//     if (spendingTime.size() != kMax) {
//       stop("Invalid length for spendingTime");
//     } else if (spendingTime[0] <= 0) {
//       stop("Elements of spendingTime must be positive");
//     } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
//       stop("Elements of spendingTime must be increasing");
//     } else if (spendingTime[kMax-1] != 1) {
//       stop("spendingTime must end with 1");
//     }
//   } else {
//     spendingTime1 = clone(informationRates1);
//   }
//
//   if (fixedFollowup && !std::isnan(studyDuration) &&
//       studyDuration < accrualDuration) {
//     stop("studyDuration must be greater than or equal to accrualDuration");
//   }
//
//   if (fixedFollowup && !std::isnan(studyDuration) &&
//       studyDuration > accrualDuration + followupTime) {
//     stop("studyDuration cannot exceed accrualDuration + followupTime");
//   }
//
//
//   // obtain criticalValues
//   if (is_true(any(is_na(criticalValues)))) {
//     if (kMax > 1 && criticalValues.size() == kMax &&
//         is_false(any(is_na(head(criticalValues, kMax-1)))) &&
//         std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto
//
//       auto f = [kMax, informationRates1,
//                 criticalValues, alpha](double aval)->double {
//                   NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
//                   for (int i=0; i<kMax-1; ++i) {
//                     u[i] = criticalValues[i];
//                   }
//                   u[kMax-1] = aval;
//
//                   List probs = exitprobcpp(u, l, zero, informationRates1);
//                   double cpu = sum(NumericVector(probs[0]));
//                   return cpu - alpha;
//                 };
//
//       criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
//     } else {
//       LogicalVector efficacyStopping1(kMax, 1);
//       criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
//                                     asf, asfpar, userAlphaSpending,
//                                     spendingTime1, efficacyStopping1);
//     }
//   }
//
//   NumericVector efficacyP(kMax);
//   for (int i=0; i<kMax; ++i) {
//     efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
//   }
//
//   NumericVector li(kMax, -6.0), ui(kMax, 6.0), zero(kMax);
//   List probs = exitprobcpp(criticalValues1, li, zero, informationRates1);
//   NumericVector cumAlphaSpent = cumsum(NumericVector(probs[0]));
//
//   // obtain the study duration
//   double studyDuration1 = studyDuration;
//   if (!fixedFollowup || std::isnan(studyDuration)) {
//     studyDuration1 = accrualDuration + followupTime;
//   }
//
//   // obtain the timing of interim analysis
//   NumericVector e0(kMax), time(kMax);
//   NumericVector u0(1, studyDuration1);
//   DataFrame lr = lrstat(u0, 1, allocationRatioPlanned,
//                         accrualTime, accrualIntensity,
//                         piecewiseSurvivalTime, stratumFraction,
//                         lambda1, lambda2, gamma1, gamma2,
//                         accrualDuration, followupTime, fixedFollowup,
//                         0, 0, 1);
//
//   e0 = sum(NumericVector(lr[2]))*informationRates1;
//   if (kMax > 1) {
//     NumericVector nevents1 = e0[Range(0, kMax-2)];
//     NumericVector time1 = caltime(
//       nevents1, allocationRatioPlanned,
//       accrualTime, accrualIntensity,
//       piecewiseSurvivalTime, stratumFraction,
//       lambda1, lambda2, gamma1, gamma2,
//       accrualDuration, followupTime, fixedFollowup);
//     time[Range(0, kMax-2)] = time1;
//   }
//   time[kMax-1] = studyDuration1;
//
//   double phi = allocationRatioPlanned/(1.0 + allocationRatioPlanned);
//   NumericVector HR(kMax), theta(kMax), I(kMax);
//
//   if (su1 == 's') {
//     HR = rep(hazardRatio, kMax);
//     theta = log(HR);
//     I = phi*(1-phi)*e0;
//
//     lr = lrstat(time, 1, allocationRatioPlanned,
//                 accrualTime, accrualIntensity,
//                 piecewiseSurvivalTime, stratumFraction,
//                 lambda1, lambda2, gamma1, gamma2,
//                 accrualDuration, followupTime, fixedFollowup, 0, 0, 1);
//   } else {
//     lr = lrstat(time, 1, allocationRatioPlanned,
//                 accrualTime, accrualIntensity,
//                 piecewiseSurvivalTime, stratumFraction,
//                 lambda1, lambda2, gamma1, gamma2,
//                 accrualDuration, followupTime, fixedFollowup, 0, 0, 3);
//
//     HR = NumericVector(lr[15]);
//     theta = log(HR);
//     I = 1.0/NumericVector(lr[16]);
//   }
//
//   NumericVector nsubjects = NumericVector(lr[1]);
//   NumericVector nsubjects1 = phi*nsubjects;
//   NumericVector nsubjects2 = (1-phi)*nsubjects;
//   NumericVector nevents = NumericVector(lr[2]);
//   NumericVector nevents1 = NumericVector(lr[3]);
//   NumericVector nevents2 = NumericVector(lr[4]);
//   NumericVector ndropouts = NumericVector(lr[5]);
//   NumericVector ndropouts1 = NumericVector(lr[6]);
//   NumericVector ndropouts2 = NumericVector(lr[7]);
//
//
//   // calculate cumulative rejection probability under H1
//   NumericVector theta10 = rep(log(hazardRatioLower), kMax);
//   NumericVector theta20 = rep(log(hazardRatioUpper), kMax);
//   NumericVector b = criticalValues1;
//   NumericVector l = b + theta10*sqrt(I);
//   NumericVector u = -b + theta20*sqrt(I);
//
//   List probs1 = exitprobcpp(pmax(l, li), li, theta, I);
//   List probs2 = exitprobcpp(ui, pmin(u, ui), theta, I);
//
//   NumericVector cpl = cumsum(NumericVector(probs1[0]));
//   NumericVector cpu = cumsum(NumericVector(probs2[1]));
//
//   // identify the last look with l[k] >= u[k] if it exists
//   IntegerVector k = which(l >= u);
//   NumericVector cp(kMax);
//   if (k.size() == 0) {
//     cp = cpl + cpu - 1;
//   } else {
//     int K = max(k);
//     IntegerVector idx = Range(0, K);
//     List a = exitprobcpp(l[idx], u[idx], theta[idx], I[idx]);
//     NumericVector ca = cumsum(NumericVector(a[0]) +
//       NumericVector(a[1]));
//
//     for (int i=0; i<kMax; ++i) {
//       if (i <= K) {
//         cp[i] = cpl[i] + cpu[i] - ca[i];
//       } else {
//         cp[i] = cpl[i] + cpu[i] - 1;
//       }
//     }
//   }
//
//   // incremental exit probabilities under H1
//   NumericVector q(kMax);
//   for (int i=0; i<kMax; ++i) {
//     if (i==0) {
//       q[i] = cp[i];
//     } else if (i<kMax-1) {
//       q[i] = cp[i] - cp[i-1];
//     } else {
//       q[i] = 1 - cp[i-1];
//     }
//   }
//
//   NumericVector rejectPerStage(kMax);
//   for (int i=0; i<kMax; ++i) {
//     if (i==0) {
//       rejectPerStage[i] = cp[i];
//     } else {
//       rejectPerStage[i] = cp[i] - cp[i-1];
//     }
//   }
//
//   NumericVector efficacyHRLower = exp(theta10 + b/sqrt(I));
//   NumericVector efficacyHRUpper = exp(theta20 - b/sqrt(I));
//
//   // calculate cumulative rejection under H10
//   List probs2H10 = exitprobcpp(ui, pmin(u, ui), theta10, I);
//
//   NumericVector cplH10 = cumAlphaSpent;
//   NumericVector cpuH10 = cumsum(NumericVector(probs2H10[1]));
//
//   NumericVector cpH10(kMax);
//   if (k.size() == 0) {
//     cpH10 = cplH10 + cpuH10 - 1;
//   } else {
//     int K = max(k);
//     IntegerVector idx = Range(0, K);
//     List a = exitprobcpp(l[idx], u[idx], theta10[idx], I[idx]);
//     NumericVector ca = cumsum(NumericVector(a[0]) +
//       NumericVector(a[1]));
//
//     for (int i=0; i<kMax; ++i) {
//       if (i <= K) {
//         cpH10[i] = cplH10[i] + cpuH10[i] - ca[i];
//       } else {
//         cpH10[i] = cplH10[i] + cpuH10[i] - 1;
//       }
//     }
//   }
//
//   // calculate cumulative rejection under H20
//   List probs1H20 = exitprobcpp(pmax(l, li), li, theta20, I);
//
//   NumericVector cplH20 = cumsum(NumericVector(probs1H20[0]));
//   NumericVector cpuH20 = cumAlphaSpent;
//
//   NumericVector cpH20(kMax);
//   if (k.size() == 0) {
//     cpH20 = cplH20 + cpuH20 - 1;
//   } else {
//     int K = max(k);
//     IntegerVector idx = Range(0, K);
//     List a = exitprobcpp(l[idx], u[idx], theta20[idx], I[idx]);
//     NumericVector ca = cumsum(NumericVector(a[0]) +
//       NumericVector(a[1]));
//
//     for (int i=0; i<kMax; ++i) {
//       if (i <= K) {
//         cpH20[i] = cplH20[i] + cpuH20[i] - ca[i];
//       } else {
//         cpH20[i] = cplH20[i] + cpuH20[i] - 1;
//       }
//     }
//   }
//
//   double overallReject = cp[kMax-1];
//   double expectedNumberOfEvents = sum(q*nevents);
//   double expectedNumberOfDropouts = sum(q*ndropouts);
//   double expectedNumberOfSubjects = sum(q*nsubjects);
//   double expectedNumberOfEvents1 = sum(q*nevents1);
//   double expectedNumberOfDropouts1 = sum(q*ndropouts1);
//   double expectedNumberOfSubjects1 = sum(q*nsubjects1);
//   double expectedNumberOfEvents2 = sum(q*nevents2);
//   double expectedNumberOfDropouts2 = sum(q*ndropouts2);
//   double expectedNumberOfSubjects2 = sum(q*nsubjects2);
//   double expectedStudyDuration = sum(q*time);
//   double expectedInformation = sum(q*I);
//
//   DataFrame overallResults = DataFrame::create(
//     _["overallReject"] = overallReject,
//     _["alpha"] = alpha,
//     _["numberOfEvents"] = (nevents[kMax-1]),
//     _["numberOfDropouts"] = (ndropouts[kMax-1]),
//     _["numberOfSubjects"] = (nsubjects[kMax-1]),
//     _["studyDuration"] = (time[kMax-1]),
//     _["information"] = (I[kMax-1]),
//     _["expectedNumberOfEvents"] = expectedNumberOfEvents,
//     _["expectedNumberOfDropouts"] = expectedNumberOfDropouts,
//     _["expectedNumberOfSubjects"] = expectedNumberOfSubjects,
//     _["expectedStudyDuration"] = expectedStudyDuration,
//     _["expectedInformation"] = expectedInformation,
//     _["kMax"] = kMax,
//     _["hazardRatioLower"] = hazardRatioLower,
//     _["hazardRatioUpper"] = hazardRatioUpper,
//     _["accrualDuration"] = accrualDuration,
//     _["followupTime"] = followupTime,
//     _["fixedFollowup"] = fixedFollowup);
//
//   DataFrame byStageResults = DataFrame::create(
//     _["informationRates"] = informationRates1,
//     _["efficacyBounds"] = criticalValues1,
//     _["rejectPerStage"] = rejectPerStage,
//     _["cumulativeRejection"] = cp,
//     _["cumulativeAlphaSpent"] = cumAlphaSpent,
//     _["cumulativeAttainedAlphaH10"] = cpH10,
//     _["cumulativeAttainedAlphaH20"] = cpH20,
//     _["numberOfEvents"] = nevents,
//     _["numberOfDropouts"] = ndropouts,
//     _["numberOfSubjects"] = nsubjects,
//     _["analysisTime"] = time,
//     _["efficacyHRLower"] = efficacyHRLower,
//     _["efficacyHRUpper"] = efficacyHRUpper,
//     _["efficacyP"] = efficacyP,
//     _["information"] = I,
//     _["HR"] = hazardRatio);
//
//   List settings = List::create(
//     _["typeAlphaSpending"] = typeAlphaSpending,
//     _["parameterAlphaSpending"] = parameterAlphaSpending,
//     _["userAlphaSpending"] = userAlphaSpending,
//     _["allocationRatioPlanned"] = allocationRatioPlanned,
//     _["accrualTime"] = accrualTime,
//     _["accrualIntensity"] = accrualIntensity,
//     _["piecewiseSurvivalTime"] = piecewiseSurvivalTime,
//     _["stratumFraction"] = stratumFraction,
//     _["lambda1"] = lambda1,
//     _["lambda2"] = lambda2,
//     _["gamma1"] = gamma1,
//     _["gamma2"] = gamma2,
//     _["typeOfComputation"] = typeOfComputation,
//     _["spendingTime"] = spendingTime);
//
//   List byTreatmentCounts = List::create(
//     _["numberOfEvents1"] = nevents1,
//     _["numberOfDropouts1"] = ndropouts1,
//     _["numberOfSubjects1"] = nsubjects1,
//     _["numberOfEvents2"] = nevents2,
//     _["numberOfDropouts2"] = ndropouts2,
//     _["numberOfSubjects2"] = nsubjects2,
//     _["expectedNumberOfEvents1"] = expectedNumberOfEvents1,
//     _["expectedNumberOfDropouts1"] = expectedNumberOfDropouts1,
//     _["expectedNumberOfSubjects1"] = expectedNumberOfSubjects1,
//     _["expectedNumberOfEvents2"] = expectedNumberOfEvents2,
//     _["expectedNumberOfDropouts2"] = expectedNumberOfDropouts2,
//     _["expectedNumberOfSubjects2"] = expectedNumberOfSubjects2);
//
//   List result = List::create(
//     _["byStageResults"] = byStageResults,
//     _["overallResults"] = overallResults,
//     _["settings"] = settings,
//     _["byTreatmentCounts"] = byTreatmentCounts);
//
//   result.attr("class") = "lrpowerequiv";
//
//   return result;
// }
//
//
// //' @title Sample Size for Equivalence in Hazard Ratio
// //' @description Obtains the sample size for equivalence in hazard ratio.
// //'
// //' @param beta The type II error.
// //' @inheritParams param_kMax
// //' @param informationRates The information rates.
// //'   Defaults to \code{(1:kMax) / kMax} if left unspecified.
// //' @inheritParams param_criticalValues
// //' @param alpha The significance level for each of the two one-sided
// //'   tests. Defaults to 0.05.
// //' @inheritParams param_typeAlphaSpending
// //' @inheritParams param_parameterAlphaSpending
// //' @inheritParams param_userAlphaSpending
// //' @param hazardRatioLower The lower equivalence limit of hazard ratio.
// //' @param hazardRatioUpper The upper equivalence limit of hazard ratio.
// //' @inheritParams param_allocationRatioPlanned
// //' @inheritParams param_accrualTime
// //' @inheritParams param_accrualIntensity
// //' @inheritParams param_piecewiseSurvivalTime
// //' @inheritParams param_stratumFraction
// //' @inheritParams param_lambda1_stratified
// //' @inheritParams param_lambda2_stratified
// //' @inheritParams param_gamma1_stratified
// //' @inheritParams param_gamma2_stratified
// //' @inheritParams param_accrualDuration
// //' @inheritParams param_followupTime
// //' @inheritParams param_fixedFollowup
// //' @inheritParams param_typeOfComputation
// //' @param interval The interval to search for the solution of
// //'   accrualDuration, followupDuration, or the proportionality constant
// //'   of accrualIntensity. Defaults to \code{c(0.001, 240)}.
// //' @param spendingTime A vector of length \code{kMax} for the error spending
// //'   time at each analysis. Defaults to missing, in which case, it is the
// //'   same as \code{informationRates}.
// //' @param rounding Whether to round up sample size.
// //'   Defaults to 1 for sample size rounding.
// //'
// //' @return An S3 class \code{lrpowerequiv} object
// //'
// //' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
// //'
// //' @seealso \code{\link{lrpowerequiv}}
// //'
// //' @examples
// //'
// //' lrsamplesizeequiv(kMax = 2, informationRates = c(0.5, 1),
// //'                   alpha = 0.05, typeAlphaSpending = "sfOF",
// //'                   hazardRatioLower = 0.71, hazardRatioUpper = 1.4,
// //'                   allocationRatioPlanned = 1, accrualTime = seq(0, 8),
// //'                   accrualIntensity = 26/9*seq(1, 9),
// //'                   piecewiseSurvivalTime = c(0, 6),
// //'                   lambda1 = c(0.0533, 0.0533),
// //'                   lambda2 = c(0.0533, 0.0533),
// //'                   gamma1 = -log(1-0.05)/12,
// //'                   gamma2 = -log(1-0.05)/12, accrualDuration = NA,
// //'                   followupTime = 18, fixedFollowup = FALSE)
// //'
// //' @export
// // [[Rcpp::export]]
// List lrsamplesizeequiv(const double beta = 0.2,
//                        const int kMax = 1,
//                        const NumericVector& informationRates = NA_REAL,
//                        const NumericVector& criticalValues = NA_REAL,
//                        const double alpha = 0.05,
//                        const std::string typeAlphaSpending = "sfOF",
//                        const double parameterAlphaSpending = NA_REAL,
//                        const NumericVector& userAlphaSpending = NA_REAL,
//                        const double hazardRatioLower = NA_REAL,
//                        const double hazardRatioUpper = NA_REAL,
//                        const double allocationRatioPlanned = 1,
//                        const NumericVector& accrualTime = 0,
//                        const NumericVector& accrualIntensity = NA_REAL,
//                        const NumericVector& piecewiseSurvivalTime = 0,
//                        const NumericVector& stratumFraction = 1,
//                        const NumericVector& lambda1 = NA_REAL,
//                        const NumericVector& lambda2 = NA_REAL,
//                        const NumericVector& gamma1 = 0,
//                        const NumericVector& gamma2 = 0,
//                        double accrualDuration = NA_REAL,
//                        double followupTime = NA_REAL,
//                        const bool fixedFollowup = 0,
//                        const std::string typeOfComputation = "direct",
//                        const NumericVector& interval =
//                          NumericVector::create(0.001, 240),
//                          const NumericVector& spendingTime = NA_REAL,
//                          const bool rounding = 1) {
//
//   NumericVector informationRates1 = clone(informationRates);
//   NumericVector criticalValues1 = clone(criticalValues);
//   NumericVector accrualIntensity1 = clone(accrualIntensity);
//   NumericVector spendingTime1 = clone(spendingTime);
//
//   std::string asf = typeAlphaSpending;
//   std::for_each(asf.begin(), asf.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//
//   double asfpar = parameterAlphaSpending;
//
//   int nstrata = static_cast<int>(stratumFraction.size());
//   int nintervals = static_cast<int>(piecewiseSurvivalTime.size());
//   int nsi = nstrata*nintervals;
//
//
//   if (beta >= 1-alpha || beta < 0.0001) {
//     stop("beta must lie in [0.0001, 1-alpha)");
//   }
//
//   if (kMax < 1) {
//     stop("kMax must be a positive integer");
//   }
//
//   if (is_false(any(is_na(informationRates)))) {
//     if (informationRates.size() != kMax) {
//       stop("Invalid length for informationRates");
//     } else if (informationRates[0] <= 0) {
//       stop("Elements of informationRates must be positive");
//     } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
//       stop("Elements of informationRates must be increasing");
//     } else if (informationRates[kMax-1] != 1) {
//       stop("informationRates must end with 1");
//     }
//   } else {
//     IntegerVector tem = seq_len(kMax);
//     informationRates1 = NumericVector(tem)/(kMax+0.0);
//   }
//
//   if (is_false(any(is_na(criticalValues)))) {
//     if (criticalValues.size() != kMax) {
//       stop("Invalid length for criticalValues");
//     }
//   }
//
//   if (!std::isnan(alpha)) {
//     if (alpha < 0.00001 || alpha >= 1) {
//       stop("alpha must lie in [0.00001, 1)");
//     }
//   }
//
//   if (is_true(any(is_na(criticalValues))) && std::isnan(alpha)) {
//     stop("alpha must be provided when criticalValues is missing");
//   }
//
//   if (is_true(any(is_na(criticalValues))) && !(asf=="of" || asf=="p" ||
//       asf=="wt" || asf=="sfof" || asf=="sfp" ||
//       asf=="sfkd" || asf=="sfhsd" || asf=="user" || asf=="none")) {
//     stop("Invalid value for typeAlphaSpending");
//   }
//
//   if ((asf=="wt" || asf=="sfkd" || asf=="sfhsd") && std::isnan(asfpar)) {
//     stop("Missing value for parameterAlphaSpending");
//   }
//
//   if (asf=="sfkd" && asfpar <= 0) {
//     stop ("parameterAlphaSpending must be positive for sfKD");
//   }
//
//   if (is_true(any(is_na(criticalValues))) && asf=="user") {
//     if (is_true(any(is_na(userAlphaSpending)))) {
//       stop("userAlphaSpending must be specified");
//     } else if (userAlphaSpending.size() < kMax) {
//       stop("Insufficient length of userAlphaSpending");
//     } else if (userAlphaSpending[0] < 0) {
//       stop("Elements of userAlphaSpending must be nonnegative");
//     } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
//       stop("Elements of userAlphaSpending must be nondecreasing");
//     } else if (userAlphaSpending[kMax-1] != alpha) {
//       stop("userAlphaSpending must end with specified alpha");
//     }
//   }
//
//   if (std::isnan(hazardRatioLower)) {
//     stop("hazardRatioLower must be provided");
//   }
//
//   if (std::isnan(hazardRatioUpper)) {
//     stop("hazardRatioUpper must be provided");
//   }
//
//   if (hazardRatioLower <= 0) {
//     stop("hazardRatioLower must be positive");
//   }
//
//   if (hazardRatioLower >= hazardRatioUpper) {
//     stop("hazardRatioLower must be less than hazardRatioUpper");
//   }
//
//   if (allocationRatioPlanned <= 0) {
//     stop("allocationRatioPlanned must be positive");
//   }
//
//   if (accrualTime[0] != 0) {
//     stop("accrualTime must start with 0");
//   }
//
//   if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
//     stop("accrualTime should be increasing");
//   }
//
//   if (is_true(any(is_na(accrualIntensity)))) {
//     stop("accrualIntensity must be provided");
//   }
//
//   if (accrualTime.size() != accrualIntensity.size()) {
//     stop("accrualTime must have the same length as accrualIntensity");
//   }
//
//   if (is_true(any(accrualIntensity < 0))) {
//     stop("accrualIntensity must be non-negative");
//   }
//
//   if (piecewiseSurvivalTime[0] != 0) {
//     stop("piecewiseSurvivalTime must start with 0");
//   }
//
//   if (nintervals > 1 && is_true(any(diff(piecewiseSurvivalTime) <= 0))) {
//     stop("piecewiseSurvivalTime should be increasing");
//   }
//
//   if (is_true(any(stratumFraction <= 0))) {
//     stop("stratumFraction must be positive");
//   }
//
//   if (sum(stratumFraction) != 1) {
//     stop("stratumFraction must sum to 1");
//   }
//
//   if (is_true(any(is_na(lambda1)))) {
//     stop("lambda1 must be provided");
//   }
//
//   if (is_true(any(is_na(lambda2)))) {
//     stop("lambda2 must be provided");
//   }
//
//   if (is_true(any(lambda1 < 0))) {
//     stop("lambda1 must be non-negative");
//   }
//
//   if (is_true(any(lambda2 < 0))) {
//     stop("lambda2 must be non-negative");
//   }
//
//   if (is_true(any(gamma1 < 0))) {
//     stop("gamma1 must be non-negative");
//   }
//
//   if (is_true(any(gamma2 < 0))) {
//     stop("gamma2 must be non-negative");
//   }
//
//   if (lambda1.size() != 1 && lambda1.size() != nintervals &&
//       lambda1.size() != nsi) {
//     stop("Invalid length for lambda1");
//   }
//
//   if (lambda2.size() != 1 && lambda2.size() != nintervals &&
//       lambda2.size() != nsi) {
//     stop("Invalid length for lambda2");
//   }
//
//   if (gamma1.size() != 1 && gamma1.size() != nintervals &&
//       gamma1.size() != nsi) {
//     stop("Invalid length for gamma1");
//   }
//
//   if (gamma2.size() != 1 && gamma2.size() != nintervals &&
//       gamma2.size() != nsi) {
//     stop("Invalid length for gamma2");
//   }
//
//
//   if (!std::isnan(accrualDuration)) {
//     if (accrualDuration <= 0) {
//       stop("accrualDuration must be positive");
//     }
//   }
//
//   if (!std::isnan(followupTime)) {
//     if (fixedFollowup && followupTime <= 0) {
//       stop("followupTime must be positive for fixed follow-up");
//     }
//
//     if (!fixedFollowup && followupTime < 0) {
//       stop("followupTime must be non-negative for variable follow-up");
//     }
//   }
//
//   if (fixedFollowup && std::isnan(followupTime)) {
//     stop("followupTime must be provided for fixed follow-up");
//   }
//
//   std::string su = typeOfComputation;
//   std::for_each(su.begin(), su.end(), [](char & c) {
//     c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
//   });
//
//   char su1 = su[0];
//
//   if (su1 != 'd' && su1 != 's') {
//     stop("typeOfComputation must be direct or Schoenfeld");
//   }
//
//   double hazardRatio = 1.0;
//   if (su1 == 's') {
//     NumericVector lambda1x = rep(lambda1, nsi/lambda1.size());
//     NumericVector lambda2x = rep(lambda2, nsi/lambda2.size());
//     NumericVector hrx = lambda1x / lambda2x;
//
//     bool proportionalHazards = 1;
//     for (int i=1; i<nsi; ++i) {
//       if (fabs(hrx[i] - hrx[0]) > 1e-8) {
//         proportionalHazards = 0;
//         break;
//       }
//     }
//
//     if (!proportionalHazards) {
//       stop("Schoenfeld method can only be used for proportional hazards");
//     } else {
//       hazardRatio = hrx[0];
//     }
//   }
//
//   if (interval.size() != 2) {
//     stop("interval must have 2 elements");
//   }
//
//   if (interval[0] < 0) {
//     stop("lower limit of interval must be positive");
//   }
//
//   if (interval[0] >= interval[1]) {
//     stop("upper limit must be greater than lower limit for interval");
//   }
//
//   if (is_false(any(is_na(spendingTime)))) {
//     if (spendingTime.size() != kMax) {
//       stop("Invalid length for spendingTime");
//     } else if (spendingTime[0] <= 0) {
//       stop("Elements of spendingTime must be positive");
//     } else if (kMax > 1 && is_true(any(diff(spendingTime) <= 0))) {
//       stop("Elements of spendingTime must be increasing");
//     } else if (spendingTime[kMax-1] != 1) {
//       stop("spendingTime must end with 1");
//     }
//   } else {
//     spendingTime1 = clone(informationRates1);
//   }
//
//
//   // obtain criticalValues
//   if (is_true(any(is_na(criticalValues)))) {
//     if (kMax > 1 && criticalValues.size() == kMax &&
//         is_false(any(is_na(head(criticalValues, kMax-1)))) &&
//         std::isnan(criticalValues[kMax-1])) { // Haybittle & Peto
//
//       auto f = [kMax, informationRates1,
//                 criticalValues, alpha](double aval)->double {
//                   NumericVector u(kMax), l(kMax, -6.0), zero(kMax);
//                   for (int i=0; i<kMax-1; ++i) {
//                     u[i] = criticalValues[i];
//                   }
//                   u[kMax-1] = aval;
//
//                   List probs = exitprobcpp(u, l, zero, informationRates1);
//                   double cpu = sum(NumericVector(probs[0]));
//                   return cpu - alpha;
//                 };
//
//       criticalValues1[kMax-1] = brent(f, -5.0, 6.0, 1.0e-6);
//     } else {
//       LogicalVector efficacyStopping1(kMax, 1);
//       criticalValues1 = getBoundcpp(kMax, informationRates1, alpha,
//                                     asf, asfpar, userAlphaSpending,
//                                     spendingTime1, efficacyStopping1);
//     }
//   }
//
//
//   std::string unknown;
//   // search for the solution according to the input
//   if (std::isnan(accrualDuration) && !std::isnan(followupTime)) {
//     unknown = "accrualDuration";
//   } else if (!std::isnan(accrualDuration) && std::isnan(followupTime)) {
//     unknown = "followupTime";
//   } else if (!std::isnan(accrualDuration) && !std::isnan(followupTime)) {
//     unknown = "accrualIntensity";
//   } else {
//     stop("accrualDuration and followupTime cannot be both missing");
//   }
//
//
//   double phi = allocationRatioPlanned/(1.0 + allocationRatioPlanned);
//   NumericVector b = criticalValues1;
//   NumericVector li(kMax, -6.0), ui(kMax, 6.0), zero(kMax);
//   double theta10 = log(hazardRatioLower), theta20 = log(hazardRatioUpper);
//
//   if (su1 == 's') {
//     List design = getDesignEquiv(
//       beta, NA_REAL, theta10, theta20, log(hazardRatio),
//       kMax, informationRates1, criticalValues1,
//       alpha, asf, asfpar, userAlphaSpending, spendingTime1);
//
//     DataFrame overallResults = DataFrame(design["overallResults"]);
//     double maxInformation = overallResults["information"];
//     double D = maxInformation/(phi*(1-phi));
//
//     auto f = [allocationRatioPlanned,
//               accrualTime, accrualIntensity,
//               piecewiseSurvivalTime, stratumFraction,
//               lambda1, lambda2, gamma1, gamma2,
//               accrualDuration, followupTime, fixedFollowup,
//               unknown, D](double aval)-> double{
//                 NumericVector accrualIntensity1 = clone(accrualIntensity);
//                 double dur1=0, dur2=0;
//
//                 if (unknown == "accrualDuration") {
//                   dur1 = aval;
//                   dur2 = followupTime;
//                 } else if (unknown == "followupTime") {
//                   dur1 = accrualDuration;
//                   dur2 = aval;
//                 } else if (unknown == "accrualIntensity") {
//                   dur1 = accrualDuration;
//                   dur2 = followupTime;
//                   accrualIntensity1 = aval*accrualIntensity;
//                 }
//
//                 // obtain the total number of events at study end
//                 NumericVector u0(1, dur1 + dur2);
//                 DataFrame lr = lrstat(
//                   u0, 1, allocationRatioPlanned,
//                   accrualTime, accrualIntensity1,
//                   piecewiseSurvivalTime, stratumFraction,
//                   lambda1, lambda2, gamma1, gamma2,
//                   dur1, dur2, fixedFollowup, 0, 0, 1);
//
//                 return sum(NumericVector(lr[2])) - D;
//               };
//
//     if (unknown == "accrualDuration") {
//       accrualDuration = brent(f, interval[0], interval[1], 1.0e-6);
//     } else if (unknown == "followupTime") {
//       followupTime = brent(f, interval[0], interval[1], 1.0e-6);
//     } else if (unknown == "accrualIntensity") {
//       double aval = brent(f, interval[0], interval[1], 1.0e-6);
//       accrualIntensity1 = aval*accrualIntensity;
//     }
//   } else {
//     auto f = [beta, kMax, informationRates1, criticalValues1,
//               theta10, theta20, allocationRatioPlanned,
//               accrualTime, accrualIntensity,
//               piecewiseSurvivalTime, stratumFraction,
//               lambda1, lambda2, gamma1, gamma2,
//               accrualDuration, followupTime, fixedFollowup,
//               unknown, li, ui](double aval)->double {
//                 NumericVector accrualIntensity1 = clone(accrualIntensity);
//                 double dur1=0, dur2=0;
//
//                 if (unknown == "accrualDuration") {
//                   dur1 = aval;
//                   dur2 = followupTime;
//                 } else if (unknown == "followupTime") {
//                   dur1 = accrualDuration;
//                   dur2 = aval;
//                 } else if (unknown == "accrualIntensity") {
//                   dur1 = accrualDuration;
//                   dur2 = followupTime;
//                   accrualIntensity1 = aval*accrualIntensity;
//                 }
//
//                 DataFrame lr;
//                 NumericVector e0(kMax), time(kMax);
//
//                 double studyDuration1 = dur1 + dur2;
//                 NumericVector u0(1, studyDuration1);
//
//                 // obtain the timing of interim analysis
//                 lr = lrstat(u0, 1, allocationRatioPlanned,
//                             accrualTime, accrualIntensity1,
//                             piecewiseSurvivalTime, stratumFraction,
//                             lambda1, lambda2, gamma1, gamma2,
//                             dur1, dur2, fixedFollowup, 0, 0, 1);
//
//                 e0 = sum(NumericVector(lr[2]))*informationRates1;
//                 if (kMax > 1) {
//                   NumericVector nevents1 = e0[Range(0, kMax-2)];
//                   NumericVector time1 = caltime(
//                     nevents1, allocationRatioPlanned,
//                     accrualTime, accrualIntensity1,
//                     piecewiseSurvivalTime, stratumFraction,
//                     lambda1, lambda2, gamma1, gamma2,
//                     dur1, dur2, fixedFollowup);
//                   time[Range(0, kMax-2)] = time1;
//                 }
//                 time[kMax-1] = studyDuration1;
//
//                 // obtain the mean and information at each stage
//                 lr = lrstat(time, 1, allocationRatioPlanned,
//                             accrualTime, accrualIntensity1,
//                             piecewiseSurvivalTime, stratumFraction,
//                             lambda1, lambda2, gamma1, gamma2,
//                             dur1, dur2, fixedFollowup, 0, 0, 3);
//
//                 NumericVector HR = NumericVector(lr[15]);
//                 NumericVector theta = log(HR);
//                 NumericVector I = 1.0/NumericVector(lr[16]);
//
//                 NumericVector b = criticalValues1;
//                 NumericVector l = b + theta10*sqrt(I);
//                 NumericVector u = -b + theta20*sqrt(I);
//
//                 List probs1 = exitprobcpp(pmax(l, li), li, theta, I);
//                 List probs2 = exitprobcpp(ui, pmin(u, ui), theta, I);
//
//                 double cpl = sum(NumericVector(probs1[0]));
//                 double cpu = sum(NumericVector(probs2[1]));
//
//                 double power;
//                 if (is_true(any(l <= u))) {
//                   power = cpl + cpu - 1;
//                 } else {
//                   List a = exitprobcpp(l, u, theta, I);
//                   double p = sum(NumericVector(a[0]) + NumericVector(a[1]));
//                   power = cpl + cpu - p;
//                 }
//
//                 return power - (1-beta);
//               };
//
//     if (unknown == "accrualDuration") {
//       accrualDuration = brent(f, interval[0], interval[1], 1.0e-6);
//     } else if (unknown == "followupTime") {
//       followupTime = brent(f, interval[0], interval[1], 1.0e-6);
//     } else if (unknown == "accrualIntensity") {
//       double aval = brent(f, interval[0], interval[1], 1.0e-6);
//       accrualIntensity1 = aval*accrualIntensity;
//     }
//   }
//
//
//   // output the results
//   List result;
//   if (rounding) {
//     NumericVector u0(1, accrualDuration + followupTime);
//     DataFrame lr = lrstat(u0, 1, allocationRatioPlanned,
//                           accrualTime, accrualIntensity1,
//                           piecewiseSurvivalTime, stratumFraction,
//                           lambda1, lambda2, gamma1, gamma2,
//                           accrualDuration, followupTime, fixedFollowup,
//                           0, 0, 1);
//
//     // round up the total number of events
//     double D0 = sum(NumericVector(lr[2]));
//     double D = std::ceil(D0 - 1.0e-12);
//
//     // adjust design parameters to obtain integer number of events
//     double n0, n, studyDuration;
//     if (!fixedFollowup) {
//       n0 = sum(NumericVector(lr[1]));
//       n = std::ceil(n0 - 1.0e-12);
//
//       if (n - n0 > 1e-6) {
//         // adjust accrual intensity or duration to obtain int # of subjects
//         if (unknown == "accrualIntensity") {
//           accrualIntensity1 = (n/n0)*accrualIntensity1;
//         } else {
//           NumericVector ns(1, n);
//           accrualDuration = getAccrualDurationFromN(ns, accrualTime,
//                                                     accrualIntensity1)[0];
//         }
//       }
//
//       // adjust follow-up time to obtain integer number of events
//       auto h = [allocationRatioPlanned,
//                 accrualTime, accrualIntensity1,
//                 piecewiseSurvivalTime, stratumFraction,
//                 lambda1, lambda2, gamma1, gamma2,
//                 accrualDuration, fixedFollowup,
//                 D](double aval)->double {
//                   NumericVector u0(1, accrualDuration + aval);
//                   DataFrame lr = lrstat(
//                     u0, 1, allocationRatioPlanned,
//                     accrualTime, accrualIntensity1,
//                     piecewiseSurvivalTime, stratumFraction,
//                     lambda1, lambda2, gamma1, gamma2,
//                     accrualDuration, aval, fixedFollowup, 0, 0, 1);
//                   return sum(NumericVector(lr[2])) - D;
//                 };
//
//       double lower = 0.0, upper = 1.1*followupTime;
//       while (h(upper) < 0) {
//         lower = upper;
//         upper = 2.0*upper;
//       }
//       followupTime = brent(h, lower, upper, 1.0e-6);
//       studyDuration = accrualDuration + followupTime;
//     } else {
//       // adjust accrual intensity or duration to obtain int number of events
//       if (unknown == "accrualIntensity") {
//         accrualIntensity1 = (D/D0)*accrualIntensity1;
//       } else {
//         auto h = [allocationRatioPlanned,
//                   accrualTime, accrualIntensity1,
//                   piecewiseSurvivalTime, stratumFraction,
//                   lambda1, lambda2, gamma1, gamma2,
//                   followupTime, fixedFollowup,
//                   D](double aval)->double {
//                     NumericVector u0(1, aval + followupTime);
//                     DataFrame lr = lrstat(
//                       u0, 1, allocationRatioPlanned,
//                       accrualTime, accrualIntensity1,
//                       piecewiseSurvivalTime, stratumFraction,
//                       lambda1, lambda2, gamma1, gamma2,
//                       aval, followupTime, fixedFollowup, 0, 0, 1);
//                     return sum(NumericVector(lr[2])) - D;
//                   };
//
//         double lower = accrualDuration, upper = 1.1*accrualDuration;
//         while (h(upper) < 0) {
//           lower = upper;
//           upper = 2.0*upper;
//         }
//         accrualDuration = brent(h, lower, upper, 1.0e-6);
//       }
//
//       NumericVector u0(1, accrualDuration);
//       n0 = accrual(u0, accrualTime, accrualIntensity1, accrualDuration)[0];
//
//       // round up the sample size
//       n = std::ceil(n0 - 1.0e-12);
//
//       if (n - n0 > 1e-6) {
//         if (unknown == "accrualIntensity") {
//           accrualIntensity1 = (n/n0)*accrualIntensity1;
//         } else {
//           NumericVector ns(1, n);
//           accrualDuration = getAccrualDurationFromN(ns, accrualTime,
//                                                     accrualIntensity1)[0];
//         }
//       }
//
//       // adjust study duration to obtain integer number of events
//       auto h = [allocationRatioPlanned,
//                 accrualTime, accrualIntensity1,
//                 piecewiseSurvivalTime, stratumFraction,
//                 lambda1, lambda2, gamma1, gamma2,
//                 accrualDuration, followupTime, fixedFollowup,
//                 D](double aval)->double {
//                   NumericVector u0(1, accrualDuration + aval);
//                   DataFrame lr = lrstat(
//                     u0, 1, allocationRatioPlanned,
//                     accrualTime, accrualIntensity1,
//                     piecewiseSurvivalTime, stratumFraction,
//                     lambda1, lambda2, gamma1, gamma2,
//                     accrualDuration, followupTime, fixedFollowup, 0, 0, 1);
//                   return sum(NumericVector(lr[2])) - D;
//                 };
//
//       double aval = brent(h, 0.0, followupTime, 1.0e-6);
//       studyDuration = accrualDuration + aval;
//     }
//
//     // update information rates to calculate new boundaries
//     NumericVector nevents = floor(D*informationRates1 + 0.5);
//     informationRates1 = nevents/nevents[kMax-1];
//
//     // recalculate boundaries
//     result = lrpowerequiv(
//       kMax, informationRates1, criticalValues,
//       alpha, typeAlphaSpending, parameterAlphaSpending,
//       userAlphaSpending, hazardRatioLower, hazardRatioUpper,
//       allocationRatioPlanned, accrualTime, accrualIntensity1,
//       piecewiseSurvivalTime, stratumFraction,
//       lambda1, lambda2, gamma1, gamma2,
//       accrualDuration, followupTime, fixedFollowup,
//       typeOfComputation, spendingTime, studyDuration);
//   } else {
//     double studyDuration = accrualDuration + followupTime;
//
//     result = lrpowerequiv(
//       kMax, informationRates1, criticalValues1,
//       alpha, typeAlphaSpending, parameterAlphaSpending,
//       userAlphaSpending, hazardRatioLower, hazardRatioUpper,
//       allocationRatioPlanned, accrualTime, accrualIntensity1,
//       piecewiseSurvivalTime, stratumFraction,
//       lambda1, lambda2, gamma1, gamma2,
//       accrualDuration, followupTime, fixedFollowup,
//       typeOfComputation, spendingTime, studyDuration);
//   }
//
//   return result;
// }
