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

std::vector<double> nb_make_breaks(
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& accrualTime,
    const double accrualDuration,
    const double maxFollowupTime,
    const double time) {

  double upper = std::min(time, maxFollowupTime);

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


// Mean exposure E(t), where t = min(time - W, C, maxFollowupTime) and
// W is the enrollment time, C is dropout time, and the expectation is taken
// with respect to the distribution of W and C: E(t) = integrate(S(t) dt)
double nb_ex_integrand(
    const double x, // analysis time
    const double time, // calendar time of analysis
    const double phi,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const double kappa,
    const double lambda,
    const std::vector<double>& zero,
    const std::vector<double>& gamma,
    const double accrualDuration) {

  // p(x) = P(at risk at time x since randomization) under dropout hazards
  double p = patrisk1(x, piecewiseSurvivalTime, zero, gamma);

  // N(time - x): number enrolled by calendar time (time-x)
  double N = accrual1(time - x, accrualTime, accrualIntensity, accrualDuration);

  return phi * N * p;
}


// Expected (Fisher) information for log rate ratio (beta) under negative
// binomial (NB) model with piecewise exponential dropouts.
// E(-d^2 logL/d beta^2) = n E(lambda * t / (1 + kappa*lambda*t))
// where t is the exposure time observed at calendar time x
// and the expectation is taken with respect to the distribution of t,
// which is determined by the dropout hazards and the accrual pattern:
// E(lambda * t/(1 + kappa * lambda *t)) =
//   integrate( lambda /(1 + kappa * lambda * t)^2 * H(time - t) G(t) dt).
// where n H(time - t) is the number of subjects enrolled by calendar
// time (time-t), and G(t) is the probability of being at risk of dropout
// at time t since randomization under the piecewise exponential dropout model.
double nb_info_integrand(
    const double x, // analysis time
    const double time, // calendar time of analysis
    const double phi,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const double kappa,
    const double lambda,
    const std::vector<double>& zero,
    const std::vector<double>& gamma,
    const double accrualDuration) {

  // p(x) = P(at risk at time x since randomization) under dropout hazards
  double p = patrisk1(x, piecewiseSurvivalTime, zero, gamma);

  // N(time - x): number enrolled by calendar time (time-x)
  double N = accrual1(time - x, accrualTime, accrualIntensity, accrualDuration);

  // key ingredient to integrate: u(x) = lambda / (1 + kappa * lambda * x)^2
  double u = lambda / sq(1.0 + kappa * lambda * x);

  return u * phi * N * p;
}


ListCpp nbstat1cpp(
    const double time,
    const double rateRatioH0,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& kappa1,
    const std::vector<double>& kappa2,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<std::vector<double>>& gamma1,
    const std::vector<std::vector<double>>& gamma2,
    const double accrualDuration,
    const double followupTime,
    const bool fixedFollowup,
    const bool nullVariance) {

  const std::size_t nstrata = stratumFraction.size();
  const std::size_t nintv = piecewiseSurvivalTime.size();
  const std::vector<double> zero(nintv, 0.0);

  // phi = P(randomized to group 1)
  const double phi = allocationRatioPlanned / (1.0 + allocationRatioPlanned);

  // max follow-up time for first enrolled subject
  const double maxFollowupTime = fixedFollowup ? followupTime :
    (accrualDuration + followupTime);

  // total enrolled by calendar time 'time'
  const double a = accrual1(time, accrualTime, accrualIntensity, accrualDuration);

  // enrolled by (time - maxFollowupTime)
  const double a2 = accrual1(time - maxFollowupTime, accrualTime,
                             accrualIntensity, accrualDuration);

  // ---- outputs ----
  std::vector<int> stratum(nstrata);
  std::vector<double> calTime(nstrata, time);
  std::vector<double> nsubjects(nstrata);
  std::vector<double> nevents(nstrata), nevents1(nstrata), nevents2(nstrata);
  std::vector<double> ndropouts(nstrata), ndropouts1(nstrata), ndropouts2(nstrata);
  std::vector<double> nfmax(nstrata), nfmax1(nstrata), nfmax2(nstrata);
  std::vector<double> exposure(nstrata), exposure1(nstrata), exposure2(nstrata);
  std::vector<double> rateRatio(nstrata);
  std::vector<double> vlogRate1(nstrata), vlogRate2(nstrata), vlogRR(nstrata);
  std::vector<double> lambda1H0(nstrata), lambda2H0(nstrata);
  std::vector<double> vlogRate1H0(nstrata), vlogRate2H0(nstrata), vlogRRH0(nstrata);

  auto breaks = nb_make_breaks(piecewiseSurvivalTime, accrualTime,
                               accrualDuration, maxFollowupTime, time);
  double tol = 1.0e-6;

  for (std::size_t h = 0; h < nstrata; ++h) {
    stratum[h] = static_cast<int>(h + 1);
    const double frac = stratumFraction[h];

    // subset per-stratum parameters
    const double k1 = kappa1[h];
    const double k2 = kappa2[h];
    const double lam1 = lambda1[h];
    const double lam2 = lambda2[h];
    const std::vector<double>& gam1 = gamma1[h];
    const std::vector<double>& gam2 = gamma2[h];

    // scale accrualIntensity by stratum fraction
    std::vector<double> accrualIntensity_frac = accrualIntensity;
    for (double& v : accrualIntensity_frac) v *= frac;

    // number of enrolled subjects
    nsubjects[h] = frac*a;

    // number of dropouts (won't be censored by recurrent events)
    bool hasDropout = false;
    for (size_t i = 0; i < nintv; ++i) {
      if (gam1[i] > 0.0 || gam2[i] > 0.0) {
        hasDropout = true;
        break;
      }
    }
    if (hasDropout) {
      auto nd_row = nevent21cpp(
        time, allocationRatioPlanned, accrualTime, accrualIntensity_frac,
        piecewiseSurvivalTime, gam1, gam2, zero, zero,
        accrualDuration, followupTime, maxFollowupTime);

      ndropouts1[h] = nd_row.first;
      ndropouts2[h] = nd_row.second;
      ndropouts[h] = ndropouts1[h] + ndropouts2[h];
    }

    // number of subjects reaching maximum follow-up
    double ncom = frac * a2;
    double p1 = patrisk1(maxFollowupTime, piecewiseSurvivalTime, zero, gam1);
    double p2 = patrisk1(maxFollowupTime, piecewiseSurvivalTime, zero, gam2);
    nfmax1[h] = phi * ncom * p1;
    nfmax2[h] = (1.0 - phi) * ncom * p2;
    nfmax[h] = nfmax1[h] + nfmax2[h];

    // number of events and information for log rate ratio
    auto f1_ex = [&](double x)->double {
      return nb_ex_integrand(
        x, time, phi, accrualTime, accrualIntensity_frac,
        piecewiseSurvivalTime, k1, lam1, zero, gam1,
        accrualDuration);
    };
    auto f2_ex = [&](double x)->double {
      return nb_ex_integrand(
        x, time, 1.0 - phi, accrualTime, accrualIntensity_frac,
        piecewiseSurvivalTime, k2, lam2, zero, gam2,
        accrualDuration);
    };

    exposure1[h] = integrate3(f1_ex, breaks, tol);
    exposure2[h] = integrate3(f2_ex, breaks, tol);
    exposure[h] = exposure1[h] + exposure2[h];
    nevents1[h] = lam1 * exposure1[h];
    nevents2[h] = lam2 * exposure2[h];
    nevents[h] = nevents1[h] + nevents2[h];

    auto f1_info = [&](double x)->double {
      return nb_info_integrand(
        x, time, phi, accrualTime, accrualIntensity_frac,
        piecewiseSurvivalTime, k1, lam1, zero, gam1,
        accrualDuration);
    };
    auto f2_info = [&](double x)->double {
      return nb_info_integrand(
        x, time, 1.0 - phi, accrualTime, accrualIntensity_frac,
        piecewiseSurvivalTime, k2, lam2, zero, gam2,
        accrualDuration);
    };

    double information1 = integrate3(f1_info, breaks, tol);
    double information2 = integrate3(f2_info, breaks, tol);
    rateRatio[h] = lam1 / lam2;
    vlogRate1[h] = 1.0 / information1;
    vlogRate2[h] = 1.0 / information2;
    vlogRR[h] = vlogRate1[h] + vlogRate2[h];

    if (nullVariance) {
      if (rateRatioH0 == 1.0) {
        lambda2H0[h] = phi * lam1 + (1.0 - phi) * lam2;
        lambda1H0[h] = lambda2H0[h];
      } else {
        // coefficients for the quadratic equation to solve for lambda2H0[h]
        // under null hypothesis lambda1H0[h] = rateRatioH0 * lambda2H0[h]
        // and the expected information for log rate ratio is 1/vlogRR[h]
        double t1 = exposure1[h] / (frac * a * phi);
        double t2 = exposure2[h] / (frac * a * (1.0 - phi));
        double a = (phi * k2 + (1.0 - phi) * k1) * rateRatioH0 * t1 * t2;
        double b = -(phi * t1 * (k2 * lam1 * t2 - rateRatioH0) +
                     (1.0 - phi) * t2 * (k1 * lam2 * rateRatioH0 * t1 - 1.0));
        double c = -(phi * lam1 * t1 + (1.0 - phi) * lam2 * t2);
        double init;
        if (k1 == 0 && k2 == 0) {
          init = -c / b;
        } else {
          init = (-b + std::sqrt(b * b - 4.0 * a * c))/(2.0 * a);
        }

        auto f = [&](double aval)->double {
          auto f1_info_H0 = [&](double x)->double {
            return nb_info_integrand(
              x, time, phi, accrualTime, accrualIntensity_frac,
              piecewiseSurvivalTime, k1, aval * rateRatioH0, zero, gam1,
              accrualDuration);
          };
          auto f2_info_H0 = [&](double x)->double {
            return nb_info_integrand(
              x, time, 1.0 - phi, accrualTime, accrualIntensity_frac,
              piecewiseSurvivalTime, k2, aval, zero, gam2,
              accrualDuration);
          };

          double a1 = integrate3(f1_info_H0, breaks, tol);
          double a2 = integrate3(f2_info_H0, breaks, tol);

          return phi * (lam1 / (aval * rateRatioH0) - 1.0) * a1 +
            (1.0 - phi)*(lam2 / aval - 1.0) * a2;
        };

        lambda2H0[h] = brent(f, 0.5 * init, 2.0 * init, tol);
        lambda1H0[h] = rateRatioH0 * lambda2H0[h];
      }

      auto f1_info_H0 = [&](double x)->double {
        return nb_info_integrand(
          x, time, phi, accrualTime, accrualIntensity_frac,
          piecewiseSurvivalTime, k1, lambda1H0[h], zero, gam1,
          accrualDuration);
      };
      auto f2_info_H0 = [&](double x)->double {
        return nb_info_integrand(
          x, time, 1.0 - phi, accrualTime, accrualIntensity_frac,
          piecewiseSurvivalTime, k2, lambda2H0[h], zero, gam2,
          accrualDuration);
      };

      double information1H0 = integrate3(f1_info_H0, breaks, tol);
      double information2H0 = integrate3(f2_info_H0, breaks, tol);
      vlogRate1H0[h] = 1.0 / information1H0;
      vlogRate2H0[h] = 1.0 / information2H0;
      vlogRRH0[h] = vlogRate1H0[h] + vlogRate2H0[h];
    }
  }

  DataFrameCpp resultsUnderH1;
  resultsUnderH1.push_back(stratum, "stratum");
  resultsUnderH1.push_back(calTime, "time");
  resultsUnderH1.push_back(nsubjects, "subjects");
  resultsUnderH1.push_back(nevents, "nevents");
  resultsUnderH1.push_back(nevents1, "nevents1");
  resultsUnderH1.push_back(nevents2, "nevents2");
  resultsUnderH1.push_back(ndropouts, "ndropouts");
  resultsUnderH1.push_back(ndropouts1, "ndropouts1");
  resultsUnderH1.push_back(ndropouts2, "ndropouts2");
  resultsUnderH1.push_back(nfmax, "nfmax");
  resultsUnderH1.push_back(nfmax1, "nfmax1");
  resultsUnderH1.push_back(nfmax2, "nfmax2");
  resultsUnderH1.push_back(exposure, "exposure");
  resultsUnderH1.push_back(exposure1, "exposure1");
  resultsUnderH1.push_back(exposure2, "exposure2");
  resultsUnderH1.push_back(rateRatio, "rateRatio");
  resultsUnderH1.push_back(vlogRate1, "vlogRate1");
  resultsUnderH1.push_back(vlogRate2, "vlogRate2");
  resultsUnderH1.push_back(vlogRR, "vlogRR");

  DataFrameCpp resultsUnderH0;
  if (!nullVariance) {
    resultsUnderH0.push_back(stratum, "stratum");
    resultsUnderH0.push_back(calTime, "time");
    resultsUnderH0.push_back(rateRatioH0, "rateRatioH0");
    resultsUnderH0.push_back(lambda1, "lambda1");
    resultsUnderH0.push_back(lambda2, "lambda2");
    resultsUnderH0.push_back(rateRatio, "rateRatio");
  } else {
    resultsUnderH0.push_back(stratum, "stratum");
    resultsUnderH0.push_back(calTime, "time");
    resultsUnderH0.push_back(rateRatioH0, "rateRatioH0");
    resultsUnderH0.push_back(lambda1H0, "lambda1H0");
    resultsUnderH0.push_back(lambda2H0, "lambda2H0");
    resultsUnderH0.push_back(vlogRate1H0, "vlogRate1H0");
    resultsUnderH0.push_back(vlogRate2H0, "vlogRate2H0");
    resultsUnderH0.push_back(vlogRRH0, "vlogRRH0");
    resultsUnderH0.push_back(lambda1, "lambda1");
    resultsUnderH0.push_back(lambda2, "lambda2");
    resultsUnderH0.push_back(rateRatio, "rateRatio");
  }

  ListCpp result;
  result.push_back(resultsUnderH1, "resultsUnderH1");
  result.push_back(resultsUnderH0, "resultsUnderH0");
  return result;
}


ListCpp nbstatcpp(
    const std::vector<double>& time,
    const double rateRatioH0,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& kappa1,
    const std::vector<double>& kappa2,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double followupTime,
    const bool fixedFollowup,
    const bool nullVariance) {

  if (!none_na(time)) throw std::invalid_argument("time must be provided");
  for (double t : time) {
    if (t < 0.0) throw std::invalid_argument("time must be non-negative");
  }
  if (rateRatioH0 <= 0.0)
    throw std::invalid_argument("rateRatioH0 must be positive");
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
  if (!none_na(kappa1)) throw std::invalid_argument("kappa1 must be provided");
  if (!none_na(kappa2)) throw std::invalid_argument("kappa2 must be provided");
  if (!none_na(lambda1)) throw std::invalid_argument("lambda1 must be provided");
  if (!none_na(lambda2)) throw std::invalid_argument("lambda2 must be provided");
  for (double v : kappa1) {
    if (v < 0.0) throw std::invalid_argument("kappa1 must be non-negative");
  }
  for (double v : kappa2) {
    if (v < 0.0) throw std::invalid_argument("kappa2 must be non-negative");
  }
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

  // Expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();

  auto kappa1x = expand1(kappa1, nstrata, "kappa1");
  auto kappa2x = expand1(kappa2, nstrata, "kappa2");
  auto lambda1x = expand1(lambda1, nstrata, "lambda1");
  auto lambda2x = expand1(lambda2, nstrata, "lambda2");
  auto gamma1x = expand_stratified(gamma1, nstrata, nintv, "gamma1");
  auto gamma2x = expand_stratified(gamma2, nstrata, nintv, "gamma2");

  // Prepare outputs (vectors length = times.size())
  size_t k = time.size();
  std::vector<double> subjects(k), nevents(k), nevents1(k), nevents2(k);
  std::vector<double> ndropouts(k), ndropouts1(k), ndropouts2(k);
  std::vector<double> nfmax(k), nfmax1(k), nfmax2(k);
  std::vector<double> exposure(k), exposure1(k), exposure2(k);
  std::vector<double> rateRatio(k);
  std::vector<double> vlogRate1(k), vlogRate2(k), vlogRR(k);
  std::vector<double> information(k), zlogRR(k);
  std::vector<double> lam1(k), lam2(k), varianceRatio(k);
  std::vector<double> lam1H0(k), lam2H0(k);
  std::vector<double> vlogRate1H0(k), vlogRate2H0(k), vlogRRH0(k);
  std::vector<double> informationH0(k), zlogRRH0(k);

  // iterate across requested times and call nbstat1cpp for each
  for (size_t i = 0; i < k; ++i) {

    ListCpp out = nbstat1cpp(
      time[i], rateRatioH0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      kappa1x, kappa2x, lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup, nullVariance);

    auto dfH1 = out.get<DataFrameCpp>("resultsUnderH1");
    auto dfH0 = out.get<DataFrameCpp>("resultsUnderH0");

    subjects[i] = extract_sum(dfH1, "subjects");
    nevents[i] = extract_sum(dfH1, "nevents");
    nevents1[i] = extract_sum(dfH1, "nevents1");
    nevents2[i] = extract_sum(dfH1, "nevents2");
    ndropouts[i] = extract_sum(dfH1, "ndropouts");
    ndropouts1[i] = extract_sum(dfH1, "ndropouts1");
    ndropouts2[i] = extract_sum(dfH1, "ndropouts2");
    nfmax[i] = extract_sum(dfH1, "nfmax");
    nfmax1[i] = extract_sum(dfH1, "nfmax1");
    nfmax2[i] = extract_sum(dfH1, "nfmax2");
    exposure[i] = extract_sum(dfH1, "exposure");
    exposure1[i] = extract_sum(dfH1, "exposure1");
    exposure2[i] = extract_sum(dfH1, "exposure2");

    // weighted geometric mean for rateRatio: exp(sum(w * log(rateRatio_per_stratum)))
    std::vector<double>& rateRatio_v = dfH1.get<double>("rateRatio");
    std::vector<double>& vlogRate1_v = dfH1.get<double>("vlogRate1");
    std::vector<double>& vlogRate2_v = dfH1.get<double>("vlogRate2");
    std::vector<double>& vlogRR_v = dfH1.get<double>("vlogRR");

    double logrr = 0.0, vv1 = 0.0, vv2 = 0.0, vvr = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      double w = stratumFraction[h];
      logrr += w * std::log(rateRatio_v[h]);
      vv1 += w * w * vlogRate1_v[h];
      vv2 += w * w * vlogRate2_v[h];
      vvr += w * w * vlogRR_v[h];
    }
    rateRatio[i] = std::exp(logrr);
    vlogRate1[i] = vv1;
    vlogRate2[i] = vv2;
    vlogRR[i] = vvr;
    information[i] = 1.0 / vvr;
    zlogRR[i] = (logrr - std::log(rateRatioH0)) * std::sqrt(information[i]);

    if (!nullVariance) {
      // dfH0 expected to contain "lambda1" and "lambda2" (per-stratum)
      auto lam1_v = dfH0.get<double>("lambda1");
      auto lam2_v = dfH0.get<double>("lambda2");
      double wg_log_l1 = 0.0, wg_log_l2 = 0.0;
      for (size_t h = 0; h < nstrata; ++h) {
        wg_log_l1 += stratumFraction[h] * std::log(lam1_v[h]);
        wg_log_l2 += stratumFraction[h] * std::log(lam2_v[h]);
      }
      lam1[i] = std::exp(wg_log_l1);
      lam2[i] = std::exp(wg_log_l2);
      varianceRatio[i] = 1.0;
    } else {
      // dfH0 contains H0-specific quantities
      auto lambda1H0_v = dfH0.get<double>("lambda1H0");
      auto lambda2H0_v = dfH0.get<double>("lambda2H0");
      auto vlogRate1H0_v = dfH0.get<double>("vlogRate1H0");
      auto vlogRate2H0_v = dfH0.get<double>("vlogRate2H0");
      auto vlogRRH0_v = dfH0.get<double>("vlogRRH0");
      auto lambda1_v = dfH0.get<double>("lambda1");
      auto lambda2_v = dfH0.get<double>("lambda2");

      double wg_log_l1H0 = 0.0, wg_log_l2H0 = 0.0, wg_log_l1 = 0.0, wg_log_l2 = 0.0;
      double vv1H0 = 0.0, vv2H0 = 0.0, vvrH0 = 0.0;
      for (size_t h = 0; h < nstrata; ++h) {
        double w = stratumFraction[h];
        wg_log_l1H0 += w * std::log(lambda1H0_v[h]);
        wg_log_l2H0 += w * std::log(lambda2H0_v[h]);
        wg_log_l1 += w * std::log(lambda1_v[h]);
        wg_log_l2 += w * std::log(lambda2_v[h]);
        vv1H0 += w * w * vlogRate1H0_v[h];
        vv2H0 += w * w * vlogRate2H0_v[h];
        vvrH0 += w * w * vlogRRH0_v[h];
      }
      lam1H0[i] = std::exp(wg_log_l1H0);
      lam2H0[i] = std::exp(wg_log_l2H0);
      vlogRate1H0[i] = vv1H0;
      vlogRate2H0[i] = vv2H0;
      vlogRRH0[i] = vvrH0;
      informationH0[i] = 1.0 / vvrH0;
      zlogRRH0[i] = (logrr - std::log(rateRatioH0)) * std::sqrt(informationH0[i]);

      lam1[i] = std::exp(wg_log_l1);
      lam2[i] = std::exp(wg_log_l2);
      varianceRatio[i] = vlogRRH0[i] / vlogRR[i];
    }
  }

  // Build DataFrameCpp resultsUnderH1
  DataFrameCpp resultsUnderH1;
  resultsUnderH1.push_back(time, "time");
  resultsUnderH1.push_back(subjects, "subjects");
  resultsUnderH1.push_back(nevents, "nevents");
  resultsUnderH1.push_back(nevents1, "nevents1");
  resultsUnderH1.push_back(nevents2, "nevents2");
  resultsUnderH1.push_back(ndropouts, "ndropouts");
  resultsUnderH1.push_back(ndropouts1, "ndropouts1");
  resultsUnderH1.push_back(ndropouts2, "ndropouts2");
  resultsUnderH1.push_back(nfmax, "nfmax");
  resultsUnderH1.push_back(nfmax1, "nfmax1");
  resultsUnderH1.push_back(nfmax2, "nfmax2");
  resultsUnderH1.push_back(exposure, "exposure");
  resultsUnderH1.push_back(exposure1, "exposure1");
  resultsUnderH1.push_back(exposure2, "exposure2");
  resultsUnderH1.push_back(rateRatio, "rateRatio");
  resultsUnderH1.push_back(vlogRate1, "vlogRate1");
  resultsUnderH1.push_back(vlogRate2, "vlogRate2");
  resultsUnderH1.push_back(vlogRR, "vlogRR");
  resultsUnderH1.push_back(information, "information");
  resultsUnderH1.push_back(zlogRR, "zlogRR");

  // Build resultsUnderH0
  DataFrameCpp resultsUnderH0;
  if (!nullVariance) {
    resultsUnderH0.push_back(time, "time");
    resultsUnderH0.push_back(rateRatioH0, "rateRatioH0");
    resultsUnderH0.push_back(varianceRatio, "varianceRatio");
    resultsUnderH0.push_back(lam1, "lambda1");
    resultsUnderH0.push_back(lam2, "lambda2");
    resultsUnderH0.push_back(rateRatio, "rateRatio");
  } else {
    resultsUnderH0.push_back(time, "time");
    resultsUnderH0.push_back(lam1H0, "lambda1H0");
    resultsUnderH0.push_back(lam2H0, "lambda2H0");
    resultsUnderH0.push_back(rateRatioH0, "rateRatioH0");
    resultsUnderH0.push_back(vlogRate1H0, "vlogRate1H0");
    resultsUnderH0.push_back(vlogRate2H0, "vlogRate2H0");
    resultsUnderH0.push_back(vlogRRH0, "vlogRRH0");
    resultsUnderH0.push_back(informationH0, "informationH0");
    resultsUnderH0.push_back(zlogRRH0, "zlogRRH0");
    resultsUnderH0.push_back(varianceRatio, "varianceRatio");
    resultsUnderH0.push_back(lam1, "lambda1");
    resultsUnderH0.push_back(lam2, "lambda2");
    resultsUnderH0.push_back(rateRatio, "rateRatio");
  }

  ListCpp out;
  out.push_back(resultsUnderH1, "resultsUnderH1");
  out.push_back(resultsUnderH0, "resultsUnderH0");
  return out;
}


//' @title Negative Binomial Rate Ratio
//' @description Obtains the number of subjects accrued, number of events,
//' number of dropouts, number of subjects reaching the maximum
//' follow-up, total exposure, and variance for log rate in each group,
//' rate ratio, variance, and Wald test statistic of
//' log rate ratio at given calendar times.
//'
//' @param time A vector of calendar times for data cut.
//' @param rateRatioH0 Rate ratio under the null hypothesis.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param kappa1 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the active treatment
//'   group by stratum.
//' @param kappa2 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the control group by
//'   stratum.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group by stratum.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group by stratum.
//' @inheritParams param_gamma1_stratified
//' @inheritParams param_gamma2_stratified
//' @inheritParams param_accrualDuration
//' @inheritParams param_followupTime
//' @inheritParams param_fixedFollowup
//' @param nullVariance Whether to calculate the variance for log rate ratio
//'   under the null hypothesis.
//'
//' @details
//' The probability mass function for a negative binomial distribution with
//' dispersion parameter \eqn{\kappa_i} and rate parameter \eqn{\lambda_i}
//' is given by
//' \deqn{P(Y_{ij} = y) = \frac{\Gamma(y+1/\kappa_i)}{\Gamma(1/\kappa_i) y!}
//' \left(\frac{1}{1 + \kappa_i \lambda_i t_{ij}}\right)^{1/\kappa_i}
//' \left(\frac{\kappa_i \lambda_i t_{ij}}
//' {1 + \kappa_i \lambda_i t_{ij}}\right)^{y},}
//' where \eqn{Y_{ij}} is the event count for subject \eqn{j} in
//' treatment group \eqn{i}, and \eqn{t_{ij}} is the exposure time for
//' the subject. If \eqn{\kappa_i=0}, the negative binomial distribution
//' reduces to the Poisson distribution.
//'
//' For treatment group \eqn{i}, let \eqn{\beta_i = \log(\lambda_i)}.
//' The log-likelihood for \eqn{\{(\kappa_i, \beta_i):i=1,2\}}
//' can be written as
//' \deqn{l = \sum_{i=1}^{2}\sum_{j=1}^{n_{i}}
//' \{\log \Gamma(y_{ij} + 1/\kappa_i) - \log \Gamma(1/\kappa_i) + y_{ij}
//' (\log(\kappa_i) + \beta_i) - (y_{ij} + 1/\kappa_i)
//' \log(1+ \kappa_i \exp(\beta_i) t_{ij})\}.}
//' It follows that
//' \deqn{\frac{\partial l}{\partial \beta_i} = \sum_{j=1}^{n_i}
//' \left\{y_{ij} - (y_{ij} + 1/\kappa_i)
//' \frac{\kappa_i \exp(\beta_i) t_{ij}}
//' {1 + \kappa_i \exp(\beta_i)t_{ij}}\right\},}
//' and
//' \deqn{-\frac{\partial^2 l}{\partial \beta_i^2} =
//' \sum_{j=1}^{n_i} (y_{ij} + 1/\kappa_i) \frac{\kappa_i \lambda_i t_{ij}}
//' {(1 + \kappa_i \lambda_i t_{ij})^2}.}
//' The Fisher information for \eqn{\beta_i} is
//' \deqn{E\left(-\frac{\partial^2 l}{\partial \beta_i^2}\right)
//' = n_i E\left(\frac{\lambda_i t_{ij}}
//' {1 + \kappa_i \lambda_i t_{ij}}\right).}
//' In addition, we can show that
//' \deqn{E\left(-\frac{\partial^2 l}
//' {\partial \beta_i \partial \kappa_i}\right) = 0.}
//' Therefore, the variance of \eqn{\hat{\beta}_i} is
//' \deqn{Var(\hat{\beta}_i) = \frac{1}{n_i} \left\{
//' E\left(\frac{\lambda_i t_{ij}}{1 + \kappa_i \lambda_i t_{ij}}\right)
//' \right\}^{-1}.}
//'
//' To evaluate the integral, we need to obtain the distribution of the
//' exposure time,
//' \deqn{t_{ij} = \min(\tau - W_{ij}, C_{ij}, T_{fmax}),}
//' where \eqn{\tau} denotes the calendar time since trial start,
//' \eqn{W_{ij}} denotes the enrollment time for subject \eqn{j}
//' in treatment group \eqn{i}, \eqn{C_{ij}} denotes the time to dropout
//' after enrollment for subject \eqn{j} in treatment group \eqn{i}, and
//' \eqn{T_{fmax}} denotes the maximum follow-up time for
//' all subjects. Therefore,
//' \deqn{P(t_{ij} \geq t) = P(W_{ij} \leq \tau - t)P(C_{ij} \geq t)
//' I(t\leq T_{fmax}).}
//' Let \eqn{H} denote the distribution function of the enrollment time,
//' and \eqn{G_i} denote the survival function of the dropout time for
//' treatment group \eqn{i}. By the change of variables, we have
//' \deqn{E\left(\frac{\lambda_i t_{ij}}{1 + \kappa_i \lambda_i t_{ij}}
//' \right) = \int_{0}^{\tau \wedge T_{fmax}}
//' \frac{\lambda_i}{(1 + \kappa_i \lambda_i t)^2} H(\tau - t) G_i(t) dt.}
//' A numerical integration algorithm for a univariate function can be
//' used to evaluate the above integral.
//'
//' For the restricted maximum likelihood (reml) estimate of
//' \eqn{(\beta_1,\beta_2)} subject to the
//' constraint that \eqn{\beta_1 - \beta_2 = \Delta}, we express the
//' log-likelihood in terms of \eqn{(\beta_2,\Delta,\kappa_1,\kappa_2)},
//' and takes the derivative of the log-likelihood function with respect
//' to \eqn{\beta_2}. The resulting score equation has asymptotic limit
//' \deqn{E\left(\frac{\partial l}{\partial \beta_2}\right) = s_1 + s_2,}
//' where
//' \deqn{s_1 = n r E\left\{\lambda_1 t_{1j} - \left(\lambda_1t_{1j}
//' + \frac{1}{\kappa_1}\right) \frac{\kappa_1 e^{\tilde{\beta}_2 +
//' \Delta}t_{1j}}{1 + \kappa_1 e^{\tilde{\beta}_2 +\Delta}t_{1j}}\right\},}
//' and
//' \deqn{s_2 = n (1-r) E\left\{\lambda_2 t_{2j} -
//' \left(\lambda_2 t_{2j} + \frac{1}{\kappa_2}\right)
//' \frac{\kappa_2 e^{\tilde{\beta}_2} t_{2j}}
//' {1 + \kappa_2 e^{\tilde{\beta}_2}t_{2j}}\right\}.}
//' Here \eqn{r} is the randomization probability for the active
//' treatment group. The asymptotic limit of the reml of \eqn{\beta_2}
//' is the solution \eqn{\tilde{\beta}_2} to
//' \eqn{E\left(\frac{\partial l}{\partial \beta_2}\right) = 0.}
//'
//' @return A list with two components:
//'
//' * \code{resultsUnderH1}: A data frame containing the following variables:
//'
//'     - \code{time}: The analysis time since trial start.
//'
//'     - \code{subjects}: The number of enrolled subjects.
//'
//'     - \code{nevents}: The total number of events.
//'
//'     - \code{nevents1}: The number of events in the active treatment
//'       group.
//'
//'     - \code{nevents2}: The number of events in the control group.
//'
//'     - \code{ndropouts}: The total number of dropouts.
//'
//'     - \code{ndropouts1}: The number of dropouts in the active treatment
//'       group.
//'
//'     - \code{ndropouts2}: The number of dropouts in the control group.
//'
//'     - \code{nfmax}: The total number of subjects reaching maximum
//'       follow-up.
//'
//'     - \code{nfmax1}: The number of subjects reaching maximum follow-up
//'       in the active treatment group.
//'
//'     - \code{nfmax2}: The number of subjects reaching maximum follow-up
//'       in the control group.
//'
//'     - \code{exposure}: The total exposure time.
//'
//'     - \code{exposure1}: The exposure time for the active treatment group.
//'
//'     - \code{exposure2}: The exposure time for the control group.
//'
//'     - \code{rateRatio}: The rate ratio of the active treatment group
//'       versus the control group.
//'
//'     - \code{vlogRate1}: The variance for the log rate
//'       parameter for the active treatment group.
//'
//'     - \code{vlogRate2}: The variance for the log rate
//'       parameter for the control group.
//'
//'     - \code{vlogRR}: The variance of log rate ratio.
//'
//'     - \code{information}: The information of log rate ratio.
//'
//'     - \code{zlogRR}: The Z-statistic for log rate ratio.
//'
//' * \code{resultsUnderH0} when \code{nullVariance = TRUE}: A data frame
//'   with the following variables:
//'
//'     - \code{time}: The analysis time since trial start.
//'
//'     - \code{lambda1H0}: The restricted maximum likelihood estimate
//'       of the event rate for the active treatment group.
//'
//'     - \code{lambda2H0}: The restricted maximum likelihood estimate
//'       of the event rate for the control group.
//'
//'     - \code{rateRatioH0}: The rate ratio under H0.
//'
//'     - \code{vlogRate1H0}: The variance for the log rate
//'       parameter for the active treatment group under H0.
//'
//'     - \code{vlogRate2H0}: The variance for the log rate
//'       parameter for the control group under H0.
//'
//'     - \code{vlogRRH0}: The variance of log rate ratio under H0.
//'
//'     - \code{informationH0}: The information of log rate ratio under H0.
//'
//'     - \code{zlogRRH0}: The Z-statistic for log rate ratio with variance
//'       evaluated under H0.
//'
//'     - \code{varianceRatio}: The ratio of the variance under H0 versus
//'       the variance under H1.
//'
//'     - \code{lambda1}: The true event rate for the active treatment group.
//'
//'     - \code{lambda2}: The true event rate for the control group.
//'
//'     - \code{rateRatio}: The true rate ratio.
//'
//' * \code{resultsUnderH0} when \code{nullVariance = FALSE}: A data frame
//'   with the following variables:
//'
//'     - \code{time}: The analysis time since trial start.
//'
//'     - \code{rateRatioH0}: The rate ratio under H0.
//'
//'     - \code{varianceRatio}: Equal to 1.
//'
//'     - \code{lambda1}: The true event rate for the active treatment group.
//'
//'     - \code{lambda2}: The true event rate for the control group.
//'
//'     - \code{rateRatio}: The true rate ratio.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//' # Example 1: Variable follow-up design
//'
//' nbstat(time = c(1, 1.25, 2, 3, 4),
//'        accrualIntensity = 1956/1.25,
//'        kappa1 = 5,
//'        kappa2 = 5,
//'        lambda1 = 0.7*0.125,
//'        lambda2 = 0.125,
//'        gamma1 = 0,
//'        gamma2 = 0,
//'        accrualDuration = 1.25,
//'        followupTime = 2.75)
//'
//' # Example 2: Fixed follow-up design
//'
//' nbstat(time = c(0.5, 1, 1.5, 2),
//'        accrualIntensity = 220/1.5,
//'        stratumFraction = c(0.2, 0.8),
//'        kappa1 = 3,
//'        kappa2 = 3,
//'        lambda1 = c(0.5*8.4, 0.6*10.5),
//'        lambda2 = c(8.4, 10.5),
//'        gamma1 = 0,
//'        gamma2 = 0,
//'        accrualDuration = 1.5,
//'        followupTime = 0.5,
//'        fixedFollowup = 1,
//'        nullVariance = 1)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List nbstat(
    const Rcpp::NumericVector& time = NA_REAL,
    const double rateRatioH0 = 1,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& kappa1 = NA_REAL,
    const Rcpp::NumericVector& kappa2 = NA_REAL,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double accrualDuration = NA_REAL,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = false,
    const bool nullVariance = false) {

  auto time1 = Rcpp::as<std::vector<double>>(time);
  auto accrualT = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualInt = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto pwSurvT = Rcpp::as<std::vector<double>>(piecewiseSurvivalTime);
  auto stratumFrac = Rcpp::as<std::vector<double>>(stratumFraction);
  auto kap1 = Rcpp::as<std::vector<double>>(kappa1);
  auto kap2 = Rcpp::as<std::vector<double>>(kappa2);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);

  auto out = nbstatcpp(
    time1, rateRatioH0, allocationRatioPlanned,
    accrualT, accrualInt, pwSurvT, stratumFrac,
    kap1, kap2, lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup, nullVariance);

  return Rcpp::wrap(out);
}


ListCpp nbpowercpp(
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
    const double rateRatioH0,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& stratumFraction,
    const std::vector<double>& kappa1,
    const std::vector<double>& kappa2,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double followupTime,
    const bool fixedFollowup,
    const std::vector<double>& spendingTime,
    const double studyDuration,
    const bool nullVariance) {

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
  if (rateRatioH0 <= 0) throw std::invalid_argument("rateRatioH0 must be positive");
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
  if (!none_na(kappa1)) throw std::invalid_argument("kappa1 must be provided");
  if (!none_na(kappa2)) throw std::invalid_argument("kappa2 must be provided");
  if (!none_na(lambda1)) throw std::invalid_argument("lambda1 must be provided");
  if (!none_na(lambda2)) throw std::invalid_argument("lambda2 must be provided");
  for (double v : kappa1) {
    if (v < 0.0) throw std::invalid_argument("kappa1 must be non-negative");
  }
  for (double v : kappa2) {
    if (v < 0.0) throw std::invalid_argument("kappa2 must be non-negative");
  }
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

  // Expand stratified rates
  size_t nstrata = stratumFraction.size();
  size_t nintv = piecewiseSurvivalTime.size();

  auto kappa1x = expand1(kappa1, nstrata, "kappa1");
  auto kappa2x = expand1(kappa2, nstrata, "kappa2");
  auto lambda1x = expand1(lambda1, nstrata, "lambda1");
  auto lambda2x = expand1(lambda2, nstrata, "lambda2");
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
  std::vector<double> exposure(K), exposure1(K), exposure2(K);
  std::vector<double> I(K);
  std::vector<double> w(K, 1.0); // square root of variance ratio

  // ---- compute maxInformation and theta via kmstat at study end ----
  ListCpp nb = nbstat1cpp(
    studyDuration1, rateRatioH0, allocationRatioPlanned,
    accrualTime, accrualIntensity,
    piecewiseSurvivalTime, stratumFraction,
    kappa1x, kappa2x, lambda1x, lambda2x, gamma1x, gamma2x,
    accrualDuration, followupTime, fixedFollowup, nullVariance);

  DataFrameCpp dfH1 = nb.get<DataFrameCpp>("resultsUnderH1");
  auto rateRatio_v = dfH1.get<double>("rateRatio");
  auto vlogRR_v = dfH1.get<double>("vlogRR");

  double wg_log = 0.0, vvr = 0.0;
  for (size_t h = 0; h < nstrata; ++h) {
    wg_log += stratumFraction[h] * std::log(rateRatio_v[h]);
    vvr += stratumFraction[h] * stratumFraction[h] * vlogRR_v[h];
  }
  double rateRatio = std::exp(wg_log);
  double maxInformation = 1.0 / vvr;
  double theta1 = -std::log(rateRatio / rateRatioH0);
  std::vector<double> theta(K, theta1);

  I[K - 1]           = maxInformation;
  time[K - 1]        = studyDuration1;
  nsubjects[K - 1]   = extract_sum(dfH1, "subjects");
  nsubjects1[K - 1]  = phi * nsubjects[K - 1];
  nsubjects2[K - 1]  = (1.0 - phi) * nsubjects[K - 1];
  nevents[K - 1]     = extract_sum(dfH1, "nevents");
  nevents1[K - 1]    = extract_sum(dfH1, "nevents1");
  nevents2[K - 1]    = extract_sum(dfH1, "nevents2");
  ndropouts[K - 1]   = extract_sum(dfH1, "ndropouts");
  ndropouts1[K - 1]  = extract_sum(dfH1, "ndropouts1");
  ndropouts2[K - 1]  = extract_sum(dfH1, "ndropouts2");
  exposure[K - 1]    = extract_sum(dfH1, "exposure");
  exposure1[K - 1]   = extract_sum(dfH1, "exposure1");
  exposure2[K - 1]   = extract_sum(dfH1, "exposure2");

  if (nullVariance) {
    DataFrameCpp dfH0 = nb.get<DataFrameCpp>("resultsUnderH0");
    auto vlogRRH0_v = dfH0.get<double>("vlogRRH0");
    double vvrH0 = 0.0;
    for (size_t h = 0; h < nstrata; ++h) {
      vvrH0 += stratumFraction[h] * stratumFraction[h] * vlogRRH0_v[h];
    }
    w[K-1] = std::sqrt(vvrH0 / vvr);
  }


  // ---- compute information, time, and other stats at interim analyses ----
  for (size_t i = 0; i < K - 1; ++i) {
    double information1 = maxInformation * infoRates[i];
    I[i] = information1;

    // solve for analysis time where total information equals information1
    auto g = [&](double t)->double {
      auto nb1 = nbstat1cpp(
        t, rateRatioH0, allocationRatioPlanned,
        accrualTime, accrualIntensity,
        piecewiseSurvivalTime, stratumFraction,
        kappa1x, kappa2x, lambda1x, lambda2x, gamma1x, gamma2x,
        accrualDuration, followupTime, fixedFollowup, false);

      DataFrameCpp dfH1 = nb1.get<DataFrameCpp>("resultsUnderH1");
      auto vlogRR_v = dfH1.get<double>("vlogRR");

      double vvr = 0.0;
      for (size_t h = 0; h < nstrata; ++h) {
        vvr += stratumFraction[h] * stratumFraction[h] * vlogRR_v[h];
      }
      return 1.0 / vvr - information1;
    };

    time[i] = brent(g, 0.001, studyDuration1, 1e-6);

    auto nb_i = nbstat1cpp(
      time[i], rateRatioH0, allocationRatioPlanned,
      accrualTime, accrualIntensity,
      piecewiseSurvivalTime, stratumFraction,
      kappa1x, kappa2x, lambda1x, lambda2x, gamma1x, gamma2x,
      accrualDuration, followupTime, fixedFollowup, nullVariance);

    DataFrameCpp dfH1_i = nb_i.get<DataFrameCpp>("resultsUnderH1");

    nsubjects[i]   = extract_sum(dfH1_i, "subjects");
    nsubjects1[i]  = phi * nsubjects[i];
    nsubjects2[i]  = (1.0 - phi) * nsubjects[i];
    nevents[i]     = extract_sum(dfH1_i, "nevents");
    nevents1[i]    = extract_sum(dfH1_i, "nevents1");
    nevents2[i]    = extract_sum(dfH1_i, "nevents2");
    ndropouts[i]   = extract_sum(dfH1_i, "ndropouts");
    ndropouts1[i]  = extract_sum(dfH1_i, "ndropouts1");
    ndropouts2[i]  = extract_sum(dfH1_i, "ndropouts2");
    exposure[i]    = extract_sum(dfH1_i, "exposure");
    exposure1[i]   = extract_sum(dfH1_i, "exposure1");
    exposure2[i]   = extract_sum(dfH1_i, "exposure2");

    if (nullVariance) {
      DataFrameCpp dfH0_i = nb_i.get<DataFrameCpp>("resultsUnderH0");
      auto vlogRRH0_v = dfH0_i.get<double>("vlogRRH0");
      double vvrH0 = 0.0;
      for (size_t h = 0; h < nstrata; ++h) {
        vvrH0 += stratumFraction[h] * stratumFraction[h] * vlogRRH0_v[h];
      }
      w[i] = std::sqrt(vvrH0 * information1);
    }
  }

  // ---- compute exit probabilities ----
  ListCpp exit_probs;
  if (!missingFutilityBounds || bsf == "none" || K == 1) {
    std::vector<double> critValues1 = critValues;
    std::vector<double> futBounds1 = futBounds;
    for (size_t i = 0; i < K; ++i) {
      critValues1[i] *= w[i];
      futBounds1[i] *= w[i];
    }
    exit_probs = exitprobcpp(critValues1, futBounds1, theta, I);
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
  double expectedExposure = 0.0;
  double expectedNumberOfEvents1 = 0.0;
  double expectedNumberOfDropouts1 = 0.0;
  double expectedNumberOfSubjects1 = 0.0;
  double expectedExposure1 = 0.0;
  double expectedNumberOfEvents2 = 0.0;
  double expectedNumberOfDropouts2 = 0.0;
  double expectedNumberOfSubjects2 = 0.0;
  double expectedExposure2 = 0.0;
  double expectedStudyDuration = 0.0;
  double expectedInformation = 0.0;
  for (size_t i = 0; i < K; ++i) {
    expectedNumberOfEvents += ptotal[i] * nevents[i];
    expectedNumberOfDropouts += ptotal[i] * ndropouts[i];
    expectedNumberOfSubjects += ptotal[i] * nsubjects[i];
    expectedExposure += ptotal[i] * exposure[i];
    expectedNumberOfEvents1 += ptotal[i] * nevents1[i];
    expectedNumberOfDropouts1 += ptotal[i] * ndropouts1[i];
    expectedNumberOfSubjects1 += ptotal[i] * nsubjects1[i];
    expectedExposure1 += ptotal[i] * exposure1[i];
    expectedNumberOfEvents2 += ptotal[i] * nevents2[i];
    expectedNumberOfDropouts2 += ptotal[i] * ndropouts2[i];
    expectedNumberOfSubjects2 += ptotal[i] * nsubjects2[i];
    expectedExposure2 += ptotal[i] * exposure2[i];
    expectedStudyDuration += ptotal[i] * time[i];
    expectedInformation += ptotal[i] * I[i];
  }

  std::vector<double> cpu(K), cpl(K);
  std::partial_sum(pu.begin(), pu.end(), cpu.begin());
  std::partial_sum(pl.begin(), pl.end(), cpl.begin());
  double overallReject = cpu[K-1];

  // efficacy/futility implied survDiff boundaries at each stage
  std::vector<double> rru(K), rrl(K);
  for (size_t i = 0; i < K; ++i) {
    rru[i] = rateRatioH0 * std::exp(-critValues[i] / std::sqrt(I[i]) * w[i]);
    rrl[i] = rateRatioH0 * std::exp(-futBounds[i] / std::sqrt(I[i]) * w[i]);
    if (critValues[i] == 6.0) { rru[i] = NaN; effStopping[i] = 0; }
    if (futBounds[i] == -6.0) { rrl[i] = NaN; futStopping[i] = 0; }
  }

  // --- Build output DataFrames and Lists ------------------------------------
  DataFrameCpp overallResults;
  overallResults.push_back(overallReject, "overallReject");
  overallResults.push_back(alpha1, "alpha");
  overallResults.push_back(nevents.back(), "numberOfEvents");
  overallResults.push_back(ndropouts.back(), "numberOfDropouts");
  overallResults.push_back(nsubjects.back(), "numberOfSubjects");
  overallResults.push_back(exposure.back(), "exposure");
  overallResults.push_back(time.back(), "studyDuration");
  overallResults.push_back(maxInformation, "information");
  overallResults.push_back(expectedNumberOfEvents, "expectedNumberOfEvents");
  overallResults.push_back(expectedNumberOfDropouts, "expectedNumberOfDropouts");
  overallResults.push_back(expectedNumberOfSubjects, "expectedNumberOfSubjects");
  overallResults.push_back(expectedExposure, "expectedExposure");
  overallResults.push_back(expectedStudyDuration, "expectedStudyDuration");
  overallResults.push_back(expectedInformation, "expectedInformation");
  overallResults.push_back(accrualDuration, "accrualDuration");
  overallResults.push_back(followupTime, "followupTime");
  overallResults.push_back(fixedFollowup, "fixedFollowup");
  overallResults.push_back(kMax, "kMax");
  overallResults.push_back(rateRatioH0, "rateRatioH0");
  overallResults.push_back(rateRatio, "rateRatio");

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
  byStageResults.push_back(std::move(exposure), "exposure");
  byStageResults.push_back(std::move(time), "analysisTime");
  byStageResults.push_back(std::move(rru), "efficacyRateRatio");
  byStageResults.push_back(std::move(rrl), "futilityRateRatio");
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
  settings.push_back(kappa1, "kappa1");
  settings.push_back(kappa2, "kappa2");
  settings.push_back(lambda1, "lambda1");
  settings.push_back(lambda2, "lambda2");
  settings.push_back(gamma1, "gamma1");
  settings.push_back(gamma2, "gamma2");
  settings.push_back(spendingTime, "spendingTime");
  settings.push_back(nullVariance, "nullVariance");

  ListCpp byTreatmentCounts;
  byTreatmentCounts.push_back(std::move(nevents1), "numberOfEvents1");
  byTreatmentCounts.push_back(std::move(ndropouts1), "numberOfDropouts1");
  byTreatmentCounts.push_back(std::move(nsubjects1), "numberOfSubjects1");
  byTreatmentCounts.push_back(std::move(exposure1), "exposure1");
  byTreatmentCounts.push_back(std::move(nevents2), "numberOfEvents2");
  byTreatmentCounts.push_back(std::move(ndropouts2), "numberOfDropouts2");
  byTreatmentCounts.push_back(std::move(nsubjects2), "numberOfSubjects2");
  byTreatmentCounts.push_back(std::move(exposure2), "exposure2");
  byTreatmentCounts.push_back(expectedNumberOfEvents1, "expectedNumberOfEvents1");
  byTreatmentCounts.push_back(expectedNumberOfDropouts1, "expectedNumberOfDropouts1");
  byTreatmentCounts.push_back(expectedNumberOfSubjects1, "expectedNumberOfSubjects1");
  byTreatmentCounts.push_back(expectedExposure1, "expectedExposure1");
  byTreatmentCounts.push_back(expectedNumberOfEvents2, "expectedNumberOfEvents2");
  byTreatmentCounts.push_back(expectedNumberOfDropouts2, "expectedNumberOfDropouts2");
  byTreatmentCounts.push_back(expectedNumberOfSubjects2, "expectedNumberOfSubjects2");
  byTreatmentCounts.push_back(expectedExposure2, "expectedExposure2");

  ListCpp result;
  result.push_back(std::move(byStageResults), "byStageResults");
  result.push_back(std::move(overallResults), "overallResults");
  result.push_back(std::move(settings), "settings");
  result.push_back(std::move(byTreatmentCounts), "byTreatmentCounts");

  return result;
}


//' @title Power for Negative Binomial Rate Ratio
//' @description Estimates the power for negative binomial rate ratio test.
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
//' @param rateRatioH0 Rate ratio under the null hypothesis.
//' @inheritParams param_allocationRatioPlanned
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_piecewiseSurvivalTime
//' @inheritParams param_stratumFraction
//' @param kappa1 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the active treatment
//'   group by stratum.
//' @param kappa2 The dispersion parameter (reciprocal of the shape
//'   parameter of the gamma mixing distribution) for the control group by
//'   stratum.
//' @param lambda1 The rate parameter of the negative binomial distribution
//'   for the active treatment group by stratum.
//' @param lambda2 The rate parameter of the negative binomial distribution
//'   for the control group by stratum.
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
//' @param nullVariance Whether to calculate the variance for log rate ratio
//'   under the null hypothesis.
//'
//' @return An S3 class \code{nbpower} object with 4 components:
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
//'     - \code{exposure}: The total exposure.
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
//'     - \code{expectedExposure}: The expected exposure.
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
//'     - \code{rateRatioH0}: The rate ratio under the null hypothesis.
//'
//'     - \code{rateRatio}: The rate ratio.
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
//'     - \code{exposure}: The exposure.
//'
//'     - \code{analysisTime}: The average time since trial start.
//'
//'     - \code{efficacyRateRatio}: The efficacy boundaries on the rate
//'       ratio scale.
//'
//'     - \code{futilityRateRatio}: The futility boundaries on the rate
//'       ratio scale.
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
//'   \code{piecewiseSurvivalTime}, \code{kappa1}, \code{kappa2},
//'   \code{lambda1}, \code{lambda2}, \code{gamma1}, \code{gamma2},
//'   \code{spendingTime}, and \code{nullVariance}.
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
//'     - \code{exposure1}: The exposure by stage for the treatment group.
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
//'     - \code{exposure2}: The exposure by stage for the control group.
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
//'     - \code{expectedExposure1}: The expected exposure for the treatment
//'       group.
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
//'     - \code{expectedExposure2}: The expected exposure for the control
//'       group.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @seealso \code{\link{nbstat}}
//'
//' @examples
//' # Example 1: Variable follow-up design
//'
//' nbpower(kMax = 2, informationRates = c(0.5, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         accrualIntensity = 1956/1.25,
//'         stratumFraction = c(0.2, 0.8),
//'         kappa1 = 5, kappa2 = 5,
//'         lambda1 = c(0.7*0.125, 0.75*0.25),
//'         lambda2 = c(0.125, 0.25),
//'         gamma1 = 0, gamma2 = 0,
//'         accrualDuration = 1.25,
//'         followupTime = 2.75, fixedFollowup = FALSE,
//'         nullVariance = 1)
//'
//' # Example 2: Fixed follow-up design
//'
//' nbpower(kMax = 2, informationRates = c(0.5, 1),
//'         alpha = 0.025, typeAlphaSpending = "sfOF",
//'         accrualIntensity = 220/1.5,
//'         kappa1 = 3, kappa2 = 3,
//'         lambda1 = 0.5*8.4, lambda2 = 8.4,
//'         gamma1 = 0, gamma2 = 0,
//'         accrualDuration = 1.5,
//'         followupTime = 0.5, fixedFollowup = TRUE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List nbpower(
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
    const double rateRatioH0 = 1,
    const double allocationRatioPlanned = 1,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& piecewiseSurvivalTime = 0,
    const Rcpp::NumericVector& stratumFraction = 1,
    const Rcpp::NumericVector& kappa1 = NA_REAL,
    const Rcpp::NumericVector& kappa2 = NA_REAL,
    const Rcpp::NumericVector& lambda1 = NA_REAL,
    const Rcpp::NumericVector& lambda2 = NA_REAL,
    const Rcpp::NumericVector& gamma1 = 0,
    const Rcpp::NumericVector& gamma2 = 0,
    const double accrualDuration = NA_REAL,
    const double followupTime = NA_REAL,
    const bool fixedFollowup = false,
    const Rcpp::NumericVector& spendingTime = NA_REAL,
    const double studyDuration = NA_REAL,
    const bool nullVariance = false) {

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
  auto kap1 = Rcpp::as<std::vector<double>>(kappa1);
  auto kap2 = Rcpp::as<std::vector<double>>(kappa2);
  auto lam1 = Rcpp::as<std::vector<double>>(lambda1);
  auto lam2 = Rcpp::as<std::vector<double>>(lambda2);
  auto gam1 = Rcpp::as<std::vector<double>>(gamma1);
  auto gam2 = Rcpp::as<std::vector<double>>(gamma2);
  auto spendTime = Rcpp::as<std::vector<double>>(spendingTime);

  ListCpp out = nbpowercpp(
    kMax, infoRates, effStopping, futStopping,
    critValues, alpha, typeAlphaSpending,
    parameterAlphaSpending, userAlpha,
    futBounds, typeBetaSpending, parameterBetaSpending,
    rateRatioH0, allocationRatioPlanned,
    accrualT, accrualInt, piecewiseTime, stratumFrac,
    kap1, kap2, lam1, lam2, gam1, gam2,
    accrualDuration, followupTime, fixedFollowup,
    spendTime, studyDuration, nullVariance);

  Rcpp::List result = Rcpp::wrap(out);
  result.attr("class") = "nbpower";
  return result;
}


