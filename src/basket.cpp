#include <Rcpp.h>
#include <boost/random.hpp>

#include "utilities.h"
#include "dataframe_list.h"

#include <algorithm>  // fill, lower_bound, max_element, memmove,
// min_element, sort, swap, upper_bound
#include <cmath>      // copysign, exp, fabs, isinf, isnan, log, pow, sqrt
#include <cstddef>    // size_t
#include <cstring>    // memcpy
#include <functional> // function
#include <limits>     // numeric_limits
#include <memory>     // make_shared, shared_ptr
#include <numeric>    // accumulate, inner_product, iota
#include <queue>      // priority_queue
#include <sstream>    // ostringstream
#include <stdexcept>  // invalid_argument, runtime_error
#include <string>     // string
#include <utility>    // pair, swap
#include <vector>     // vector


ListCpp simonBayesAnalysiscpp(
    const int nstrata,
    const std::vector<double>& r,
    const std::vector<double>& n,
    const double lambda,
    const double gamma,
    const double phi,
    const double plo) {

  if (nstrata == INT_MIN) {
    throw std::invalid_argument("nstrata must be provided.");
  }
  if (nstrata <= 0) {
    throw std::invalid_argument("nstrata must be a positive integer.");
  }
  if (!none_na(r)) {
    throw std::invalid_argument("r must be provided.");
  }
  if (!none_na(n)) {
    throw std::invalid_argument("n must be provided.");
  }
  if (static_cast<int>(r.size()) != nstrata) {
    throw std::invalid_argument("Invalid length for r.");
  }
  if (static_cast<int>(n.size()) != nstrata) {
    throw std::invalid_argument("Invalid length for n.");
  }
  if (std::any_of(r.begin(), r.end(), [](double val) { return val < 0.0; })) {
    throw std::invalid_argument("r must be nonnegative.");
  }

  for (size_t i = 0; i < r.size(); ++i) {
    if (r[i] > n[i]) {
      throw std::invalid_argument("r must be less than or equal to n.");
    }
  }

  if (std::isnan(lambda)) {
    throw std::invalid_argument("lambda must be provided.");
  }
  if (std::isnan(gamma)) {
    throw std::invalid_argument("gamma must be provided.");
  }
  if (std::isnan(phi)) {
    throw std::invalid_argument("phi must be provided.");
  }
  if (std::isnan(plo)) {
    throw std::invalid_argument("plo must be provided.");
  }
  if (lambda <= 0 || lambda >= 1) {
    throw std::invalid_argument("lambda must lie between 0 and 1.");
  }
  if (gamma <= 0 || gamma >= 1) {
    throw std::invalid_argument("gamma must lie between 0 and 1.");
  }
  if (phi <= 0 || phi >= 1) {
    throw std::invalid_argument("phi must lie between 0 and 1.");
  }
  if (plo <= 0 || plo >= 1) {
    throw std::invalid_argument("plo must lie between 0 and 1.");
  }
  if (plo >= phi) {
    throw std::invalid_argument("plo must be less than phi.");
  }

  int ncases = static_cast<int>(std::pow(2, nstrata));
  FlatMatrix incid(ncases, nstrata);
  std::vector<double> prior(ncases), like(ncases);
  for (int i=0; i<ncases; ++i) {
    int number = ncases - i - 1;
    std::vector<double> cc(nstrata);
    for (int j=0; j<nstrata; ++j) {
      cc[j] = (number / static_cast<int>(std::pow(2, nstrata - 1 - j))) % 2;
    }

    bool all_ones = std::all_of(cc.begin(), cc.end(), [](int x) { return x == 1; });
    bool all_zeros = std::all_of(cc.begin(), cc.end(), [](int x) { return x == 0; });

    if (all_ones || all_zeros) {
      prior[i] = lambda * (gamma * (cc[0] == 1) + (1 - gamma) * (cc[0] == 0)) +
        (1 - lambda) * (std::pow(gamma, nstrata) * (cc[0] == 1) +
        std::pow(1 - gamma, nstrata) * (cc[0]==0));
    } else {
      double y = std::accumulate(cc.begin(), cc.end(), 0.0);
      prior[i] = (1 - lambda) * std::pow(gamma, y) *
        std::pow(1 - gamma, nstrata - y);
    }

    int m = static_cast<int>(cc.size());
    std::vector<double> x(m);
    for (int j = 0; j < m; ++j) {
      x[j] = phi * cc[j] + plo * (1 - cc[j]);
    }

    // Calculate sum of r*log(x) + (n-r)*log(1-x)
    double log_sum = 0.0;
    for (int j = 0; j < m; ++j) {
      log_sum += r[j] * std::log(x[j]) + (n[j] - r[j]) * std::log(1 - x[j]);
    }
    like[i] = std::exp(log_sum);

    // Copy cc to row i of incid matrix
    for (int j = 0; j < m; ++j) {
      incid(i, j) = cc[j];
    }
  }

  std::vector<double> post(ncases);
  for (int i = 0; i < ncases; ++i) {
    post[i] = prior[i] * like[i];
  }

  // Normalize: post = post / sum(post)
  double post_sum = std::accumulate(post.begin(), post.end(), 0.0);
  for (int i = 0; i < ncases; ++i) {
    post[i] /= post_sum;
  }

  // Calculate q_prior and q_post for each stratum
  std::vector<double> q_prior(nstrata, 0.0);
  std::vector<double> q_post(nstrata, 0.0);
  for (int j = 0; j < nstrata; ++j) {
    for (int i = 0; i < ncases; ++i) {
      q_prior[j] += incid(i, j) * prior[i];
      q_post[j] += incid(i, j) * post[i];
    }
  }

  ListCpp result;
  result.push_back(std::move(incid), "case");
  result.push_back(std::move(prior), "prior_case");
  result.push_back(std::move(q_prior), "prior_stratum");
  result.push_back(std::move(post), "post_case");
  result.push_back(std::move(q_post), "post_stratum");
  return result;
}

//' @title Analysis of Simon's Bayesian Basket Trials
//' @description Obtains the prior and posterior probabilities for
//' Simon's Bayesian basket discovery trials.
//'
//' @param nstrata The number of strata.
//' @param r The vector of number of responders across strata.
//' @param n The vector of number of subjects across strata.
//' @param lambda The prior probability that the drug activity is
//'   homogeneous across strata.
//' @param gamma The prior probability that the drug is active in a
//'   stratum.
//' @param phi The response probability for an active drug.
//' @param plo The response probability for an inactive drug.
//'
//' @return A list containing the following five components:
//'
//' * \code{case}: The matrix with each row corresponding to a combination
//'   of drug activity over strata represented by the columns.
//'
//' * \code{prior_case}: The vector of joint prior probabilities
//'   for the stratum-specific response rates.
//'
//' * \code{prior_stratum}: The vector of marginal prior probabilities
//'   for the stratum-specific response rates.
//'
//' * \code{post_case}: The vector of joint posterior probabilities
//'   for the stratum-specific response rates.
//'
//' * \code{post_stratum}: The vector of marginal posterior probabilities
//'   for the stratum-specific response rates.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' a = simonBayesAnalysis(
//'   nstrata = 10,
//'   r = c(8,0,1,1,6,2,0,0,3,3),
//'   n = c(19,10,26,8,14,7,8,5,4,14),
//'   lambda = 0.5, gamma = 0.33,
//'   phi = 0.35, plo = 0.15)
//'
//' a$post_stratum
//'
//' @export
// [[Rcpp::export]]
Rcpp::List simonBayesAnalysis(
    const int nstrata = NA_INTEGER,
    const Rcpp::NumericVector& r = NA_REAL,
    const Rcpp::NumericVector& n = NA_REAL,
    const double lambda = NA_REAL,
    const double gamma = NA_REAL,
    const double phi = NA_REAL,
    const double plo = NA_REAL) {
  auto r1 = Rcpp::as<std::vector<double>>(r);
  auto n1 = Rcpp::as<std::vector<double>>(n);
  ListCpp result = simonBayesAnalysiscpp(nstrata, r1, n1, lambda, gamma, phi, plo);
  return Rcpp::wrap(result);
}


ListCpp simonBayesSimcpp(
    const std::vector<double>& p,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& stratumFraction,
    const double lambda,
    const double gamma,
    const double phi,
    const double plo,
    const double T,
    const int maxSubjects,
    const std::vector<int>& plannedSubjects,
    const int maxNumberOfIterations,
    const int maxNumberOfRawDatasets,
    const int seed) {

  int nstrata = static_cast<int>(stratumFraction.size());

  if (!none_na(p)) {
    throw std::invalid_argument("p must be provided.");
  }
  if (std::any_of(p.begin(), p.end(), [](double x) { return x <= 0 || x >= 1; })) {
    throw std::invalid_argument("p must lie between 0 and 1.");
  }

  if (accrualTime[0] != 0) {
    throw std::invalid_argument("accrualTime must start with 0.");
  }
  if (any_nonincreasing(accrualTime)) {
    throw std::invalid_argument("accrualTime should be increasing.");
  }
  if (!none_na(accrualIntensity)) {
    throw std::invalid_argument("accrualIntensity must be provided.");
  }
  if (accrualTime.size() != accrualIntensity.size()) {
    throw std::invalid_argument(
        "accrualTime must have the same length as accrualIntensity.");
  }
  if (std::any_of(accrualIntensity.begin(), accrualIntensity.end(),
                  [](double val) { return val < 0.0; })) {
    throw std::invalid_argument("accrualIntensity must be non-negative.");
  }
  if (std::any_of(stratumFraction.begin(), stratumFraction.end(),
                  [](double val) { return val <= 0.0; })) {
    throw std::invalid_argument("stratumFraction must be positive.");
  }

  double sum_stratumFraction = std::accumulate(stratumFraction.begin(),
                                               stratumFraction.end(), 0.0);
  if (std::fabs(sum_stratumFraction - 1.0) > 1.0e-8) {
    throw std::invalid_argument("stratumFraction must sum to 1.");
  }

  if (std::isnan(lambda)) {
    throw std::invalid_argument("lambda must be provided.");
  }
  if (std::isnan(gamma)) {
    throw std::invalid_argument("gamma must be provided.");
  }
  if (std::isnan(phi)) {
    throw std::invalid_argument("phi must be provided.");
  }
  if (std::isnan(plo)) {
    throw std::invalid_argument("plo must be provided.");
  }
  if (std::isnan(T)) {
    throw std::invalid_argument("T must be provided.");
  }
  if (lambda <= 0 || lambda >= 1) {
    throw std::invalid_argument("lambda must lie between 0 and 1.");
  }
  if (gamma <= 0 || gamma >= 1) {
    throw std::invalid_argument("gamma must lie between 0 and 1.");
  }
  if (phi <= 0 || phi >= 1) {
    throw std::invalid_argument("phi must lie between 0 and 1.");
  }
  if (plo <= 0 || plo >= 1) {
    throw std::invalid_argument("plo must lie between 0 and 1.");
  }
  if (plo >= phi) {
    throw std::invalid_argument("plo must be less than phi.");
  }
  if (T <= 0 || T >= 1) {
    throw std::invalid_argument("T must lie between 0 and 1.");
  }
  if (maxSubjects == INT_MIN) {
    throw std::invalid_argument("maxSubjects must be provided.");
  }
  if (maxSubjects < 1) {
    throw std::invalid_argument("maxSubjects must be a positive integer.");
  }
  if (!none_na(plannedSubjects)) {
    throw std::invalid_argument("plannedSubjects must be provided.");
  }
  if (plannedSubjects[0] <= 0) {
    throw std::invalid_argument("plannedSubjects must be positive.");
  }
  if (any_nonincreasing(plannedSubjects)) {
    throw std::invalid_argument("plannedSubjects must be increasing.");
  }
  if (maxNumberOfIterations < 1) {
    throw std::invalid_argument("maxNumberOfIterations must be positive.");
  }
  if (maxNumberOfRawDatasets < 0) {
    throw std::invalid_argument("maxNumberOfRawDatasets must be non-negative.");
  }


  std::vector<unsigned char> act(p.size());
  for (size_t i = 0; i < p.size(); ++i) {
    act[i] = (p[i] == phi) ? 1 : 0;
  }
  int nactive = std::accumulate(act.begin(), act.end(), 0);

  std::vector<double> cumStratumFraction(nstrata);
  std::partial_sum(stratumFraction.begin(), stratumFraction.end(),
                   cumStratumFraction.begin());
  std::vector<double> post_stratum(nstrata), n(nstrata), r(nstrata);
  std::vector<unsigned char> open(nstrata), pos(nstrata), neg(nstrata);


  std::vector<double> arrivalTime(maxSubjects);
  std::vector<int> stratum(maxSubjects), y(maxSubjects);
  std::vector<int> iterationNumber(maxNumberOfIterations);
  std::vector<double> N(maxNumberOfIterations);
  std::vector<double> nact(maxNumberOfIterations), nopen(maxNumberOfIterations);
  std::vector<double> tpos(maxNumberOfIterations), fneg(maxNumberOfIterations);
  std::vector<double> fpos(maxNumberOfIterations), tneg(maxNumberOfIterations);
  std::vector<int> numberOfStrata(maxNumberOfIterations, nstrata);

  // cache for the patient-level raw data to extract
  int kMax = static_cast<int>(plannedSubjects.size());
  int nrow1 = kMax * maxNumberOfRawDatasets * maxSubjects;
  std::vector<int> iterationNumberx(nrow1);
  std::vector<int> stageNumberx(nrow1);
  std::vector<int> subjectIdx(nrow1);
  std::vector<double> arrivalTimex(nrow1);
  std::vector<int> stratumx(nrow1);
  std::vector<int> yx(nrow1);

  // cache for the summary data to extract
  int nrow2 = kMax * maxNumberOfIterations * nstrata;
  std::vector<int> iterationNumbery(nrow2);
  std::vector<int> stageNumbery(nrow2);
  std::vector<int> stratumy(nrow2);
  std::vector<unsigned char> activey(nrow2);
  std::vector<int> ny(nrow2);
  std::vector<int> ry(nrow2);
  std::vector<double> posty(nrow2);
  std::vector<unsigned char> openy(nrow2);
  std::vector<unsigned char> posy(nrow2);
  std::vector<unsigned char> negy(nrow2);

  // random number generator
  boost::random::mt19937_64 rng(seed);
  boost::random::uniform_real_distribution<double> unif(0.0, 1.0);

  int index1 = 0, index2 = 0;
  for (int iter = 0; iter < maxNumberOfIterations; ++iter) {
    // initialize the contents in each stratum
    std::fill(n.begin(), n.end(), 0.0);
    std::fill(r.begin(), r.end(), 0.0);
    std::fill(open.begin(), open.end(), 1);
    std::fill(pos.begin(), pos.end(), 0);
    std::fill(neg.begin(), neg.end(), 0);

    int k = 0;      // index of the number of subjects included in analysis
    int stage = 0;
    double enrollt = 0.0;
    for (int i = 0; i < 100000; ++i) {
      // generate accrual time
      double u = unif(rng);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);

      // generate stratum information
      u = unif(rng);
      int j = 0;
      for (; j < nstrata; ++j) {
        if (cumStratumFraction[j] > u) break;
      }

      // if the stratum is open, generate the response for the subject
      if (open[j]) {
        arrivalTime[k] = enrollt;
        stratum[k] = j+1;
        y[k] = (unif(rng) < p[j]) ? 1 : 0;

        // update the number of subjects and responders in the stratum
        n[j] = n[j] + 1;
        r[j] = r[j] + y[k];

        ++k;

        // interim analysis
        if (std::any_of(plannedSubjects.begin(), plannedSubjects.end(),
                        [k](double val) { return val == k; })) {
          // output raw data
          if (iter < maxNumberOfRawDatasets) {
            for (int idx = 0; idx < k; ++idx) {
              iterationNumberx[index1] = iter + 1;
              stageNumberx[index1] = stage + 1;
              subjectIdx[index1] = idx + 1;
              arrivalTimex[index1] = arrivalTime[idx];
              stratumx[index1] = stratum[idx];
              yx[index1] = y[idx];

              ++index1;
            }
          }

          // calculate the posterior probabilities
          ListCpp a = simonBayesAnalysiscpp(nstrata, r, n, lambda, gamma, phi, plo);
          post_stratum = a.get<std::vector<double>>("post_stratum");

          // whether to close the stratum due to positive or negative results
          for (int l = 0; l < nstrata; ++l) {
            if (open[l]) {
              if (post_stratum[l] > T) {
                pos[l] = 1; open[l] = 0;
              } else if (post_stratum[l] < 1 - T) {
                neg[l] = 1; open[l] = 0;
              }
            }

            // output strata-level data
            iterationNumbery[index2] = iter + 1;
            stageNumbery[index2] = stage + 1;
            stratumy[index2] = l + 1;
            activey[index2] = act[l];
            ny[index2] = n[l];
            ry[index2] = r[l];
            posty[index2] = post_stratum[l];
            openy[index2] = open[l];
            posy[index2] = pos[l];
            negy[index2] = neg[l];

            ++index2;
          }

          ++stage;
        }


        // stop the trial if all strata are closed or max subjects reached
        // Check if all elements in open are 0, or k reached maxSubjects
        bool all_closed = std::all_of(open.begin(), open.end(),
                                      [](unsigned char val) { return val == 0; });

        if (all_closed || (k == maxSubjects)) {
          iterationNumber[iter] = iter + 1;
          N[iter] = k;
          nact[iter] = nactive;

          // Calculate true positives: sum(pos & act)
          int sum_pos_act = 0;
          for (int i = 0; i < nstrata; ++i) {
            sum_pos_act += pos[i] & act[i];
          }
          tpos[iter] = sum_pos_act;

          // Calculate false negatives: sum(neg & act)
          int sum_neg_act = 0;
          for (int i = 0; i < nstrata; ++i) {
            sum_neg_act += neg[i] & act[i];
          }
          fneg[iter] = sum_neg_act;

          // Calculate false positives: sum(pos & !act)
          int sum_pos_notact = 0;
          for (int i = 0; i < nstrata; ++i) {
            sum_pos_notact += pos[i] & !act[i];
          }
          fpos[iter] = sum_pos_notact;

          // Calculate true negatives: sum(neg & !act)
          int sum_neg_notact = 0;
          for (int i = 0; i < nstrata; ++i) {
            sum_neg_notact += neg[i] & !act[i];
          }
          tneg[iter] = sum_neg_notact;

          // Calculate number of open: sum(open)
          int sum_open = std::accumulate(open.begin(), open.end(), 0);
          nopen[iter] = sum_open;

          break;
        }
      }
    }
  }

  // subject-level raw data set
  DataFrameCpp rawdata;
  if (maxNumberOfRawDatasets > 0) {
    subset_in_place(iterationNumberx, 0, index1);
    subset_in_place(stageNumberx, 0, index1);
    subset_in_place(subjectIdx, 0, index1);
    subset_in_place(arrivalTimex, 0, index1);
    subset_in_place(stratumx, 0, index1);
    subset_in_place(yx, 0, index1);

    rawdata.push_back(std::move(iterationNumberx), "iterationNumber");
    rawdata.push_back(std::move(stageNumberx), "stageNumber");
    rawdata.push_back(std::move(subjectIdx), "subjectId");
    rawdata.push_back(std::move(arrivalTimex), "arrivalTime");
    rawdata.push_back(std::move(stratumx), "stratum");
    rawdata.push_back(std::move(yx), "y");
  }

  // simulation summary data set
  subset_in_place(iterationNumbery, 0, index2);
  subset_in_place(stageNumbery, 0, index2);
  subset_in_place(stratumy, 0, index2);
  subset_in_place(activey, 0, index2);
  subset_in_place(ny, 0, index2);
  subset_in_place(ry, 0, index2);
  subset_in_place(posty, 0, index2);
  subset_in_place(openy, 0, index2);
  subset_in_place(posy, 0, index2);
  subset_in_place(negy, 0, index2);

  DataFrameCpp sumdata1;
  sumdata1.push_back(std::move(iterationNumbery), "iterationNumber");
  sumdata1.push_back(std::move(stageNumbery), "stageNumber");
  sumdata1.push_back(std::move(stratumy), "stratum");
  sumdata1.push_back(std::move(activey), "active");
  sumdata1.push_back(std::move(ny), "n");
  sumdata1.push_back(std::move(ry), "r");
  sumdata1.push_back(std::move(posty), "posterior");
  sumdata1.push_back(std::move(openy), "open");
  sumdata1.push_back(std::move(posy), "positive");
  sumdata1.push_back(std::move(negy), "negative");

  DataFrameCpp sumdata2;
  sumdata2.push_back(std::move(iterationNumber), "iterationNumber");
  sumdata2.push_back(std::move(numberOfStrata), "numberOfStrata");
  sumdata2.push_back(std::move(nact), "n_active_strata");
  sumdata2.push_back(std::move(tpos), "true_positive");
  sumdata2.push_back(std::move(fneg), "false_negative");
  sumdata2.push_back(std::move(fpos), "false_positive");
  sumdata2.push_back(std::move(tneg), "true_negative");
  sumdata2.push_back(std::move(nopen), "n_indet_strata");
  sumdata2.push_back(std::move(N), "numberOfSubjects");

  double mn_nact = mean_kahan(nact);
  double mn_tpos = mean_kahan(tpos);
  double mn_fneg = mean_kahan(fneg);
  double mn_fpos = mean_kahan(fpos);
  double mn_tneg = mean_kahan(tneg);
  double mn_inde = mean_kahan(nopen);
  double mn_N = mean_kahan(N);

  DataFrameCpp overview;
  overview.push_back(nstrata, "numberOfStrata");
  overview.push_back(mn_nact, "n_active_strata");
  overview.push_back(mn_tpos, "true_positive");
  overview.push_back(mn_fneg, "false_negative");
  overview.push_back(mn_fpos, "false_positive");
  overview.push_back(mn_tneg, "true_negative");
  overview.push_back(mn_inde, "n_indet_strata");
  overview.push_back(mn_N, "numberOfSubjects");

  ListCpp result;
  if (maxNumberOfRawDatasets > 0) {
    result.push_back(std::move(rawdata), "rawdata");
    result.push_back(std::move(sumdata1), "sumdata1");
    result.push_back(std::move(sumdata2), "sumdata2");
    result.push_back(std::move(overview), "overview");
  } else {
    result.push_back(std::move(sumdata1), "sumdata1");
    result.push_back(std::move(sumdata2), "sumdata2");
    result.push_back(std::move(overview), "overview");
  }
  return result;
}


//' @title Simulation of Simon's Bayesian Basket Trials
//' @description Obtains the simulated raw and summary data for Simon's
//' Bayesian basket discovery trials.
//'
//' @param p The vector of true response probabilities across strata.
//' @inheritParams param_accrualTime
//' @inheritParams param_accrualIntensity
//' @inheritParams param_stratumFraction
//' @param lambda The prior probability that the drug activity is
//'   homogeneous across strata.
//' @param gamma The prior probability that the drug is active in a
//'   stratum.
//' @param phi The response probability for an active drug.
//' @param plo The response probability for an inactive drug.
//' @param T The threshold for a conclusive posterior probability to
//'   stop enrollment.
//' @param maxSubjects The maximum total sample size.
//' @param plannedSubjects The planned cumulative number of subjects
//'   at each stage.
//' @param maxNumberOfIterations The number of simulation iterations.
//'   Defaults to 1000.
//' @param maxNumberOfRawDatasets The number of raw datasets to extract.
//' @param seed The seed to reproduce the simulation results.
//'   The seed from the environment will be used if left unspecified,
//'
//' @return A list containing the following four components:
//'
//' * \code{rawdata}: A data frame for subject-level data, containing
//'   the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage number.
//'
//'     - \code{subjectId}: The subject ID.
//'
//'     - \code{arrivalTime}: The enrollment time for the subject.
//'
//'     - \code{stratum}: The stratum for the subject.
//'
//'     - \code{y}: Whether the subject was a responder (1) or
//'       nonresponder (0).
//'
//' * \code{sumdata1}: A data frame for simulation and stratum-level
//'   summary data, containing the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{stageNumber}: The stage number.
//'
//'     - \code{stratum}: The stratum number.
//'
//'     - \code{active}: Whether the drug is active in the stratum.
//'
//'     - \code{n}: The number of subjects in the stratum.
//'
//'     - \code{r}: The number of responders in the stratum.
//'
//'     - \code{posterior}: The posterior probability that the drug is
//'       active in the stratum.
//'
//'     - \code{open}: Whether the stratum is still open for enrollment.
//'
//'     - \code{positive}: Whether the stratum has been determined to be
//'       a positive stratum.
//'
//'     - \code{negative}: Whether the stratum has been determined to be
//'       a negative stratum.
//'
//' * \code{sumdata2}: A data frame for the simulation level summary data,
//'   containing the following variables:
//'
//'     - \code{iterationNumber}: The iteration number.
//'
//'     - \code{numberOfStrata}: The total number of strata.
//'
//'     - \code{n_active_strata}: The number of active strata.
//'
//'     - \code{true_positive}: The number of true positive strata.
//'
//'     - \code{false_negative}: The number of false negative strata.
//'
//'     - \code{false_positive}: The number of false positive strata.
//'
//'     - \code{true_negative}: The number of true negative strata.
//'
//'     - \code{n_indet_strata}: The number of indeterminate strata.
//'
//'     - \code{numberOfSubjects}: The number of subjects.
//'
//' * \code{overview}: A data frame for the summary across simulations,
//'   containing the following variables:
//'
//'     - \code{numberOfStrata}: The total number of strata.
//'
//'     - \code{n_active_strata}: The average number of active strata.
//'
//'     - \code{true_positive}: The average number of true positive strata.
//'
//'     - \code{false_negative}: The average number of false negative strata.
//'
//'     - \code{false_positive}: The average number of false positive strata.
//'
//'     - \code{true_negative}: The average number of true negative strata.
//'
//'     - \code{n_indet_strata}: The average number of indeterminate strata.
//'
//'     - \code{numberOfSubjects}: The average number of subjects.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' sim1 = simonBayesSim(
//'   p = c(0.25, 0.25, 0.05),
//'   accrualIntensity = 5,
//'   stratumFraction = c(1/3, 1/3, 1/3),
//'   lambda = 0.33, gamma = 0.5,
//'   phi = 0.25, plo = 0.05,
//'   T = 0.8, maxSubjects = 50,
//'   plannedSubjects = seq(5, 50, 5),
//'   maxNumberOfIterations = 1000,
//'   maxNumberOfRawDatasets = 1,
//'   seed = 314159)
//'
//' sim1$overview
//'
//' @export
// [[Rcpp::export]]
Rcpp::List simonBayesSim(
    const Rcpp::NumericVector& p = NA_REAL,
    const Rcpp::NumericVector& accrualTime = 0,
    const Rcpp::NumericVector& accrualIntensity = NA_REAL,
    const Rcpp::NumericVector& stratumFraction = 1,
    const double lambda = NA_REAL,
    const double gamma = NA_REAL,
    const double phi = NA_REAL,
    const double plo = NA_REAL,
    const double T = NA_REAL,
    const int maxSubjects = NA_INTEGER,
    const Rcpp::IntegerVector& plannedSubjects = NA_INTEGER,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasets = 1,
    const int seed = 0) {

  auto p1 = Rcpp::as<std::vector<double>>(p);
  auto accrualTime1 = Rcpp::as<std::vector<double>>(accrualTime);
  auto accrualIntensity1 = Rcpp::as<std::vector<double>>(accrualIntensity);
  auto stratumFraction1 = Rcpp::as<std::vector<double>>(stratumFraction);
  auto plannedSubjects1 = Rcpp::as<std::vector<int>>(plannedSubjects);

  ListCpp result = simonBayesSimcpp(
    p1, accrualTime1, accrualIntensity1, stratumFraction1,
    lambda, gamma, phi, plo, T, maxSubjects, plannedSubjects1,
    maxNumberOfIterations, maxNumberOfRawDatasets, seed);

  return Rcpp::wrap(result);
}
