#include "survival_analysis.h"
#include "utilities.h"
#include "dataframe_list.h"
#include "thread_utils.h"

#include <Rcpp.h>
#include <boost/random.hpp>

#include <algorithm> // accumulate, any_of, max_element, min_element, none_of, sort
#include <cctype>    // isalnum, isalpha, tolower, toupper
#include <climits>   // INT_MIN
#include <cmath>     // exp, fabs, isinf, isnan, log, sqrt
#include <cstring>   // memcpy
#include <limits>    // numeric_limits
#include <numeric>   // iota, inner_product
#include <random>    // mt19937
#include <stdexcept> // invalid_argument, runtime_error
#include <string>    // string
#include <vector>    // vector

// Helper function to compute confidence interval for survival probability
std::vector<double> fsurvci(double surv, double sesurv, std::string& ct, double z) {
  double grad, hw, lower = NaN, upper = NaN;
  if (surv == 1.0 && sesurv == 0.0) {
    lower = upper = 1.0;
  } else if (surv == 0.0 && sesurv == 0.0) {
    lower = upper = 0.0;
  } else if (ct == "plain" || ct == "linear") {
    lower = std::max(surv - z * sesurv, 0.0);
    upper = std::min(surv + z * sesurv, 1.0);
  } else if (ct == "log") {
    grad = 1.0 / surv;
    hw = z * grad * sesurv;
    lower = std::exp(std::log(surv) - hw);
    upper = std::min(std::exp(std::log(surv) + hw), 1.0);
  } else if (ct == "log-log" || ct == "loglog" || ct == "cloglog") {
    grad = 1.0 / (surv * std::log(surv));
    hw = z * grad * sesurv;
    lower = std::exp(-std::exp(std::log(-std::log(surv)) - hw));
    upper = std::exp(-std::exp(std::log(-std::log(surv)) + hw));
  } else if (ct == "logit") {
    grad = 1.0 / (surv * (1.0 - surv));
    hw = z * grad * sesurv;
    lower = boost_plogis(boost_qlogis(surv) - hw);
    upper = boost_plogis(boost_qlogis(surv) + hw);
  } else if (ct == "arcsin" || ct == "asin" || ct == "asinsqrt") {
    grad = 1.0 / (2.0 * std::sqrt(surv * (1.0 - surv)));
    hw = z * grad * sesurv;
    lower = sq(std::sin(std::asin(std::sqrt(surv)) - hw));
    upper = sq(std::sin(std::asin(std::sqrt(surv)) + hw));
  } else {
    throw std::invalid_argument("Unknown confidence type: " + ct);
  }

  return {lower, upper};
}


// Compute survival quantiles and confidence intervals
DataFrameCpp survQuantilecpp(const std::vector<double>& time,
                             const std::vector<int>& event,
                             const double cilevel,
                             const std::string& transform,
                             const std::vector<double>& probs) {
  
  // Basic input checks
  if (time.size() != event.size()) {
    throw std::invalid_argument("time and event must have the same length");
  }
  int n = static_cast<int>(time.size());
  for (int i = 0; i < n; ++i) {
    if (std::isnan(time[i])) throw std::invalid_argument("time must be provided");
  }
  for (int i = 0; i < n; ++i) {
    // event is integer; check not NaN is irrelevant; but check allowed values:
    if (!(event[i] == 0 || event[i] == 1)) {
      throw std::invalid_argument("event must be 1 or 0");
    }
  }
  
  int total_events = 0;
  for (int v : event) if (v == 1) ++total_events;
  if (total_events == 0) throw std::invalid_argument("at least 1 event is needed");
  
  if (!(cilevel > 0.0 && cilevel < 1.0)) {
    throw std::invalid_argument("cilevel must lie between 0 and 1");
  }
  
  std::string ct = transform;
  std::for_each(ct.begin(), ct.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (!(ct == "linear" || ct == "plain" || ct == "log" ||
      ct == "loglog" || ct == "log-log" || ct == "cloglog" ||
      ct == "asinsqrt" || ct == "arcsin"|| ct == "asin" ||
      ct == "logit")) {
    throw std::invalid_argument("Invalid value for transform");
  }
  
  int code;
  if (ct == "plain" || ct == "linear") code = 1;
  else if (ct == "log") code = 2;
  else if (ct == "log-log" || ct == "loglog" || ct == "cloglog") code = 3;
  else if (ct == "logit") code = 4;
  else if (ct == "arcsin" || ct == "asin" || ct == "asinsqrt") code = 5;
  else throw std::invalid_argument("Unknown confidence type: " + ct);
  
  // probs checks
  for (double p : probs) {
    if (std::isnan(p)) throw std::invalid_argument("probs must be provided");
    if (!(p > 0.0 && p < 1.0)) 
      throw std::invalid_argument("Elements of probs must lie between 0 and 1");
  }

  // sort by time, and event with event in descending order
  std::vector<int> order = seqcpp(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (time[i] != time[j]) return time[i] < time[j];
    return event[i] > event[j];
  });
  
  std::vector<double> time2(n);
  std::vector<int> event2(n);
  for (int i = 0; i < n; ++i) {
    time2[i] = time[order[i]];
    event2[i] = event[order[i]];
  }
  
  std::vector<double> time0, nrisk0, nevent0, surv0, sesurv0;
  time0.reserve(total_events);
  nrisk0.reserve(total_events);
  nevent0.reserve(total_events);
  surv0.reserve(total_events);
  sesurv0.reserve(total_events);
  
  double t = 0, nrisk = n, nevent = 0, surv = 1, vcumhaz = 0, sesurv;
  
  bool cache = false; // buffer for the current event time
  for (int j = 0; j < n; ++j) {
    if ((j == 0 && event2[j] == 1) ||
        (j >= 1 && event2[j] == 1 && time2[j] > time2[j-1])) {
      // new event, add the info for the previous event
      if (cache) {
        surv *= (1.0 - nevent / nrisk);
        if (nrisk > nevent) {
          vcumhaz += nevent / (nrisk * (nrisk - nevent));
        } else {
          vcumhaz = NaN;
        }
        sesurv = surv * std::sqrt(vcumhaz);
        
        time0.push_back(t);
        nrisk0.push_back(nrisk);
        nevent0.push_back(nevent);
        surv0.push_back(surv);
        sesurv0.push_back(sesurv);
      }
      
      // update the buffer for the current event time
      t = time2[j];
      nrisk = n-j;
      nevent = 1;
      cache = true;
    } else if (j >= 1 && event2[j] == 1 && event2[j-1] == 1 && 
      time2[j] == time2[j-1]) {
      // tied event
      ++nevent;
    } else if (j >= 1 && event2[j] == 0 && event2[j-1] == 1) {
      // new censoring, add the info for the previous event
      surv *= (1.0 - nevent / nrisk);
      if (nrisk > nevent) {
        vcumhaz += nevent / (nrisk * (nrisk - nevent));
      } else {
        vcumhaz = NaN;
      }
      sesurv = surv * std::sqrt(vcumhaz);
      
      time0.push_back(t);
      nrisk0.push_back(nrisk);
      nevent0.push_back(nevent);
      surv0.push_back(surv);
      sesurv0.push_back(sesurv);
      
      // empty the buffer for the current event time
      cache = false;
    }
  }
  
  // add the info for the last event
  if (cache) {
    surv *= (1.0 - nevent / nrisk);
    if (nrisk > nevent) {
      vcumhaz += nevent / (nrisk * (nrisk - nevent));
    } else {
      vcumhaz = NaN;
    }
    sesurv = surv * std::sqrt(vcumhaz);
    
    time0.push_back(t);
    nrisk0.push_back(nrisk);
    nevent0.push_back(nevent);
    surv0.push_back(surv);
    sesurv0.push_back(sesurv);
  }
  
  int n0 = time0.size();

  std::vector<double> z(n0, NaN), grad(n0, NaN);
  double zcrit = boost_qnorm((1.0 + cilevel) / 2.0);
  
  int m = probs.size();
  std::vector<double> quantile(m), lower(m), upper(m);
  for (int j = 0; j < m; ++j) {
    double q = 1.0 - probs[j];
    for (int i = 0; i < n0; ++i) {
      if (nrisk0[i] > nevent0[i]) {
        switch (code) {
        case 1: // linear or plain
          z[i] = (surv0[i] - q) / sesurv0[i];
          break;
          
        case 2: // log
          grad[i] = 1.0 / surv0[i];
          z[i] = (std::log(surv0[i]) - std::log(q)) / (grad[i] * sesurv0[i]);
          break;
          
        case 3: // loglog / log-log / cloglog
          grad[i] = 1.0 / (surv0[i] * std::log(surv0[i]));
          z[i] = (std::log(-std::log(surv0[i])) - std::log(-std::log(q))) / 
            (grad[i] * sesurv0[i]);
          break;
          
        case 4: // logit
          grad[i] = 1.0 / (surv0[i] * (1.0 - surv0[i]));
          z[i] = (boost_qlogis(surv0[i]) - boost_qlogis(q)) / 
            (grad[i] * sesurv0[i]);
          break;
          
        case 5: // arcsin / asin / asinsqrt
          grad[i] = 1.0 / (2.0 * std::sqrt(surv0[i] * (1.0 - surv0[i])));
          z[i] = (std::asin(std::sqrt(surv0[i])) - std::asin(std::sqrt(q))) / 
            (grad[i] * sesurv0[i]);
          break;
        }
      }
    }
    
   // find indices with |z| <= zcrit (we only need min and max)
    int imin = -1, imax = -1;
    for (int i = 0; i < n0; ++i) {
      if (!std::isnan(z[i]) && std::fabs(z[i]) <= zcrit) {
        if (imin == -1) imin = i;
        imax = i;
      }
    }
    if (imin == -1) {
      lower[j] = upper[j] = NaN;
    } else {
      lower[j] = time0[imin];
      if (imax < n0 - 1) upper[j] = time0[imax + 1];
      else upper[j] = NaN;
    }
    
    // quantile: first time where surv0 < q. surv0 is nonincreasing, do binary search
    int first_lt = -1;
    int lo = 0, hi = n0 - 1;
    while (lo <= hi) {
      int mid = (lo + hi) >> 1;
      if (surv0[mid] < q) {
        first_lt = mid;
        hi = mid - 1;
      } else {
        lo = mid + 1;
      }
    }
    if (first_lt >= 0) quantile[j] = time0[first_lt];
    else quantile[j] = NaN;
  }
  
  DataFrameCpp result;
  result.push_back(probs, "prob");
  result.push_back(std::move(quantile), "quantile");
  result.push_back(std::move(lower), "lower");
  result.push_back(std::move(upper), "upper");
  result.push_back(cilevel, "cilevel");
  result.push_back(ct, "transform");
  return result;
}

//' @title Brookmeyer-Crowley Confidence Interval for Quantiles of
//' Right-Censored Time-to-Event Data
//' @description Obtains the Brookmeyer-Crowley confidence
//' interval for quantiles of right-censored time-to-event data.
//'
//' @param time The vector of possibly right-censored survival times.
//' @param event The vector of event indicators.
//' @param cilevel The confidence interval level. Defaults to 0.95.
//' @param transform The transformation of the survival function to use
//'   to construct the confidence interval. Options include
//'   "linear" (alternatively "plain"), "log",
//'   "loglog" (alternatively "log-log" or "cloglog"),
//'   "asinsqrt" (alternatively "asin" or "arcsin"), and "logit".
//'   Defaults to "loglog".
//'
//' @param probs The vector of probabilities to calculate the quantiles.
//'   Defaults to c(0.25, 0.5, 0.75).
//'
//' @return A data frame containing the estimated quantile and
//' confidence interval corresponding to each specified probability.
//' It includes the following variables:
//'
//' * \code{prob}: The probability to calculate the quantile.
//'
//' * \code{quantile}: The estimated quantile.
//'
//' * \code{lower}: The lower limit of the confidence interval.
//'
//' * \code{upper}: The upper limit of the confidence interval.
//'
//' * \code{cilevel}: The confidence interval level.
//'
//' * \code{transform}: The transformation of the survival function to use
//'   to construct the confidence interval.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' survQuantile(
//'   time = c(33.7, 3.9, 10.5, 5.4, 19.5, 23.8, 7.9, 16.9, 16.6,
//'            33.7, 17.1, 7.9, 10.5, 38),
//'   event = c(0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1),
//'   probs = c(0.25, 0.5, 0.75))
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame survQuantile(
    const Rcpp::NumericVector& time,
    const Rcpp::NumericVector& event,
    const double cilevel = 0.95,
    const std::string& transform = "loglog",
    const Rcpp::NumericVector& probs = Rcpp::NumericVector::create(0.25, 0.5, 0.75)) {
  
  auto timev = Rcpp::as<std::vector<double>>(time);
  auto eventv = Rcpp::as<std::vector<int>>(event);
  auto probsv = Rcpp::as<std::vector<double>>(probs);
  
  auto cpp_result = survQuantilecpp(timev, eventv, cilevel, transform, probsv);
  return Rcpp::wrap(cpp_result);
}


// Compute Kaplan-Meier estimator
DataFrameCpp kmestcpp(const DataFrameCpp& data,
                      const std::vector<std::string>& stratum,
                      const std::string& time,
                      const std::string& time2,
                      const std::string& event,
                      const std::string& weight,
                      const std::string& conftype,
                      const double conflev,
                      const bool keep_censor) {
  
  int n = data.nrows();
  
  bool has_stratum = false;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
  }
  
  if (!data.containElementNamed(time))
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (double v : timen) if (v < 0.0)
    throw std::invalid_argument("time must be nonnegative for each subject");
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    if (data.int_cols.count(time2)) {
      const std::vector<int>& vi = data.get<int>(time2);
      for (int i = 0; i < n; ++i) time2n[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(time2)) {
      time2n = data.get<double>(time2);
    } else {
      throw std::invalid_argument("time2 variable must be integer or numeric");
    }
    for (int i = 0; i < n; ++i) {
      if (time2n[i] <= timen[i]) {
        throw std::invalid_argument(
            "time2 must be greater than time for each observation");
      }
    }
  }
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = std::move(timen);
    // tstartn remains zero-initialized
  } else {
    tstartn = std::move(timen);
    tstopn = std::move(time2n);
  }
  
  if (!data.containElementNamed(event))
    throw std::invalid_argument("data must contain the event variable");
  std::vector<int> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
  } else if (data.int_cols.count(event)) {
    eventn = data.get<int>(event);
  } else if (data.numeric_cols.count(event)) {
    const std::vector<double>& vd = data.get<double>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1)
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const std::vector<int>& vi = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else {
      throw std::invalid_argument("weight variable must be integer or numeric");
    }
    for (double v : weightn) if (v <= 0.0)
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::string ct = conftype;
  std::for_each(ct.begin(), ct.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (!(ct == "none" || ct == "plain" || ct == "log" || ct == "log-log" || 
      ct == "logit" || ct == "arcsin")) {
    throw std::invalid_argument(
        "conftype must be none, plain, log, log-log, logit, or arcsin");
  }
  
  if (conflev <= 0.0 || conflev >= 1.0) {
    throw std::invalid_argument("conflev must lie between 0 and 1");
  }
  
  // confidence interval for survival probability
  double z = boost_qnorm((1.0 + conflev) / 2.0);
  
  // sort by stopping time in descending order within each stratum
  std::vector<int> order = seqcpp(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (stratumn[i] != stratumn[j]) return stratumn[i] > stratumn[j];
    if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
    return eventn[i] < eventn[j];
  });
  
  subset_in_place(stratumn, order);
  subset_in_place(tstartn, order);
  subset_in_place(tstopn, order);
  subset_in_place(eventn, order);
  subset_in_place(weightn, order);
  
  // sort by starting time in descending order within each stratum
  std::vector<int> order1 = seqcpp(0, n-1);
  std::sort(order1.begin(), order1.end(), [&](int i, int j) {
    if (stratumn[i] != stratumn[j]) return stratumn[i] > stratumn[j];
    return tstartn[i] > tstartn[j];
  });
  
  double nrisk = 0;   // number at risk
  double nriskw = 0;  // weighted # at risk
  double nriskw2 = 0; // weight^2 # at risk
  double nevent = 0;  // number of events
  double neventw = 0; // weighted # events
  double ncensor = 0; // number of censored
  double surv = 1;    // survival probability
  double sesurv = 0;  // standard error of survival probability
  double vcumhaz = 0; // cumulative hazard variance
  
  int istratum = stratumn.empty() ? 0 : stratumn[0]; // current stratum
  int i1 = 0;                 // index for removing out-of-risk subjects
  
  int stratum_size = 0;
  std::vector<int> stratum0, size00(n);
  std::vector<double> time0, nrisk0, nriskw0, nriskw20;
  std::vector<double> nevent0, neventw0, ncensor0;
  stratum0.reserve(n); time0.reserve(n); 
  nrisk0.reserve(n); nriskw0.reserve(n); nriskw20.reserve(n);
  nevent0.reserve(n); neventw0.reserve(n); ncensor0.reserve(n);
  
  for (int i = 0; i < n; ) {
    // Reset when entering a new stratum
    if (stratumn[i] != istratum) {
      istratum = stratumn[i];
      i1 = i;
      stratum_size = 0;
      nrisk = nriskw = nriskw2 = 0;
    }
    
    const double dtime = tstopn[i];
    
    // Process all persons tied at this dtime
    while (i < n && tstopn[i] == dtime) {
      ++nrisk;
      nriskw += weightn[i];
      nriskw2 += weightn[i] * weightn[i];
      
      if (eventn[i] == 1) {
        ++nevent;
        neventw += weightn[i];
      } else {
        ++ncensor;
      }
      
      if (tstartn[i] == 0) ++stratum_size; // unique subjects in the stratum
      
      ++i;
      
      if (i < n && stratumn[i] != istratum) {
        size00[istratum] = stratum_size; // update size of the stratum
        break;
      }
    }
    
    // remove subjects no longer at risk
    for (; i1 < n; ++i1) {
      const int p1 = order1[i1];
      if (tstartn[p1] < dtime || stratumn[p1] != istratum) break;
      
      --nrisk;
      nriskw -= weightn[p1];
      nriskw2 -= weightn[p1] * weightn[p1];
    }
    
    if (nevent > 0 || keep_censor) {
      stratum0.push_back(istratum);
      time0.push_back(dtime);
      nrisk0.push_back(nrisk);
      nriskw0.push_back(nriskw);
      nriskw20.push_back(nriskw2);
      nevent0.push_back(nevent);
      neventw0.push_back(neventw);
      ncensor0.push_back(ncensor);
      
      nevent = neventw = 0;
    }
    
    ncensor = 0;
  }
  
  
  size00[istratum] = stratum_size; // update size of the last stratum
  
  int m = stratum0.size();
  std::vector<int> size0(m);
  std::vector<double> surv0(m), sesurv0(m), lower0(m), upper0(m);
  
  if (m > 0) {
    istratum = stratum0[m - 1];
    for (int i = m - 1; i >= 0; --i) {
      double nevent_l = nevent0[i];
      double neventw_l = neventw0[i];
      double nriskw_l = nriskw0[i];
      double nriskw2_l = nriskw20[i];
      
      if (stratum0[i] != istratum) { // hit a new stratum
        // reset temporary variables
        istratum = stratum0[i];
        surv = 1.0; vcumhaz = 0.0; sesurv = 0.0;
      }
      
      if (nevent_l > 0) {
        double p = neventw_l / nriskw_l;
        surv *= (1.0 - p);
        double a = nriskw_l * nriskw_l / nriskw2_l;
        vcumhaz += p / (a * (1.0 - p));
        sesurv = surv * std::sqrt(vcumhaz);
      }
      
      size0[i] = size00[istratum];
      surv0[i] = surv;
      sesurv0[i] = sesurv;
      if (ct != "none") {
        std::vector<double> ci = fsurvci(surv, sesurv, ct, z);
        lower0[i] = ci[0];
        upper0[i] = ci[1];
      }
    }
  }
  
  // update the global m (reverse)
  std::vector<int> stratum1(m), size1(m);
  std::vector<double> time1(m), nrisk1(m), nevent1(m), ncensor1(m);
  std::vector<double> surv1(m), sesurv1(m), lower1(m), upper1(m);
  
  for (int i = 0; i < m; ++i) {
    int k = m - 1 - i;               // source index (reversed)
    stratum1[i] = stratum0[k];
    size1[i]    = size0[k];
    time1[i]    = time0[k];
    nrisk1[i]   = nrisk0[k];
    nevent1[i]  = nevent0[k];
    ncensor1[i] = ncensor0[k];
    surv1[i]    = surv0[k];
    sesurv1[i]  = sesurv0[k];
    if (ct != "none") {
      lower1[i] = lower0[k];
      upper1[i] = upper0[k];
    }
  }
  
  DataFrameCpp result;
  result.push_back(std::move(size1), "size");
  result.push_back(std::move(time1), "time");
  result.push_back(std::move(nrisk1), "nrisk");
  result.push_back(std::move(nevent1), "nevent");
  result.push_back(std::move(ncensor1), "ncensor");
  result.push_back(std::move(surv1), "surv");
  result.push_back(std::move(sesurv1), "sesurv");
  
  if (ct != "none") {
    result.push_back(std::move(lower1), "lower");
    result.push_back(std::move(upper1), "upper");
    result.push_back(conflev, "conflev");
    result.push_back(ct, "conftype");
  }
  
  // add stratum lookup columns (if any)
  if (has_stratum) {
    for (int i = 0; i < p_stratum; ++i) {
      std::string s = stratum[i];
      if (u_stratum.int_cols.count(s)) {
        auto v = u_stratum.get<int>(s);
        result.push_back(subset(v, stratum1), s);
      } else if (u_stratum.numeric_cols.count(s)) {
        auto v = u_stratum.get<double>(s);
        result.push_back(subset(v, stratum1), s);
      } else if (u_stratum.string_cols.count(s)) {
        auto v = u_stratum.get<std::string>(s);
        result.push_back(subset(v, stratum1), s);
      } else {
        throw std::invalid_argument("unsupported type for stratum variable " + s);
      }
    }
  }
  
  return result;
}


//' @title Kaplan-Meier Estimates of Survival Curve
//' @description Obtains the Kaplan-Meier estimates of the survival curve.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{time}: The follow-up time for right censored data, or
//'     the left end of each interval for counting process data.
//'
//'   * \code{time2}: The right end of each interval for counting process
//'     data. Intervals are assumed to be open on the left
//'     and closed on the right, and event indicates whether an event
//'     occurred at the right end of each interval.
//'
//'   * \code{event}: The event indicator, 1=event, 0=no event.
//'
//'   * \code{weight}: The weight for each observation.
//'
//' @param stratum The name(s) of the stratum variable(s) in the input data.
//' @param time The name of the time variable or the left end of each
//'   interval for counting process data in the input data.
//' @param time2 The name of the right end of each interval for counting
//'   process data in the input data.
//' @param event The name of the event variable in the input data.
//' @param weight The name of the weight variable in the input data.
//' @param conftype The type of the confidence interval. One of "none",
//'   "plain", "log", "log-log" (the default), or "arcsin".
//'   The arcsin option bases the intervals on asin(sqrt(survival)).
//' @param conflev The level of the two-sided confidence interval for
//'   the survival probabilities. Defaults to 0.95.
//' @param keep_censor Whether to retain the censoring time in the output
//'   data frame.
//'
//' @return A data frame with the following variables:
//'
//' * \code{size}: The number of subjects in the stratum.
//'
//' * \code{time}: The event time.
//'
//' * \code{nrisk}: The number of subjects at risk.
//'
//' * \code{nevent}: The number of subjects having the event.
//'
//' * \code{ncensor}: The number of censored subjects.
//'
//' * \code{surv}: The Kaplan-Meier estimate of the survival probability.
//'
//' * \code{sesurv}: The standard error of the estimated survival
//'   probability based on the Greendwood formula.
//'
//' * \code{lower}: The lower bound of confidence interval if requested.
//'
//' * \code{upper}: The upper bound of confidence interval if requested.
//'
//' * \code{conflev}: The level of confidence interval if requested.
//'
//' * \code{conftype}: The type of confidence interval if requested.
//'
//' * \code{stratum}: The stratum.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' kmest(data = aml, stratum = "x", time = "time", event = "status")
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame kmest(const Rcpp::DataFrame& data,
                      const Rcpp::StringVector& stratum = "",
                      const std::string& time = "time",
                      const std::string& time2 = "",
                      const std::string& event = "event",
                      const std::string& weight = "",
                      const std::string& conftype = "log-log",
                      const double conflev = 0.95,
                      const bool keep_censor = false) {
  
  auto dfcpp = convertRDataFrameToCpp(data);
  auto stratumcpp = Rcpp::as<std::vector<std::string>>(stratum);
  
  auto cpp_result = kmestcpp(
    dfcpp, stratumcpp, time, time2, event, weight, 
    conftype, conflev, keep_censor
  );
  
  return Rcpp::wrap(cpp_result);
}


// Compute difference in KM estimates at a milestone time
DataFrameCpp kmdiffcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat,
                       const std::string& time,
                       const std::string& time2,
                       const std::string& event,
                       const std::string& weight,
                       const double milestone,
                       const double survDiffH0,
                       const double conflev) {
  int n = data.nrows();

  bool has_stratum = false;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
  }
  
  if (!data.containElementNamed(treat)) {
    throw std::invalid_argument("data must contain the treat variable");
  }
  
  // create the numeric treat variable
  std::vector<int> treatn(n);
  std::vector<int> treatwi;
  std::vector<double> treatwn;
  std::vector<std::string> treatwc;
  if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
    std::vector<int> treatv(n);
    if (data.bool_cols.count(treat)) {
      const std::vector<unsigned char>& treatvb = data.get<unsigned char>(treat);
      for (int i = 0; i < n; ++i) treatv[i] = treatvb[i] ? 1 : 0;
    } else treatv = data.get<int>(treat);
    treatwi = unique_sorted(treatv); // obtain unique treatment values
    if (treatwi.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwi.begin(), treatwi.end(), [](int v) { 
      return v == 0 || v == 1; })) {
      treatwi = {1, 0}; // special handling for 1 / 0 treatment coding
      for (int i = 0; i < n; ++i) treatn[i] = 2 - treatv[i];
    } else {
      treatn = matchcpp(treatv, treatwi, 1);
    }
  } else if (data.numeric_cols.count(treat)) {
    const std::vector<double>& treatv = data.get<double>(treat);
    treatwn = unique_sorted(treatv);
    if (treatwn.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwn.begin(), treatwn.end(), [](double v) { 
      return v == 0.0 || v == 1.0; })) {
      treatwn = {1.0, 0.0};
      for (int i = 0; i < n; ++i) treatn[i] = 2 - static_cast<int>(treatv[i]);
    } else {
      treatn = matchcpp(treatv, treatwn, 1);
    }
  } else if (data.string_cols.count(treat)) {
    const std::vector<std::string>& treatv = data.get<std::string>(treat);
    treatwc = unique_sorted(treatv);
    if (treatwc.size() != 2) 
      throw std::invalid_argument("treat must have two and only two distinct values");
    treatn = matchcpp(treatv, treatwc, 1);
  } else {
    throw std::invalid_argument(
        "incorrect type for the treat variable in the input data");
  }
  
  if (!data.containElementNamed(time)) {
    throw std::invalid_argument("data must contain the time variable");
  }
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (double v : timen) if (v < 0.0)
    throw std::invalid_argument("time must be nonnegative for each subject");
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    if (data.int_cols.count(time2)) {
      const std::vector<int>& vi = data.get<int>(time2);
      for (int i = 0; i < n; ++i) time2n[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(time2)) {
      time2n = data.get<double>(time2);
    } else {
      throw std::invalid_argument("time2 variable must be integer or numeric");
    }
    for (int i = 0; i < n; ++i) {
      if (time2n[i] <= timen[i]) {
        throw std::invalid_argument(
            "time2 must be greater than time for each observation");
      }
    }
  }
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = std::move(timen);
    // tstartn remains zero-initialized
  } else {
    tstartn = std::move(timen);
    tstopn = std::move(time2n);
  }
  
  if (!data.containElementNamed(event))
    throw std::invalid_argument("data must contain the event variable");
  std::vector<int> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
  } else if (data.int_cols.count(event)) {
    eventn = data.get<int>(event);
  } else if (data.numeric_cols.count(event)) {
    const std::vector<double>& vd = data.get<double>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1)
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const std::vector<int>& vi = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else {
      throw std::invalid_argument("weight variable must be integer or numeric");
    }
    for (double v : weightn) if (v <= 0.0)
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  if (std::isnan(milestone)) {
    throw std::invalid_argument("milestone must be provided");
  }
  
  if (milestone <= 0) {
    throw std::invalid_argument("milestone must be positive");
  }
  
  if (survDiffH0 <= -1 || survDiffH0 >= 1) {
    throw std::invalid_argument("survDiffH0 must lie between -1 and 1");
  }
  
  if (conflev <= 0 || conflev >= 1) {
    throw std::invalid_argument("conflev must lie between 0 and 1");
  }
  
  bool noerr = true;
  
  // sort by stratum
  std::vector<int> order = seqcpp(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return stratumn[i] < stratumn[j];
  });
  
  subset_in_place(stratumn, order);
  subset_in_place(treatn, order);
  subset_in_place(tstartn, order);
  subset_in_place(tstopn, order);
  subset_in_place(eventn, order);
  subset_in_place(weightn, order);
  
  // identify the locations of the unique values of stratum
  std::vector<int> idx(1,0);
  if (has_stratum) {
    for (int i = 1; i < n; ++i) {
      if (stratumn[i] != stratumn[i-1]) {
        idx.push_back(i);
      }
    }
  }

  int nstrata = idx.size();
  idx.push_back(n);
  
  // whether the milestone exceeds the largest observed time
  for (int i = 0; i < nstrata; ++i) {
    int start = idx[i], end = idx[i+1];
    int n1 = end - start;
    std::vector<int> treat1 = subset(treatn, start, end);
    std::vector<double> tstop1 = subset(tstopn, start, end);
    std::vector<double> tstop11, tstop12;
    for (int i = 0; i < n1; ++i) {
      if (treat1[i] == 1) {
        tstop11.push_back(tstop1[i]);
      } else {
        tstop12.push_back(tstop1[i]);
      }
    }
    
    double max_time11 = *std::max_element(tstop11.begin(), tstop11.end());
    double max_time12 = *std::max_element(tstop12.begin(), tstop12.end());
    if (milestone > std::min(max_time11, max_time12)) {
      std::string stratumerr;
      if (!has_stratum) {
        stratumerr = "";
      } else {
        for (int j = 0; j < p_stratum; ++j) {
          std::string s = stratum[j];
          if (u_stratum.int_cols.count(s)) {
            auto v = u_stratum.get<int>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.numeric_cols.count(s)) {
            auto v = u_stratum.get<double>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.string_cols.count(s)) {
            auto v = u_stratum.get<std::string>(s);
            stratumerr += " " + s + " = " + v[i];
          } else {
            throw std::invalid_argument("unsupported type for stratum variable " + s);
          }
        }
      }
      
      int k = milestone > max_time11 ? 0 : 1;
      std::string treaterr;
      if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
        treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
      } else if (data.numeric_cols.count(treat)) {
        treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
      } else if (data.string_cols.count(treat)) {
        treaterr = " " + treat + " = " + treatwc[k];
      } else {
        throw std::invalid_argument("unsupported type for treat variable " + treat);
      }
      
      std::string str1 = "The milestone is larger than the largest observed time for";
      std::string errmsg = str1 + treaterr;
      if (!stratumerr.empty()) errmsg = errmsg + " " + stratumerr;
      
      if (noerr) {
        thread_utils::push_thread_warning(
          errmsg + "\nAdditional warning messages are suppressed.");
        noerr = false;
      }
      
      continue;
    }
  }

  DataFrameCpp dfin;
  dfin.push_back(std::move(stratumn), "stratum");
  dfin.push_back(std::move(treatn), "treat");
  dfin.push_back(std::move(tstartn), "tstart");
  dfin.push_back(std::move(tstopn), "tstop");
  dfin.push_back(std::move(eventn), "event");
  dfin.push_back(std::move(weightn), "weight");
  
  DataFrameCpp dfout = kmestcpp(dfin, {"stratum", "treat"}, "tstart", "tstop", 
                                "event", "weight", "none", 0.95, false);
  
  std::vector<int> stratum2 = dfout.get<int>("stratum");
  std::vector<int> treat2 = dfout.get<int>("treat");
  std::vector<int> treatsize = dfout.get<int>("size");
  std::vector<double> time20 = dfout.get<double>("time");
  std::vector<double> survival2 = dfout.get<double>("surv");
  std::vector<double> stderr2 = dfout.get<double>("sesurv");
  
  int n2 = stratum2.size();
  
  // identify the locations of the unique values of stratum
  std::vector<int> idx2(1,0);
  if (has_stratum) {
    for (int i = 1; i < n2; ++i) {
      if (stratum2[i] != stratum2[i-1]) {
        idx2.push_back(i);
      }
    } 
  }
  
  nstrata = idx2.size();
  idx2.push_back(n2);
  
  std::vector<double> m(nstrata); // number of subjects in each stratum
  for (int i = 0; i < nstrata; ++i) {
    int j1 = idx2[i], j2 = idx2[i+1] - 1;
    if (treat2[j1] != 1 || treat2[j2] != 2) {
      std::string stratumerr;
      if (!has_stratum) {
        stratumerr = "";
      } else {
        for (int j = 0; j < p_stratum; ++j) {
          std::string s = stratum[j];
          if (u_stratum.int_cols.count(s)) {
            auto v = u_stratum.get<int>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.numeric_cols.count(s)) {
            auto v = u_stratum.get<double>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.string_cols.count(s)) {
            auto v = u_stratum.get<std::string>(s);
            stratumerr += " " + s + " = " + v[i];
          } else {
            throw std::invalid_argument("unsupported type for stratum variable " + s);
          }
        }
      }
      
      int k = treat2[j1] != 1 ? 0 : 1;
      std::string treaterr;
      if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
        treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
      } else if (data.numeric_cols.count(treat)) {
        treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
      } else if (data.string_cols.count(treat)) {
        treaterr = " " + treat + " = " + treatwc[k];
      } else {
        throw std::invalid_argument("unsupported type for treat variable " + treat);
      }
      
      std::string str1 = "The data set does not contain";
      std::string errmsg = str1 + treaterr;
      if (!stratumerr.empty()) errmsg = errmsg + " " + stratumerr;
      
      if (noerr) {
        thread_utils::push_thread_warning(
          errmsg + "\nAdditional warning messages are suppressed.");
        noerr = false;
      }
      
      continue;
    }
    
    m[i] = treatsize[j1] + treatsize[j2];
  }
  
  double M = std::accumulate(m.begin(), m.end(), 0.0);
  std::vector<double> p(nstrata);
  
  double surv1 = 0.0, surv2 = 0.0, vsurv1 = 0.0, vsurv2 = 0.0;
  for (int i = 0; i < nstrata; ++i) {
    p[i] = m[i] / M; // fraction of subjects in the stratum
    int start = idx2[i], end = idx2[i+1];
    int nx = end - start;
    std::vector<int> treatx = subset(treat2, start, end);
    std::vector<double> timex = subset(time20, start, end);
    std::vector<double> survivalx = subset(survival2, start, end);
    std::vector<double> stderrx = subset(stderr2, start, end);
    
    std::vector<double> surv(2), vsurv(2);
    for (int j = 0; j < 2; ++j) {
      
      std::vector<double> time0, survival0, stderr0;
      for (int k = 0; k < nx; ++k) {
        if (treatx[k] == j+1) {
          time0.push_back(timex[k]);
          survival0.push_back(survivalx[k]);
          stderr0.push_back(stderrx[k]);
        }
      }
      int K = time0.size();
      
      // find the latest event time before milestone for each treat
      int k = 0;
      for (; k < K; ++k) {
        if (time0[k] > milestone) break;
      }
      
      if (k == 0) {
        surv[j] = 1;
        vsurv[j] = 0;
      } else {
        --k;
        surv[j] = survival0[k];
        vsurv[j] = stderr0[k] * stderr0[k];
      }
    }
    
    surv1 += p[i] * surv[0];
    surv2 += p[i] * surv[1];
    vsurv1 += p[i] * p[i] * vsurv[0];
    vsurv2 += p[i] * p[i] * vsurv[1];
  }
  
  double z = boost_qnorm((1.0 + conflev) / 2.0);
  
  double survDiff = surv1 - surv2;
  double sesurvDiff = std::sqrt(vsurv1 + vsurv2);
  double survDiffZ = (survDiff - survDiffH0) / sesurvDiff;
  double survDiffPValue = 2.0 * boost_pnorm(-std::fabs(survDiffZ));
  double lower = survDiff - z * sesurvDiff;
  double upper = survDiff + z * sesurvDiff;
  
  DataFrameCpp result;
  result.push_back(milestone, "milestone");
  result.push_back(survDiffH0, "survDiffH0");
  result.push_back(surv1, "surv1");
  result.push_back(surv2, "surv2");
  result.push_back(survDiff, "survDiff");
  result.push_back(vsurv1, "vsurv1");
  result.push_back(vsurv2, "vsurv2");
  result.push_back(sesurvDiff, "sesurvDiff");
  result.push_back(survDiffZ, "survDiffZ");
  result.push_back(survDiffPValue, "survDiffPValue");
  result.push_back(lower, "lower");
  result.push_back(upper, "upper");
  result.push_back(conflev, "conflev");
  
  return result;
}



//' @title Estimate of Milestone Survival Difference
//' @description Obtains the estimate of milestone survival difference
//' between two treatment groups.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{treat}: The treatment.
//'
//'   * \code{time}: The follow-up time for right censored data, or
//'     the left end of each interval for counting process data.
//'
//'   * \code{time2}: The right end of each interval for counting process
//'     data. Intervals are assumed to be open on the left
//'     and closed on the right, and event indicates whether an event
//'     occurred at the right end of each interval.
//'
//'   * \code{event}: The event indicator, 1=event, 0=no event.
//'
//'   * \code{weight}: The weight for each observation.
//'
//' @param stratum The name of the stratum variable in the input data.
//' @param treat The name of the treatment variable in the input data.
//' @param time The name of the time variable or the left end of each
//'   interval for counting process data in the input data.
//' @param time2 The name of the right end of each interval for counting
//'   process data in the input data.
//' @param event The name of the event variable in the input data.
//' @param weight The name of the weight variable in the input data.
//' @param milestone The milestone time at which to calculate the
//'   survival probability.
//' @param survDiffH0 The difference in milestone survival probabilities
//'   under the null hypothesis. Defaults to 0 for superiority test.
//' @param conflev The level of the two-sided confidence interval for
//'   the difference in milestone survival probabilities. Defaults to 0.95.
//'
//' @return A data frame with the following variables:
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{survDiffH0}: The difference in milestone survival probabilities
//'   under the null hypothesis.
//'
//' * \code{surv1}: The estimated milestone survival probability for
//'   the treatment group.
//'
//' * \code{surv2}: The estimated milestone survival probability for
//'   the control group.
//'
//' * \code{survDiff}: The estimated difference in milestone survival
//'   probabilities.
//'
//' * \code{vsurv1}: The variance for surv1.
//'
//' * \code{vsurv2}: The variance for surv2.
//'
//' * \code{sesurvDiff}: The standard error for survDiff.
//'
//' * \code{survDiffZ}: The Z-statistic value.
//'
//' * \code{survDiffPValue}: The two-sided p-value.
//'
//' * \code{lower}: The lower bound of confidence interval.
//'
//' * \code{upper}: The upper bound of confidence interval.
//'
//' * \code{conflev}: The level of confidence interval.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' kmdiff(data = rawdata[rawdata$iterationNumber == 1, ],
//'        stratum = "stratum", treat = "treatmentGroup",
//'        time = "timeUnderObservation", event = "event",
//'        milestone = 12)
//'              
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame kmdiff(const Rcpp::DataFrame& data,
                       const Rcpp::StringVector& stratum = "",
                       const std::string& treat = "treat",
                       const std::string& time = "time",
                       const std::string& time2 = "",
                       const std::string& event = "event",
                       const std::string& weight = "",
                       const double milestone = 0,
                       const double survDiffH0 = 0,
                       const double conflev = 0.95) {
  
  auto dfcpp = convertRDataFrameToCpp(data);
  auto stratumcpp = Rcpp::as<std::vector<std::string>>(stratum);
  
  auto cpp_result = kmdiffcpp(
    dfcpp, stratumcpp, treat, time, time2, event, weight, 
    milestone, survDiffH0, conflev
  );
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(cpp_result);
}


// Compute log-rank test statistic
DataFrameCpp lrtestcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat,
                       const std::string& time,
                       const std::string& time2,
                       const std::string& event,
                       const std::string& weight,
                       const bool weight_readj,
                       const double rho1,
                       const double rho2) {
  int n = data.nrows();
  
  std::vector<int> stratumn(n);
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    stratumn = out.get<std::vector<int>>("index");
  }
  
  // ---- treatment -> integer coding 1 / 2 ----
  if (!data.containElementNamed(treat)) {
    throw std::invalid_argument("data must contain the treat variable");
  }
  std::vector<int> treatn(n);
  if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
    std::vector<int> treatv(n);
    if (data.bool_cols.count(treat)) {
      const auto& vb = data.get<unsigned char>(treat);
      for (int i = 0; i < n; ++i) treatv[i] = vb[i] ? 1 : 0;
    } else treatv = data.get<int>(treat);
    auto treatwi = unique_sorted(treatv); // obtain unique treatment values
    if (treatwi.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwi.begin(), treatwi.end(), [](int v) { 
      return v == 0 || v == 1; })) {
      // special handling for 0 / 1 => map to 2,1 (so treated = 1, control = 2)
      for (int i = 0; i < n; ++i) treatn[i] = 2 - treatv[i];
    } else {
      treatn = matchcpp(treatv, treatwi, 1);
    }
  } else if (data.numeric_cols.count(treat)) {
    const auto& tv = data.get<double>(treat);
    auto treatwn = unique_sorted(tv);
    if (treatwn.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwn.begin(), treatwn.end(), [](double v) { 
      return v == 0.0 || v == 1.0; })) {
      for (int i = 0; i < n; ++i) treatn[i] = 2 - static_cast<int>(tv[i]);
    } else {
      treatn = matchcpp(tv, treatwn, 1);
    }
  } else if (data.string_cols.count(treat)) {
    const auto& tv = data.get<std::string>(treat);
    auto treatwc = unique_sorted(tv);
    if (treatwc.size() != 2) 
      throw std::invalid_argument("treat must have two and only two distinct values");
    treatn = matchcpp(tv, treatwc, 1);
  } else {
    throw std::invalid_argument(
        "incorrect type for the treat variable in the input data");
  }
  
  // ---- time / time2 ----
  if (!data.containElementNamed(time)) {
    throw std::invalid_argument("data must contain the time variable");
  }
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const auto& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (double v : timen) if (v < 0.0)
    throw std::invalid_argument("time must be nonnegative for each subject");
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    if (data.int_cols.count(time2)) {
      const auto& vi = data.get<int>(time2);
      for (int i = 0; i < n; ++i) time2n[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(time2)) {
      time2n = data.get<double>(time2);
    } else {
      throw std::invalid_argument("time2 variable must be integer or numeric");
    }
    for (int i = 0; i < n; ++i) {
      if (time2n[i] <= timen[i]) {
        throw std::invalid_argument(
            "time2 must be greater than time for each observation");
      }
    }
  }
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = std::move(timen);
    // tstartn remains zero-initialized
  } else {
    tstartn = std::move(timen);
    tstopn = std::move(time2n);
  }
  
  if (!data.containElementNamed(event))
    throw std::invalid_argument("data must contain the event variable");
  std::vector<int> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
  } else if (data.int_cols.count(event)) {
    eventn = data.get<int>(event);
  } else if (data.numeric_cols.count(event)) {
    const std::vector<double>& vd = data.get<double>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1)
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const auto& vi = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else {
      throw std::invalid_argument("weight variable must be integer or numeric");
    }
    for (double v : weightn) if (v <= 0.0)
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  if (rho1 < 0) throw std::invalid_argument("rho1 must be non-negative");
  if (rho2 < 0) throw std::invalid_argument("rho2 must be non-negative");
  
  // sort by stopping time in descending order within each stratum
  std::vector<int> order = seqcpp(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
    if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
    return eventn[i] < eventn[j];
  });
  
  subset_in_place(stratumn, order);
  subset_in_place(treatn, order);
  subset_in_place(tstartn, order);
  subset_in_place(tstopn, order);
  subset_in_place(eventn, order);
  subset_in_place(weightn, order);
  
  // sort by starting time in descending order within each stratum
  std::vector<int> order1 = seqcpp(0, n-1);
  std::sort(order1.begin(), order1.end(), [&](int i, int j) {
    if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
    return tstartn[i] > tstartn[j];
  });
  
  double u = 0; // score
  double v = 0; // information (variance for score)
  
  double nrisk = 0, nrisk_1 = 0, nrisk_2 = 0;    // number at risk
  double nriskw = 0, nriskw_1 = 0, nriskw_2 = 0; // weighted # at risk
  double nriskw2_1 = 0, nriskw2_2 = 0;           // weight^2 # at risk
  double nevent = 0;                             // number of events
  double neventw = 0, neventw_1 = 0;             // weighted # events

  int istratum = stratumn[0]; // current stratum
  int i1 = 0;                 // index for removing out-of-risk subjects
  
  if (rho1 == 0 && rho2 == 0) {
    for (int i = 0; i < n; ) {
      // Reset when entering a new stratum
      if (stratumn[i] != istratum) {
        istratum = stratumn[i];
        i1 = i;
        
        nrisk = nrisk_1 = nrisk_2 = 0;
        nriskw = nriskw_1 = nriskw_2 = 0;
        nriskw2_1 = nriskw2_2 = 0;
      }
      
      const double dtime = tstopn[i];
      
      // Process all persons tied at this dtime
      for (; i < n && tstopn[i] == dtime && stratumn[i] == istratum; ++i) {
        ++nrisk;
        double w = weightn[i];
        nriskw += w;
        if (treatn[i] == 1) {
          ++nrisk_1; nriskw_1 += w; nriskw2_1 += w * w;
        } else {
          ++nrisk_2; nriskw_2 += w; nriskw2_2 += w * w;
        }
        
        if (eventn[i] == 1) {
          ++nevent; neventw += w;
          if (treatn[i] == 1) neventw_1 += w;
        }
      }
      
      // Remove subjects leaving risk set
      for (; i1 < n; ++i1) {
        const int p1 = order1[i1];
        if (tstartn[p1] < dtime || stratumn[p1] != istratum) break;
        
        --nrisk;
        double w = weightn[p1]; 
        nriskw -= w;
        if (treatn[p1] == 1) {
          --nrisk_1; nriskw_1 -= w; nriskw2_1 -= w * w;
        } else {
          --nrisk_2; nriskw_2 -= w; nriskw2_2 -= w * w;
        }
      }
      
      // Add contributions for deaths at this time
      if (nevent > 0) {
        if (nrisk_1 > 0 && nrisk_2 > 0) {
          if (!weight_readj) {
            u += neventw_1 - neventw * (nriskw_1 / nriskw);
            double v1 = sq(nriskw_1 / nriskw) * nriskw2_2;
            double v2 = sq(nriskw_2 / nriskw) * nriskw2_1;
            v += nevent * (nrisk - nevent) / (nrisk * (nrisk - 1)) * (v1 + v2);
          } else {
            double neventw_2 = neventw - neventw_1;
            double neventwa_1 = neventw_1 * nrisk_1 / nriskw_1;
            double neventwa_2 = neventw_2 * nrisk_2 / nriskw_2;
            double neventwa = neventwa_1 + neventwa_2;
            u += neventwa_1 - neventwa * (nrisk_1 / nrisk);
            double v1 = sq(nrisk_1 / nrisk) * nriskw2_2 * sq(nrisk_2 / nriskw_2);
            double v2 = sq(nrisk_2 / nrisk) * nriskw2_1 * sq(nrisk_1 / nriskw_1);
            v += nevent * (nrisk - nevent) / (nrisk * (nrisk - 1)) * (v1 + v2);
          }
        }
        
        // Reset after processing deaths
        nevent = neventw = neventw_1 = 0;
      }
    }
  } else { // need to collect summaries for each event time and then iterate backward
    std::vector<int> stratum0;
    std::vector<double> time0;
    std::vector<double> nrisk0, nrisk_10, nrisk_20;
    std::vector<double> nriskw0, nriskw_10, nriskw_20;
    std::vector<double> nriskw2_10, nriskw2_20;
    std::vector<double> nevent0, neventw0, neventw_10;
    int d = n / 4; // upper bound on number of unique event times
    stratum0.reserve(d); time0.reserve(d);
    nrisk0.reserve(d); nrisk_10.reserve(d); nrisk_20.reserve(d);
    nriskw0.reserve(d); nriskw_10.reserve(d); nriskw_20.reserve(d);
    nriskw2_10.reserve(d); nriskw2_20.reserve(d);
    nevent0.reserve(d); neventw0.reserve(d); neventw_10.reserve(d);
    
    for (int i = 0; i < n; ) {
      // Reset when entering a new stratum
      if (stratumn[i] != istratum) {
        istratum = stratumn[i];
        i1 = i;
        
        nrisk = nrisk_1 = nrisk_2 = 0;
        nriskw = nriskw_1 = nriskw_2 = 0;
        nriskw2_1 = nriskw2_2 = 0;
      }
      
      const double dtime = tstopn[i];
      
      // Process all persons tied at this dtime
      for (; i < n && tstopn[i] == dtime && stratumn[i] == istratum; ++i) {
        ++nrisk;
        double w = weightn[i]; 
        nriskw += w;
        if (treatn[i] == 1) {
          ++nrisk_1; nriskw_1 += w; nriskw2_1 += w * w;
        } else {
          ++nrisk_2; nriskw_2 += w; nriskw2_2 += w * w;
        }
        
        if (eventn[i] == 1) {
          ++nevent; neventw += w;
          if (treatn[i] == 1) neventw_1 += w;
        }
      }
      
      // Remove subjects leaving risk set
      for (; i1 < n; ++i1) {
        const int p1 = order1[i1];
        if (tstartn[p1] < dtime || stratumn[p1] != istratum) break;
        
        --nrisk;
        double w = weightn[p1]; 
        nriskw -= w;
        if (treatn[p1] == 1) {
          --nrisk_1; nriskw_1 -= w; nriskw2_1 -= w * w;
        } else {
          --nrisk_2; nriskw_2 -= w; nriskw2_2 -= w * w;
        }
      }
      
      // add to the temporary storage
      if (nevent > 0) {
        stratum0.push_back(istratum);
        time0.push_back(dtime);
        nrisk0.push_back(nrisk);
        nrisk_10.push_back(nrisk_1);
        nrisk_20.push_back(nrisk_2);
        nriskw0.push_back(nriskw);
        nriskw_10.push_back(nriskw_1);
        nriskw_20.push_back(nriskw_2);
        nriskw2_10.push_back(nriskw2_1);
        nriskw2_20.push_back(nriskw2_2);
        nevent0.push_back(nevent);
        neventw0.push_back(neventw);
        neventw_10.push_back(neventw_1);

        // Reset after processing deaths
        nevent = neventw = neventw_1 = 0;
      }
    }
    
    int m = time0.size();
    if (m > 0) {
      double surv = 1.0;
      istratum = stratum0[m - 1];
      for (int i = m-1; i >= 0; --i) {
        nrisk = nrisk0[i];
        nrisk_1 = nrisk_10[i];
        nrisk_2 = nrisk_20[i];
        nriskw = nriskw0[i];
        nriskw_1 = nriskw_10[i];
        nriskw_2 = nriskw_20[i];
        nriskw2_1 = nriskw2_10[i];
        nriskw2_2 = nriskw2_20[i];
        nevent = nevent0[i];
        neventw = neventw0[i];
        neventw_1 = neventw_10[i];
        
        // Reset when entering a new stratum
        if (stratum0[i] != istratum) {
          istratum = stratum0[i];
          surv = 1;
        }
        
        // add to the main terms with S(t-)^rho1 * (1 - S(t-))^rho2 weights
        double w = std::pow(surv, rho1) * std::pow(1 - surv, rho2);
        
        double neventwa = 0, neventwa_1 = 0, neventwa_2 = 0;
        if (weight_readj) {
          double neventw_2 = neventw - neventw_1;
          neventwa_1 = neventw_1 * nrisk_1 / nriskw_1;
          neventwa_2 = neventw_2 * nrisk_2 / nriskw_2;
          neventwa = neventwa_1 + neventwa_2;
        }
        
        if (nrisk_1 > 0 && nrisk_2 > 0) {
          if (!weight_readj) {
            u += w * (neventw_1 - neventw * (nriskw_1 / nriskw));
            double v1 = sq(nriskw_1 / nriskw) * nriskw2_2;
            double v2 = sq(nriskw_2 / nriskw) * nriskw2_1;
            v += w * w * nevent * (nrisk - nevent) / (nrisk * (nrisk - 1)) * (v1 + v2);
          } else {
            u += w * (neventwa_1 - neventwa * (nrisk_1 / nrisk));
            double v1 = sq(nrisk_1 / nrisk) * nriskw2_2 * sq(nrisk_2 / nriskw_2);
            double v2 = sq(nrisk_2 / nrisk) * nriskw2_1 * sq(nrisk_1 / nriskw_1);
            v += w * w * nevent * (nrisk - nevent) / (nrisk * (nrisk - 1)) * (v1 + v2);
          }
        }
        
        // update survival probability
        if (!weight_readj) surv *= (1.0 - neventw / nriskw);
        else surv *= (1.0 - neventwa / nrisk);
      }
    }
  }
  
  double z, p;
  if (v <= 0.0) {
    z = NaN; p = NaN;
  } else {
    z = u / std::sqrt(v);
    p = 2.0 * boost_pnorm(-std::fabs(z));
  }
  
  DataFrameCpp result;
  result.push_back(u, "uscore");
  result.push_back(v, "vscore");
  result.push_back(z, "logRankZ");
  result.push_back(p, "logRankPValue");
  result.push_back(weight_readj, "weight_readj");
  result.push_back(rho1, "rho1");
  result.push_back(rho2, "rho2");
  
  return result;
}


//' @title Log-Rank Test of Survival Curve Difference
//' @description Obtains the log-rank test using the Fleming-Harrington
//' family of weights.
//'
//' @param data The input data frame or list of data frames that contains 
//' the following variables:
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{treat}: The treatment.
//'
//'   * \code{time}: The follow-up time for right censored data, or
//'     the left end of each interval for counting process data.
//'
//'   * \code{time2}: The right end of each interval for counting process
//'     data. Intervals are assumed to be open on the left
//'     and closed on the right, and event indicates whether an event
//'     occurred at the right end of each interval.
//'
//'   * \code{event}: The event indicator, 1=event, 0=no event.
//'
//'   * \code{weight}: The weight for each observation.
//'
//' @param stratum The name(s) of the stratum variable(s) in the input data.
//' @param treat The name of the treatment variable in the input data.
//' @param time The name of the time variable or the left end of each
//'   interval for counting process data in the input data.
//' @param time2 The name of the right end of each interval for counting
//'   process data in the input data.
//' @param event The name of the event variable in the input data.
//' @param weight The name of the weight variable in the input data.
//' @param weight_readj Whether the weight variable at each event time
//'   will be readjusted to be proportional to the number at risk by
//'   treatment group. Defaults to `FALSE`.
//' @param rho1 The first parameter of the Fleming-Harrington family of
//'   weighted log-rank test. Defaults to 0 for conventional log-rank test.
//' @param rho2 The second parameter of the Fleming-Harrington family of
//'   weighted log-rank test. Defaults to 0 for conventional log-rank test.
//'   
//' @return A data frame with the following variables:
//'
//' * \code{uscore}: The numerator of the log-rank test statistic.
//'
//' * \code{vscore}: The variance of the log-rank score test statistic.
//'
//' * \code{logRankZ}: The Z-statistic value.
//'
//' * \code{logRankPValue}: The two-sided p-value.
//'
//' * \code{weight_readj}: Whether the weight variable will be readjusted.
//'
//' * \code{rho1}: The first parameter of the Fleming-Harrington weights.
//'
//' * \code{rho2}: The second parameter of the Fleming-Harrington weights.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' lrtest(rawdata[rawdata$iterationNumber == 1, ], 
//'        stratum = "stratum", treat = "treatmentGroup", 
//'        time = "timeUnderObservation", event = "event", 
//'        rho1 = 0.5, rho2 = 0)
//'              
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame lrtest(const Rcpp::DataFrame data,
                       const Rcpp::StringVector& stratum = "",
                       const std::string& treat = "treat",
                       const std::string& time = "time",
                       const std::string& time2 = "",
                       const std::string& event = "event",
                       const std::string& weight = "",
                       const bool weight_readj = false,
                       const double rho1 = 0,
                       const double rho2 = 0) {
  
  auto dfcpp = convertRDataFrameToCpp(data);
  auto stratumcpp = Rcpp::as<std::vector<std::string>>(stratum);
  
  auto cpp_result = lrtestcpp(
    dfcpp, stratumcpp, treat, time, time2, event, weight, 
    weight_readj, rho1, rho2
  );
  
  return Rcpp::wrap(cpp_result);
}


// Compute Restricted Mean Survival Time
DataFrameCpp rmestcpp(const DataFrameCpp& data,
                      const std::vector<std::string>& stratum,
                      const std::string& time,
                      const std::string& event,
                      const double milestone,
                      const double conflev,
                      const bool biascorrection) {
  int n = data.nrows();
  
  bool has_stratum = false;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
  }
  
  if (!data.containElementNamed(time))
    throw std::invalid_argument("data must contain the time variable");
  
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (double v : timen) if (v < 0.0)
    throw std::invalid_argument("time must be nonnegative for each subject");
  
  if (!data.containElementNamed(event))
    throw std::invalid_argument("data must contain the event variable");
  std::vector<int> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
  } else if (data.int_cols.count(event)) {
    eventn = data.get<int>(event);
  } else if (data.numeric_cols.count(event)) {
    const std::vector<double>& vd = data.get<double>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1)
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  if (std::isnan(milestone)) {
    throw std::invalid_argument("milestone must be provided");
  }
  
  if (milestone <= 0) {
    throw std::invalid_argument("milestone must be positive");
  }
  
  if (conflev <= 0.0 || conflev >= 1.0) {
    throw std::invalid_argument("conflev must lie between 0 and 1");
  }
  
  double z = boost_qnorm((1.0 + conflev) / 2.0);
  
  std::vector<int> stratum0, size0;
  std::vector<double> rmst0, stderr0, lower0, upper0;
  stratum0.reserve(n); size0.reserve(n);
  rmst0.reserve(n); stderr0.reserve(n); lower0.reserve(n); upper0.reserve(n);
  
  bool noerr = true;
  
  // sort by stratum, time, and event with event in descending order
  std::vector<int> order = seqcpp(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
    if (timen[i] != timen[j]) return timen[i] < timen[j];
    return eventn[i] > eventn[j];
  });
  
  subset_in_place(stratumn, order);
  subset_in_place(timen, order);
  subset_in_place(eventn, order);
  
  // identify the locations of the unique values of stratum
  std::vector<int> idx1(1,0);
  if (has_stratum) {
    for (int i = 1; i < n; ++i) {
      if (stratumn[i] != stratumn[i-1]) {
        idx1.push_back(i);
      }
    }
  }
  
  int nstrata = idx1.size();
  idx1.push_back(n);
  
  for (int i = 0; i < nstrata; ++i) {
    int start = idx1[i], end = idx1[i+1];
    int n2 = end - start;
    std::vector<double> time2 = subset(timen, start, end);
    std::vector<int> event2 = subset(eventn, start, end);
    
    double max_time = *std::max_element(time2.begin(), time2.end());
    if (milestone > max_time) {
      std::string stratumerr;
      if (!has_stratum) {
        stratumerr = "";
      } else {
        for (int j = 0; j < p_stratum; ++j) {
          std::string s = stratum[j];
          if (u_stratum.int_cols.count(s)) {
            auto v = u_stratum.get<int>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.numeric_cols.count(s)) {
            auto v = u_stratum.get<double>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.string_cols.count(s)) {
            auto v = u_stratum.get<std::string>(s);
            stratumerr += " " + s + " = " + v[i];
          } else {
            throw std::invalid_argument("unsupported type for stratum variable " + s);
          }
        }
      }
      
      std::string str1 = "The milestone is larger than the largest observed time";
      std::string errmsg = str1;
      if (!stratumerr.empty()) {
        errmsg = errmsg + " " + stratumerr;
      }
      
      if (noerr) {
        thread_utils::push_thread_warning(
          errmsg + "\nAdditional warning messages are suppressed.");
        noerr = false;
      }
      
      continue;
    }
    
    std::vector<double> time0, nrisk0, nevent0, surv0;
    time0.reserve(n2); nrisk0.reserve(n2); nevent0.reserve(n2); surv0.reserve(n2);
    
    double t = 0, nrisk = n2, nevent = 0, surv = 1.0;
    bool cache = false;
    for (int j = 0; j < n2; ++j) {
      if ((j == 0 && event2[j] == 1) || (j >= 1 && event2[j] == 1 && 
          time2[j] > time2[j-1])) {
        // new event
        // add the info for the previous event
        if (cache) {
          surv *= (1.0 - nevent / nrisk);
          
          time0.push_back(t);
          nrisk0.push_back(nrisk);
          nevent0.push_back(nevent);
          surv0.push_back(surv);
        }
        
        // update the buffer for the current event time
        t = time2[j];
        nrisk = n2-j;
        nevent = 1;
        
        cache = true;
      } else if (j >= 1 && event2[j] == 1 && event2[j-1] == 1 && 
        time2[j] == time2[j-1]) { 
        // tied event
        ++nevent;
      } else if (j >= 1 && event2[j] == 0 && event2[j-1] == 1) {
        // new censoring
        // add the info for the previous event
        surv *= (1.0 - nevent / nrisk);
        
        time0.push_back(t);
        nrisk0.push_back(nrisk);
        nevent0.push_back(nevent);
        surv0.push_back(surv);
        
        // empty the cache for the current event time
        cache = false;
      }
    }
    
    // add the info for the last event
    if (cache) {
      surv *= (1.0 - nevent / nrisk);
      
      time0.push_back(t);
      nrisk0.push_back(nrisk);
      nevent0.push_back(nevent);
      surv0.push_back(surv);
    }
    
    // locate the latest event time before milestone
    std::vector<double> milestone1(1, milestone);
    int N = findInterval3(milestone1, time0)[0];
    
    // prepend time zero information
    time0.insert(time0.begin(), 0.0);
    nrisk0.insert(nrisk0.begin(), n2);
    nevent0.insert(nevent0.begin(), 0.0);
    surv0.insert(surv0.begin(), 1.0);
    
    // replace the last time of interest with milestone
    if (N == static_cast<int>(time0.size()) - 1) {
      time0.push_back(milestone);
    } else {
      time0[N+1] = milestone;
    }
    
    // calculate the partial sum of the trapezoid integration
    std::vector<double> rmstx(N+1);
    rmstx[0] = surv0[0] * (time0[1] - time0[0]);
    for (int k = 1; k <= N; ++k) {
      rmstx[k] = rmstx[k-1] + surv0[k] * (time0[k+1] - time0[k]);
    }
    
    // calculate rmst and its variance
    double u = rmstx[N];
    double v = 0.0;
    for (int k = 1; k <= N; ++k) {
      // rmst from the kth event time to milestone
      double a = u - rmstx[k-1];
      // do not add variance if the largest observed time is an event time
      if (nrisk0[k] > nevent0[k]) {
        v += nevent0[k] * a * a / (nrisk0[k] * (nrisk0[k] - nevent0[k]));
      }
    }
    
    // apply bias correction if requested
    if (biascorrection) {
      double m1 = 0;
      for (int k = 1; k <= N; ++k) {
        m1 += nevent0[k];
      }
      
      if (m1 <= 1.0) {
        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else {
          for (int j = 0; j < p_stratum; ++j) {
            std::string s = stratum[j];
            if (u_stratum.int_cols.count(s)) {
              auto v = u_stratum.get<int>(s);
              stratumerr += " " + s + " = " + std::to_string(v[i]);
            } else if (u_stratum.numeric_cols.count(s)) {
              auto v = u_stratum.get<double>(s);
              stratumerr += " " + s + " = " + std::to_string(v[i]);
            } else if (u_stratum.string_cols.count(s)) {
              auto v = u_stratum.get<std::string>(s);
              stratumerr += " " + s + " = " + v[i];
            } else {
              throw std::invalid_argument(
                  "unsupported type for stratum variable " + s);
            }
          }
        }
        
        std::string str1 = "Bias correction is not done due to no or";
        std::string str2 = "only 1 event before the milestone time:";
        std::string errmsg = str1 + " " + str2;
        if (!stratumerr.empty()) {
          errmsg = errmsg + " " + stratumerr;
        }
        
        if (noerr) {
          thread_utils::push_thread_warning(
            errmsg + "\nAdditional warning messages are suppressed.");
          noerr = false;
        }
        
        continue;
      } else {
        v = m1 / (m1 - 1.0) * v;
      }
    }
    
    double se = std::sqrt(v);
    stratum0.push_back(i);
    size0.push_back(n2);
    rmst0.push_back(u);
    stderr0.push_back(se);
    lower0.push_back(u - z * se);
    upper0.push_back(u + z * se);
  }
  
  DataFrameCpp result;
  result.push_back(std::move(size0), "size");
  result.push_back(milestone, "milestone");
  result.push_back(std::move(rmst0), "rmst");
  result.push_back(std::move(stderr0), "stderr");
  result.push_back(std::move(lower0), "lower");
  result.push_back(std::move(upper0), "upper");
  result.push_back(conflev, "conflev");
  result.push_back(biascorrection, "biascorrection");
  
  if (has_stratum) {
    for (int i = 0; i < p_stratum; ++i) {
      std::string s = stratum[i];
      if (u_stratum.int_cols.count(s)) {
        auto v = u_stratum.get<int>(s);
        result.push_back(subset(v, stratum0), s);
      } else if (u_stratum.numeric_cols.count(s)) {
        auto v = u_stratum.get<double>(s);
        result.push_back(subset(v, stratum0), s);
      } else if (u_stratum.string_cols.count(s)) {
        auto v = u_stratum.get<std::string>(s);
        result.push_back(subset(v, stratum0), s);
      } else {
        throw std::invalid_argument("unsupported type for stratum variable " + s);
      }
    }
  }
  
  return result;
}


//' @title Estimate of Restricted Mean Survival Time
//' @description Obtains the estimate of restricted means survival time
//' for each stratum.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{time}: The possibly right-censored survival time.
//'
//'   * \code{event}: The event indicator.
//'
//' @param stratum The name of the stratum variable in the input data.
//' @param time The name of the time variable in the input data.
//' @param event The name of the event variable in the input data.
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param conflev The level of the two-sided confidence interval for
//'   the survival probabilities. Defaults to 0.95.
//' @param biascorrection Whether to apply bias correction for the
//'   variance estimate. Defaults to no bias correction.
//'
//' @return A data frame with the following variables:
//'
//' * \code{stratum}: The stratum variable.
//'
//' * \code{size}: The number of subjects in the stratum.
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{rmst}: The estimate of restricted mean survival time.
//'
//' * \code{stderr}: The standard error of the estimated rmst.
//'
//' * \code{lower}: The lower bound of confidence interval if requested.
//'
//' * \code{upper}: The upper bound of confidence interval if requested.
//'
//' * \code{conflev}: The level of confidence interval if requested.
//'
//' * \code{biascorrection}: Whether to apply bias correction for the
//'   variance estimate.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' rmest(data = aml, stratum = "x",
//'       time = "time", event = "status", milestone = 24)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame rmest(const Rcpp::DataFrame& data,
                      const Rcpp::StringVector& stratum = "",
                      const std::string& time = "time",
                      const std::string& event = "event",
                      const double milestone = 0,
                      const double conflev = 0.95,
                      const bool biascorrection = false) {
  
  auto dfcpp = convertRDataFrameToCpp(data);
  auto stratumcpp = Rcpp::as<std::vector<std::string>>(stratum);
  
  auto cpp_result = rmestcpp(
    dfcpp, stratumcpp, time, event, milestone, conflev, biascorrection
  );
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(cpp_result);
}



DataFrameCpp rmdiffcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat,
                       const std::string& time,
                       const std::string& event,
                       const double milestone,
                       const double rmstDiffH0,
                       const double conflev,
                       const bool biascorrection) {
  int n = data.nrows();
  
  bool has_stratum = false;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
  }
  
  if (!data.containElementNamed(treat)) {
    throw std::invalid_argument("data must contain the treat variable");
  }
  
  // create the numeric treat variable
  std::vector<int> treatn(n);
  std::vector<int> treatwi;
  std::vector<double> treatwn;
  std::vector<std::string> treatwc;
  if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
    std::vector<int> treatv(n);
    if (data.bool_cols.count(treat)) {
      const std::vector<unsigned char>& treatvb = data.get<unsigned char>(treat);
      for (int i = 0; i < n; ++i) treatv[i] = treatvb[i] ? 1 : 0;
    } else treatv = data.get<int>(treat);
    treatwi = unique_sorted(treatv); // obtain unique treatment values
    if (treatwi.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwi.begin(), treatwi.end(), [](int v) { 
      return v == 0 || v == 1; })) {
      treatwi = {1, 0}; // special handling for 1 / 0 treatment coding
      for (int i = 0; i < n; ++i) treatn[i] = 2 - treatv[i];
    } else {
      treatn = matchcpp(treatv, treatwi, 1);
    }
  } else if (data.numeric_cols.count(treat)) {
    const std::vector<double>& treatv = data.get<double>(treat);
    treatwn = unique_sorted(treatv);
    if (treatwn.size() != 2)
      throw std::invalid_argument("treat must have two and only two distinct values");
    if (std::all_of(treatwn.begin(), treatwn.end(), [](double v) { 
      return v == 0.0 || v == 1.0; })) {
      treatwn = {1.0, 0.0};
      for (int i = 0; i < n; ++i) treatn[i] = 2 - static_cast<int>(treatv[i]);
    } else {
      treatn = matchcpp(treatv, treatwn, 1);
    }
  } else if (data.string_cols.count(treat)) {
    const std::vector<std::string>& treatv = data.get<std::string>(treat);
    treatwc = unique_sorted(treatv);
    if (treatwc.size() != 2) 
      throw std::invalid_argument("treat must have two and only two distinct values");
    treatn = matchcpp(treatv, treatwc, 1);
  } else {
    throw std::invalid_argument(
        "incorrect type for the treat variable in the input data");
  }
  
  if (!data.containElementNamed(time))
    throw std::invalid_argument("data must contain the time variable");
  
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (double v : timen) if (v < 0.0)
    throw std::invalid_argument("time must be nonnegative for each subject");
  
  if (!data.containElementNamed(event))
    throw std::invalid_argument("data must contain the event variable");
  std::vector<int> eventn(n);
  if (data.bool_cols.count(event)) {
    const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
    for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
  } else if (data.int_cols.count(event)) {
    eventn = data.get<int>(event);
  } else if (data.numeric_cols.count(event)) {
    const std::vector<double>& vd = data.get<double>(event);
    for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
  } else {
    throw std::invalid_argument("event variable must be bool, integer or numeric");
  }
  for (double val : eventn) if (val != 0 && val != 1)
    throw std::invalid_argument("event must be 1 or 0 for each observation");
  
  
  if (std::isnan(milestone)) {
    throw std::invalid_argument("milestone must be provided");
  }
  
  if (milestone <= 0) {
    throw std::invalid_argument("milestone must be positive");
  }
  
  if (conflev <= 0.0 || conflev >= 1.0) {
    throw std::invalid_argument("conflev must lie between 0 and 1");
  }
  
  bool noerr = true;
  
  DataFrameCpp dfin; 
  dfin.push_back(std::move(stratumn), "stratum");
  dfin.push_back(std::move(treatn), "treat");
  dfin.push_back(std::move(timen), "time");
  dfin.push_back(std::move(eventn), "event");
  
  DataFrameCpp dfout = rmestcpp(dfin, {"stratum", "treat"}, "time", "event",
                                milestone, 0.95, biascorrection);
  
  std::vector<int> stratum1 = dfout.get<int>("stratum");
  std::vector<int> treat1 = dfout.get<int>("treat");
  std::vector<int> treatsize = dfout.get<int>("size");
  std::vector<double> rmstime1 = dfout.get<double>("rmst");
  std::vector<double> stderr1 = dfout.get<double>("stderr");
  int n1 = stratum1.size();
  
  // identify the locations of the unique values of stratum
  std::vector<int> idx(1,0);
  if (has_stratum) {
    for (int i = 1; i < n1; ++i) {
      if (stratum1[i] != stratum1[i-1]) {
        idx.push_back(i);
      }
    }
  }
  
  int nstrata = idx.size();
  idx.push_back(n1);
  
  std::vector<int> m(nstrata, 0); // number of subjects in each stratum
  for (int i = 0; i < nstrata; ++i) {
    int j1 = idx[i], j2 = idx[i+1] - 1;
    if (treat1[j1] != 1 || treat1[j2] != 2) {
      std::string stratumerr;
      if (!has_stratum) {
        stratumerr = "";
      } else {
        for (int j = 0; j < p_stratum; ++j) {
          std::string s = stratum[j];
          if (u_stratum.int_cols.count(s)) {
            auto v = u_stratum.get<int>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.numeric_cols.count(s)) {
            auto v = u_stratum.get<double>(s);
            stratumerr += " " + s + " = " + std::to_string(v[i]);
          } else if (u_stratum.string_cols.count(s)) {
            auto v = u_stratum.get<std::string>(s);
            stratumerr += " " + s + " = " + v[i];
          } else {
            throw std::invalid_argument("unsupported type for stratum variable " + s);
          }
        }
      }
      
      int k = treat1[j1] != 1 ? 0 : 1;
      std::string treaterr;
      if (data.bool_cols.count(treat) || data.int_cols.count(treat)) {
        treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
      } else if (data.numeric_cols.count(treat)) {
        treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
      } else if (data.string_cols.count(treat)) {
        treaterr = " " + treat + " = " + treatwc[k];
      } else {
        throw std::invalid_argument("unsupported type for treat variable " + treat);
      }
      
      std::string str1 = "The data set does not contain";
      std::string errmsg = str1 + treaterr;
      if (!stratumerr.empty()) errmsg = errmsg + " " + stratumerr;
      
      if (noerr) {
        thread_utils::push_thread_warning(
          errmsg + "\nAdditional warning messages are suppressed.");
        noerr = false;
      }
      
      continue;
    }
    
    m[i] += treatsize[j1] + treatsize[j2];
  }
  
  double M = std::accumulate(m.begin(), m.end(), 0.0);
  std::vector<double> p(nstrata);
  
  double rmst1 = 0.0, rmst2 = 0.0, vrmst1 = 0.0, vrmst2 = 0.0;
  for (int i = 0; i < nstrata; ++i) {
    p[i] = m[i] / M; // fraction of subjects in the stratum
    int start = idx[i], end = idx[i+1];
    std::vector<double> rmst = subset(rmstime1, start, end);
    std::vector<double> stderrx = subset(stderr1, start, end);
    std::vector<double> vrmst(2); 
    for (int j = 0; j < 2; ++j) vrmst[j] = stderrx[j] * stderrx[j];
    rmst1 += p[i] * rmst[0];
    rmst2 += p[i] * rmst[1];
    vrmst1 += p[i] * p[i] * vrmst[0];
    vrmst2 += p[i] * p[i] * vrmst[1];
  }
  
  double z = boost_qnorm((1.0 + conflev) / 2.0);
  
  double rmstDiff = rmst1 - rmst2;
  double sermstDiff= std::sqrt(vrmst1 + vrmst2);
  double rmstDiffZ = (rmstDiff - rmstDiffH0) / sermstDiff;
  double rmstDiffPValue = 2.0 * boost_pnorm(-std::fabs(rmstDiffZ));
  double lower = rmstDiff - z * sermstDiff;
  double upper = rmstDiff + z * sermstDiff;
  
  DataFrameCpp result;
  result.push_back(milestone, "milestone");
  result.push_back(rmstDiffH0, "rmstDiffH0");
  result.push_back(rmst1, "rmst1");
  result.push_back(rmst2, "rmst2");
  result.push_back(rmstDiff, "rmstDiff");
  result.push_back(vrmst1, "vrmst1");
  result.push_back(vrmst2, "vrmst2");
  result.push_back(sermstDiff, "sermstDiff");
  result.push_back(rmstDiffZ, "rmstDiffZ");
  result.push_back(rmstDiffPValue, "rmstDiffPValue");
  result.push_back(lower, "lower");
  result.push_back(upper, "upper");
  result.push_back(conflev, "conflev");
  result.push_back(biascorrection, "biascorrection");
  
  return result;
}


//' @title Estimate of Restricted Mean Survival Time Difference
//' @description Obtains the estimate of restricted mean survival time
//' difference between two treatment groups.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{treat}: The treatment.
//'
//'   * \code{time}: The possibly right-censored survival time.
//'
//'   * \code{event}: The event indicator.
//'
//' @param stratum The name of the stratum variable in the input data.
//' @param treat The name of the treatment variable in the input data.
//' @param time The name of the time variable in the input data.
//' @param event The name of the event variable in the input data.
//' @param milestone The milestone time at which to calculate the
//'   restricted mean survival time.
//' @param rmstDiffH0 The difference in restricted mean survival times
//'   under the null hypothesis. Defaults to 0 for superiority test.
//' @param conflev The level of the two-sided confidence interval for
//'   the difference in restricted mean survival times. Defaults to 0.95.
//' @param biascorrection Whether to apply bias correction for the
//'   variance estimate of individual restricted mean survival times.
//'   Defaults to no bias correction.
//'
//' @return A data frame with the following variables:
//'
//' * \code{milestone}: The milestone time relative to randomization.
//'
//' * \code{rmstDiffH0}: The difference in restricted mean survival times
//'   under the null hypothesis.
//'
//' * \code{rmst1}: The estimated restricted mean survival time for
//'   the treatment group.
//'
//' * \code{rmst2}: The estimated restricted mean survival time for
//'   the control group.
//'
//' * \code{rmstDiff}: The estimated difference in restricted mean
//'   survival times.
//'
//' * \code{vrmst1}: The variance for rmst1.
//'
//' * \code{vrmst2}: The variance for rmst2.
//'
//' * \code{sermstDiff}: The standard error for rmstDiff.
//'
//' * \code{rmstDiffZ}: The Z-statistic value.
//'
//' * \code{rmstDiffPValue}: The two-sided p-value.
//'
//' * \code{lower}: The lower bound of confidence interval.
//'
//' * \code{upper}: The upper bound of confidence interval.
//'
//' * \code{conflev}: The level of confidence interval.
//'
//' * \code{biascorrection}: Whether to apply bias correction for the
//'   variance estimate of individual restricted mean survival times.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' rmdiff(data = rawdata[rawdata$iterationNumber == 1, ], 
//'        stratum = "stratum", treat = "treatmentGroup",
//'        time = "timeUnderObservation", event = "event",
//'        milestone = 12)
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame rmdiff(const Rcpp::DataFrame& data,
                       const Rcpp::StringVector& stratum = "",
                       const std::string& treat = "treat",
                       const std::string& time = "time",
                       const std::string& event = "event",
                       const double milestone = 0,
                       const double rmstDiffH0 = 0,
                       const double conflev = 0.95,
                       const bool biascorrection = false) {
  
  auto dfcpp = convertRDataFrameToCpp(data);
  auto stratumcpp = Rcpp::as<std::vector<std::string>>(stratum);
  
  auto cpp_result = rmdiffcpp(
    dfcpp, stratumcpp, treat, time, event, milestone,
    rmstDiffH0, conflev,  biascorrection
  );
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(cpp_result);
}


struct aftparams {
  int dist_code; // 1: exponential, 2: weibull, 3: lognormal, 4: normal, 
  // 5: loglogistic, 6: logistic
  std::vector<int> strata;
  std::vector<double> tstart;
  std::vector<double> tstop;
  std::vector<int> status;
  std::vector<double> weight;
  std::vector<double> offset;
  FlatMatrix z;
  int nstrata;
};


// function for log-likelihood, score, and information matrix for the AFT model
ListCpp f_der_1(int p, const std::vector<double>& par, void* ex) {
  aftparams *param = (aftparams *) ex;
  const int n = param->z.nrow;
  const int nvar = param->z.ncol;
  const int dist_code = param->dist_code;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& status = param->status; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const double* zptr = param->z.data_ptr(); // column-major: zptr[col * n + row]
  
  // compute linear predictor eta efficiently using column-major storage:
  std::vector<double> eta = offset; // initialize with offset
  // add contributions of each coefficient times column
  for (int i = 0; i < nvar; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    const double* zcol = zptr + i * n;
    for (int r = 0; r < n; ++r) {
      eta[r] += beta * zcol[r];
    }
  }
  
  // Precompute sigma and logsigma per person (exponential has sigma = 1)
  std::vector<double> sigma(n, 1.0); 
  std::vector<double> logsigma(n, 0.0); 
  if (dist_code != 1) { // all except exponential 
    for (int person = 0; person < n; ++person) { 
      int k = strata[person] + nvar; // index in par for log(sigma) 
      sigma[person] = std::exp(par[k]); 
      logsigma[person] = par[k]; // par[k] == log(sigma) 
    } 
  }
  
  // Initialize accumulators
  double loglik = 0.0;
  std::vector<double> score(p);
  FlatMatrix imat(p,p);
  
  // Main loop over persons 
  for (int person = 0; person < n; ++person) {
    const double wt = weight[person]; 
    const double s = sigma[person]; 
    const double inv_s = 1.0 / s; 
    const double logsig = logsigma[person]; 
    const double eta_p = eta[person]; 
    const double tstart_p = tstart[person]; 
    const double tstop_p = tstop[person]; 
    const int st = status[person]; 
    const int k = strata[person] + nvar;
    
    // helper to get z_{j}(person) / sigma without allocating
    auto z = [&](int j)->double {
      return zptr[j * n + person] * inv_s;
    };
    
    switch (st) {
    case 1: // event
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
      double u = (std::log(tstop_p) - eta_p) * inv_s;
      double eu = std::exp(u);
      loglik += wt * (u - eu - logsig);
      
      double c1 = -wt * (1.0 - eu);
      for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
      if (dist_code == 2)
        score[k] += wt * ((1.0 - eu) * (-u) - 1.0);
      
      c1 = wt * eu;
      for (int j = 0; j < nvar; ++j) {
        double zj = z(j);
        for (int i = 0; i <= j; ++i) {
          imat(i, j) += c1 * z(i) * zj;
        }
      }
      if (dist_code == 2) { // weibull
        double c2 = wt * (eu * u - (1.0 - eu));
        for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
        imat(k, k) += c2 * u;
      }
      break;
    }
      case 3: case 4: { // lognormal / normal
        double u;
        if (dist_code == 3) u = (std::log(tstop_p) - eta_p) * inv_s;
        else u = (tstop_p - eta_p) * inv_s;
        loglik += wt * (std::log(boost_dnorm(u)) - logsig);
        
        double c1 = wt * u;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += wt * (u * u - 1.0);
        
        // information: beta-beta
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += wt * z(i) * zj;
          }
        }
        double c2 = wt * 2.0 * u;
        for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
        imat(k, k) += c2 * u;
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u;
        if (dist_code == 5) u = (std::log(tstop_p) - eta_p) * inv_s;
        else u = (tstop_p - eta_p) * inv_s;
        loglik += wt * (std::log(boost_dlogis(u)) - logsig);
        
        double c = 1.0 - 2.0 * boost_plogis(u, 0.0, 1.0, 0);
        double c1 = wt * c;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += wt * (c * u - 1.0);
        
        c1 = wt * 2.0 * boost_dlogis(u);
        double c2 = wt * (2.0 * boost_dlogis(u) * u + 
                          1.0 - 2.0 * boost_plogis(u, 0.0, 1.0, 0));
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
        imat(k, k) += c2 * u;
        break;
      }
      }
      break;
      
    case 3: // interval censoring
      switch (dist_code) {
      case 1: case 2: {
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double u2 = (std::log(tstop_p)  - eta_p) * inv_s;
        double e_u1 = std::exp(u1);
        double e_u2 = std::exp(u2);
        double q1 = std::exp(-e_u1);
        double q2 = std::exp(-e_u2);
        double d1 = e_u1 * q1;
        double d2 = e_u2 * q2;
        double num = d1 - d2;
        double den = q1 - q2;
        double tmp = num / den;
        double ddu = d1 * u1 - d2 * u2;
        double term = ddu / den;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        if (dist_code == 2) {
          score[k] += wt * term;
        }
        
        c1 = wt * (tmp * tmp + (d1 * (1.0 - e_u1) - d2 * (1.0 - e_u2)) / den);
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        if (dist_code == 2) {
          double du1 = d1 * (1.0 + (1.0 - e_u1) * u1);
          double du2 = d2 * (1.0 + (1.0 - e_u2) * u2);
          double c2 = wt * (tmp * term +  (du1 - du2) / den);
          for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
          imat(k, k) += wt * (term * term + (du1 * u1 - du2 * u2) / den);
        }
        break;
      }
      case 3: case 4: {
        double u1, u2;
        if (dist_code == 3) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p)  - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p  - eta_p) * inv_s;
        }
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0);
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0);
        double d1 = boost_dnorm(u1);
        double d2 = boost_dnorm(u2);
        double num = d1 - d2;
        double den = q1 - q2;
        double tmp = num / den;
        double ddu = d1 * u1 - d2 * u2;
        double term = ddu / den;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += wt * term;
        
        c1 = wt * ( tmp * tmp - term );
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double du1 = d1 * (1.0 - u1 * u1);
        double du2 = d2 * (1.0 - u2 * u2);
        double c2 = wt * ( tmp * term + (du1 - du2) / den );
        for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
        imat(k, k) += wt * ( term * term + (du1 * u1 - du2 * u2) / den );
        break;
      }
      case 5: case 6: {
        double u1, u2;
        if (dist_code == 5) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p)  - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p  - eta_p) * inv_s;
        }
        double q1 = boost_plogis(u1, 0.0, 1.0, 0);
        double q2 = boost_plogis(u2, 0.0, 1.0, 0);
        double d1 = boost_dlogis(u1);
        double d2 = boost_dlogis(u2);
        double num = d1 - d2;
        double den = q1 - q2;
        double tmp = num / den;
        double ddu = d1 * u1 - d2 * u2;
        double term = ddu / den;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += wt * term;
        
        c1 = wt * (tmp * tmp + (d1 * (2.0 * q1 - 1.0) - d2 * (2.0 * q2 - 1.0)) / den);
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double du1 = d1 * (1.0 + (2.0 * q1 - 1.0) * u1);
        double du2 = d2 * (1.0 + (2.0 * q2 - 1.0) * u2);
        double c2 = wt * (tmp * term + (du1 - du2) / den);
        for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
        imat(k, k) += wt * ( term * term + (du1 * u1 - du2 * u2) / den );
        break;
      }
      } // dist_code
      break;
      
    case 2: // left censoring
      switch (dist_code) {
      case 1: case 2: {
        double u2 = (std::log(tstop_p) - eta_p) * inv_s;
        double e_u2 = std::exp(u2);
        double q2 = std::exp(-e_u2);
        double d2 = e_u2 * q2;
        double num = -d2;
        double den = 1.0 - q2;
        double tmp = num / den;
        double term = tmp * u2;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        if (dist_code == 2) score[k] += c1 * u2;
        
        c1 = wt * ( tmp * tmp - d2 * (1.0 - e_u2) / den );
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        if (dist_code == 2) {
          double du2 = d2 * (1.0 + (1.0 - e_u2) * u2);
          double c2 = wt * ( tmp * term - du2 / den );
          for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
          imat(k, k) += c2 * u2;
        }
        break;
      }
      case 3: case 4: {
        double u2;
        if (dist_code == 3) u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0); 
        double d2 = boost_dnorm(u2);
        double num = -d2;
        double den = 1.0 - q2;
        double tmp = num / den;
        double term = tmp * u2;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += c1 * u2;
        
        c1 = wt * ( tmp * tmp - term );
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double du2 = d2 * (1.0 - u2 * u2);
        double c2 = wt * ( tmp * term - du2 / den );
        for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
        imat(k, k) += c2 * u2;
        break;
      }
      case 5: case 6: {
        double u2;
        if (dist_code == 5) u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_plogis(u2, 0.0, 1.0, 0); 
        double d2 = boost_dlogis(u2);
        double num = -d2;
        double den = 1.0 - q2;
        double tmp = num / den;
        
        loglik += wt * std::log(den);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += c1 * u2;
        
        c1 = wt * d2;
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double c2 = wt * (d2 * u2 - q2);
        for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
        imat(k, k) += c2 * u2;
        break;
      }
      }
      break;
      
    case 0: // right censoring
      switch (dist_code) {
      case 1: case 2: {
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double e_u1 = std::exp(u1);
        loglik += wt * (-e_u1);
        
        double c1 = wt * e_u1;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        if (dist_code == 2) score[k] += c1 * u1;
        
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        if (dist_code == 2) {
          double c2 = wt * e_u1 * (1.0 + u1);
          for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
          imat(k, k) += c2 * u1;
        }
        break;
      }
      case 3: case 4: {
        double u1;
        if (dist_code == 3) u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else u1 = (tstart_p - eta_p) * inv_s;
        double d1 = boost_dnorm(u1);
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0);
        double tmp = d1 / q1;
        double term = tmp * u1;
        
        loglik += wt * std::log(q1);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += c1 * u1;
        
        c1 = wt * ( tmp * tmp - term );
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double c2 = wt * ( tmp * term + tmp * (1.0 - u1 * u1) );
        for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
        imat(k, k) += c2 * u1;
        break;
      }
      case 5: case 6: {
        double u1;
        if (dist_code == 5) u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else u1 = (tstart_p - eta_p) * inv_s;
        double q1 = boost_plogis(u1, 0.0, 1.0, 0);
        double d1 = boost_dlogis(u1); 
        double tmp  = d1 / q1;
        
        loglik += wt * std::log(q1);
        
        double c1 = wt * tmp;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z(i);
        score[k] += c1 * u1;
        
        c1 = wt * d1;
        for (int j = 0; j < nvar; ++j) {
          double zj = z(j);
          for (int i = 0; i <= j; ++i) {
            imat(i, j) += c1 * z(i) * zj;
          }
        }
        
        double c2 = wt * (1.0 - q1 + d1 * u1);
        for (int j = 0; j < nvar; ++j) imat(j, k) += c2 * z(j);
        imat(k, k) += c2 * u1;
        break;
      }
      }
      break;
      
    default:
      throw std::runtime_error("Unknown status: " + std::to_string(st));
    } // switch(status)
  } // person loop
  
  // mirror lower triangle to upper triangle (imat is symmetric) 
  for (int j = 0; j < p - 1; ++j) 
    for (int i = j + 1; i < p; ++i) 
      imat(i, j) = imat(j, i);
  
  // Build result 
  ListCpp result; 
  result.push_back(loglik, "loglik"); 
  result.push_back(std::move(score), "score"); 
  result.push_back(std::move(imat), "imat"); 
  return result; 
}


// score residual matrix
FlatMatrix f_ressco_1(int p, const std::vector<double>& par, void *ex) {
  aftparams *param = (aftparams *) ex;
  const int n = param->z.nrow;
  const int nvar = param->z.ncol;
  const int dist_code = param->dist_code;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& status = param->status; 
  const std::vector<double>& offset = param->offset; 
  double* zptr = param->z.data_ptr(); // column-major: zptr[col * n + row]
  
  // compute linear predictor eta efficiently using column-major storage:
  std::vector<double> eta = offset; // initialize with offset
  // add contributions of each coefficient times column
  for (int i = 0; i < nvar; ++i) {
    double beta = par[i];
    if (beta == 0.0) continue;
    const double* zcol = zptr + i * n;
    for (int r = 0; r < n; ++r) {
      eta[r] += beta * zcol[r];
    }
  }
  
  // Precompute sigma per person (exponential has sigma = 1)
  std::vector<double> sigma(n, 1.0); 
  if (dist_code != 1) { // all except exponential 
    for (int person = 0; person < n; ++person) { 
      int k = strata[person] + nvar; // index in par for log(sigma) 
      sigma[person] = std::exp(par[k]); 
    } 
  }
  
  // Main loop to compute residuals
  FlatMatrix resid(n, p);
  double* rptr = resid.data_ptr(); // column-major: rptr[col * n + row]
  for (int person = 0; person < n; ++person) {
    const double s = sigma[person];
    const double inv_s = 1.0 / s;
    const double eta_p = eta[person];
    const double tstart_p = tstart[person];
    const double tstop_p = tstop[person];
    const int st = status[person];
    const int k = strata[person] + nvar;
    
    // helper to get z_{j}(person) / sigma without allocating
    auto z = [&](int j)->double {
      return zptr[j * n + person] * inv_s;
    };
    
    // helper to get resid_{j}(person) without allocating
    auto set_resid = [&](int person, int j, double val) {
      rptr[j * n + person] = val;
    };
    
    switch (st) {
    
    case 1: // event
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
      double u = (std::log(tstop_p) - eta_p) * inv_s;
      double c1 = -(1.0 - std::exp(u)); // -f' / f
      for (int j = 0; j < nvar; ++j) {
        set_resid(person, j, c1 * z(j));
      }
      if (dist_code == 2) {
        set_resid(person, k, c1 * u - 1.0);
      }
      break;
    }
      case 3: case 4: { // lognormal / normal
        double u;
        if (dist_code == 3)
          u = (std::log(tstop_p) - eta_p) * inv_s;
        else
          u = (tstop_p - eta_p) * inv_s;
        double c1 = u;
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u - 1.0);
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u;
        if (dist_code == 5)
          u = (std::log(tstop_p) - eta_p) * inv_s;
        else
          u = (tstop_p - eta_p) * inv_s;
        double c1 = 1.0 - 2.0 * boost_plogis(u, 0.0, 1.0, 0);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u - 1.0);
        break;
      }
      }
      break;
      
    case 3: // interval censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double u2 = (std::log(tstop_p) - eta_p) * inv_s;
        double e_u1 = std::exp(u1);
        double e_u2 = std::exp(u2);
        double q1 = std::exp(-e_u1);
        double q2 = std::exp(-e_u2);
        double d1 = e_u1 * q1;
        double d2 = e_u2 * q2;
        double c1 = (d1 - d2) / (q1 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        if (dist_code == 2) {
          set_resid(person, k, (d1 * u1 - d2 * u2) / (q1 - q2));
        }
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u1, u2;
        if (dist_code == 3) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p - eta_p) * inv_s;
        }
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0); 
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0);
        double d1 = boost_dnorm(u1); 
        double d2 = boost_dnorm(u2);
        double c1 = (d1 - d2) / (q1 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, (d1 * u1 - d2 * u2) / (q1 - q2));
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u1, u2;
        if (dist_code == 5) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p - eta_p) * inv_s;
        }
        double q1 = boost_plogis(u1, 0.0, 1.0, 0); 
        double q2 = boost_plogis(u2, 0.0, 1.0, 0);
        double d1 = boost_dlogis(u1); 
        double d2 = boost_dlogis(u2);
        double c1 = (d1 - d2) / (q1 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, (d1 * u1 - d2 * u2) / (q1 - q2));
        break;
      }
      }
      break;
      
    case 2: // left censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u2 = (std::log(tstop_p) - eta_p) * inv_s;
        double e_u2 = std::exp(u2);
        double q2 = std::exp(-e_u2);
        double d2 = e_u2 * q2;
        double c1 = -d2 / (1.0 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        if (dist_code == 2) set_resid(person, k, c1 * u2);
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u2;
        if (dist_code == 3)
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else
          u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0);
        double d2 = boost_dnorm(u2);
        double c1 = -d2 / (1.0 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u2);
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u2;
        if (dist_code == 5)
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else
          u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_plogis(u2, 0.0, 1.0, 0);
        double d2 = boost_dlogis(u2);
        double c1 = -d2 / (1.0 - q2);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u2);
        break;
      }
      }
      break;
      
    case 0: // right censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double c1 = std::exp(u1);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        if (dist_code == 2) set_resid(person, k, c1 * u1);
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u1;
        if (dist_code == 3)
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else
          u1 = (tstart_p - eta_p) * inv_s;
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0);
        double d1 = boost_dnorm(u1);
        double c1 = d1 / q1;
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u1);
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u1;
        if (dist_code == 5)
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else
          u1 = (tstart_p - eta_p) * inv_s;
        double c1 = boost_plogis(u1);
        for (int j = 0; j < nvar; ++j) {
          set_resid(person, j, c1 * z(j));
        }
        set_resid(person, k, c1 * u1);
        break;
      }
      }
      break;
      
    default:
      throw std::runtime_error("Unknown status: " + std::to_string(st));
    }
  }
  
  return resid;
}


// substitute information matrix guaranteed to be positive definite
FlatMatrix f_jj_1(int p, const std::vector<double>& par, void *ex) {
  aftparams *param = (aftparams *) ex;
  const int n = param->z.nrow;
  
  FlatMatrix resid = f_ressco_1(p, par, param);
  FlatMatrix jj(p,p);
  
  // Fast access pointers
  const double* rptr = resid.data_ptr();      // length n * p, column-major
  const double* wptr  = param->weight.data(); // length n
  double* jjptr       = jj.data_ptr();        // length p * p, column-major
  
  // Compute jj(a,b) = sum_{person=0..n-1} w[person] *resid(person,a) *resid(person,b)
  // We compute only lower triangle and mirror to upper triangle for efficiency.
  for (int a = 0; a < p; ++a) {
    const double* colA = rptr + a * n; // resid[:, a]
    for (int b = 0; b <= a; ++b) {
      const double* colB = rptr + b * n; // resid[:, b]
      double sum = 0.0;
      // accumulate dot product of colA and colB, weighted by wptr
      for (int person = 0; person < n; ++person) {
        sum += wptr[person] * colA[person] * colB[person];
      }
      // store at (row=a, col=b) and (row=b, col=a)
      // column-major index: data[col * nrow + row], here nrow == p for jj
      jjptr[ b * p + a ] = sum; // (a, b)
      if (a != b) jjptr[ a * p + b ] = sum; // (b, a) mirror
    }
  }
  
  return jj;
}


// underlying optimization algorithm for lifereg
ListCpp liferegloop(int p, const std::vector<double>& par, void *ex,
                    int maxiter, double eps,
                    const std::vector<int>& colfit, int ncolfit) {
  aftparams *param = (aftparams *) ex;
  
  int iter = 0, halving = 0;
  bool fail = false;
  
  int nstrata = param->nstrata;
  int nvar = param->z.ncol;
  int nsub = param->z.nrow;
  
  FlatMatrix z1 = param->z;
  std::vector<double> mu(nvar), sigma(nvar);
  FlatMatrix z2(nsub, nvar);
  
  // --- standardize z once --- (work column-by-column using column-major layout)
  for (int c = 0; c < nvar; ++c) {
    const double* colptr = z1.data_ptr() + c * nsub;
    // compute mean and sd
    double m, s;
    mean_sd(colptr, nsub, m, s);
    
    // check if column is indicator (only 0 or 1) -> keep m=0, s=1
    bool all_zero_or_one = true;
    for (int r = 0; r < nsub; ++r) {
      double v = colptr[r];
      if (!(v == 0.0 || v == 1.0)) { all_zero_or_one = false; break; }
    }
    if (all_zero_or_one) { m = 0.0; s = 1.0; }
    
    mu[c] = m;
    sigma[c] = s;
    
    // fill standardized column into z2
    for (int r = 0; r < nsub; ++r) {
      z2(r, c) = (colptr[r] - m) / s;
    }
  }
  
  // --- initial beta ---
  std::vector<double> beta(p), newbeta(p);
  beta[0] = par[0];
  for (int i = 1; i < nvar; ++i) {
    beta[i] = par[i] * sigma[i];
    beta[0] += par[i] * mu[i];
  }
  if (param->dist_code != 1)
    std::copy(par.begin() + nvar, par.end(), beta.begin() + nvar);
  
  // local aftparams using standardized covariates z2  
  aftparams para = {param->dist_code, param->strata, param->tstart,
                    param->tstop, param->status, param->weight,
                    param->offset, z2, nstrata};
  
  ListCpp der = f_der_1(p, beta, &para);
  double loglik = der.get<double>("loglik");
  double newlk = 0;
  std::vector<double> u = der.get<std::vector<double>>("score");
  FlatMatrix imat = der.get<FlatMatrix>("imat");
  FlatMatrix jj; // will be used if needed
  std::vector<double> uu1(ncolfit);
  double* u1 = uu1.data();
  FlatMatrix imat1(ncolfit, ncolfit);
  FlatMatrix jj1(ncolfit, ncolfit);
  
  // fill u1 with selected components
  for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
  
  // fill imat1 from imat using colfit indices
  for (int j = 0; j < ncolfit; ++j) {
    for (int i = 0; i < ncolfit; ++i) {
      imat1(i, j) = imat(colfit[i], colfit[j]);
    }
  }
  
  // --- first step: solve system using imat1 (cholesky) or fallback to jj1 ---
  if (cholesky2(imat1, ncolfit) < 0) {
    jj = f_jj_1(p, beta, &para); // substitute information matrix
    for (int j = 0; j < ncolfit; ++j)
      for (int i = 0; i < ncolfit; ++i)
        jj1(i, j) = jj(colfit[i], colfit[j]);
    cholesky2(jj1, ncolfit);
    chsolve2(jj1, ncolfit, u1);
  } else {
    chsolve2(imat1, ncolfit, u1);
  }
  
  // construct update vector u (length p) from solved u1
  std::fill(u.begin(), u.end(), 0.0);
  for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
  
  // newbeta = beta + u
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + u[i];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_1(p, newbeta, &para);
    newlk = der.get<double>("loglik");
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && halving == 0 && std::fabs(1.0 - loglik / newlk) < eps) break;
    
    if (fail || newlk < loglik) {
      ++halving;
      for (int i = 0; i < p; ++i) newbeta[i] = 0.5 * (beta[i] + newbeta[i]);
      
      // special handling of sigmas
      if (halving == 1 && param->dist_code != 1) {
        for (int i = 0; i < nstrata; ++i) {
          int idx = nvar + i;
          if (beta[idx] - newbeta[idx] > 1.1) newbeta[idx] = beta[idx] - 1.1;
        }
      }
      continue;
    }
    
    // --- update step: accept newbeta and compute next increment ---
    halving = 0;
    beta = newbeta;         // copy accepted parameters
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<FlatMatrix>("imat");
    
    // extract relevant components for solving
    for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
    
    for (int j = 0; j < ncolfit; ++j)
      for (int i = 0; i < ncolfit; ++i)
        imat1(i, j) = imat(colfit[i], colfit[j]);
    
    if (cholesky2(imat1, ncolfit) < 0) {
      jj = f_jj_1(p, beta, &para);
      for (int j = 0; j < ncolfit; ++j)
        for (int i = 0; i < ncolfit; ++i)
          jj1(i, j) = jj(colfit[i], colfit[j]);
      cholesky2(jj1, ncolfit);
      chsolve2(jj1, ncolfit, u1);
    } else {
      chsolve2(imat1, ncolfit, u1);
    }
    
    std::fill(u.begin(), u.end(), 0.0);
    for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
    for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + u[i];
  }
  
  if (iter == maxiter) fail = true;
  
  // --- rescale back (undo standardization) ---
  for (int i = 1; i < nvar; ++i) {
    newbeta[i] /= sigma[i];
    newbeta[0] -= newbeta[i] * mu[i];
  }
  
  // --- rescale the information matrix accordingly ---
  imat = der.get<FlatMatrix>("imat");
  FlatMatrix jmat = imat; // copy
  
  // adjust the top-left nvar x nvar block
  for (int j = 0; j < nvar; ++j) {
    for (int i = j; i < nvar; ++i) {
      imat(i, j) = jmat(0,0) * mu[i] * mu[j]
      + jmat(j,0) * mu[i] * sigma[j]
      + jmat(i,0) * mu[j] * sigma[i]
      + jmat(i,j) * sigma[i] * sigma[j];
      if (i != j) imat(j, i) = imat(i, j); // symmetric 
    }
  }
  
  // adjust remaining rows / cols that involve shape / scale parameters
  for (int j = 0; j < nvar; ++j) {
    for (int i = nvar; i < p; ++i) {
      imat(i,j) = jmat(i,0) * mu[j] + jmat(i,j) * sigma[j];
      imat(j,i) = imat(i, j); // symmetric
    }
  }
  
  // compute variance matrix for the fitted parameters (subset colfit)
  for (int j = 0; j < ncolfit; ++j)
    for (int i = 0; i < ncolfit; ++i)
      imat1(i, j) = imat(colfit[i], colfit[j]);
  
  FlatMatrix var1 = invsympd(imat1, ncolfit); // inverse of submatrix
  FlatMatrix var(p, p); // zero-initialized
  // place var1 into the appropriate locations in the full variance matrix
  for (int j = 0; j < ncolfit; ++j)
    for (int i = 0; i < ncolfit; ++i)
      var(colfit[i], colfit[j]) = var1(i, j);
  
  // Build and return result as ListCpp
  ListCpp result;
  result.push_back(std::move(newbeta), "coef");
  result.push_back(iter, "iter");
  result.push_back(std::move(var), "var");
  result.push_back(newlk, "loglik");
  result.push_back(fail, "fail");
  return result;
}


// confidence limit of profile likelihood method
double liferegplloop(int p, const std::vector<double>& par, void *ex,
                     int maxiter, double eps,
                     int k, int direction, double l0) {
  aftparams *param = (aftparams *) ex;
  int iter = 0;
  bool fail = false;
  
  // use std::vector for parameter vectors
  std::vector<double> beta = par;
  std::vector<double> newbeta(p, 0.0);
  double loglik = 0.0;
  double newlk = 0.0;
  
  // containers for score, delta and matrices (FlatMatrix)
  std::vector<double> u(p), delta(p);
  FlatMatrix imat(p, p), jj(p, p), v(p, p);
  
  // --- first step ---
  // Evaluate derivatives / information at initial beta
  ListCpp der = f_der_1(p, beta, param);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  imat = der.get<FlatMatrix>("imat");
  
  // compute inverse of imat (with cholesky fallback to jj)
  // test cholesky on a copy
  jj = imat; // copy
  if (cholesky2(jj, p) < 0) {
    // fallback: compute substitute information jj and invert
    jj = f_jj_1(p, beta, param);
    v = invsympd(jj, p); // inv of jj
  } else {
    v = invsympd(imat, p); // inv of imat
  }
  
  // Lagrange multiplier method used in SAS PROC LOGISTIC
  double w = -quadsym(u, v);
  double underroot = -2.0 * (l0 - loglik + 0.5 * w) / v(k, k);
  double lambda = (underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot));
  u[k] += lambda;
  delta = mat_vec_mult(v, u);
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_1(p, newbeta, param);
    newlk = der.get<double>("loglik");
    
    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && std::fabs(newlk - l0) < eps && w < eps) break;
    
    // update step: accept newbeta and recompute
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<FlatMatrix>("imat");
    jj = imat; // copy for cholesky test
    
    if (cholesky2(jj, p) < 0) {
      jj = f_jj_1(p, beta, param);
      v = invsympd(jj, p);
    } else {
      v = invsympd(imat, p);
    }
    
    // Lagrange multiplier step again
    w = -quadsym(u, v);
    underroot = -2.0 * (l0 - loglik + 0.5 * w) / v(k, k);
    lambda = underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot);
    u[k] += lambda;
    delta = mat_vec_mult(v, u);
    for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  }
  
  if (iter == maxiter) fail = true;
  if (fail) thread_utils::push_thread_warning("liferegplloop did not converge.");
  
  return newbeta[k];
}


// main liferegcpp function
ListCpp liferegcpp(const DataFrameCpp& data,
                   const std::vector<std::string>& stratum,
                   const std::string time,
                   const std::string time2,
                   const std::string event,
                   const std::vector<std::string>& covariates,
                   const std::string weight,
                   const std::string offset,
                   const std::string id,
                   const std::string dist,
                   const std::vector<double>& init,
                   const bool robust,
                   const bool plci,
                   const double alpha,
                   const int maxiter,
                   const double eps) {
  
  // --- sizes and distribution normalization ---
  int n = data.nrows();
  int nvar = static_cast<int>(covariates.size()) + 1;
  if (nvar == 2 && covariates[0] == "") nvar = 1;
  
  std::string dist1 = dist;
  std::for_each(dist1.begin(), dist1.end(), [](char &c){ 
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  if (dist1 == "log-logistic" || dist1 == "llogistic") dist1 = "loglogistic";
  else if (dist1 == "log-normal" || dist1 == "lnormal") dist1 = "lognormal";
  else if (dist1 == "gaussian") dist1 = "normal";
  
  int dist_code = 0;
  if (dist1 == "exponential") dist_code = 1;
  else if (dist1 == "weibull") dist_code = 2;
  else if (dist1 == "lognormal") dist_code = 3;
  else if (dist1 == "normal") dist_code = 4;
  else if (dist1 == "loglogistic") dist_code = 5;
  else if (dist1 == "logistic") dist_code = 6;
  else throw std::invalid_argument("invalid distribution: " + dist1);
  
  // --- handle strata (bygroup) ---
  std::vector<int> stratumn(n);
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    stratumn = out.get<std::vector<int>>("index");
  }
  
  std::vector<int> stratumn1 = unique_sorted(stratumn);
  int nstrata = static_cast<int>(stratumn1.size());
  int p = (dist_code == 1) ? nvar : (nvar + nstrata);
  if (dist_code == 1 && nstrata > 1) {
    throw std::invalid_argument(
        "Stratification is not valid with the exponential distribution");
  }
  
  // --- time / time2 existence and checks ---
  if (!data.containElementNamed(time)) 
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  if ((dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)) {
    for (int i = 0; i < n; ++i) if (!std::isnan(timen[i]) && timen[i] <= 0.0)
      throw std::invalid_argument("time must be positive for the " + 
                                  dist1 + " distribution");
  }
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    if (data.int_cols.count(time2)) {
      const std::vector<int>& vi = data.get<int>(time2);
      for (int i = 0; i < n; ++i) time2n[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(time2)) {
      time2n = data.get<double>(time2);
    } else {
      throw std::invalid_argument("time2 variable must be integer or numeric");
    }    
    if ((dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)) {
      for (int i = 0; i < n; ++i) if (!std::isnan(time2n[i]) && time2n[i] <= 0.0)
        throw std::invalid_argument("time2 must be positive for the " + 
                                    dist1 + " distribution");
    }
  }
  
  // --- event variable ---
  bool has_event = !event.empty() && data.containElementNamed(event);
  if (!has_time2 && !has_event) {
    throw std::invalid_argument(
        "data must contain the event variable for right censored data"); 
  }
  std::vector<int> eventn(n);
  if (has_event) {
    if (data.bool_cols.count(event)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
      for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
    } else if (data.int_cols.count(event)) {
      eventn = data.get<int>(event);
    } else if (data.numeric_cols.count(event)) {
      const std::vector<double>& vd = data.get<double>(event);
      for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
    } else {
      throw std::invalid_argument("event variable must be bool, integer or numeric");
    }
    for (double val : eventn) if (val != 0 && val != 1)
      throw std::invalid_argument("event must be 1 or 0 for each observation");
    
    for (double val : eventn) if (val != 0 && val != 1)
      throw std::invalid_argument("event must be 1 or 0 for each observation");
  }
  
  // --- build design matrix zn (n x nvar) column-major FlatMatrix ---
  FlatMatrix zn(n, nvar);
  // intercept
  for (int i = 0; i < n; ++i) zn.data[i] = 1.0;
  // covariates
  for (int j = 0; j < nvar - 1; ++j) {
    const std::string& zj = covariates[j];
    if (!data.containElementNamed(zj)) 
      throw std::invalid_argument("data must contain the variables in covariates");
    double* zn_col = zn.data_ptr() + (j + 1) * n;
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
      for (int i = 0; i < n; ++i) zn_col[i] = vb[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& vi = data.get<int>(zj);
      for (int i = 0; i < n; ++i) zn_col[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(zj)) {
      const std::vector<double>& vd = data.get<double>(zj);
      std::memcpy(zn_col, vd.data(), n * sizeof(double));
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }
  
  // --- weight and offset ---
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const std::vector<int>& vi = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else {
      throw std::invalid_argument("weight variable must be integer or numeric");
    }
    for (double w : weightn) if (std::isnan(w) || w <= 0.0) 
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    if (data.int_cols.count(offset)) {
      const std::vector<int>& vi = data.get<int>(offset);
      for (int i = 0; i < n; ++i) offsetn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(offset)) {
      offsetn = data.get<double>(offset);
    } else {
      throw std::invalid_argument("offset variable must be integer or numeric");
    }
  }
  
  // --- id mapping ---
  bool has_id = !id.empty() && data.containElementNamed(id);
  std::vector<int> idn(n);
  if (!has_id) {
    std::iota(idn.begin(), idn.end(), 0);
  } else {
    if (data.int_cols.count(id)) {
      auto v = data.get<int>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.numeric_cols.count(id)) {
      auto v = data.get<double>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.string_cols.count(id)) {
      auto v = data.get<std::string>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else throw std::invalid_argument(
        "incorrect type for the id variable in the input data");
  }
  
  // --- exclude observations with missing values ---
  std::vector<unsigned char> sub(n, 1);
  for (int i = 0; i < n; ++i) {
    if (stratumn[i] == INT_MIN || idn[i] == INT_MIN ||
        (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
        std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
      sub[i] = 0;
      continue;
    }
    // check covariates columns
    for (int j = 0; j < nvar - 1; ++j) {
      if (std::isnan(zn(i, j+1))) { sub[i] = 0; break; }
    }
  }
  std::vector<int> keep = which(sub);
  if (keep.empty()) throw std::invalid_argument(
      "no observations without missing values");
  
  subset_in_place(stratumn, keep);
  subset_in_place(timen, keep);
  subset_in_place(time2n, keep);
  subset_in_place(eventn, keep);
  subset_in_place(weightn, keep);
  subset_in_place(offsetn, keep);
  subset_in_place(idn, keep);
  subset_in_place_flatmatrix(zn, keep);
  n = keep.size();
  
  // sumstat data set
  int nobs, nevents;
  double loglik0, loglik1;
  int niter;
  bool fail;
  
  // parest data set
  std::vector<std::string> par(p);
  std::vector<double> b(p), seb(p), rseb(p);
  std::vector<double> z(p), expbeta(p);
  FlatMatrix vb(p, p), rvb(p, p);
  std::vector<double> lb(p), ub(p), prob(p);
  std::vector<std::string> clparm(p);
  
  // linear predictor and fitted values for all observations
  std::vector<double> linear_predictors(n);
  
  double zcrit = boost_qnorm(1.0 - alpha / 2.0);
  double xcrit = zcrit * zcrit;
  
  // unify right censored data with interval censored data
  std::vector<double> tstart(n), tstop(n);
  if (!has_time2) {
    for (int i = 0; i < n; ++i) {
      tstart[i] = timen[i];
      tstop[i] = eventn[i] == 1 ? tstart[i] : NaN;
    }
  } else {
    tstart = timen;
    tstop = time2n;
  }
  
  std::vector<int> status(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) && tstart[i] == tstop[i]) {
      status[i] = 1; // event
    } else if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) && 
      tstart[i] < tstop[i]) {
      status[i] = 3; // interval censoring
    } else if (std::isnan(tstart[i]) && !std::isnan(tstop[i])) {
      status[i] = 2; // left censoring
    } else if (!std::isnan(tstart[i]) && std::isnan(tstop[i])) {
      status[i] = 0; // right censoring
    } else {
      status[i] = -1; // exclude the observation
    }
  }
  
  nobs = n;
  nevents = 0;
  for (int i = 0; i < n; ++i) if (status[i] == 1) ++nevents;
  if (nevents == 0) {
    for (int i = 0; i < p; ++i) {
      if (i == 0) par[i] = "(Intercept)";
      else if (i < nvar) par[i] = covariates[i-1];
      else {
        if (nstrata == 1) par[i] = "Log(scale)";
        else par[i] = std::string("Log(scale ") + std::to_string(i - nvar + 1) + ")";
      }
      
      b[i] = NaN;
      seb[i] = 0;
      rseb[i] = 0;
      z[i] = NaN;
      expbeta[i] = NaN;
      lb[i] = NaN;
      ub[i] = NaN;
      prob[i] = NaN;
      clparm[i] = "Wald";
    }
    
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < p; ++i) {
        vb(i,j) = 0;
        rvb(i,j) = 0;
      }
    }
    
    for (int i = 0; i < n; ++i) {
      linear_predictors[i] = offsetn[i];
    }
    
    loglik0 = NaN;
    loglik1 = NaN;
    niter = 0;
    fail = true;
  } else { // nevents > 0
    // exclude invalid status
    std::vector<unsigned char> good(n, 1);
    for (int i = 0; i < n; ++i) if (status[i] == -1) good[i] = 0;
    std::vector<int> q = which(good);
    int n1 = q.size();
    
    if (n1 < n) {
      subset_in_place(stratumn, q);
      subset_in_place(tstart, q);
      subset_in_place(tstop, q);
      subset_in_place(status, q);
      subset_in_place(weightn, q);
      subset_in_place(offsetn, q);
      subset_in_place(idn, q);
      subset_in_place_flatmatrix(zn, q);
    }
    
    // intercept only model
    std::vector<double> time0(n1);
    for (int i = 0; i < n1; ++i) {
      if (status[i] == 0 || status[i] == 1) { // right censoring or event
        time0[i] = tstart[i];
      } else if (status[i] == 2) { // left censoring
        time0[i] = tstop[i];
      } else if (status[i] == 3) { // interval censoring
        time0[i] = 0.5 * (tstart[i] + tstop[i]);
      }
    }
    
    std::vector<double> y0 = time0;
    if (dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5) {
      for (int i = 0; i < n1; ++i) y0[i] = std::log(y0[i]);
    }
    
    double int0, sig0;
    mean_sd(y0.data(), n1, int0, sig0);
    double logsig0 = std::log(sig0);
    
    std::vector<double> bint0(p);
    int ncolfit0 = (dist_code == 1) ? 1 : (nstrata + 1);
    std::vector<int> colfit0(ncolfit0); // indices of parameters to be fitted
    if (dist_code == 1) {
      bint0[0] = int0;
      colfit0[0] = 0;
    } else {
      bint0[0] = int0;
      for (int i = 0; i < nstrata; ++i) bint0[nvar + i] = logsig0;
      colfit0[0] = 0;
      for (int i = 0; i < nstrata; ++i) colfit0[i + 1] = nvar + i;
    }
    
    // parameter estimates and standard errors for the null model
    aftparams param = {dist_code, stratumn, tstart, tstop, status, 
                       weightn, offsetn, zn, nstrata};
    ListCpp outint = liferegloop(p, bint0, &param, maxiter, eps, colfit0, ncolfit0);
    
    std::vector<double> bint = outint.get<std::vector<double>>("coef");
    FlatMatrix vbint = outint.get<FlatMatrix>("var");
    
    ListCpp out;
    
    if (nvar > 1) {
      std::vector<int> colfit = seqcpp(0, p-1);
      if (!init.empty() && static_cast<int>(init.size()) == p && 
          std::none_of(init.begin(), init.end(), [](double val){ 
            return std::isnan(val); })) {
        out = liferegloop(p, init, &param, maxiter, eps, colfit, p);
      } else {
        out = liferegloop(p, bint, &param, maxiter, eps, colfit, p);
      }
      
      fail = out.get<bool>("fail");
      
      if (fail) {
        // obtain initial values for model parameters using OLS
        std::vector<double> y1(n1);
        for (int i = 0; i < n1; ++i) y1[i] = y0[i] - offsetn[i];
        
        FlatMatrix v1(nvar, nvar); // X'WX
        const double *zptr = zn.data_ptr();
        double *vptr = v1.data_ptr();
        for (int j = 0; j < nvar; ++j) {
          const double* zj = zptr + j * n1;          // pointer to Z(:,j)
          for (int k = j; k < nvar; ++k) {
            const double* zk = zptr + k * n1;        // pointer to Z(:,k)
            double sum = 0.0;
            // inner loop reads zj[i] and zk[i] contiguously
            for (int i = 0; i < n1; ++i) {
              sum += weightn[i] * zj[i] * zk[i];
            }
            // write into v1(j,k) and mirror
            vptr[k * nvar + j] = sum;                
            if (k != j) {
              vptr[j * nvar + k] = sum;            
            }
          }
        }
        
        std::vector<double> uu1(nvar);               // X'Wy
        double* u1 = uu1.data();
        for (int j = 0; j < nvar; ++j) {
          const double* zj = zptr + j * n1;          // pointer to Z(:,j)
          double sum = 0.0;
          // inner loop reads zj[i] contiguously
          for (int i = 0; i < n1; ++i) {
            sum += weightn[i] * zj[i] * y1[i];
          }
          u1[j] = sum;
        }
        
        cholesky2(v1, nvar);
        chsolve2(v1, nvar, u1);
        
        std::vector<double> binit(p);
        for (int j = 0; j < nvar; ++j) binit[j] = u1[j];
        
        if (dist_code != 1) {
          double ssum = 0.0;
          double wsum = 0.0;
          for (int i = 0; i < n1; ++i) {
            double pred = 0.0;
            for (int j = 0; j < nvar; ++j) pred += zn(i,j) * u1[j];
            double r = y1[i] - pred;
            ssum += weightn[i] * r * r;
            wsum += weightn[i];
          }
          double s = 0.5 * std::log(ssum / wsum * n1 / std::max(1, n1 - nvar));  
          // log(sigma)
          for (int j = nvar; j < p; ++j) binit[j] = s;
        }
        
        // fit the model using the initial values
        out = liferegloop(p, binit, &param, maxiter, eps, colfit, p);
        fail = out.get<bool>("fail");
      }
      
      if (fail) thread_utils::push_thread_warning(
          "The algorithm in liferegr did not converge");
      
      b = out.get<std::vector<double>>("coef");
      vb = out.get<FlatMatrix>("var");
    } else {
      // intercept-only
      b = bint;
      vb = vbint;
      out = outint;
    }
    
    // compute standard errors
    for (int j = 0; j < p; ++j) {
      seb[j] = std::sqrt(vb(j,j));
    }
    
    // fill parest outputs
    for (int i = 0; i < p; ++i) {
      if (i == 0) par[i] = "(Intercept)";
      else if (i < nvar) par[i] = covariates[i-1];
      else {
        if (nstrata == 1) par[i] = "Log(scale)";
        else par[i] = std::string("Log(scale ") + std::to_string(i - nvar + 1) + ")";
      }
    }
    
    // fill linear predictors
    for (int i = 0; i < nvar; ++i) {
      double beta = b[i];
      if (beta == 0.0) continue;
      const double* zn_col = zn.data_ptr() + i * n1;
      for (int r = 0; r < n1; ++r) {
        linear_predictors[q[r]] += beta * zn_col[r];
      }
    }
    
    niter = out.get<int>("iter");
    fail = out.get<bool>("fail");
    
    // robust variance estimates
    if (robust) {
      FlatMatrix ressco = f_ressco_1(p, b, &param);
      
      int nr; // number of rows in the score residual matrix
      if (!has_id) { // no clustering, just weight the score residuals
        for (int j = 0; j < p; ++j) {
          double* rcol = ressco.data_ptr() + j * n1;
          for (int i = 0; i < n1; ++i) {
            rcol[i] *= weightn[i];
          }
        }
        nr = n1;
      } else { // sum up the score residuals by id
        std::vector<int> order = seqcpp(0, n1-1);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return idn[i] < idn[j];
        });
        
        std::vector<int> id1 = subset(idn, order);
        std::vector<int> idx(1,0);
        for (int i = 1; i < n1; ++i) {
          if (id1[i] != id1[i-1]) {
            idx.push_back(i);
          }
        }
        int nids = idx.size();
        idx.push_back(n1);
        
        FlatMatrix ressco1(nids, p); // score residuals summed by id
        for (int j = 0; j < p; ++j) {
          const double* rcol = ressco.data_ptr() + j * n1;
          double* rcol1 = ressco1.data_ptr() + j * nids;
          for (int i = 0; i < nids; ++i) {
            double sum = 0.0;
            for (int k = idx[i]; k < idx[i+1]; ++k) {
              int row = order[k];
              sum  += weightn[row] * rcol[row];
            }
            rcol1[i] = sum;
          }
        }
        
        ressco = std::move(ressco1);  // update the score residuals
        nr = nids;
      }
      
      FlatMatrix D = mat_mat_mult(ressco, vb); // DFBETA
      
      const double* Dptr = D.data_ptr();
      double* rvbptr = rvb.data_ptr();
      for (int j = 0; j < p; ++j) {
        const double* Dj = Dptr + j * nr; // pointer to D(:,j)
        for (int k = 0; k <= j; ++k) {
          const double* Dk = Dptr + k * nr; // pointer to D(:,k)
          double sum = 0.0;
          for (int i = 0; i < nr; ++i) {
            sum += Dj[i] * Dk[i];
          }
          rvbptr[k * p + j] = sum;
          if (j != k) rvbptr[j * p + k] = sum;
        }
      }
      
      for (int i = 0; i < p; ++i) {
        rseb[i] = std::sqrt(rvb(i,i));
      }
    }
    
    // profile likelihood confidence interval for regression coefficients
    if (plci) {
      double lmax = out.get<double>("loglik");
      double l0 = lmax - 0.5 * xcrit;
      
      for (int k = 0; k < p; ++k) {
        lb[k] = liferegplloop(p, b, &param, maxiter, eps, k, -1, l0);
        ub[k] = liferegplloop(p, b, &param, maxiter, eps, k, 1, l0);
        
        std::vector<int> colfit1(p-1);
        for (int i = 0, j = 0; i < p; ++i) {
          if (i == k) continue;
          colfit1[j++] = i;
        }
        
        std::vector<double> b0(p);
        ListCpp out0 = liferegloop(p, b0, &param, maxiter, eps, colfit1, p-1);
        double lmax0 = out0.get<double>("loglik");
        prob[k] = boost_pchisq(-2.0 * (lmax0 - lmax), 1, 0);
        clparm[k] = "PL";
      }
    } else {
      for (int k = 0; k < p; ++k) {
        if (!robust) {
          lb[k] = b[k] - zcrit * seb[k];
          ub[k] = b[k] + zcrit * seb[k];
          prob[k] = boost_pchisq(sq(b[k] / seb[k]), 1, 0);
        } else {
          lb[k] = b[k] - zcrit * rseb[k];
          ub[k] = b[k] + zcrit * rseb[k];
          prob[k] = boost_pchisq(sq(b[k] / rseb[k]), 1, 0);
        }
        clparm[k] = "Wald";
      }
    }
    
    // log-likelihoods
    loglik0 = outint.get<double>("loglik");
    loglik1 = out.get<double>("loglik");
    
    // compute exp(beta)
    for (int i = 0; i < p; ++i) {
      expbeta[i] = std::exp(b[i]);
    }
    
    // compute z statistics
    if (robust) {
      for (int i = 0; i < p; ++i) {
        if (rseb[i] == 0) {
          z[i] = NaN;
        } else {
          z[i] = b[i] / rseb[i];
        }
      }
    } else {
      for (int i = 0; i < p; ++i) {
        if (seb[i] == 0) {
          z[i] = NaN;
        } else {
          z[i] = b[i] / seb[i];
        }
      }
    }
  }
  
  DataFrameCpp sumstat;
  sumstat.push_back(nobs, "n");
  sumstat.push_back(nevents, "nevents");
  sumstat.push_back(loglik0, "loglik0");
  sumstat.push_back(loglik1, "loglik1");
  sumstat.push_back(niter, "niter");
  sumstat.push_back(dist1, "dist");
  sumstat.push_back(p, "p");
  sumstat.push_back(nvar - 1, "nvar");
  sumstat.push_back(robust, "robust");
  sumstat.push_back(fail, "fail");
  
  std::vector sebeta = robust ? rseb : seb;
  FlatMatrix vbeta = robust ? rvb : vb;
  
  DataFrameCpp parest;
  parest.push_back(std::move(par), "param");
  parest.push_back(std::move(b), "beta");
  parest.push_back(std::move(sebeta), "sebeta");
  parest.push_back(std::move(z), "z");
  parest.push_back(std::move(expbeta), "expbeta");
  parest.push_back(std::move(lb), "lower");
  parest.push_back(std::move(ub), "upper");
  parest.push_back(std::move(prob), "p");
  parest.push_back(std::move(clparm), "method");
  if (robust) parest.push_back(std::move(seb), "sebeta_naive");

  ListCpp result;
  result.push_back(std::move(sumstat), "sumstat");
  result.push_back(std::move(parest), "parest");
  result.push_back(std::move(vbeta), "vbeta");
  if (robust) result.push_back(std::move(vb), "vbeta_naive");
  result.push_back(std::move(linear_predictors), "linear_predictors");
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List liferegRcpp(const Rcpp::DataFrame& data,
                       const std::vector<std::string>& stratum,
                       const std::string time,
                       const std::string time2,
                       const std::string event,
                       const std::vector<std::string>& covariates,
                       const std::string weight,
                       const std::string offset,
                       const std::string id,
                       const std::string dist,
                       const std::vector<double>& init,
                       const bool robust,
                       const bool plci,
                       const double alpha,
                       const int maxiter,
                       const double eps) {
  
  auto dfcpp = convertRDataFrameToCpp(data);
  
  auto cpp_result = liferegcpp(
    dfcpp, stratum, time, time2, event, covariates, weight, offset, id, 
    dist, init, robust, plci, alpha, maxiter, eps
  );
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(cpp_result);
}


// first and second derivatives of log likelihood with respect to eta and tau
ListCpp f_ld_1(std::vector<double>& eta, std::vector<double>& sigma, void *ex) {
  aftparams *param = (aftparams *) ex;
  const int n = param->tstart.size();
  const int dist_code = param->dist_code;
  
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& status = param->status; 
  
  std::vector<double> g(n), dg(n), ddg(n), ds(n), dds(n), dsg(n);
  
  // Main loop to compute  derivatives
  for (int person = 0; person < n; ++person) {
    const double s = sigma[person];
    const double inv_s = 1.0 / s;
    const double logsig = std::log(s);
    const double eta_p = eta[person];
    const double tstart_p = tstart[person]; 
    const double tstop_p = tstop[person];
    const int st = status[person];
    
    double vg = NaN, vdg = NaN, vddg = NaN, vds = NaN, vdds = NaN, vdsg = NaN;
    
    switch (st) {
    
    case 1: // event
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
      double u = (std::log(tstop_p) - eta_p) * inv_s;
      double eu = std::exp(u);
      double su = s * u;
      vg = u - eu - logsig;        // log f(u) - log(s)
      vdg = -(1.0 - eu) * inv_s;   // -f'(u) / f(u) * 1 / s
      vddg = -eu * inv_s * inv_s;  // (f''(u) / f(u) - (f'(u) / f(u))^2) * 1 / s^2
      vds = -1.0 + vdg * su;
      vdds = vddg * su * su - vdg * su;
      vdsg = vddg * su - vdg;
      break;
    }
      case 3: case 4: { // lognormal / normal
        double u;
        if (dist_code == 3) u = (std::log(tstop_p) - eta_p) * inv_s;
        else u = (tstop_p - eta_p) * inv_s;
        double su = s * u;
        double d = boost_dnorm(u);
        vg = std::log(d) - logsig;
        vdg = u * inv_s;
        vddg = -inv_s * inv_s;
        vds = -1.0 + vdg * su;
        vdds = vddg * su * su - vdg * su;
        vdsg = vddg * su - vdg;
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u;
        if (dist_code == 5) u = (std::log(tstop_p) - eta_p) * inv_s;
        else u = (tstop_p - eta_p) * inv_s;
        double su = s * u;
        double q = boost_plogis(u, 0.0, 1.0, 0);
        double d = boost_dlogis(u);
        vg = std::log(d) - logsig;
        vdg = (1.0 - 2.0 * q) * inv_s;
        vddg = -2.0 * d * inv_s * inv_s;
        vds = -1.0 + vdg * su;
        vdds = vddg * su * su - vdg * su;
        vdsg = vddg * su - vdg;
        break;
      }
      }
      break;
      
    case 3: // interval censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double u2 = (std::log(tstop_p) - eta_p) * inv_s;
        double e_u1 = std::exp(u1); 
        double e_u2 = std::exp(u2);
        double q1 = std::exp(-e_u1); 
        double q2 = std::exp(-e_u2);
        double d1 = e_u1 * q1; 
        double d2 = e_u2 * q2;
        double r1 = 1.0 - e_u1;
        double r2 = 1.0 - e_u2;
        double den = q1 - q2;
        vg = std::log(den);
        vdg = (d1 - d2) / den * inv_s;
        vddg = -(d1 * r1 - d2 * r2) / den * inv_s * inv_s - vdg * vdg;
        vds = (u1 * d1 - u2 * d2) / den;
        vdds = (u2 * u2 * d2 * r2 - u1 * u1 * d1 * r1) / den - vds * (1.0 + vds);
        vdsg = (u2 * d2 * r2 - u1 * d1 * r1) / den * inv_s - vdg * (1.0 + vds);
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u1, u2;
        if (dist_code == 3) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p - eta_p) * inv_s;
        }
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0); 
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0);
        double d1 = boost_dnorm(u1); 
        double d2 = boost_dnorm(u2);
        double r1 = -u1;
        double r2 = -u2;
        double den = q1 - q2;
        vg = std::log(den);
        vdg = (d1 - d2) / den * inv_s;
        vddg = -(d1 * r1 - d2 * r2) / den * inv_s * inv_s - vdg * vdg;
        vds = (u1 * d1 - u2 * d2) / den;
        vdds = (u2 * u2 * d2 * r2 - u1 * u1 * d1 * r1) / den - vds * (1.0 + vds);
        vdsg = (u2 * d2 * r2 - u1 * d1 * r1) / den * inv_s - vdg * (1.0 + vds);
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u1, u2;
        if (dist_code == 5) {
          u1 = (std::log(tstart_p) - eta_p) * inv_s;
          u2 = (std::log(tstop_p) - eta_p) * inv_s;
        } else {
          u1 = (tstart_p - eta_p) * inv_s;
          u2 = (tstop_p - eta_p) * inv_s;
        }
        double d1 = boost_dlogis(u1); 
        double d2 = boost_dlogis(u2);
        double q1 = boost_plogis(u1, 0.0, 1.0, 0); 
        double q2 = boost_plogis(u2, 0.0, 1.0, 0);
        double r1 = 2.0 * q1 - 1.0;
        double r2 = 2.0 * q2 - 1.0;
        double den = q1 - q2;
        vg = std::log(den);
        vdg = (d1 - d2) / den * inv_s;
        vddg = -(d1 * r1 - d2 * r2) / den * inv_s * inv_s - vdg * vdg;
        vds = (u1 * d1 - u2 * d2) / den;
        vdds = (u2 * u2 * d2 * r2 - u1 * u1 * d1 * r1) / den - vds * (1.0 + vds);
        vdsg = (u2 * d2 * r2 - u1 * d1 * r1) / den * inv_s - vdg * (1.0 + vds);
        break;
      }
      }
      break;
      
    case 2: // left censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u2 = (std::log(tstop_p) - eta_p) * inv_s;
        double e_u2 = std::exp(u2);
        double q2 = std::exp(-e_u2);
        double d2 = e_u2 * q2;
        double r2 = 1.0 - e_u2;
        double den = 1.0 - q2;
        double su = s * u2;
        vg = std::log(den);
        vdg = -d2 / den * inv_s;
        vddg = -vdg * r2 * inv_s - vdg * vdg;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u2;
        if (dist_code == 3) u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_pnorm(u2, 0.0, 1.0, 0);
        double d2 = boost_dnorm(u2);
        double r2 = -u2;
        double den = 1.0 - q2;
        double su = s * u2;
        vg = std::log(den);
        vdg = -d2 / den * inv_s;
        vddg = -vdg * r2 * inv_s - vdg * vdg;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u2;
        if (param->dist_code==5) u2 = (std::log(tstop_p) - eta_p) * inv_s;
        else u2 = (tstop_p - eta_p) * inv_s;
        double q2 = boost_plogis(u2, 0.0, 1.0, 0);
        double d2 = boost_dlogis(u2);
        double den = 1.0 - q2;
        double su = s * u2;
        vg = std::log(den);
        vdg = -q2 * inv_s;
        vddg = -d2 * inv_s * inv_s;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      }
      break;
      
    case 0: // right censoring
      switch (dist_code) {
      case 1: case 2: { // exponential / weibull
        double u1 = (std::log(tstart_p) - eta_p) * inv_s;
        double e_u1 = std::exp(u1);
        double su = s * u1;
        vg = -e_u1;
        vdg = e_u1 * inv_s;
        vddg = -e_u1 * inv_s * inv_s;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      case 3: case 4: { // lognormal / normal
        double u1;
        if (dist_code==3) u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else u1 = (tstart_p - eta_p) * inv_s;
        double q1 = boost_pnorm(u1, 0.0, 1.0, 0);
        double d1 = boost_dnorm(u1);
        double r1 = -u1;
        double su = s * u1;
        vg = std::log(q1);
        vdg = d1 / q1 * inv_s;
        vddg = -vdg * r1 * inv_s - vdg * vdg;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      case 5: case 6: { // loglogistic / logistic
        double u1;
        if (dist_code == 5) u1 = (std::log(tstart_p) - eta_p) * inv_s;
        else u1 = (tstart_p - eta_p) * inv_s;
        double q1 = boost_plogis(u1, 0.0, 1.0, 0);
        double d1 = boost_dlogis(u1);
        double su = s * u1;
        vg = std::log(q1);
        vdg = (1.0 - q1) * inv_s;
        vddg = -d1 * inv_s * inv_s;
        vds = vdg * su;
        vdds = vddg * su * su - vds;
        vdsg = vddg * su - vdg;
        break;
      }
      }
      break;
      
    default:
      throw std::runtime_error("Unknown status: " + std::to_string(st));
    }
    
    g[person] = vg;
    dg[person] = vdg;
    ddg[person] = vddg;
    ds[person] = vds;
    dds[person] = vdds;
    dsg[person] = vdsg;
  }
  
  ListCpp result;
  result.push_back(std::move(g), "g");
  result.push_back(std::move(dg), "dg");
  result.push_back(std::move(ddg), "ddg");
  result.push_back(std::move(ds), "ds");
  result.push_back(std::move(dds), "dds");
  result.push_back(std::move(dsg), "dsg");
  
  return result;
}


// residuals of the AFT model
FlatMatrix residuals_liferegcpp(const std::vector<double>& beta,
                                const FlatMatrix& vbeta,
                                const DataFrameCpp& data,
                                const std::vector<std::string>& stratum,
                                const std::string& time,
                                const std::string& time2,
                                const std::string& event,
                                const std::vector<std::string>& covariates,
                                const std::string& weight,
                                const std::string& offset,
                                const std::string& id,
                                const std::string& dist,
                                const std::string& type,
                                bool collapse,
                                bool weighted) {
  // --- sizes and distribution normalization ---
  int n = data.nrows();
  int nvar = static_cast<int>(covariates.size()) + 1;
  if (nvar == 2 && covariates[0] == "") nvar = 1;
  
  std::string dist1 = dist;
  std::for_each(dist1.begin(), dist1.end(), [](char &c){
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  if (dist1 == "log-logistic" || dist1 == "llogistic") dist1 = "loglogistic";
  else if (dist1 == "log-normal" || dist1 == "lnormal") dist1 = "lognormal";
  else if (dist1 == "gaussian") dist1 = "normal";
  
  int dist_code = 0;
  if (dist1 == "exponential") dist_code = 1;
  else if (dist1 == "weibull") dist_code = 2;
  else if (dist1 == "lognormal") dist_code = 3;
  else if (dist1 == "normal") dist_code = 4;
  else if (dist1 == "loglogistic") dist_code = 5;
  else if (dist1 == "logistic") dist_code = 6;
  else throw std::invalid_argument("invalid distribution: " + dist1);
  
  // --- handle strata (bygroup) ---
  std::vector<int> stratumn(n);
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    stratumn = out.get<std::vector<int>>("index");
  }
  
  std::vector<int> stratumn1 = unique_sorted(stratumn);
  int nstrata = static_cast<int>(stratumn1.size());
  int p = (dist_code == 1) ? nvar : (nvar + nstrata);
  if (dist_code == 1 && nstrata > 1) {
    throw std::invalid_argument(
        "Stratification is not valid with the exponential distribution");
  }
  
  // --- time / time2 existence and checks ---
  if (!data.containElementNamed(time)) 
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  if ((dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)) {
    for (int i = 0; i < n; ++i) if (!std::isnan(timen[i]) && timen[i] <= 0.0)
      throw std::invalid_argument("time must be positive for the " + 
                                  dist1 + " distribution");
  }
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    if (data.int_cols.count(time2)) {
      const std::vector<int>& vi = data.get<int>(time2);
      for (int i = 0; i < n; ++i) time2n[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(time2)) {
      time2n = data.get<double>(time2);
    } else {
      throw std::invalid_argument("time2 variable must be integer or numeric");
    }    
    if ((dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)) {
      for (int i = 0; i < n; ++i) if (!std::isnan(time2n[i]) && time2n[i] <= 0.0)
        throw std::invalid_argument("time2 must be positive for the " + 
                                    dist1 + " distribution");
    }
  }
  
  // --- event variable ---
  bool has_event = !event.empty() && data.containElementNamed(event);
  if (!has_time2 && !has_event) {
    throw std::invalid_argument(
        "data must contain the event variable for right censored data");
  }
  std::vector<int> eventn(n);
  if (has_event) {
    if (data.bool_cols.count(event)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
      for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
    } else if (data.int_cols.count(event)) {
      eventn = data.get<int>(event);
    } else if (data.numeric_cols.count(event)) {
      const std::vector<double>& vd = data.get<double>(event);
      for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
    } else {
      throw std::invalid_argument("event variable must be bool, integer or numeric");
    }
    for (double val : eventn) if (val != 0 && val != 1)
      throw std::invalid_argument("event must be 1 or 0 for each observation");
  }
  
  // --- build design matrix zn (n x nvar) column-major FlatMatrix ---
  FlatMatrix zn(n, nvar);
  // intercept
  for (int i = 0; i < n; ++i) zn.data[i] = 1.0;
  // covariates
  for (int j = 0; j < nvar - 1; ++j) {
    const std::string& zj = covariates[j];
    if (!data.containElementNamed(zj)) 
      throw std::invalid_argument("data must contain the variables in covariates");
    double* zn_col = zn.data_ptr() + (j + 1) * n;
    if (data.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
      for (int i = 0; i < n; ++i) zn_col[i] = vb[i] ? 1.0 : 0.0;
    } else if (data.int_cols.count(zj)) {
      const std::vector<int>& vi = data.get<int>(zj);
      for (int i = 0; i < n; ++i) zn_col[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(zj)) {
      const std::vector<double>& vd = data.get<double>(zj);
      std::memcpy(zn_col, vd.data(), n * sizeof(double));
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }
  
  // --- weight and offset ---
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const std::vector<int>& vi = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else {
      throw std::invalid_argument("weight variable must be integer or numeric");
    }
    for (double w : weightn) if (std::isnan(w) || w <= 0.0) 
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    if (data.int_cols.count(offset)) {
      const std::vector<int>& vi = data.get<int>(offset);
      for (int i = 0; i < n; ++i) offsetn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(offset)) {
      offsetn = data.get<double>(offset);
    } else {
      throw std::invalid_argument("offset variable must be integer or numeric");
    }
  }
  
  // --- id mapping ---
  bool has_id = !id.empty() && data.containElementNamed(id);
  std::vector<int> idn(n);
  if (!has_id) {
    std::iota(idn.begin(), idn.end(), 0);
  } else {
    if (data.int_cols.count(id)) {
      auto v = data.get<int>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.numeric_cols.count(id)) {
      auto v = data.get<double>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.string_cols.count(id)) {
      auto v = data.get<std::string>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else throw std::invalid_argument(
        "incorrect type for the id variable in the input data");
  }
  
  // --- exclude observations with missing values ---
  std::vector<unsigned char> sub(n, 1);
  for (int i = 0; i < n; ++i) {
    if (stratumn[i] == INT_MIN || idn[i] == INT_MIN ||
        (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
        std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
      sub[i] = 0;
      continue;
    }
    // check covariates columns
    for (int j = 0; j < nvar - 1; ++j) {
      double v = zn.data[(j + 1) * n + i];
      if (std::isnan(v)) { sub[i] = 0; break; }
    }
  }
  std::vector<int> keep = which(sub);
  if (keep.empty()) throw std::invalid_argument(
      "no observations without missing values");
  
  subset_in_place(stratumn, keep);
  subset_in_place(timen, keep);
  subset_in_place(time2n, keep);
  subset_in_place(eventn, keep);
  subset_in_place(weightn, keep);
  subset_in_place(offsetn, keep);
  subset_in_place(idn, keep);
  subset_in_place_flatmatrix(zn, keep);
  n = static_cast<int>(keep.size());
  
  // --- tstart / tstop and status ---
  std::vector<double> tstart(n), tstop(n);
  if (!has_time2) {
    for (int i = 0; i < n; ++i) {
      tstart[i] = timen[i];
      tstop[i] = (eventn[i] == 1) ? tstart[i] : NaN;
    }
  } else {
    tstart = timen;
    tstop = time2n;
  }
  std::vector<int> status(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) && tstart[i] == tstop[i]) 
      status[i] = 1;
    else if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) && tstart[i] < tstop[i]) 
      status[i] = 3;
    else if (std::isnan(tstart[i]) && !std::isnan(tstop[i])) status[i] = 2;
    else if (!std::isnan(tstart[i]) && std::isnan(tstop[i])) status[i] = 0;
    else status[i] = -1;
  }
  
  // exclude invalid status
  std::vector<unsigned char> good(n, 1);
  for (int i = 0; i < n; ++i) if (status[i] == -1) good[i] = 0;
  std::vector<int> q = which(good);
  int n1 = static_cast<int>(q.size());
  if (n1 == 0) throw std::invalid_argument("no valid records after status filtering");
  
  if (n1 < n) {
    subset_in_place(stratumn, q);
    subset_in_place(tstart, q);
    subset_in_place(tstop, q);
    subset_in_place(status, q);
    subset_in_place(weightn, q);
    subset_in_place(offsetn, q);
    subset_in_place(idn, q);
    subset_in_place_flatmatrix(zn, q);
  }
  
  // --- compute eta (linear predictor) ---
  std::vector<double> eta = offsetn; // initialize with offset
  // add contributions of each coefficient times column
  const double* zptr = zn.data_ptr();
  for (int j = 0; j < nvar; ++j) {
    double b = beta[j];
    if (b == 0.0) continue;
    const double* zcol = zptr + j * n1;
    for (int i = 0; i < n1; ++i) eta[i] += b * zcol[i];
  }
  
  // --- compute sigma per observation ---
  std::vector<double> s(n1, 1.0);
  if (dist_code != 1) {
    for (int i = 0; i < n1; ++i) {
      int k = stratumn[i] + nvar;
      s[i] = std::exp(beta[k]);
    }
  }
  
  // --- map type string to code ---
  int K = 1;
  if (type == "dfbeta" || type == "dfbetas") K = p;
  else if (type == "matrix") K = 6;
  
  int type_code = 0;
  if (type == "response") type_code = 1;
  else if (type == "martingale") type_code = 2;
  else if (type == "deviance") type_code = 3;
  else if (type == "working") type_code = 4;
  else if (type == "dfbeta") type_code = 5;
  else if (type == "dfbetas") type_code = 6;
  else if (type == "ldcase") type_code = 7;
  else if (type == "ldresp") type_code = 8;
  else if (type == "ldshape") type_code = 9;
  else if (type == "matrix") type_code = 10;
  else throw std::invalid_argument("invalid type of residuals: " + type);
  
  // --- prepare result matrix rr (n1 x K) ---
  FlatMatrix rr(n1, K);
  double* rrptr = rr.data_ptr();
  
  // Helper lambdas for rr indexing: rr(i,k) -> rrptr[k * n1 + i]
  
  // --- simple types ---
  if (type_code == 1) { // response
    std::vector<double> yhat0(n1);
    for (int i = 0; i < n1; ++i) {
      if (status[i] == 0 || status[i] == 1) {
        if (dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)
          yhat0[i] = std::log(tstart[i]);
        else yhat0[i] = tstart[i];
      } else if (status[i] == 2) {
        if (dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)
          yhat0[i] = std::log(tstop[i]);
        else yhat0[i] = tstop[i];
      } else { // interval
        if (dist_code == 1 || dist_code == 2) {
          double width = (std::log(tstop[i]) - std::log(tstart[i])) / s[i];
          yhat0[i] = std::log(tstart[i]) - 
            s[i] * std::log(width / (std::exp(width) - 1.0));
        } else if (dist_code == 3 || dist_code == 5) {
          yhat0[i] = 0.5 * (std::log(tstart[i]) + std::log(tstop[i]));
        } else {
          yhat0[i] = 0.5 * (tstart[i] + tstop[i]);
        }
      }
      double val;
      if (dist_code == 1 || dist_code == 2 || dist_code == 3 || dist_code == 5)
        val = std::exp(yhat0[i]) - std::exp(eta[i]);
      else val = yhat0[i] - eta[i];
      rrptr[i] = val;
    }
  } else if (type_code == 2) { // martingale
    if (dist_code == 4 || dist_code == 6)
      throw std::invalid_argument(
          "incorrect distribution for martingale residuals: " + dist1);
    for (int i = 0; i < n1; ++i) {
      if (status[i] == 0 || status[i] == 1) {
        double y = (std::log(tstart[i]) - eta[i]) / s[i];
        double val = 0.0;
        double evt = status[i] == 1 ? 1.0 : 0.0;
        if (dist_code == 1 || dist_code == 2) val = evt - std::exp(y);
        else if (dist_code == 3) val = evt + std::log(boost_pnorm(y,0,1,0));
        else if (dist_code == 5) val = evt + std::log(boost_plogis(y,0,1,0));
        rrptr[i] = val;
      } else {
        rrptr[i] = NaN;
      }
    }
  } else {
    // --- complex types using f_ld_1 ---
    aftparams param = { dist_code, stratumn, tstart, tstop, status, weightn, 
                        offsetn, zn, nstrata };
    ListCpp der = f_ld_1(eta, s, &param);
    std::vector<double> g   = der.get<std::vector<double>>("g");
    std::vector<double> dg  = der.get<std::vector<double>>("dg");
    std::vector<double> ddg = der.get<std::vector<double>>("ddg");
    std::vector<double> ds = der.get<std::vector<double>>("ds");
    std::vector<double> dds = der.get<std::vector<double>>("dds");
    std::vector<double> dsg = der.get<std::vector<double>>("dsg");
    
    if (type_code == 3) { // deviance
      std::vector<double> loglik(n1);
      for (int i = 0; i < n1; ++i) {
        switch (status[i]) {
        case 0: case 2:
          loglik[i] = 0.0;
          break;
        case 1:
          if (dist_code == 1 || dist_code == 2) loglik[i] = -std::log(s[i]) - 1.0;
          else if (dist_code == 3 || dist_code == 4)
            loglik[i] = -std::log(std::sqrt(2.0 * M_PI) * s[i]);
          else loglik[i] = -std::log(4.0 * s[i]);
          break;
        default: { // interval censored
            double width;
            if (dist_code == 1 || dist_code == 2) {
              width = (std::log(tstop[i]) - std::log(tstart[i])) / s[i];
              loglik[i] = - width / (std::exp(width) - 1.0) + 
                std::log(1.0 - std::exp(-width));
            } else if (dist_code == 3) {
              width = (std::log(tstop[i]) - std::log(tstart[i])) / s[i];
              loglik[i] = std::log(2.0 * boost_pnorm(width / 2.0) - 1.0);
            } else if (dist_code == 4) {
              width = (tstop[i] - tstart[i]) / s[i];
              loglik[i] = std::log(2.0 * boost_pnorm(width / 2.0) - 1.0);
            } else if (dist_code == 5) {
              width = (std::log(tstop[i]) - std::log(tstart[i])) / s[i];
              loglik[i] = std::log((std::exp(width / 2.0) - 1.0) / 
                (std::exp(width / 2.0) + 1.0));
            } else {
              width = (tstop[i] - tstart[i]) / s[i];
              loglik[i] = std::log((std::exp(width / 2.0) - 1.0) / 
                (std::exp(width / 2.0) + 1.0));
            }
          }
        }
        double val = -dg[i] / ddg[i];
        double dev = 0.0;
        if (std::isfinite(loglik[i]) && std::isfinite(g[i])) {
          double tmp = 2.0 * (loglik[i] - g[i]);
          dev = (tmp > 0.0) ? std::sqrt(tmp) : 0.0;
          if (val < 0) dev = -dev;
        } else dev = NaN;
        rrptr[i] = dev;
      }
    } else if (type_code == 4) { // working
      for (int i = 0; i < n1; ++i) rrptr[i] = -dg[i] / ddg[i];
    } else if (type_code == 5 || type_code == 6 || type_code == 7) { 
      // dfbeta, dfbetas, ldcase
      // vbeta is p x p FlatMatrix (column-major)
      const double* vptr = vbeta.data_ptr();
      for (int i = 0; i < n1; ++i) {
        // compute score vector
        std::vector<double> score(p, 0.0);
        for (int j = 0; j < nvar; ++j) score[j] = dg[i] * zn(i,j);
        for (int j = nvar; j < p; ++j) 
          score[j] = (stratumn[i] == j - nvar) ? ds[i] : 0.0;
        // resid = score * vbeta  (1 x p) * (p x p) -> vector length p
        std::vector<double> resid(p, 0.0);
        for (int k = 0; k < p; ++k) {
          double sum = 0.0;
          const double* colvk = vptr + k * p; // vbeta(:,k)
          for (int j = 0; j < p; ++j) sum += score[j] * colvk[j];
          resid[k] = sum;
        }
        if (type_code == 6) {
          for (int k = 0; k < p; ++k) {
            double denom = vbeta(k,k);
            if (denom <= 0.0) rrptr[k * n1 + i] = NaN;
            else rrptr[k * n1 + i] = resid[k] / std::sqrt(denom);
          }
        } else if (type_code == 7) {
          double acc = 0.0;
          for (int k = 0; k < p; ++k) acc += resid[k] * score[k];
          rrptr[i] = acc;
        } else {
          for (int k = 0; k < p; ++k) rrptr[k * n1 + i] = resid[k];
        }
      }
    } else if (type_code == 8) { // ldresp
      const double* vptr = vbeta.data_ptr();
      for (int i = 0; i < n1; ++i) {
        std::vector<double> rscore(p, 0.0);
        for (int j = 0; j < nvar; ++j) rscore[j] = -ddg[i] * zn(i,j) * s[i];
        for (int j = nvar; j < p; ++j) 
          rscore[j] = (stratumn[i] == j - nvar) ? -dsg[i] * s[i] : 0.0;
        std::vector<double> temp(p, 0.0);
        for (int k = 0; k < p; ++k) {
          const double* colvk = vptr + k * p;
          double sum = 0.0;
          for (int j = 0; j < p; ++j) sum += rscore[j] * colvk[j];
          temp[k] = sum;
        }
        double acc = 0.0;
        for (int k = 0; k < p; ++k) acc += temp[k] * rscore[k];
        rrptr[i] = acc;
      }
    } else if (type_code == 9) { // ldshape
      const double* vptr = vbeta.data_ptr();
      for (int i = 0; i < n1; ++i) {
        std::vector<double> sscore(p, 0.0);
        for (int j = 0; j < nvar; ++j) sscore[j] = dsg[i] * zn(i,j);
        for (int j = nvar; j < p; ++j) 
          sscore[j] = (stratumn[i] == j - nvar) ? dds[i] : 0.0;
        std::vector<double> temp(p, 0.0);
        for (int k = 0; k < p; ++k) {
          const double* colvk = vptr + k * p;
          double sum = 0.0;
          for (int j = 0; j < p; ++j) sum += sscore[j] * colvk[j];
          temp[k] = sum;
        }
        double acc = 0.0;
        for (int k = 0; k < p; ++k) acc += temp[k] * sscore[k];
        rrptr[i] = acc;
      }
    } else if (type_code == 10) { // matrix
      for (int i = 0; i < n1; ++i) {
        rrptr[i] = g[i];
        rrptr[1 * n1 + i] = dg[i];
        rrptr[2 * n1 + i] = ddg[i];
        rrptr[3 * n1 + i] = ds[i];
        rrptr[4 * n1 + i] = dds[i];
        rrptr[5 * n1 + i] = dsg[i];
      }
    }
  } // end complex types
  
  // --- apply case weights if requested ---
  if (weighted) {
    for (int k = 0; k < K; ++k) {
      double* rrcol = rrptr + k * n1;
      for (int i = 0; i < n1; ++i) {
        rrcol[i] *= weightn[i];
      }
    }
  }
  
  // --- collapse by id if requested ---
  if (collapse) {
    // order by id
    std::vector<int> order = seqcpp(0, n1 - 1);
    std::sort(order.begin(), order.end(), [&](int i, int j){ 
      return idn[i] < idn[j]; });
    std::vector<int> id1 = subset(idn, order);
    std::vector<int> idx(1,0);
    for (int i = 1; i < n1; ++i) if (id1[i] != id1[i-1]) idx.push_back(i);
    int nids = static_cast<int>(idx.size());
    idx.push_back(n1);
    
    FlatMatrix rr1(nids, K);
    double* rr1ptr = rr1.data_ptr();
    for (int i = 0; i < nids; ++i) {
      for (int k = 0; k < K; ++k) {
        double acc = 0.0;
        const double* rrcol = rrptr + k * n1;
        for (int j = idx[i]; j < idx[i+1]; ++j) {
          acc += rrcol[order[j]];
        }
        rr1ptr[ k * nids + i ] = acc;
      }
    }
    return rr1;
  }
  
  return rr;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix residuals_liferegRcpp(
    const std::vector<double>& beta,
    const Rcpp::NumericMatrix& vbeta,
    const Rcpp::DataFrame& data,
    const std::vector<std::string>& stratum,
    const std::string& time,
    const std::string& time2,
    const std::string& event,
    const std::vector<std::string>& covariates,
    const std::string& weight,
    const std::string& offset,
    const std::string& id,
    const std::string& dist,
    const std::string& type,
    const bool collapse,
    const bool weighted) {
  
  auto dfcpp = convertRDataFrameToCpp(data);
  auto vbetacpp = flatmatrix_from_Rmatrix(vbeta);
  
  auto rrcpp = residuals_liferegcpp(
    beta, vbetacpp, dfcpp, stratum, time, time2, event, covariates, 
    weight, offset, id, dist, type, collapse, weighted
  );
  
  return Rcpp::wrap(rrcpp);
}

struct coxparams {
  int nused;
  std::vector<int> strata;
  std::vector<double> tstart;
  std::vector<double> tstop;
  std::vector<int> event;
  std::vector<double> weight;
  std::vector<double> offset;
  FlatMatrix z;
  std::vector<int> order1;
  int method; // 1: breslow, 2: efron
};

// all-in-one function for log-likelihood, score, and information matrix
// for the Cox model with or without Firth's correction
ListCpp f_der_2(int p, const std::vector<double>& par, void* ex, bool firth) {
  coxparams* param = (coxparams*) ex;
  const int n = param->tstop.size();
  const int nused = param->nused;
  const int method = param->method;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& event = param->event; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const std::vector<int>& order1 = param->order1;
  const FlatMatrix& z = param->z;
  
  // Precompute eta and risk = exp(eta)
  std::vector<double> eta(nused);
  for (int person = 0; person < nused; ++person) {
    eta[person] = offset[person];
  }
  if (p > 0) {
    for (int i = 0; i < p; ++i) {
      double beta = par[i];
      if (beta == 0.0) continue;
      const double* zcol = z.data_ptr() + i * n;
      for (int person = 0; person < nused; ++person) {
        eta[person] += beta * zcol[person];
      }
    }
  }
  std::vector<double> risk(nused);
  for (int person = 0; person < nused; ++person) {
    risk[person] = std::exp(eta[person]);
  }
  
  double loglik = 0.0;        // log-likelihood
  std::vector<double> u(p);   // score vector
  FlatMatrix imat(p,p);       // information matrix
  FlatArray dimat(p,p,p);     // tensor for third order derivatives
  std::vector<double> a(p);   // s1(beta,k,t)
  std::vector<double> a2(p);  // sum of w*exp(zbeta)*z for the deaths
  FlatMatrix cmat(p,p);       // s2(beta,k,t)
  FlatMatrix cmat2(p,p);      // sum of w*exp(zbeta)*z*z' for the deaths
  FlatArray dmat(p,p,p);      // q2(beta,k,t)
  FlatArray dmat2(p,p,p);     // sum of w*exp(zbeta)*z*z*z' for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  double ndead = 0.0;         // number of deaths at this time point
  
  int istrata = strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    // Reset when entering a new stratum
    if (strata[person] != istrata) {
      istrata = strata[person];
      i1 = person;
      denom = 0.0;
      
      for (int i = 0; i < p; ++i) {
        a[i] = 0.0;
      }
      
      for (int j = 0; j < p; ++j) {
        for (int i = j; i < p; ++i) {
          cmat(i,j) = 0.0;
        }
      }
      
      if (firth) {
        for (int k = 0; k < p; ++k) {
          for (int j = k; j < p; ++j) {
            for (int i = j; i < p; ++i) {
              dmat(i,j,k) = 0.0;
            }
          }
        }
      }
    }
    
    const double dtime = tstop[person];
    
    // Process all persons tied at this dtime
    for (; person < nused && tstop[person] == dtime && 
         strata[person] == istrata; ++person) {
      
      const double w = weight[person];
      const double r = w * risk[person];
      
      if (event[person] == 0) {
        denom += r;
        
        for (int i = 0; i < p; ++i) {
          a[i] += r * z(person,i);
        }
        
        for (int j = 0; j < p; ++j) {
          const double zj = z(person,j);
          for (int i = j; i < p; ++i) {
            cmat(i,j) += r * z(person,i) * zj;  
          }
        }
        
        if (firth) {
          for (int k = 0; k < p; ++k) {
            const double zk = z(person,k);
            for (int j = k; j < p; ++j) {
              const double zj = z(person,j);
              for (int i = j; i < p; ++i) {
                dmat(i,j,k) += r * z(person,i) * zj * zk;
              }
            }
          }
        }
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
        loglik += w * eta[person];
        
        for (int i = 0; i < p; ++i) {
          const double zi = z(person,i);
          a2[i] += r * zi;
          u[i] += w * zi;
        }
        
        for (int j = 0; j < p; ++j) {
          const double zj = z(person,j);
          for (int i = j; i < p; ++i) {
            cmat2(i,j) += r * z(person,i) * zj;  
          }
        }
        
        if (firth) {
          for (int k = 0; k < p; ++k) {
            const double zk = z(person,k);
            for (int j = k; j < p; ++j) {
              const double zj = z(person,j);
              for (int i = j; i < p; ++i) {
                dmat2(i,j,k) += r * z(person,i) * zj * zk;
              }
            }
          }
        }
      }
    }
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = order1[i1];
      if (tstart[p1] < dtime || strata[p1] != istrata) break;
      
      const double r = weight[p1] * risk[p1];
      denom -= r;
      
      for (int i = 0; i < p; ++i) {
        a[i] -= r * z(p1,i);
      }
      
      for (int j = 0; j < p; ++j) {
        const double zj = z(p1,j);
        for (int i = j; i < p; ++i) {
          cmat(i,j) -= r * z(p1,i) * zj;  
        }
      }
      
      if (firth) {
        for (int k = 0; k < p; ++k) {
          const double zk = z(p1,k);
          for (int j = k; j < p; ++j) {
            const double zj = z(p1,j);
            for (int i = j; i < p; ++i) {
              dmat(i,j,k) -= r * z(p1,i) * zj * zk;
            }
          }
        }
      }
    }
    
    // Add contributions for deaths at this time
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        loglik -= deadwt * std::log(denom);
        
        std::vector<double> xbar(p);
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          xbar[i] = a[i] / denom;
          u[i] -= deadwt * xbar[i];
        }
        
        for (int j = 0; j < p; ++j) {
          for (int i = j; i < p; ++i) {
            cmat(i,j) += cmat2(i,j);
            imat(i,j) += deadwt * (cmat(i,j) - xbar[i] * a[j]) / denom;
          }
        }
        
        if (firth) {
          for (int k = 0; k < p; ++k) {
            for (int j = k; j < p; ++j) {
              for (int i = j; i < p; ++i) {
                dmat(i,j,k) += dmat2(i,j,k);
                dimat(i,j,k) += deadwt * (dmat(i,j,k) -
                  (cmat(i,j) * a[k] + cmat(i,k) * a[j] + cmat(j,k) * a[i]) / denom +
                  2.0 * a[i] * a[j] * a[k] / (denom * denom)) / denom;
              }
            }
          }
        }
      } else { // Efron method
        const double meanwt = deadwt / ndead;
        const double increment = denom2 / ndead;
        for (int l = 0; l < ndead; ++l) {
          denom += increment;
          loglik -= meanwt * std::log(denom);
          
          std::vector<double> xbar(p);
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            xbar[i] = a[i] / denom;
            u[i] -= meanwt * xbar[i];
          }
          
          for (int j = 0; j < p; ++j) {
            for (int i = j; i < p; ++i) {
              cmat(i,j) += cmat2(i,j) / ndead;
              imat(i,j) += meanwt * (cmat(i,j) - xbar[i] * a[j]) / denom;
            }
          }
          
          if (firth) {
            for (int k = 0; k < p; ++k) {
              for (int j = k; j < p; ++j) {
                for (int i = j; i < p; ++i) {
                  dmat(i,j,k) += dmat2(i,j,k) / ndead;
                  dimat(i,j,k) += meanwt * (dmat(i,j,k) -
                    (cmat(i,j) * a[k] + cmat(i,k) * a[j] + cmat(j,k) * a[i]) / denom +
                    2.0 * a[i] * a[j] * a[k] / (denom * denom)) / denom;
                }
              }
            }
          }
        }
      }
      
      // Reset after processing deaths
      ndead = deadwt = denom2 = 0.0;
      
      for (int i = 0; i < p; ++i) {
        a2[i] = 0.0;
      }
      
      for (int j = 0; j < p; ++j) {
        for (int i = j; i < p; ++i) {
          cmat2(i,j) = 0.0;
        }
      }
      
      if (firth) {
        for (int k = 0; k < p; ++k) {
          for (int j = k; j < p; ++j) {
            for (int i = j; i < p; ++i) {
              dmat2(i,j,k) = 0.0;
            }
          }
        }
      }
    }
  }
  
  // fill the symmetric elements of the information matrix
  for (int j = 1; j < p; ++j)
    for (int i = 0; i < j; ++i)
      imat(i,j) = imat(j,i);
  
  // fill the symmetric elements of the tensor array
  if (firth) {
    for (int k = 0; k < p-1; ++k)
      for (int j = k+1; j < p; ++j) 
        for (int i = k; i < j; ++i) 
          dimat(i,j,k) = dimat(j,i,k);
    
    for (int k = 1; k < p; ++k) 
      for (int j = 0; j < k; ++j)
        for (int i = j; i < p; ++i) 
          dimat(i,j,k) = dimat(i,k,j);
    
    for (int k = 1; k < p; ++k) {
      for (int j = 1; j < p; ++j) {
        int l = std::min(j,k);
        for (int i = 0; i < l; ++i) 
          dimat(i,j,k) = dimat(k,j,i);
      }
    }
  }
  
  ListCpp result;
  
  // Firth adjustment
  if (p > 0 && firth) {
    // obtain the determinant of information matrix
    FlatMatrix imat0 = imat;
    cholesky2(imat0, p);
    double* base = imat0.data_ptr();
    
    double v = 0.0;
    for (int i = 0; i < p; ++i) {
      v += std::log(imat0(i,i));
    }
    
    // penalized log-likelihood adjustment
    double penloglik = loglik + 0.5 * v;
    
    // compute the bias adjustment to the score function
    FlatMatrix y(p,p);
    std::vector<double> g(p);
    double* yptr = y.data_ptr();
    
    for (int k = 0; k < p; ++k) {
      // partial derivative of the information matrix w.r.t. beta[k]
      for (int j = 0; j < p; ++j) {
        for (int i = 0; i < p; ++i) {
          y(i,j) = dimat(i,j,k);
        }
      }
      
      // solve(imat, y)
      for (int h = 0; h < p; ++h) {
        double* yh = yptr + h * p;
        
        for (int j = 0; j < p - 1; ++j) {
          double yjh = yh[j];
          if (yjh == 0.0) continue;
          double* col_j = base + j * p;
          for (int i = j + 1; i < p; ++i) {
            yh[i] -= yjh * col_j[i];
          }
        }
        for (int i = p - 1; i >= 0; --i) {
          double* col_i = base + i * p;
          double diag = col_i[i];
          if (diag == 0.0) yh[i] = 0.0;
          else {
            double temp = yh[i] / diag;
            for (int j = i + 1; j < p; ++j)
              temp -= yh[j] * col_i[j];
            yh[i] = temp;
          }
        }
      }
      
      // trace
      for (int i = 0; i < p; ++i) g[k] += y(i,i);
      
      g[k] = u[k] + 0.5 * g[k];
    }
    
    result.push_back(penloglik, "loglik");
    result.push_back(std::move(g), "score");
    result.push_back(std::move(imat), "imat");
    result.push_back(loglik, "regloglik");
    result.push_back(std::move(u), "regscore");
  } else {
    result.push_back(loglik, "loglik");
    if (p > 0) {
      result.push_back(std::move(u), "score");
      result.push_back(std::move(imat), "imat");
    }
  }
  
  return result;
}


// underlying optimization algorithm for phreg
ListCpp phregloop(int p, const std::vector<double>& par, void *ex,
                  int maxiter, double eps, bool firth,
                  const std::vector<int>& colfit, int ncolfit) {
  coxparams *param = (coxparams *) ex;
  
  int iter = 0, halving = 0;
  bool fail = false;
  
  std::vector<double> beta = par;
  std::vector<double> newbeta(p);
  double loglik = 0.0, newlk = 0.0;
  std::vector<double> u(p);
  std::vector<double> uu1(ncolfit);
  double* u1 = uu1.data();
  FlatMatrix imat(p,p);
  FlatMatrix imat1(ncolfit, ncolfit);
  
  // --- first step ---
  ListCpp der = f_der_2(p, beta, param, firth);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  imat = der.get<FlatMatrix>("imat");
  
  for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
  
  for (int j = 0; j < ncolfit; ++j)
    for (int i = 0; i < ncolfit; ++i)
      imat1(i,j) = imat(colfit[i], colfit[j]);
  
  cholesky2(imat1, ncolfit);
  chsolve2(imat1, ncolfit, u1);
  
  std::fill(u.begin(), u.end(), 0.0);
  for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + u[i];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_2(p, newbeta, param, firth);
    newlk = der.get<double>("loglik");
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && halving == 0 && std::fabs(1 - (loglik / newlk)) < eps) break;
    
    if (fail || newlk < loglik) {
      ++halving; // adjust step size if likelihood decreases
      for (int i = 0; i < p; ++i) {
        newbeta[i] = 0.5 * (beta[i] + newbeta[i]);
      }
      continue;
    }
    
    // --- update ---
    halving = 0;
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<FlatMatrix>("imat");
    
    for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
    
    for (int j = 0; j < ncolfit; ++j)
      for (int i = 0; i < ncolfit; ++i)
        imat1(i,j) = imat(colfit[i], colfit[j]);
    
    cholesky2(imat1, ncolfit);
    chsolve2(imat1, ncolfit, u1);
    
    std::fill(u.begin(), u.end(), 0.0);
    for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
    for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + u[i];
  }
  
  if (iter == maxiter) fail = true;
  
  // --- final variance calculation ---
  imat = der.get<FlatMatrix>("imat");
  for (int j = 0; j < ncolfit; ++j)
    for (int i = 0; i < ncolfit; ++i)
      imat1(i,j) = imat(colfit[i], colfit[j]);
  
  FlatMatrix var1 = invsympd(imat1, ncolfit);
  FlatMatrix var(p, p);
  for (int j = 0; j < ncolfit; ++j)
    for (int i = 0; i < ncolfit; ++i)
      var(colfit[i], colfit[j]) = var1(i,j);
  
  ListCpp result;
  result.push_back(std::move(newbeta), "coef");
  result.push_back(iter, "iter");
  result.push_back(std::move(var), "var");
  result.push_back(newlk, "loglik");
  result.push_back(fail, "fail");
  
  if (firth) {
    double regloglik = der.get<double>("regloglik");
    result.push_back(regloglik, "regloglik");
  }
  
  return result;
}


// confidence limit of profile likelihood method
double phregplloop(int p, const std::vector<double>& par, void *ex,
                   int maxiter, double eps, bool firth,
                   int k, int direction, double l0) {
  coxparams *param = (coxparams *) ex;
  int iter;
  bool fail = false;
  
  std::vector<double> beta = par;
  std::vector<double> newbeta(p);
  double loglik = 0.0, newlk = 0.0;
  std::vector<double> u(p);
  std::vector<double> delta(p);
  FlatMatrix imat(p, p);
  FlatMatrix v(p, p);
  
  ListCpp der = f_der_2(p, beta, param, firth);
  loglik = der.get<double>("loglik");
  u = der.get<std::vector<double>>("score");
  imat = der.get<FlatMatrix>("imat");
  v = invsympd(imat, p);
  
  // compute w = - u^T v u
  double w = -quadsym(u, v);
  double underroot = -2 * (l0 - loglik + 0.5 * w) / v(k, k);
  double lambda = underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot);
  u[k] += lambda;
  delta = mat_vec_mult(v, u);
  for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_2(p, newbeta, param, firth);
    newlk = der.get<double>("loglik");
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && std::fabs(newlk - l0) < eps && w < eps) break;
    beta = newbeta;
    loglik = newlk;
    u = der.get<std::vector<double>>("score");
    imat = der.get<FlatMatrix>("imat");
    v = invsympd(imat, p);
    w = -quadsym(u, v);
    underroot = -2 * (l0 - newlk + 0.5 * w) / v(k, k);
    lambda = underroot < 0.0 ? 0.0 : direction * std::sqrt(underroot);
    u[k] += lambda;
    delta = mat_vec_mult(v, u);
    for (int i = 0; i < p; ++i) newbeta[i] = beta[i] + delta[i];
  }
  
  if (iter == maxiter) fail = true;
  if (fail) thread_utils::push_thread_warning("phregplloop did not converge.");
  
  return newbeta[k];
}


// baseline hazard estimates
ListCpp f_basehaz(int p, const std::vector<double>& par, void *ex) {
  coxparams *param = (coxparams *) ex;
  const int n = param->tstop.size();
  const int nused = param->nused;
  const int method = param->method;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& event = param->event; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const std::vector<int>& order1 = param->order1;
  const FlatMatrix& z = param->z;
  
  // Precompute eta and risk = exp(eta)
  std::vector<double> eta(nused);
  for (int person = 0; person < nused; ++person) {
    eta[person] = offset[person];
  }
  if (p > 0) {
    for (int i = 0; i < p; ++i) {
      double beta = par[i];
      if (beta == 0.0) continue;
      const double* zcol = z.data_ptr() + i * n;
      for (int person = 0; person < nused; ++person) {
        eta[person] += beta * zcol[person];
      }
    }
  }
  std::vector<double> risk(nused);
  for (int person = 0; person < nused; ++person) {
    risk[person] = std::exp(eta[person]);
  }
  
  std::vector<double> a(p);   // s1(beta,k,t)
  std::vector<double> a2(p);  // sum of w*exp(zbeta)*z for the deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for the deaths
  double natrisk = 0;         // number at risk at this time point
  double ndead = 0;           // number of deaths at this time point
  double ncens = 0;           // number of censored at this time point
  
  // locate the first observation within each stratum
  std::vector<int> istratum(1,0);
  for (int i = 1; i < nused; ++i) {
    if (strata[i] != strata[i-1]) {
      istratum.push_back(i);
    }
  }
  
  int nstrata = static_cast<int>(istratum.size());
  istratum.push_back(nused);
  
  // add time 0 to each stratum
  int J = nstrata;
  for (int i = 0; i < nstrata; ++i) {
    std::vector<double> utime = subset(tstop, istratum[i], istratum[i+1]);
    utime = unique_sorted(utime);
    J += static_cast<int>(utime.size());
  }
  
  std::vector<int> stratum(J);
  std::vector<double> time(J), nrisk(J), nevent(J), ncensor(J), haz(J), varhaz(J);
  FlatMatrix gradhaz(J,p);
  
  int istrata = strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  int j = J;  // index the unique time in ascending order
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    if (strata[person] != istrata) { // hit a new stratum
      // add time 0 at the start of a new stratum
      j--;
      stratum[j] = istrata;
      time[j] = 0.0;
      nrisk[j] = natrisk;
      
      istrata = strata[person]; // reset temporary variables
      i1 = person;
      natrisk = 0;
      denom = 0.0;
      std::fill(a.begin(), a.end(), 0.0);
    }
    
    const double dtime = tstop[person];
    
    // Process all persons tied at this dtime
    bool first = true;
    for (; person < nused && tstop[person] == dtime && 
         strata[person] == istrata; ++person) {
      
      if (first) { // first incidence at this time
        j--;
        stratum[j] = strata[person];
        time[j] = dtime;
        first = false;
      }
      
      const double w = weight[person];
      const double r = w * risk[person];
      
      ++natrisk;
      if (event[person] == 0) {
        ++ncens;
        denom += r;
        for (int i = 0; i < p; ++i) {
          a[i] += r * z(person,i);
        }
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
        for (int i = 0; i < p; ++i) {
          a2[i] += r * z(person,i);
        }
      }
    }
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = order1[i1];
      if (tstart[p1] < dtime || strata[p1] != istrata) break;
      
      const double r = weight[p1] * risk[p1];
      
      natrisk--;
      denom -= r;
      for (int i = 0; i < p; ++i) {
        a[i] -= r * z(p1,i);
      }
    }
    
    // Add contributions for deaths at this time
    nrisk[j] = natrisk;
    nevent[j] = ndead;
    ncensor[j] = ncens;
    ncens = 0; // reset for the next time point
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        const double temp = deadwt / denom;
        haz[j] = temp;
        varhaz[j] = temp / denom;
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          gradhaz(j,i) = temp * a[i] / denom;
        }
      } else { // Efron method
        const double meanwt = deadwt / ndead;
        for (int k = 0; k < ndead; ++k) {
          denom += denom2 / ndead;
          const double temp = meanwt / denom;
          haz[j] += temp;
          varhaz[j] += temp / denom;
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            gradhaz(j,i) += temp * a[i] / denom;
          }
        }
      }
      
      // reset for the next death time
      ndead = deadwt = denom2 = 0.0;
      std::fill(a2.begin(), a2.end(), 0.0);
    }
  }
  
  // add time 0 for the first stratum
  stratum[0] = istrata;
  time[0] = 0.0;
  nrisk[0] = natrisk;
  
  ListCpp result;
  result.push_back(std::move(stratum), "stratum");
  result.push_back(std::move(time), "time");
  result.push_back(std::move(nrisk), "nrisk");
  result.push_back(std::move(nevent), "nevent");
  result.push_back(std::move(ncensor), "ncensor");
  result.push_back(std::move(haz), "haz");
  result.push_back(std::move(varhaz), "varhaz");
  
  if (p > 0) result.push_back(std::move(gradhaz), "gradhaz");
  
  return result;
}


// martingale residuals
std::vector<double> f_resmart(int p, const std::vector<double>& par, void *ex) {
  coxparams *param = (coxparams *) ex;
  const int n = param->tstop.size();
  const int nused = param->nused;
  const int method = param->method;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& event = param->event; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const std::vector<int>& order1 = param->order1;
  const FlatMatrix& z = param->z;
  
  // Precompute eta and risk = exp(eta)
  std::vector<double> eta(nused);
  for (int person = 0; person < nused; ++person) {
    eta[person] = offset[person];
  }
  if (p > 0) {
    for (int i = 0; i < p; ++i) {
      double beta = par[i];
      if (beta == 0.0) continue;
      const double* zcol = z.data_ptr() + i * n;
      for (int person = 0; person < nused; ++person) {
        eta[person] += beta * zcol[person];
      }
    }
  }
  std::vector<double> risk(nused);
  for (int person = 0; person < nused; ++person) {
    risk[person] = std::exp(eta[person]);
  }
  
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  double ndead = 0.0;         // number of deaths at this time point
  
  // initialize the residuals to the event indicators
  std::vector<double> resid(n);
  for (int person = 0; person < nused; ++person) {
    resid[person] = event[person];
  }
  
  int istrata = strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  int j0 = 0; // first person in the stratum
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    if (strata[person] != istrata) { // hit a new stratum
      istrata = strata[person]; // reset temporary variables
      i1 = person;
      j0 = person;
      denom = 0.0;
    }
    
    const double dtime = tstop[person];
    
    // process all persons tied at this dtime
    int j1 = person;   // first person in the stratum with the tied time
    for (; person < nused && tstop[person] == dtime && 
         strata[person] == istrata; ++person) {
      
      const double w = weight[person];
      const double r = w * risk[person];
      
      if (event[person] == 0) {
        denom += r;
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
      }
    }
    
    int j2 = person - 1; // last person in the stratum with the tied time
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = order1[i1];
      if (tstart[p1] < dtime || strata[p1] != istrata) break;
      denom -= weight[p1] * risk[p1];
    }
    
    // Add contributions for deaths at this time
    if (ndead > 0) {
      denom += denom2;
      
      for (int j = j0; j <= j2; ++j) {
        if (tstart[j] < dtime) {
          double hazard;
          if (method == 0 || ndead == 1) {
            hazard = deadwt / denom;
          } else {
            hazard = 0.0;
            const double meanwt = deadwt / ndead;
            if (j < j1 || event[j] == 0) {
              for (int i = 0; i < ndead; ++i) {
                hazard += meanwt / (denom - i / ndead * denom2);
              }
            } else {
              for (int i = 0; i < ndead; ++i) {
                hazard += (1 - i / ndead) * meanwt / (denom - i / ndead * denom2);
              }
            }
          }
          resid[j] -= hazard * risk[j];
        }
      }
      
      // reset for the next death time
      ndead = deadwt = denom2 = 0.0;
    }
  }
  
  return resid;
}


// score residual matrix
FlatMatrix f_ressco_2(int p, const std::vector<double>& par, void *ex) {
  coxparams *param = (coxparams *) ex;
  const int n = param->tstop.size();
  const int nused = param->nused;
  const int method = param->method;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& event = param->event; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const std::vector<int>& order1 = param->order1;
  const FlatMatrix& z = param->z;
  
  // Precompute eta and risk = exp(eta)
  std::vector<double> eta(nused);
  for (int person = 0; person < nused; ++person) {
    eta[person] = offset[person];
  }
  if (p > 0) {
    for (int i = 0; i < p; ++i) {
      double beta = par[i];
      if (beta == 0.0) continue;
      const double* zcol = z.data_ptr() + i * n;
      for (int person = 0; person < nused; ++person) {
        eta[person] += beta * zcol[person];
      }
    }
  }
  std::vector<double> risk(nused);
  for (int person = 0; person < nused; ++person) {
    risk[person] = std::exp(eta[person]);
  }
  
  FlatMatrix resid(n,p);      // residual matrix
  std::vector<double> a(p);   // s1(beta,k,t)
  std::vector<double> a2(p);  // sum of w*exp(zbeta)*z for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  double ndead = 0.0;         // number of deaths at this time point
  double cumhaz = 0.0;        // cumulative hazard
  
  std::vector<double> xhaz(p), mh1(p), mh2(p), mh3(p); // temp vectors
  
  int istrata = strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    // Reset when entering a new stratum
    if (strata[person] != istrata) {
      istrata = strata[person];
      
      // first obs of a new stratum, finish off the prior stratum
      for (; i1 < nused && order1[i1] < person; ++i1) {
        const int p1 = order1[i1];
        for (int i = 0; i < p; ++i) {
          resid(p1,i) -= risk[p1] * (z(p1,i) * cumhaz - xhaz[i]);
        }
      }
      
      denom = 0.0; // reset temporary variables
      cumhaz = 0.0;
      std::fill(a.begin(), a.end(), 0.0);
      std::fill(xhaz.begin(), xhaz.end(), 0.0);
    }
    
    const double dtime = tstop[person];
    
    // process all persons tied at this dtime
    for (; person < nused && tstop[person] == dtime && 
         strata[person] == istrata; ++person) {
      
      // initialize residuals to score[i] * (x[i] * cumhaz - xhaz), before
      // updating cumhaz and xhaz
      for (int i = 0; i < p; ++i) {
        resid(person,i) = risk[person] * (z(person,i) * cumhaz - xhaz[i]);
      }
      
      const double w = weight[person];
      const double r = w * risk[person];
      
      if (event[person] == 0) {
        denom += r;
        for (int i = 0; i < p; ++i) {
          a[i] += r * z(person,i);
        }
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
        for (int i = 0; i < p; ++i) {
          a2[i] += r * z(person,i);
        }
      }
    }
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = order1[i1];
      if (tstart[p1] < dtime || strata[p1] != istrata) break;
      
      const double r = weight[p1] * risk[p1];
      denom -= r;
      for (int i = 0; i < p; ++i) {
        // finish the residual by subtracting score[i] * (x[i] * cumhaz - xhaz)
        resid(p1,i) -= risk[p1] * (z(p1,i) * cumhaz - xhaz[i]);
        a[i] -= r * z(p1,i);
      }
    }
    
    // Add contributions for deaths at this time
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        const double hazard = deadwt / denom;
        cumhaz += hazard;
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          const double xbar = a[i] / denom;
          xhaz[i] += xbar * hazard;
          for (int j = person - 1; j >= person - ndead; j--) {
            resid(j,i) += z(j,i) - xbar;
          }
        }
      } else {  // Efron method
        for (int i = 0; i < p; ++i) {
          mh1[i] = 0.0;
          mh2[i] = 0.0;
          mh3[i] = 0.0;
        }
        
        const double meanwt = deadwt / ndead;
        const double increment = denom2 / ndead;
        
        for (int k = 0; k < ndead; ++k) {
          denom += increment;
          const double hazard = meanwt / denom;
          cumhaz += hazard;
          const double downwt = (ndead - k - 1.0) / ndead;
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            const double xbar = a[i] / denom;
            xhaz[i] += xbar * hazard;
            mh1[i]  += hazard * downwt;
            mh2[i]  += xbar * hazard * downwt;
            mh3[i]  += xbar / ndead;
          }
        }
        
        for (int i = 0; i < p; ++i) {
          for (int j = person-1; j >= person - ndead; j--) {
            resid(j,i) += (z(j,i) - mh3[i]) +
              risk[j] * (z(j,i) * mh1[i] - mh2[i]);
          }
        }
      }
      
      // Reset after processing deaths
      ndead = deadwt = denom2 = 0.0;
      std::fill(a2.begin(), a2.end(), 0.0);
    }
  }
  
  // finish those remaining in the final stratum
  for (; i1 < nused; ++i1) {
    const int p1 = order1[i1];
    for (int i = 0; i < p; ++i)
      resid(p1,i) -= risk[p1] * (z(p1,i) * cumhaz - xhaz[i]);
  }
  
  return resid;
}


// main function for phreg
ListCpp phregcpp(const DataFrameCpp& data,
                 const std::vector<std::string>& stratum,
                 const std::string& time,
                 const std::string& time2,
                 const std::string& event,
                 const std::vector<std::string>& covariates,
                 const std::string& weight,
                 const std::string& offset,
                 const std::string& id,
                 const std::string& ties,
                 const std::vector<double>& init,
                 const bool robust,
                 const bool est_basehaz,
                 const bool est_resid,
                 const bool firth,
                 const bool plci,
                 const double alpha,
                 const int maxiter,
                 const double eps) {
  
  int n = data.nrows();
  int p = static_cast<int>(covariates.size());
  if (p == 1 && covariates[0] == "") p = 0;
  
  // --- handle strata (bygroup) ---
  bool has_stratum = false;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
  }
  
  // --- time / time2 existence and checks ---
  if (!data.containElementNamed(time)) 
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(timen[i]) && timen[i] < 0.0)
      throw std::invalid_argument("time must be nonnegative");
  }
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    if (data.int_cols.count(time2)) {
      const std::vector<int>& vi = data.get<int>(time2);
      for (int i = 0; i < n; ++i) time2n[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(time2)) {
      time2n = data.get<double>(time2);
    } else {
      throw std::invalid_argument("time2 variable must be integer or numeric");
    }    
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(timen[i]) && !std::isnan(time2n[i]) && time2n[i] <= timen[i])
        throw std::invalid_argument("time2 must be greater than time");
    }
  }
  
  // --- event variable ---
  bool has_event = !event.empty() && data.containElementNamed(event);
  if (!has_time2 && !has_event) {
    throw std::invalid_argument(
        "data must contain the event variable for right censored data"); 
  }
  std::vector<int> eventn(n);
  if (has_event) {
    if (data.bool_cols.count(event)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
      for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
    } else if (data.int_cols.count(event)) {
      eventn = data.get<int>(event);
    } else if (data.numeric_cols.count(event)) {
      const std::vector<double>& vd = data.get<double>(event);
      for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
    } else {
      throw std::invalid_argument("event variable must be bool, integer or numeric");
    }
    for (double val : eventn) if (val != 0 && val != 1)
      throw std::invalid_argument("event must be 1 or 0 for each observation");
  }
  
  // --- build design matrix zn (n x p) column-major FlatMatrix ---
  FlatMatrix zn(n, p);
  if (p > 0) {
    for (int j = 0; j < p; ++j) {
      const std::string& zj = covariates[j];
      if (!data.containElementNamed(zj)) 
        throw std::invalid_argument("data must contain the variables in covariates");
      double* zn_col = zn.data_ptr() + j * n;
      if (data.bool_cols.count(zj)) {
        const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
        for (int i = 0; i < n; ++i) zn_col[i] = vb[i] ? 1.0 : 0.0;
      } else if (data.int_cols.count(zj)) {
        const std::vector<int>& vi = data.get<int>(zj);
        for (int i = 0; i < n; ++i) zn_col[i] = static_cast<double>(vi[i]);
      } else if (data.numeric_cols.count(zj)) {
        const std::vector<double>& vd = data.get<double>(zj);
        std::memcpy(zn_col, vd.data(), n * sizeof(double));
      } else {
        throw std::invalid_argument("covariates must be bool, integer or numeric");
      }
    }
  }
  
  // --- weight and offset ---
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const std::vector<int>& vi = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else {
      throw std::invalid_argument("weight variable must be integer or numeric");
    }
    for (double w : weightn) if (std::isnan(w) || w <= 0.0) 
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    if (data.int_cols.count(offset)) {
      const std::vector<int>& vi = data.get<int>(offset);
      for (int i = 0; i < n; ++i) offsetn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(offset)) {
      offsetn = data.get<double>(offset);
    } else {
      throw std::invalid_argument("offset variable must be integer or numeric");
    }
  }
  
  // --- id mapping ---
  bool has_id = !id.empty() && data.containElementNamed(id);
  std::vector<int> idn(n);
  if (!has_id) {
    std::iota(idn.begin(), idn.end(), 0);
  } else {
    if (data.int_cols.count(id)) {
      auto v = data.get<int>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.numeric_cols.count(id)) {
      auto v = data.get<double>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.string_cols.count(id)) {
      auto v = data.get<std::string>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else throw std::invalid_argument(
        "incorrect type for the id variable in the input data");
  }
  
  if (robust && has_time2 && !has_id) {
    throw std::invalid_argument(
        "id is needed for counting process data with robust variance");
  }
  
  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  int method = meth == "efron" ? 1 : 0;
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  
  // exclude observations with missing values
  std::vector<unsigned char> sub(n,1);
  for (int i = 0; i < n; ++i) {
    if (stratumn[i] == INT_MIN || idn[i] == INT_MIN ||
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) ||
        eventn[i] == INT_MIN || std::isnan(weightn[i]) ||
        std::isnan(offsetn[i])) {
      sub[i] = 0;
      continue;
    }
    for (int j = 0; j < p; ++j) {
      if (std::isnan(zn(i,j))) { sub[i] = 0; break; }
    }
  }
  std::vector<int> keep = which(sub);
  if (keep.empty()) 
    throw std::invalid_argument("no observations without missing values");
  
  subset_in_place(stratumn, keep);
  subset_in_place(tstartn, keep);
  subset_in_place(tstopn, keep);
  subset_in_place(eventn, keep);
  subset_in_place(weightn, keep);
  subset_in_place(offsetn, keep);
  subset_in_place(idn, keep);
  if (p > 0) subset_in_place_flatmatrix(zn, keep);
  n = keep.size();
  
  // sumstat data set
  double loglik0, loglik1;
  double regloglik0, regloglik1;
  double scoretest;
  int niter;
  bool fail;
  
  // parest data set
  std::vector<std::string> par(p);
  std::vector<double> b(p), seb(p), rseb(p);
  std::vector<double> z(p), expbeta(p);
  FlatMatrix vb(p, p), rvb(p, p);
  std::vector<double> lb(p), ub(p), prob(p);
  std::vector<std::string> clparm(p);
  
  // baseline hazards data set
  std::vector<int> stratumn1 = unique_sorted(stratumn);
  int nstrata = static_cast<int>(stratumn1.size());
  int N = n + nstrata; // add a time 0 row for each stratum
  std::vector<int> dstratum;
  std::vector<double> dtime, dnrisk, dnevent, dncensor;
  std::vector<double> dhaz, dvarhaz;
  FlatMatrix dgradhaz(N,p);
  dstratum.reserve(N); 
  dtime.reserve(N); dnrisk.reserve(N); dnevent.reserve(N); dncensor.reserve(N);
  dhaz.reserve(N); dvarhaz.reserve(N);
  
  // martingale residuals
  std::vector<double> resmart(n);
  
  // linear predictors
  std::vector<double> linear_predictors(n);
  
  double zcrit = boost_qnorm(1.0 - alpha / 2.0);
  double xcrit = zcrit * zcrit;
  
  int nobs = n;
  int nevents = std::accumulate(eventn.begin(), eventn.end(), 0); 
  
  if (nevents == 0) {
    if (p > 0) {
      for (int i = 0; i < p; ++i) {
        par[i] = covariates[i];
        b[i] = NaN;
        seb[i] = 0;
        rseb[i] = 0;
        z[i] = NaN;
        expbeta[i] = NaN;
        lb[i] = NaN;
        ub[i] = NaN;
        prob[i] = NaN;
        clparm[i] = "Wald";
      }
      
      for (int j = 0; j < p; ++j) {
        for (int i = 0; i < p; ++i) {
          vb(i,j) = 0;
          rvb(i,j) = 0;
        }
      }
    }
    
    // baseline hazard
    if (est_basehaz) {
      dstratum[0] = 0;
      dtime[0] = 0;
      dnrisk[0] = n;
      dnevent[0] = 0;
      dncensor[0] = 0;
      dhaz[0] = 0;
      dvarhaz[0] = 0;
      if (p > 0) {
        for (int i = 0; i < p; ++i) {
          dgradhaz(0,i) = 0;
        }
      }
    }
    
    // martingale residuals
    if (est_resid) {
      for (int i = 0; i < n; ++i) {
        resmart[i] = 0;
      }
    }
    
    // linear predictors
    for (int i = 0; i < n; ++i) {
      linear_predictors[i] = offsetn[i];
    }
    
    loglik0 = NaN;
    loglik1 = NaN;
    regloglik0 = NaN;
    regloglik1 = NaN;
    scoretest = NaN;
    niter = 0;
    fail = true;
  } else {
    // sort by stratum
    std::vector<int> order0 = seqcpp(0, n-1);
    std::sort(order0.begin(), order0.end(), [&](int i, int j) {
      return stratumn[i] < stratumn[j];
    });
    
    std::vector<int> stratumnz = subset(stratumn, order0);
    std::vector<double> tstartnz = subset(tstartn, order0);
    std::vector<double> tstopnz = subset(tstopn, order0);
    std::vector<int> eventnz = subset(eventn, order0);
    
    // locate the first observation within each stratum
    std::vector<int> istratum(1,0);
    for (int i = 1; i < n; ++i) {
      if (stratumnz[i] != stratumnz[i-1]) {
        istratum.push_back(i);
      }
    }
    
    istratum.push_back(n);
    
    // ignore subjects not at risk for any event time
    std::vector<int> ignorenz(n);
    for (int i = 0; i < nstrata; ++i) {
      int start = istratum[i], end = istratum[i+1];
      int n0 = end - start;
      std::vector<double> tstart0 = subset(tstartnz, start, end);
      std::vector<double> tstop0 = subset(tstopnz, start, end);
      std::vector<int> event0 = subset(eventnz, start, end);
      
      // unique event times
      std::vector<double> etime;
      etime.reserve(n0);
      for (int j = 0; j < n0; ++j) {
        if (event0[j] == 1) etime.push_back(tstop0[j]);
      }
      etime = unique_sorted(etime);
      
      std::vector<int> index1 = findInterval3(tstart0, etime);
      std::vector<int> index2 = findInterval3(tstop0, etime);
      for (int j = istratum[i]; j < istratum[i+1]; ++j) {
        int j0 = j - istratum[i];
        if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
          ignorenz[j] = 1;
        } else {
          ignorenz[j] = 0;
        }
      }
    }
    
    std::vector<int> ignoren(n); // back to the original order
    for (int i = 0; i < n; ++i) {
      ignoren[order0[i]] = ignorenz[i];
    }
    
    int nused = n - std::accumulate(ignoren.begin(), ignoren.end(), 0);
    
    // sort by stopping time in descending order within each stratum
    std::vector<int> order1 = seqcpp(0, n-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      if (ignoren[i] != ignoren[j]) return ignoren[i] < ignoren[j];
      if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
      if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
      return eventn[i] < eventn[j];
    });
    
    std::vector<int> stratumna = subset(stratumn, order1);
    std::vector<double> tstartna = subset(tstartn, order1);
    std::vector<double> tstopna = subset(tstopn, order1);
    std::vector<int> eventna = subset(eventn, order1);
    std::vector<double> weightna = subset(weightn, order1);
    std::vector<double> offsetna = subset(offsetn, order1);
    std::vector<int> idna = subset(idn, order1);
    std::vector<int> ignorena = subset(ignoren, order1);
    FlatMatrix zna;
    if (p > 0) zna = subset_flatmatrix(zn, order1);
    
    // sort by starting time in descending order within each stratum
    std::vector<int> orderna = seqcpp(0, n-1);
    std::sort(orderna.begin(), orderna.end(), [&](int i, int j) {
      if (ignorena[i] != ignorena[j]) return ignorena[i] < ignorena[j];
      if (stratumna[i] != stratumna[j]) return stratumna[i] < stratumna[j];
      return tstartna[i] > tstartna[j];
    });
    
    coxparams param = {nused, stratumna, tstartna, tstopna, eventna,
                       weightna, offsetna, zna, orderna, method};
    
    std::vector<double> bint(p);
    ListCpp derint = f_der_2(p, bint, &param, firth);
    
    ListCpp out;
    
    if (p > 0) {
      std::vector<int> colfit = seqcpp(0, p - 1);
      if  (!init.empty() && static_cast<int>(init.size()) == p &&
           std::none_of(init.begin(), init.end(), [](double val){
             return std::isnan(val); })) {
        out = phregloop(p, init, &param, maxiter, eps, firth, colfit, p);
      } else {
        out = phregloop(p, bint, &param, maxiter, eps, firth, colfit, p);
      }
      
      niter = out.get<int>("iter");
      fail = out.get<bool>("fail");
      if (fail) {
        thread_utils::push_thread_warning(
          "phregloop failed to converge for the full model; "
          "continuing with current results.");
      }
      
      b = out.get<std::vector<double>>("coef");
      vb = out.get<FlatMatrix>("var");
      
      for (int j = 0; j < p; ++j) {
        seb[j] = std::sqrt(vb(j,j));
      }
      
      for (int i = 0; i < p; ++i) {
        par[i] = covariates[i];
      }
      
      // score statistic
      std::vector<double> scorebint;
      if (firth) scorebint = derint.get<std::vector<double>>("regscore");
      else scorebint = derint.get<std::vector<double>>("score");
      FlatMatrix infobint = derint.get<FlatMatrix>("imat");
      FlatMatrix vbint = invsympd(infobint, p);
      scoretest = quadsym(scorebint, vbint);
      
      
      // robust variance estimates
      if (robust) {
        FlatMatrix ressco = f_ressco_2(p, b, &param);
        
        int nr; // number of rows in the score residual matrix
        if (!has_id) {
          for (int j = 0; j < p; ++j) {
            double* rcol = ressco.data_ptr() + j * n;
            for (int i = 0; i < n; ++i) {
              rcol[i] *= weightna[i];
            }
          }
          nr = n;
        } else { // need to sum up score residuals by id
          std::vector<int> order = seqcpp(0, n-1);
          std::sort(order.begin(), order.end(), [&](int i, int j) {
            return idna[i] < idna[j];
          });
          
          std::vector<int> id1 = subset(idna, order);
          std::vector<int> idx(1,0);
          for (int i = 1; i < n; ++i) {
            if (id1[i] != id1[i-1]) {
              idx.push_back(i);
            }
          }
          
          int nids = static_cast<int>(idx.size());
          idx.push_back(n);
          
          FlatMatrix ressco1(nids,p);
          for (int j = 0; j < p; ++j) {
            const double* rcol = ressco.data_ptr() + j * n;
            double* rcol1 = ressco1.data_ptr() + j * nids;
            for (int i = 0; i < nids; ++i) {
              double sum = 0.0;
              for (int k = idx[i]; k < idx[i+1]; ++k) {
                int row = order[k];
                sum  += weightna[row] * rcol[row];
              }
              rcol1[i] = sum;
            }
          }
          
          ressco = std::move(ressco1);  // update the score residuals
          nr = nids;
        }
        
        FlatMatrix D = mat_mat_mult(ressco, vb); // DFBETA
        
        const double* Dptr = D.data_ptr();
        double* rvbptr = rvb.data_ptr();
        for (int j = 0; j < p; ++j) {
          const double* Dj = Dptr + j * nr; // pointer to D(:,j)
          for (int k = 0; k <= j; ++k) {
            const double* Dk = Dptr + k * nr; // pointer to D(:,k)
            double sum = 0.0;
            for (int i = 0; i < nr; ++i) {
              sum += Dj[i] * Dk[i];
            }
            rvbptr[k * p + j] = sum;
            if (j != k) rvbptr[j * p + k] = sum;
          }
        }
        
        for (int i = 0; i < p; ++i) {
          rseb[i] = std::sqrt(rvb(i,i));
        }
      }
      
      // profile likelihood confidence interval for regression coefficients
      if (plci) {
        double lmax = out.get<double>("loglik");
        double l0 = lmax - 0.5 * xcrit;
        
        for (int k = 0; k < p; ++k) {
          lb[k] = phregplloop(p, b, &param, maxiter, eps, firth, k, -1, l0);
          ub[k] = phregplloop(p, b, &param, maxiter, eps, firth, k, 1, l0);
          
          std::vector<int> colfit1(p-1);
          for (int i = 0, j = 0; i < p; ++i) {
            if (i == k) continue;
            colfit1[j++] = i;
          }
          
          std::vector<double> b0(p);
          ListCpp out0 = phregloop(p, b0, &param, maxiter, eps, firth, colfit1, p-1);
          double lmax0 = out0.get<double>("loglik");
          prob[k] = boost_pchisq(-2.0 * (lmax0 - lmax), 1, 0);
          clparm[k] = "PL";
        }
      } else {
        for (int k = 0; k < p; ++k) {
          if (!robust) {
            lb[k] = b[k] - zcrit * seb[k];
            ub[k] = b[k] + zcrit * seb[k];
            prob[k] = boost_pchisq(sq(b[k] / seb[k]), 1, 0);
          } else {
            lb[k] = b[k] - zcrit * rseb[k];
            ub[k] = b[k] + zcrit * rseb[k];
            prob[k] = boost_pchisq(sq(b[k] / rseb[k]), 1, 0);
          }
          clparm[k] = "Wald";
        }
      }
      
      if (firth) {
        loglik0 = derint.get<double>("loglik");
        loglik1 = out.get<double>("loglik");
        regloglik0 = derint.get<double>("regloglik");
        regloglik1 = out.get<double>("regloglik");
      } else {
        loglik0 = derint.get<double>("loglik");
        loglik1 = out.get<double>("loglik");
        regloglik0 = loglik0;
        regloglik1 = loglik1;
      }
    } else {
      fail = false;
      loglik0 = derint.get<double>("loglik");
      loglik1 = derint.get<double>("loglik");
      regloglik0 = loglik0;
      regloglik1 = loglik1;
      scoretest = 0.0;
      niter = 0;
    }
    
    // estimate baseline hazard
    if (est_basehaz) {
      // prepare the data for estimating baseline hazards at all time points
      
      // sort by stopping time in descending order within each stratum
      std::vector<int> order2 = seqcpp(0, n-1);
      std::sort(order2.begin(), order2.end(), [&](int i, int j) {
        if (stratumn[i] != stratumn[j]) return stratumn[i] > stratumn[j];
        return tstopn[i] > tstopn[j];
      });
      
      std::vector<int> stratumnb = subset(stratumn, order2);
      std::vector<double> tstartnb = subset(tstartn, order2);
      std::vector<double> tstopnb = subset(tstopn, order2);
      std::vector<int> eventnb = subset(eventn, order2);
      std::vector<double> weightnb = subset(weightn, order2);
      std::vector<double> offsetnb = subset(offsetn, order2);
      FlatMatrix znb;
      if (p > 0) znb = subset_flatmatrix(zn, order2);
      
      // sort by starting time in descending order within each stratum
      std::vector<int> ordernb = seqcpp(0, n-1);
      std::sort(ordernb.begin(), ordernb.end(), [&](int i, int j) {
        if (stratumnb[i] != stratumnb[j]) return stratumnb[i] > stratumnb[j];
        return tstartnb[i] > tstartnb[j];
      });
      
      coxparams paramb = {n, stratumnb, tstartnb, tstopnb, eventnb,
                          weightnb, offsetnb, znb, ordernb, method};
      
      ListCpp basehazn = f_basehaz(p, b, &paramb);
      
      dstratum = basehazn.get<std::vector<int>>("stratum");
      dtime = basehazn.get<std::vector<double>>("time");
      dnrisk = basehazn.get<std::vector<double>>("nrisk");
      dnevent = basehazn.get<std::vector<double>>("nevent");
      dncensor = basehazn.get<std::vector<double>>("ncensor");
      dhaz = basehazn.get<std::vector<double>>("haz");
      dvarhaz = basehazn.get<std::vector<double>>("varhaz");
      if (p > 0) dgradhaz = basehazn.get<FlatMatrix>("gradhaz");
    }
    
    // martingale residuals
    if (est_resid) {
      std::vector<double> resid = f_resmart(p, b, &param);
      for (int i = 0; i < n; ++i) {
        resmart[order1[i]] = resid[i];
      }
    }
    
    // linear predictors
    for (int i = 0; i < n; ++i) {
      linear_predictors[order1[i]] = offsetna[i];
    }
    
    
    if (p > 0) {
      for (int j = 0; j < p; ++j) {
        double beta = b[j];
        if (beta == 0.0) continue;
        const double* zna_col = zna.data_ptr() + j * n;
        for (int i = 0; i < n; ++i) {
          linear_predictors[order1[i]] += beta * zna_col[i];
        }
      }
    }
    
    // compute exp(beta)
    for (int i = 0; i < p; ++i) {
      expbeta[i] = std::exp(b[i]);
    }
    
    // compute z statistics
    if (robust) {
      for (int i = 0; i < p; ++i) {
        if (rseb[i] == 0) {
          z[i] = NaN;
        } else {
          z[i] = b[i] / rseb[i];
        }
      }
    } else {
      for (int i = 0; i < p; ++i) {
        if (seb[i] == 0) {
          z[i] = NaN;
        } else {
          z[i] = b[i] / seb[i];
        }
      }
    }
  }
  
  // prepare the output data sets
  DataFrameCpp sumstat;
  sumstat.push_back(nobs, "n");
  sumstat.push_back(nevents, "nevents");
  sumstat.push_back(loglik0, "loglik0");
  sumstat.push_back(loglik1, "loglik1");
  sumstat.push_back(scoretest, "scoretest");
  sumstat.push_back(niter, "niter");
  sumstat.push_back(meth, "ties");
  sumstat.push_back(p, "p");
  sumstat.push_back(robust, "robust");
  sumstat.push_back(firth, "firth");
  sumstat.push_back(fail, "fail");
  
  if (p > 0 && firth) {
    sumstat.push_back(regloglik0, "loglik0_unpenalized");
    sumstat.push_back(regloglik1, "loglik1_unpenalized");
  }
  
  ListCpp result;
  result.push_back(std::move(sumstat), "sumstat");
  
  if (p > 0) {
    std::vector<double> sebeta = robust ? rseb : seb;
    FlatMatrix vbeta = robust ? rvb : vb;
    DataFrameCpp parest;
    parest.push_back(std::move(par), "param");
    parest.push_back(std::move(b), "beta");
    parest.push_back(std::move(sebeta), "sebeta");
    parest.push_back(std::move(z), "z");
    parest.push_back(std::move(expbeta), "expbeta");
    parest.push_back(std::move(lb), "lower");
    parest.push_back(std::move(ub), "upper");
    parest.push_back(std::move(prob), "p");
    parest.push_back(std::move(clparm), "method");
    if (robust) parest.push_back(std::move(seb), "sebeta_naive");
    
    result.push_back(std::move(parest), "parest");
    result.push_back(std::move(vbeta), "vbeta");
    if (robust) result.push_back(std::move(vb), "vbeta_naive");
  }
  
  if (est_basehaz) {
    DataFrameCpp basehaz;
    basehaz.push_back(std::move(dtime), "time");
    basehaz.push_back(std::move(dnrisk), "nrisk");
    basehaz.push_back(std::move(dnevent), "nevent");
    basehaz.push_back(std::move(dncensor), "ncensor");
    basehaz.push_back(std::move(dhaz), "haz");
    basehaz.push_back(std::move(dvarhaz), "varhaz");
    
    if (p > 0) basehaz.push_back(std::move(dgradhaz), "gradhaz");
    
    if (has_stratum) {
      for (int i = 0; i < p_stratum; ++i) {
        std::string s = stratum[i];
        if (u_stratum.int_cols.count(s)) {
          auto v = u_stratum.get<int>(s);
          basehaz.push_back(subset(v, dstratum), s);
        } else if (u_stratum.numeric_cols.count(s)) {
          auto v = u_stratum.get<double>(s);
          basehaz.push_back(subset(v, dstratum), s);
        } else if (u_stratum.string_cols.count(s)) {
          auto v = u_stratum.get<std::string>(s);
          basehaz.push_back(subset(v, dstratum), s);
        } else {
          throw std::invalid_argument("unsupported type for stratum variable " + s);
        }
      }
    }
    
    result.push_back(std::move(basehaz), "basehaz");
  }
  
  if (est_resid) result.push_back(std::move(resmart), "residuals");
  result.push_back(std::move(linear_predictors), "linear_predictors");
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List phregRcpp(const Rcpp::DataFrame& data,
                     const std::vector<std::string>& stratum,
                     const std::string& time,
                     const std::string& time2,
                     const std::string& event,
                     const std::vector<std::string>& covariates,
                     const std::string& weight,
                     const std::string& offset,
                     const std::string& id,
                     const std::string& ties,
                     const std::vector<double>& init,
                     const bool robust,
                     const bool est_basehaz,
                     const bool est_resid,
                     const bool firth,
                     const bool plci,
                     const double alpha,
                     const int maxiter,
                     const double eps) {
  
  auto dfcpp = convertRDataFrameToCpp(data);
  
  auto cpp_result = phregcpp(
    dfcpp, stratum, time, time2, event, covariates, weight, offset, id, 
    ties, init, robust, est_basehaz, est_resid, firth, plci, alpha, 
    maxiter, eps
  );
  
  thread_utils::drain_thread_warnings_to_R();
  return Rcpp::wrap(cpp_result);
}


// survival function estimation based on Cox model
DataFrameCpp survfit_phregcpp(const int p,
                              const std::vector<double>& beta,
                              const FlatMatrix& vbeta,
                              const DataFrameCpp& basehaz,
                              const DataFrameCpp& newdata,
                              const std::vector<std::string>& covariates,
                              const std::vector<std::string>& stratum,
                              const std::string& offset,
                              const std::string& id,
                              const std::string& tstart,
                              const std::string& tstop,
                              const bool sefit,
                              const std::string& conftype,
                              const double conflev) {
  
  int n0 = basehaz.nrows();
  int n = newdata.nrows();
  int nvar = static_cast<int>(covariates.size());
  
  std::string ct = conftype;
  std::for_each(ct.begin(), ct.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (!(ct=="none" || ct=="plain" || ct=="log" || ct=="log-log" || 
      ct=="logit" || ct=="arcsin")) {
    throw std::invalid_argument(
        "conftype must be none, plain, log, log-log, logit, or arcsin");
  }
  
  if (conflev <= 0.0 || conflev >= 1.0) {
    throw std::invalid_argument("conflev must lie between 0 and 1");
  }
  
  double zcrit = boost_qnorm((1.0 + conflev) / 2.0);
  
  std::vector<int> stratumn0(n0);
  DataFrameCpp u_stratum0;
  std::vector<int> nlevels;
  ListPtr lookups_ptr;
  int p_stratum = static_cast<int>(stratum.size());
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(basehaz, stratum);
    stratumn0 = out.get<std::vector<int>>("index");
    u_stratum0 = out.get<DataFrameCpp>("lookup");
    nlevels = out.get<std::vector<int>>("nlevels");
    lookups_ptr = out.get<ListPtr>("lookups_per_variable");
  }
  
  bool nullmodel = (p == 0 || (nvar == 1 && covariates[0] == ""));
  
  if (!nullmodel && nvar != p) {
    throw std::invalid_argument("incorrect number of covariates for the Cox model");
  }
  
  FlatMatrix zn(n,p);
  for (int j = 0; j < p; ++j) {
    const std::string& zj = covariates[j];
    if (!newdata.containElementNamed(zj))
      throw std::invalid_argument(
          "newdata must contain the variables in covariates");
    double* zn_col = zn.data_ptr() + j * n;
    if (newdata.bool_cols.count(zj)) {
      const std::vector<unsigned char>& vb = newdata.get<unsigned char>(zj);
      for (int i = 0; i < n; ++i) zn_col[i] = vb[i] ? 1.0 : 0.0;
    } else if (newdata.int_cols.count(zj)) {
      const std::vector<int>& vi = newdata.get<int>(zj);
      for (int i = 0; i < n; ++i) zn_col[i] = static_cast<double>(vi[i]);
    } else if (newdata.numeric_cols.count(zj)) {
      const std::vector<double>& vd = newdata.get<double>(zj);
      std::memcpy(zn_col, vd.data(), n * sizeof(double));
    } else {
      throw std::invalid_argument("covariates must be bool, integer or numeric");
    }
  }
  
  bool has_stratum = false;
  std::vector<int> stratumn(n);
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    has_stratum = true;
    const ListCpp& lookups = *lookups_ptr;
    
    // match stratum in newdata to stratum in basehaz
    int orep = u_stratum0.nrows();
    for (int i = 0; i < p_stratum; ++i) {
      orep /= nlevels[i];
      std::string s = stratum[i];
      std::vector<int> idx;
      
      if (newdata.bool_cols.count(s) || newdata.int_cols.count(s)) {
        std::vector<int> v;
        std::vector<int> w;
        if (newdata.bool_cols.count(s)) {
          auto vb = newdata.get<unsigned char>(s);
          auto wb = lookups.get<std::vector<unsigned char>>(s);
          v.resize(n);
          w.resize(n);
          for (int j = 0; j < n; ++j) {
            v[j] = vb[j] ? 1 : 0;
            w[j] = wb[j] ? 1 : 0;
          }
        } else {
          v = newdata.get<int>(s);
          w = lookups.get<std::vector<int>>(s);
        }
        idx = matchcpp(v, w);
      } else if (newdata.numeric_cols.count(s)) {
        auto v = newdata.get<double>(s);
        auto w = lookups.get<std::vector<double>>(s);
        idx = matchcpp(v, w);
      } else if (newdata.string_cols.count(s)) {
        auto v = newdata.get<std::string>(s);
        auto w = lookups.get<std::vector<std::string>>(s);
        idx = matchcpp(v, w);
      } else {
        throw std::invalid_argument("Unsupported type for stratum variable: " + s);
      }
      
      for (int person = 0; person < n; ++person) {
        stratumn[person] += idx[person] * orep;
      }
    }
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && newdata.containElementNamed(offset)) {
    if (newdata.int_cols.count(offset)) {
      const std::vector<int>& vi = newdata.get<int>(offset);
      for (int i = 0; i < n; ++i) offsetn[i] = static_cast<double>(vi[i]);
    } else if (newdata.numeric_cols.count(offset)) {
      offsetn = newdata.get<double>(offset);
    } else {
      throw std::invalid_argument("offset variable must be integer or numeric");
    }
  }
  
  std::vector<double> time0 = basehaz.get<double>("time");
  std::vector<double> nrisk0 = basehaz.get<double>("nrisk");
  std::vector<double> nevent0 = basehaz.get<double>("nevent");
  std::vector<double> ncensor0 = basehaz.get<double>("ncensor");
  std::vector<double> haz0 = basehaz.get<double>("haz");
  std::vector<double> vhaz0 = basehaz.get<double>("varhaz");
  FlatMatrix ghaz0(n0,p);
  for (int j = 0; j < p; ++j) {
    std::string col_name = "gradhaz";
    if (p>1) col_name += "." + std::to_string(j+1);
    std::vector<double> u = basehaz.get<double>(col_name);
    for (int i = 0; i < n0; ++i) {
      ghaz0(i,j) = u[i];
    }
  }
  
  // create the numeric id variable
  bool has_id = !id.empty() && newdata.containElementNamed(id);
  std::vector<int> idn(n);
  std::vector<int> idwi;
  std::vector<double> idwn;
  std::vector<std::string> idwc;
  if (!has_id) {
    idn = seqcpp(0, n-1);
  } else { // input data has the counting process style of input
    if (newdata.int_cols.count(id)) {
      auto v = newdata.get<int>(id);
      idwi = unique_sorted(v);
      idn = matchcpp(v, idwi);
    } else if (newdata.numeric_cols.count(id)) {
      auto v = newdata.get<double>(id);
      idwn = unique_sorted(v);
      idn = matchcpp(v, idwn);
    } else if (newdata.string_cols.count(id)) {
      auto v = newdata.get<std::string>(id);
      idwc = unique_sorted(v);
      idn = matchcpp(v, idwc);
    } else throw std::invalid_argument(
        "incorrect type for the id variable in newdata");
  }
  
  // unify right-censoring data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_id) { // right-censored data
    double maxt0 = *std::max_element(time0.begin(), time0.end()) + 1.0;
    std::fill(tstopn.begin(), tstopn.end(), maxt0);
  } else {
    if (!newdata.containElementNamed(tstart))
      throw std::invalid_argument("newdata must contain the tstart variable");
    if (newdata.int_cols.count(tstart)) {
      const std::vector<int>& vi = newdata.get<int>(tstart);
      for (int i = 0; i < n; ++i) tstartn[i] = static_cast<double>(vi[i]);
    } else if (newdata.numeric_cols.count(tstart)) {
      tstartn = newdata.get<double>(tstart);
    } else {
      throw std::invalid_argument("tstart variable must be integer or numeric");
    }
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(tstartn[i]) && tstartn[i] < 0.0)
        throw std::invalid_argument("tstart must be nonnegative");
    }
    
    if (!newdata.containElementNamed(tstop))
      throw std::invalid_argument("newdata must contain the tstop variable");
    if (newdata.int_cols.count(tstop)) {
      const std::vector<int>& vi = newdata.get<int>(tstop);
      for (int i = 0; i < n; ++i) tstopn[i] = static_cast<double>(vi[i]);
    } else if (newdata.numeric_cols.count(tstop)) {
      tstopn = newdata.get<double>(tstop);
    } else {
      throw std::invalid_argument("tstop variable must be integer or numeric");
    }
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(tstartn[i]) && !std::isnan(tstopn[i]) && 
          tstopn[i] <= tstartn[i])
        throw std::invalid_argument("tstop must be greater than tstart");
    }
  }
  
  // order data by id and tstop, assuming consecutive intervals
  std::vector<int> order = seqcpp(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (idn[i] != idn[j]) return idn[i] < idn[j];
    return tstopn[i] < tstopn[j];
  });
  
  subset_in_place(idn, order);
  subset_in_place(stratumn, order);
  subset_in_place(tstartn, order);
  subset_in_place(tstopn, order);
  subset_in_place(offsetn, order);
  if (p > 0) subset_in_place_flatmatrix(zn, order);
  
  // exclude observations with missing values
  std::vector<unsigned char> sub(n,1);
  for (int i = 0; i < n; ++i) {
    if (stratumn[i] == INT_MIN || idn[i] == INT_MIN ||
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) ||
        std::isnan(offsetn[i])) {
      sub[i] = 0; continue;
    }
    for (int j = 0; j < p; ++j) {
      if (std::isnan(zn(i,j))) { sub[i] = 0; break; }
    }
  }
  
  std::vector<int> keep = which(sub);
  if (keep.empty()) 
    throw std::invalid_argument("no observations without missing values");
  subset_in_place(stratumn, keep);
  subset_in_place(offsetn, keep);
  subset_in_place(idn, keep);
  subset_in_place(tstartn, keep);
  subset_in_place(tstopn, keep);
  if (p > 0) subset_in_place_flatmatrix(zn, keep);
  n = static_cast<int>(keep.size());
  
  // risk score
  std::vector<double> eta = offsetn;
  for (int j = 0; j < p; ++j) {
    double b = beta[j];
    if (b == 0.0) continue;
    const double* zn_col = zn.data_ptr() + j * n;
    for (int i = 0; i < n; ++i) {
      eta[i] += b * zn_col[i];
    }
  }
  std::vector<double> risk(n);
  for (int i = 0; i < n; ++i) {
    risk[i] = std::exp(eta[i]);
  }
  
  // count number of observations for each id
  std::vector<int> idx(1,0);
  for (int i = 1; i < n; ++i) {
    if (idn[i] != idn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nids = static_cast<int>(idx.size());
  idx.push_back(n);
  
  int N = nids * n0; // upper bound on the number of rows in the output
  std::vector<double> time(N), nrisk(N), nevent(N), ncensor(N);
  std::vector<double> cumhaz(N), vcumhaz(N), secumhaz(N);
  std::vector<int> strata(N), ids(N);
  FlatMatrix z(N,p);
  
  // process by id
  int l = 0;
  for (int h = 0; h < nids; ++h) {
    int start = idx[h], end = idx[h+1];
    int n1 = end - start;
    std::vector<int> id1 = subset(idn, start, end);
    std::vector<int> stratum1 = subset(stratumn, start, end);
    std::vector<double> tstart1 = subset(tstartn, start, end);
    std::vector<double> tstop1 = subset(tstopn, start, end);
    std::vector<double> risk1 = subset(risk, start, end);
    FlatMatrix z1 = subset_flatmatrix(zn, start, end);
    
    std::vector<double> tstop2(n1 + 1);
    tstop2[0] = tstart1[0];
    std::memcpy(tstop2.data() + 1, tstop1.data(), n1 * sizeof(double));
    
    // match the stratum in basehaz
    std::vector<int> idx1;
    for (int i = 0; i < n0; ++i) {
      if (stratumn0[i] == stratum1[0]) idx1.push_back(i);
    }
    std::vector<double> time01 = subset(time0, idx1);
    
    // left-open and right-closed intervals containing the event time
    std::vector<int> idx2 = findInterval3(time01, tstop2, 0, 0, 1);
    std::vector<int> sub;
    for (int i = 0; i < static_cast<int>(idx2.size()); ++i) {
      if (idx2[i] >= 1 && idx2[i] <= n1) sub.push_back(i);
    }
    int m1 = sub.size();
    
    if (m1 != 0) {
      std::vector<int> idx3 = subset(idx1, sub);
      std::vector<double> time1 = subset(time0, idx3);
      std::vector<double> nrisk1 = subset(nrisk0, idx3);
      std::vector<double> nevent1 = subset(nevent0, idx3);
      std::vector<double> ncensor1 = subset(ncensor0, idx3);
      std::vector<double> haz1 = subset(haz0, idx3);
      
      std::vector<int> idx4 = subset(idx2, sub);
      for (int i = 0; i < m1; ++i) idx4[i] -= 1; // change to 0-1 indexing
      
      // cumulative hazards
      for (int i = 0; i < m1; ++i) {
        int r = l + i;
        time[r] = time1[i];
        nrisk[r] = nrisk1[i];
        nevent[r] = nevent1[i];
        ncensor[r] = ncensor1[i];
        
        int k = idx4[i];
        ids[r] = id1[k];
        strata[r] = stratum1[k];
        for (int j = 0; j < p; ++j) {
          z(r,j) = z1(k,j);
        }
        
        if (i==0) {
          cumhaz[r] = haz1[i] * risk1[k];
        } else {
          cumhaz[r] = cumhaz[r-1] + haz1[i] * risk1[k];
        }
      }
      
      if (sefit) {
        std::vector<double> vhaz1 = subset(vhaz0, idx3);
        FlatMatrix ghaz1(m1,p);
        for (int j = 0; j < p; ++j) {
          for (int i = 0; i < m1; ++i) {
            ghaz1(i,j) = ghaz0(idx3[i],j);
          }
        }
        
        FlatMatrix a(m1,p);
        for (int j = 0; j < p; ++j) {
          for (int i = 0; i < m1; ++i) {
            int k = idx4[i];
            if (i == 0) {
              a(i,j) = (haz1[i] * z1(k,j) - ghaz1(i,j)) * risk1[k];
            } else {
              a(i,j) = a(i-1,j) + (haz1[i] * z1(k,j) - ghaz1(i,j)) * risk1[k];
            }
          }
        }
        
        // calculate the first component of variance
        for (int i = 0; i < m1; ++i) {
          int r = l + i;
          int k = idx4[i];
          if (i == 0) {
            vcumhaz[r] = vhaz1[i] * risk1[k] * risk1[k];
          } else {
            vcumhaz[r] = vcumhaz[r-1] + vhaz1[i] * risk1[k] * risk1[k];
          }
        }
        
        // add the second component of variance
        for (int k = 0; k < p; ++k) {
          for (int j = 0; j < p; ++j) {
            for (int i = 0; i < m1; ++i) {
              int r = l + i;
              vcumhaz[r] += a(i,j) * vbeta(j,k) * a(i,k);
            }
          }
        }
        
        for (int i = 0; i < m1; ++i) {
          int r = l + i;
          secumhaz[r] = std::sqrt(vcumhaz[r]);
        }
      }
      
      l += m1;
    }
  }
  
  subset_in_place(time, 0, l);
  subset_in_place(nrisk, 0, l);
  subset_in_place(nevent, 0, l);
  subset_in_place(ncensor, 0, l);
  subset_in_place(cumhaz, 0, l);
  subset_in_place(secumhaz, 0, l);
  subset_in_place(strata, 0, l);
  subset_in_place(ids, 0, l);
  subset_in_place_flatmatrix(z, 0, l);
  
  std::vector<double> surv(l);
  for (int i = 0; i < l; ++i) {
    surv[i] = std::exp(-cumhaz[i]);
  }
  
  DataFrameCpp result;
  result.push_back(std::move(time), "time");
  result.push_back(std::move(nrisk), "nrisk");
  result.push_back(std::move(nevent), "nevent");
  result.push_back(std::move(ncensor), "ncensor");
  result.push_back(std::move(cumhaz), "cumhaz");
  result.push_back(surv, "surv");
  
  if (sefit) {
    std::vector<double> sesurv(l);
    for (int i = 0; i < l; ++i) {
      sesurv[i] = surv[i] * secumhaz[i];
    }
    
    std::vector<double> lower(l), upper(l);
    for (int i = 0; i < l; ++i) {
      std::vector<double> ci = fsurvci(surv[i], sesurv[i], ct, zcrit);
      lower[i] = ci[0];
      upper[i] = ci[1];
    }
    
    result.push_back(std::move(sesurv), "sesurv");
    result.push_back(std::move(lower), "lower");
    result.push_back(std::move(upper), "upper");
    result.push_back(conflev, "conflev");
    result.push_back(ct, "conftype");
  }
  
  for (int j = 0; j < p; ++j) {
    std::string zj = covariates[j];
    std::vector<double> u(l);
    const double* zcol = z.data_ptr() + j * l;
    for (int i = 0; i < l; ++i) {
      u[i] = zcol[i];
    }
    result.push_back(std::move(u), zj);
  }
  
  if (has_stratum) {
    for (int i = 0; i < p_stratum; ++i) {
      std::string s = stratum[i];
      if (u_stratum0.int_cols.count(s)) {
        auto v = u_stratum0.get<int>(s);
        result.push_back(subset(v, strata), s);
      } else if (u_stratum0.numeric_cols.count(s)) {
        auto v = u_stratum0.get<double>(s);
        result.push_back(subset(v, strata), s);
      } else if (u_stratum0.string_cols.count(s)) {
        auto v = u_stratum0.get<std::string>(s);
        result.push_back(subset(v, strata), s);
      } else {
        throw std::invalid_argument("Unsupported type for stratum variable: " + s);
      }
    }
  }
  
  if (has_id) {
    if (newdata.int_cols.count(id)) {
      result.push_back(subset(idwi, ids), id);
    } else if (newdata.numeric_cols.count(id)) {
      result.push_back(subset(idwn, ids), id);
    } else if (newdata.string_cols.count(id)) {
      result.push_back(subset(idwc, ids), id);
    } else {
      throw std::invalid_argument("incorrect type for the id variable in newdata");
    }
  }
  
  return result;
}


// [[Rcpp::export]]
Rcpp::DataFrame survfit_phregRcpp(const int p,
                                  const std::vector<double>& beta,
                                  const Rcpp::NumericMatrix& vbeta,
                                  const Rcpp::DataFrame& basehaz,
                                  const Rcpp::DataFrame& newdata,
                                  const std::vector<std::string>& covariates,
                                  const std::vector<std::string>& stratum,
                                  const std::string& offset,
                                  const std::string& id,
                                  const std::string& tstart,
                                  const std::string& tstop,
                                  const bool sefit,
                                  const std::string& conftype,
                                  const double conflev) {
  
  auto vbetacpp = flatmatrix_from_Rmatrix(vbeta);
  auto basehcpp = convertRDataFrameToCpp(basehaz);
  auto newdfcpp = convertRDataFrameToCpp(newdata);
  
  auto cpp_result = survfit_phregcpp(
    p, beta, vbetacpp, basehcpp, newdfcpp, covariates, stratum, 
    offset, id, tstart, tstop, sefit, conftype, conflev
  );
  
  return Rcpp::wrap(cpp_result);
}


// schoenfeld residuals
ListCpp f_ressch(int p, const std::vector<double>& par, void *ex) {
  coxparams *param = (coxparams *) ex;
  const int n = param->tstop.size();
  const int nused = param->nused;
  const int method = param->method;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& event = param->event; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const std::vector<int>& order1 = param->order1;
  const FlatMatrix& z = param->z;
  
  // Precompute eta and risk = exp(eta)
  std::vector<double> eta(nused);
  for (int person = 0; person < nused; ++person) {
    eta[person] = offset[person];
  }
  if (p > 0) {
    for (int i = 0; i < p; ++i) {
      double beta = par[i];
      if (beta == 0.0) continue;
      const double* zcol = z.data_ptr() + i * n;
      for (int person = 0; person < nused; ++person) {
        eta[person] += beta * zcol[person];
      }
    }
  }
  std::vector<double> risk(nused);
  for (int person = 0; person < nused; ++person) {
    risk[person] = std::exp(eta[person]);
  }
  
  int nevent = 0;
  for (int i = 0; i < nused; ++i) if (event[i] != 0) ++nevent;
  
  FlatMatrix resid(nevent, p);     // residual matrix
  std::vector<int> index(nevent);  // index of residuals
  
  std::vector<double> xbar(p);     // weighted mean covariate at this time
  std::vector<double> a(p);        // s1(beta,k,t)
  std::vector<double> a2(p);       // sum of w*exp(zbeta)*z for the deaths
  double denom = 0.0;              // s0(beta,k,t)
  double denom2 = 0.0;             // sum of weighted risks for the deaths
  double ndead = 0.0;              // number of deaths at this time point
  
  int istrata = strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  int j = nevent;  // index the events in descending order
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    if (strata[person] != istrata) { // hit a new stratum
      istrata = strata[person]; // reset temporary variables
      i1 = person;
      denom = 0.0;
      std::fill(a.begin(), a.end(), 0.0);
    }
    
    const double dtime = tstop[person];
    
    // process all persons tied at this dtime
    for (; person < nused && tstop[person] == dtime && 
         strata[person] == istrata; ++person) {
      
      const double r = weight[person] * risk[person];
      
      if (event[person] == 0) {
        denom += r;
        for (int i = 0; i < p; ++i) {
          a[i] += r * z(person,i);
        }
      } else {
        --j;
        for (int i = 0; i < p; ++i) {
          resid(j,i) = z(person,i);  
        }
        index[j] = person;
        
        ++ndead;
        denom2 += r;
        for (int i = 0; i < p; ++i) {
          a2[i] += r * z(person,i);
        }
      }
    }
    
    // remove subjects no longer at risk
    for (; i1 < nused; ++i1) {
      const int p1 = order1[i1];
      if (tstart[p1] < dtime || strata[p1] != istrata) break;
      const double r = weight[p1] * risk[p1];
      denom -= r;
      for (int i = 0; i < p; ++i) {
        a[i] -= r * z(p1,i);
      }
    }
    
    // add to the main terms
    if (ndead > 0) {
      if (method == 0 || ndead == 1) {
        denom += denom2;
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          xbar[i] = a[i] / denom;
        }
      } else {
        std::fill(xbar.begin(), xbar.end(), 0.0);
        for (int k = 0; k < ndead; ++k) {
          denom += denom2 / ndead;
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            xbar[i] += a[i] / denom;
          }
        }
        for (int i = 0; i < p; ++i) {
          xbar[i] /= ndead;  
        }
      }
      
      for (int i = 0; i < p; ++i) {
        double xc = xbar[i];
        double* rptr = resid.data_ptr() + i * nevent;
        for (int k = 0; k < ndead; ++k) {
          rptr[j + k] -= xc;
        }
      }
      
      // reset for the next death time
      ndead = denom2 = 0.0;
      std::fill(a2.begin(), a2.end(), 0.0);
    }
  }
  
  ListCpp result;
  result.push_back(std::move(resid), "resid");
  result.push_back(std::move(index), "index");
  return result;
}


// residuals for phreg
ListCpp residuals_phregcpp(const int p,
                           const std::vector<double>& beta,
                           const FlatMatrix& vbeta,
                           const std::vector<double>& resmart,
                           const DataFrameCpp& data,
                           const std::vector<std::string>& stratum,
                           const std::string& time,
                           const std::string& time2,
                           const std::string& event,
                           const std::vector<std::string>& covariates,
                           const std::string& weight,
                           const std::string& offset,
                           const std::string& id,
                           const std::string& ties,
                           const std::string& type,
                           const bool collapse,
                           const bool weighted) {
  
  int n = data.nrows();
  
  // --- handle strata (bygroup) ---
  bool has_stratum = false;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  int p_stratum = static_cast<int>(stratum.size());
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
  }
  
  // --- time / time2 existence and checks ---
  if (!data.containElementNamed(time)) 
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(timen[i]) && timen[i] < 0.0)
      throw std::invalid_argument("time must be nonnegative");
  }
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    if (data.int_cols.count(time2)) {
      const std::vector<int>& vi = data.get<int>(time2);
      for (int i = 0; i < n; ++i) time2n[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(time2)) {
      time2n = data.get<double>(time2);
    } else {
      throw std::invalid_argument("time2 variable must be integer or numeric");
    }    
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(timen[i]) && !std::isnan(time2n[i]) && time2n[i] <= timen[i])
        throw std::invalid_argument("time2 must be greater than time");
    }
  }
  
  // --- event variable ---
  bool has_event = !event.empty() && data.containElementNamed(event);
  if (!has_time2 && !has_event) {
    throw std::invalid_argument(
        "data must contain the event variable for right censored data");
  }
  std::vector<int> eventn(n);
  if (has_event) {
    if (data.bool_cols.count(event)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
      for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
    } else if (data.int_cols.count(event)) {
      eventn = data.get<int>(event);
    } else if (data.numeric_cols.count(event)) {
      const std::vector<double>& vd = data.get<double>(event);
      for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
    } else {
      throw std::invalid_argument("event variable must be bool, integer or numeric");
    }
    for (double val : eventn) if (val != 0 && val != 1)
      throw std::invalid_argument("event must be 1 or 0 for each observation");
  }
  
  // --- build design matrix zn (n x p) column-major FlatMatrix ---
  FlatMatrix zn(n, p);
  if (p > 0) {
    for (int j = 0; j < p; ++j) {
      const std::string& zj = covariates[j];
      if (!data.containElementNamed(zj)) 
        throw std::invalid_argument("data must contain the variables in covariates");
      double* zn_col = zn.data_ptr() + j * n;
      if (data.bool_cols.count(zj)) {
        const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
        for (int i = 0; i < n; ++i) zn_col[i] = vb[i] ? 1.0 : 0.0;
      } else if (data.int_cols.count(zj)) {
        const std::vector<int>& vi = data.get<int>(zj);
        for (int i = 0; i < n; ++i) zn_col[i] = static_cast<double>(vi[i]);
      } else if (data.numeric_cols.count(zj)) {
        const std::vector<double>& vd = data.get<double>(zj);
        std::memcpy(zn_col, vd.data(), n * sizeof(double));
      } else {
        throw std::invalid_argument("covariates must be bool, integer or numeric");
      }
    }
  }
  
  // --- weight and offset ---
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const std::vector<int>& vi = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else {
      throw std::invalid_argument("weight variable must be integer or numeric");
    }
    for (double w : weightn) if (std::isnan(w) || w <= 0.0) 
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    if (data.int_cols.count(offset)) {
      const std::vector<int>& vi = data.get<int>(offset);
      for (int i = 0; i < n; ++i) offsetn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(offset)) {
      offsetn = data.get<double>(offset);
    } else {
      throw std::invalid_argument("offset variable must be integer or numeric");
    }
  }
  
  // --- id mapping ---
  bool has_id = !id.empty() && data.containElementNamed(id);
  std::vector<int> idn(n);
  if (!has_id) {
    std::iota(idn.begin(), idn.end(), 0);
  } else {
    if (data.int_cols.count(id)) {
      auto v = data.get<int>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.numeric_cols.count(id)) {
      auto v = data.get<double>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else if (data.string_cols.count(id)) {
      auto v = data.get<std::string>(id);
      auto w = unique_sorted(v);
      idn = matchcpp(v, w);
    } else throw std::invalid_argument(
        "incorrect type for the id variable in the input data");
  }
  
  if (type != "martingale" && type != "deviance" &&
      type != "score" && type != "dfbeta" && type != "dfbetas" &&
      type != "schoenfeld" && type != "scaledsch") {
    throw std::invalid_argument("unknown type of residuals");
  }
  
  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  int method = meth == "efron" ? 1 : 0;
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  // exclude observations with missing values
  std::vector<unsigned char> sub(n, 1);
  for (int i = 0; i < n; ++i) {
    if (stratumn[i] == INT_MIN || idn[i] == INT_MIN ||
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) ||
        eventn[i] == INT_MIN || std::isnan(weightn[i]) ||
        std::isnan(offsetn[i])) {
      sub[i] = 0; continue;
    }
    for (int j = 0; j < p; ++j) {
      if (std::isnan(zn(i,j))) { sub[i] = 0; break; }
    }
  }
  
  std::vector<int> keep = which(sub);
  if (keep.empty()) throw std::invalid_argument(
      "no observations left after removing missing values");
  
  subset_in_place(stratumn, keep);
  subset_in_place(tstartn, keep);
  subset_in_place(tstopn, keep);
  subset_in_place(eventn, keep);
  subset_in_place(weightn, keep);
  subset_in_place(offsetn, keep);
  subset_in_place(idn, keep);
  if (p > 0) subset_in_place_flatmatrix(zn, keep);
  n = static_cast<int>(keep.size());
  
  // sort by stratum
  std::vector<int> ordern = seqcpp(0, n-1);
  std::sort(ordern.begin(), ordern.end(), [&](int i, int j) {
    return stratumn[i] < stratumn[j];
  });
  
  std::vector<int> stratumnz = subset(stratumn, ordern);
  std::vector<double> tstartnz = subset(tstartn, ordern);
  std::vector<double> tstopnz = subset(tstopn, ordern);
  std::vector<int> eventnz = subset(eventn, ordern);
  
  // locate the first observation within each stratum
  std::vector<int> istratum(1,0);
  for (int i = 1; i < n; ++i) {
    if (stratumnz[i] != stratumnz[i-1]) {
      istratum.push_back(i);
    }
  }
  int nstrata = static_cast<int>(istratum.size());
  istratum.push_back(n);
  
  // ignore subjects not at risk for any event time
  std::vector<int> ignore1z(n, 0);
  for (int i = 0; i < nstrata; ++i) {
    int start = istratum[i], end = istratum[i+1];
    int m = end - start;
    std::vector<double> tstart0 = subset(tstartnz, start, end);
    std::vector<double> tstop0 = subset(tstopnz, start, end);
    std::vector<int> event0 = subset(eventnz, start, end);
    std::vector<double> etime;
    etime.reserve(m);
    for (int j = 0; j < m; ++j) {
      if (event0[j] == 1) etime.push_back(tstop0[j]); 
    }
    etime = unique_sorted(etime);
    std::vector<int> index1 = findInterval3(tstart0, etime);
    std::vector<int> index2 = findInterval3(tstop0, etime);
    for (int j = istratum[i]; j < istratum[i+1]; ++j) {
      int j0 = j - istratum[i];
      if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
        ignore1z[j] = 1;
      } else {
        ignore1z[j] = 0;
      }
    }
  }
  
  std::vector<int> ignore(n, 0);
  for (int i = 0; i < n; ++i) ignore[ordern[i]] = ignore1z[i];
  
  int nused = n - std::accumulate(ignore.begin(), ignore.end(), 0);
  
  std::vector<int> order = seqcpp(0, n-1);
  std::vector<int> idx(1,0);
  int nids = n;
  if (has_id) {
    std::sort(order.begin(), order.end(), [&](int i, int j){ 
      return idn[i] < idn[j]; 
    });
    std::vector<int> id1 = subset(idn, order);
    for (int i = 1; i < n; ++i) if (id1[i] != id1[i-1]) idx.push_back(i);
    nids = static_cast<int>(idx.size());
    idx.push_back(n);
  }
  
  ListCpp result;
  
  if (type == "martingale") {
    std::vector<double> rr = resmart;
    if (weighted) {
      for (int i = 0; i < n; ++i) rr[i] *= weightn[i];
    }
    if (collapse && has_id) { // collapse over id
      std::vector<double> rr1(nids, 0.0);
      for (int i = 0; i < nids; ++i) {
        for (int j = idx[i]; j < idx[i+1]; ++j) {
          rr1[i] += rr[order[j]];
        }
      }
      rr = std::move(rr1);
    }
    result.push_back(std::move(rr), "resid");
    return result;
  } 
  
  if (type == "deviance") {
    std::vector<double> rr = resmart;
    std::vector<int> status = eventn;
    int m = n;
    if (weighted) {
      for (int i = 0; i < n; ++i) rr[i] *= weightn[i];
    }
    if (collapse && has_id) {
      std::vector<double> rr1(nids, 0.0);
      std::vector<int> status1(nids, 0);
      for (int i = 0; i < nids; ++i) {
        for (int j = idx[i]; j < idx[i+1]; ++j) {
          int k = order[j];
          rr1[i] += rr[k];
          status1[i] += eventn[k];
        }
      }
      rr = std::move(rr1);
      status = std::move(status1);
      m = nids;
    }
    for (int i = 0; i < m; ++i) {
      double temp = status[i] == 0 ? 0 : status[i] * std::log(status[i] - rr[i]);
      rr[i] = ((rr[i] > 0) - (rr[i] < 0)) * std::sqrt(-2.0 * (rr[i] + temp));
    }
    result.push_back(std::move(rr), "resid");
    return result;
  } 
  
  if (p == 0) throw std::invalid_argument(
    "covariates must be present for score and schoenfeld residuals");
  
  if (type == "score" || type == "dfbeta" || type == "dfbetas") {
    // sort by stopping time in descending order within each stratum
    std::vector<int> order0 = seqcpp(0, n-1);
    std::sort(order0.begin(), order0.end(), [&](int i, int j){
      if (ignore[i] != ignore[j]) return ignore[i] < ignore[j];
      if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
      if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
      return eventn[i] < eventn[j];
    });
    
    std::vector<int> stratum1 = subset(stratumn, order0);
    std::vector<double> tstart1 = subset(tstartn, order0);
    std::vector<double> tstop1 = subset(tstopn, order0);
    std::vector<int> event1 = subset(eventn, order0);
    std::vector<double> weight1 = subset(weightn, order0);
    std::vector<double> offset1 = subset(offsetn, order0);
    std::vector<int> ignore1 = subset(ignore, order0);
    FlatMatrix z1 = subset_flatmatrix(zn, order0);
    
    // sort by starting time in descending order within each stratum
    std::vector<int> order1 = seqcpp(0, n-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j){
      if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
      if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
      return tstart1[i] > tstart1[j];
    });
    
    coxparams param = {nused, stratum1, tstart1, tstop1, event1, 
                       weight1, offset1, z1, order1, method};
    FlatMatrix ressco = f_ressco_2(p, beta, &param);
    
    // re-order ressco back to original order
    FlatMatrix score(n, p);
    for (int j = 0; j < p; ++j) {
      const double* rcol = ressco.data_ptr() + j * n;
      double* scol = score.data_ptr() + j * n;
      for (int i = 0; i < n; ++i) {
        scol[order0[i]] = rcol[i]; // original order
      }
    }
    
    // compute rr depending on type
    FlatMatrix rr(n, p);
    if (type == "dfbeta" || type == "dfbetas") {
      rr = mat_mat_mult(score, vbeta);
      
      if (type == "dfbetas") {
        for (int k = 0; k < p; ++k) {
          double* rcol = rr.data_ptr() + k * n;
          const double seb = std::sqrt(vbeta(k,k));
          for (int i = 0; i < n; ++i) {
            rcol[i] /= seb;
          }
        }
      }
    } else {
      rr = score;
    }
    
    if (weighted) {
      for (int k = 0; k < p; ++k) {
        double* rcol = rr.data_ptr() + k * n;
        for (int i = 0; i < n; ++i) {
          rcol[i] *= weightn[i];
        }
      }
    }
    
    if (collapse && has_id) { // collapse over id
      FlatMatrix rr1(nids,p);
      for (int k = 0; k < p; ++k) {
        const double* rcol = rr.data_ptr() + k * n;
        double* rcol1 = rr1.data_ptr() + k * nids;
        for (int i = 0; i < nids; ++i) {
          double sum = 0.0;
          for (int j = idx[i]; j < idx[i+1]; ++j) {
            sum += rcol[order[j]];
          }
          rcol1[i] = sum;
        }
      }
      rr = std::move(rr1);
    }
    
    result.push_back(std::move(rr), "resid");
    return result;
  }
  
  if (type == "schoenfeld" || type == "scaledsch") {
    // sort by stopping time in descending order within each stratum
    std::vector<int> order0 = seqcpp(0, n-1);
    std::sort(order0.begin(), order0.end(), [&](int i, int j){
      if (ignore[i] != ignore[j]) return ignore[i] < ignore[j];
      if (stratumn[i] != stratumn[j]) return stratumn[i] > stratumn[j];
      if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
      if (eventn[i] != eventn[j]) return eventn[i] < eventn[j];
      return idn[i] > idn[j];
    });
    
    std::vector<int> stratum1 = subset(stratumn, order0);
    std::vector<double> tstart1 = subset(tstartn, order0);
    std::vector<double> tstop1 = subset(tstopn, order0);
    std::vector<int> event1 = subset(eventn, order0);
    std::vector<double> weight1 = subset(weightn, order0);
    std::vector<double> offset1 = subset(offsetn, order0);
    std::vector<int> id1 = subset(idn, order0);
    std::vector<int> ignore1 = subset(ignore, order0);
    FlatMatrix z1;
    if (p > 0) z1 = subset_flatmatrix(zn, order0);
    
    // sort by starting time in descending order within each stratum
    std::vector<int> order1 = seqcpp(0, n-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j){
      if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
      if (stratum1[i] != stratum1[j]) return stratum1[i] > stratum1[j];
      return tstart1[i] > tstart1[j];
    });
    
    coxparams param = {nused, stratum1, tstart1, tstop1, event1,
                       weight1, offset1, z1, order1, method};
    
    ListCpp out = f_ressch(p, beta, &param);
    FlatMatrix rr = out.get<FlatMatrix>("resid");
    std::vector<int> index = out.get<std::vector<int>>("index");
    
    int ndead = static_cast<int>(index.size());
    std::vector<int> stratum2 = subset(stratum1, index);
    std::vector<double> time1 = subset(tstop1, index);
    
    if (weighted) {
      for (int k = 0; k < p; ++k) {
        double* rrcol = rr.data_ptr() + k * ndead;
        for (int i = 0; i < ndead; ++i) {
          rrcol[i] *= weightn[i];
        }
      }
    }
    
    if (type == "scaledsch") {
      FlatMatrix rr1 = mat_mat_mult(rr, vbeta);
      for (int k = 0; k < p; ++k) {
        double* rr1_col = rr1.data_ptr() + k * ndead;
        const double b = beta[k];
        for (int i = 0; i < ndead; ++i) {
          rr1_col[i] = rr1_col[i] * ndead + b;
        }
      }
      rr = std::move(rr1);
    }
    
    result.push_back(std::move(rr), "resid");
    result.push_back(std::move(time1), "time");
    
    if (has_stratum) {
      DataFrameCpp strata_df;
      for (int i = 0; i < p_stratum; ++i) {
        std::string s = stratum[i];
        if (u_stratum.int_cols.count(s)) {
          auto v = u_stratum.get<int>(s);
          strata_df.push_back(subset(v, stratum2), s);
        } else if (u_stratum.numeric_cols.count(s)) {
          auto v = u_stratum.get<double>(s);
          strata_df.push_back(subset(v, stratum2), s);
        } else if (u_stratum.string_cols.count(s)) {
          auto v = u_stratum.get<std::string>(s);
          strata_df.push_back(subset(v, stratum2), s);
        } else {
          throw std::invalid_argument("Unsupported type for stratum variable: " + s);
        }
      }
      result.push_back(std::move(strata_df), "strata");
    }
    
    return result;
  }
  
  throw std::invalid_argument("unknown type of residuals");
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List residuals_phregRcpp(const int p,
                               const std::vector<double>& beta,
                               const Rcpp::NumericMatrix& vbeta,
                               const std::vector<double>& resmart,
                               const Rcpp::DataFrame& data,
                               const std::vector<std::string>& stratum,
                               const std::string& time,
                               const std::string& time2,
                               const std::string& event,
                               const std::vector<std::string>& covariates,
                               const std::string& weight,
                               const std::string& offset,
                               const std::string& id,
                               const std::string& ties,
                               const std::string& type,
                               const bool collapse,
                               const bool weighted) {
  
  auto vbetacpp = flatmatrix_from_Rmatrix(vbeta);
  auto dfcpp = convertRDataFrameToCpp(data);
  
  auto cpp_result = residuals_phregcpp(
    p, beta, vbetacpp, resmart, dfcpp, stratum, time, time2, event, 
    covariates, weight, offset, id, ties, type, collapse, weighted
  );
  
  return Rcpp::wrap(cpp_result);
}


// function for individual contributions to score and info matrix for Cox model
ListCpp f_der_i_2(int p, const std::vector<double>& par, void* ex) {
  coxparams* param = (coxparams*) ex;
  const int n = param->tstop.size();
  const int nused = param->nused;
  const int method = param->method;
  
  const std::vector<int>& strata = param->strata; 
  const std::vector<double>& tstart = param->tstart; 
  const std::vector<double>& tstop = param->tstop; 
  const std::vector<int>& event = param->event; 
  const std::vector<double>& weight = param->weight; 
  const std::vector<double>& offset = param->offset; 
  const std::vector<int>& order1 = param->order1;
  const FlatMatrix& z = param->z;
  
  // Precompute eta and risk = exp(eta)
  std::vector<double> eta(nused);
  for (int person = 0; person < nused; ++person) {
    eta[person] = offset[person];
  }
  if (p > 0) {
    for (int i = 0; i < p; ++i) {
      double beta = par[i];
      if (beta == 0.0) continue;
      const double* zcol = z.data_ptr() + i * n;
      for (int person = 0; person < nused; ++person) {
        eta[person] += beta * zcol[person];
      }
    }
  }
  std::vector<double> risk(nused);
  for (int person = 0; person < nused; ++person) {
    risk[person] = std::exp(eta[person]);
  }
  
  FlatMatrix u(nused,p);         // score vector for each individual
  FlatArray imat(nused,p,p);     // information matrix for each individual
  std::vector<double> a(p);      // s1(beta,k,t)
  std::vector<double> a2(p);     // sum of w*exp(zbeta)*z for the deaths
  FlatMatrix cmat(p,p);          // s2(beta,k,t)
  FlatMatrix cmat2(p,p);         // sum of w*exp(zbeta)*z*z' for the deaths
  double denom = 0.0;            // s0(beta,k,t)
  double denom2 = 0.0;           // sum of weighted risks for deaths
  double ndead = 0.0;            // number of deaths at this time point
  
  int istrata = strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    // Reset when entering a new stratum
    if (strata[person] != istrata) {
      istrata = strata[person];
      i1 = person;
      denom = 0.0;
      
      for (int i = 0; i < p; ++i) {
        a[i] = 0.0;
      }
      
      for (int j = 0; j < p; ++j) {
        for (int i = j; i < p; ++i) {
          cmat(i,j) = 0.0;
        }
      }
    }
    
    const double dtime = tstop[person];
    
    // Process all persons tied at this dtime
    int person1 = person;
    
    for (; person < nused && tstop[person] == dtime &&
         strata[person] == istrata; ++person) {
      
      const double w = weight[person];
      const double r = w * risk[person];
      
      if (event[person] == 0) {
        denom += r;
        
        for (int i = 0; i < p; ++i) {
          a[i] += r * z(person,i);
        }
        
        for (int j = 0; j < p; ++j) {
          const double zj = z(person,j);
          for (int i = j; i < p; ++i) {
            cmat(i,j) += r * z(person,i) * zj;
          }
        }
      } else {
        ++ndead;
        denom2 += r;
        
        for (int i = 0; i < p; ++i) {
          const double zi = z(person,i);
          a2[i] += r * zi;
          u(person, i) = w * zi;
        }
        
        for (int j = 0; j < p; ++j) {
          const double zj = z(person,j);
          for (int i = j; i < p; ++i) {
            cmat2(i,j) += r * z(person,i) * zj;  
          }
        }
      }
    }
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = order1[i1];
      if (tstart[p1] < dtime || strata[p1] != istrata) break;
      
      const double r = weight[p1] * risk[p1];
      denom -= r;
      
      for (int i = 0; i < p; ++i) {
        a[i] -= r * z(p1,i);
      }
      
      for (int j = 0; j < p; ++j) {
        const double zj = z(p1,j);
        for (int i = j; i < p; ++i) {
          cmat(i,j) -= r * z(p1,i) * zj;
        }
      }
    }
    
    // Add contributions for deaths at this time
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        
        std::vector<double> xbar(p);
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          xbar[i] = a[i] / denom;
          for (int human = person1; human < person; ++human) {
            if (event[human] == 1) u(human,i) -= weight[human] * xbar[i];
          }
        }
        
        for (int j = 0; j < p; ++j) {
          for (int i = j; i < p; ++i) {
            cmat(i,j) += cmat2(i,j);
            for (int human = person1; human < person; ++human) {
              if (event[human] == 1) imat(human,i,j) = weight[human] *
                (cmat(i,j) - xbar[i] * a[j]) / denom;
            }
          }
        }
      } else { // Efron method
        const double increment = denom2 / ndead;
        for (int l = 0; l < ndead; ++l) {
          denom += increment;
          
          std::vector<double> xbar(p);
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            xbar[i] = a[i] / denom;
            for (int human = person1; human < person; ++human) {
              if (event[human] == 1) u(human,i) -= weight[human] / ndead * xbar[i];
            }
          }
          
          for (int j = 0; j < p; ++j) {
            for (int i = j; i < p; ++i) {
              cmat(i,j) += cmat2(i,j) / ndead;
              for (int human = person1; human < person; ++human) {
                if (event[human] == 1) imat(human,i,j) += weight[human] / ndead *
                  (cmat(i,j) - xbar[i] * a[j]) / denom;
              }
            }
          }
        }
      }
      
      // Reset after processing deaths
      ndead = denom2 = 0.0;
      for (int i = 0; i < p; ++i) {
        a2[i] = 0;
      }
      
      for (int j = 0; j < p; ++j) {
        for (int i = j; i < p; ++i) {
          cmat2(i,j) = 0;
        }
      }
    }
  }
  
  // fill the symmetric elements of the information matrix
  for (int j = 1; j < p; ++j)
    for (int i = 0; i < j; ++i)
      for (int person = 0; person < nused; ++person)
        imat(person,i,j) = imat(person,j,i);
  
  ListCpp result;
  result.push_back(std::move(u), "score_i");
  result.push_back(std::move(imat), "imat_i");
  return result;
}


// assess proportional hazards assumption via score process
ListCpp assess_phregcpp(const int p,
                        const std::vector<double>& beta,
                        const FlatMatrix& vbeta,
                        const DataFrameCpp& data,
                        const std::vector<std::string>& stratum ,
                        const std::string& time,
                        const std::string& time2,
                        const std::string& event,
                        const std::vector<std::string>& covariates,
                        const std::string& weight,
                        const std::string& offset,
                        const std::string& ties,
                        const int resample,
                        const int seed) {
  
  boost::random::mt19937_64 rng(seed);
  boost::random::normal_distribution<double> dist(0.0, 1.0);
  
  if (p <= 0) {
    throw std::invalid_argument(
        "covariates must be present to test proportional hazards");
  }
  
  int n = data.nrows();
  
  std::vector<int> stratumn(n);
  int p_stratum = static_cast<int>(stratum.size());
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    stratumn = out.get<std::vector<int>>("index");
  }
  
  // --- time / time2 existence and checks ---
  if (!data.containElementNamed(time))
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(timen[i]) && timen[i] < 0.0)
      throw std::invalid_argument("time must be nonnegative");
  }
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    if (data.int_cols.count(time2)) {
      const std::vector<int>& vi = data.get<int>(time2);
      for (int i = 0; i < n; ++i) time2n[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(time2)) {
      time2n = data.get<double>(time2);
    } else {
      throw std::invalid_argument("time2 variable must be integer or numeric");
    }
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(timen[i]) && !std::isnan(time2n[i]) && time2n[i] <= timen[i])
        throw std::invalid_argument("time2 must be greater than time");
    }
  }
  
  // --- event variable ---
  bool has_event = !event.empty() && data.containElementNamed(event);
  if (!has_time2 && !has_event) {
    throw std::invalid_argument(
        "data must contain the event variable for right censored data");
  }
  std::vector<int> eventn(n);
  if (has_event) {
    if (data.bool_cols.count(event)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
      for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
    } else if (data.int_cols.count(event)) {
      eventn = data.get<int>(event);
    } else if (data.numeric_cols.count(event)) {
      const std::vector<double>& vd = data.get<double>(event);
      for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
    } else {
      throw std::invalid_argument("event variable must be bool, integer or numeric");
    }
    for (double val : eventn) if (val != 0 && val != 1)
      throw std::invalid_argument("event must be 1 or 0 for each observation");
  }
  
  // --- build design matrix zn (n x p) column-major FlatMatrix ---
  FlatMatrix zn(n, p);
  if (p > 0) {
    for (int j = 0; j < p; ++j) {
      const std::string& zj = covariates[j];
      if (!data.containElementNamed(zj))
        throw std::invalid_argument("data must contain the variables in covariates");
      double* zn_col = zn.data_ptr() + j * n;
      if (data.bool_cols.count(zj)) {
        const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
        for (int i = 0; i < n; ++i) zn_col[i] = vb[i] ? 1.0 : 0.0;
      } else if (data.int_cols.count(zj)) {
        const std::vector<int>& vi = data.get<int>(zj);
        for (int i = 0; i < n; ++i) zn_col[i] = static_cast<double>(vi[i]);
      } else if (data.numeric_cols.count(zj)) {
        const std::vector<double>& vd = data.get<double>(zj);
        std::memcpy(zn_col, vd.data(), n * sizeof(double));
      } else {
        throw std::invalid_argument("covariates must be bool, integer or numeric");
      }
    }
  }
  
  // --- weight and offset ---
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const std::vector<int>& vi = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else {
      throw std::invalid_argument("weight variable must be integer or numeric");
    }
    for (double w : weightn) if (std::isnan(w) || w <= 0.0)
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    if (data.int_cols.count(offset)) {
      const std::vector<int>& vi = data.get<int>(offset);
      for (int i = 0; i < n; ++i) offsetn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(offset)) {
      offsetn = data.get<double>(offset);
    } else {
      throw std::invalid_argument("offset variable must be integer or numeric");
    }
  }
  
  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  int method = meth == "efron" ? 1 : 0;
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  // exclude observations with missing values
  std::vector<unsigned char> sub(n,1);
  for (int i = 0; i < n; ++i) {
    if (stratumn[i] == INT_MIN || std::isnan(tstartn[i]) || 
        std::isnan(tstopn[i]) || eventn[i] == INT_MIN || 
        std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
      sub[i] = 0; continue;
    }
    for (int j = 0; j < p; ++j) {
      if (std::isnan(zn(i,j))) { sub[i] = 0; break; }
    }
  }
  std::vector<int> keep = which(sub);
  if (keep.empty())
    throw std::invalid_argument("no observations without missing values");
  
  subset_in_place(stratumn, keep);
  subset_in_place(tstartn, keep);
  subset_in_place(tstopn, keep);
  subset_in_place(eventn, keep);
  subset_in_place(weightn, keep);
  subset_in_place(offsetn, keep);
  if (p > 0) subset_in_place_flatmatrix(zn, keep);
  n = keep.size();
  
  // sort by stratum
  std::vector<int> order0 = seqcpp(0, n-1);
  std::sort(order0.begin(), order0.end(), [&](int i, int j) {
    return stratumn[i] < stratumn[j];
  });
  
  std::vector<int> stratum1z = subset(stratumn, order0);
  std::vector<double> tstart1z = subset(tstartn, order0);
  std::vector<double> tstop1z = subset(tstopn, order0);
  std::vector<int> event1z = subset(eventn, order0);
  
  // locate the first observation within each stratum
  std::vector<int> istratum(1,0);
  for (int i = 1; i < n; ++i) {
    if (stratum1z[i] != stratum1z[i-1]) {
      istratum.push_back(i);
    }
  }
  int nstrata = static_cast<int>(istratum.size());
  istratum.push_back(n);
  
  // ignore subjects not at risk for any event time
  std::vector<int> ignore1z(n);
  for (int i = 0; i < nstrata; ++i) {
    int start = istratum[i], end = istratum[i+1];
    int n0 = end - start;
    std::vector<double> tstart0 = subset(tstart1z, start, end);
    std::vector<double> tstop0 = subset(tstop1z, start, end);
    std::vector<int> event0 = subset(event1z, start, end);
    
    // unique event times
    std::vector<double> etime;
    etime.reserve(n0);
    for (int j = 0; j < n0; ++j) {
      if (event0[j] == 1) etime.push_back(tstop0[j]);
    }
    etime = unique_sorted(etime);
    
    std::vector<int> index1 = findInterval3(tstart0, etime);
    std::vector<int> index2 = findInterval3(tstop0, etime);
    for (int j = istratum[i]; j < istratum[i+1]; ++j) {
      int j0 = j - istratum[i];
      if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
        ignore1z[j] = 1;
      } else {
        ignore1z[j] = 0;
      }
    }
  }
  
  std::vector<int> ignoren(n); // back to the original order
  for (int i = 0; i < n; ++i) {
    ignoren[order0[i]] = ignore1z[i];
  }
  
  int nused = n - std::accumulate(ignoren.begin(), ignoren.end(), 0);
  
  // sort by stopping time in descending order within each stratum
  std::vector<int> order1 = seqcpp(0, n-1);
  std::sort(order1.begin(), order1.end(), [&](int i, int j) {
    if (ignoren[i] != ignoren[j]) return ignoren[i] < ignoren[j];
    if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
    if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
    return eventn[i] < eventn[j];
  });
  
  std::vector<int> stratum1 = subset(stratumn, order1);
  std::vector<double> tstart1 = subset(tstartn, order1);
  std::vector<double> tstop1 = subset(tstopn, order1);
  std::vector<int> event1 = subset(eventn, order1);
  std::vector<double> weight1 = subset(weightn, order1);
  std::vector<double> offset1 = subset(offsetn, order1);
  std::vector<int> ignore1 = subset(ignoren, order1);
  FlatMatrix z1;
  if (p > 0) z1 = subset_flatmatrix(zn, order1);
  
  // sort by starting time in descending order within each stratum
  std::vector<int> order10 = seqcpp(0, n-1);
  std::sort(order10.begin(), order10.end(), [&](int i, int j) {
    if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
    if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
    return tstart1[i] > tstart1[j];
  });
  
  coxparams param = {nused, stratum1, tstart1, tstop1, event1,
                     weight1, offset1, z1, order10, method};
  
  ListCpp out = f_der_i_2(p, beta, &param);
  FlatMatrix score_i = out.get<FlatMatrix>("score_i");
  FlatArray imat_i = out.get<FlatArray>("imat_i");
  
  // order data by ascending tstop
  std::vector<int> order2 = seqcpp(0, n-1);
  std::sort(order2.begin(), order2.end(), [&](int i, int j) {
    return tstop1[i] < tstop1[j];
  });
  
  std::vector<double> tstop2 = subset(tstop1, order2);
  std::vector<int> event2 = subset(event1, order2);
  FlatMatrix score_i2 = subset_flatmatrix(score_i, order2);
  FlatArray imat_i2 = subset_flatarray(imat_i, order2);
  
  // only consider event times
  std::vector<int> q;
  q.reserve(n);
  for (int i = 0; i < n; ++i) if (event2[i] == 1) q.push_back(i);
  
  std::vector<double> tstop3 = subset(tstop2, q);
  FlatMatrix score_i3 = subset_flatmatrix(score_i2, q);
  FlatArray imat_i3 = subset_flatarray(imat_i2, q);
  int d = static_cast<int>(q.size());
  
  // identify unique event times
  std::vector<double> t = unique_sorted(tstop3);
  
  int nt = static_cast<int>(t.size());
  FlatMatrix score_i4(nt, p);
  FlatArray imat_i4(nt, p, p);
  
  std::vector<double> tstop4(d+1);
  tstop4[0] = -1.0;
  std::memcpy(tstop4.data() + 1, tstop3.data(), d * sizeof(double));
  
  // sums over unique event times
  for (int j = 0; j < p; ++j) {
    int h = -1;
    for (int i = 1; i <= d; ++i) {
      if (tstop4[i] != tstop4[i-1]) ++h;
      score_i4(h,j) += score_i3(i-1,j);
    }
  }
  
  for (int k = 0; k < p; ++k) {
    for (int j = 0; j < p; ++j) {
      int h = -1;
      for (int i = 1; i <= d; ++i) {
        if (tstop4[i] != tstop4[i-1]) ++h;
        imat_i4(h,j,k) += imat_i3(i-1,j,k);
      }
    }
  }
  
  // cumulative sums of score and information processes
  for (int j = 0; j < p; ++j) {
    for (int i = 1; i < nt; ++i) {
      score_i4(i,j) += score_i4(i-1,j);
    }
  }
  
  for (int k = 0; k < p; ++k) {
    for (int j = 0; j < p; ++j) {
      for (int i = 1; i < nt; ++i) {
        imat_i4(i,j,k) += imat_i4(i-1,j,k);
      }
    }
  }
  
  // standardize the score processes
  for (int i = 0; i < nt; ++i) {
    for (int j = 0; j < p; ++j) {
      score_i4(i,j) *= std::sqrt(vbeta(j,j));
    }
  }
  
  // resampling to approximate the null distribution
  FlatMatrix G(d, resample);
  for (int r = 0; r < resample; ++r) {
    for (int i = 0; i < d; ++i) {
      G(i,r) = dist(rng);
    }
  }
  
  // sum(U * G)
  FlatMatrix U(p, resample);
  for (int r = 0; r < resample; ++r) {
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < d; ++i) {
        U(j,r) += score_i3(i,j) * G(i,r);
      }
    }
  }
  
  // sum(I(X <= t) * U * G)
  FlatArray score_i_vec(nt, p, resample);
  for (int r = 0; r < resample; ++r) {
    // sums of unique event times
    for (int j = 0; j < p; ++j) {
      int h = -1;
      for (int i = 1; i <= d; ++i) {
        if (tstop4[i] != tstop4[i-1]) ++h;
        score_i_vec(h,j,r) += score_i3(i-1,j) * G(i-1,r);
      }
    }
    
    // cumulative sums of score processes
    for (int j = 0; j < p; ++j) {
      for (int i = 1; i < nt; ++i) {
        score_i_vec(i,j,r) += score_i_vec(i-1,j,r);
      }
    }
  }
  
  // piece together the result
  FlatMatrix vU = mat_mat_mult(vbeta, U);
  FlatArray ivU(nt, p, resample);
  for (int r = 0; r < resample; ++r) {
    for (int k = 0; k < p; ++k) {
      for (int j = 0; j < p; ++j) {
        for (int i = 0; i < nt; ++i) {
          ivU(i,j,r) += imat_i4(i,j,k) * vU(k,r);
        }
      }
    }
  }
  FlatArray score_i_vec_2(nt, p, resample);
  for (int r = 0; r < resample; ++r) {
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < nt; ++i) {
        score_i_vec_2(i,j,r) = score_i_vec(i,j,r) - ivU(i,j,r);
      }
    }
  }
  
  // standardize the resampled score processes
  for (int r = 0;  r <resample; ++r) {
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < nt; ++i) {
        score_i_vec_2(i,j,r) *= std::sqrt(vbeta(j,j));
      }
    }
  }
  
  // calculate individual p-values for proportional hazards assumption test
  std::vector<double> p_values(p+1);
  
  // observed test statistics for each covariate and overall
  std::vector<double> tobs(p+1);
  for (int j = 0; j < p; ++j) {
    std::vector<double> obs(nt);
    for (int i = 0; i < nt; ++i) {
      obs[i] = std::fabs(score_i4(i,j));
    }
    tobs[j] = *std::max_element(obs.begin(), obs.end());
  }
  
  std::vector<double> obs(nt);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < nt; ++i) {
      obs[i] += std::fabs(score_i4(i,j));
    }
  }
  tobs[p] = *std::max_element(obs.begin(), obs.end());
  
  // count exceedances in resampled test statistics
  std::vector<double> count(p+1);
  for (int r = 0; r < resample; ++r) {
    for (int j = 0; j < p; ++j) {
      std::vector<double> sim(nt);
      for (int i = 0; i < nt; ++i) {
        sim[i] = std::fabs(score_i_vec_2(i,j,r));
      }
      double tsim = *std::max_element(sim.begin(), sim.end());
      if (tsim >= tobs[j]) ++count[j];
    }
    
    std::vector<double> sim(nt);
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < nt; ++i) {
        sim[i] += std::fabs(score_i_vec_2(i,j,r));
      }
    }
    double tsim = *std::max_element(sim.begin(), sim.end());
    if (tsim >= tobs[p]) ++count[p];
  }
  
  // calculate p-values
  for (int j = 0; j <= p; ++j) {
    p_values[j] = count[j] / resample;
  }
  
  ListCpp result;
  result.push_back(std::move(t), "time");
  result.push_back(std::move(score_i4), "score_t");
  result.push_back(std::move(score_i_vec_2), "score_t_list");
  result.push_back(std::move(tobs), "max_abs_value");
  result.push_back(std::move(p_values), "p_value");
  result.push_back(resample, "resample");
  result.push_back(seed, "seed");
  return result;
}

// [[Rcpp::export]]
Rcpp::List assess_phregRcpp(const int p,
                            const std::vector<double>& beta,
                            const Rcpp::NumericMatrix& vbeta,
                            const Rcpp::DataFrame& data,
                            const std::vector<std::string>& stratum ,
                            const std::string& time,
                            const std::string& time2,
                            const std::string& event,
                            const std::vector<std::string>& covariates,
                            const std::string& weight,
                            const std::string& offset,
                            const std::string& ties,
                            const int resample,
                            const std::uint32_t seed) {
  
  auto vbetacpp = flatmatrix_from_Rmatrix(vbeta);
  auto dfcpp = convertRDataFrameToCpp(data);
  
  auto cpp_result = assess_phregcpp(
    p, beta, vbetacpp, dfcpp, stratum, time, time2, event,
    covariates, weight, offset, ties, resample, seed
  );
  
  return Rcpp::wrap(cpp_result);
}


// test proportional hazards assumption via scaled Schoenfeld residuals
ListCpp zph_phregcpp(int p,
                     const std::vector<double>& beta,
                     const FlatMatrix& vbeta,
                     const std::vector<double>& resmart,
                     const DataFrameCpp& data,
                     const std::vector<std::string>& stratum,
                     const std::string& time,
                     const std::string& time2,
                     const std::string& event,
                     const std::vector<std::string>& covariates,
                     const std::string& weight,
                     const std::string& offset,
                     const std::string& ties,
                     const std::string& transform) {
  
  if (p <= 0) {
    throw std::invalid_argument(
        "covariates must be present to test proportional hazards");
  }
  
  int n = data.nrows();
  
  bool has_stratum = false;
  std::vector<int> stratumn(n);
  DataFrameCpp u_stratum;
  int p_stratum = stratum.size();
  if (!(p_stratum == 0 || (p_stratum == 1 && stratum[0] == ""))) {
    ListCpp out = bygroup(data, stratum);
    has_stratum = true;
    stratumn = out.get<std::vector<int>>("index");
    u_stratum = out.get<DataFrameCpp>("lookup");
  }
  
  // --- time / time2 existence and checks ---
  if (!data.containElementNamed(time))
    throw std::invalid_argument("data must contain the time variable");
  std::vector<double> timen(n);
  if (data.int_cols.count(time)) {
    const std::vector<int>& vi = data.get<int>(time);
    for (int i = 0; i < n; ++i) timen[i] = static_cast<double>(vi[i]);
  } else if (data.numeric_cols.count(time)) {
    timen = data.get<double>(time);
  } else {
    throw std::invalid_argument("time variable must be integer or numeric");
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isnan(timen[i]) && timen[i] < 0.0)
      throw std::invalid_argument("time must be nonnegative");
  }
  
  bool has_time2 = !time2.empty() && data.containElementNamed(time2);
  std::vector<double> time2n(n);
  if (has_time2) {
    if (data.int_cols.count(time2)) {
      const std::vector<int>& vi = data.get<int>(time2);
      for (int i = 0; i < n; ++i) time2n[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(time2)) {
      time2n = data.get<double>(time2);
    } else {
      throw std::invalid_argument("time2 variable must be integer or numeric");
    }
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(timen[i]) && !std::isnan(time2n[i]) && time2n[i] <= timen[i])
        throw std::invalid_argument("time2 must be greater than time");
    }
  }
  
  // --- event variable ---
  bool has_event = !event.empty() && data.containElementNamed(event);
  if (!has_time2 && !has_event) {
    throw std::invalid_argument(
        "data must contain the event variable for right censored data");
  }
  std::vector<int> eventn(n);
  if (has_event) {
    if (data.bool_cols.count(event)) {
      const std::vector<unsigned char>& vb = data.get<unsigned char>(event);
      for (int i = 0; i < n; ++i) eventn[i] = vb[i] ? 1 : 0;
    } else if (data.int_cols.count(event)) {
      eventn = data.get<int>(event);
    } else if (data.numeric_cols.count(event)) {
      const std::vector<double>& vd = data.get<double>(event);
      for (int i = 0; i < n; ++i) eventn[i] = static_cast<int>(vd[i]);
    } else {
      throw std::invalid_argument("event variable must be bool, integer or numeric");
    }
    for (double val : eventn) if (val != 0 && val != 1)
      throw std::invalid_argument("event must be 1 or 0 for each observation");
  }
  
  // --- build design matrix zn (n x p) column-major FlatMatrix ---
  FlatMatrix zn(n, p);
  if (p > 0) {
    for (int j = 0; j < p; ++j) {
      const std::string& zj = covariates[j];
      if (!data.containElementNamed(zj))
        throw std::invalid_argument("data must contain the variables in covariates");
      double* zn_col = zn.data_ptr() + j * n;
      if (data.bool_cols.count(zj)) {
        const std::vector<unsigned char>& vb = data.get<unsigned char>(zj);
        for (int i = 0; i < n; ++i) zn_col[i] = vb[i] ? 1.0 : 0.0;
      } else if (data.int_cols.count(zj)) {
        const std::vector<int>& vi = data.get<int>(zj);
        for (int i = 0; i < n; ++i) zn_col[i] = static_cast<double>(vi[i]);
      } else if (data.numeric_cols.count(zj)) {
        const std::vector<double>& vd = data.get<double>(zj);
        std::memcpy(zn_col, vd.data(), n * sizeof(double));
      } else {
        throw std::invalid_argument("covariates must be bool, integer or numeric");
      }
    }
  }
  
  // --- weight and offset ---
  std::vector<double> weightn(n, 1.0);
  if (!weight.empty() && data.containElementNamed(weight)) {
    if (data.int_cols.count(weight)) {
      const std::vector<int>& vi = data.get<int>(weight);
      for (int i = 0; i < n; ++i) weightn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(weight)) {
      weightn = data.get<double>(weight);
    } else {
      throw std::invalid_argument("weight variable must be integer or numeric");
    }
    for (double w : weightn) if (std::isnan(w) || w <= 0.0)
      throw std::invalid_argument("weight must be greater than 0");
  }
  
  std::vector<double> offsetn(n, 0.0);
  if (!offset.empty() && data.containElementNamed(offset)) {
    if (data.int_cols.count(offset)) {
      const std::vector<int>& vi = data.get<int>(offset);
      for (int i = 0; i < n; ++i) offsetn[i] = static_cast<double>(vi[i]);
    } else if (data.numeric_cols.count(offset)) {
      offsetn = data.get<double>(offset);
    } else {
      throw std::invalid_argument("offset variable must be integer or numeric");
    }
  }
  
  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  int method = meth == "efron" ? 1 : 0;
  
  // unify right censored data with counting process data
  std::vector<double> tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  // exclude observations with missing values
  std::vector<unsigned char> sub(n,1);
  for (int i = 0; i < n; ++i) {
    if (stratumn[i] == INT_MIN || std::isnan(tstartn[i]) || 
        std::isnan(tstopn[i]) || eventn[i] == INT_MIN || 
        std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
      sub[i] = 0; continue;
    }
    for (int j = 0; j < p; ++j) {
      if (std::isnan(zn(i,j))) { sub[i] = 0; break; }
    }
  }
  std::vector<int> keep = which(sub);
  if (keep.empty())
    throw std::invalid_argument("no observations without missing values");
  
  subset_in_place(stratumn, keep);
  subset_in_place(tstartn, keep);
  subset_in_place(tstopn, keep);
  subset_in_place(eventn, keep);
  subset_in_place(weightn, keep);
  subset_in_place(offsetn, keep);
  if (p > 0) subset_in_place_flatmatrix(zn, keep);
  n = keep.size();
  
  // sort by stratum
  std::vector<int> order0 = seqcpp(0, n-1);
  std::sort(order0.begin(), order0.end(), [&](int i, int j) {
    return stratumn[i] < stratumn[j];
  });
  
  std::vector<int> stratum1z = subset(stratumn, order0);
  std::vector<double> tstart1z = subset(tstartn, order0);
  std::vector<double> tstop1z = subset(tstopn, order0);
  std::vector<int> event1z = subset(eventn, order0);
  
  // locate the first observation within each stratum
  std::vector<int> istratum(1,0);
  for (int i = 1; i < n; ++i) {
    if (stratum1z[i] != stratum1z[i-1]) {
      istratum.push_back(i);
    }
  }
  int nstrata = static_cast<int>(istratum.size());
  istratum.push_back(n);
  
  // ignore subjects not at risk for any event time
  std::vector<int> ignore1z(n);
  for (int i = 0; i < nstrata; ++i) {
    int start = istratum[i], end = istratum[i+1];
    int n0 = end - start;
    std::vector<double> tstart0 = subset(tstart1z, start, end);
    std::vector<double> tstop0 = subset(tstop1z, start, end);
    std::vector<int> event0 = subset(event1z, start, end);
    
    // unique event times
    std::vector<double> etime;
    etime.reserve(n0);
    for (int j = 0; j < n0; ++j) {
      if (event0[j] == 1) etime.push_back(tstop0[j]);
    }
    etime = unique_sorted(etime);
    
    std::vector<int> index1 = findInterval3(tstart0, etime);
    std::vector<int> index2 = findInterval3(tstop0, etime);
    for (int j = istratum[i]; j < istratum[i+1]; ++j) {
      int j0 = j - istratum[i];
      if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
        ignore1z[j] = 1;
      } else {
        ignore1z[j] = 0;
      }
    }
  }
  
  std::vector<int> ignoren(n); // back to the original order
  for (int i = 0; i < n; ++i) {
    ignoren[order0[i]] = ignore1z[i];
  }
  
  int nused = n - std::accumulate(ignoren.begin(), ignoren.end(), 0);
  
  // sort by stopping time in descending order within each stratum
  std::vector<int> order1 = seqcpp(0, n-1);
  std::sort(order1.begin(), order1.end(), [&](int i, int j) {
    if (ignoren[i] != ignoren[j]) return ignoren[i] < ignoren[j];
    if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
    if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
    return eventn[i] < eventn[j];
  });
  
  std::vector<int> stratum1 = subset(stratumn, order1);
  std::vector<double> tstart1 = subset(tstartn, order1);
  std::vector<double> tstop1 = subset(tstopn, order1);
  std::vector<int> event1 = subset(eventn, order1);
  std::vector<double> weight1 = subset(weightn, order1);
  std::vector<double> offset1 = subset(offsetn, order1);
  std::vector<int> ignore1 = subset(ignoren, order1);
  FlatMatrix z1;
  if (p > 0) z1 = subset_flatmatrix(zn, order1);
  
  // sort by starting time in descending order within each stratum
  std::vector<int> order10 = seqcpp(0, n-1);
  std::sort(order10.begin(), order10.end(), [&](int i, int j) {
    if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
    if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
    return tstart1[i] > tstart1[j];
  });
  
  coxparams param = {nused, stratum1, tstart1, tstop1, event1,
                     weight1, offset1, z1, order10, method};
  
  ListCpp out = f_der_i_2(p, beta, &param);
  FlatMatrix score_i = out.get<FlatMatrix>("score_i");
  FlatArray imat_i = out.get<FlatArray>("imat_i");
  
  // order data by ascending tstop
  std::vector<int> order2 = seqcpp(0, n-1);
  std::sort(order2.begin(), order2.end(), [&](int i, int j) {
    return tstop1[i] < tstop1[j];
  });
  
  std::vector<int> stratum2 = subset(stratum1, order2);
  std::vector<double> tstop2 = subset(tstop1, order2);
  std::vector<int> event2 = subset(event1, order2);
  FlatMatrix score_i2 = subset_flatmatrix(score_i, order2);
  FlatArray imat_i2 = subset_flatarray(imat_i, order2);
  
  // only consider event times
  std::vector<int> q;
  q.reserve(n);
  for (int i = 0; i < n; ++i) if (event2[i] == 1) q.push_back(i);
  
  std::vector<int> stratum3 = subset(stratum2, q);
  std::vector<double> tstop3 = subset(tstop2, q);
  FlatMatrix score_i3 = subset_flatmatrix(score_i2, q);
  FlatArray imat_i3 = subset_flatarray(imat_i2, q);
  int d = static_cast<int>(q.size());
  
  // identify unique event times
  std::vector<double> t = unique_sorted(tstop3);
  
  int nt = static_cast<int>(t.size());
  FlatMatrix score_i4(nt, p);
  FlatArray imat_i4(nt, p, p);
  
  std::vector<double> tstop4(d+1);
  tstop4[0] = -1.0;
  std::memcpy(tstop4.data() + 1, tstop3.data(), d * sizeof(double));
  
  // sums over unique event times
  for (int j = 0; j < p; ++j) {
    int h = -1;
    for (int i = 1; i <= d; ++i) {
      if (tstop4[i] != tstop4[i-1]) ++h;
      score_i4(h,j) += score_i3(i-1,j);
    }
  }
  
  for (int k = 0; k < p; ++k) {
    for (int j = 0; j < p; ++j) {
      int h = -1;
      for (int i = 1; i <= d; ++i) {
        if (tstop4[i] != tstop4[i-1]) ++h;
        imat_i4(h,j,k) += imat_i3(i-1,j,k);
      }
    }
  }
  
  // transformed time points
  std::vector<double> g(nt);
  if (transform == "identity") {
    for (int i = 0; i < nt; ++i) g[i] = t[i];
  } else if (transform == "log") {
    for (int i = 0; i < nt; ++i) g[i] = std::log(t[i]);
  } else if (transform == "rank") {
    for (int i = 0; i < nt; ++i) g[i] = static_cast<double>(i+1);
  } else if (transform == "km") {
    DataFrameCpp temp;
    temp.push_back(std::move(tstart1), "tstart");
    temp.push_back(std::move(tstop1), "tstop");
    temp.push_back(std::move(event1), "event");
    
    DataFrameCpp df_km = kmestcpp(temp, {""}, "tstart", "tstop", "event", 
                                  "", "none", 0.95, false);
    
    std::vector<double> surv = df_km.get<double>("surv");
    surv.insert(surv.begin(), 1.0); // for left-continuous step function
    for (int i = 0; i < nt; ++i) g[i] = 1.0 - surv[i];
  } else {
    throw std::invalid_argument("Unsupported transform: " + transform);
  }
  
  // score for theta
  std::vector<double> u_theta(p);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < nt; ++i) {
      u_theta[j] += score_i4(i,j) * g[i];
    }
  }
  
  // covariance between u_theta and u_beta
  FlatMatrix imat_theta_beta(p, p);
  for (int k = 0; k < p; ++k) {
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < nt; ++i) {
        imat_theta_beta(j,k) += imat_i4(i,j,k) * g[i];
      }
    }
  }
  
  // covariance for u_theta
  FlatMatrix imat_theta(p, p);
  for (int k = 0; k < p; ++k) {
    for (int j = 0; j < p; ++j) {
      for (int i = 0; i < nt; ++i) {
        imat_theta(j,k) += imat_i4(i,j,k) * g[i] * g[i];
      }
    }
  }
  
  // conditional variance for u_theta given u_beta
  FlatMatrix ivbeta = mat_mat_mult(imat_theta_beta, vbeta);
  
  FlatMatrix vtheta(p, p);
  for (int s = 0; s < p; ++s) {
    for (int k = 0; k < p; ++k) {
      for (int j = 0; j < p; ++j) {
        vtheta(j,k) -= ivbeta(j,s) * imat_theta_beta(k,s);
      }
    }
  }
  
  for (int k = 0; k < p; ++k) {
    for (int j = 0; j < p; ++j) {
      vtheta(j,k) += imat_theta(j,k);
    }
  }
  
  // individual score test for theta = 0
  std::vector<double> score_test(p);
  for (int j = 0; j < p; ++j) {
    score_test[j] = u_theta[j] * u_theta[j] / vtheta(j,j);
  }
  
  // global score test for theta = 0
  FlatMatrix vtheta_inv = invsympd(vtheta, p);
  double score_test_global = quadsym(u_theta, vtheta_inv);
  
  FlatMatrix table(p+1,3);
  for (int j = 0; j < p; ++j) {
    table(j,0) = score_test[j];
    table(j,1) = 1.0;
    table(j,2) = 1.0 - boost_pchisq(score_test[j], 1.0);
  }
  table(p,0) = score_test_global;
  table(p,1) = static_cast<double>(p);
  table(p,2) = 1.0 - boost_pchisq(score_test_global, static_cast<double>(p));
  
  // obtain scaled schoenfeld residuals
  ListCpp resid_list = residuals_phregcpp(
    p, beta, vbeta, resmart, data, stratum, time, time2, event,
    covariates, weight, offset, "", ties, "scaledsch", false, true);
  
  FlatMatrix sresid = resid_list.get<FlatMatrix>("resid");
  
  // repeat the g values according to the number of events
  std::vector<double> g_rep; 
  g_rep.reserve(d);
  for (int i = 0; i < nt; ++i) {
    int count = 0;
    for (int j = 0; j < d; ++j) if (tstop3[j] == t[i]) ++count;
    for (int k = 0; k < count; ++k) g_rep.push_back(g[i]);
  }
  FlatMatrix var = vbeta;
  for (double &x : var.data) x *= static_cast<double>(d);
  
  ListCpp result;
  result.push_back(std::move(table), "table");
  result.push_back(std::move(g_rep), "x");
  result.push_back(std::move(tstop3), "time");
  result.push_back(std::move(sresid), "y");
  result.push_back(std::move(var), "var");
  result.push_back(transform, "transform");
  
  if (has_stratum) {
    for (int i = 0; i < p_stratum; ++i) {
      std::string s = stratum[i];
      if (u_stratum.int_cols.count(s)) {
        auto v = u_stratum.get<int>(s);
        result.push_back(subset(v, stratum3), s);
      } else if (u_stratum.numeric_cols.count(s)) {
        auto v = u_stratum.get<double>(s);
        result.push_back(subset(v, stratum3), s);
      } else if (u_stratum.string_cols.count(s)) {
        auto v = u_stratum.get<std::string>(s);
        result.push_back(subset(v, stratum3), s);
      } else {
        throw std::invalid_argument("unsupported type for stratum variable " + s);
      }
    }
  }
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List zph_phregRcpp(int p,
                         const std::vector<double>& beta,
                         const Rcpp::NumericMatrix& vbeta,
                         const std::vector<double>& resmart,
                         const Rcpp::DataFrame& data,
                         const std::vector<std::string>& stratum,
                         const std::string& time,
                         const std::string& time2,
                         const std::string& event,
                         const std::vector<std::string>& covariates,
                         const std::string& weight,
                         const std::string& offset,
                         const std::string& ties,
                         const std::string& transform) {
  
  auto vbetacpp = flatmatrix_from_Rmatrix(vbeta);
  auto dfcpp = convertRDataFrameToCpp(data);
  
  auto cpp_result = zph_phregcpp(
    p, beta, vbetacpp, resmart, dfcpp, stratum, time, time2, event,
    covariates, weight, offset, ties, transform
  );
  
  return Rcpp::wrap(cpp_result);
}
