#include <Rcpp.h>
#include <R_ext/Applic.h>
#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


NumericVector fsurvci(double surv, double sesurv, std::string ct, double z) {
  double grad, hw, lower = NA_REAL, upper = NA_REAL;
  if (surv == 1.0 && sesurv == 0.0) {
    lower = upper = 1.0;
  } else if (ct == "plain" || ct == "linear") {
    lower = std::max(surv - z * sesurv, 0.0);
    upper = std::min(surv + z * sesurv, 1.0);
  } else if (ct == "log") {
    grad = 1.0 / surv;
    hw = z * grad * sesurv;
    lower = exp(log(surv) - hw);
    upper = std::min(exp(log(surv) + hw), 1.0);
  } else if (ct == "log-log" || ct == "loglog" || ct == "cloglog") {
    grad = 1.0 / (surv * log(surv));
    hw = z * grad * sesurv;
    lower = exp(-exp(log(-log(surv)) - hw));
    upper = exp(-exp(log(-log(surv)) + hw));
  } else if (ct == "logit") {
    grad = 1.0 / (surv * (1.0 - surv));
    hw = z * grad * sesurv;
    lower = R::plogis(R::qlogis(surv, 0, 1, 1, 0) - hw, 0, 1, 1, 0);
    upper = R::plogis(R::qlogis(surv, 0, 1, 1, 0) + hw, 0, 1, 1, 0);
  } else if (ct == "arcsin" || ct == "asin" || ct == "asinsqrt") {
    grad = 1.0 / (2.0 * sqrt(surv * (1.0 - surv)));
    hw = z * grad * sesurv;
    lower = pow(sin(asin(sqrt(surv)) - hw), 2);
    upper = pow(sin(asin(sqrt(surv)) + hw), 2);
  } else {
    stop("Unknown confidence type: " + ct);
  }
  
  return NumericVector::create(lower, upper);
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
DataFrame survQuantile(const NumericVector& time = NA_REAL,
                       const IntegerVector& event = NA_INTEGER,
                       const double cilevel = 0.95,
                       const std::string transform = "loglog",
                       const NumericVector& probs = NA_REAL) {
  
  if (is_true(any(is_na(time)))) {
    stop("time must be provided");
  }
  
  if (is_true(any(is_na(event)))) {
    stop("event must be provided");
  }
  
  if (is_true(any((event != 1) & (event != 0)))) {
    stop("event must be 1 or 0");
  }
  
  if (is_true(all(event == 0))) {
    stop("at least 1 event is needed");
  }
  
  if (cilevel <= 0 || cilevel >= 1) {
    stop("cilevel must lie between 0 and 1");
  }
  
  std::string ct = transform;
  std::for_each(ct.begin(), ct.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (!(ct == "linear" || ct == "plain" || ct == "log" ||
      ct == "loglog" || ct == "log-log" || ct == "cloglog" ||
      ct == "asinsqrt" || ct == "arcsin"|| ct == "asin" ||
      ct == "logit")) {
    stop("Invalid value for transform");
  }
  
  int code;
  if (ct == "plain" || ct == "linear") {
    code = 1;
  } else if (ct == "log") {
    code = 2;
  } else if (ct == "log-log" || ct == "loglog" || ct == "cloglog") {
    code = 3;
  } else if (ct == "logit") {
    code = 4;
  } else if (ct == "arcsin" || ct == "asin" || ct == "asinsqrt") {
    code = 5;
  } else {
    Rcpp::stop("Unknown confidence type: " + ct);
  }
  
  if (is_true(any((probs <= 0) | (probs >= 1)))) {
    stop("Elements of probs must lie between 0 and 1");
  }
  
  int n = static_cast<int>(time.size());
  
  // sort by time, and event with event in descending order
  IntegerVector order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (time[i] != time[j]) return time[i] < time[j];
    return event[i] > event[j];
  });
  
  NumericVector time2 = time[order];
  IntegerVector event2 = event[order];
  
  NumericVector time0(n, NA_REAL), nrisk0(n), nevent0(n);
  NumericVector surv0(n), sesurv0(n);
  double t = 0, nrisk = n, nevent = 0, surv = 1, vcumhaz = 0, sesurv;
  
  int i = 0;
  bool cache = false; // buffer for the current event time
  for (int j=0; j<n; ++j) {
    if ((j == 0 && event2[j] == 1) ||
        (j >= 1 && event2[j] == 1 && time2[j] > time2[j-1])) {
      // new event, add the info for the previous event
      if (cache) {
        surv *= (1.0 - nevent/nrisk);
        if (nrisk > nevent) {
          vcumhaz += nevent/(nrisk*(nrisk - nevent));
        } else {
          vcumhaz = NA_REAL;
        }
        sesurv = surv*sqrt(vcumhaz);
        
        time0[i] = t;
        nrisk0[i] = nrisk;
        nevent0[i] = nevent;
        surv0[i] = surv;
        sesurv0[i] = sesurv;
        
        ++i;
      }
      
      // update the buffer for the current event time
      t = time2[j];
      nrisk = n-j;
      nevent = 1;
      cache = true;
    } else if (j >= 1 && event2[j] == 1 && event2[j-1] == 1 &&
      time2[j] == time2[j-1]) { // tied event
      ++nevent;
    } else if (j >= 1 && event2[j] == 0 && event2[j-1] == 1) {
      // new censoring, add the info for the previous event
      surv *= (1.0 - nevent/nrisk);
      if (nrisk > nevent) {
        vcumhaz += nevent/(nrisk*(nrisk - nevent));
      } else {
        vcumhaz = NA_REAL;
      }
      sesurv = surv*sqrt(vcumhaz);
      
      time0[i] = t;
      nrisk0[i] = nrisk;
      nevent0[i] = nevent;
      surv0[i] = surv;
      sesurv0[i] = sesurv;
      
      ++i;
      
      // empty the buffer for the current event time
      cache = false;
    }
  }
  
  // add the info for the last event
  if (cache) {
    surv *= (1.0 - nevent/nrisk);
    if (nrisk > nevent) {
      vcumhaz += nevent/(nrisk*(nrisk - nevent));
    } else {
      vcumhaz = NA_REAL;
    }
    sesurv = surv*sqrt(vcumhaz);
    
    time0[i] = t;
    nrisk0[i] = nrisk;
    nevent0[i] = nevent;
    surv0[i] = surv;
    sesurv0[i] = sesurv;
    
    ++i;
  }
  
  // only keep nonmissing records
  LogicalVector sub = !is_na(time0);
  time0 = time0[sub];
  nrisk0 = nrisk0[sub];
  nevent0 = nevent0[sub];
  surv0 = surv0[sub];
  sesurv0 = sesurv0[sub];
  
  int n0 = sum(sub);
  NumericVector z(n0, NA_REAL), grad(n0, NA_REAL);
  double zcrit = R::qnorm((1+cilevel)/2,0,1,1,0);
  
  int m = static_cast<int>(probs.size());
  NumericVector quantile(m), lower(m), upper(m);
  for (int j=0; j<m; ++j) {
    double q = 1.0 - probs[j];
    for (int i=0; i<n0; ++i) {
      if (nrisk0[i] > nevent0[i]) {
        switch (code) {
        case 1: // linear or plain
          z[i] = (surv0[i] - q) / sesurv0[i];
          break;
        
        case 2: // log
          grad[i] = 1.0 / surv0[i];
          z[i] = (log(surv0[i]) - log(q)) / (grad[i] * sesurv0[i]);
          break;
        
        case 3: // loglog / log-log / cloglog
          grad[i] = 1.0 / (surv0[i] * log(surv0[i]));
          z[i] = (log(-log(surv0[i])) - log(-log(q))) / (grad[i] * sesurv0[i]);
          break;
          
        case 4: // logit
          grad[i] = 1.0 / (surv0[i] * (1.0 - surv0[i]));
          z[i] = (R::qlogis(surv0[i], 0, 1, 1, 0) - R::qlogis(q, 0, 1, 1, 0)) /
            (grad[i] * sesurv0[i]);
          break;
        
        case 5: // arcsin / asin / asinsqrt
          grad[i] = 1.0 / (2.0 * sqrt(surv0[i] * (1.0 - surv0[i])));
          z[i] = (asin(sqrt(surv0[i])) - asin(sqrt(q))) / 
            (grad[i] * sesurv0[i]);
          break;
        }
      }
    }
    
    IntegerVector index;
    for (int i = 0; i < z.size(); ++i) {
      if (!std::isnan(z[i]) && std::abs(z[i]) <= zcrit) {
        index.push_back(i);
      }
    }
    
    if (index.size() == 0) {
      lower[j] = upper[j] = NA_REAL;
    } else {
      lower[j] = time0[min(index)];
      upper[j] = max(index) < n0-1 ? time0[max(index)+1] : NA_REAL;
    }
    
    if (is_true(any(surv0 < q))) {
      quantile[j] = time0[min(which(surv0 < q))];
    } else {
      quantile[j] = NA_REAL;
    }
  }
  
  DataFrame result = DataFrame::create(
    Named("prob") = probs,
    Named("quantile") = quantile,
    Named("lower") = lower,
    Named("upper") = upper,
    Named("cilevel") = cilevel,
    Named("transform") = ct
  );
  
  return result;
}


//' @title Kaplan-Meier Estimates of Survival Curve
//' @description Obtains the Kaplan-Meier estimates of the survival curve.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{rep}: The replication for by-group processing.
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
//' @param rep The name(s) of the replication variable(s) in the input data.
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
//' * \code{rep}: The replication.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' kmest(data = aml, stratum = "x", time = "time", event = "status")
//'
//' @export
// [[Rcpp::export]]
DataFrame kmest(const DataFrame data,
                const StringVector& rep = "",
                const StringVector& stratum = "",
                const std::string time = "time",
                const std::string time2 = "",
                const std::string event = "event",
                const std::string weight = "",
                const std::string conftype = "log-log",
                const double conflev = 0.95,
                const bool keep_censor = false) {
  
  int n = data.nrows();
  
  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(0);
  } else {
    List out = bygroup(data, rep);
    has_rep = 1;
    repn = out["index"];
    u_rep = DataFrame(out["lookup"]);
  }
  
  bool has_stratum;
  IntegerVector stratumn(n);
  DataFrame u_stratum;
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(0);
  } else {
    List out = bygroup(data, stratum);
    has_stratum = 1;
    stratumn = out["index"];
    u_stratum = DataFrame(out["lookup"]);
  }
  
  bool has_time = hasVariable(data, time);
  if (!has_time) stop("data must contain the time variable");
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each subject");
  }
  
  bool has_time2 = hasVariable(data, time2);
  NumericVector time2n(n);
  if (has_time2) {
    NumericVector time2nz = data[time2];
    time2n = clone(time2nz);
    if (is_true(any(time2n <= timen))) {
      stop("time2 must be greater than time for each observation");
    }
  }
  
  // unify right censored data with counting process data
  NumericVector tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  bool has_event = hasVariable(data, event);
  if (!has_event) stop("data must contain the event variable");
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }
  
  bool has_weight = hasVariable(data, weight);
  NumericVector weightn(n, 1.0);
  if (has_weight) {
    NumericVector weightnz = data[weight];
    weightn = clone(weightnz);
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }
  
  std::string ct = conftype;
  std::for_each(ct.begin(), ct.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (!(ct=="none" || ct=="plain" || ct=="log" || ct=="log-log" ||
      ct=="logit" || ct=="arcsin")) {
    stop("conftype must be none, plain, log, log-log, logit, or arcsin");
  }
  
  if (conflev <= 0 || conflev >= 1) {
    stop("conflev must lie between 0 and 1");
  }
  
  // confidence interval for survival probability
  double z = R::qnorm((1+conflev)/2,0,1,1,0);
  
  // sort the data by rep
  if (has_rep) {
    IntegerVector order = seq(0, n-1);
    std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
      return repn[i] < repn[j];
    });
    
    repn = repn[order];
    stratumn = stratumn[order];
    tstartn = tstartn[order];
    tstopn = tstopn[order];
    eventn = eventn[order];
    weightn = weightn[order];
  }
  
  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (int i=1; i<n; ++i) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);
  
  IntegerVector rep8(n);
  IntegerVector stratum8(n), size8(n);
  NumericVector time8(n), nrisk8(n), nevent8(n), ncensor8(n);
  NumericVector surv8(n), sesurv8(n);
  NumericVector lower8(n), upper8(n);
  
  int index = 0;
  for (int h = 0; h < nreps; ++h) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());
    
    IntegerVector stratum1 = stratumn[q1];
    NumericVector tstart1 = tstartn[q1];
    NumericVector tstop1 = tstopn[q1];
    IntegerVector event1 = eventn[q1];
    NumericVector weight1 = weightn[q1];
    
    // sort by stopping time in descending order within each stratum
    IntegerVector order = seq(0, n1-1);
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      if (stratum1[i] != stratum1[j]) return stratum1[i] > stratum1[j];
      if (tstop1[i] != tstop1[j]) return tstop1[i] > tstop1[j];
      return event1[i] < event1[j];
    });
    
    stratum1 = stratum1[order];
    tstart1 = tstart1[order];
    tstop1 = tstop1[order];
    event1 = event1[order];
    weight1 = weight1[order];
    
    // sort by starting time in descending order within each stratum
    IntegerVector order1 = seq(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      if (stratum1[i] != stratum1[j]) return stratum1[i] > stratum1[j];
      return tstart1[i] > tstart1[j];
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
    
    int istratum = stratum1[0]; // current stratum
    int i1 = 0;                 // index for removing out-of-risk subjects

    int size = 0;
    IntegerVector stratum0(n1), size1(n1);
    NumericVector time0(n1), nrisk0(n1), nriskw0(n1), nriskw20(n1);
    NumericVector nevent0(n1), neventw0(n1), ncensor0(n1);
    
    int index1 = 0;
    for (int i = 0; i < n1; ) {
      // Reset when entering a new stratum
      if (stratum1[i] != istratum) {
        istratum = stratum1[i];
        i1 = i;
        size = 0; 
        nrisk = nriskw = nriskw2 = 0;
      }
      
      const double dtime = tstop1[i];
      
      // Process all persons tied at this dtime
      while (i < n1 && tstop1[i] == dtime) {
        ++nrisk;
        nriskw += weight1[i];
        nriskw2 += weight1[i]*weight1[i];
        
        if (event1[i] == 1) {
          ++nevent;
          neventw += weight1[i];
        } else {
          ++ncensor;
        }
        
        if (tstart1[i] == 0) ++size; // unique subjects in the stratum
        
        ++i;
        
        if (i < n1 && stratum1[i] != istratum) {
          size1[istratum] = size; // update size of the stratum
          break; 
        }
      }
      
      // remove subjects no longer at risk
      for (; i1<n1; ++i1) {
        const int p1 = order1[i1];
        if (tstart1[p1] < dtime || stratum1[p1] != istratum) break;
        
        nrisk--;
        nriskw -= weight1[p1];
        nriskw2 -= weight1[p1]*weight1[p1];
      }
      
      if (nevent > 0 || keep_censor) {
        stratum0[index1] = istratum;
        time0[index1] = dtime;
        nrisk0[index1] = nrisk;
        nriskw0[index1] = nriskw;
        nriskw20[index1] = nriskw2;
        nevent0[index1] = nevent;
        neventw0[index1] = neventw;
        ncensor0[index1] = ncensor;
        
        nevent = neventw = 0;
        ++index1;
      }
      
      ncensor = 0;
    }
    
    size1[istratum] = size; // update size of the last stratum
    
    IntegerVector size0(index1);
    NumericVector surv0(index1), sesurv0(index1);
    NumericVector lower0(index1), upper0(index1);
    
    istratum = stratum0[index1-1];
    for (int i = index1-1; i >= 0; i--) {
      double nevent = nevent0[i];
      double neventw = neventw0[i];
      double nriskw = nriskw0[i];
      double nriskw2 = nriskw20[i];
      
      if (stratum0[i] != istratum) { // hit a new stratum
        // reset temporary variables
        istratum = stratum0[i];
        surv = 1; vcumhaz = 0; sesurv = 0;
      }
      
      if (nevent > 0) {
        double p = neventw/nriskw;
        surv *= (1.0 - p);
        double m = nriskw*nriskw/nriskw2;
        vcumhaz += p/(m*(1.0 - p));
        sesurv = surv*sqrt(vcumhaz);
      }
      
      size0[i] = size1[istratum];
      surv0[i] = surv;
      sesurv0[i] = sesurv;
      if (ct != "none") {
        NumericVector ci = fsurvci(surv, sesurv, ct, z);
        lower0[i] = ci[0];
        upper0[i] = ci[1];
      }
    }
    
    // update the global index
    for (int i = 0; i < index1; ++i) {
      int j = index + i;
      int k = index1 - 1 - i;
      rep8[j] = h;
      stratum8[j] = stratum0[k];
      size8[j] = size0[k];
      time8[j] = time0[k];
      nrisk8[j] = nrisk0[k];
      nevent8[j] = nevent0[k];
      ncensor8[j] = ncensor0[k];
      surv8[j] = surv0[k];
      sesurv8[j] = sesurv0[k];
      if (ct != "none") {
        lower8[j] = lower0[k];
        upper8[j] = upper0[k];
      }
    }
    
    index += index1;
  } // end of rep loop
  
  // subset to the used portion
  IntegerVector sub = Range(0, index-1);
  rep8 = rep8[sub];
  stratum8 = stratum8[sub];
  size8 = size8[sub];
  time8 = time8[sub];
  nrisk8 = nrisk8[sub];
  nevent8 = nevent8[sub];
  ncensor8 = ncensor8[sub];
  surv8 = surv8[sub];
  sesurv8 = sesurv8[sub];
  if (ct != "none") {
    lower8 = lower8[sub];
    upper8 = upper8[sub];
  }
  
  List result = List::create(
    Named("size") = size8,
    Named("time") = time8,
    Named("nrisk") = nrisk8,
    Named("nevent") = nevent8,
    Named("ncensor") = ncensor8,
    Named("surv") = surv8,
    Named("sesurv") = sesurv8
  );
  
  if (ct != "none") {
    result.push_back(lower8, "lower");
    result.push_back(upper8, "upper");
    result.push_back(conflev, "conflev");
    result.push_back(ct, "conftype");
  }
  
  // add stratum and rep variables
  if (has_stratum) {
    for (int i = 0; i < p_stratum; ++i) {
      std::string s = as<std::string>(stratum[i]);
      SEXP col = u_stratum[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == INTSXP) {
        IntegerVector v = col;
        result.push_back(v[stratum8], s);
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        result.push_back(v[stratum8], s);
      } else if (col_type == STRSXP) {
        StringVector v = col;
        result.push_back(v[stratum8], s);
      } else {
        stop("Unsupported type for stratum variable" + s);
      }
    }
  }
  
  if (has_rep) {
    for (int i = 0; i < p_rep; ++i) {
      std::string s = as<std::string>(rep[i]);
      SEXP col = u_rep[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == INTSXP) {
        IntegerVector v = col;
        result.push_back(v[rep8], s);
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        result.push_back(v[rep8], s);
      } else if (col_type == STRSXP) {
        StringVector v = col;
        result.push_back(v[rep8], s);
      } else {
        stop("Unsupported type for rep variable" + s);
      }
    }
  }
  
  return result;
}


//' @title Estimate of Milestone Survival Difference
//' @description Obtains the estimate of milestone survival difference
//' between two treatment groups.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{rep}: The replication for by-group processing.
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{treat}: The treatment.
//'
//'   * \code{time}: The possibly right-censored survival time.
//'
//'   * \code{event}: The event indicator.
//'
//' @param rep The name of the replication variable in the input data.
//' @param stratum The name of the stratum variable in the input data.
//' @param treat The name of the treatment variable in the input data.
//' @param time The name of the time variable in the input data.
//' @param event The name of the event variable in the input data.
//' @param milestone The milestone time at which to calculate the
//'   survival probability.
//' @param survDiffH0 The difference in milestone survival probabilities
//'   under the null hypothesis. Defaults to 0 for superiority test.
//' @param conflev The level of the two-sided confidence interval for
//'   the difference in milestone survival probabilities. Defaults to 0.95.
//'
//' @return A data frame with the following variables:
//'
//' * \code{rep}: The replication.
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
//' df <- kmdiff(data = rawdata, rep = "iterationNumber",
//'              stratum = "stratum", treat = "treatmentGroup",
//'              time = "timeUnderObservation", event = "event",
//'              milestone = 12)
//' head(df)
//'
//' @export
// [[Rcpp::export]]
DataFrame kmdiff(const DataFrame data,
                 const StringVector& rep = "",
                 const StringVector& stratum = "",
                 const std::string treat = "treat",
                 const std::string time = "time",
                 const std::string event = "event",
                 const double milestone = NA_REAL,
                 const double survDiffH0 = 0,
                 const double conflev = 0.95) {
  int k, n = data.nrows();
  
  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(0);
  } else {
    List out = bygroup(data, rep);
    has_rep = 1;
    repn = out["index"];
    u_rep = DataFrame(out["lookup"]);
  }
  
  bool has_stratum;
  IntegerVector stratumn(n);
  DataFrame u_stratum;
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(0);
  } else {
    List out = bygroup(data, stratum);
    has_stratum = 1;
    stratumn = out["index"];
    u_stratum = DataFrame(out["lookup"]);
  }
  
  
  bool has_treat = hasVariable(data, treat);
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);
  
  if (!has_treat) {
    stop("data must contain the treat variable");
  }
  
  if (!has_time) {
    stop("data must contain the time variable");
  }
  
  if (!has_event) {
    stop("data must contain the event variable");
  }
  
  // create the numeric treat variable
  IntegerVector treatn(n);
  IntegerVector treatwi;
  NumericVector treatwn;
  StringVector treatwc;
  SEXP col = data[treat];
  SEXPTYPE col_type = TYPEOF(col);
  if (col_type == LGLSXP || col_type == INTSXP) {
    IntegerVector treatv = col;
    treatwi = unique(treatv);
    if (treatwi.size() != 2)
      stop("treat must have two and only two distinct values");
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwi == 0) | (treatwi == 1)))) {
      treatwi = IntegerVector::create(1, 0);
      treatn = 2 - treatv;
    } else {
      treatwi.sort();
      treatn = match(treatv, treatwi);
    }
  } else if (col_type == REALSXP) {
    NumericVector treatv = col;
    treatwn = unique(treatv);
    if (treatwn.size() != 2)
      stop("treat must have two and only two distinct values");
    if (is_true(all((treatwn == 0) | (treatwn == 1)))) {
      treatwn = NumericVector::create(1, 0);
      treatn = 2 - as<IntegerVector>(treatv);
    } else {
      treatwn.sort();
      treatn = match(treatv, treatwn);
    }
  } else if (col_type == STRSXP) {
    StringVector treatv = col;
    treatwc = unique(treatv);
    if (treatwc.size() != 2)
      stop("treat must have two and only two distinct values");
    treatwc.sort();
    treatn = match(treatv, treatwc);
  } else {
    stop("incorrect type for the treat variable in the input data");
  }
  
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each subject");
  }
  
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }
  
  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }
  
  if (milestone <= 0) {
    stop("milestone must be positive");
  }
  
  if (survDiffH0 <= -1 || survDiffH0 >= 1) {
    stop("survDiffH0 must lie between -1 and 1");
  }
  
  if (conflev <= 0 || conflev >= 1) {
    stop("conflev must lie between 0 and 1");
  }
  
  // sort the data by rep
  if (has_rep) {
    IntegerVector order = seq(0, n-1);
    std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
      return repn[i] < repn[j];
    });
    
    repn = repn[order];
    stratumn = stratumn[order];
    treatn = treatn[order];
    timen = timen[order];
    eventn = eventn[order];
  }
  
  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (int i=1; i<n; ++i) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);
  
  IntegerVector rep0(nreps);
  NumericVector surv10(nreps), surv20(nreps), survDiff0(nreps);
  NumericVector vsurv10(nreps), vsurv20(nreps), sesurvDiff0(nreps);
  NumericVector survDiffZ0(nreps), survDiffPValue0(nreps);
  NumericVector lower0(nreps), upper0(nreps);
  
  double z = R::qnorm((1+conflev)/2,0,1,1,0);
  
  bool noerr = true;
  int index = 0;
  for (int h=0; h<nreps; ++h) {
    bool skip = false;
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());
    
    IntegerVector stratum1 = stratumn[q1];
    IntegerVector treat1 = treatn[q1];
    NumericVector time1 = timen[q1];
    IntegerVector event1 = eventn[q1];
    
    // sort by stratum in descending order
    IntegerVector order1 = seq(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      return stratum1[i] < stratum1[j];
    });
    
    stratum1 = stratum1[order1];
    treat1 = treat1[order1];
    time1 = time1[order1];
    event1 = event1[order1];
    
    // identify the locations of the unique values of stratum
    IntegerVector idx1(1,0);
    for (int i=1; i<n1; ++i) {
      if (stratum1[i] != stratum1[i-1]) {
        idx1.push_back(i);
      }
    }
    
    int nstrata = static_cast<int>(idx1.size());
    idx1.push_back(n1);
    
    // whether the milestone exceeds the largest observed time
    for (int i=0; i<nstrata; ++i) {
      IntegerVector q2 = Range(idx1[i], idx1[i+1]-1);
      IntegerVector treat2 = treat1[q2];
      NumericVector time2 = time1[q2];
      NumericVector time21 = time2[treat2==1];
      NumericVector time22 = time2[treat2==2];
      
      if (milestone > std::min(max(time21), max(time22))) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else {
          for (int j=0; j<p_rep; ++j) {
            std::string s = as<std::string>(rep[j]);
            SEXP col = u_rep[s];
            SEXPTYPE col_type = TYPEOF(col);
            if (col_type == INTSXP) {
              IntegerVector repwi = col;
              reperr += " " + s + " = " + std::to_string(repwi[h]);
            } else if (col_type == REALSXP) {
              NumericVector repwn = col;
              reperr += " " + s + " = " + std::to_string(repwn[h]);
            } else if (col_type == STRSXP) {
              StringVector repwc = col;
              reperr += " " + s + " = " + repwc[h];
            } else {
              stop("unsupported type in rep variable: " + s);
            }
          }
        }
        
        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else {
          for (int j=0; j<p_stratum; ++j) {
            std::string s = as<std::string>(stratum[j]);
            SEXP col = u_stratum[s];
            SEXPTYPE col_type = TYPEOF(col);
            if (col_type == INTSXP) {
              IntegerVector stratumwi = col;
              stratumerr += " " + s + " = " + std::to_string(stratumwi[i]);
            } else if (col_type == REALSXP) {
              NumericVector stratumwn = col;
              stratumerr += " " + s + " = " + std::to_string(stratumwn[i]);
            } else if (col_type == STRSXP) {
              StringVector stratumwc = col;
              stratumerr += " " + s + " = " + stratumwc[i];
            } else {
              stop("unsupported type in stratum variable: " + s);
            }
          }
        }

        int k = milestone > max(time21) ? 0 : 1;
        std::string treaterr;
        SEXP col = data[treat];
        SEXPTYPE col_type = TYPEOF(col);
        if (col_type == LGLSXP || col_type == INTSXP) {
          IntegerVector treatwi = col;
          treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
        } else if (col_type == REALSXP) {
          NumericVector treatwn = col;
          treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
        } else if (col_type == STRSXP) {
          StringVector treatwc = col;
          treaterr = " " + treat + " = " + treatwc[k];
        } else {
          stop("unsupported type for treat variable: " + treat);
        }

        std::string str1 = "The milestone is larger than";
        std::string str2 = "the largest observed time for";
        std::string errmsg = str1 + " " + str2 + treaterr;
        if (!reperr.empty() || !stratumerr.empty()) {
          errmsg = errmsg + ":" + reperr + stratumerr;
        }
        
        if (noerr) {
          Rcout << errmsg << "\n";
          Rcout << "Additional warning messages are suppressed" << "\n";
          noerr = false;
        }
        
        skip = true;
        break;
      }
    }
    
    // skip the replication if there is a stratum with max time < milestone
    if (skip) continue;
    
    DataFrame dfin = DataFrame::create(
      _["stratum"] = stratum1,
      _["treat"] = treat1,
      _["time"] = time1,
      _["event"] = event1);
    
    DataFrame dfout = kmest(dfin, "stratum", "treat", "time", "", "event",
                            "", "none", 0.95, 0);
    
    IntegerVector stratum2 = dfout["stratum"];
    IntegerVector treat2 = dfout["treat"];
    IntegerVector treatsize = dfout["size"];
    NumericVector time2 = dfout["time"];
    NumericVector survival2 = dfout["surv"];
    NumericVector stderr2 = dfout["sesurv"];
    
    int n2 = static_cast<int>(stratum2.size());
    
    // identify the locations of the unique values of stratum
    IntegerVector idx2(1,0);
    for (int i=1; i<n2; ++i) {
      if (stratum2[i] != stratum2[i-1]) {
        idx2.push_back(i);
      }
    }
    
    idx2.push_back(n2);
    
    IntegerVector m(nstrata, 0); // number of subjects in each stratum
    for (int i=0; i<nstrata; ++i) {
      int j1 = idx2[i], j2 = idx2[i+1] - 1;
      if (treat2[j1] != 1 || treat2[j2] != 2) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else {
          for (int j=0; j<p_rep; ++j) {
            std::string s = as<std::string>(rep[j]);
            SEXP col = u_rep[s];
            SEXPTYPE col_type = TYPEOF(col);
            if (col_type == INTSXP) {
              IntegerVector repwi = col;
              reperr += " " + s + " = " + std::to_string(repwi[h]);
            } else if (col_type == REALSXP) {
              NumericVector repwn = col;
              reperr += " " + s + " = " + std::to_string(repwn[h]);
            } else if (col_type == STRSXP) {
              StringVector repwc = col;
              reperr += " " + s + " = " + repwc[h];
            } else {
              stop("unsupported type in rep variable: " + s);
            }
          }
        }
        
        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else {
          for (int j=0; j<p_stratum; ++j) {
            std::string s = as<std::string>(stratum[j]);
            SEXP col = u_stratum[s];
            SEXPTYPE col_type = TYPEOF(col);
            if (col_type == INTSXP) {
              IntegerVector stratumwi = col;
              stratumerr += " " + s + " = " + std::to_string(stratumwi[i]);
            } else if (col_type == REALSXP) {
              NumericVector stratumwn = col;
              stratumerr += " " + s + " = " + std::to_string(stratumwn[i]);
            } else if (col_type == STRSXP) {
              StringVector stratumwc = col;
              stratumerr += " " + s + " = " + stratumwc[i];
            } else {
              stop("unsupported type in stratum variable: " + s);
            }
          }
        }
        
        int k = treat2[j1] != 1 ? 0 : 1;
        std::string treaterr;
        SEXP col = data[treat];
        SEXPTYPE col_type = TYPEOF(col);
        if (col_type == LGLSXP || col_type == INTSXP) {
          IntegerVector treatwi = col;
          treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
        } else if (col_type == REALSXP) {
          NumericVector treatwn = col;
          treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
        } else if (col_type == STRSXP) {
          StringVector treatwc = col;
          treaterr = " " + treat + " = " + treatwc[k];
        } else {
          stop("unsupported type for treat variable: " + treat);
        }
        
        std::string str1 = "The data set does not contain";
        std::string errmsg = str1 + treaterr;
        if (!reperr.empty() || !stratumerr.empty()) {
          errmsg = errmsg + ":" + reperr + stratumerr;
        }
        
        if (noerr) {
          Rcout << errmsg << "\n";
          Rcout << "Additional warning messages are suppressed" << "\n";
          noerr = false;
        }
        
        skip = true;
        break;
      }
      
      m[i] += treatsize[j1] + treatsize[j2];
    }
    
    // skip the replication if there is a stratum without both treatments
    if (skip) continue;
    
    double M = sum(m);
    NumericVector p(nstrata);
    
    double surv1 = 0.0, surv2 = 0.0, vsurv1 = 0.0, vsurv2 = 0.0;
    for (int i=0; i<nstrata; ++i) {
      p[i] = m[i]/M; // fraction of subjects in the stratum
      IntegerVector q = Range(idx2[i], idx2[i+1]-1);
      IntegerVector treatx = treat2[q];
      NumericVector timex = time2[q];
      NumericVector survivalx = survival2[q];
      NumericVector stderrx = stderr2[q];
      
      NumericVector surv(2), vsurv(2);
      for (int j=0; j<2; ++j) {
        LogicalVector sub = (treatx == j+1);
        NumericVector time0 = timex[sub];
        NumericVector survival0 = survivalx[sub];
        NumericVector stderr0 = stderrx[sub];
        int K = sum(sub);
        
        // find the latest event time before milestone for each treat
        for (k = 0; k < K; ++k) {
          if (time0[k] > milestone) break;
        }
        
        if (k == 0) {
          surv[j] = 1;
          vsurv[j] = 0;
        } else {
          k--;
          surv[j] = survival0[k];
          vsurv[j] = stderr0[k]*stderr0[k];
        }
      }
      
      surv1 += p[i]*surv[0];
      surv2 += p[i]*surv[1];
      vsurv1 += p[i]*p[i]*vsurv[0];
      vsurv2 += p[i]*p[i]*vsurv[1];
    }
    
    rep0[index] = h;
    surv10[index] = surv1;
    surv20[index] = surv2;
    vsurv10[index] = vsurv1;
    vsurv20[index] = vsurv2;
    survDiff0[index] = surv1 - surv2;
    sesurvDiff0[index] = sqrt(vsurv1 + vsurv2);
    survDiffZ0[index] = (survDiff0[index] - survDiffH0)/sesurvDiff0[index];
    survDiffPValue0[index] = 2*R::pnorm(-fabs(survDiffZ0[index]),0,1,1,0);
    lower0[index] = survDiff0[index] - z*sesurvDiff0[index];
    upper0[index] = survDiff0[index] + z*sesurvDiff0[index];
    
    ++index;
  }
  
  // only keep nonmissing records
  IntegerVector sub = Range(0, index-1);
  rep0 = rep0[sub];
  surv10 = surv10[sub];
  surv20 = surv20[sub];
  survDiff0 = survDiff0[sub];
  vsurv10 = vsurv10[sub];
  vsurv20 = vsurv20[sub];
  sesurvDiff0 = sesurvDiff0[sub];
  survDiffZ0 = survDiffZ0[sub];
  survDiffPValue0 = survDiffPValue0[sub];
  lower0 = lower0[sub];
  upper0 = upper0[sub];
  
  List result = List::create(
    _["milestone"] = milestone,
    _["survDiffH0"] = survDiffH0,
    _["surv1"] = surv10,
    _["surv2"] = surv20,
    _["survDiff"] = survDiff0,
    _["vsurv1"] = vsurv10,
    _["vsurv2"] = vsurv20,
    _["sesurvDiff"] = sesurvDiff0,
    _["survDiffZ"] = survDiffZ0,
    _["survDiffPValue"] = survDiffPValue0,
    _["lower"] = lower0,
    _["upper"] = upper0,
    _["conflev"] = conflev);
  
  if (has_rep) {
    for (int i = 0; i < p_rep; ++i) {
      std::string s = as<std::string>(rep[i]);
      SEXP col = u_rep[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == INTSXP) {
        IntegerVector v = col;
        result.push_back(v[rep0], s);
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        result.push_back(v[rep0], s);
      } else if (col_type == STRSXP) {
        StringVector v = col;
        result.push_back(v[rep0], s);
      } else {
        stop("Unsupported type for rep variable" + s);
      }
    }
  }

  return result;
}


//' @title Log-Rank Test of Survival Curve Difference
//' @description Obtains the log-rank test using the Fleming-Harrington
//' family of weights.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{rep}: The replication for by-group processing.
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
//' @param rep The name(s) of the replication variable(s) in the input data.
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
//' * \code{rep}: The replication.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' df <- lrtest(data = rawdata, rep = "iterationNumber",
//'              stratum = "stratum", treat = "treatmentGroup",
//'              time = "timeUnderObservation", event = "event",
//'              rho1 = 0.5, rho2 = 0)
//' head(df)
//'
//' @export
// [[Rcpp::export]]
DataFrame lrtest(const DataFrame data,
                 const StringVector& rep = "",
                 const StringVector& stratum = "",
                 const std::string treat = "treat",
                 const std::string time = "time",
                 const std::string time2 = "",
                 const std::string event = "event",
                 const std::string weight = "",
                 const bool weight_readj = false,
                 const double rho1 = 0,
                 const double rho2 = 0) {
  int n = data.nrows();
  
  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(0);
  } else {
    List out = bygroup(data, rep);
    has_rep = 1;
    repn = out["index"];
    u_rep = DataFrame(out["lookup"]);
  }
  
  IntegerVector stratumn(n);
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    stratumn.fill(0);
  } else {
    List out = bygroup(data, stratum);
    stratumn = out["index"];
  }

  bool has_treat = hasVariable(data, treat);
  if (!has_treat) stop("data must contain the treat variable");
  
  // create the numeric treat variable
  IntegerVector treatn(n);
  IntegerVector treatwi;
  NumericVector treatwn;
  StringVector treatwc;
  SEXP col = data[treat];
  SEXPTYPE col_type = TYPEOF(col);
  if (col_type == LGLSXP || col_type == INTSXP) {
    IntegerVector treatv = col;
    treatwi = unique(treatv);
    if (treatwi.size() != 2)
      stop("treat must have two and only two distinct values");
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwi == 0) | (treatwi == 1)))) {
      treatwi = IntegerVector::create(1, 0);
      treatn = 2 - treatv;
    } else {
      treatwi.sort();
      treatn = match(treatv, treatwi);
    }
  } else if (col_type == REALSXP) {
    NumericVector treatv = col;
    treatwn = unique(treatv);
    if (treatwn.size() != 2)
      stop("treat must have two and only two distinct values");
    if (is_true(all((treatwn == 0) | (treatwn == 1)))) {
      treatwn = NumericVector::create(1, 0);
      treatn = 2 - as<IntegerVector>(treatv);
    } else {
      treatwn.sort();
      treatn = match(treatv, treatwn);
    }
  } else if (col_type == STRSXP) {
    StringVector treatv = col;
    treatwc = unique(treatv);
    if (treatwc.size() != 2)
      stop("treat must have two and only two distinct values");
    treatwc.sort();
    treatn = match(treatv, treatwc);
  } else {
    stop("incorrect type for the treat variable in the input data");
  }
  
  bool has_time = hasVariable(data, time);
  if (!has_time) stop("data must contain the time variable");
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each subject");
  }
  
  bool has_time2 = hasVariable(data, time2);
  NumericVector time2n(n);
  if (has_time2) {
    NumericVector time2nz = data[time2];
    time2n = clone(time2nz);
    if (is_true(any(time2n <= timen))) {
      stop("time2 must be greater than time for each observation");
    }
  }
  
  // unify right censored data with counting process data
  NumericVector tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  bool has_event = hasVariable(data, event);
  if (!has_event) stop("data must contain the event variable");
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }
  
  bool has_weight = hasVariable(data, weight);
  NumericVector weightn(n, 1.0);
  if (has_weight) {
    NumericVector weightnz = data[weight];
    weightn = clone(weightnz);
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }
  
  if (rho1 < 0) stop("rho1 must be non-negative");
  if (rho2 < 0) stop("rho2 must be non-negative");
  
  // sort the data by rep
  if (has_rep) {
    IntegerVector order = seq(0, n-1);
    std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
      return repn[i] < repn[j];
    });
    
    repn = repn[order];
    stratumn = stratumn[order];
    treatn = treatn[order];
    tstartn = tstartn[order];
    tstopn = tstopn[order];
    eventn = eventn[order];
    weightn = weightn[order];
  }
  
  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (int i=1; i<n; ++i) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);
  
  IntegerVector rep0(nreps);
  NumericVector uscore0(nreps), vscore0(nreps);
  NumericVector logRankZ0(nreps), logRankPValue0(nreps);

  for (int h=0; h<nreps; ++h) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());
    
    IntegerVector stratum1 = stratumn[q1];
    IntegerVector treat1 = treatn[q1];
    NumericVector tstart1 = tstartn[q1];
    NumericVector tstop1 = tstopn[q1];
    IntegerVector event1 = eventn[q1];
    NumericVector weight1 = weightn[q1];
    
    // sort by stopping time in descending order within each stratum
    IntegerVector order = seq(0, n1-1);
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
      if (tstop1[i] != tstop1[j]) return tstop1[i] > tstop1[j];
      return event1[i] < event1[j];
    });
    
    stratum1 = stratum1[order];
    treat1 = treat1[order];
    tstart1 = tstart1[order];
    tstop1 = tstop1[order];
    event1 = event1[order];
    weight1 = weight1[order];
    
    // sort by starting time in descending order within each stratum
    IntegerVector order1 = seq(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
      return tstart1[i] > tstart1[j];
    });
    
    double u = 0; // score
    double v = 0; // information (variance for score)
    
    double nrisk = 0, nrisk_1 = 0, nrisk_2 = 0;    // number at risk
    double nriskw = 0, nriskw_1 = 0, nriskw_2 = 0; // weighted # at risk
    double nriskw2_1 = 0, nriskw2_2 = 0;           // weight^2 # at risk
    double nevent = 0;                             // number of events
    double neventw = 0, neventw_1 = 0, neventw_2;  // weighted # events
    double neventwa, neventwa_1, neventwa_2;  // adjusted weighted # events
    
    int istratum = stratum1[0]; // current stratum
    int i1 = 0;                 // index for removing out-of-risk subjects

    if (rho1 == 0 && rho2 == 0) {
      for (int i = 0; i < n1; ) {
        // Reset when entering a new stratum
        if (stratum1[i] != istratum) {
          istratum = stratum1[i];
          i1 = i;
          
          nrisk = nrisk_1 = nrisk_2 = 0;
          nriskw = nriskw_1 = nriskw_2 = 0; 
          nriskw2_1 = nriskw2_2 = 0;
        }
        
        const double dtime = tstop1[i];
        
        // Process all persons tied at this dtime
        for (; i < n1 && tstop1[i] == dtime && stratum1[i] == istratum; ++i) {
          ++nrisk;
          nriskw += weight1[i];
          if (treat1[i] == 1) {
            ++nrisk_1;
            nriskw_1 += weight1[i];
            nriskw2_1 += weight1[i]*weight1[i];
          } else {
            ++nrisk_2;
            nriskw_2 += weight1[i];
            nriskw2_2 += weight1[i]*weight1[i];
          }
          
          if (event1[i] == 1) {
            ++nevent;
            neventw += weight1[i];
            if (treat1[i] == 1) neventw_1 += weight1[i];
          }
        }
        
        // Remove subjects leaving risk set
        for (; i1<n1; ++i1) {
          const int p1 = order1[i1];
          if (tstart1[p1] < dtime || stratum1[p1] != istratum) break;
          
          nrisk--;
          nriskw -= weight1[p1];
          if (treat1[p1] == 1) {
            nrisk_1--;
            nriskw_1 -= weight1[p1];
            nriskw2_1 -= weight1[p1]*weight1[p1];
          } else {
            nrisk_2--;
            nriskw_2 -= weight1[p1];
            nriskw2_2 -= weight1[p1]*weight1[p1];
          }
        }
        
        // Add contributions for deaths at this time
        if (nevent > 0) {
          if (nrisk_1 > 0 && nrisk_2 > 0) {
            if (!weight_readj) {
              u += neventw_1 - neventw*(nriskw_1/nriskw);
              double v1 = pow(nriskw_1/nriskw, 2)*nriskw2_2;
              double v2 = pow(nriskw_2/nriskw, 2)*nriskw2_1;
              v += nevent*(nrisk - nevent)/(nrisk*(nrisk - 1))*(v1 + v2);
            } else {
              neventw_2 = neventw - neventw_1;
              neventwa_1 = neventw_1*nrisk_1/nriskw_1;
              neventwa_2 = neventw_2*nrisk_2/nriskw_2;
              neventwa = neventwa_1 + neventwa_2;
              u += neventwa_1 - neventwa*(nrisk_1/nrisk);
              double v1 = pow(nrisk_1/nrisk, 2)*nriskw2_2*
                pow(nrisk_2/nriskw_2, 2);
              double v2 = pow(nrisk_2/nrisk, 2)*nriskw2_1*
                pow(nrisk_1/nriskw_1, 2);
              v += nevent*(nrisk - nevent)/(nrisk*(nrisk - 1))*(v1 + v2);
            }
          }
          
          // Reset after processing deaths
          nevent = neventw = neventw_1 = 0;
        }
      }
    } else {
      IntegerVector stratum0(n1);
      NumericVector time0(n1);
      NumericVector nrisk0(n1), nrisk_10(n1), nrisk_20(n1);
      NumericVector nriskw0(n1), nriskw_10(n1), nriskw_20(n1);
      NumericVector nriskw2_10(n1), nriskw2_20(n1);
      NumericVector nevent0(n1);
      NumericVector neventw0(n1), neventw_10(n1);
      
      int index1 = 0;
      for (int i = 0; i < n1; ) {
        // Reset when entering a new stratum
        if (stratum1[i] != istratum) {
          istratum = stratum1[i];
          i1 = i;
          
          nrisk = nrisk_1 = nrisk_2 = 0;
          nriskw = nriskw_1 = nriskw_2 = 0;
          nriskw2_1 = nriskw2_2 = 0;
        }
        
        const double dtime = tstop1[i];
        
        // Process all persons tied at this dtime
        for (; i < n1 && tstop1[i] == dtime && stratum1[i] == istratum; ++i) {
          ++nrisk;
          nriskw += weight1[i];
          if (treat1[i] == 1) {
            ++nrisk_1;
            nriskw_1 += weight1[i];
            nriskw2_1 += weight1[i]*weight1[i];
          } else {
            ++nrisk_2;
            nriskw_2 += weight1[i];
            nriskw2_2 += weight1[i]*weight1[i];
          }
          
          if (event1[i] == 1) {
            ++nevent;
            neventw += weight1[i];
            if (treat1[i] == 1) neventw_1 += weight1[i];
          }
        }
        
        // Remove subjects leaving risk set
        for (; i1<n1; ++i1) {
          const int p1 = order1[i1];
          if (tstart1[p1] < dtime || stratum1[p1] != istratum) break;
          
          nrisk--;
          nriskw -= weight1[p1];
          if (treat1[p1] == 1) {
            nrisk_1--;
            nriskw_1 -= weight1[p1];
            nriskw2_1 -= weight1[p1]*weight1[p1];
          } else {
            nrisk_2--;
            nriskw_2 -= weight1[p1];
            nriskw2_2 -= weight1[p1]*weight1[p1];
          }
        }
        
        // add to the temporary storage
        if (nevent > 0) {
          stratum0[index1] = istratum;
          time0[index1] = dtime;
          nrisk0[index1] = nrisk;
          nrisk_10[index1] = nrisk_1;
          nrisk_20[index1] = nrisk_2;
          nriskw0[index1] = nriskw;
          nriskw_10[index1] = nriskw_1;
          nriskw_20[index1] = nriskw_2;
          nriskw2_10[index1] = nriskw2_1;
          nriskw2_20[index1] = nriskw2_2;
          nevent0[index1] = nevent;
          neventw0[index1] = neventw;
          neventw_10[index1] = neventw_1;
          
          ++index1;
          
          // Reset after processing deaths
          nevent = neventw = neventw_1 = 0;
        }
      }
      
      double surv = 1.0;
      istratum = stratum0[index1 - 1];
      for (int i = index1-1; i >= 0; i--) {
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
        
        // add to the main terms with S(t-)^rho1*(1 - S(t-))^rho2 weights
        double w = pow(surv, rho1)*pow(1 - surv, rho2);
        
        if (weight_readj) {
          neventw_2 = neventw - neventw_1;
          neventwa_1 = neventw_1*nrisk_1/nriskw_1;
          neventwa_2 = neventw_2*nrisk_2/nriskw_2;
          neventwa = neventwa_1 + neventwa_2;
        }
        
        if (nrisk_1 > 0 && nrisk_2 > 0) {
          if (!weight_readj) {
            u += w*(neventw_1 - neventw*(nriskw_1/nriskw));
            double v1 = pow(nriskw_1/nriskw, 2)*nriskw2_2;
            double v2 = pow(nriskw_2/nriskw, 2)*nriskw2_1;
            v += w*w*nevent*(nrisk - nevent)/(nrisk*(nrisk - 1))*(v1 + v2);
          } else {
            u += w*(neventwa_1 - neventwa*(nrisk_1/nrisk));
            double v1 = pow(nrisk_1/nrisk, 2)*nriskw2_2*
              pow(nrisk_2/nriskw_2, 2);
            double v2 = pow(nrisk_2/nrisk, 2)*nriskw2_1*
              pow(nrisk_1/nriskw_1, 2);
            v += w*w*nevent*(nrisk - nevent)/(nrisk*(nrisk - 1))*(v1 + v2);
          }
        }
        
        // update survival probability
        if (!weight_readj) {
          surv *= (1.0 - neventw/nriskw);
        } else {
          surv *= (1.0 - neventwa/nrisk);
        }
      }
    }

    double z = u/sqrt(v);
    double p = 2*R::pnorm(-fabs(z), 0.0, 1.0, 1, 0);
    
    rep0[h] = h;
    uscore0[h] = u;
    vscore0[h] = v;
    logRankZ0[h] = z;
    logRankPValue0[h] = p;
  }
  
  List result = List::create(
    Named("uscore") = uscore0,
    Named("vscore") = vscore0,
    Named("logRankZ") = logRankZ0,
    Named("logRankPValue") = logRankPValue0,
    Named("weight_readj") = weight_readj,
    Named("rho1") = rho1,
    Named("rho2") = rho2
  );
  
  // add back rep variables
  if (has_rep) {
    for (int i = 0; i < p_rep; ++i) {
      std::string s = as<std::string>(rep[i]);
      SEXP col = u_rep[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == INTSXP) {
        IntegerVector v = col;
        result.push_back(v[rep0], s);
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        result.push_back(v[rep0], s);
      } else if (col_type == STRSXP) {
        StringVector v = col;
        result.push_back(v[rep0], s);
      } else {
        stop("Unsupported type for rep variable" + s);
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
//'   * \code{rep}: The replication for by-group processing.
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{time}: The possibly right-censored survival time.
//'
//'   * \code{event}: The event indicator.
//'
//' @param rep The name of the replication variable in the input data.
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
//' * \code{rep}: The replication.
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
DataFrame rmest(const DataFrame data,
                const StringVector& rep = "",
                const StringVector& stratum = "",
                const std::string time = "time",
                const std::string event = "event",
                const double milestone = NA_REAL,
                const double conflev = 0.95,
                const bool biascorrection = false) {
  int n = data.nrows();
  
  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(0);
  } else {
    List out = bygroup(data, rep);
    has_rep = 1;
    repn = out["index"];
    u_rep = DataFrame(out["lookup"]);
  }
  
  bool has_stratum;
  IntegerVector stratumn(n);
  DataFrame u_stratum;
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(0);
  } else {
    List out = bygroup(data, stratum);
    has_stratum = 1;
    stratumn = out["index"];
    u_stratum = DataFrame(out["lookup"]);
  }
  
  
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);
  
  if (!has_time) {
    stop("data must contain the time variable");
  }
  
  if (!has_event) {
    stop("data must contain the event variable");
  }
  
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each subject");
  }
  
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }
  
  
  if (std::isnan(milestone)) {
    stop("milestone must be provided");
  }
  
  if (milestone <= 0) {
    stop("milestone must be positive");
  }
  
  if (conflev <= 0 || conflev >= 1) {
    stop("conflev must lie between 0 and 1");
  }
  
  // sort the data by rep
  if (has_rep) {
    IntegerVector order = seq(0, n-1);
    std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
      return repn[i] < repn[j];
    });
    
    repn = repn[order];
    stratumn = stratumn[order];
    timen = timen[order];
    eventn = eventn[order];
  }
  
  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (int i=1; i<n; ++i) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);
  
  IntegerVector rep0(n), stratum0(n), size0(n);
  NumericVector rmst0(n), stderr0(n), lower0(n), upper0(n);
  
  double z = R::qnorm((1+conflev)/2,0,1,1,0);
  
  bool noerr = true;
  int index = 0;
  for (int h=0; h<nreps; ++h) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());
    
    IntegerVector stratum1 = stratumn[q1];
    NumericVector time1 = timen[q1];
    IntegerVector event1 = eventn[q1];
    
    // sort by stratum, time, and event with event in descending order
    IntegerVector order1 = seq(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
      if (time1[i] != time1[j]) return time1[i] < time1[j];
      return event1[i] > event1[j];
    });
    
    stratum1 = stratum1[order1];
    time1 = time1[order1];
    event1 = event1[order1];
    
    // identify the locations of the unique values of stratum
    IntegerVector idx1(1,0);
    for (int i=1; i<n1; ++i) {
      if (stratum1[i] != stratum1[i-1]) {
        idx1.push_back(i);
      }
    }
    
    int nstrata = static_cast<int>(idx1.size());
    idx1.push_back(n1);
    
    for (int i=0; i<nstrata; ++i) {
      IntegerVector q2 = Range(idx1[i], idx1[i+1]-1);
      NumericVector time2 = time1[q2];
      IntegerVector event2 = event1[q2];
      int n2 = static_cast<int>(q2.size());
      
      if (milestone > max(time2)) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else {
          for (int j=0; j<p_rep; ++j) {
            std::string s = as<std::string>(rep[j]);
            SEXP col = u_rep[s];
            SEXPTYPE col_type = TYPEOF(col);
            if (col_type == INTSXP) {
              IntegerVector repwi = col;
              reperr += " " + s + " = " + std::to_string(repwi[h]);
            } else if (col_type == REALSXP) {
              NumericVector repwn = col;
              reperr += " " + s + " = " + std::to_string(repwn[h]);
            } else if (col_type == STRSXP) {
              StringVector repwc = col;
              reperr += " " + s + " = " + repwc[h];
            } else {
              stop("unsupported type in rep variable: " + s);
            }
          }
        }
        
        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else {
          for (int j=0; j<p_stratum; ++j) {
            std::string s = as<std::string>(stratum[j]);
            SEXP col = u_stratum[s];
            SEXPTYPE col_type = TYPEOF(col);
            if (col_type == INTSXP) {
              IntegerVector stratumwi = col;
              stratumerr += " " + s + " = " + std::to_string(stratumwi[i]);
            } else if (col_type == REALSXP) {
              NumericVector stratumwn = col;
              stratumerr += " " + s + " = " + std::to_string(stratumwn[i]);
            } else if (col_type == STRSXP) {
              StringVector stratumwc = col;
              stratumerr += " " + s + " = " + stratumwc[i];
            } else {
              stop("unsupported type in stratum variable: " + s);
            }
          }
        }
        
        std::string str1 = "The milestone is larger than";
        std::string str2 = "the largest observed time";
        std::string errmsg = str1 + " " + str2;
        if (!reperr.empty() || !stratumerr.empty()) {
          errmsg = errmsg + ":" + reperr + stratumerr;
        }
        
        if (noerr) {
          Rcout << errmsg << "\n";
          Rcout << "Additional warning messages are suppressed" << "\n";
          noerr = false;
        }
        
        continue;
      }
      
      NumericVector time0(n2), nrisk0(n2), nevent0(n2), surv0(n2);
      int index1 = 0;
      double t = 0, nrisk = n2, nevent = 0, surv = 1.0;
      bool cache = false;
      for (int j=0; j<n2; ++j) {
        if ((j == 0 && event2[j] == 1) ||
            (j >= 1 && event2[j] == 1 && time2[j] > time2[j-1])) {
          // new event
          // add the info for the previous event
          if (cache) {
            surv *= (1.0 - nevent/nrisk);
            
            time0[index1] = t;
            nrisk0[index1] = nrisk;
            nevent0[index1] = nevent;
            surv0[index1] = surv;
            
            ++index1;
          }
          
          // update the buffer for the current event time
          t = time2[j];
          nrisk = n2-j;
          nevent = 1;
          
          cache = true;
        } else if (j >= 1 && event2[j] == 1 && event2[j-1] == 1 &&
          time2[j] == time2[j-1]) { // tied event
          ++nevent;
        } else if (j >= 1 && event2[j] == 0 && event2[j-1] == 1) {
          // new censoring
          // add the info for the previous event
          surv *= (1.0 - nevent/nrisk);
          
          time0[index1] = t;
          nrisk0[index1] = nrisk;
          nevent0[index1] = nevent;
          surv0[index1] = surv;
          
          ++index1;
          
          // empty the cache for the current event time
          cache = false;
        }
      }
      
      // add the info for the last event
      if (cache) {
        surv *= (1.0 - nevent/nrisk);
        
        time0[index1] = t;
        nrisk0[index1] = nrisk;
        nevent0[index1] = nevent;
        surv0[index1] = surv;
        
        ++index1;
      }
      
      // only keep nonmissing records
      int N;
      if (index1 == 0) { // no event
        N = 0;
        time0 = NumericVector::create(0.0);
        nrisk0 = NumericVector::create(n2);
        nevent0 = NumericVector::create(0.0);
        surv0 = NumericVector::create(1.0);
      } else { // at least 1 event
        IntegerVector sub = seq(0, index1 - 1);
        time0 = time0[sub];
        nrisk0 = nrisk0[sub];
        nevent0 = nevent0[sub];
        surv0 = surv0[sub];
        
        // locate the latest event time before milestone
        NumericVector milestone1(1, milestone);
        N = findInterval3(milestone1, time0, 0, 0, 0)[0];
        
        // prepend time zero information
        time0.push_front(0.0);
        nrisk0.push_front(n2);
        nevent0.push_front(0.0);
        surv0.push_front(1.0);
      }
      
      // replace the last time of interest with milestone
      if (N == time0.size() - 1) {
        time0.push_back(milestone);
      } else {
        time0[N+1] = milestone;
      }
      
      // calculate the partial sum of the trapezoid integration
      NumericVector rmstx(N+1);
      rmstx[0] = surv0[0]*(time0[1] - time0[0]);
      for (int k=1; k<=N; ++k) {
        rmstx[k] = rmstx[k-1] + surv0[k]*(time0[k+1] - time0[k]);
      }
      
      // calculate rmst and its variance
      double u = rmstx[N];
      double v = 0.0;
      for (int k=1; k<=N; ++k) {
        // rmst from the kth event time to milestone
        double a = u - rmstx[k-1];
        // do not add variance if the largest observed time is an event time
        if (nrisk0[k] > nevent0[k]) {
          v += nevent0[k]*a*a/(nrisk0[k]*(nrisk0[k] - nevent0[k]));
        }
      }
      
      // apply bias correction if requested
      if (biascorrection) {
        double m1 = 0;
        for (int k=1; k<=N; ++k) {
          m1 += nevent0[k];
        }
        
        if (m1 <= 1.0) {
          std::string reperr;
          if (!has_rep) {
            reperr = "";
          } else {
            for (int j=0; j<p_rep; ++j) {
              std::string s = as<std::string>(rep[j]);
              SEXP col = u_rep[s];
              SEXPTYPE col_type = TYPEOF(col);
              if (col_type == INTSXP) {
                IntegerVector repwi = col;
                reperr += " " + s + " = " + std::to_string(repwi[h]);
              } else if (col_type == REALSXP) {
                NumericVector repwn = col;
                reperr += " " + s + " = " + std::to_string(repwn[h]);
              } else if (col_type == STRSXP) {
                StringVector repwc = col;
                reperr += " " + s + " = " + repwc[h];
              } else {
                stop("unsupported type in rep variable: " + s);
              }
            }
          }
          
          std::string stratumerr;
          if (!has_stratum) {
            stratumerr = "";
          } else {
            for (int j=0; j<p_stratum; ++j) {
              std::string s = as<std::string>(stratum[j]);
              SEXP col = u_stratum[s];
              SEXPTYPE col_type = TYPEOF(col);
              if (col_type == INTSXP) {
                IntegerVector stratumwi = col;
                stratumerr += " " + s + " = " + std::to_string(stratumwi[i]);
              } else if (col_type == REALSXP) {
                NumericVector stratumwn = col;
                stratumerr += " " + s + " = " + std::to_string(stratumwn[i]);
              } else if (col_type == STRSXP) {
                StringVector stratumwc = col;
                stratumerr += " " + s + " = " + stratumwc[i];
              } else {
                stop("unsupported type in stratum variable: " + s);
              }
            }
          }
          
          std::string str1 = "Bias correction is not done due to no or";
          std::string str2 = "only 1 event before the milestone time:";
          std::string errmsg = str1 + " " + str2;
          if (!reperr.empty() || !stratumerr.empty()) {
            errmsg = errmsg + ":" + reperr + stratumerr;
          }
          
          if (noerr) {
            Rcout << errmsg << "\n";
            Rcout << "Additional warning messages are suppressed" << "\n";
            noerr = false;
          }
        } else {
          v = m1/(m1 - 1.0)*v;
        }
      }
      
      rep0[index] = h;
      stratum0[index] = i;
      size0[index] = n2;
      rmst0[index] = u;
      stderr0[index] = sqrt(v);
      lower0[index] = u - z*stderr0[index];
      upper0[index] = u + z*stderr0[index];
      
      ++index;
    }
  }
  
  // only keep nonmissing records
  IntegerVector sub = Range(0, index - 1);
  if (index == 0) {
    stop("no replication enables valid inference");
  }
  
  rep0 = rep0[sub];
  stratum0 = stratum0[sub];
  size0 = size0[sub];
  rmst0 = rmst0[sub];
  stderr0 = stderr0[sub];
  lower0 = lower0[sub];
  upper0 = upper0[sub];
  
  List result = List::create(
    _["size"] = size0,
    _["milestone"] = milestone,
    _["rmst"] = rmst0,
    _["stderr"] = stderr0,
    _["lower"] = lower0,
    _["upper"] = upper0,
    _["conflev"] = conflev,
    _["biascorrection"] = biascorrection);
  
  // add stratum and rep variables
  if (has_stratum) {
    for (int i = 0; i < p_stratum; ++i) {
      std::string s = as<std::string>(stratum[i]);
      SEXP col = u_stratum[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == INTSXP) {
        IntegerVector v = col;
        result.push_back(v[stratum0], s);
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        result.push_back(v[stratum0], s);
      } else if (col_type == STRSXP) {
        StringVector v = col;
        result.push_back(v[stratum0], s);
      } else {
        stop("Unsupported type for stratum variable" + s);
      }
    }
  }

  if (has_rep) {
    for (int i = 0; i < p_rep; ++i) {
      std::string s = as<std::string>(rep[i]);
      SEXP col = u_rep[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == INTSXP) {
        IntegerVector v = col;
        result.push_back(v[rep0], s);
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        result.push_back(v[rep0], s);
      } else if (col_type == STRSXP) {
        StringVector v = col;
        result.push_back(v[rep0], s);
      } else {
        stop("Unsupported type for rep variable" + s);
      }
    }
  }

  return result;
}


//' @title Estimate of Restricted Mean Survival Time Difference
//' @description Obtains the estimate of restricted mean survival time
//' difference between two treatment groups.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{rep}: The replication for by-group processing.
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{treat}: The treatment.
//'
//'   * \code{time}: The possibly right-censored survival time.
//'
//'   * \code{event}: The event indicator.
//'
//' @param rep The name of the replication variable in the input data.
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
//' * \code{rep}: The replication number.
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
//' df <- rmdiff(data = rawdata, rep = "iterationNumber",
//'              stratum = "stratum", treat = "treatmentGroup",
//'              time = "timeUnderObservation", event = "event",
//'              milestone = 12)
//' head(df)
//'
//' @export
// [[Rcpp::export]]
DataFrame rmdiff(const DataFrame data,
                 const StringVector& rep = "",
                 const StringVector& stratum = "",
                 const std::string treat = "treat",
                 const std::string time = "time",
                 const std::string event = "event",
                 const double milestone = NA_REAL,
                 const double rmstDiffH0 = 0,
                 const double conflev = 0.95,
                 const bool biascorrection = false) {
  int n = data.nrows();
  
  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(0);
  } else {
    List out = bygroup(data, rep);
    has_rep = 1;
    repn = out["index"];
    u_rep = DataFrame(out["lookup"]);
  }
  
  bool has_stratum;
  IntegerVector stratumn(n);
  DataFrame u_stratum;
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(0);
  } else {
    List out = bygroup(data, stratum);
    has_stratum = 1;
    stratumn = out["index"];
    u_stratum = DataFrame(out["lookup"]);
  }
  
  bool has_treat = hasVariable(data, treat);
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);
  if (!has_treat) stop("data must contain the treat variable");
  if (!has_time) stop("data must contain the time variable");
  if (!has_event) stop("data must contain the event variable");

  // create the numeric treat variable
  IntegerVector treatn(n);
  IntegerVector treatwi;
  NumericVector treatwn;
  StringVector treatwc;
  SEXP col = data[treat];
  SEXPTYPE col_type = TYPEOF(col);
  if (col_type == LGLSXP || col_type == INTSXP) {
    IntegerVector treatv = col;
    treatwi = unique(treatv);
    if (treatwi.size() != 2)
      stop("treat must have two and only two distinct values");
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwi == 0) | (treatwi == 1)))) {
      treatwi = IntegerVector::create(1, 0);
      treatn = 2 - treatv;
    } else {
      treatwi.sort();
      treatn = match(treatv, treatwi);
    }
  } else if (col_type == REALSXP) {
    NumericVector treatv = col;
    treatwn = unique(treatv);
    if (treatwn.size() != 2)
      stop("treat must have two and only two distinct values");
    if (is_true(all((treatwn == 0) | (treatwn == 1)))) {
      treatwn = NumericVector::create(1, 0);
      treatn = 2 - as<IntegerVector>(treatv);
    } else {
      treatwn.sort();
      treatn = match(treatv, treatwn);
    }
  } else if (col_type == STRSXP) {
    StringVector treatv = col;
    treatwc = unique(treatv);
    if (treatwc.size() != 2)
      stop("treat must have two and only two distinct values");
    treatwc.sort();
    treatn = match(treatv, treatwc);
  } else {
    stop("incorrect type for the treat variable in the input data");
  }
  
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each subject");
  }
  
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }
  
  
  if (std::isnan(milestone)) stop("milestone must be provided");
  if (milestone <= 0) stop("milestone must be positive");
  
  if (conflev <= 0 || conflev >= 1) {
    stop("conflev must lie between 0 and 1");
  }
  
  // sort the data by rep
  if (has_rep) {
    IntegerVector order = seq(0, n-1);
    std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
      return repn[i] < repn[j];
    });
    
    repn = repn[order];
    stratumn = stratumn[order];
    treatn = treatn[order];
    timen = timen[order];
    eventn = eventn[order];
  }
  
  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (int i=1; i<n; ++i) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);
  
  IntegerVector rep0(nreps);
  NumericVector rmst10(nreps), rmst20(nreps), rmstDiff0(nreps);
  NumericVector vrmst10(nreps), vrmst20(nreps), sermstDiff0(nreps);
  NumericVector rmstDiffZ0(nreps), rmstDiffPValue0(nreps);
  NumericVector lower0(nreps), upper0(nreps);
  
  double z = R::qnorm((1.0 + conflev)/2.0,0,1,1,0);
  
  bool noerr = true;
  int index = 0;
  for (int h=0; h<nreps; ++h) {
    bool skip = false;
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    
    IntegerVector stratum1 = stratumn[q1];
    IntegerVector treat1 = treatn[q1];
    NumericVector time1 = timen[q1];
    IntegerVector event1 = eventn[q1];
    
    DataFrame dfin = DataFrame::create(
      _["stratum"] = stratum1,
      _["treat"] = treat1,
      _["time"] = time1,
      _["event"] = event1);
    
    DataFrame dfout = rmest(dfin, "stratum", "treat", "time", "event",
                            milestone, 0.95, biascorrection);
    
    IntegerVector stratum2 = dfout["stratum"];
    IntegerVector treat2 = dfout["treat"];
    IntegerVector treatsize = dfout["size"];
    NumericVector rmstime2 = dfout["rmst"];
    NumericVector stderr2 = dfout["stderr"];
    int n2 = static_cast<int>(stratum2.size());
    
    // identify the locations of the unique values of stratum
    IntegerVector idx2(1,0);
    for (int i=1; i<n2; ++i) {
      if (stratum2[i] != stratum2[i-1]) {
        idx2.push_back(i);
      }
    }
    
    int nstrata = static_cast<int>(idx2.size());
    idx2.push_back(n2);
    
    IntegerVector m(nstrata, 0); // number of subjects in each stratum
    for (int i=0; i<nstrata; ++i) {
      int j1 = idx2[i], j2 = idx2[i+1] - 1;
      if (treat2[j1] != 1 || treat2[j2] != 2) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else {
          for (int j=0; j<p_rep; ++j) {
            std::string s = as<std::string>(rep[j]);
            SEXP col = u_rep[s];
            SEXPTYPE col_type = TYPEOF(col);
            if (col_type == INTSXP) {
              IntegerVector repwi = col;
              reperr += " " + s + " = " + std::to_string(repwi[h]);
            } else if (col_type == REALSXP) {
              NumericVector repwn = col;
              reperr += " " + s + " = " + std::to_string(repwn[h]);
            } else if (col_type == STRSXP) {
              StringVector repwc = col;
              reperr += " " + s + " = " + repwc[h];
            } else {
              stop("unsupported type in rep variable: " + s);
            }
          }
        }
        
        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else {
          for (int j=0; j<p_stratum; ++j) {
            std::string s = as<std::string>(stratum[j]);
            SEXP col = u_stratum[s];
            SEXPTYPE col_type = TYPEOF(col);
            if (col_type == INTSXP) {
              IntegerVector stratumwi = col;
              stratumerr += " " + s + " = " + std::to_string(stratumwi[i]);
            } else if (col_type == REALSXP) {
              NumericVector stratumwn = col;
              stratumerr += " " + s + " = " + std::to_string(stratumwn[i]);
            } else if (col_type == STRSXP) {
              StringVector stratumwc = col;
              stratumerr += " " + s + " = " + stratumwc[i];
            } else {
              stop("unsupported type in stratum variable: " + s);
            }
          }
        }
        
        int k = treat2[j1] != 1 ? 0 : 1;
        std::string treaterr;
        SEXP col = data[treat];
        SEXPTYPE col_type = TYPEOF(col);
        if (col_type == LGLSXP || col_type == INTSXP) {
          IntegerVector treatwi = col;
          treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
        } else if (col_type == REALSXP) {
          NumericVector treatwn = col;
          treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
        } else if (col_type == STRSXP) {
          StringVector treatwc = col;
          treaterr = " " + treat + " = " + treatwc[k];
        } else {
          stop("unsupported type for treat variable: " + treat);
        }
        
        std::string str1 = "The data set does not contain";
        std::string errmsg = str1 + treaterr;
        if (!reperr.empty() || !stratumerr.empty()) {
          errmsg = errmsg + ":" + reperr + stratumerr;
        }
        
        if (noerr) {
          Rcout << errmsg << "\n";
          Rcout << "Additional warning messages are suppressed" << "\n";
          noerr = false;
        }
        
        skip = true;
        break;
      }
      
      m[i] += treatsize[j1] + treatsize[j2];
    }
    
    // skip the replication if there is a stratum without both treatments
    if (skip) continue;
    
    double M = sum(m);
    NumericVector p(nstrata);
    
    double rmst1 = 0.0, rmst2 = 0.0, vrmst1 = 0.0, vrmst2 = 0.0;
    for (int i=0; i<nstrata; ++i) {
      p[i] = m[i]/M; // fraction of subjects in the stratum
      IntegerVector q = Range(idx2[i], idx2[i+1]-1);
      NumericVector rmst = rmstime2[q];
      NumericVector stderrx = stderr2[q];
      NumericVector vrmst = stderrx*stderrx;
      
      rmst1 += p[i]*rmst[0];
      rmst2 += p[i]*rmst[1];
      vrmst1 += p[i]*p[i]*vrmst[0];
      vrmst2 += p[i]*p[i]*vrmst[1];
    }
    
    rep0[index] = h;
    rmst10[index] = rmst1;
    rmst20[index] = rmst2;
    vrmst10[index] = vrmst1;
    vrmst20[index] = vrmst2;
    rmstDiff0[index] = rmst1 - rmst2;
    sermstDiff0[index] = sqrt(vrmst1 + vrmst2);
    rmstDiffZ0[index] = (rmstDiff0[index] - rmstDiffH0)/sermstDiff0[index];
    rmstDiffPValue0[index] = 2*R::pnorm(-fabs(rmstDiffZ0[index]),0,1,1,0);
    lower0[index] = rmstDiff0[index] - z*sermstDiff0[index];
    upper0[index] = rmstDiff0[index] + z*sermstDiff0[index];
    
    ++index;
  }
  
  if (index == 0) stop("no replication enables valid inference");
  
  IntegerVector sub = Range(0, index-1);
  rep0 = rep0[sub];
  rmst10 = rmst10[sub];
  rmst20 = rmst20[sub];
  rmstDiff0 = rmstDiff0[sub];
  vrmst10 = vrmst10[sub];
  vrmst20 = vrmst20[sub];
  sermstDiff0 = sermstDiff0[sub];
  rmstDiffZ0 = rmstDiffZ0[sub];
  rmstDiffPValue0 = rmstDiffPValue0[sub];
  lower0 = lower0[sub];
  upper0 = upper0[sub];
  
  List result = List::create(
    _["milestone"] = milestone,
    _["rmstDiffH0"] = rmstDiffH0,
    _["rmst1"] = rmst10,
    _["rmst2"] = rmst20,
    _["rmstDiff"] = rmstDiff0,
    _["vrmst1"] = vrmst10,
    _["vrmst2"] = vrmst20,
    _["sermstDiff"] = sermstDiff0,
    _["rmstDiffZ"] = rmstDiffZ0,
    _["rmstDiffPValue"] = rmstDiffPValue0,
    _["lower"] = lower0,
    _["upper"] = upper0,
    _["conflev"] = conflev,
    _["biascorrection"] = biascorrection);
  
  if (has_rep) {
    for (int i = 0; i < p_rep; ++i) {
      std::string s = as<std::string>(rep[i]);
      SEXP col = u_rep[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == INTSXP) {
        IntegerVector v = col;
        result.push_back(v[rep0], s);
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        result.push_back(v[rep0], s);
      } else if (col_type == STRSXP) {
        StringVector v = col;
        result.push_back(v[rep0], s);
      } else {
        stop("Unsupported type for rep variable" + s);
      }
    }
  }

  return result;
}


// all-in-one function for log-likelihood, score, and information matrix
// for the AFT model
List f_der_1(int p, const NumericVector& par, void* ex) {
  aftparams *param = (aftparams *) ex;
  const int n = param->z.nrow();
  const int nvar = param->z.ncol();
  const int dist_code = param->dist_code;
  
  NumericVector eta(n);
  for (int person = 0; person < n; ++person) {
    double val = param->offset[person];
    for (int i = 0; i < nvar; ++i) {
      val += par[i] * param->z(person,i);
    }
    eta[person] = val;
  }
  
  NumericVector sig(n, 1.0);
  if (dist_code != 1) { // all except exponential
    for (int person = 0; person < n; ++person) {
      int k = param->strata[person] + nvar;
      sig[person] = exp(par[k]);
    }
  }
  
  // Main loop 
  double loglik = 0.0;
  NumericVector score(p);
  NumericMatrix imat(p,p);
  for (int person = 0; person < n; ++person) {
    int k = param->strata[person] + nvar;
    double wt = param->weight[person];
    double sigma = sig[person];
    double logsig = log(sigma);
    double eta_p = eta[person];
    double tstart_p = param->tstart[person];
    double tstop_p = param->tstop[person];
    NumericVector z = param->z(person, _)/sigma;
    
    double u, u1, u2, c0, c1, c2, d, d1, d2;
    double q, q1, q2, w1, w2, num, den, term;
    
    switch (param->status[person]) {
    
    case 1: // event
      switch (dist_code) {
      
      case 1: case 2: // exponential / weibull
        u = (log(tstop_p) - eta_p)/sigma;
        loglik += wt * (u - exp(u) - logsig);
        
        c1 = -wt * (1 - exp(u));
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        if (dist_code == 2)
          score[k] += wt * ((1 - exp(u)) * (-u) - 1);
        
        c1 = wt*exp(u);
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        if (dist_code==2) { // weibull
          c2 = wt*(exp(u)*u - (1 - exp(u)));
          for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
          imat(k,k) += c2*u;
        } 
        
        break;
        
      case 3: case 4: // lognormal / normal
        u = (dist_code == 3) ? (log(tstop_p) - eta_p) / sigma
        : (tstop_p - eta_p) / sigma;
        loglik += wt * (R::dnorm(u,0,1,1) - logsig);
        
        c1 = wt * u;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        score[k] += wt * (u * u - 1);
        
        c2 = wt*2*u;
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += wt*z[i]*z[j];
        for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
        imat(k,k) += c2*u;
        
        break;
        
      case 5: case 6: // loglogistic / logistic
        u = (dist_code == 5) ? (log(tstop_p) - eta_p) / sigma
        : (tstop_p - eta_p) / sigma;
        loglik += wt * (R::dlogis(u,0,1,1) - logsig);
        
        c0 = 1 - 2 * R::plogis(u, 0, 1, 0, 0);
        c1 = wt * c0;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        score[k] += wt * (c0 * u - 1);
        
        c1 = wt*2*R::dlogis(u,0,1,0);
        c2 = wt*(2*R::dlogis(u,0,1,0)*u + 1 - 2*R::plogis(u,0,1,0,0));
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
        imat(k,k) += c2*u;
        
        break;
        
      }
      
      break;
      
      
    case 3: // interval censoring
      switch (dist_code) {
      
      case 1: case 2: // exponential / weibull
        u1 = (log(tstart_p) - eta_p)/sigma;
        u2 = (log(tstop_p) - eta_p)/sigma;
        loglik += wt * log(exp(-exp(u1)) - exp(-exp(u2)));
        
        num = exp(u1 - exp(u1)) - exp(u2 - exp(u2));
        den = exp(-exp(u1)) - exp(-exp(u2));
        c1 = wt * num / den;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        if (dist_code == 2) {
          term = (exp(u1 - exp(u1)) * u1 - exp(u2 - exp(u2)) * u2) / den;
          score[k] += wt * term;
        }
        
        w1 = exp(u1); w2 = exp(u2);
        q1 = exp(-w1); q2 = exp(-w2);
        d1 = w1*q1; d2 = w2*q2;
        c1 = wt*(pow((d1 - d2)/(q1 - q2), 2) + 
          (d1*(1-w1) - d2*(1-w2))/(q1 - q2));
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        if (dist_code==2) { // weibull
          c2 = wt*((d1 - d2)*(d1*u1 - d2*u2)/pow(q1 - q2, 2) + 
            (d1*(1 + (1-w1)*u1) - d2*(1 + (1-w2)*u2))/(q1 - q2));
          for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
          imat(k,k) += wt*(pow((d1*u1 - d2*u2)/(q1 - q2), 2) +
            (d1*(1 + (1-w1)*u1)*u1 - d2*(1 + (1-w2)*u2)*u2)/(q1 - q2));
        }
        
        break;
        
      case 3: case 4: // lognormal / normal
        u1 = (dist_code == 3) ? (log(tstart_p) - eta_p) / sigma
        : (tstart_p - eta_p) / sigma;
        u2 = (dist_code == 3)  ? (log(tstop_p) - eta_p) / sigma
        : (tstop_p - eta_p) / sigma;
        loglik += wt * log(R::pnorm(u1,0,1,0,0) - R::pnorm(u2,0,1,0,0));
        
        d1 = R::dnorm(u1, 0, 1, 0); d2 = R::dnorm(u2, 0, 1, 0);
        q1 = R::pnorm(u1, 0, 1, 0, 0); q2 = R::pnorm(u2, 0, 1, 0, 0);
        c1 = wt * (d1 - d2) / (q1 - q2);
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        score[k] += wt * (d1 * u1 - d2 * u2) / (q1 - q2);
        
        c1 = wt*(pow((d1 - d2)/(q1 -  q2), 2) + 
          (-d1*u1 + d2*u2)/(q1 - q2));
        c2 = wt*((d1 - d2)*(d1*u1 - d2*u2)/pow(q1 - q2, 2) + 
          (d1*(1 - u1*u1) - d2*(1 - u2*u2))/(q1 - q2));
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
        imat(k,k) += wt*(pow((d1*u1 - d2*u2)/(q1 - q2), 2) + 
          (d1*(1 - u1*u1)*u1 - d2*(1 - u2*u2)*u2)/(q1 - q2));
        break;
        
      case 5: case 6: // loglogistic / logistic
        u1 = (dist_code == 5) ? (log(tstart_p) - eta_p) / sigma
        : (tstart_p - eta_p) / sigma;
        u2 = (dist_code == 5) ? (log(tstop_p) - eta_p) / sigma
        : (tstop_p - eta_p) / sigma;
        loglik += wt * log(R::plogis(u1,0,1,0,0) - R::plogis(u2,0,1,0,0));
        
        d1 = R::dlogis(u1, 0, 1, 0); d2 = R::dlogis(u2, 0, 1, 0);
        q1 = R::plogis(u1, 0, 1, 0, 0); q2 = R::plogis(u2, 0, 1, 0, 0);
        c1 = wt * (d1 - d2) / (q1 - q2);
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        score[k] += wt * (d1 * u1 - d2 * u2) / (q1 - q2);
        
        c1 = wt*(pow((d1 - d2)/(q1 - q2), 2) + 
          (d1*(2*q1-1) - d2*(2*q2-1))/(q1 - q2));
        c2 = wt*((d1 - d2)*(d1*u1 - d2*u2)/pow(q1 - q2, 2) + 
          (d1*(1+(2*q1-1)*u1) - d2*(1+(2*q2-1)*u2))/(q1 - q2));
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
        imat(k,k) += wt*(pow((d1*u1 - d2*u2)/(q1 - q2), 2) + 
          (d1*(1+(2*q1-1)*u1)*u1 - d2*(1+(2*q2-1)*u2)*u2)/(q1 - q2));
        
        break;
        
      }
      
      break;
      
      
    case 2: // left censoring
      switch (dist_code) {
      
      case 1: case 2: // exponential / weibull
        u = (log(tstop_p) - eta_p)/sigma;
        loglik += wt * log(1.0 - exp(-exp(u)));
        
        num = -exp(u - exp(u));
        den = 1 - exp(-exp(u));
        c1 = wt * num / den;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        if (dist_code == 2) score[k] += c1 * u;
        
        w2 = exp(u); q2 = exp(-w2); d2 = w2*q2;
        c1 = wt*(pow(d2/(1 - q2), 2) - d2*(1-w2)/(1 - q2));
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        if (dist_code==2) {
          c2 = wt*(pow(d2/(1-q2),2)*u - d2*(1+(1-w2)*u)/(1-q2));
          for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
          imat(k,k) += c2*u;
        }
        
        break;
        
      case 3: case 4: // lognormal / normal
        u = (dist_code == 3) ? (log(tstop_p) - eta_p) / sigma
        : (tstop_p - eta_p) / sigma;
        loglik += wt * R::pnorm(u,0,1,1,1);
        
        d = R::dnorm(u, 0, 1, 0);
        q = R::pnorm(u, 0, 1, 1, 0);
        c1 = -wt * d / q;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        score[k] += c1 * u;
        
        c1 = wt*(pow(d/q, 2) + d*u/q);
        c2 = wt*(pow(d/q, 2)*u - d*(1-u*u)/q);
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
        imat(k,k) += c2*u;
        
        break;
        
      case 5: case 6: // loglogistic / logistic
        u = (dist_code == 5) ? (log(tstop_p) - eta_p) / sigma
        : (tstop_p - eta_p) / sigma;
        loglik += wt * R::plogis(u,0,1,1,1);
        
        q = R::plogis(u, 0, 1, 0, 0);
        c1 = -wt * q;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        score[k] += c1 * u;
        
        d2 = R::dlogis(u,0,1,0);
        c1 = wt*d2;
        c2 = wt*(d2*u - q);
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
        imat(k,k) += c2*u;
        
        break;
        
      }
      
      break;
      
      
    case 0: // right censoring
      switch (dist_code) {
      
      case 1: case 2: // exponential / weibull
        u = (log(tstart_p) - eta_p)/sigma;
        loglik += wt * (-exp(u));
        
        c1 = wt * exp(u);
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        if (dist_code == 2) score[k] += c1 * u;
        
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        if (dist_code==2) {
          c2 = wt*exp(u)*(1+u);
          for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
          imat(k,k) += c2*u;
        }
        
        break;
        
      case 3: case 4: // lognormal / normal
        u = (dist_code == 3) ? (log(tstart_p) - eta_p) / sigma
        : (tstart_p - eta_p) / sigma;
        loglik += wt * R::pnorm(u,0,1,0,1);
        
        d = R::dnorm(u, 0, 1, 0);
        q = R::pnorm(u, 0, 1, 0, 0);
        c1 = wt * d / q;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        score[k] += c1 * u;
        
        c1 = wt*(pow(d/q,2) - d*u/q);
        c2 = wt*(pow(d/q,2)*u + d*(1-u*u)/q);
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
        imat(k,k) += c2*u;
        
        break;
        
      case 5: case 6: // loglogistic / logistic
        u = (dist_code == 5) ? (log(tstart_p) - eta_p) / sigma
        : (tstart_p - eta_p) / sigma;
        loglik += wt * R::plogis(u,0,1,0,1);
        
        q = R::plogis(u, 0, 1, 1, 0);
        c1 = wt * q;
        for (int i = 0; i < nvar; ++i) score[i] += c1 * z[i];
        score[k] += c1 * u;
        
        d1 = R::dlogis(u,0,1,0); q1 = R::plogis(u,0,1,0,0);
        c1 = wt*d1;
        c2 = wt*(1-q1 + d1*u);
        for (int i=0; i<nvar; ++i)
          for (int j=0; j<=i; ++j)
            imat(i,j) += c1*z[i]*z[j];
        for (int j=0; j<nvar; ++j) imat(k,j) += c2*z[j];
        imat(k,k) += c2*u;
        
        break;
        
      }
      
      break;
      
    default:
      stop("Unknown status: " + std::to_string(param->status[person]));
    }
  }
  
  // fill upper triangle
  for (int i=0; i<p-1; ++i)
    for (int j=i+1; j<p; ++j)
      imat(i,j) = imat(j,i);
  
  return List::create(
    Named("loglik") = loglik,
    Named("score") = score,
    Named("imat") = imat
  );
}


// score residual matrix
NumericMatrix f_ressco_1(int p, const NumericVector& par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  
  // compute eta
  NumericVector eta(n);
  for (int person = 0; person < n; ++person) {
    double val = param->offset[person];
    for (int i = 0; i < nvar; ++i) val += par[i] * param->z(person, i);
    eta[person] = val;
  }
  
  // compute sigma
  NumericVector sig(n, 1.0);
  if (param->dist_code != 1) { // not exponential
    for (int person = 0; person < n; ++person) {
      int k = param->strata[person] + nvar;
      sig[person] = exp(par[k]);
    }
  }
  
  // Main loop to compute residuals
  NumericMatrix resid(n, p);
  for (int person = 0; person < n; ++person) {
    double sigma = sig[person];
    NumericVector z = param->z(person, _) / sigma;
    int k = param->strata[person] + nvar;
    double eta_p = eta[person];
    
    double u, u1, u2, c1, w1, w2, q1, q2, d1, d2;
    
    switch (param->status[person]) {
    
    case 1: // event
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u = (log(param->tstop[person]) - eta_p) / sigma;
        c1 = -(1 - exp(u));
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        if (param->dist_code == 2) resid(person, k) = (1 - exp(u)) * (-u) - 1;
        break;
      case 3: case 4: // lognormal / normal
        u = (param->dist_code == 3) ?
        (log(param->tstop[person]) - eta_p) / sigma :
        (param->tstop[person] - eta_p) / sigma;
        for (int i = 0; i < nvar; ++i) resid(person, i) = u * z[i];
        resid(person, k) = u * u - 1;
        break;
      case 5: case 6: // loglogistic / logistic
        u = (param->dist_code == 5) ?
        (log(param->tstop[person]) - eta_p) / sigma :
        (param->tstop[person] - eta_p) / sigma;
        c1 = 1 - 2 * R::plogis(u, 0, 1, 0, 0);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = c1 * u - 1;
        break;
      }
      break;
      
    case 3: // interval censoring
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u1 = (log(param->tstart[person]) - eta_p) / sigma;
        u2 = (log(param->tstop[person]) - eta_p) / sigma;
        w1 = exp(u1); w2 = exp(u2);
        q1 = exp(-w1); q2 = exp(-w2);
        d1 = w1 * q1; d2 = w2 * q2;
        c1 = (d1 - d2) / (q1 - q2);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        if (param->dist_code == 2)
          resid(person, k) = (d1 * u1 - d2 * u2) / (q1 - q2);
        break;
      case 3: case 4: // lognormal / normal
        u1 = (param->dist_code == 3) ?
        (log(param->tstart[person]) - eta_p) / sigma :
        (param->tstart[person] - eta_p) / sigma;
        u2 = (param->dist_code == 3) ?
        (log(param->tstop[person]) - eta_p) / sigma :
          (param->tstop[person] - eta_p) / sigma;
        d1 = R::dnorm(u1, 0, 1, 0); d2 = R::dnorm(u2, 0, 1, 0);
        q1 = R::pnorm(u1, 0, 1, 0, 0); q2 = R::pnorm(u2, 0, 1, 0, 0);
        c1 = (d1 - d2) / (q1 - q2);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = (d1 * u1 - d2 * u2) / (q1 - q2);
        break;
      case 5: case 6: // loglogistic / logistic
        u1 = (param->dist_code == 5) ?
        (log(param->tstart[person]) - eta_p) / sigma :
        (param->tstart[person] - eta_p) / sigma;
        u2 = (param->dist_code == 5) ?
        (log(param->tstop[person]) - eta_p) / sigma :
          (param->tstop[person] - eta_p) / sigma;
        d1 = R::dlogis(u1, 0, 1, 0); d2 = R::dlogis(u2, 0, 1, 0);
        q1 = R::plogis(u1, 0, 1, 0, 0); q2 = R::plogis(u2, 0, 1, 0, 0);
        c1 = (d1 - d2) / (q1 - q2);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = (d1 * u1 - d2 * u2) / (q1 - q2);
        break;
      }
      break;
      
    case 2: // left censoring
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u = (log(param->tstop[person]) - eta_p) / sigma;
        w2 = exp(u); q2 = exp(-w2); d2 = w2 * q2;
        c1 = -d2 / (1 - q2);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        if (param->dist_code == 2) resid(person, k) = c1 * u;
        break;
      case 3: case 4: // lognormal / normal
        u = (param->dist_code == 3) ?
        (log(param->tstop[person]) - eta_p) / sigma :
        (param->tstop[person] - eta_p) / sigma;
        c1 = -R::dnorm(u, 0, 1, 0) / R::pnorm(u, 0, 1, 1, 0);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = c1 * u;
        break;
      case 5: case 6: // loglogistic / logistic
        u = (param->dist_code == 5) ?
        (log(param->tstop[person]) - eta_p) / sigma :
        (param->tstop[person] - eta_p) / sigma;
        c1 = -R::plogis(u, 0, 1, 0, 0);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = c1 * u;
        break;
      }
      break;
      
    case 0: // right censoring
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u = (log(param->tstart[person]) - eta_p) / sigma;
        c1 = exp(u);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        if (param->dist_code == 2) resid(person, k) = c1 * u;
        break;
      case 3: case 4: // lognormal / normal
        u = (param->dist_code == 3) ?
        (log(param->tstart[person]) - eta_p) / sigma :
        (param->tstart[person] - eta_p) / sigma;
        c1 = R::dnorm(u, 0, 1, 0) / R::pnorm(u, 0, 1, 0, 0);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = c1 * u;
        break;
      case 5: case 6: // loglogistic / logistic
        u = (param->dist_code == 5) ?
        (log(param->tstart[person]) - eta_p) / sigma :
        (param->tstart[person] - eta_p) / sigma;
        c1 = R::plogis(u, 0, 1, 1, 0);
        for (int i = 0; i < nvar; ++i) resid(person, i) = c1 * z[i];
        resid(person, k) = c1 * u;
        break;
      }
      break;
      
    default:
      stop("Unknown status: " + std::to_string(param->status[person]));
    }
  }
  
  return resid;
}


// substitute information matrix guaranteed to be positive definite
NumericMatrix f_jj_1(int p, const NumericVector& par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();

  NumericMatrix resid = f_ressco_1(p, par, param);
  NumericMatrix jj(p,p);
  for (int person = 0; person < n; ++person) {
    double w = param->weight[person];
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        jj(i,j) += w * resid(person,i) * resid(person,j);
      }
    }
  }
  
  return jj;
}


// underlying optimization algorithm for lifereg
List liferegloop(int p, const NumericVector& par, void *ex,
                 int maxiter, double eps,
                 const IntegerVector& colfit, int ncolfit) {
  aftparams *param = (aftparams *) ex;
  
  int iter, halving = 0;
  bool fail = false;
  double toler = 1e-12;
  
  int nstrata = param->nstrata;
  int nvar = param->z.ncol();
  int nsub = param->z.nrow();
  
  NumericMatrix z1 = param->z;
  NumericVector mu(nvar), sigma(nvar);
  NumericMatrix z2(nsub, nvar);
  
  // --- standardize z once ---
  for (int i = 0; i < nvar; ++i) {
    NumericVector u = z1(_, i);
    double s = sd(u), m = mean(u);
    if (is_true(all((u == 0) | (u == 1)))) {
      m = 0; s = 1;
    }
    mu[i] = m;
    sigma[i] = s;
    for (int k = 0; k < nsub; ++k) z2(k, i) = (u[k] - m) / s;
  }
  
  // --- initial beta ---
  NumericVector beta(p), newbeta(p);
  beta[0] = par[0];
  for (int i = 1; i < nvar; ++i) {
    beta[i] = par[i] * sigma[i];
    beta[0] += par[i] * mu[i];
  }
  if (param->dist_code != 1)
    std::copy(par.begin() + nvar, par.end(), beta.begin() + nvar);
  
  aftparams para = {param->dist_code, param->strata, param->tstart, 
                    param->tstop, param->status, param->weight, 
                    param->offset, z2, nstrata};
  
  
  List der = f_der_1(p, beta, &para);
  double loglik = der["loglik"];
  double newlk = 0;
  NumericVector u = der["score"];
  NumericMatrix imat = as<NumericMatrix>(der["imat"]);
  NumericMatrix jj;
  NumericVector u1(ncolfit);
  NumericMatrix imat1(ncolfit, ncolfit);
  NumericMatrix jj1(ncolfit, ncolfit);

  for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
  
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      imat1(i, j) = imat(colfit[i], colfit[j]);
  
  // --- first step ---
  if (cholesky2(imat1, ncolfit, toler) < 0) {
    jj = f_jj_1(p, beta, &para);
    for (int i = 0; i < ncolfit; ++i)
      for (int j = 0; j < ncolfit; ++j)
        jj1(i, j) = jj(colfit[i], colfit[j]);
    cholesky2(jj1, ncolfit, toler);
    chsolve2(jj1, ncolfit, u1);
  } else {
    chsolve2(imat1, ncolfit, u1);
  }
  
  u.fill(0.0);
  for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
  newbeta = beta + u;
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_1(p, newbeta, &para);
    newlk = der["loglik"];
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && halving == 0 && fabs(1 - loglik / newlk) < eps) break;
    
    if (fail || newlk < loglik) {
      ++halving;
      for (int i = 0; i < p; ++i)
        newbeta[i] = 0.5 * (beta[i] + newbeta[i]);
      
      // special handling of sigmas
      if (halving == 1 && param->dist_code != 1) {
        for (int i = 0; i < nstrata; ++i)
          if (beta[nvar + i] - newbeta[nvar + i] > 1.1)
            newbeta[nvar + i] = beta[nvar + i] - 1.1;
      }
      continue;
    }
    
    // --- update ---
    halving = 0;
    beta = clone(newbeta);
    loglik = newlk;
    u = der["score"];
    imat = as<NumericMatrix>(der["imat"]);
    
    for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
    
    for (int i = 0; i < ncolfit; ++i)
      for (int j = 0; j < ncolfit; ++j)
        imat1(i, j) = imat(colfit[i], colfit[j]);
    
    if (cholesky2(imat1, ncolfit, toler) < 0) {
      jj = f_jj_1(p, beta, &para);
      for (int i = 0; i < ncolfit; ++i)
        for (int j = 0; j < ncolfit; ++j)
          jj1(i, j) = jj(colfit[i], colfit[j]);
      cholesky2(jj1, ncolfit, toler);
      chsolve2(jj1, ncolfit, u1);
    } else {
      chsolve2(imat1, ncolfit, u1);
    }
    
    u.fill(0.0);
    for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
    newbeta = beta + u;
  }
  
  if (iter == maxiter) fail = true;
  
  // --- rescale back ---
  for (int i = 1; i < nvar; ++i) {
    newbeta[i] /= sigma[i];
    newbeta[0] -= newbeta[i] * mu[i];
  }
  
  
  // rescale the information matrix accordingly
  imat = as<NumericMatrix>(der["imat"]);
  NumericMatrix jmat = clone(imat);
  for (int i = 0; i < nvar; ++i)
    for (int j = 0; j < nvar; ++j)
      imat(i,j) = jmat(0,0)*mu[i]*mu[j] + jmat(0,j)*mu[i]*sigma[j] +
        jmat(i,0)*mu[j]*sigma[i] + jmat(i,j)*sigma[i]*sigma[j];
  
  for (int i = nvar; i < p; ++i) {
    for (int j = 0; j < nvar; ++j) {
      imat(i,j) = jmat(i,0)*mu[j] + jmat(i,j)*sigma[j];
      imat(j,i) = imat(i,j);
    }
  }
  
  // compute variance matrix
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      imat1(i, j) = imat(colfit[i], colfit[j]);
  
  NumericMatrix var1 = invsympd(imat1, ncolfit, toler);
  NumericMatrix var(p, p);
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      var(colfit[i], colfit[j]) = var1(i, j);
  
  return List::create(
    Named("coef") = newbeta,
    Named("iter") = iter,
    Named("var") = var,
    Named("loglik") = newlk,
    Named("fail") = fail);
}


// confidence limit of profile likelihood method
double liferegplloop(int p, const NumericVector& par, void *ex,
                     int maxiter, double eps,
                     int k, int which, double l0) {
  aftparams *param = (aftparams *) ex;

  int iter;
  bool fail = false;
  double toler = 1e-12;
  
  NumericVector beta(p), newbeta(p);
  double loglik, newlk;
  NumericVector u(p), delta(p);
  NumericMatrix imat(p, p), jj(p, p), v(p, p);
  
  // --- first step ---
  beta = clone(par);
  
  List der = f_der_1(p, beta, param);
  loglik = der["loglik"];
  u = der["score"];
  imat = as<NumericMatrix>(der["imat"]);
  
  // compute inverse of inv (with cholesky fallback to jj)
  jj = clone(imat); // test cholesky on a copy (cholesky2 overwrites)
  
  if (cholesky2(jj, p, toler) < 0) {
    jj = f_jj_1(p, beta, param);
    v = invsympd(jj, p, toler); // inv of jj
  } else {
    v = invsympd(imat, p, toler); // inv of imat
  }

  // Lagrange multiplier method as used in SAS PROC LOGISTIC
  double w = 0;
  for (int i=0; i<p; ++i) 
    for (int j=0; j<p; ++j) 
      w -= u[i]*v(i,j)*u[j];
  
  double underroot = -2*(l0 - loglik + 0.5*w)/v(k,k);
  double lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
  u[k] += lambda;
  
  delta.fill(0.0);
  for (int i=0; i<p; ++i) 
    for (int j=0; j<p; ++j) 
      delta[i] += v(i,j)*u[j];

  // update beta
  newbeta = beta + delta;
  
  // --- main iteration ---
  for (iter=0; iter<maxiter; ++iter) {
    der = f_der_1(p, newbeta, param);
    newlk = der["loglik"];
    
    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && fabs(newlk - l0) < eps && w < eps) break;
    
    // update step
    beta = clone(newbeta);
    loglik = newlk;
    u = der["score"];
    imat = as<NumericMatrix>(der["imat"]);
    jj = clone(imat);
    
    if (cholesky2(jj, p, toler) < 0) {
      jj = f_jj_1(p, beta, param);
      v = invsympd(jj, p, toler);
    } else {
      v = invsympd(imat, p, toler);
    }
    
    // Lagrange multiplier method as used in SAS PROC LOGISTIC
    double w = 0;
    for (int i=0; i<p; ++i) 
      for (int j=0; j<p; ++j) 
        w -= u[i]*v(i,j)*u[j];
    
    double underroot = -2*(l0 - loglik + 0.5*w)/v(k,k);
    double lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
    u[k] += lambda;
    
    delta.fill(0.0);
    for (int i=0; i<p; ++i) 
      for (int j=0; j<p; ++j) 
        delta[i] += v(i,j)*u[j];
    
    // update beta
    newbeta = beta + delta;
  }
  
  if (iter == maxiter) fail = true;
  if (fail) warning("The algorithm in liferegplloop did not converge");

  return newbeta[k];
}



// first and second derivatives of log likelihood with respect to eta and tau
List f_ld_1(NumericVector eta, NumericVector sig, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->tstart.size();
  
  NumericVector g(n), dg(n), ddg(n), ds(n), dds(n), dsg(n);
  
  // Main loop to compute  derivatives
  for (int i = 0; i < n; ++i) {
    double sigma = sig[i];
    double u, u1, u2, d1, d2, q1, q2, w1, w2;
    
    switch (param->status[i]) {
    
    case 1: // event
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u = (log(param->tstop[i]) - eta[i])/sigma;
        g[i] = u - exp(u) - log(sigma);
        dg[i] = -(1-exp(u))/sigma;
        ddg[i] = -exp(u)/(sigma*sigma);
        ds[i] = -1 + dg[i]*(sigma*u);
        dds[i] = ddg[i]*pow(sigma*u,2) - dg[i]*(sigma*u);
        dsg[i] = ddg[i]*(sigma*u) - dg[i];
        break;
      case 3: case 4: // lognormal / normal
        u = (param->dist_code==3) ? 
        (log(param->tstop[i]) - eta[i])/sigma 
        : (param->tstop[i] - eta[i])/sigma;
        g[i] = R::dnorm(u,0,1,1) - log(sigma);
        dg[i] = u/sigma;
        ddg[i] = -1/(sigma*sigma);
        ds[i] = -1 + dg[i]*(sigma*u);
        dds[i] = ddg[i]*pow(sigma*u,2) - dg[i]*(sigma*u);
        dsg[i] = ddg[i]*(sigma*u) - dg[i];
        break;
      case 5: case 6: // loglogistic / logistic
        u = (param->dist_code==5) 
        ? (log(param->tstop[i]) - eta[i])/sigma 
        : (param->tstop[i] - eta[i])/sigma;
        g[i] = R::dlogis(u,0,1,1) - log(sigma);
        dg[i] = (1 - 2*R::plogis(u,0,1,0,0))/sigma;
        ddg[i] = -2*R::dlogis(u,0,1,0)/(sigma*sigma);
        ds[i] = -1 + dg[i]*(sigma*u);
        dds[i] = ddg[i]*pow(sigma*u,2) - dg[i]*(sigma*u);
        dsg[i] = ddg[i]*(sigma*u) - dg[i];
        break;
      }
      break;
      
    case 3: // interval censoring
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u1 = (log(param->tstart[i]) - eta[i])/sigma;
        u2 = (log(param->tstop[i]) - eta[i])/sigma;
        w1 = exp(u1); w2 = exp(u2);
        q1 = exp(-w1); q2 = exp(-w2);
        d1 = w1*q1; d2 = w2*q2;
        g[i] = log(q1 - q2);
        dg[i] = (d1 - d2)/(q1 - q2)/sigma;
        ddg[i] = -(d1*(1-w1) - d2*(1-w2))/(q1 - q2)/(sigma*sigma) - 
          pow(dg[i],2);
        ds[i] = (u1*d1 - u2*d2)/(q1 - q2);
        dds[i] = (u2*u2*(1-w2)*d2 - u1*u1*(1-w1)*d1)/(q1 - q2) - 
          ds[i]*(1+ds[i]);
        dsg[i] = (u2*(1-w2)*d2 - u1*(1-w1)*d1)/((q1 - q2)*sigma) - 
          dg[i]*(1+ds[i]);
        break;
      case 3: case 4: // lognormal / normal
        u1 = (param->dist_code==3) 
        ? (log(param->tstart[i]) - eta[i])/sigma 
        : (param->tstart[i] - eta[i])/sigma;
        u2 = (param->dist_code==3) 
          ? (log(param->tstop[i]) - eta[i])/sigma 
        : (param->tstop[i] - eta[i])/sigma;
        d1 = R::dnorm(u1,0,1,0); d2 = R::dnorm(u2,0,1,0);
        q1 = R::pnorm(u1,0,1,0,0); q2 = R::pnorm(u2,0,1,0,0);
        g[i] = log(q1 - q2);
        dg[i] = (d1 - d2)/(q1 - q2)/sigma;
        ddg[i] = (d1*u1 - d2*u2)/(q1 - q2)/(sigma*sigma) - pow(dg[i],2);
        ds[i] = (u1*d1 - u2*d2)/(q1 - q2);
        dds[i] = (-u2*u2*u2*d2 + u1*u1*u1*d1)/(q1 - q2) - ds[i]*(1+ds[i]);
        dsg[i] = (-u2*u2*d2 + u1*u1*d1)/((q1 - q2)*sigma) - dg[i]*(1+ds[i]);
        break;
      case 5: case 6: // loglogistic / logistic
        u1 = (param->dist_code==5) 
        ? (log(param->tstart[i]) - eta[i])/sigma 
        : (param->tstart[i] - eta[i])/sigma;
        u2 = (param->dist_code==5) ? 
        (log(param->tstop[i]) - eta[i])/sigma 
        : (param->tstop[i] - eta[i])/sigma;
        d1 = R::dlogis(u1,0,1,0); d2 = R::dlogis(u2,0,1,0);
        q1 = R::plogis(u1,0,1,0,0); q2 = R::plogis(u2,0,1,0,0);
        g[i] = log(q1 - q2);
        dg[i] = (d1 - d2)/(q1 - q2)/sigma;
        ddg[i] = -(d1*(2*q1-1) - d2*(2*q2-1))/(q1 - q2)/(sigma*sigma) - 
          pow(dg[i],2);
        ds[i] = (u1*d1 - u2*d2)/(q1 - q2);
        dds[i] = (u2*u2*d2*(2*q2-1) - u1*u1*d1*(2*q1-1))/(q1 - q2) - 
          ds[i]*(1+ds[i]);
        dsg[i] = (u2*d2*(2*q2-1) - u1*d1*(2*q1-1))/((q1 - q2)*sigma) - 
          dg[i]*(1+ds[i]);
        break;
      }
      break;
      
    case 2: // left censoring
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u = (log(param->tstop[i]) - eta[i])/sigma;
        g[i] = log(1 - exp(-exp(u)));
        dg[i] = -exp(u - exp(u))/(1 - exp(-exp(u)))/sigma;
        ddg[i] = (1 - exp(u) - exp(-exp(u)))*exp(u - exp(u))/
          pow((1 - exp(-exp(u)))*sigma, 2);
        ds[i] = dg[i]*(sigma*u);
        dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
        dsg[i] = ddg[i]*(sigma*u) - dg[i];
        break;
      case 3: case 4: // lognormal / normal
        u = (param->dist_code==3) 
        ? (log(param->tstop[i]) - eta[i])/sigma 
        : (param->tstop[i] - eta[i])/sigma;
        g[i] = R::pnorm(u,0,1,1,1);
        dg[i] = -R::dnorm(u,0,1,0)/R::pnorm(u,0,1,1,0)/sigma;
        ddg[i] = -u*R::dnorm(u,0,1,0)/(R::pnorm(u,0,1,1,0)*sigma*sigma) - 
          pow(dg[i],2);
        ds[i] = dg[i]*(sigma*u);
        dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
        dsg[i] = ddg[i]*(sigma*u) - dg[i];
        break;
      case 5: case 6: // loglogistic / logistic
        u = (param->dist_code==5) 
        ? (log(param->tstop[i]) - eta[i])/sigma 
        : (param->tstop[i] - eta[i])/sigma;
        g[i] = R::plogis(u,0,1,1,1);
        dg[i] = -R::plogis(u,0,1,0,0)/sigma;
        ddg[i] = -R::dlogis(u,0,1,0)/(sigma*sigma);
        ds[i] = dg[i]*(sigma*u);
        dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
        dsg[i] = ddg[i]*(sigma*u) - dg[i];
        break;
      }
      break;
      
    case 0: // right censoring
      switch (param->dist_code) {
      case 1: case 2: // exponential / weibull
        u = (log(param->tstart[i]) - eta[i])/sigma;
        g[i] = -exp(u);
        dg[i] = exp(u)/sigma;
        ddg[i] = -exp(u)/(sigma*sigma);
        ds[i] = dg[i]*(sigma*u);
        dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
        dsg[i] = ddg[i]*(sigma*u) - dg[i];
        break;
      case 3: case 4: // lognormal / normal
        u = (param->dist_code==3) 
        ? (log(param->tstart[i]) - eta[i])/sigma 
        : (param->tstart[i] - eta[i])/sigma;
        g[i] = R::pnorm(u,0,1,0,1);
        dg[i] = R::dnorm(u,0,1,0)/R::pnorm(u,0,1,0,0)/sigma;
        ddg[i] = u*R::dnorm(u,0,1,0)/(R::pnorm(u,0,1,0,0)*sigma*sigma) - 
          pow(dg[i],2);
        ds[i] = dg[i]*(sigma*u);
        dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
        dsg[i] = ddg[i]*(sigma*u) - dg[i];
        break;
      case 5: case 6: // loglogistic / logistic
        u = (param->dist_code==5) 
        ? (log(param->tstart[i]) - eta[i])/sigma 
        : (param->tstart[i] - eta[i])/sigma;
        g[i] = R::plogis(u,0,1,0,1);
        dg[i] = R::plogis(u,0,1,1,0)/sigma;
        ddg[i] = -R::dlogis(u,0,1,0)/(sigma*sigma);
        ds[i] = dg[i]*(sigma*u);
        dds[i] = ddg[i]*pow(sigma*u,2) - ds[i];
        dsg[i] = ddg[i]*(sigma*u) - dg[i];
        break;
      }
      break;
      
    default:
      stop("Unknown status: " + std::to_string(param->status[i]));
    }
  }
  
  return List::create(
    Named("g") = g,
    Named("dg") = dg,
    Named("ddg") = ddg,
    Named("ds") = ds,
    Named("dds") = dds,
    Named("dsg") = dsg
  );
}


// [[Rcpp::export]]
List liferegcpp(const DataFrame data,
                const StringVector& rep = "",
                const StringVector& stratum = "",
                const std::string time = "time",
                const std::string time2 = "",
                const std::string event = "event",
                const StringVector& covariates = "",
                const std::string weight = "",
                const std::string offset = "",
                const std::string id = "",
                const std::string dist = "weibull",
                const NumericVector& init = NA_REAL,
                const bool robust = false,
                const bool plci = false,
                const double alpha = 0.05,
                const int maxiter = 50,
                const double eps = 1.0e-9) {
  
  int n = data.nrows();
  int nvar = static_cast<int>(covariates.size()) + 1;
  if (nvar == 2 && (covariates[0] == "" || covariates[0] == "none")) {
    nvar = 1;
  }
  
  std::string dist1 = dist;
  std::for_each(dist1.begin(), dist1.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (dist1 == "log-logistic" || dist1 == "llogistic") {
    dist1 = "loglogistic";
  } else if (dist1 == "log-normal" || dist1 == "lnormal") {
    dist1 = "lognormal";
  } else if (dist1 == "gaussian") {
    dist1 = "normal";
  }
  
  int dist_code; 
  if (dist1 == "exponential") dist_code = 1;
  else if (dist1 == "weibull") dist_code = 2;
  else if (dist1 == "lognormal") dist_code = 3;
  else if (dist1 == "normal") dist_code = 4;
  else if (dist1 == "loglogistic") dist_code = 5;
  else if (dist1 == "logistic") dist_code = 6;
  else stop("invalid distribution: " + dist1);
  
  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(0);
  } else {
    List out = bygroup(data, rep);
    has_rep = 1;
    repn = out["index"];
    u_rep = DataFrame(out["lookup"]);
  }
  
  IntegerVector stratumn(n);
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    stratumn.fill(0);
  } else {
    List out = bygroup(data, stratum);
    stratumn = out["index"];
  }
  
  IntegerVector stratumn1 = unique(stratumn);
  int nstrata = static_cast<int>(stratumn1.size());
  int p = dist_code == 1 ? nvar : (nvar+nstrata);
  
  if (dist_code == 1 && nstrata > 1) {
    stop("Stratification is not valid with the exponential distribution");
  }
  
  bool has_time = hasVariable(data, time);
  if (!has_time) stop("data must contain the time variable");
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  for (int i=0; i<n; ++i) {
    if (!std::isnan(timen[i]) && (dist_code == 1 || dist_code == 2 || 
        dist_code == 3 || dist_code == 5) && timen[i] <= 0) {
      std::string str1 = "time must be positive for each subject for the";
      std::string str2 = "distribution";
      std::string errmsg = str1 + " " + dist1 + " " + str2;
      stop(errmsg);
    }
  }
  
  bool has_time2 = hasVariable(data, time2);
  NumericVector time2n(n);
  if (has_time2) {
    NumericVector time2nz = data[time2];
    time2n = clone(time2nz);
    for (int i=0; i<n; ++i) {
      if (!std::isnan(time2n[i]) && (dist_code == 1 || dist_code == 2 || 
          dist_code == 3 || dist_code == 5) && time2n[i] <= 0) {
        std::string str1 = "time2 must be positive for each subject for the";
        std::string str2 = "distribution";
        std::string errmsg = str1 + " " + dist1 + " " + str2;
        stop(errmsg);
      }
    }
  }
  
  bool has_event = hasVariable(data, event);
  if (!has_time2 && !has_event) {
    stop("data must contain the event variable for right censored data");
  }
  
  IntegerVector eventn(n);
  if (has_event) {
    IntegerVector eventnz = data[event];
    eventn = clone(eventnz);
    if (is_true(any((eventn != 1) & (eventn != 0)))) {
      stop("event must be 1 or 0 for each subject");
    }
  }
  
  NumericMatrix zn(n,nvar);
  for (int i=0; i<n; ++i) {
    zn(i,0) = 1; // intercept
  }
  
  for (int j=0; j<nvar-1; ++j) {
    String zj = covariates[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in covariates");
    }
    NumericVector u = data[zj];
    for (int i=0; i<n; ++i) {
      zn(i,j+1) = u[i];
    }
  }
  
  bool has_weight = hasVariable(data, weight);
  NumericVector weightn(n, 1.0);
  if (has_weight) {
    NumericVector weightnz = data[weight];
    weightn = clone(weightnz);
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }
  
  bool has_offset = hasVariable(data, offset);
  NumericVector offsetn(n);
  if (has_offset) {
    NumericVector offsetnz = data[offset];
    offsetn = clone(offsetnz);
  }
  
  // create the numeric id variable
  bool has_id = hasVariable(data, id);
  IntegerVector idn(n);
  if (!has_id) {
    idn = seq(0, n - 1);
  } else {
    SEXP col = data[id];
    SEXPTYPE col_type = TYPEOF(col);
    if (col_type == INTSXP) {
      IntegerVector v = col;
      IntegerVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else if (col_type == REALSXP) {
      NumericVector v = col;
      NumericVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else if (col_type == STRSXP) {
      StringVector v = col;
      StringVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else {
      stop("incorrect type for the id variable in the input data");
    }
  }
  
  // sort the data by rep
  if (has_rep) {
    IntegerVector order = seq(0, n-1);
    std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
      return repn[i] < repn[j];
    });
    
    repn = repn[order];
    stratumn = stratumn[order];
    timen = timen[order];
    time2n = time2n[order];
    eventn = eventn[order];
    weightn = weightn[order];
    offsetn = offsetn[order];
    idn = idn[order];
    zn = subset_matrix_by_row(zn, order);
  }
  
  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (int i=0; i<n; ++i) {
    if (repn[i] == NA_INTEGER || stratumn[i] == NA_INTEGER ||
        (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
        std::isnan(weightn[i]) || std::isnan(offsetn[i]) ||
        idn[i] == NA_INTEGER) {
      sub[i] = 0;
    }
    for (int j=0; j<nvar-1; ++j) {
      if (std::isnan(zn(i,j+1))) sub[i] = 0;
    }
  }
  
  IntegerVector order = which(sub);
  repn = repn[order];
  stratumn = stratumn[order];
  timen = timen[order];
  time2n = time2n[order];
  eventn = eventn[order];
  weightn = weightn[order];
  offsetn = offsetn[order];
  idn = idn[order];
  zn = subset_matrix_by_row(zn, order);
  n = sum(sub);
  if (n == 0) stop("no observation is left after removing missing values");
  
  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (int i=1; i<n; ++i) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);
  
  // variables in the output data sets
  // sumstat data set
  IntegerVector rep01 = seq(0,nreps-1);
  IntegerVector nobs(nreps), nevents(nreps);
  NumericVector loglik0(nreps), loglik1(nreps);
  IntegerVector niter(nreps);
  LogicalVector fails(nreps);
  
  // parest data set
  IntegerVector rep0(nreps*p);
  StringVector par0(nreps*p);
  NumericVector beta0(nreps*p), sebeta0(nreps*p), rsebeta0(nreps*p);
  NumericVector z0(nreps*p), expbeta0(nreps*p);
  NumericMatrix vbeta0(nreps*p,p), rvbeta0(nreps*p,p);
  NumericVector lb0(nreps*p), ub0(nreps*p), prob0(nreps*p);
  StringVector clparm0(nreps*p);
  
  NumericVector linear_predictors(n);
  
  double zcrit = R::qnorm(1-alpha/2,0,1,1,0);
  double xcrit = zcrit * zcrit;
  int bign0 = 0;
  for (int h=0; h<nreps; ++h) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());
    
    IntegerVector stratum1 = stratumn[q1];
    NumericVector time1 = timen[q1];
    NumericVector time21 = time2n[q1];
    IntegerVector event1 = eventn[q1];
    NumericVector weight1 = weightn[q1];
    NumericVector offset1 = offsetn[q1];
    IntegerVector id1 = idn[q1];
    NumericMatrix z1 = subset_matrix_by_row(zn, q1);
    
    // unify right censored data with interval censored data
    NumericVector tstart(n1), tstop(n1);
    if (!has_time2) {
      tstart = time1;
      for (int i=0; i<n1; ++i) {
        tstop[i] = event1[i] == 1 ? tstart[i] : NA_REAL;
      }
    } else {
      tstart = time1;
      tstop = time21;
    }
    
    IntegerVector status(n1);
    for (int i=0; i<n1; ++i) {
      if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
          tstart[i] == tstop[i]) {
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
    
    nobs[h] = n1;
    nevents[h] = sum(status == 1);
    
    if (nevents[h] == 0) {
      for (int i=0; i<p; ++i) {
        int k = h*p+i;
        rep0[k] = h;
        
        if (i==0) {
          par0[k] = "(Intercept)";
        } else if (i < nvar) {
          par0[k] = covariates[i-1];
        } else {
          if (nstrata == 1) {
            par0[k] = "Log(scale)";
          } else {
            std::string str1 = "Log(scale ";
            std::string str2 = ")";
            par0[k] = str1 + std::to_string(i-nvar+1) + str2;
          }
        }
        
        beta0[k] = NA_REAL;
        sebeta0[k] = 0;
        rsebeta0[k] = 0;
        z0[k] = NA_REAL;
        expbeta0[k] = NA_REAL;
        for (int j=0; j<p; ++j) {
          vbeta0(k,j) = 0;
          rvbeta0(k,j) = 0;
        }
        lb0[k] = NA_REAL;
        ub0[k] = NA_REAL;
        prob0[k] = NA_REAL;
        clparm0[k] = "Wald";
      }

      for (int i=0; i<n1; ++i) {
        linear_predictors[bign0+i] = offset1[i];
      }
      
      bign0 += n1;
      
      continue;
    }
    
    // exclude records with invalid status
    IntegerVector q2 = which(status != -1);
    int n2 = static_cast<int>(q2.size());
    
    if (n2 < n1) {
      stratum1 = stratum1[q2];
      tstart = tstart[q2];
      tstop = tstop[q2];
      status = status[q2];
      weight1 = weight1[q2];
      offset1 = offset1[q2];
      id1 = id1[q2];
      z1 = subset_matrix_by_row(z1,q2);
    }
    
    // intercept only model
    NumericVector time0(n2);
    for (int i=0; i<n2; ++i) {
      if (status[i] == 0 || status[i] == 1) { // right censoring or event
        time0[i] = tstart[i];
      } else if (status[i] == 2) { // left censoring
        time0[i] = tstop[i];
      } else if (status[i] == 3) { // interval censoring
        time0[i] = (tstart[i] + tstop[i])/2;
      }
    }
    
    NumericVector y0 = clone(time0);
    if (dist_code == 1 || dist_code == 2 || 
        dist_code == 3 || dist_code == 5) {
      y0 = log(y0);
    }
    
    double int0 = mean(y0);
    double logsig0 = log(sd(y0));
    
    NumericVector bint0(p);
    int ncolfit0 = dist_code == 1 ? 1 : nstrata + 1;
    IntegerVector colfit0(ncolfit0);
    if (dist_code == 1) {
      bint0[0] = int0;
      ncolfit0 = 1;
      colfit0[0] = 0;
    } else {
      bint0[0] = int0;
      for (int i=0; i<nstrata; ++i) {
        bint0[nvar+i] = logsig0;
      }
      
      colfit0[0] = 0;
      for (int i=0; i<nstrata; ++i) {
        colfit0[i+1] = nvar+i;
      }
    }
    
    // parameter estimates and standard errors for the null model
    aftparams param = {dist_code, stratum1, tstart, tstop, status, weight1,
                       offset1, z1, nstrata};
    
    List outint = liferegloop(p, bint0, &param, maxiter, eps,
                              colfit0, ncolfit0);
    NumericVector bint = outint["coef"];
    NumericMatrix vbint = as<NumericMatrix>(outint["var"]);
    
    NumericVector b(p);
    NumericMatrix vb(p,p);
    List out;
    
    if (nvar > 1) {
      IntegerVector colfit = seq(0,p-1);
      if (is_false(any(is_na(init))) && init.size() == p) {
        out = liferegloop(p, init, &param, maxiter, eps, colfit, p);
      } else {
        out = liferegloop(p, bint, &param, maxiter, eps, colfit, p);
      }
      
      bool fail = out["fail"];
      if (fail) {
        // obtain initial values for model parameters using OLS
        NumericVector y1 = y0 - offset1;
        NumericMatrix v1(nvar,nvar);  // XWX matrix
        NumericVector u1(nvar);       // XWY vector
        for (int i=0; i<n2; ++i) {
          for (int j=0; j<nvar; ++j) {
            for (int k=0; k<nvar; ++k) {
              v1(j,k) += weight1[i]*z1(i,j)*z1(i,k);
            }
            u1[j] += weight1[i]*z1(i,j)*y1[i];
          }
        }
        
        double toler = 1e-12;
        cholesky2(v1, nvar, toler);
        chsolve2(v1, nvar, u1);
        
        NumericVector binit(p);
        for (int j=0; j<nvar; ++j) {
          binit[j] = u1[j];
        }
        
        if (dist_code != 1) {
          double s = 0;
          for (int i=0; i<n2; ++i) {
            double r = y1[i] - std::inner_product(
              z1(i, _).begin(), z1(i, _).end(), u1.begin(), 0.0);
            s += weight1[i]*r*r;
          }
          s = 0.5*log(s/sum(weight1)*n2/(n2-nvar));  // log(sigma)
          
          for (int j=nvar; j<p; ++j) {
            binit[j] = s;
          }
        }
        
        // fit the model using the initial values
        out = liferegloop(p, binit, &param, maxiter, eps, colfit, p);
        fail = out["fail"];
      }
      
      if (fail) warning("The algorithm in liferegr did not converge");
      
      b = out["coef"];
      vb = as<NumericMatrix>(out["var"]);
    } else {
      b = bint;
      vb = vbint;
      out = outint;
    }
    
    NumericVector seb(p);
    for (int j=0; j<p; ++j) {
      seb[j] = sqrt(vb(j,j));
    }
    
    for (int i=0; i<p; ++i) {
      rep0[h*p+i] = h;
      
      if (i==0) {
        par0[h*p+i] = "(Intercept)";
      } else if (i < nvar) {
        par0[h*p+i] = covariates[i-1];
      } else {
        if (nstrata == 1) {
          par0[h*p+i] = "Log(scale)";
        } else {
          std::string str1 = "Log(scale ";
          std::string str2 = ")";
          par0[h*p+i] = str1 + std::to_string(i-nvar+1) + str2;
        }
      }
      
      beta0[h*p+i] = b[i];
      sebeta0[h*p+i] = seb[i];
      for (int j=0; j<p; ++j) {
        vbeta0(h*p+i,j) = vb(i,j);
      }
    }
    
    for (int i=0; i<n2; ++i) {
      linear_predictors[bign0+q2[i]] = offset1[i];
      for (int j=0; j<nvar; ++j) {
        linear_predictors[bign0+q2[i]] += b[j]*z1(i,j);
      }
    }
    
    bign0 += n1;
    
    niter[h] = out["iter"];
    fails[h] = out["fail"];
    
    // robust variance estimates
    NumericVector rseb(p);  // robust standard error for betahat
    if (robust) {
      NumericMatrix ressco = f_ressco_1(p, b, &param);
      
      int nr; // number of rows in the score residual matrix
      if (!has_id) {
        for (int i=0; i<n2; ++i) {
          for (int j=0; j<p; ++j) {
            ressco(i,j) = weight1[i]*ressco(i,j);
          }
        }
        nr = n2;
      } else { // need to sum up score residuals by id
        IntegerVector order = seq(0, n2-1);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return id1[i] < id1[j];
        });
        
        IntegerVector id2 = id1[order];
        IntegerVector idx(1,0);
        for (int i=1; i<n2; ++i) {
          if (id2[i] != id2[i-1]) {
            idx.push_back(i);
          }
        }
        
        int nids = static_cast<int>(idx.size());
        idx.push_back(n2);
        
        NumericVector weight2 = weight1[order];
        
        NumericMatrix ressco2(nids,p);
        for (int i=0; i<nids; ++i) {
          for (int j=0; j<p; ++j) {
            for (int k=idx[i]; k<idx[i+1]; ++k) {
              ressco2(i,j) += weight2[k]*ressco(order[k],j);
            }
          }
        }
        
        ressco = ressco2;  // update the score residuals
        nr = nids;
      }
      
      NumericMatrix D(nr,p); // DFBETA
      for (int i=0; i<nr; ++i) {
        for (int j=0; j<p; ++j) {
          for (int k=0; k<p; ++k) {
            D(i,j) += ressco(i,k)*vb(k,j);
          }
        }
      }
      
      NumericMatrix rvb(p,p); // robust variance matrix for betahat
      for (int j=0; j<p; ++j) {
        for (int k=0; k<p; ++k) {
          for (int i=0; i<nr; ++i) {
            rvb(j,k) += D(i,j)*D(i,k);
          }
        }
      }
      
      for (int i=0; i<p; ++i) {
        rseb[i] = sqrt(rvb(i,i));
      }
      
      for (int i=0; i<p; ++i) {
        rsebeta0[h*p+i] = rseb[i];
        for (int j=0; j<p; ++j) {
          rvbeta0(h*p+i,j) = rvb(i,j);
        }
      }
    }
    
    // profile likelihood confidence interval for regression coefficients
    NumericVector lb(p), ub(p), prob(p);
    StringVector clparm(p);
    
    if (plci) {
      double lmax = out["loglik"];
      double l0 = lmax - 0.5*xcrit;
      
      for (int k=0; k<p; ++k) {
        lb[k] = liferegplloop(p, b, &param, maxiter, eps, k, -1, l0);
        ub[k] = liferegplloop(p, b, &param, maxiter, eps, k, 1, l0);
        
        IntegerVector colfit1(p-1);
        for (int i=0; i<k; ++i) {
          colfit1[i] = i;
        }
        for (int i=k+1; i<p; ++i) {
          colfit1[i-1] = i;
        }
        
        NumericVector b0(p);
        List out0 = liferegloop(p, b0, &param, maxiter, eps, colfit1, p-1);
        double lmax0 = out0["loglik"];
        prob[k] = R::pchisq(-2*(lmax0 - lmax), 1, 0, 0);
        clparm[k] = "PL";
      }
    } else {
      for (int k=0; k<p; ++k) {
        if (!robust) {
          lb[k] = b[k] - zcrit*seb[k];
          ub[k] = b[k] + zcrit*seb[k];
          prob[k] = R::pchisq(pow(b[k]/seb[k], 2), 1, 0, 0);
        } else {
          lb[k] = b[k] - zcrit*rseb[k];
          ub[k] = b[k] + zcrit*rseb[k];
          prob[k] = R::pchisq(pow(b[k]/rseb[k], 2), 1, 0, 0);
        }
        clparm[k] = "Wald";
      }
    }
    
    for (int i=0; i<p; ++i) {
      lb0[h*p+i] = lb[i];
      ub0[h*p+i] = ub[i];
      prob0[h*p+i] = prob[i];
      clparm0[h*p+i] = clparm[i];
    }
    
    // log-likelihoods
    loglik0[h] = outint["loglik"];
    loglik1[h] = out["loglik"];
  }
  
  expbeta0 = exp(beta0);
  if (!robust) z0 = beta0/sebeta0;
  else z0 = beta0/rsebeta0;
  
  List sumstat = List::create(
    _["n"] = nobs,
    _["nevents"] = nevents,
    _["loglik0"] = loglik0,
    _["loglik1"] = loglik1,
    _["niter"] = niter,
    _["dist"] = dist1,
    _["p"] = p,
    _["nvar"] = nvar-1,
    _["robust"] = robust,
    _["fail"] = fails);
  
  List parest = List::create(
    _["param"] = par0,
    _["beta"] = beta0,
    _["sebeta"] = robust ? rsebeta0 : sebeta0,
    _["z"] = z0,
    _["expbeta"] = expbeta0,
    _["vbeta"] = robust ? rvbeta0 : vbeta0,
    _["lower"] = lb0,
    _["upper"] = ub0,
    _["p"] = prob0,
    _["method"] = clparm0);
  
  if (robust) {
    parest.push_back(sebeta0, "sebeta_naive");
    parest.push_back(vbeta0, "vbeta_naive");
  }
  
  if (has_rep) {
    for (int i=0; i<p_rep; ++i) {
      std::string s = as<std::string>(rep[i]);
      SEXP col = u_rep[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == INTSXP) {
        IntegerVector v = col;
        sumstat.push_back(v[rep01], s);
        parest.push_back(v[rep0], s);
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        sumstat.push_back(v[rep01], s);
        parest.push_back(v[rep0], s);
      } else if (col_type == STRSXP) {
        StringVector v = col;
        sumstat.push_back(v[rep01], s);
        parest.push_back(v[rep0], s);
      } else {
        stop("Unsupported type for rep variable" + s);
      }
    }
  }
  
  List result = List::create(
    _["sumstat"] = as<DataFrame>(sumstat),
    _["parest"] = as<DataFrame>(parest),
    _["linear_predictors"] = linear_predictors);
  
  return result;
}


// [[Rcpp::export]]
NumericMatrix residuals_liferegcpp(const NumericVector& beta,
                                   const NumericMatrix& vbeta,
                                   DataFrame data,
                                   const StringVector& stratum = "",
                                   const std::string time = "time",
                                   const std::string time2 = "",
                                   const std::string event = "event",
                                   const StringVector& covariates = "",
                                   const std::string weight = "",
                                   const std::string offset = "",
                                   const std::string id = "",
                                   const std::string dist = "weibull",
                                   const std::string type = "response",
                                   const bool collapse = false,
                                   const bool weighted = false) {
  
  int n = data.nrows();
  int nvar = static_cast<int>(covariates.size()) + 1;
  if (nvar == 2 && (covariates[0] == "" || covariates[0] == "none")) {
    nvar = 1;
  }
  
  std::string dist1 = dist;
  std::for_each(dist1.begin(), dist1.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (dist1 == "log-logistic" || dist1 == "llogistic") {
    dist1 = "loglogistic";
  } else if  (dist1 == "log-normal" || dist1 == "lnormal") {
    dist1 = "lognormal";
  } else if (dist1 == "gaussian") {
    dist1 = "normal";
  }
  
  int dist_code; 
  if (dist1 == "exponential") dist_code = 1;
  else if (dist1 == "weibull") dist_code = 2;
  else if (dist1 == "lognormal") dist_code = 3;
  else if (dist1 == "normal") dist_code = 4;
  else if (dist1 == "loglogistic") dist_code = 5;
  else if (dist1 == "logistic") dist_code = 6;
  else stop("invalid distribution: " + dist1);
  
  IntegerVector stratumn(n);
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    stratumn.fill(0);
  } else {
    List out = bygroup(data, stratum);
    stratumn = out["index"];
  }
  
  IntegerVector stratumn1 = unique(stratumn);
  int nstrata = static_cast<int>(stratumn1.size());
  int p = dist_code == 1 ? nvar : (nvar+nstrata);
  
  if (dist_code == 1 && nstrata > 1) {
    stop("Stratification is not valid with the exponential distribution");
  }
  
  bool has_time = hasVariable(data, time);
  if (!has_time) stop("data must contain the time variable");
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  for (int i=0; i<n; ++i) {
    if (!std::isnan(timen[i]) && (dist_code == 1 || dist_code == 2 || 
        dist_code == 3 || dist_code == 5) && timen[i] <= 0) {
      std::string str1 = "time must be positive for each subject for the";
      std::string str2 = "distribution";
      std::string errmsg = str1 + " " + dist1 + " " + str2;
      stop(errmsg);
    }
  }
  
  bool has_time2 = hasVariable(data, time2);
  NumericVector time2n(n);
  if (has_time2) {
    NumericVector time2nz = data[time2];
    time2n = clone(time2nz);
    for (int i=0; i<n; ++i) {
      if (!std::isnan(time2n[i]) && (dist_code == 1 || dist_code == 2 || 
          dist_code == 3 || dist_code == 5) && time2n[i] <= 0) {
        std::string str1 = "time2 must be positive for each subject for the";
        std::string str2 = "distribution";
        std::string errmsg = str1 + " " + dist1 + " " + str2;
        stop(errmsg);
      }
    }
  }
  
  bool has_event = hasVariable(data, event);
  if (!has_time2 && !has_event) {
    stop("data must contain the event variable for right censored data");
  }
  
  IntegerVector eventn(n);
  if (has_event) {
    IntegerVector eventnz = data[event];
    eventn = clone(eventnz);
    if (is_true(any((eventn != 1) & (eventn != 0)))) {
      stop("event must be 1 or 0 for each subject");
    }
  }
  
  NumericMatrix zn(n,nvar);
  for (int i=0; i<n; ++i) {
    zn(i,0) = 1; // intercept
  }
  
  for (int j=0; j<nvar-1; ++j) {
    String zj = covariates[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in covariates");
    }
    NumericVector u = data[zj];
    for (int i=0; i<n; ++i) {
      zn(i,j+1) = u[i];
    }
  }
  
  bool has_weight = hasVariable(data, weight);
  NumericVector weightn(n, 1.0);
  if (has_weight) {
    NumericVector weightnz = data[weight];
    weightn = clone(weightnz);
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }
  
  bool has_offset = hasVariable(data, offset);
  NumericVector offsetn(n);
  if (has_offset) {
    NumericVector offsetnz = data[offset];
    offsetn = clone(offsetnz);
  }
  
  // create the numeric id variable
  bool has_id = hasVariable(data, id);
  IntegerVector idn(n);
  if (!has_id) {
    idn = seq(0, n - 1);
  } else {
    SEXP col = data[id];
    SEXPTYPE col_type = TYPEOF(col);
    if (col_type == INTSXP) {
      IntegerVector v = col;
      IntegerVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else if (col_type == REALSXP) {
      NumericVector v = col;
      NumericVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else if (col_type == STRSXP) {
      StringVector v = col;
      StringVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else {
      stop("incorrect type for the id variable in the input data");
    }
  }
  
  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (int i=0; i<n; ++i) {
    if (stratumn[i] == NA_INTEGER || idn[i] == NA_INTEGER ||
        (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
        std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
      sub[i] = 0;
    }
    for (int j=0; j<nvar-1; ++j) {
      if (std::isnan(zn(i,j+1))) sub[i] = 0;
    }
  }
  
  IntegerVector q1 = which(sub);
  stratumn = stratumn[q1];
  timen = timen[q1];
  time2n = time2n[q1];
  eventn = eventn[q1];
  weightn = weightn[q1];
  offsetn = offsetn[q1];
  idn = idn[q1];
  zn = subset_matrix_by_row(zn, q1);
  int n1 = sum(sub);
  if (n1 == 0) stop("no observation is left after removing missing values");
  
  // unify right censored data with interval censored data
  NumericVector tstart(n1), tstop(n1);
  if (!has_time2) {
    tstart = timen;
    for (int i=0; i<n1; ++i) {
      tstop[i] = eventn[i] == 1 ? tstart[i] : NA_REAL;
    }
  } else {
    tstart = timen;
    tstop = time2n;
  }
  
  IntegerVector status(n1);
  for (int i=0; i<n1; ++i) {
    if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
        tstart[i] == tstop[i]) {
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
  
  // exclude records with invalid status
  IntegerVector q2 = which(status != -1);
  int n2 = static_cast<int>(q2.size());
  
  if (n2 < n1) {
    stratumn = stratumn[q2];
    tstart = tstart[q2];
    tstop = tstop[q2];
    status = status[q2];
    weightn = weightn[q2];
    offsetn = offsetn[q2];
    idn = idn[q2];
    zn = subset_matrix_by_row(zn,q2);
  }
  
  NumericVector eta(n2);
  for (int i = 0; i < n2; ++i) {
    eta[i] = offsetn[i];
    for (int j=0; j<nvar; ++j) {
      eta[i] += beta[j]*zn(i,j);
    }
  }
  
  NumericVector sig(n2, 1.0);
  if (dist_code != 1) {
    for (int i = 0; i < n2; ++i) {
      int j = stratumn[i] + nvar;
      sig[i] = exp(beta[j]);
    }
  }
  
  int K = 1;
  if (type == "dfbeta" || type == "dfbetas") {
    K = p;
  } else if (type == "matrix") {
    K = 6;
  }
  
  // Map type to integer code
  int type_code;
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
  else stop("invalid type of residuals" + type);
  
  // rr: residual matrix
  NumericMatrix rr(n2, K);
  
  switch (type_code) {
  
  case 1: { // response
    NumericVector yhat0(n2);
    for (int i = 0; i < n2; ++i) {
      switch (status[i]) {
      case 0: case 1: // right-censored or event
        switch (dist_code) {
        case 1: case 2: case 3: case 5: yhat0[i] = log(tstart[i]); break;
        default: yhat0[i] = tstart[i];
        }
        break;
      case 2: // left-censored
        switch (dist_code) {
        case 1: case 2: case 3: case 5: yhat0[i] = log(tstop[i]); break;
        default: yhat0[i] = tstop[i];
        }
        break;
      default: // interval-censored
        switch (dist_code) {
        case 1: case 2: {
          double width = (log(tstop[i]) - log(tstart[i])) / sig[i];
          yhat0[i] = log(tstart[i]) - sig[i] * log(width / (exp(width) - 1));
          break;
        }
        case 3: case 5: 
          yhat0[i] = 0.5 * (log(tstart[i]) + log(tstop[i])); break;
        default: yhat0[i] = 0.5 * (tstart[i] + tstop[i]);
        }
      }
      
      switch (dist_code) {
      case 1: case 2: case 3: case 5: 
        rr(i,0) = exp(yhat0[i]) - exp(eta[i]); break;
      default: rr(i,0) = yhat0[i] - eta[i];
      }
    }
    break;
  }
    
  case 2: { // martingale
    if (dist_code == 4 || dist_code == 6) 
      stop("incorrect type of distribution" + dist1 + 
        " for martingale residuals");
    for (int i = 0; i < n2; ++i) {
      if (status[i] == 0 || status[i] == 1) {
        double y = (log(tstart[i]) - eta[i]) / sig[i];
        switch (dist_code) {
        case 1: case 2: rr(i,0) = (status[i] == 1) - exp(y); break;
        case 3: rr(i,0) = (status[i] == 1) + R::pnorm(y,0,1,0,1); break;
        case 5: rr(i,0) = (status[i] == 1) + R::plogis(y,0,1,0,1); break;
        }
      } else rr(i,0) = NA_REAL;
    }
    break;
  }
    
  default: { // other types: deviance, working, dfbeta, ld*, matrix
    aftparams param = {dist_code, stratumn, tstart, tstop, status, weightn,
                       offsetn, zn, nstrata};
    List der = f_ld_1(eta, sig, &param);
    NumericVector g = der["g"], dg = der["dg"], ddg = der["ddg"];
    NumericVector ds = der["ds"], dds = der["dds"], dsg = der["dsg"];
    
    switch (type_code) {
    
    case 3: { // deviance
      NumericVector loglik(n2);
      for (int i = 0; i < n2; ++i) {
        switch (status[i]) {
        case 0: case 2: loglik[i] = 0; break; // right or left censored
        case 1:  // event
          switch (dist_code) {
          case 1: case 2: loglik[i] = -log(sig[i]) - 1; break;
          case 3: case 4: loglik[i] = -log(sqrt(2*M_PI)*sig[i]); break;
          default: loglik[i] = -log(4*sig[i]);
          }
          break;
        default: { // interval censored
            double width;
            switch (dist_code) {
            case 1: case 2: width = (log(tstop[i]) - log(tstart[i])) / sig[i];
              loglik[i] = - width/(exp(width)-1) + log(1 - exp(-width)); break;
            case 3: width = (log(tstop[i]) - log(tstart[i])) / sig[i];
              loglik[i] = log(2*R::pnorm(width/2,0,1,1,0) - 1); break;
            case 4: width = (tstop[i] - tstart[i]) / sig[i];
              loglik[i] = log(2*R::pnorm(width/2,0,1,1,0) - 1); break;
            case 5: width = (log(tstop[i]) - log(tstart[i])) / sig[i];
              loglik[i] = log((exp(width/2)-1)/(exp(width/2)+1)); break;
            default: width = (tstop[i] - tstart[i]) / sig[i];
            loglik[i] = log((exp(width/2)-1)/(exp(width/2)+1));
            }
          }
        }
        double val = -dg[i]/ddg[i];
        rr(i,0) = ((val>0) - (val<0))*sqrt(2*(loglik[i] - g[i]));
      }
      break;
    }
      
    case 4: { // working
      for (int i=0; i<n2; ++i) rr(i,0) = -dg[i]/ddg[i];
      break;
    }
      
    case 5: case 6: case 7: { // dfbeta, dfbetas, ldcase
      for (int i=0; i<n2; ++i) {
      NumericVector score(p), resid(p);
      for (int j=0; j<nvar; ++j) score[j] = dg[i]*zn(i,j);
      for (int j=nvar; j<p; ++j) score[j] = stratumn[i]==j-nvar ? ds[i] : 0;
      for (int k=0; k<p; ++k) for (int j=0; j<p; ++j) 
        resid[k] += score[j]*vbeta(j,k);
      if (type_code==6) for (int k=0; k<p; ++k) resid[k] /= sqrt(vbeta(k,k));
      if (type_code==7) for (int k=0; k<p; ++k) rr(i,0) += resid[k]*score[k];
      else for (int k=0; k<p; ++k) rr(i,k) = resid[k];
    }
      break;
    }
      
    case 8: { // ldresp
      for (int i=0; i<n2; ++i) {
      NumericVector rscore(p), temp(p);
      for (int j=0; j<nvar; ++j) rscore[j] = -ddg[i]*zn(i,j)*sig[i];
      for (int j=nvar; j<p; ++j) 
        rscore[j] = stratumn[i]==j-nvar ? -dsg[i]*sig[i] : 0;
      for (int k=0; k<p; ++k) for (int j=0; j<p; ++j) 
        temp[k] += rscore[j]*vbeta(j,k);
      for (int k=0; k<p; ++k) rr(i,0) += temp[k]*rscore[k];
    }
      break;
    }
      
    case 9: { // ldshape
      for (int i=0; i<n2; ++i) {
      NumericVector sscore(p), temp(p);
      for (int j=0; j<nvar; ++j) sscore[j] = dsg[i]*zn(i,j);
      for (int j=nvar; j<p; ++j) sscore[j] = stratumn[i]==j-nvar ? dds[i] : 0;
      for (int k=0; k<p; ++k) for (int j=0; j<p; ++j) 
        temp[k] += sscore[j]*vbeta(j,k);
      for (int k=0; k<p; ++k) rr(i,0) += temp[k]*sscore[k];
    }
      break;
    }
      
    case 10: { // matrix
      for (int i=0; i<n2; ++i) {
      rr(i,0)=g[i]; rr(i,1)=dg[i]; rr(i,2)=ddg[i];
      rr(i,3)=ds[i]; rr(i,4)=dds[i]; rr(i,5)=dsg[i];
    }
      break;
    }
      
    }
  } // end default
  } // end switch (type_code)
  
  // case weights
  if (weighted) {
    for (int i = 0; i < n2; ++i) {
      for (int k=0; k<K; ++k) {
        rr(i,k) *= weightn[i];
      }
    }
  }
  
  // collapse if needed
  if (collapse) {
    // order data by id
    IntegerVector order = seq(0, n2-1);
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      return idn[i] < idn[j];
    });
    
    IntegerVector id2 = idn[order];
    IntegerVector idx(1,0);
    for (int i=1; i<n2; ++i) {
      if (id2[i] != id2[i-1]) {
        idx.push_back(i);
      }
    }
    
    int nids = static_cast<int>(idx.size());
    idx.push_back(n2);
    
    // collapse over id
    NumericMatrix rr2(nids,K);
    for (int i=0; i<nids; ++i) {
      for (int k=0; k<K; ++k) {
        for (int j=idx[i]; j<idx[i+1]; ++j) {
          rr2(i,k) += rr(order[j],k);
        }
      }
    }
    
    rr = rr2;
  }
  
  return rr;
}



// all-in-one function for log-likelihood, score, and information matrix
// for the Cox model with or without Firth's correction
List f_der_2(int p, const NumericVector& par, void* ex, bool firth) {
  coxparams* param = (coxparams*) ex;
  
  const int nused = param->nused;
  const int method = param->method;
  
  // Precompute eta and exp(eta)
  NumericVector eta(nused);
  NumericVector exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    double val = param->offset[person];
    for (int i = 0; i < p; ++i) {
      val += par[i] * param->z(person,i);
    }
    eta[person] = val;
    exp_eta[person] = exp(val);
  }
  
  double loglik = 0.0;        // log-likelihood
  NumericVector u(p);         // score vector
  NumericMatrix imat(p,p);    // information matrix
  NumericMatrix dimat(p,p*p); // tensor for third order derivatives
  NumericVector a(p);         // s1(beta,k,t)
  NumericVector a2(p);        // sum of w*exp(zbeta)*z for the deaths
  NumericMatrix cmat(p,p);    // s2(beta,k,t)
  NumericMatrix cmat2(p,p);   // sum of w*exp(zbeta)*z*z' for the deaths
  NumericMatrix dmat(p,p*p);  // q2(beta,k,t)
  NumericMatrix dmat2(p,p*p); // sum of w*exp(zbeta)*z*z*z' for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  int ndead = 0;              // number of deaths at this time point
  
  
  int istrata = param->strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    // Reset when entering a new stratum
    if (param->strata[person] != istrata) {
      istrata = param->strata[person];
      i1 = person;
      denom = 0.0;
      
      for (int i = 0; i < p; ++i) {
        a[i] = 0.0;
        for (int j = 0; j <= i; ++j) {
          cmat(i,j) = 0.0;
          if (firth) {
            for (int k = 0; k <= j; ++k) {
              dmat(i,k*p+j) = 0.0;
            }
          }
        }
      }
    }
    
    const double dtime = param->tstop[person];
    
    // Process all persons tied at this dtime
    for (; person < nused && param->tstop[person] == dtime &&
         param->strata[person] == istrata; ++person) {
      
      const double w = param->weight[person];
      const double r = w * exp_eta[person];
      
      if (param->event[person] == 0) {
        denom += r;
        
        for (int i = 0; i < p; ++i) {
          const double zi = param->z(person,i);
          a[i] += r * zi;
          for (int j = 0; j <= i; ++j) {
            const double zj = param->z(person,j);
            cmat(i,j) += r * zi * zj;
            if (firth) {
              for (int k = 0; k <= j; ++k) {
                dmat(i,k*p+j) += r * zi * zj * param->z(person,k);
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
          const double zi = param->z(person,i);
          a2[i] += r * zi;
          u[i] += w * zi;
          for (int j = 0; j <= i; ++j) {
            const double zj = param->z(person,j);
            cmat2(i,j) += r * zi * zj;
            if (firth) {
              for (int k = 0; k <= j; ++k) {
                dmat2(i,k*p+j) += r * zi * zj * param->z(person,k);
              }
            }
          }
        }
      }
    }
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      
      const double r = param->weight[p1] * exp_eta[p1];
      denom -= r;
      for (int i = 0; i < p; ++i) {
        const double zi = param->z(p1,i);
        a[i] -= r * zi;
        for (int j = 0; j <= i; ++j) {
          const double zj = param->z(p1,j);
          cmat(i,j) -= r * zi * zj;
          if (firth) {
            for (int k = 0; k <= j; ++k) {
              dmat(i,k*p+j) -= r * zi * zj * param->z(p1,k);
            }
          }
        }
      }
    }
    
    // Add contributions for deaths at this time
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        loglik -= deadwt * log(denom);
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          const double xbar = a[i] / denom;
          u[i] -= deadwt * xbar;
          for (int j = 0; j <= i; ++j) {
            cmat(i,j) += cmat2(i,j);
            imat(i,j) += deadwt * (cmat(i,j) - xbar * a[j]) / denom;
            if (firth) {
              for (int k = 0; k <= j; ++k) {
                dmat(i,k*p+j) += dmat2(i,k*p+j);
                dimat(i,k*p+j) += deadwt * (dmat(i,k*p+j) - 
                  (cmat(i,j)*a[k] + cmat(i,k)*a[j] + cmat(j,k)*a[i]) / denom + 
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
          loglik -= meanwt * log(denom);
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            const double xbar = a[i] / denom;
            u[i] -= meanwt * xbar;
            for (int j = 0; j <= i; ++j) {
              cmat(i,j) += cmat2(i,j) / ndead;
              imat(i,j) += meanwt * (cmat(i,j) - xbar * a[j]) / denom;
              if (firth) {
                for (int k = 0; k <= j; ++k) {
                  dmat(i,k*p+j) += dmat2(i,k*p+j) / ndead;
                  dimat(i,k*p+j) += meanwt * (dmat(i,k*p+j) - 
                    (cmat(i,j)*a[k] + cmat(i,k)*a[j] + cmat(j,k)*a[i])/denom + 
                    2.0 * a[i] * a[j] * a[k] / (denom * denom)) / denom;
                }
              }
            }
          }
        }
      }
      
      // Reset after processing deaths
      ndead = 0;
      deadwt = denom2 = 0.0;
      for (int i = 0; i < p; ++i) {
        a2[i] = 0;
        for (int j = 0; j <= i; ++j) {
          cmat2(i,j) = 0;
          if (firth) {
            for (int k = 0; k <= j; ++k) {
              dmat2(i,k*p+j) = 0;
            }
          }
        }
      }
    }
  }
  
  
  // fill the symmetric elements of the information matrix
  for (int i = 0; i < p - 1; ++i)
    for (int j = i+1; j < p; ++j)
      imat(i,j) = imat(j,i);
  
  
  // fill the symmetric elements of the tensor array
  if (firth) {
    for (int i = 0; i < p-1; ++i) 
      for (int j = i+1; j < p; ++j) 
        for (int k = 0; k <= i; ++k) 
          dimat(i,k*p+j) = dimat(j,k*p+i);
    
    for (int j = 0; j < p-1; ++j)
      for (int i = j; i < p; ++i)
        for (int k = j+1; k < p; ++k)
          dimat(i,k*p+j) = dimat(i,j*p+k);
    
    for (int i = 0; i < p-1; ++i)
      for (int j = i+1; j < p; ++j)
        for (int k = i+1; k < p; ++k)
          dimat(i,k*p+j) = dimat(k,i*p+j);
  }
  
  
  List result;
  
  // Firth adjustment
  if (p > 0 && firth) {
    // obtain the determinant of information matrix
    NumericMatrix imat0 = clone(imat);
    double toler = 1e-12;
    cholesky2(imat0, p, toler);
    
    double v = 0;
    for (int i=0; i<p; ++i) {
      v += log(imat0(i,i));
    }
    
    // penalized log-likelihood adjustment
    double penloglik = loglik + 0.5*v;
    
    // compute the bias adjustment to the score function
    NumericMatrix y(p,p);
    NumericVector g(p);
    
    for (int k = 0; k < p; ++k) {
      // partial derivative of the information matrix w.r.t. beta[k]
      for (int i = 0; i < p; ++i) {
        for (int j = 0; j < p; ++j) {
          y(i,j) = dimat(i,k*p+j);
        }
      }
      
      // solve(imat, y)
      for (int h = 0; h < p; ++h) {
        for (int i = 0; i < p; ++i) {
          double temp = y(i,h);
          for (int j = 0; j < i; ++j)
            temp -= y(j,h) * imat0(j,i);
          y(i,h) = temp;
        }
        
        for (int i = p-1; i >= 0; --i) {
          if (imat0(i,i) == 0) y(i,h) = 0;
          else {
            double temp = y(i,h) / imat0(i,i);
            for (int j = i+1; j < p; ++j)
              temp -= y(j,h)*imat0(i,j);
            y(i,h) = temp;
          }
        }
      }
      
      // trace
      for (int i = 0; i < p; ++i) g[k] += y(i,i);
      
      g[k] = u[k] + 0.5 * g[k];
    }
    
    result = List::create(
      _["loglik"] = penloglik,
      _["score"] = g,
      _["imat"] = imat,
      _["regloglik"] = loglik,
      _["regscore"] = u
    );
  } else {
    result = List::create(
      _["loglik"] = loglik
    );
    
    if (p > 0) {
      result.push_back(u, "score");
      result.push_back(imat, "imat");
    }
  }
  
  return result;
}


// underlying optimization algorithm for phreg
List phregloop(int p, const NumericVector& par, void *ex,
               int maxiter, double eps, bool firth,
               const IntegerVector& colfit, int ncolfit) {
  coxparams *param = (coxparams *) ex;
  
  int iter, halving = 0;
  bool fail = false;
  double toler = 1e-12;
  
  NumericVector beta(p), newbeta(p);
  double loglik, newlk = 0;
  NumericVector u(p);
  NumericMatrix imat(p,p);
  NumericVector u1(ncolfit);
  NumericMatrix imat1(ncolfit, ncolfit);
  NumericMatrix z1 = param->z;
  
  // --- first step ---
  beta = clone(par);
  
  List der = f_der_2(p, beta, param, firth);
  loglik = der["loglik"];
  u = der["score"];
  imat = as<NumericMatrix>(der["imat"]);
  
  for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
  
  for (int i = 0; i < ncolfit; ++i) 
    for (int j = 0; j < ncolfit; ++j) 
      imat1(i,j) = imat(colfit[i], colfit[j]);
  
  cholesky2(imat1, ncolfit, toler);
  chsolve2(imat1, ncolfit, u1);
  
  u.fill(0.0);
  for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
  newbeta = beta + u;
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_2(p, newbeta, param, firth);
    newlk = der["loglik"];
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && halving == 0 && fabs(1 - (loglik/newlk)) < eps) break;
    
    if (fail || newlk < loglik) {
      ++halving; // adjust step size if likelihood decreases
      for (int i = 0; i < p; ++i) {
        newbeta[i] = 0.5 * (beta[i] + newbeta[i]);
      }
      continue;
    }
    
    // --- update ---
    halving = 0;
    beta = clone(newbeta);
    loglik = newlk;
    u = der["score"];
    imat = as<NumericMatrix>(der["imat"]);
    
    for (int i = 0; i < ncolfit; ++i) u1[i] = u[colfit[i]];
    
    for (int i = 0; i < ncolfit; ++i) 
      for (int j = 0; j < ncolfit; ++j) 
        imat1(i,j) = imat(colfit[i], colfit[j]);
    
    cholesky2(imat1, ncolfit, toler);
    chsolve2(imat1, ncolfit, u1);
    
    u.fill(0.0);
    for (int i = 0; i < ncolfit; ++i) u[colfit[i]] = u1[i];
    newbeta = beta + u;
  }
  
  if (iter == maxiter) fail = true;
  
  // --- final variance calculation ---
  imat = as<NumericMatrix>(der["imat"]);
  for (int i = 0; i < ncolfit; ++i) 
    for (int j = 0; j < ncolfit; ++j) 
      imat1(i,j) = imat(colfit[i], colfit[j]);
  
  NumericMatrix var1 = invsympd(imat1, ncolfit, toler);
  NumericMatrix var(p,p);
  for (int i = 0; i < ncolfit; ++i)
    for (int j = 0; j < ncolfit; ++j)
      var(colfit[i], colfit[j]) = var1(i,j);
  
  List result = List::create(
    Named("coef") = newbeta,
    Named("iter") = iter,
    Named("var") = var,
    Named("loglik") = newlk,
    Named("fail") = fail);
  
  if (firth) {
    double regloglik = as<double>(der["regloglik"]);
    result.push_back(regloglik, "regloglik");
  }
  
  return result;
}


// confidence limit of profile likelihood method
double phregplloop(int p, const NumericVector& par, void *ex,
                   int maxiter, double eps, bool firth,
                   int k, int which, double l0) {
  coxparams *param = (coxparams *) ex;
  
  int iter;
  bool fail = false;
  double toler = 1e-12;
  
  NumericVector beta(p), newbeta(p);
  double loglik, newlk;
  NumericVector u(p);
  NumericVector delta(p);
  NumericMatrix imat(p,p);
  NumericMatrix v(p,p);
  
  // --- first step ---
  beta = clone(par);
  
  List der = f_der_2(p, beta, param, firth);
  loglik = der["loglik"];
  u = der["score"];
  imat = as<NumericMatrix>(der["imat"]);
  
  // Lagrange multiplier method as used in SAS PROC LOGISTIC
  v = invsympd(imat, p, toler);
  
  double w = 0.0;
  for (int i = 0; i < p; ++i)
    for (int j = 0; j < p; ++j)
      w -= u[i] * v(i,j) * u[j];
  
  double underroot = -2*(l0 - loglik + 0.5*w)/v(k,k);
  double lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
  u[k] += lambda;
  
  delta.fill(0.0);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j < p; ++j)
      delta[i] += v(i,j) * u[j];
  
  // update beta
  newbeta = beta + delta;
  
  // --- main iteration ---
  for (iter = 0; iter < maxiter; ++iter) {
    der = f_der_2(p, newbeta, param, firth);
    newlk = der["loglik"];
    
    fail = std::isnan(newlk) || std::isinf(newlk);
    if (!fail && fabs(newlk - l0) < eps && w < eps) break;
    
    beta = clone(newbeta);
    loglik = newlk;
    u = as<NumericVector>(der["score"]);
    imat = as<NumericMatrix>(der["imat"]);
    
    // Lagrange multiplier method as used in SAS PROC LOGISTIC
    v = invsympd(imat, p, toler);
    
    double w = 0.0;
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        w -= u[i] * v(i,j) * u[j];
    
    double underroot = -2*(l0 - loglik + 0.5*w)/v(k,k);
    double lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
    u[k] += lambda;
    
    delta.fill(0.0);
    for (int i = 0; i < p; ++i)
      for (int j = 0; j < p; ++j)
        delta[i] += v(i,j) * u[j];
    
    // update beta
    newbeta = beta + delta;
  }
  
  if (iter == maxiter) fail = true;
  if (fail) warning("The algorithm in phregplloop did not converge");
  
  return newbeta[k];
}


// baseline hazard estimates
List f_basehaz(int p, const NumericVector& par, void *ex) {
  coxparams *param = (coxparams *) ex;
  
  const int nused = param->nused;
  const int method = param->method;

  // Precompute exp(eta)
  NumericVector exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    double val = param->offset[person];
    for (int i = 0; i < p; ++i) {
      val += par[i] * param->z(person,i);
    }
    exp_eta[person] = exp(val);
  }
  
  NumericVector a(p);         // s1(beta,k,t)
  NumericVector a2(p);        // sum of w*exp(zbeta)*z for the deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for the deaths
  int natrisk = 0;            // number at risk at this time point
  int ndead = 0;              // number of deaths at this time point
  int ncens = 0;              // number of censored at this time point
  
  // locate the first observation within each stratum
  IntegerVector istratum(1,0);
  for (int i = 1; i < nused; ++i) {
    if (param->strata[i] != param->strata[i-1]) {
      istratum.push_back(i);
    }
  }
  
  int nstrata = static_cast<int>(istratum.size());
  istratum.push_back(nused);
  
  // add time 0 to each stratum
  int J = nstrata;
  for (int i = 0; i < nstrata; ++i) {
    IntegerVector idx = seq(istratum[i], istratum[i+1] - 1);
    NumericVector utime = param->tstop[idx];
    utime = unique(utime);
    J += static_cast<int>(utime.size());
  }
  
  IntegerVector stratum(J);
  NumericVector time(J), nrisk(J), nevent(J), ncensor(J), haz(J), varhaz(J);
  NumericMatrix gradhaz(J,p);
  
  int istrata = param->strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  int j = J;  // index the unique time in ascending order
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      // add time 0 at the start of a new stratum
      j--;
      stratum[j] = istrata;
      time[j] = 0.0;
      nrisk[j] = natrisk;
      
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      natrisk = 0;
      denom = 0.0;
      a.fill(0.0);
    }
    
    const double dtime = param->tstop[person];
    
    // Process all persons tied at this dtime
    bool first = true;
    for (; person < nused && param->tstop[person] == dtime &&
         param->strata[person] == istrata; ++person) {
      
      if (first) { // first incidence at this time
        j--;
        stratum[j] = param->strata[person];
        time[j] = dtime;
        first = false;
      }
      
      const double w = param->weight[person];
      const double r = w * exp_eta[person];
      
      ++natrisk;
      if (param->event[person] == 0) {
        ++ncens;
        denom += r;
        for (int i = 0; i < p; ++i) {
          a[i] += r * param->z(person,i);
        }
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
        for (int i = 0; i < p; ++i) {
          a2[i] += r * param->z(person,i);
        }
      }
    }
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      
      const double r = param->weight[p1] * exp_eta[p1];
      
      natrisk--;
      denom -= r;
      for (int i = 0; i < p; ++i) {
        a[i] -= r * param->z(p1,i);
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
      ndead = 0;
      deadwt = denom2 = 0.0;
      a2.fill(0.0);
    }
  }
  
  // add time 0 for the first stratum
  stratum[0] = istrata;
  time[0] = 0.0;
  nrisk[0] = natrisk;
  
  List result = List::create(
    _["stratum"] = stratum,
    _["time"] = time,
    _["nrisk"] = nrisk,
    _["nevent"] = nevent,
    _["ncensor"] = ncensor,
    _["haz"] = haz,
    _["varhaz"] = varhaz
  );
  
  if (p > 0) result.push_back(gradhaz, "gradhaz");
  
  return result;
}



// martingale residuals
NumericVector f_resmart(int p, const NumericVector& par, void *ex) {
  coxparams *param = (coxparams *) ex;
  
  const int nused = param->nused;
  const int method = param->method;
  const int n = static_cast<int>(param->tstop.size());
  
  // Precompute exp(eta)  
  NumericVector exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    double val = param->offset[person];
    for (int i = 0; i < p; ++i) {
      val += par[i] * param->z(person,i);
    }
    exp_eta[person] = exp(val);
  }

  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  int ndead = 0;              // number of deaths at this time point
  
  // initialize the residuals to the event indicators
  NumericVector resid(n);
  for (int person = 0; person < nused; ++person) {
    resid[person] = param->event[person];
  }
  
  int istrata = param->strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  int j0 = 0; // first person in the stratum
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      j0 = person;
      denom = 0.0;
    }
    
    const double dtime = param->tstop[person];
    
    // process all persons tied at this dtime
    int j1 = person;   // first person in the stratum with the tied time
    for (; person < nused && param->tstop[person] == dtime &&
         param->strata[person] == istrata; ++person) {
      
      const double w = param->weight[person];
      const double r = w * exp_eta[person];
      
      if (param->event[person] == 0) {
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
      const int p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      denom -= param->weight[p1] * exp_eta[p1];
    }
    
    // Add contributions for deaths at this time
    if (ndead > 0) {
      denom += denom2;
      
      for (int j = j0; j <= j2; ++j) {
        if (param->tstart[j] < dtime) {
          double hazard;
          if (method == 0 || ndead == 1) {
            hazard = deadwt / denom;
          } else {
            hazard = 0.0;
            const double meanwt = deadwt / ndead;
            if (j < j1 || param->event[j] == 0) {
              for (int i = 0; i < ndead; ++i) {
                hazard += meanwt /(denom - (i + 0.0)/ndead * denom2);
              }
            } else {
              for (int i = 0; i < ndead; ++i) {
                hazard += (1 - (i + 0.0)/ndead) * meanwt / 
                  (denom - (i + 0.0)/ndead * denom2);
              }
            }
          }
          resid[j] -= hazard * exp_eta[j];
        }
      }
      
      // reset for the next death time
      ndead = 0;
      deadwt = denom2 = 0.0;
    }
  }
  
  return resid;
}


// score residual matrix
NumericMatrix f_ressco_2(int p, const NumericVector& par, void *ex) {
  coxparams *param = (coxparams *) ex;
  
  const int nused = param->nused;
  const int method = param->method;
  const int n = static_cast<int>(param->tstart.size());
  
  // Precompute exp(eta)
  NumericVector exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    double val = param->offset[person];
    for (int i = 0; i < p; ++i) {
      val += par[i] * param->z(person,i);
    }
    exp_eta[person] = exp(val);
  }
  
  NumericMatrix resid(n,p);   // residual matrix
  NumericVector a(p);         // s1(beta,k,t)
  NumericVector a2(p);        // sum of w*exp(zbeta)*z for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  double deadwt = 0.0;        // sum of weights for the deaths
  int ndead = 0;              // number of deaths at this time point
  double cumhaz = 0.0;        // cumulative hazard
  
  NumericVector xhaz(p), mh1(p), mh2(p), mh3(p); // temp vectors
  
  int istrata = param->strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    // Reset when entering a new stratum
    if (param->strata[person] != istrata) {
      istrata = param->strata[person];
      
      // first obs of a new stratum, finish off the prior stratum
      for (; i1 < nused && param->order1[i1] < person; ++i1) {
        const int p1 = param->order1[i1];
        for (int i = 0; i < p; ++i) {
          resid(p1,i) -= exp_eta[p1] * (param->z(p1,i) * cumhaz - xhaz[i]);
        }
      }
      
      denom = 0.0; // reset temporary variables
      cumhaz = 0.0;
      a.fill(0.0); xhaz.fill(0.0);
    }
    
    const double dtime = param->tstop[person];
    
    // process all persons tied at this dtime
    for (; person < nused && param->tstop[person] == dtime &&
         param->strata[person] == istrata; ++person) {
      
      // initialize residuals to score[i] * (x[i] * cumhaz - xhaz), before
      // updating cumhaz and xhaz
      for (int i = 0; i < p; ++i) {
        resid(person,i) = exp_eta[person] * 
          (param->z(person,i) * cumhaz - xhaz[i]);
      }
      
      const double w = param->weight[person];
      const double r = w * exp_eta[person];
      
      if (param->event[person] == 0) {
        denom += r;
        for (int i = 0; i < p; ++i) {
          a[i] += r * param->z(person,i);
        }
      } else {
        ++ndead;
        deadwt += w;
        denom2 += r;
        for (int i = 0; i < p; ++i) {
          a2[i] += r * param->z(person,i);
        }
      }
    }
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      
      const double r = param->weight[p1] * exp_eta[p1];
      denom -= r;
      for (int i = 0; i < p; ++i) {
        // finish the residual by subtracting score[i] * (x[i] * cumhaz - xhaz)
        resid(p1,i) -= exp_eta[p1] * (param->z(p1,i) * cumhaz - xhaz[i]);
        a[i] -= r * param->z(p1,i);
      }
    }
    
    // Add contributions for deaths at this time
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        const double hazard = deadwt/denom;
        cumhaz += hazard;
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          const double xbar = a[i]/denom;
          xhaz[i] += xbar*hazard;
          for (int j = person-1; j >= person - ndead; j--) {
            resid(j,i) += param->z(j,i) - xbar;
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
          const double hazard = meanwt/denom;
          cumhaz += hazard;
          const double downwt = (ndead - k - 1.0)/ndead;
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            const double xbar = a[i] / denom;
            xhaz[i] += xbar * hazard;
            mh1[i]  += hazard * downwt;
            mh2[i]  += xbar * hazard * downwt;
            mh3[i]  += xbar / ndead;
          }
        }
        
        for (int j = person-1; j >= person - ndead; j--) {
          for (int i = 0; i < p; ++i) {
            resid(j,i) += (param->z(j,i) - mh3[i]) +
              exp_eta[j] * (param->z(j,i) * mh1[i] - mh2[i]);
          }
        }
      }
      
      // Reset after processing deaths
      ndead = 0;
      deadwt = denom2 = 0.0;
      a2.fill(0.0);
    }
  }
  
  // finish those remaining in the final stratum
  for (; i1 < nused; ++i1) {
    const int p1 = param->order1[i1];
    for (int i = 0; i < p; ++i)
      resid(p1,i) -= exp_eta[p1] * (param->z(p1,i) * cumhaz - xhaz[i]);
  }
  
  return resid;
}


// schoenfeld residuals
List f_ressch(int p, const NumericVector& par, void *ex) {
  coxparams *param = (coxparams *) ex;

  const int nused = param->nused;
  const int method = param->method;

  // Precompute exp(eta)
  NumericVector exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    double val = param->offset[person];
    for (int i = 0; i < p; ++i) {
      val += par[i] * param->z(person, i);
    }
    exp_eta[person] = exp(val);
  }
  
  NumericVector xbar(p);      // weighted mean covariate at this time
  NumericVector a(p);         // s1(beta,k,t)
  NumericVector a2(p);        // sum of w*exp(zbeta)*z for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for the deaths
  int ndead = 0;              // number of deaths at this time point
  
  int nevent = sum(param->event);   // total number of events
  NumericMatrix resid(nevent, p);   // residual matrix
  IntegerVector index(nevent);      // index of residuals
  
  int istrata = param->strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  int j = nevent;  // index the events in descending order
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      denom = 0.0;
      a.fill(0.0);
    }
    
    const double dtime = param->tstop[person];
    
    // process all persons tied at this dtime
    for (; person < nused && param->tstop[person] == dtime &&
         param->strata[person] == istrata; ++person) {
      
      const double r = param->weight[person] * exp_eta[person];
      
      if (param->event[person] == 0) {
        denom += r;
        for (int i=0; i<p; ++i) {
          a[i] += r * param->z(person,i);
        }
      } else {
        j--;
        resid(j,_) = param->z(person,_);
        index[j] = person;
        
        ++ndead;
        denom2 += r;
        for (int i=0; i<p; ++i) {
          a2[i] += r * param->z(person,i);
        }
      }
    }
    
    // remove subjects no longer at risk
    for (; i1 < nused; ++i1) {
      const int p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      const double r = param->weight[p1] * exp_eta[p1];
      denom -= r;
      for (int i = 0; i < p; ++i) {
        a[i] -= r * param->z(p1,i);
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
        xbar.fill(0.0);
        for (int k = 0; k < ndead; ++k) {
          denom += denom2 / ndead;
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            xbar[i] += a[i] / denom;
          }
        }
        xbar = xbar / ndead;
      }
      
      for (int k = 0; k < ndead; ++k) {
        for (int i = 0; i < p; ++i) {
          resid(j+k,i) -= xbar[i];
        }
      }
      
      // reset for the next death time
      ndead = 0;
      denom2 = 0.0;
      a2.fill(0.0);
    }
  }
  
  return List::create(
    Named("resid") = resid,
    Named("index") = index);
}


// [[Rcpp::export]]
List phregcpp(const DataFrame data,
              const StringVector& rep = "",
              const StringVector& stratum = "",
              const std::string time = "time",
              const std::string time2 = "",
              const std::string event = "event",
              const StringVector& covariates = "",
              const std::string weight = "",
              const std::string offset = "",
              const std::string id = "",
              const std::string ties = "efron",
              const NumericVector& init = NA_REAL,
              const bool robust = false,
              const bool est_basehaz = true,
              const bool est_resid = true,
              const bool firth = false,
              const bool plci = false,
              const double alpha = 0.05,
              const int maxiter = 50,
              const double eps = 1.0e-9) {
  
  int n = data.nrows();
  int p = static_cast<int>(covariates.size());
  if (p == 1 && (covariates[0] == "" || covariates[0] == "none")) {
    p = 0;
  }
  
  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(0);
  } else {
    List out = bygroup(data, rep);
    has_rep = 1;
    repn = out["index"];
    u_rep = DataFrame(out["lookup"]);
  }
  
  bool has_stratum;
  IntegerVector stratumn(n);
  DataFrame u_stratum;
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(0);
  } else {
    List out = bygroup(data, stratum);
    has_stratum = 1;
    stratumn = out["index"];
    u_stratum = DataFrame(out["lookup"]);
  }
  
  bool has_time = hasVariable(data, time);
  if (!has_time) stop("data must contain the time variable");
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each observation");
  }
  
  bool has_time2 = hasVariable(data, time2);
  NumericVector time2n(n);
  if (has_time2) {
    NumericVector time2nz = data[time2];
    time2n = clone(time2nz);
    if (is_true(any(time2n <= timen))) {
      stop("time2 must be greater than time for each observation");
    }
  }
  
  bool has_event = hasVariable(data, event);
  if (!has_event) stop("data must contain the event variable");
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each observation");
  }
  
  NumericMatrix zn(n,p);
  if (p > 0) {
    for (int j=0; j<p; ++j) {
      String zj = covariates[j];
      if (!hasVariable(data, zj)) {
        stop("data must contain the variables in covariates");
      }
      NumericVector u = data[zj];
      for (int i=0; i<n; ++i) {
        zn(i,j) = u[i];
      }
    }
  }
  
  bool has_weight = hasVariable(data, weight);
  NumericVector weightn(n, 1.0);
  if (has_weight) {
    NumericVector weightnz = data[weight];
    weightn = clone(weightnz);
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }
  
  bool has_offset = hasVariable(data, offset);
  NumericVector offsetn(n);
  if (has_offset) {
    NumericVector offsetnz = data[offset];
    offsetn = clone(offsetnz);
  }
  
  // create the numeric id variable
  bool has_id = hasVariable(data, id);
  IntegerVector idn(n);
  if (!has_id) {
    idn = seq(0, n - 1);
  } else {
    SEXP col = data[id];
    SEXPTYPE col_type = TYPEOF(col);
    if (col_type == INTSXP) {
      IntegerVector v = col;
      IntegerVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else if (col_type == REALSXP) {
      NumericVector v = col;
      NumericVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else if (col_type == STRSXP) {
      StringVector v = col;
      StringVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else {
      stop("incorrect type for the id variable in the input data");
    }
  }
  
  if (robust && has_time2 && !has_id) {
    stop("id is needed for counting process data with robust variance");
  }
  
  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  int method = meth == "efron" ? 1 : 0;
  
  // unify right censored data with counting process data
  NumericVector tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  // sort the data by rep
  if (has_rep) {
    IntegerVector order = seq(0, n-1);
    std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
      return repn[i] < repn[j];
    });
    
    repn = repn[order];
    stratumn = stratumn[order];
    tstartn = tstartn[order];
    tstopn = tstopn[order];
    eventn = eventn[order];
    weightn = weightn[order];
    offsetn = offsetn[order];
    idn = idn[order];
    if (p > 0) zn = subset_matrix_by_row(zn, order);
  }
  
  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (int i=0; i<n; ++i) {
    if (repn[i] == NA_INTEGER || stratumn[i] == NA_INTEGER ||
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) || 
        eventn[i] == NA_INTEGER || std::isnan(weightn[i]) || 
        std::isnan(offsetn[i]) || idn[i] == NA_INTEGER) {
      sub[i] = 0;
    }
    for (int j=0; j<p; ++j) {
      if (std::isnan(zn(i,j))) sub[i] = 0;
    }
  }
  
  IntegerVector order = which(sub);
  repn = repn[order];
  stratumn = stratumn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  eventn = eventn[order];
  weightn = weightn[order];
  offsetn = offsetn[order];
  idn = idn[order];
  if (p > 0) zn = subset_matrix_by_row(zn, order);
  n = sum(sub);
  if (n == 0) stop("no observations without missing values");
  
  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (int i=1; i<n; ++i) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);
  
  // variables in the output data sets
  // sumstat data set
  IntegerVector rep01 = seq(0,nreps-1);
  IntegerVector nobs(nreps), nevents(nreps);
  NumericMatrix loglik(nreps,2);
  NumericMatrix regloglik(nreps,2);
  NumericVector scoretest(nreps);
  IntegerVector niter(nreps);
  LogicalVector fails(nreps);
  
  // parest data set
  IntegerVector rep0(nreps*p);
  StringVector par0(nreps*p);
  NumericVector beta0(nreps*p), sebeta0(nreps*p), rsebeta0(nreps*p);
  NumericVector z0(nreps*p), expbeta0(nreps*p);
  NumericMatrix vbeta0(nreps*p,p), rvbeta0(nreps*p,p);
  NumericVector lb0(nreps*p), ub0(nreps*p), prob0(nreps*p);
  StringVector clparm0(nreps*p);
  
  // baseline hazards data set
  int N = 2*n; // account for additional time 0 rows
  IntegerVector drep(N), dstratum(N);
  NumericVector dtime(N), dnrisk(N), dnevent(N), dncensor(N);
  NumericVector dhaz(N), dvarhaz(N);
  NumericMatrix dgradhaz(N,p);

  // martingale residuals
  NumericVector resmart(n);
  
  // linear predictors
  NumericVector linear_predictors(n);
  
  int n0 = 0; // number of rows in the baseline hazard data set
  int bign0 = 0; // number of elements in the martingale residuals vector
  double toler = 1e-12;
  double zcrit = R::qnorm(1-alpha/2,0,1,1,0);
  double xcrit = zcrit * zcrit;
  
  for (int h=0; h<nreps; ++h) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());
    
    IntegerVector stratum1 = stratumn[q1];
    NumericVector tstart1 = tstartn[q1];
    NumericVector tstop1 = tstopn[q1];
    IntegerVector event1 = eventn[q1];
    NumericVector weight1 = weightn[q1];
    NumericVector offset1 = offsetn[q1];
    IntegerVector id1 = idn[q1];
    
    NumericMatrix z1(n1,p);
    if (p > 0) z1 = subset_matrix_by_row(zn, q1);
    
    nobs[h] = n1;
    nevents[h] = sum(event1);
    
    if (nevents[h] == 0) {
      if (p > 0) {
        for (int i=0; i<p; ++i) {
          int k = h*p+i;
          rep0[k] = h;
          par0[k] = covariates[i];
          beta0[k] = NA_REAL;
          sebeta0[k] = 0;
          rsebeta0[k] = 0;
          z0[k] = NA_REAL;
          expbeta0[k] = NA_REAL;
          for (int j=0; j<p; ++j) {
            vbeta0(k,j) = 0;
            rvbeta0(k,j) = 0;
          }
          lb0[k] = NA_REAL;
          ub0[k] = NA_REAL;
          prob0[k] = NA_REAL;
          clparm0[k] = "Wald";
        }
      }
      
      // baseline hazard
      if (est_basehaz) {
        drep[n0] = h;
        dstratum[n0] = 0;
        dtime[n0] = 0;
        dnrisk[n0] = n1;
        dnevent[n0] = 0;
        dncensor[n0] = 0;
        dhaz[n0] = 0;
        dvarhaz[n0] = 0;
        if (p > 0) {
          for (int i=0; i<p; ++i) {
            dgradhaz(n0,i) = 0;
          }
        }
        ++n0;
      }
      
      // martingale residuals
      if (est_resid) {
        for (int i=0; i<n1; ++i) {
          resmart[bign0+i] = 0;
        }
      }
      
      // linear predictors
      for (int i=0; i<n1; ++i) {
        linear_predictors[bign0+i] = offset1[i];
      }
      
      bign0 += n1;
      
      continue;
    }
    
    // sort by stratum
    IntegerVector order0 = seq(0, n1-1);
    std::sort(order0.begin(), order0.end(), [&](int i, int j) {
      return stratum1[i] < stratum1[j];
    });
    
    IntegerVector stratum1z = stratum1[order0];
    NumericVector tstart1z = tstart1[order0];
    NumericVector tstop1z = tstop1[order0];
    IntegerVector event1z = event1[order0];
    
    // locate the first observation within each stratum
    IntegerVector istratum(1,0);
    for (int i=1; i<n1; ++i) {
      if (stratum1z[i] != stratum1z[i-1]) {
        istratum.push_back(i);
      }
    }
    
    int nstrata = static_cast<int>(istratum.size());
    istratum.push_back(n1);
    
    // ignore subjects not at risk for any event time
    IntegerVector ignore1z(n1);
    for (int i=0; i<nstrata; ++i) {
      IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
      NumericVector tstart0 = tstart1z[q0];
      NumericVector tstop0 = tstop1z[q0];
      IntegerVector event0 = event1z[q0];
      NumericVector etime = tstop0[event0==1];
      etime = unique(etime);
      etime.sort();
      IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
      IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
      for (int j=istratum[i]; j<istratum[i+1]; ++j) {
        int j0 = j-istratum[i];
        if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
          ignore1z[j] = 1;
        } else {
          ignore1z[j] = 0;
        }
      }
    }
    
    IntegerVector ignore1(n1); // back to the original order
    for (int i=0; i<n1; ++i) {
      ignore1[order0[i]] = ignore1z[i];
    }
    
    int nused = n1 - sum(ignore1); // number of used observations 
    
    // sort by stopping time in descending order within each stratum
    IntegerVector order2 = seq(0, n1-1);
    std::sort(order2.begin(), order2.end(), [&](int i, int j) {
      if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
      if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
      if (tstop1[i] != tstop1[j]) return tstop1[i] > tstop1[j];
      return event1[i] < event1[j];
    });
    
    IntegerVector stratum1a = stratum1[order2];
    NumericVector tstart1a = tstart1[order2];
    NumericVector tstop1a = tstop1[order2];
    IntegerVector event1a = event1[order2];
    NumericVector weight1a = weight1[order2];
    NumericVector offset1a = offset1[order2];
    IntegerVector id1a = id1[order2];
    IntegerVector ignore1a = ignore1[order2];
    NumericMatrix z1a;
    if (p > 0) z1a = subset_matrix_by_row(z1, order2);
    
    // sort by starting time in descending order within each stratum
    IntegerVector order1a = seq(0, n1-1);
    std::sort(order1a.begin(), order1a.end(), [&](int i, int j) {
      if (ignore1a[i] != ignore1a[j]) return ignore1a[i] < ignore1a[j];
      if (stratum1a[i] != stratum1a[j]) return stratum1a[i] < stratum1a[j];
      return tstart1a[i] > tstart1a[j];
    });
    
    coxparams param = {nused, stratum1a, tstart1a, tstop1a, event1a,
                       weight1a, offset1a, z1a, order1a, method};
    
    NumericVector bint(p);
    List derint = f_der_2(p, bint, &param, firth);
    
    NumericVector b(p);
    NumericMatrix vb(p,p);
    List out;
    
    if (p > 0) {
      IntegerVector colfit = seq(0,p-1);
      if (is_false(any(is_na(init))) && init.size() == p) {
        out = phregloop(p, init, &param, maxiter, eps, firth, colfit, p);
      } else {
        out = phregloop(p, bint, &param, maxiter, eps, firth, colfit, p);
      }

      bool fail = out["fail"];
      if (fail) warning("The algorithm in phregr did not converge");
      
      b = out["coef"];
      vb = as<NumericMatrix>(out["var"]);
      
      NumericVector seb(p);
      for (int j=0; j<p; ++j) {
        seb[j] = sqrt(vb(j,j));
      }
      
      for (int i=0; i<p; ++i) {
        int k = h*p+i;
        rep0[k] = h;
        par0[k] = covariates[i];
        beta0[k] = b[i];
        sebeta0[k] = seb[i];
        for (int j=0; j<p; ++j) {
          vbeta0(k,j) = vb(i,j);
        }
      }

      // score statistic
      
      NumericVector scorebint = firth ? derint["regscore"] : derint["score"];
      NumericMatrix infobint = as<NumericMatrix>(derint["imat"]); 
      NumericMatrix vbint = invsympd(infobint, p, toler);
      for (int i=0; i<p; ++i) {
        for (int j=0; j<p; ++j) {
          scoretest[h] += scorebint[i]*vbint(i,j)*scorebint[j];
        }
      }
      
      niter[h] = out["iter"];
      fails[h] = out["fail"];
      
      // robust variance estimates
      NumericVector rseb(p);  // robust standard error for betahat
      if (robust) {
        NumericMatrix ressco = f_ressco_2(p, b, &param);
        
        int nr; // number of rows in the score residual matrix
        if (!has_id) {
          for (int i=0; i<n1; ++i) {
            for (int j=0; j<p; ++j) {
              ressco(i,j) = weight1a[i]*ressco(i,j);
            }
          }
          nr = n1;
        } else { // need to sum up score residuals by id
          IntegerVector order = seq(0, n1-1);
          std::sort(order.begin(), order.end(), [&](int i, int j) {
            return id1a[i] < id1a[j];
          });
          
          IntegerVector id2 = id1a[order];
          IntegerVector idx(1,0);
          for (int i=1; i<n1; ++i) {
            if (id2[i] != id2[i-1]) {
              idx.push_back(i);
            }
          }
          
          int nids = static_cast<int>(idx.size());
          idx.push_back(n1);
          
          NumericVector weight2 = weight1a[order];
          
          NumericMatrix ressco2(nids,p);
          for (int i=0; i<nids; ++i) {
            for (int j=0; j<p; ++j) {
              for (int k=idx[i]; k<idx[i+1]; ++k) {
                ressco2(i,j) += weight2[k]*ressco(order[k],j);
              }
            }
          }
          
          ressco = ressco2;  // update the score residuals
          nr = nids;
        }
        
        NumericMatrix D(nr,p); // DFBETA
        for (int i=0; i<nr; ++i) {
          for (int j=0; j<p; ++j) {
            for (int k=0; k<p; ++k) {
              D(i,j) += ressco(i,k)*vb(k,j);
            }
          }
        }
        
        NumericMatrix rvb(p,p); // robust variance matrix for betahat
        for (int j=0; j<p; ++j) {
          for (int k=0; k<p; ++k) {
            for (int i=0; i<nr; ++i) {
              rvb(j,k) += D(i,j)*D(i,k);
            }
          }
        }
        
        for (int i=0; i<p; ++i) {
          rseb[i] = sqrt(rvb(i,i));
        }
        
        for (int i=0; i<p; ++i) {
          int k = h*p+i;
          rsebeta0[k] = rseb[i];
          for (int j=0; j<p; ++j) {
            rvbeta0(k,j) = rvb(i,j);
          }
        }
      }
      
      // profile likelihood confidence interval for regression coefficients
      NumericVector lb(p), ub(p), prob(p);
      StringVector clparm(p);
      
      if (plci) {
        double lmax = out["loglik"];
        double l0 = lmax - 0.5*xcrit;
        
        for (int k=0; k<p; ++k) {
          lb[k] = phregplloop(p, b, &param, maxiter, eps, firth, k, -1, l0);
          ub[k] = phregplloop(p, b, &param, maxiter, eps, firth, k, 1, l0);
          
          IntegerVector colfit1(p-1);
          for (int i = 0, j = 0; i < p; ++i) {
            if (i == k) continue;
            colfit1[j++] = i;
          }
          
          NumericVector b0(p);
          List out0 = phregloop(p, b0, &param, maxiter,eps,firth,colfit1,p-1);
          double lmax0 = out0["loglik"];
          prob[k] = R::pchisq(-2*(lmax0 - lmax), 1, 0, 0);
          clparm[k] = "PL";
        }
      } else {
        for (int k=0; k<p; ++k) {
          if (!robust) {
            lb[k] = b[k] - zcrit*seb[k];
            ub[k] = b[k] + zcrit*seb[k];
            prob[k] = R::pchisq(pow(b[k]/seb[k], 2), 1, 0, 0);
          } else {
            lb[k] = b[k] - zcrit*rseb[k];
            ub[k] = b[k] + zcrit*rseb[k];
            prob[k] = R::pchisq(pow(b[k]/rseb[k], 2), 1, 0, 0);
          }
          clparm[k] = "Wald";
        }
      }
      
      for (int i=0; i<p; ++i) {
        int k = h*p+i;
        lb0[k] = lb[i];
        ub0[k] = ub[i];
        prob0[k] = prob[i];
        clparm0[k] = clparm[i];
      }
    }
    
    // log-likelihoods
    if (p > 0 && firth) {
      loglik(h,0) = derint["loglik"];
      loglik(h,1) = out["loglik"];
      regloglik(h,0) = derint["regloglik"];
      regloglik(h,1) = out["regloglik"];
    } else {
      loglik(h,0) = derint["loglik"];
      loglik(h,1) = p > 0 ? out["loglik"] : derint["loglik"];
    }

    // estimate baseline hazard
    if (est_basehaz) {
      // prepare the data for estimating baseline hazards at all time points
      
      // sort by stopping time in descending order within each stratum
      IntegerVector order3 = seq(0, n1-1);
      std::sort(order3.begin(), order3.end(), [&](int i, int j) {
        if (stratum1[i] != stratum1[j]) return stratum1[i] > stratum1[j];
        return tstop1[i] > tstop1[j];
      });
      
      IntegerVector stratum1b = stratum1[order3];
      NumericVector tstart1b = tstart1[order3];
      NumericVector tstop1b = tstop1[order3];
      IntegerVector event1b = event1[order3];
      NumericVector weight1b = weight1[order3];
      NumericVector offset1b = offset1[order3];
      NumericMatrix z1b(n1,p);
      if (p > 0) z1b = subset_matrix_by_row(z1, order3);
      
      // sort by starting time in descending order within each stratum
      IntegerVector order1b = seq(0, n1-1);
      std::sort(order1b.begin(), order1b.end(), [&](int i, int j) {
        if (stratum1b[i] != stratum1b[j]) return stratum1b[i] > stratum1b[j];
        return tstart1b[i] > tstart1b[j];
      });
      
      coxparams paramb = {n1, stratum1b, tstart1b, tstop1b, event1b,
                          weight1b, offset1b, z1b, order1b, method};
      
      List basehaz1 = f_basehaz(p, b, &paramb);
      
      IntegerVector dstratum1 = basehaz1["stratum"];
      NumericVector dtime1 = basehaz1["time"];
      NumericVector dnrisk1 = basehaz1["nrisk"];
      NumericVector dnevent1 = basehaz1["nevent"];
      NumericVector dncensor1 = basehaz1["ncensor"];
      NumericVector dhaz1 = basehaz1["haz"];
      NumericVector dvarhaz1 = basehaz1["varhaz"];
      int J = static_cast<int>(dstratum1.size());
      
      // add to output data frame
      for (int j=0; j<J; ++j) {
        int k = n0 + j;
        drep[k] = h;
        dstratum[k] = dstratum1[j];
        dtime[k] = dtime1[j];
        dnrisk[k] = dnrisk1[j];
        dnevent[k] = dnevent1[j];
        dncensor[k] = dncensor1[j];
        dhaz[k] = dhaz1[j];
        dvarhaz[k] = dvarhaz1[j];
        
        if (p > 0) {
          NumericMatrix dgradhaz1 = basehaz1["gradhaz"];
          for (int i=0; i<p; ++i) {
            dgradhaz(k,i) = dgradhaz1(j,i);
          }
        }
      }
      
      n0 += J;
    }
    
    // martingale residuals
    if (est_resid) {
      NumericVector resid = f_resmart(p, b, &param);
      
      for (int i=0; i<n1; ++i) {
        resmart[bign0 + order2[i]] = resid[i];
      }
    }
    
    // linear predictors
    for (int i=0; i<n1; ++i) {
      linear_predictors[bign0 + order2[i]] = offset1a[i];
      if (p > 0) {
        for (int j=0; j<p; ++j) {
          linear_predictors[bign0 + order2[i]] += b[j]*z1a(i,j);
        }
      }
    }
    
    bign0 += n1;
  }
  
  if (est_basehaz) {
    IntegerVector sub = Range(0, n0 - 1);
    drep = drep[sub];
    dstratum = dstratum[sub];
    dtime = dtime[sub];
    dnrisk = dnrisk[sub];
    dnevent = dnevent[sub];
    dncensor = dncensor[sub];
    dhaz = dhaz[sub];
    dvarhaz = dvarhaz[sub];
    if (p > 0) dgradhaz = subset_matrix_by_row(dgradhaz, sub);
  }

  // prepare the output data sets
  List sumstat = List::create(
    _["n"] = nobs,
    _["nevents"] = nevents,
    _["loglik0"] = loglik(_,0),
    _["loglik1"] = loglik(_,1),
    _["scoretest"] = scoretest,
    _["niter"] = niter,
    _["ties"] = meth,
    _["p"] = p,
    _["robust"] = robust,
    _["firth"] = firth,
    _["fail"] = fails);
  
  if (p > 0 && firth) {
    sumstat.push_back(regloglik(_,0), "loglik0_unpenalized");
    sumstat.push_back(regloglik(_,1), "loglik1_unpenalized");
  }
  
  if (has_rep) {
    for (int i = 0; i < p_rep; ++i) {
      std::string s = as<std::string>(rep[i]);
      SEXP col = u_rep[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == INTSXP) {
        IntegerVector v = col;
        sumstat.push_back(v[rep01], s);
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        sumstat.push_back(v[rep01], s);
      } else if (col_type == STRSXP) {
        StringVector v = col;
        sumstat.push_back(v[rep01], s);
      } else {
        stop("Unsupported type for rep variable" + s);
      }
    }
  }
  
  List result = List::create(
    _["sumstat"] = as<DataFrame>(sumstat)
  );
  
  
  if (p > 0) {
    expbeta0 = exp(beta0);
    if (!robust) z0 = beta0/sebeta0;
    else z0 = beta0/rsebeta0;
    
    List parest = List::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = robust ? rsebeta0 : sebeta0,
      _["z"] = z0,
      _["expbeta"] = expbeta0,
      _["vbeta"] = robust ? rvbeta0 : vbeta0,
      _["lower"] = lb0,
      _["upper"] = ub0,
      _["p"] = prob0,
      _["method"] = clparm0);
    
    if (robust) {
      parest.push_back(sebeta0, "sebeta_naive");
      parest.push_back(vbeta0, "vbeta_naive");
    }
    
    if (has_rep) {
      for (int i=0; i<p_rep; ++i) {
        std::string s = as<std::string>(rep[i]);
        SEXP col = u_rep[s];
        SEXPTYPE col_type = TYPEOF(col);
        if (col_type == INTSXP) {
          IntegerVector v = col;
          parest.push_back(v[rep0], s);
        } else if (col_type == REALSXP) {
          NumericVector v = col;
          parest.push_back(v[rep0], s);
        } else if (col_type == STRSXP) {
          StringVector v = col;
          parest.push_back(v[rep0], s);
        } else {
          stop("Unsupported type for rep variable" + s);
        }
      }
    }
    
    result.push_back(as<DataFrame>(parest), "parest");
  }
  
  
  if (est_basehaz) {
    List basehaz = List::create(
      _["time"] = dtime,
      _["nrisk"] = dnrisk,
      _["nevent"] = dnevent,
      _["ncensor"] = dncensor,
      _["haz"] = dhaz,
      _["varhaz"] = dvarhaz
    );
    
    if (p > 0) {
      basehaz.push_back(dgradhaz, "gradhaz");
    }
    
    if (has_stratum) {
      for (int i = 0; i < p_stratum; ++i) {
        std::string s = as<std::string>(stratum[i]);
        SEXP col = u_stratum[s];
        SEXPTYPE col_type = TYPEOF(col);
        if (col_type == INTSXP) {
          IntegerVector v = col;
          basehaz.push_back(v[dstratum], s);
        } else if (col_type == REALSXP) {
          NumericVector v = col;
          basehaz.push_back(v[dstratum], s);
        } else if (col_type == STRSXP) {
          StringVector v = col;
          basehaz.push_back(v[dstratum], s);
        } else {
          stop("Unsupported type for stratum variable" + s);
        }
      }
    }
    
    if (has_rep) {
      for (int i = 0; i < p_rep; ++i) {
        std::string s = as<std::string>(rep[i]);
        SEXP col = u_rep[s];
        SEXPTYPE col_type = TYPEOF(col);
        if (col_type == INTSXP) {
          IntegerVector v = col;
          basehaz.push_back(v[drep], s);
        } else if (col_type == REALSXP) {
          NumericVector v = col;
          basehaz.push_back(v[drep], s);
        } else if (col_type == STRSXP) {
          StringVector v = col;
          basehaz.push_back(v[drep], s);
        } else {
          stop("Unsupported type for rep variable" + s);
        }
      }
    }
    
    result.push_back(as<DataFrame>(basehaz), "basehaz");
  }
  
  
  if (est_resid) {
    result.push_back(resmart, "residuals");
  }
  
  result.push_back(linear_predictors, "linear_predictors");
  
  return result;
}


// [[Rcpp::export]]
DataFrame survfit_phregcpp(const int p,
                           const NumericVector& beta,
                           const NumericMatrix& vbeta,
                           DataFrame basehaz,
                           DataFrame newdata,
                           const StringVector& covariates = "",
                           const StringVector& stratum = "",
                           const std::string offset = "",
                           const std::string id = "",
                           const std::string tstart = "",
                           const std::string tstop = "",
                           const bool sefit = true,
                           const String conftype = "log-log",
                           const double conflev = 0.95) {
  
  int n0 = basehaz.nrows();
  int n = newdata.nrows();
  int nvar = static_cast<int>(covariates.size());
  
  std::string ct = conftype;
  std::for_each(ct.begin(), ct.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  if (!(ct=="none" || ct=="plain" || ct=="log" || ct=="log-log" ||
      ct=="logit" || ct=="arcsin")) {
    stop("conftype must be none, plain, log, log-log, logit, or arcsin");
  }
  
  if (conflev <= 0 || conflev >= 1) {
    stop("conflev must lie between 0 and 1");
  }
  
  double zcrit = R::qnorm((1+conflev)/2,0,1,1,0);
  
  IntegerVector stratumn0(n0);
  DataFrame u_stratum0;
  IntegerVector nlevels;
  List lookups;
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    stratumn0.fill(0);
  } else {
    List out = bygroup(basehaz, stratum);
    stratumn0 = out["index"];
    u_stratum0 = DataFrame(out["lookup"]);
    nlevels = out["nlevels"];
    lookups = out["lookups"];
  }
  
  bool nullmodel = (p == 0 || (nvar == 1 &&
                    (covariates[0] == "" || covariates[0] == "none")));
  
  if (!nullmodel && nvar != p) {
    stop("incorrect number of covariates for the Cox model");
  }
  
  NumericMatrix zn(n,p);
  for (int j=0; j<p; ++j) {
    String zj = covariates[j];
    if (!hasVariable(newdata, zj)) {
      stop("newdata must contain the variables in covariates");
    }
    NumericVector u = newdata[zj];
    for (int i=0; i<n; ++i) {
      zn(i,j) = u[i];
    }
  }
  
  bool has_stratum;
  IntegerVector stratumn(n);
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = false;
    stratumn.fill(0);
  } else {
    has_stratum = true;
    // match stratum in newdata to stratum in basehaz
    int orep = u_stratum0.nrow();
    for (int i=0; i<p_stratum; ++i) {
      orep /= nlevels[i];
      std::string s = as<std::string>(stratum[i]);
      IntegerVector idx;
      SEXP col = newdata[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == LGLSXP || col_type == INTSXP) {
        IntegerVector v = col;
        IntegerVector w = lookups[i];
        idx = match(v, w) - 1;
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        NumericVector w = lookups[i];
        idx = match(v, w) - 1;
      } else if (col_type == STRSXP) {
        StringVector v = col;
        StringVector w = lookups[i];
        idx = match(v, w) - 1;
      } else {
        stop("Unsupported type for stratum variable: " + s);
      }
      
      stratumn = stratumn + idx * orep;
    }
  }
  
  bool has_offset = hasVariable(newdata, offset);
  NumericVector offsetn(n);
  if (has_offset) {
    NumericVector offsetnz = newdata[offset];
    offsetn = clone(offsetnz);
  }
  
  NumericVector time0 = basehaz["time"];
  NumericVector nrisk0 = basehaz["nrisk"];
  NumericVector nevent0 = basehaz["nevent"];
  NumericVector ncensor0 = basehaz["ncensor"];
  NumericVector haz0 = basehaz["haz"];
  NumericVector vhaz0 = basehaz["varhaz"];
  NumericMatrix ghaz0(n0,p);
  for (int j=0; j<p; ++j) {
    std::string col_name = "gradhaz";
    if (p>1) col_name += "." + std::to_string(j+1);
    NumericVector u = basehaz[col_name];
    ghaz0(_,j) = u;
  }
  
  // create the numeric id variable
  bool has_id = hasVariable(newdata, id);
  IntegerVector idn(n);
  IntegerVector idwi;
  NumericVector idwn;
  StringVector idwc;
  if (!has_id) {
    idn = seq(0,n-1);
  } else { // input data has the counting process style of input
    SEXP col = newdata[id];
    SEXPTYPE col_type = TYPEOF(col);
    if (col_type == INTSXP) {
      IntegerVector idv = col;
      idwi = unique(idv);
      idwi.sort();
      idn = match(idv, idwi) - 1;
    } else if (col_type == REALSXP) {
      NumericVector idv = col;
      idwn = unique(idv);
      idwn.sort();
      idn = match(idv, idwn) - 1;
    } else if (col_type == STRSXP) {
      StringVector idv = col;
      idwc = unique(idv);
      idwc.sort();
      idn = match(idv, idwc) - 1;
    } else {
      stop("incorrect type for the id variable in newdata");
    }
  }
  
  // unify right-censoring data with counting process data
  NumericVector tstartn(n), tstopn(n);
  if (!has_id) { // right-censored data
    tstartn.fill(0.0);
    double maxt0 = max(time0) + 1.0;
    tstopn.fill(maxt0);
  } else {
    bool has_tstart = hasVariable(newdata, tstart);
    if (!has_tstart) stop("newdata must contain the tstart variable");
    NumericVector tstartnz = newdata[tstart];
    tstartn = clone(tstartnz);
    if (is_true(any(tstartn < 0))) {
      stop("tstart must be nonnegative for each observation");
    }
    
    bool has_tstop = hasVariable(newdata, tstop);
    if (!has_tstop) stop("newdata must contain the tstop variable");
    NumericVector tstopnz = newdata[tstop];
    tstopn = clone(tstopnz);
    if (is_true(any(tstopn <= tstartn))) {
      stop("tstop must be greater than tstart for each observation");
    }
  }
  
  // order data by id and tstop, assuming consecutive intervals
  IntegerVector order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    if (idn[i] != idn[j]) return idn[i] < idn[j];
    return tstopn[i] < tstopn[j];
  });
  
  idn = idn[order];
  stratumn = stratumn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  offsetn = offsetn[order];
  if (p > 0) zn = subset_matrix_by_row(zn, order);
  
  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (int i=0; i<n; ++i) {
    if (stratumn[i] == NA_INTEGER || idn[i] == NA_INTEGER || 
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) || 
        std::isnan(offsetn[i])) {
      sub[i] = 0;
    }
    for (int j=0; j<p; ++j) {
      if (std::isnan(zn(i,j))) sub[i] = 0;
    }
  }
  
  order = which(sub);
  stratumn = stratumn[order];
  offsetn = offsetn[order];
  idn = idn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  if (p > 0) zn = subset_matrix_by_row(zn, order);
  n = sum(sub);
  if (n == 0) stop("no observations left after removing missing values");
  
  // risk score
  NumericVector risk(n);
  for (int i=0; i<n; ++i) {
    double val = offsetn[i];
    for (int j=0; j<p; ++j) {
      val += beta[j]*zn(i,j);
    }
    risk[i] = exp(val);
  }
  
  // count number of observations for each id
  IntegerVector idx(1,0);
  for (int i=1; i<n; ++i) {
    if (idn[i] != idn[i-1]) {
      idx.push_back(i);
    }
  }
  
  int nids = static_cast<int>(idx.size());
  idx.push_back(n);
  
  int N = nids*n0; // upper bound on the number of rows in the output
  NumericVector time(N);
  NumericVector nrisk(N), nevent(N), ncensor(N);
  NumericVector cumhaz(N), vcumhaz(N), secumhaz(N);
  IntegerVector strata(N);
  NumericMatrix z(N,p);
  IntegerVector ids(N);
  
  // process by id
  int l = 0;
  for (int h=0; h<nids; ++h) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());
    
    IntegerVector id2 = idn[q1];
    IntegerVector stratum2 = stratumn[q1];
    NumericVector tstart2 = tstartn[q1];
    NumericVector tstop2 = tstopn[q1];
    NumericVector risk2 = risk[q1];
    NumericMatrix z2 = subset_matrix_by_row(zn, q1);
    tstop2.push_front(tstart2[0]);
    
    // match the stratum in basehaz
    IntegerVector idx1 = which(stratumn0 == stratum2[0]);
    NumericVector time01 = time0[idx1];
    
    // left-open and right-closed intervals containing the event time
    IntegerVector idx2 = findInterval3(time01, tstop2, 0, 0, 1);
    IntegerVector sub = which((idx2 >= 1) & (idx2 <= n1));
    int m1 = sub.size();
    
    if (m1 != 0) {
      IntegerVector idx3 = idx1[sub];
      NumericVector time1 = time0[idx3];
      NumericVector nrisk1 = nrisk0[idx3];
      NumericVector nevent1 = nevent0[idx3];
      NumericVector ncensor1 = ncensor0[idx3];
      NumericVector haz1 = haz0[idx3];
      
      IntegerVector idx4 = idx2[sub];
      idx4 = idx4 - 1; // change to 0-1 indexing
      
      // cumulative hazards
      for (int i=0; i<m1; ++i) {
        int r = l + i;
        time[r] = time1[i];
        nrisk[r] = nrisk1[i];
        nevent[r] = nevent1[i];
        ncensor[r] = ncensor1[i];
        
        int k = idx4[i];
        ids[r] = id2[k];
        strata[r] = stratum2[k];
        for (int j=0; j<p; ++j) {
          z(r,j) = z2(k,j);
        }
        
        if (i==0) {
          cumhaz[r] = haz1[i]*risk2[k]; 
        } else {
          cumhaz[r] = cumhaz[r-1] + haz1[i]*risk2[k]; 
        }
      }
      
      if (sefit) {
        NumericVector vhaz1 = vhaz0[idx3];
        NumericMatrix ghaz1(m1,p);
        for (int j=0; j<p; ++j) {
          for (int i=0; i<m1; ++i) {
            ghaz1(i,j) = ghaz0(idx3[i],j);
          }
        }
        
        NumericMatrix a(m1,p);
        for (int j=0; j<p; ++j) {
          for (int i=0; i<m1; ++i) {
            int k = idx4[i];
            if (i==0) {
              a(i,j) = (haz1[i]*z2(k,j) - ghaz1(i,j))*risk2[k]; 
            } else {
              a(i,j) = a(i-1,j) + (haz1[i]*z2(k,j) - ghaz1(i,j))*risk2[k]; 
            }
          }
        }
        
        // calculate the first component of variance
        for (int i=0; i<m1; ++i) {
          int r = l + i;
          int k = idx4[i];
          if (i==0) {
            vcumhaz[r] = vhaz1[i]*risk2[k]*risk2[k]; 
          } else {
            vcumhaz[r] = vcumhaz[r-1] + vhaz1[i]*risk2[k]*risk2[k]; 
          }
        }
        
        // add the second component of variance
        for (int i=0; i<m1; ++i) {
          int r = l + i;
          for (int j=0; j<p; ++j) {
            for (int k=0; k<p; ++k) {
              vcumhaz[r] += a(i,j)*vbeta(j,k)*a(i,k);
            }
          }
          secumhaz[r] = sqrt(vcumhaz[r]);
        }
      }
      
      l += m1;
    }
  }
  
  IntegerVector sub2 = Range(0,l-1);
  time = time[sub2];
  nrisk = nrisk[sub2];
  nevent = nevent[sub2];
  ncensor = ncensor[sub2];
  cumhaz = cumhaz[sub2];
  secumhaz = secumhaz[sub2];
  strata = strata[sub2];
  z = subset_matrix_by_row(z, sub2);
  ids = ids[sub2];
  
  NumericVector surv = exp(-cumhaz);
  
  DataFrame result = DataFrame::create(
    _["time"] = time,
    _["nrisk"] = nrisk,
    _["nevent"] = nevent,
    _["ncensor"] = ncensor,
    _["cumhaz"] = cumhaz,
    _["surv"] = surv);
  
  if (sefit) {
    NumericVector sesurv = surv*secumhaz;
    
    NumericVector lower(l), upper(l);
    for (int i=0; i<l; ++i) {
      NumericVector ci = fsurvci(surv[i], sesurv[i], ct, zcrit);
      lower[i] = ci[0];
      upper[i] = ci[1];
    }
    
    result.push_back(sesurv, "sesurv");
    result.push_back(lower, "lower");
    result.push_back(upper, "upper");
    result.push_back(conflev, "conflev");
    result.push_back(ct, "conftype");
  }
  
  for (int j=0; j<p; ++j) {
    NumericVector u = z(_,j);
    String zj = covariates[j];
    result.push_back(u, zj);
  }
  
  if (has_stratum) {
    for (int i=0; i<p_stratum; ++i) {
      std::string s = as<std::string>(stratum[i]);
      SEXP col = u_stratum0[s];
      SEXPTYPE col_type = TYPEOF(col);
      if (col_type == INTSXP) {
        IntegerVector v = col;
        result.push_back(v[strata], s);
      } else if (col_type == REALSXP) {
        NumericVector v = col;
        result.push_back(v[strata], s);
      } else if (col_type == STRSXP) {
        StringVector v = col;
        result.push_back(v[strata], s);
      } else {
        stop("Unsupported type for stratum variable: " + s);
      }
    }
  }
  
  if (has_id) {
    SEXP col = newdata[id];
    SEXPTYPE col_type = TYPEOF(col);
    if (col_type == INTSXP) {
      result.push_back(idwi[ids], id);
    } else if (col_type == REALSXP) {
      result.push_back(idwn[ids], id);
    } else if (col_type == STRSXP) {
      result.push_back(idwc[ids], id);
    } else {
      stop("incorrect type for the id variable in newdata");
    }
  }
  
  return result;
}


// [[Rcpp::export]]
List residuals_phregcpp(const int p,
                        const NumericVector& beta,
                        const NumericMatrix& vbeta,
                        const NumericVector& resmart,
                        DataFrame data,
                        const StringVector& stratum = "",
                        const std::string time = "time",
                        const std::string time2 = "",
                        const std::string event = "event",
                        const StringVector& covariates = "",
                        const std::string weight = "",
                        const std::string offset = "",
                        const std::string id = "",
                        const std::string ties = "efron",
                        const std::string type = "schoenfeld",
                        const bool collapse = false,
                        const bool weighted = false) {
  
  int n = data.nrows();
  
  IntegerVector stratumn(n);
  DataFrame u_stratum;
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    stratumn.fill(0);
  } else {
    List out = bygroup(data, stratum);
    stratumn = out["index"];
    u_stratum = out["lookup"];
  }
  
  bool has_time = hasVariable(data, time);
  if (!has_time) stop("data must contain the time variable");
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each observation");
  }
  
  bool has_time2 = hasVariable(data, time2);
  NumericVector time2n(n);
  if (has_time2) {
    NumericVector time2nz = data[time2];
    time2n = clone(time2nz);
    if (is_true(any(time2n <= timen))) {
      stop("time2 must be greater than time for each observation");
    }
  }
  
  bool has_event = hasVariable(data, event);
  if (!has_event) stop("data must contain the event variable");
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each observation");
  }
  
  NumericMatrix zn(n,p);
  if (p > 0) {
    for (int j=0; j<p; ++j) {
      String zj = covariates[j];
      if (!hasVariable(data, zj)) {
        stop("data must contain the variables in covariates");
      }
      NumericVector u = data[zj];
      for (int i=0; i<n; ++i) {
        zn(i,j) = u[i];
      }
    }
  }
  
  bool has_weight = hasVariable(data, weight);
  NumericVector weightn(n, 1.0);
  if (has_weight) {
    NumericVector weightnz = data[weight];
    weightn = clone(weightnz);
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }
  
  bool has_offset = hasVariable(data, offset);
  NumericVector offsetn(n);
  if (has_offset) {
    NumericVector offsetnz = data[offset];
    offsetn = clone(offsetnz);
  }
  
  // create the numeric id variable
  bool has_id = hasVariable(data, id);
  IntegerVector idn(n);
  if (!has_id) {
    idn = seq(0, n - 1);
  } else {
    SEXP col = data[id];
    SEXPTYPE col_type = TYPEOF(col);
    if (col_type == INTSXP) {
      IntegerVector v = col;
      IntegerVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else if (col_type == REALSXP) {
      NumericVector v = col;
      NumericVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else if (col_type == STRSXP) {
      StringVector v = col;
      StringVector w = unique(v);
      w.sort();
      idn = match(v, w) - 1;
    } else {
      stop("incorrect type for the id variable in the input data");
    }
  }
  
  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  int method = meth == "efron" ? 1 : 0;
  
  // unify right censored data with counting process data
  NumericVector tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (int i=0; i<n; ++i) {
    if (stratumn[i] == NA_INTEGER || idn[i] == NA_INTEGER || 
        std::isnan(tstartn[i]) || std::isnan(tstopn[i]) || 
        eventn[i] == NA_INTEGER || std::isnan(weightn[i]) || 
        std::isnan(offsetn[i])) {
      sub[i] = 0;
    }
    for (int j=0; j<p; ++j) {
      if (std::isnan(zn(i,j))) sub[i] = 0;
    }
  }
  
  IntegerVector order = which(sub);
  stratumn = stratumn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  eventn = eventn[order];
  weightn = weightn[order];
  offsetn = offsetn[order];
  idn = idn[order];
  if (p > 0) zn = subset_matrix_by_row(zn, order);
  n = sum(sub);
  if (n == 0) stop("no observations left after removing missing values");
  
  // sort by stratum
  IntegerVector order0 = seq(0, n-1);
  std::sort(order0.begin(), order0.end(), [&](int i, int j) {
    return stratumn[i] < stratumn[j];
  });
  
  IntegerVector stratum1z = stratumn[order0];
  NumericVector tstart1z = tstartn[order0];
  NumericVector tstop1z = tstopn[order0];
  IntegerVector event1z = eventn[order0];
  
  // locate the first observation within each stratum
  IntegerVector istratum(1,0);
  for (int i=1; i<n; ++i) {
    if (stratum1z[i] != stratum1z[i-1]) {
      istratum.push_back(i);
    }
  }
  
  int nstrata = static_cast<int>(istratum.size());
  istratum.push_back(n);
  
  // ignore subjects not at risk for any event time
  IntegerVector ignore1z(n);
  for (int i=0; i<nstrata; ++i) {
    IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
    NumericVector tstart0 = tstart1z[q0];
    NumericVector tstop0 = tstop1z[q0];
    IntegerVector event0 = event1z[q0];
    NumericVector etime = tstop0[event0==1];
    etime = unique(etime);
    etime.sort();
    IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
    IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
    for (int j=istratum[i]; j<istratum[i+1]; ++j) {
      int j0 = j - istratum[i];
      if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
        ignore1z[j] = 1;
      } else {
        ignore1z[j] = 0;
      }
    }
  }
  
  IntegerVector ignore(n);
  for (int i=0; i<n; ++i) {
    ignore[order0[i]] = ignore1z[i];
  }
  
  int nused = n - sum(ignore);
  
  order = seq(0, n-1);
  IntegerVector idx(1,0);
  int nids = n;
  if (has_id) { // collapse over id
    std::sort(order.begin(), order.end(), [&](int i, int j) {
      return idn[i] < idn[j];
    });
    
    IntegerVector id2 = idn[order];
    for (int i=1; i<n; ++i) {
      if (id2[i] != id2[i-1]) {
        idx.push_back(i);
      }
    }
    
    nids = static_cast<int>(idx.size());
    idx.push_back(n);
  }
  
  
  List result;
  if (type == "martingale") {
    NumericVector rr = clone(resmart);
    if (weighted) rr = rr*weightn;
    
    if (collapse) { // collapse over id
      NumericVector rr2(nids);
      for (int i=0; i<nids; ++i) {
        for (int j=idx[i]; j<idx[i+1]; ++j) {
          rr2[i] += rr[order[j]];
        }
      }
      
      rr = rr2;
    }
    
    result = List::create(Named("resid") = rr);
  } else if (type == "deviance") {
    NumericVector rr = clone(resmart);
    IntegerVector status = clone(eventn);
    int m = n;
    
    if (weighted) rr = rr*weightn;
    
    if (collapse) { // collapse over id
      NumericVector rr2(nids);
      IntegerVector status2(nids);
      for (int i=0; i<nids; ++i) {
        for (int j=idx[i]; j<idx[i+1]; ++j) {
          int k = order[j];
          rr2[i] += rr[k];
          status2[i] += eventn[k];
        }
      }
      
      rr = rr2;
      status = status2;
      m = nids;
    }
    
    for (int i=0; i<m; ++i) {
      double temp = status[i] == 0 ? 0 : status[i]*log(status[i] - rr[i]);
      rr[i] = ((rr[i]>0) - (rr[i]<0))*sqrt(-2*(rr[i] + temp));
    }
    
    result = List::create(Named("resid") = rr);
  } else if (p == 0) {
    stop("covariates must be present for score and schoenfeld residuals");
  } else {
    NumericMatrix rr(n,p);
    if (type == "score" || type == "dfbeta" || type == "dfbetas") {
      // sort by stopping time in descending order within each stratum
      IntegerVector order2 = seq(0, n-1);
      std::sort(order2.begin(), order2.end(), [&](int i, int j) {
        if (ignore[i] != ignore[j]) return ignore[i] < ignore[j];
        if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
        if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
        return eventn[i] < eventn[j];
      });
      
      IntegerVector stratum1 = stratumn[order2];
      NumericVector tstart1 = tstartn[order2];
      NumericVector tstop1 = tstopn[order2];
      IntegerVector event1 = eventn[order2];
      NumericVector weight1 = weightn[order2];
      NumericVector offset1 = offsetn[order2];
      IntegerVector ignore1 = ignore[order2];
      NumericMatrix z1 = subset_matrix_by_row(zn, order2);
      
      // sort by starting time in descending order within each stratum
      IntegerVector order1 = seq(0, n-1);
      std::sort(order1.begin(), order1.end(), [&](int i, int j) {
        if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
        if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
        return tstart1[i] > tstart1[j];
      });
      
      coxparams param = {nused, stratum1, tstart1, tstop1, event1,
                         weight1, offset1, z1, order1, method};
      
      NumericMatrix ressco = f_ressco_2(p, beta, &param);
      
      NumericMatrix score(n,p);
      for (int i=0; i<n; ++i) {
        score(order2[i],_) = ressco(i,_); // original order
      }
      
      if (type == "dfbeta" || type == "dfbetas") {
        for (int i=0; i<n; ++i) {
          for (int k=0; k<p; ++k) {
            for (int j=0; j<p; ++j) {
              rr(i,k) += score(i,j)*vbeta(j,k);
            }
            if (type == "dfbetas") {
              rr(i,k) /= sqrt(vbeta(k,k));
            }
          }
        }
      } else {
        rr = score;
      }
      
      if (weighted) {
        for (int i=0; i<n; ++i) {
          for (int k=0; k<p; ++k) {
            rr(i,k) = rr(i,k)*weightn[i];
          }
        }
      }
      
      if (collapse) { // collapse over id
        NumericMatrix rr2(nids,p);
        for (int i=0; i<nids; ++i) {
          for (int k=0; k<p; ++k) {
            for (int j=idx[i]; j<idx[i+1]; ++j) {
              rr2(i,k) += rr(order[j],k);
            }
          }
        }
        
        rr = rr2;
      }
      
      result = List::create(Named("resid") = rr);
    } else if (type == "schoenfeld" || type == "scaledsch") {
      // sort by stopping time in descending order within each stratum
      IntegerVector order2 = seq(0, n-1);
      std::sort(order2.begin(), order2.end(), [&](int i, int j) {
        if (ignore[i] != ignore[j]) return ignore[i] < ignore[j];
        if (stratumn[i] != stratumn[j]) return stratumn[i] > stratumn[j];
        if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
        if (eventn[i] != eventn[j]) return eventn[i] < eventn[j];
        return idn[i] > idn[j];
      });
      
      IntegerVector stratum1 = stratumn[order2];
      NumericVector tstart1 = tstartn[order2];
      NumericVector tstop1 = tstopn[order2];
      IntegerVector event1 = eventn[order2];
      NumericVector weight1 = weightn[order2];
      NumericVector offset1 = offsetn[order2];
      IntegerVector id1 = idn[order2];
      IntegerVector ignore1 = ignore[order2];
      NumericMatrix z1 = subset_matrix_by_row(zn, order2);
      
      // sort by starting time in descending order within each stratum
      IntegerVector order1 = seq(0, n-1);
      std::sort(order1.begin(), order1.end(), [&](int i, int j) {
        if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
        if (stratum1[i] != stratum1[j]) return stratum1[i] > stratum1[j];
        return tstart1[i] > tstart1[j];
      });
      
      coxparams param = {nused, stratum1, tstart1, tstop1, event1,
                         weight1, offset1, z1, order1, method};
      
      List out = f_ressch(p, beta, &param);
      rr = as<NumericMatrix>(out["resid"]);
      IntegerVector index = out["index"];
      IntegerVector stratum2 = stratum1[index];
      NumericVector time3 = tstop1[index];
      int ndead = static_cast<int>(index.size());
      
      if (weighted) {
        for (int i=0; i<ndead; ++i) {
          for (int k=0; k<p; ++k) {
            rr(i,k) = rr(i,k)*weightn[i];
          }
        }
      }
      
      if (type == "scaledsch") {
        NumericMatrix rr2(ndead,p);
        for (int i=0; i<ndead; ++i) {
          for (int k=0; k<p; ++k) {
            for (int j=0; j<p; ++j) {
              rr2(i,k) += rr(i,j)*vbeta(j,k);
            }
            rr2(i,k) = rr2(i,k)*ndead + beta[k];
          }
        }
        rr = rr2;
      }
      
      result = List::create(
        Named("resid") = rr,
        Named("time") = time3);
      
      IntegerVector stratum3 = unique(stratum2);
      if (stratum3.size() > 1) {
        List strata(p_stratum);
        for (int i = 0; i < p_stratum; ++i) {
          std::string s = as<std::string>(stratum[i]);
          SEXP col = u_stratum[s];
          SEXPTYPE col_type = TYPEOF(col);
          if (col_type == INTSXP) {
            IntegerVector v = col;
            strata[i] = v[stratum2];
          } else if (col_type == REALSXP) {
            NumericVector v = col;
            strata[i] = v[stratum2];
          } else if (col_type == STRSXP) {
            StringVector v = col;
            strata[i] = v[stratum2];
          } else {
            stop("Unsupported type for stratum variable: " + s);
          }
        }
        strata.attr("names") = stratum;
        result.push_back(as<DataFrame>(strata), "strata");
      }
    } else {
      stop("unknown type of residuals");
    }
  }
  
  return result;
}


// function for individual contributions to score and information matrix
// for the Cox model
List f_der_i_2(int p, const NumericVector& par, void* ex) {
  coxparams* param = (coxparams*) ex;
  
  const int nused = param->nused;
  const int method = param->method;
  
  // Precompute eta and exp(eta)
  NumericVector eta(nused);
  NumericVector exp_eta(nused);
  for (int person = 0; person < nused; ++person) {
    double val = param->offset[person];
    for (int i = 0; i < p; ++i) {
      val += par[i] * param->z(person,i);
    }
    eta[person] = val;
    exp_eta[person] = exp(val);
  }
  
  NumericMatrix u(nused,p);         // score vector for each individual
  NumericMatrix imat(nused,p*p);    // information matrix for each individual
  NumericVector a(p);         // s1(beta,k,t)
  NumericVector a2(p);        // sum of w*exp(zbeta)*z for the deaths
  NumericMatrix cmat(p,p);    // s2(beta,k,t)
  NumericMatrix cmat2(p,p);   // sum of w*exp(zbeta)*z*z' for the deaths
  double denom = 0.0;         // s0(beta,k,t)
  double denom2 = 0.0;        // sum of weighted risks for deaths
  int ndead = 0;              // number of deaths at this time point
  
  
  int istrata = param->strata[0];
  int i1 = 0; // index for removing out-of-risk subjects
  
  // Loop through subjects
  for (int person = 0; person < nused; ) {
    // Reset when entering a new stratum
    if (param->strata[person] != istrata) {
      istrata = param->strata[person];
      i1 = person;
      denom = 0.0;
      
      for (int i = 0; i < p; ++i) {
        a[i] = 0.0;
        for (int j = 0; j <= i; ++j) {
          cmat(i,j) = 0.0;
        }
      }
    }
    
    const double dtime = param->tstop[person];
    
    // Process all persons tied at this dtime
    int person1 = person;
    
    for (; person < nused && param->tstop[person] == dtime &&
         param->strata[person] == istrata; ++person) {
      
      const double w = param->weight[person];
      const double r = w * exp_eta[person];
      
      if (param->event[person] == 0) {
        denom += r;
        
        for (int i = 0; i < p; ++i) {
          const double zi = param->z(person,i);
          a[i] += r * zi;
          for (int j = 0; j <= i; ++j) {
            const double zj = param->z(person,j);
            cmat(i,j) += r * zi * zj;
          }
        }
      } else {
        ++ndead;
        denom2 += r;
        for (int i = 0; i < p; ++i) {
          const double zi = param->z(person,i);
          a2[i] += r * zi;
          u(person,i) = w * zi;
          for (int j = 0; j <= i; ++j) {
            const double zj = param->z(person,j);
            cmat2(i,j) += r * zi * zj;
          }
        }
      }
    }
    
    // Remove subjects leaving risk set
    for (; i1 < nused; ++i1) {
      const int p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      
      const double r = param->weight[p1] * exp_eta[p1];
      denom -= r;
      for (int i = 0; i < p; ++i) {
        const double zi = param->z(p1,i);
        a[i] -= r * zi;
        for (int j = 0; j <= i; ++j) {
          const double zj = param->z(p1,j);
          cmat(i,j) -= r * zi * zj;
        }
      }
    }
    
    // Add contributions for deaths at this time
    if (ndead > 0) {
      if (method == 0 || ndead == 1) { // Breslow or single event
        denom += denom2;
        for (int i = 0; i < p; ++i) {
          a[i] += a2[i];
          const double xbar = a[i] / denom;
          for (int human = person1; human < person; ++human) {
            if (param->event[human] == 1) 
              u(human,i) -= param->weight[human] * xbar;
          }

          for (int j = 0; j <= i; ++j) {
            cmat(i,j) += cmat2(i,j);
            
            for (int human = person1; human < person; ++human) {
              if (param->event[human] == 1) 
                imat(human, i*p+j) = param->weight[human] * 
                  (cmat(i,j) - xbar * a[j]) / denom;
            }
          }
        }
      } else { // Efron method
        const double increment = denom2 / ndead;
        for (int l = 0; l < ndead; ++l) {
          denom += increment;
          for (int i = 0; i < p; ++i) {
            a[i] += a2[i] / ndead;
            const double xbar = a[i] / denom;
            for (int human = person1; human < person; ++human) {
              if (param->event[human] == 1)
                u(human,i) -= param->weight[human]/ndead * xbar;
            }
            
            for (int j = 0; j <= i; ++j) {
              cmat(i,j) += cmat2(i,j) / ndead;
              for (int human = person1; human < person; ++human) {
                if (param->event[human] == 1)
                  imat(human,i*p+j) +=  param->weight[human]/ndead * 
                    (cmat(i,j) - xbar * a[j]) / denom;
              }
            }
          }
        }
      }
      
      // Reset after processing deaths
      ndead = 0;
      denom2 = 0.0;
      for (int i = 0; i < p; ++i) {
        a2[i] = 0;
        for (int j = 0; j <= i; ++j) {
          cmat2(i,j) = 0;
        }
      }
    }
  }
  
  // fill the symmetric elements of the information matrix
  for (int person = 0; person < nused; ++person) 
    for (int i = 0; i < p - 1; ++i)
      for (int j = i+1; j < p; ++j)
        imat(person,i*p+j) = imat(person,j*p+i);

  List result = List::create(
    _["score_i"] = u,
    _["imat_i"] = imat
  );

  return result;
}


// [[Rcpp::export]]
List assess_phregcpp(const int p,
                     const NumericVector& beta,
                     const NumericMatrix& vbeta,
                     DataFrame data,
                     const StringVector& stratum = "",
                     const std::string time = "time",
                     const std::string time2 = "",
                     const std::string event = "event",
                     const StringVector& covariates = "",
                     const std::string weight = "",
                     const std::string offset = "",
                     const std::string ties = "efron",
                     const int resample = 1000,
                     const int seed = NA_INTEGER) {
  
  if (seed != NA_INTEGER) set_seed(seed);
  
  if (p <= 0) {
    stop("covariates must be present to test proportional hazards");
  }
  
  int n = data.nrows();
  
  IntegerVector stratumn(n);
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    stratumn.fill(0);
  } else {
    List out = bygroup(data, stratum);
    stratumn = out["index"];
  }
  
  bool has_time = hasVariable(data, time);
  if (!has_time) stop("data must contain the time variable");
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each observation");
  }
  
  bool has_time2 = hasVariable(data, time2);
  NumericVector time2n(n);
  if (has_time2) {
    NumericVector time2nz = data[time2];
    time2n = clone(time2nz);
    if (is_true(any(time2n <= timen))) {
      stop("time2 must be greater than time for each observation");
    }
  }
  
  bool has_event = hasVariable(data, event);
  if (!has_event) stop("data must contain the event variable");
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each observation");
  }
  
  NumericMatrix zn(n,p);
  if (p > 0) {
    for (int j=0; j<p; ++j) {
      String zj = covariates[j];
      if (!hasVariable(data, zj)) {
        stop("data must contain the variables in covariates");
      }
      NumericVector u = data[zj];
      for (int i=0; i<n; ++i) {
        zn(i,j) = u[i];
      }
    }
  }
  
  bool has_weight = hasVariable(data, weight);
  NumericVector weightn(n, 1.0);
  if (has_weight) {
    NumericVector weightnz = data[weight];
    weightn = clone(weightnz);
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }
  
  bool has_offset = hasVariable(data, offset);
  NumericVector offsetn(n);
  if (has_offset) {
    NumericVector offsetnz = data[offset];
    offsetn = clone(offsetnz);
  }
  
  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  int method = meth == "efron" ? 1 : 0;
  
  // unify right censored data with counting process data
  NumericVector tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (int i=0; i<n; ++i) {
    if (stratumn[i] == NA_INTEGER || std::isnan(tstartn[i]) || 
        std::isnan(tstopn[i]) || eventn[i] == NA_INTEGER || 
        std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
      sub[i] = 0;
    }
    for (int j=0; j<p; ++j) {
      if (std::isnan(zn(i,j))) sub[i] = 0;
    }
  }
  
  IntegerVector order = which(sub);
  stratumn = stratumn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  eventn = eventn[order];
  weightn = weightn[order];
  offsetn = offsetn[order];
  zn = subset_matrix_by_row(zn, order);
  n = sum(sub);
  if (n == 0) stop("no observations left after removing missing values");
  
  // sort by stratum
  IntegerVector order0 = seq(0, n-1);
  std::sort(order0.begin(), order0.end(), [&](int i, int j) {
    return stratumn[i] < stratumn[j];
  });
  
  IntegerVector stratum1z = stratumn[order0];
  NumericVector tstart1z = tstartn[order0];
  NumericVector tstop1z = tstopn[order0];
  IntegerVector event1z = eventn[order0];
  
  // locate the first observation within each stratum
  IntegerVector istratum(1,0);
  for (int i=1; i<n; ++i) {
    if (stratum1z[i] != stratum1z[i-1]) {
      istratum.push_back(i);
    }
  }
  
  int nstrata = static_cast<int>(istratum.size());
  istratum.push_back(n);
  
  // ignore subjects not at risk for any event time
  IntegerVector ignore1z(n);
  for (int i=0; i<nstrata; ++i) {
    IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
    NumericVector tstart0 = tstart1z[q0];
    NumericVector tstop0 = tstop1z[q0];
    IntegerVector event0 = event1z[q0];
    NumericVector etime = tstop0[event0==1];
    etime = unique(etime);
    etime.sort();
    IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
    IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
    for (int j=istratum[i]; j<istratum[i+1]; ++j) {
      int j0 = j - istratum[i];
      if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
        ignore1z[j] = 1;
      } else {
        ignore1z[j] = 0;
      }
    }
  }
  
  IntegerVector ignoren(n); // final ignore vector
  for (int i=0; i<n; ++i) {
    ignoren[order0[i]] = ignore1z[i];
  }
  
  int nused = n - sum(ignoren);
  
  // sort by stopping time in descending order within each stratum
  IntegerVector order1 = seq(0, n-1);
  std::sort(order1.begin(), order1.end(), [&](int i, int j) {
    if (ignoren[i] != ignoren[j]) return ignoren[i] < ignoren[j];
    if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
    if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
    return eventn[i] < eventn[j];
  });
  
  IntegerVector stratum1 = stratumn[order1];
  NumericVector tstart1 = tstartn[order1];
  NumericVector tstop1 = tstopn[order1];
  IntegerVector event1 = eventn[order1];
  NumericVector weight1 = weightn[order1];
  NumericVector offset1 = offsetn[order1];
  IntegerVector ignore1 = ignoren[order1];
  NumericMatrix z1 = subset_matrix_by_row(zn, order1);
  
  // sort by starting time in descending order within each stratum
  IntegerVector order10 = seq(0, n-1);
  std::sort(order10.begin(), order10.end(), [&](int i, int j) {
    if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
    if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
    return tstart1[i] > tstart1[j];
  });
  
  coxparams param = {nused, stratum1, tstart1, tstop1, event1,
                     weight1, offset1, z1, order10, method};
  
  List out = f_der_i_2(p, beta, &param);
  NumericMatrix score_i = as<NumericMatrix>(out["score_i"]);
  NumericMatrix imat_i = as<NumericMatrix>(out["imat_i"]);
  
  // order data by ascending tstop
  IntegerVector order2 = seq(0, n-1);
  std::sort(order2.begin(), order2.end(), [&](int i, int j) {
    return tstop1[i] < tstop1[j];
  });
  
  NumericVector tstop2 = tstop1[order2];
  IntegerVector event2 = event1[order2];
  NumericMatrix score_i2 = subset_matrix_by_row(score_i, order2);
  NumericMatrix imat_i2 = subset_matrix_by_row(imat_i, order2);
  
  // only consider event times
  IntegerVector q = which(event2 == 1);
  NumericVector tstop3 = tstop2[q];
  NumericMatrix score_i3 = subset_matrix_by_row(score_i2, q);
  NumericMatrix imat_i3 = subset_matrix_by_row(imat_i2, q);
  int nq = static_cast<int>(q.size());
  
  // identify unique event times
  NumericVector t = unique(tstop3);
  t.sort();
  
  int nt = static_cast<int>(t.size());
  NumericMatrix score_i4(nt, p);
  NumericMatrix imat_i4(nt, p*p);
  
  tstop3.push_front(-1.0); // to facilitate the loop
  int k = -1;
  for (int i=1; i<=nq; ++i) {
    if (tstop3[i] != tstop3[i-1]) k++;
    
    for (int j=0; j<p; ++j) {
      score_i4(k,j) += score_i3(i-1,j);
      for (int l=0; l<p; ++l) {
        imat_i4(k,j*p+l) += imat_i3(i-1,j*p+l);
      }
    }
  }
  
  // cumulative sums of score and information processes
  for (int i=1; i<nt; ++i) {
    for (int j=0; j<p; ++j) {
      score_i4(i,j) += score_i4(i-1,j);
      for (int l=0; l<p; ++l) {
        imat_i4(i,j*p+l) += imat_i4(i-1,j*p+l);
      }
    }
  }
  
  // standardize the score processes
  for (int i=0; i<nt; ++i) {
    for (int j=0; j<p; ++j) {
      score_i4(i,j) *= sqrt(vbeta(j,j));
    }
  }
  
  // resampling to approximate the null distribution
  List G_list(resample);
  for (int r=0; r<resample; ++r) {
    NumericVector G(nq);
    for (int i=0; i<nq; ++i) {
      G[i] = R::rnorm(0.0, 1.0);
    }
    G_list[r] = G;
  }

  // sum(u*G)
  List U_list(resample);
  for (int r=0; r<resample; ++r) {
    NumericVector U(p);
    NumericVector G = G_list[r];
    for (int i=0; i<nq; ++i) {
      for (int j=0; j<p; ++j) {
        U[j] += score_i3(i,j)*G[i];
      }
    }
    U_list[r] = U;
  }
 
  // sum(I(X <= t)*U*G)  
  List score_i_list(resample);
  for (int r=0; r<resample; ++r) {
    NumericMatrix score_i4(nt, p);
    NumericVector G = G_list[r];
    
    int k = -1;
    for (int i=1; i<=nq; ++i) {
      if (tstop3[i] != tstop3[i-1]) k++;
      
      for (int j=0; j<p; ++j) {
        score_i4(k,j) += score_i3(i-1,j)*G[i-1];
      }
    }
    
    // cumulative sums of score processes
    for (int i=1; i<nt; ++i) {
      for (int j=0; j<p; ++j) {
        score_i4(i,j) += score_i4(i-1,j);
      }
    }
    
    score_i_list[r] = score_i4;
  }
  
  // piece together the result
  List score_i_list_2(resample);
  for (int r=0; r<resample; ++r) {
    NumericMatrix score_i4 = score_i_list[r];
    NumericMatrix score_i4_2(nt, p);
    NumericVector U = U_list[r];
    for (int i=0; i<nt; ++i) {
      for (int j=0; j<p; ++j) {
        score_i4_2(i,j) = score_i4(i,j);
        for (int l=0; l<p; ++l) {
          for (int s=0; s<p; ++s) {
            score_i4_2(i,j) -= imat_i4(i,j*p+l)*vbeta(l,s)*U[s];  
          }
        }
      }
    }
    score_i_list_2[r] = score_i4_2;
  }
  
  // standardize the resampled score processes
  for (int r=0; r<resample; ++r) {
    NumericMatrix score_i4 = score_i_list_2[r];
    for (int i=0; i<nt; ++i) {
      for (int j=0; j<p; ++j) {
        score_i4(i,j) *= sqrt(vbeta(j,j));
      }
    }
    score_i_list_2[r] = score_i4;
  }
  
  // calculate individual p-values for proportional hazards assumption test
  NumericVector p_values(p+1);
  NumericVector maxabs(p+1);
  for (int j=0; j<p; ++j) {
    NumericVector obs(nt);
    for (int i=0; i<nt; ++i) {
      obs[i] = fabs(score_i4(i,j));
    }
    double tobs = max(obs);
    maxabs[j] = tobs;
    
    int count = 0;
    for (int r=0; r<resample; ++r) {
      NumericMatrix score_i4 = score_i_list_2[r];
      NumericVector sim(nt);
      for (int i=0; i<nt; ++i) {
        sim[i] = fabs(score_i4(i,j));
      }
      double tsim = max(sim);
      if (tsim >= tobs) count++;
    }
    p_values[j] = static_cast<double>(count) / static_cast<double>(resample);
  }
  
  // p-value for the global test
  NumericVector obs(nt);
  for (int i=0; i<nt; ++i) {
    for (int j=0; j<p; ++j) {
      obs[i] += fabs(score_i4(i,j));
    }
  }
  double tobs = max(obs);
  
  int count = 0;
  for (int r=0; r<resample; ++r) {
    NumericMatrix score_i4 = score_i_list_2[r];
    NumericVector sim(nt);
    for (int i=0; i<nt; ++i) {
      for (int j=0; j<p; ++j) {
        sim[i] += fabs(score_i4(i,j));
      }
    }
    double tsim = max(sim);
    if (tsim >= tobs) count++;
  }
  p_values[p] = static_cast<double>(count) / static_cast<double>(resample);
  
  maxabs[p] = tobs;
  List result = List::create(
    Named("time") = t,
    Named("score_t") = score_i4,
    Named("score_t_list") = score_i_list_2,
    Named("max_abs_value") = maxabs,
    Named("p_value") = p_values
  );
  
  return result;
}


// [[Rcpp::export]]
List zph_phregcpp(const int p,
                  const NumericVector& beta,
                  const NumericMatrix& vbeta,
                  const NumericVector& resmart,
                  DataFrame data,
                  const StringVector& stratum = "",
                  const std::string time = "time",
                  const std::string time2 = "",
                  const std::string event = "event",
                  const StringVector& covariates = "",
                  const std::string weight = "",
                  const std::string offset = "",
                  const std::string ties = "efron",
                  const std::string transform= "km") {
  
  if (p <= 0) {
    stop("covariates must be present to test proportional hazards");
  }
  
  int n = data.nrows();
  
  bool has_stratum;
  IntegerVector stratumn(n);
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = false;
    stratumn.fill(0);
  } else {
    has_stratum = true;
    List out = bygroup(data, stratum);
    stratumn = out["index"];
  }
  
  bool has_time = hasVariable(data, time);
  if (!has_time) stop("data must contain the time variable");
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each observation");
  }
  
  bool has_time2 = hasVariable(data, time2);
  NumericVector time2n(n);
  if (has_time2) {
    NumericVector time2nz = data[time2];
    time2n = clone(time2nz);
    if (is_true(any(time2n <= timen))) {
      stop("time2 must be greater than time for each observation");
    }
  }
  
  bool has_event = hasVariable(data, event);
  if (!has_event) stop("data must contain the event variable");
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each observation");
  }
  
  NumericMatrix zn(n,p);
  if (p > 0) {
    for (int j=0; j<p; ++j) {
      String zj = covariates[j];
      if (!hasVariable(data, zj)) {
        stop("data must contain the variables in covariates");
      }
      NumericVector u = data[zj];
      for (int i=0; i<n; ++i) {
        zn(i,j) = u[i];
      }
    }
  }
  
  bool has_weight = hasVariable(data, weight);
  NumericVector weightn(n, 1.0);
  if (has_weight) {
    NumericVector weightnz = data[weight];
    weightn = clone(weightnz);
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }
  
  bool has_offset = hasVariable(data, offset);
  NumericVector offsetn(n);
  if (has_offset) {
    NumericVector offsetnz = data[offset];
    offsetn = clone(offsetnz);
  }
  
  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  int method = meth == "efron" ? 1 : 0;
  
  // unify right censored data with counting process data
  NumericVector tstartn(n), tstopn(n);
  if (!has_time2) {
    tstopn = timen;
  } else {
    tstartn = timen;
    tstopn = time2n;
  }
  
  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (int i=0; i<n; ++i) {
    if (stratumn[i] == NA_INTEGER || std::isnan(tstartn[i]) || 
        std::isnan(tstopn[i]) || eventn[i] == NA_INTEGER || 
        std::isnan(weightn[i]) || std::isnan(offsetn[i])) {
      sub[i] = 0;
    }
    for (int j=0; j<p; ++j) {
      if (std::isnan(zn(i,j))) sub[i] = 0;
    }
  }
  
  IntegerVector order = which(sub);
  stratumn = stratumn[order];
  tstartn = tstartn[order];
  tstopn = tstopn[order];
  eventn = eventn[order];
  weightn = weightn[order];
  offsetn = offsetn[order];
  zn = subset_matrix_by_row(zn, order);
  n = sum(sub);
  if (n == 0) stop("no observations left after removing missing values");
  
  // sort by stratum
  IntegerVector order0 = seq(0, n-1);
  std::sort(order0.begin(), order0.end(), [&](int i, int j) {
    return stratumn[i] < stratumn[j];
  });
  
  IntegerVector stratum1z = stratumn[order0];
  NumericVector tstart1z = tstartn[order0];
  NumericVector tstop1z = tstopn[order0];
  IntegerVector event1z = eventn[order0];
  
  // locate the first observation within each stratum
  IntegerVector istratum(1,0);
  for (int i=1; i<n; ++i) {
    if (stratum1z[i] != stratum1z[i-1]) {
      istratum.push_back(i);
    }
  }
  
  int nstrata = static_cast<int>(istratum.size());
  istratum.push_back(n);
  
  // ignore subjects not at risk for any event time
  IntegerVector ignore1z(n);
  for (int i=0; i<nstrata; ++i) {
    IntegerVector q0 = Range(istratum[i], istratum[i+1]-1);
    NumericVector tstart0 = tstart1z[q0];
    NumericVector tstop0 = tstop1z[q0];
    IntegerVector event0 = event1z[q0];
    NumericVector etime = tstop0[event0==1];
    etime = unique(etime);
    etime.sort();
    IntegerVector index1 = findInterval3(tstart0, etime, 0, 0, 0);
    IntegerVector index2 = findInterval3(tstop0, etime, 0, 0, 0);
    for (int j=istratum[i]; j<istratum[i+1]; ++j) {
      int j0 = j - istratum[i];
      if (index1[j0] == index2[j0]) { // no event in (tstart, tstop]
        ignore1z[j] = 1;
      } else {
        ignore1z[j] = 0;
      }
    }
  }
  
  IntegerVector ignoren(n); // final ignore vector
  for (int i=0; i<n; ++i) {
    ignoren[order0[i]] = ignore1z[i];
  }
  
  int nused = n - sum(ignoren);
  
  // sort by stopping time in descending order within each stratum
  IntegerVector order1 = seq(0, n-1);
  std::sort(order1.begin(), order1.end(), [&](int i, int j) {
    if (ignoren[i] != ignoren[j]) return ignoren[i] < ignoren[j];
    if (stratumn[i] != stratumn[j]) return stratumn[i] < stratumn[j];
    if (tstopn[i] != tstopn[j]) return tstopn[i] > tstopn[j];
    return eventn[i] < eventn[j];
  });
  
  IntegerVector stratum1 = stratumn[order1];
  NumericVector tstart1 = tstartn[order1];
  NumericVector tstop1 = tstopn[order1];
  IntegerVector event1 = eventn[order1];
  NumericVector weight1 = weightn[order1];
  NumericVector offset1 = offsetn[order1];
  IntegerVector ignore1 = ignoren[order1];
  NumericMatrix z1 = subset_matrix_by_row(zn, order1);
  
  // sort by starting time in descending order within each stratum
  IntegerVector order10 = seq(0, n-1);
  std::sort(order10.begin(), order10.end(), [&](int i, int j) {
    if (ignore1[i] != ignore1[j]) return ignore1[i] < ignore1[j];
    if (stratum1[i] != stratum1[j]) return stratum1[i] < stratum1[j];
    return tstart1[i] > tstart1[j];
  });
  
  coxparams param = {nused, stratum1, tstart1, tstop1, event1,
                     weight1, offset1, z1, order10, method};
  
  List out = f_der_i_2(p, beta, &param);
  NumericMatrix score_i = as<NumericMatrix>(out["score_i"]);
  NumericMatrix imat_i = as<NumericMatrix>(out["imat_i"]);
  
  // order data by ascending tstop
  IntegerVector order2 = seq(0, n-1);
  std::sort(order2.begin(), order2.end(), [&](int i, int j) {
    return tstop1[i] < tstop1[j];
  });
  
  IntegerVector stratum2 = stratum1[order2];
  NumericVector tstop2 = tstop1[order2];
  IntegerVector event2 = event1[order2];
  NumericMatrix score_i2 = subset_matrix_by_row(score_i, order2);
  NumericMatrix imat_i2 = subset_matrix_by_row(imat_i, order2);
  
  // only consider event times
  IntegerVector q = which(event2 == 1);
  IntegerVector stratum3 = stratum2[q];
  NumericVector tstop3 = tstop2[q];
  NumericMatrix score_i3 = subset_matrix_by_row(score_i2, q);
  NumericMatrix imat_i3 = subset_matrix_by_row(imat_i2, q);
  int nq = static_cast<int>(q.size());
  
  // identify unique event times
  NumericVector t = unique(tstop3);
  t.sort();
  
  int nt = static_cast<int>(t.size());
  NumericMatrix score_i4(nt, p);
  NumericMatrix imat_i4(nt, p*p);
  
  NumericVector tstop3a = clone(tstop3);
  tstop3a.push_front(-1.0); // to facilitate the loop
  int k = -1;
  for (int i=1; i<=nq; ++i) {
    if (tstop3a[i] != tstop3a[i-1]) k++;
    
    for (int j=0; j<p; ++j) {
      score_i4(k,j) += score_i3(i-1,j);
      for (int l=0; l<p; ++l) {
        imat_i4(k,j*p+l) += imat_i3(i-1,j*p+l);
      }
    }
  }
  
  // transformed time points
  NumericVector g(nt);
  if (transform == "identity") {
    for (int i=0; i<nt; ++i) {
      g[i] = t[i];
    }
  } else if (transform == "log") {
    for (int i=0; i<nt; ++i) {
      g[i] = log(t[i]);
    }
  } else if (transform == "rank") {
    for (int i=0; i<nt; ++i) {
      g[i] = static_cast<double>(i+1);
    }
  } else if (transform == "km") {
    DataFrame temp = DataFrame::create(
      Named("tstart") = tstart1,
      Named("tstop") = tstop1,
      Named("event") = event1
    );
    
    DataFrame df_km = kmest(temp, "", "", "tstart", "tstop", "event",
                            "", "none", 0.95, false);
    
    NumericVector surv = df_km["surv"];
    surv.push_front(1.0); // for left-continuous step function
    for (int i=0; i<nt; ++i) {
      g[i] = 1.0 - surv[i];
    }
  } else {
    Function f(transform);
    for (int i=0; i<nt; ++i) {
      g[i] = as<double>(f(t[i]));
    }
  }
  
  // score for theta
  NumericVector u_theta(p);
  for (int i=0; i<nt; ++i) {
    for (int j=0; j<p; ++j) {
      u_theta[j] += score_i4(i,j)*g[i];
    }
  }
  
  // covariance between u_theta and u_beta
  NumericMatrix imat_theta_beta(p, p);
  for (int i=0; i<nt; ++i) {
    for (int j=0; j<p; ++j) {
      for (int l=0; l<p; ++l) {
        imat_theta_beta(j,l) += imat_i4(i,j*p+l)*g[i];
      }
    }
  }
  
  // covariance for u_theta
  NumericMatrix imat_theta(p, p);
  for (int i=0; i<nt; ++i) {
    for (int j=0; j<p; ++j) {
      for (int l=0; l<p; ++l) {
        imat_theta(j,l) += imat_i4(i,j*p+l)*g[i]*g[i];
      }
    }
  }
  
  // conditional variance for u_theta given u_beta
  NumericMatrix vtheta(p, p);
  for (int j=0; j<p; ++j) {
    for (int l=0; l<p; ++l) {
      vtheta(j,l) = imat_theta(j,l);
      for (int s=0; s<p; ++s) {
        for (int r=0; r<p; ++r) {
          vtheta(j,l) -= imat_theta_beta(j,s)*vbeta(s,r)*imat_theta_beta(l,r);
        }
      }
    }
  }
  
  // individual score test for theta = 0
  NumericVector score_test(p);
  for (int j=0; j<p; ++j) {
    score_test[j] = u_theta[j]*u_theta[j]/vtheta(j,j);
  }
  
  // global score test for theta = 0
  double score_test_global = 0.0;
  NumericMatrix vtheta_inv = invsympd(vtheta, p, 1e-12);
  for (int j=0; j<p; ++j) {
    for (int l=0; l<p; ++l) {
      score_test_global += u_theta[j]*vtheta_inv(j,l)*u_theta[l];
    }
  }
  
  score_test.push_back(score_test_global);
  NumericVector df(p+1,1.0);
  df[p] = p; // degrees of freedom for the global test
  
  // p-values
  NumericVector p_values(p+1);
  for (int j=0; j<=p; ++j) {
    p_values[j] = 1.0 - R::pchisq(score_test[j], df[j], true, false);
  }

  // assemble the table
  NumericMatrix table(p+1,3);
  for (int j=0; j<=p; ++j) {
    table(j,0) = score_test[j];
    table(j,1) = df[j];
    table(j,2) = p_values[j];
  }
  
  StringVector covariates_vec(p+1);
  for (int j=0; j<p; ++j) covariates_vec[j] = covariates[j];
  covariates_vec[p] = "GLOBAL";
  rownames(table) = covariates_vec;
  colnames(table) = StringVector::create("chisq", "df", "p");
  
  // obtain scaled schoenfeld residuals
  List resid_list = residuals_phregcpp(
    p, beta, vbeta, resmart, data, stratum, time, time2, event, 
    covariates, weight, offset, "", ties, "scaledsch", false, true);
  
  NumericMatrix sresid = as<NumericMatrix>(resid_list["resid"]);
  colnames(sresid) = covariates;
    
  // repeat the g values according to the number of events
  NumericVector g_rep(nq);
  int idx = 0;
  for (int i=0; i<nt; ++i) {
    int count = 0;
    for (int j=0; j<nq; ++j) {
      if (tstop3[j] == t[i]) {
        count++;
      }
    }
    for (int j=0; j<count; ++j) {
      g_rep[idx] = g[i];
      idx++;
    }
  }
  
  List result = List::create(
    Named("table") = table,
    Named("x") = g_rep,
    Named("time") = tstop3,
    Named("y") = sresid,
    Named("var") = vbeta*nq,
    Named("transform") = transform
  );
  
  if (has_stratum) {
    result.push_back(stratum3, "strata");
  }
  
  return result;
}
