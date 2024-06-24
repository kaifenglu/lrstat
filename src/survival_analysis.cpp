#include <Rcpp.h>
#include <R_ext/Applic.h>
#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


//' @title Kaplan-Meier estimates of the survival curve
//' @description Obtains the Kaplan-Meier estimates of the survival curve.
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
//' @param conftype The type of confidence interval. One of "none",
//'   "plain", "log", "log-log" (the default), or "arcsin".
//'   The arcsin option bases the intervals on asin(sqrt(survival)).
//' @param confint The level of the two-sided confidence interval for
//'   the survival probabilities. Defaults to 0.95.
//'
//' @return A data frame with the following variables:
//'
//' * \code{rep}: The replication.
//'
//' * \code{stratum}: The stratum.
//'
//' * \code{size}: The number of subjects in the stratum.
//'
//' * \code{time}: The event time.
//'
//' * \code{nrisk}: The number of subjects at risk.
//'
//' * \code{nevent}: The number of subjects having the event.
//'
//' * \code{survival}: The Kaplan-Meier estimate of the survival probability.
//'
//' * \code{stderr}: The standard error of the estimated survival
//'   probability based on the Greendwood formula.
//'
//' * \code{lower}: The lower bound of confidence interval if requested.
//'
//' * \code{upper}: The upper bound of confidence interval if requested.
//'
//' * \code{confint}: The level of confidence interval if requested.
//'
//' * \code{conftype}: The type of confidence interval if requested.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' kmest(data = aml, stratum = "x",
//'       time = "time", event = "status")
//'
//' @export
// [[Rcpp::export]]
DataFrame kmest(const DataFrame data,
                const std::string rep = "rep",
                const std::string stratum = "stratum",
                const std::string time = "time",
                const std::string event = "event",
                const std::string conftype = "log-log",
                const double confint = 0.95) {
  int h, i, j, n = data.nrows();

  bool has_rep = hasVariable(data, rep);
  bool has_stratum = hasVariable(data, stratum);
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);

  if (!has_time) {
    stop("data must contain the time variable");
  }

  if (!has_event) {
    stop("data must contain the event variable");
  }

  NumericVector timen = data[time];
  NumericVector eventn = data[event];

  if (is_true(any(timen <= 0))) {
    stop("time must be positive for each subject");
  }

  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }


  // create the numeric rep variable
  IntegerVector repn(n);
  IntegerVector repwi;
  NumericVector repwn;
  StringVector repwc;
  if (!has_rep) {
    repn.fill(1);
  } else {
    if (TYPEOF(data[rep]) == INTSXP) {
      IntegerVector repv = data[rep];
      repwi = unique(repv);
      repwi.sort();
      repn = match(repv, repwi);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      NumericVector repv = data[rep];
      repwn = unique(repv);
      repwn.sort();
      repn = match(repv, repwn);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      StringVector repv = data[rep];
      repwc = unique(repv);
      repwc.sort();
      repn = match(repv, repwc);
    } else {
      stop("incorrect type for the rep variable in the input data");
    }
  }


  // create the numeric stratum variable
  IntegerVector stratumn(n);
  IntegerVector stratumwi;
  NumericVector stratumwn;
  StringVector stratumwc;
  if (!has_stratum) {
    stratumn.fill(1);
  } else {
    if (TYPEOF(data[stratum]) == INTSXP) {
      IntegerVector stratumv = data[stratum];
      stratumwi = unique(stratumv);
      stratumwi.sort();
      stratumn = match(stratumv, stratumwi);
    } else if (TYPEOF(data[stratum]) == REALSXP) {
      NumericVector stratumv = data[stratum];
      stratumwn = unique(stratumv);
      stratumwn.sort();
      stratumn = match(stratumv, stratumwn);
    } else if (TYPEOF(data[stratum]) == STRSXP) {
      StringVector stratumv = data[stratum];
      stratumwc = unique(stratumv);
      stratumwc.sort();
      stratumn = match(stratumv, stratumwc);
    } else {
      stop("incorrect type for the stratum variable in the input data");
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

  if (confint <= 0 || confint >= 1) {
    stop("confint must lie between 0 and 1");
  }


  // confidence interval for survival probability
  double z = R::qnorm((1.0 + confint)/2.0, 0, 1, 1, 0);

  auto f = [ct, z](double surv, double sesurv)->NumericVector {
    double grad, hw, lower = NA_REAL, upper = NA_REAL;
    if (ct == "plain") {
      lower = std::max(surv - z*sesurv, 0.0);
      upper = std::min(surv + z*sesurv, 1.0);
    } else if (ct == "log") {
      grad = 1.0/surv;
      hw = z*grad*sesurv;
      lower = exp(log(surv) - hw);
      upper = std::min(exp(log(surv) + hw), 1.0);
    } else if (ct == "log-log") {
      grad = 1.0/(surv*log(surv));
      hw = z*grad*sesurv;
      lower = exp(-exp(log(-log(surv)) - hw));
      upper = exp(-exp(log(-log(surv)) + hw));
    } else if (ct == "logit") {
      grad = 1.0/(surv*(1.0-surv));
      hw = z*grad*sesurv;
      lower = R::plogis(R::qlogis(surv, 0, 1, 1, 0) - hw, 0, 1, 1, 0);
      upper = R::plogis(R::qlogis(surv, 0, 1, 1, 0) + hw, 0, 1, 1, 0);
    } else if (ct == "arcsin") {
      grad = 1.0/(2.0*sqrt(surv*(1.0 - surv)));
      hw = z*grad*sesurv;
      lower = pow(sin(asin(sqrt(surv)) - hw), 2);
      upper = pow(sin(asin(sqrt(surv)) + hw), 2);
    }

    return NumericVector::create(lower, upper);
  };

  // sort the data by rep
  IntegerVector order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return repn[i] < repn[j];
  });

  repn = repn[order];
  stratumn = stratumn[order];
  timen = timen[order];
  eventn = eventn[order];

  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (i=1; i<n; i++) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }

  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);

  IntegerVector rep0(n, NA_INTEGER);
  IntegerVector stratum0(n), size0(n);
  NumericVector time0(n), nrisk0(n), nevent0(n);
  NumericVector surv0(n), sesurv0(n);
  NumericVector lower0(n), upper0(n);

  int index = 0;
  for (h=0; h<nreps; h++) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());
    int iter = repn[q1[0]];

    IntegerVector stratum1 = stratumn[q1];
    NumericVector time1 = timen[q1];
    NumericVector event1 = eventn[q1];

    // sort by stratum, time, and event with event in descending order
    IntegerVector order1 = seq(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      return (stratum1[i] < stratum1[j]) ||
        ((stratum1[i] == stratum1[j]) && (time1[i] < time1[j])) ||
        ((stratum1[i] == stratum1[j]) && (time1[i] == time1[j]) &&
        (event1[i] > event1[j]));
    });

    stratum1 = stratum1[order1];
    time1 = time1[order1];
    event1 = event1[order1];

    // identify the locations of the unique values of stratum
    IntegerVector idx1(1,0);
    for (i=1; i<n1; i++) {
      if (stratum1[i] != stratum1[i-1]) {
        idx1.push_back(i);
      }
    }

    int nstrata = static_cast<int>(idx1.size());
    idx1.push_back(n1);

    for (i=0; i<nstrata; i++) {
      IntegerVector q2 = Range(idx1[i], idx1[i+1]-1);
      NumericVector time2 = time1[q2];
      NumericVector event2 = event1[q2];

      int s = stratum1[q2[0]], n2 = static_cast<int>(q2.size());
      double t, nrisk, nevent, surv = 1, vcumhaz = 0, sesurv;
      bool cache = 0;
      for (j=0; j<n2; j++) {
        if (((j == 0) && (event2[j] == 1)) ||
            ((j >= 1) && (event2[j] == 1) && (time2[j] > time2[j-1]))) {
          // new event
          // add the info for the previous event
          if (cache) {
            surv = surv*(1.0 - nevent/nrisk);
            if (nrisk > nevent) {
              vcumhaz = vcumhaz + nevent/(nrisk*(nrisk - nevent));
            } else {
              vcumhaz = NA_REAL;
            }
            sesurv = surv*sqrt(vcumhaz);

            rep0[index] = iter;
            stratum0[index] = s;
            size0[index] = n2;
            time0[index] = t;
            nrisk0[index] = nrisk;
            nevent0[index] = nevent;
            surv0[index] = surv;
            sesurv0[index] = sesurv;

            if (ct != "none") {
              NumericVector ci = f(surv, sesurv);
              lower0[index] = ci[0];
              upper0[index] = ci[1];
            }

            index++;
          }

          // update the buffer for the current event time
          t = time2[j];
          nrisk = n2-j;
          nevent = 1;

          cache = 1;
        } else if ((j >= 1) && (event2[j] == 1) && (event2[j-1] == 1) &&
          (time2[j] == time2[j-1])) { // tied event
          nevent = nevent + 1;
        } else if ((j >= 1) && (event2[j] == 0) && (event2[j-1] == 1)) {
          // new censoring
          // add the info for the previous event
          surv = surv*(1.0 - nevent/nrisk);
          if (nrisk > nevent) {
            vcumhaz = vcumhaz + nevent/(nrisk*(nrisk - nevent));
          } else {
            vcumhaz = NA_REAL;
          }
          sesurv = surv*sqrt(vcumhaz);

          rep0[index] = iter;
          stratum0[index] = s;
          size0[index] = n2;
          time0[index] = t;
          nrisk0[index] = nrisk;
          nevent0[index] = nevent;
          surv0[index] = surv;
          sesurv0[index] = sesurv;

          if (ct != "none") {
            NumericVector ci = f(surv, sesurv);
            lower0[index] = ci[0];
            upper0[index] = ci[1];
          }

          index++;

          // empty the cache for the current event time
          cache = 0;
        }
      }

      // add the info for the last event
      if (cache) {
        surv = surv*(1.0 - nevent/nrisk);
        if (nrisk > nevent) {
          vcumhaz = vcumhaz + nevent/(nrisk*(nrisk - nevent));
        } else {
          vcumhaz = NA_REAL;
        }
        sesurv = surv*sqrt(vcumhaz);

        rep0[index] = iter;
        stratum0[index] = s;
        size0[index] = n2;
        time0[index] = t;
        nrisk0[index] = nrisk;
        nevent0[index] = nevent;
        surv0[index] = surv;
        sesurv0[index] = sesurv;

        if (ct != "none") {
          NumericVector ci = f(surv, sesurv);
          lower0[index] = ci[0];
          upper0[index] = ci[1];
        }

        index++;
      }
    }
  }

  // only keep nonmissing records
  LogicalVector sub = !is_na(rep0);
  if (is_false(any(sub))) {
    stop("no replication enables valid inference");
  }

  rep0 = rep0[sub];
  stratum0 = stratum0[sub];
  size0 = size0[sub];
  time0 = time0[sub];
  nrisk0 = nrisk0[sub];
  nevent0 = nevent0[sub];
  surv0 = surv0[sub];
  sesurv0 = sesurv0[sub];

  DataFrame result = DataFrame::create(
    _["size"] = size0,
    _["time"] = time0,
    _["nrisk"] = nrisk0,
    _["nevent"] = nevent0,
    _["survival"] = surv0,
    _["stderr"] = sesurv0);


  if (ct != "none") {
    result.push_back(lower0[sub], "lower");
    result.push_back(upper0[sub], "upper");
    result.push_back(confint, "confint");
    result.push_back(conftype, "conftype");
  }

  if (has_rep) {
    if (TYPEOF(data[rep]) == INTSXP) {
      result.push_back(repwi[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      result.push_back(repwn[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      result.push_back(repwc[rep0-1], rep);
    }
  }

  if (has_stratum) {
    if (TYPEOF(data[stratum]) == INTSXP) {
      result.push_back(stratumwi[stratum0-1], stratum);
    } else if (TYPEOF(data[stratum]) == REALSXP) {
      result.push_back(stratumwn[stratum0-1], stratum);
    } else if (TYPEOF(data[stratum]) == STRSXP) {
      result.push_back(stratumwc[stratum0-1], stratum);
    }
  }

  return result;
}


//' @title Log-rank test of survival curve difference
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
//'   * \code{time}: The possibly right-censored survival time.
//'
//'   * \code{event}: The event indicator.
//'
//' @param rep The name of the replication variable in the input data.
//' @param stratum The name of the stratum variable in the input data.
//' @param treat The name of the treatment variable in the input data.
//' @param time The name of the time variable in the input data.
//' @param event The name of the event variable in the input data.
//' @inheritParams param_rho1
//' @inheritParams param_rho2
//'
//' @return A data frame with the following variables:
//'
//' * \code{rep}: The replication.
//'
//' * \code{uscore}: The numerator of the log-rank test statistic.
//'
//' * \code{vscore}: The variance of the log-rank score test statistic.
//'
//' * \code{logRankZ}: The Z-statistic value.
//'
//' * \code{logRankPValue}: The one-sided p-value.
//'
//' * \code{rho1}: The first parameter of the Fleming-Harrington weights.
//'
//' * \code{rho2}: The second parameter of the Fleming-Harrington weights.
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
                 const std::string rep = "rep",
                 const std::string stratum = "stratum",
                 const std::string treat = "treat",
                 const std::string time = "time",
                 const std::string event = "event",
                 const double rho1 = 0,
                 const double rho2 = 0) {

  int h, i, j, k, n = data.nrows();

  bool has_rep = hasVariable(data, rep);
  bool has_stratum = hasVariable(data, stratum);
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
  if (TYPEOF(data[treat]) == LGLSXP) {
    IntegerVector treatv = data[treat];
    treatwi = unique(treatv);
    if (treatwi.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    treatwi = IntegerVector::create(1,0);
    treatn = 2 - treatv;
  } else if (TYPEOF(data[treat]) == INTSXP) {
    IntegerVector treatv = data[treat];
    treatwi = unique(treatv);
    if (treatwi.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    treatwi.sort();
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwi == 0) | (treatwi == 1)))) {
      treatwi = IntegerVector::create(1,0);
      treatn = 2 - treatv;
    } else {
      treatn = match(treatv, treatwi);
    }
  } else if (TYPEOF(data[treat]) == REALSXP) {
    NumericVector treatv = data[treat];
    treatwn = unique(treatv);
    if (treatwn.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    treatwn.sort();
    treatn = match(treatv, treatwn);
  } else if (TYPEOF(data[treat]) == STRSXP) {
    StringVector treatv = data[treat];
    treatwc = unique(treatv);
    if (treatwc.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    treatwc.sort();
    treatn = match(treatv, treatwc);
  } else {
    stop("incorrect type for the treat variable in the input data");
  }


  NumericVector timen = data[time];
  NumericVector eventn = data[event];

  if (is_true(any(timen <= 0))) {
    stop("time must be positive for each subject");
  }

  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }


  // create the numeric rep variable
  IntegerVector repn(n);
  IntegerVector repwi;
  NumericVector repwn;
  StringVector repwc;
  if (!has_rep) {
    repn.fill(1);
  } else {
    if (TYPEOF(data[rep]) == INTSXP) {
      IntegerVector repv = data[rep];
      repwi = unique(repv);
      repwi.sort();
      repn = match(repv, repwi);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      NumericVector repv = data[rep];
      repwn = unique(repv);
      repwn.sort();
      repn = match(repv, repwn);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      StringVector repv = data[rep];
      repwc = unique(repv);
      repwc.sort();
      repn = match(repv, repwc);
    } else {
      stop("incorrect type for the rep variable in the input data");
    }
  }


  // create the numeric stratum variable
  IntegerVector stratumn(n);
  IntegerVector stratumwi;
  NumericVector stratumwn;
  StringVector stratumwc;
  if (!has_stratum) {
    stratumn.fill(1);
  } else {
    if (TYPEOF(data[stratum]) == INTSXP) {
      IntegerVector stratumv = data[stratum];
      stratumwi = unique(stratumv);
      stratumwi.sort();
      stratumn = match(stratumv, stratumwi);
    } else if (TYPEOF(data[stratum]) == REALSXP) {
      NumericVector stratumv = data[stratum];
      stratumwn = unique(stratumv);
      stratumwn.sort();
      stratumn = match(stratumv, stratumwn);
    } else if (TYPEOF(data[stratum]) == STRSXP) {
      StringVector stratumv = data[stratum];
      stratumwc = unique(stratumv);
      stratumwc.sort();
      stratumn = match(stratumv, stratumwc);
    } else {
      stop("incorrect type for the stratum variable in the input data");
    }
  }


  if (rho1 < 0) {
    stop("rho1 must be non-negative");
  }

  if (rho2 < 0) {
    stop("rho2 must be non-negative");
  }


  // sort the data by rep
  IntegerVector order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return repn[i] < repn[j];
  });

  repn = repn[order];
  stratumn = stratumn[order];
  treatn = treatn[order];
  timen = timen[order];
  eventn = eventn[order];

  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (i=1; i<n; i++) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }

  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);

  IntegerVector rep0(nreps, NA_INTEGER);
  NumericVector uscore0(nreps), vscore0(nreps);
  NumericVector logRankZ0(nreps), logRankPValue0(nreps);

  bool noerr = 1;
  int index = 0;
  for (h=0; h<nreps; h++) {
    bool skip = 0;
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());

    IntegerVector stratum1 = stratumn[q1];
    IntegerVector treat1 = treatn[q1];
    NumericVector time1 = timen[q1];
    NumericVector event1 = eventn[q1];

    // sort by stratum, time, and event with event in descending order
    IntegerVector order1 = seq(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      return (stratum1[i] < stratum1[j]) ||
        ((stratum1[i] == stratum1[j]) && (time1[i] < time1[j])) ||
        ((stratum1[i] == stratum1[j]) && (time1[i] == time1[j]) &&
        (event1[i] > event1[j]));
    });

    stratum1 = stratum1[order1];
    treat1 = treat1[order1];
    time1 = time1[order1];
    event1 = event1[order1];

    // identify the locations of the unique values of stratum
    IntegerVector idx1(1,0);
    for (i=1; i<n1; i++) {
      if (stratum1[i] != stratum1[i-1]) {
        idx1.push_back(i);
      }
    }

    int nstrata = static_cast<int>(idx1.size());
    idx1.push_back(n1);

    double uscore = 0.0, vscore = 0.0;
    for (i=0; i<nstrata; i++) {
      IntegerVector q = Range(idx1[i], idx1[i+1]-1);
      IntegerVector treat2 = treat1[q];
      NumericVector time2 = time1[q];
      NumericVector event2 = event1[q];
      int n2 = static_cast<int>(q.size());

      // running index of number of subjects at risk
      int nriskx = n2;
      int nrisk1x = sum(treat2 == 1), nrisk2x = sum(treat2 == 2);

      if ((nrisk1x == 0) || (nrisk2x == 0)) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else if (TYPEOF(data[rep]) == INTSXP) {
          reperr = " " + rep + " = " + std::to_string(repwi[repn[idx[h]]-1]);
        } else if (TYPEOF(data[rep]) == REALSXP) {
          reperr = " " + rep + " = " + std::to_string(repwn[repn[idx[h]]-1]);
        } else {
          reperr = " " + rep + " = " + repwc[repn[idx[h]]-1];
        }

        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else if (TYPEOF(data[stratum]) == INTSXP) {
          stratumerr = " " + stratum + " = " +
            std::to_string(stratumwi[stratum1[idx1[i]]-1]);
        } else if (TYPEOF(data[stratum]) == REALSXP) {
          stratumerr = " " + stratum + " = " +
            std::to_string(stratumwn[stratum1[idx1[i]]-1]);
        } else {
          stratumerr = " " + stratum + " = " +
            stratumwc[stratum1[idx1[i]]-1];
        }


        k = nrisk1x == 0 ? 0 : 1;
        std::string treaterr;
        if ((TYPEOF(data[treat]) == LGLSXP) ||
            (TYPEOF(data[treat]) == INTSXP)) {
          treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
        } else if (TYPEOF(data[treat]) == REALSXP) {
          treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
        } else {
          treaterr = " " + treat + " = " + treatwc[k];
        }

        std::string str1 = "Warning: The data set does not contain";
        std::string errmsg = str1 + treaterr;
        if (!reperr.empty() || !stratumerr.empty()) {
          errmsg = errmsg + ":" + reperr + stratumerr;
        }

        if (noerr) {
          Rcout << errmsg << "\n";
          Rcout << "Additional warning messages are suppressed" << "\n";
          noerr = 0;
        }

        skip = 1;
        break;
      }


      double nrisk = nriskx, nrisk1 = nrisk1x, nrisk2 = nrisk2x;
      double nevent, nevent1, nevent2;
      double surv = 1.0;

      NumericVector nrisk0(n2, NA_REAL), nevent0(n2);
      NumericVector nrisk10(n2), nevent10(n2);
      NumericVector nrisk20(n2), nevent20(n2);
      NumericVector surv0(n2);

      int index1 = 0;
      bool cache = 0;
      for (j=0; j<n2; j++) {
        if (((j == 0) && (event2[j] == 1)) ||
            ((j >= 1) && (event2[j] == 1) && (time2[j] > time2[j-1]))) {
          // new event
          // add the info for the previous event
          if (cache) {
            surv = surv*(1.0 - nevent/nrisk);

            nrisk0[index1] = nrisk;
            nevent0[index1] = nevent;
            nrisk10[index1] = nrisk1;
            nevent10[index1] = nevent1;
            nrisk20[index1] = nrisk2;
            nevent20[index1] = nevent2;
            surv0[index1] = surv;

            index1++;
          }

          // update the cache for the current event time
          nrisk = nriskx;
          nrisk1 = nrisk1x;
          nrisk2 = nrisk2x;
          nevent = 1;
          nevent1 = (treat2[j] == 1);
          nevent2 = (treat2[j] == 2);

          cache = 1;
        } else if ((j >= 1) && (event2[j] == 1) && (event2[j-1] == 1) &&
          (time2[j] == time2[j-1])) { // tied event
          nevent = nevent + 1;
          nevent1 = nevent1 + (treat2[j] == 1);
          nevent2 = nevent2 + (treat2[j] == 2);
        } else if ((j >= 1) && (event2[j] == 0) && (event2[j-1] == 1)) {
          // new censoring
          // add the info for the previous event
          surv = surv*(1.0 - nevent/nrisk);

          nrisk0[index1] = nrisk;
          nevent0[index1] = nevent;
          nrisk10[index1] = nrisk1;
          nevent10[index1] = nevent1;
          nrisk20[index1] = nrisk2;
          nevent20[index1] = nevent2;
          surv0[index1] = surv;

          index1++;

          // empty the cache for the current event time
          cache = 0;
        }

        nriskx--;
        if (treat2[j] == 1) nrisk1x--;
        if (treat2[j] == 2) nrisk2x--;
      }

      // add the info for the last event
      if (cache) {
        surv = surv*(1.0 - nevent/nrisk);

        nrisk0[index1] = nrisk;
        nevent0[index1] = nevent;
        nrisk10[index1] = nrisk1;
        nevent10[index1] = nevent1;
        nrisk20[index1] = nrisk2;
        nevent20[index1] = nevent2;
        surv0[index1] = surv;

        index1++;
      }

      // only keep nonmissing records
      double uscore1 = 0.0, vscore1 = 0.0;
      LogicalVector sub = !is_na(nrisk0);
      if (is_true(any(sub))) { // at least 1 event
        nrisk0 = nrisk0[sub];
        nevent0 = nevent0[sub];
        nrisk10 = nrisk10[sub];
        nevent10 = nevent10[sub];
        nrisk20 = nrisk20[sub];
        nevent20 = nevent20[sub];
        surv0 = surv0[sub];

        int K = sum(sub);
        NumericVector w(K);
        surv0.push_front(1.0);
        for (k=0; k<K; k++) {
          double w1, w2;
          if (surv0[k] == 1.0) {
            w1 = 1.0;
            w2 = rho2 > 0.0 ? 0.0 : 1.0;
          } else if (surv0[k] == 0.0) {
            w1 = rho1 > 0.0 ? 0.0 : 1.0;
            w2 = 1.0;
          } else {
            w1 = pow(surv0[k], rho1);
            w2 = pow(1.0 - surv0[k], rho2);
          }
          w[k] = w1*w2;
        }

        for (k=0; k<K; k++) {
          uscore1 += w[k]*(nevent10[k] - nevent0[k]*nrisk10[k]/nrisk0[k]);
          if (nrisk0[k] > 1.0) {
            vscore1 += pow(w[k],2)*nevent0[k]*(nrisk0[k]-nevent0[k])*
              nrisk10[k]*nrisk20[k]/(pow(nrisk0[k],2)*(nrisk0[k]-1.0));
          }
        }
      }

      uscore += uscore1;
      vscore += vscore1;
    }

    // skip the replication if there is a stratum without both treatments
    if (skip) continue;

    rep0[index] = repn[idx[h]];
    uscore0[index] = uscore;
    vscore0[index] = vscore;
    logRankZ0[index] = uscore/sqrt(vscore);
    logRankPValue0[index] = R::pnorm(logRankZ0[index], 0, 1, 1, 0);

    index++;
  }

  // only keep nonmissing records
  LogicalVector sub = !is_na(rep0);
  if (is_false(any(sub))) {
    stop("no replication enables valid inference");
  }

  rep0 = rep0[sub];
  uscore0 = uscore0[sub];
  vscore0 = vscore0[sub];
  logRankZ0 = logRankZ0[sub];
  logRankPValue0 = logRankPValue0[sub];


  DataFrame result = DataFrame::create(
    _["uscore"] = uscore0,
    _["vscore"] = vscore0,
    _["logRankZ"] = logRankZ0,
    _["logRankPValue"] = logRankPValue0,
    _["rho1"] = rho1,
    _["rho2"] = rho2);

  if (has_rep) {
    if (TYPEOF(data[rep]) == INTSXP) {
      result.push_back(repwi[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      result.push_back(repwn[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      result.push_back(repwc[rep0-1], rep);
    }
  }

  return result;
}



// define functions in likelihood inference

// negative log likelihood
double f_nllik_1(int p, double *par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int person, i, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = exp(par[k]);
    }
  }

  double loglik = 0;
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];

    if (param->status[person] == 1) { // event
      double logsig = log(sigma);
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*(u - exp(u) - logsig);
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*(R::dnorm(u, 0, 1, 1) - logsig);
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*(R::dlogis(u, 0, 1, 1) - logsig);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*(R::dnorm(u, 0, 1, 1) - logsig);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*(R::dlogis(u, 0, 1, 1) - logsig);
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*log(exp(-exp(v)) - exp(-exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*log(R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*log(R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*log(R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*log(R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*log(1.0 - exp(-exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*R::pnorm(u, 0, 1, 1, 1);
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*R::plogis(u, 0, 1, 1, 1);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*R::pnorm(u, 0, 1, 1, 1);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*R::plogis(u, 0, 1, 1, 1);
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*(-exp(v));
      } else if (param->dist == "lognormal") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*R::pnorm(v, 0, 1, 0, 1);
      } else if (param->dist == "loglogistic") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*R::plogis(v, 0, 1, 0, 1);
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*R::pnorm(v, 0, 1, 0, 1);
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*R::plogis(v, 0, 1, 0, 1);
      }
    }
  }

  return -loglik;
}


// negative score vector
void f_nscore_1(int p, double *par, double *gr, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int person, i, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = exp(par[k]);
    }
  }

  NumericVector score(p);
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];
    NumericVector z = param->z(person, _)/sigma;
    k = param->strata[person] + nvar - 1;

    if (param->status[person] == 1) { // event
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -wt*(1 - exp(u));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -wt*(1 - exp(u));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*((1 - exp(u))*(-u) - 1);
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*u;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(u*u - 1);
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c0 = 1 - 2*R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*c0;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(c0*u - 1);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*u;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(u*u - 1);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c0 = 1 - 2*R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*c0;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(c0*u - 1);
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(exp(v - exp(v)) - exp(u - exp(u)))/
          (exp(-exp(v)) - exp(-exp(u)));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(exp(v - exp(v)) - exp(u - exp(u)))/
          (exp(-exp(v)) - exp(-exp(u)));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(exp(v - exp(v))*v - exp(u - exp(u))*u)/
          (exp(-exp(v)) - exp(-exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(R::dnorm(v, 0, 1, 0)*v - R::dnorm(u, 0, 1, 0)*u)/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(R::dlogis(v, 0, 1, 0) - R::dlogis(u, 0, 1, 0))/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(R::dlogis(v, 0, 1, 0)*v - R::dlogis(u, 0, 1, 0)*u)/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(R::dnorm(v, 0, 1, 0)*v - R::dnorm(u, 0, 1, 0)*u)/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(R::dlogis(v, 0, 1, 0) - R::dlogis(u, 0, 1, 0))/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(R::dlogis(v, 0, 1, 0)*v - R::dlogis(u, 0, 1, 0)*u)/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-exp(u - exp(u))/(1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-exp(u - exp(u))/(1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-R::dnorm(u, 0, 1, 0)/R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(-R::dnorm(u, 0, 1, 0)/R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(-R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "lognormal") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*R::dnorm(v, 0, 1, 0)/R::pnorm(v, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "loglogistic") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*R::plogis(v, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*R::dnorm(v, 0, 1, 0)/R::pnorm(v, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*R::plogis(v, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      }
    }
  }

  for (i=0; i<p; i++) gr[i] = -score[i];
}


// observed information matrix
NumericMatrix f_info_1(int p, double *par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int person, i, j, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = exp(par[k]);
    }
  }

  NumericMatrix imat(p,p);
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];
    NumericVector z = param->z(person, _)/sigma;
    k = param->strata[person] + nvar - 1;

    if (param->status[person] == 1) { // event
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*exp(u);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*exp(u);
        double c2 = wt*(exp(u)*u - (1 - exp(u)));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c2 = wt*2*u;
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += wt*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*2*R::dlogis(u, 0, 1, 0);
        double c2 = wt*(2*R::dlogis(u, 0, 1, 0)*u +
                        1 - 2*R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c2 = wt*2*u;
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += wt*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*2*R::dlogis(u, 0, 1, 0);
        double c2 = wt*(2*R::dlogis(u, 0, 1, 0)*u +
                        1 - 2*R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(pow((exp(v - exp(v)) - exp(u - exp(u)))/
                        (exp(-exp(v)) - exp(-exp(u))), 2) +
                          (exp(v - exp(v))*(1 - exp(v)) -
                          exp(u - exp(u))*(1 - exp(u)))/
                            (exp(-exp(v)) - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(pow((exp(v - exp(v)) - exp(u - exp(u)))/
                        (exp(-exp(v)) - exp(-exp(u))), 2) +
                          (exp(v - exp(v))*(1 - exp(v)) -
                          exp(u - exp(u))*(1 - exp(u)))/
                            (exp(-exp(v)) - exp(-exp(u))));
        double c2 = wt*((exp(v - exp(v)) - exp(u - exp(u)))*
                        (exp(v - exp(v))*v - exp(u - exp(u))*u)/
                          pow(exp(-exp(v)) - exp(-exp(u)), 2) +
                            (exp(v - exp(v))*(1 + (1 - exp(v))*v) -
                            exp(u - exp(u))*(1 + (1 - exp(u))*u))/
                              (exp(-exp(v)) - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((exp(v - exp(v))*v - exp(u - exp(u))*u)/
          (exp(-exp(v)) - exp(-exp(u))), 2) +
            (exp(v - exp(v))*(1 + (1 - exp(v))*v)*v
            - exp(u - exp(u))*(1 + (1 - exp(u))*u)*u)/
            (exp(-exp(v)) - exp(-exp(u))));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(pow((R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
                        (R::pnorm(v, 0, 1, 0, 0) -
                          R::pnorm(u, 0, 1, 0, 0)), 2) +
                          (-R::dnorm(v, 0, 1, 0)*v +
                          R::dnorm(u, 0, 1, 0)*u)/
                            (R::pnorm(v, 0, 1, 0, 0) -
                              R::pnorm(u, 0, 1, 0, 0)));
        double c2 = wt*((R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))*
                        (R::dnorm(v, 0, 1, 0)*v - R::dnorm(u, 0, 1, 0)*u)/
                          pow(R::pnorm(v, 0, 1, 0, 0) -
                            R::pnorm(u, 0, 1, 0, 0), 2) +
                            (R::dnorm(v, 0, 1, 0)*(1 - v*v) -
                            R::dnorm(u, 0, 1, 0)*(1 - u*u))/
                              (R::pnorm(v, 0, 1, 0, 0) -
                                R::pnorm(u, 0, 1, 0, 0)));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((R::dnorm(v, 0, 1, 0)*v -
          R::dnorm(u, 0, 1, 0)*u)/
            (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0)), 2) +
              (R::dnorm(v, 0, 1, 0)*(1 - v*v)*v -
              R::dnorm(u, 0, 1, 0)*(1 - u*u)*u)/
                (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0)));
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double q1 = R::plogis(v, 0, 1, 0, 0);
        double q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*(pow((q1*(1-q1) - q2*(1-q2))/(q1-q2), 2) +
                        (q1*(1-q1)*(2*q1-1) - q2*(1-q2)*(2*q2-1))/(q1-q2));
        double c2 = wt*((q1*(1-q1) - q2*(1-q2))*
                        (q1*(1-q1)*v - q2*(1-q2)*u)/pow(q1-q2, 2) +
                        (q1*(1-q1)*(1+(2*q1-1)*v) -
                        q2*(1-q2)*(1+(2*q2-1)*u))/(q1-q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((q1*(1-q1)*v - q2*(1-q2)*u)/(q1-q2), 2) +
          (q1*(1-q1)*(1+(2*q1-1)*v)*v - q2*(1-q2)*(1+(2*q2-1)*u)*u)/(q1-q2));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(pow((R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
                        (R::pnorm(v, 0, 1, 0, 0) -
                          R::pnorm(u, 0, 1, 0, 0)), 2) +
                          (-R::dnorm(v, 0, 1, 0)*v +
                          R::dnorm(u, 0, 1, 0)*u)/
                            (R::pnorm(v, 0, 1, 0, 0) -
                              R::pnorm(u, 0, 1, 0, 0)));
        double c2 = wt*((R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))*
                        (R::dnorm(v, 0, 1, 0)*v - R::dnorm(u, 0, 1, 0)*u)/
                          pow(R::pnorm(v, 0, 1, 0, 0) -
                            R::pnorm(u, 0, 1, 0, 0), 2) +
                            (R::dnorm(v, 0, 1, 0)*(1 - v*v) -
                            R::dnorm(u, 0, 1, 0)*(1 - u*u))/
                              (R::pnorm(v, 0, 1, 0, 0) -
                                R::pnorm(u, 0, 1, 0, 0)));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((R::dnorm(v, 0, 1, 0)*v -
          R::dnorm(u, 0, 1, 0)*u)/
            (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0)), 2) +
              (R::dnorm(v, 0, 1, 0)*(1 - v*v)*v -
              R::dnorm(u, 0, 1, 0)*(1 - u*u)*u)/
                (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0)));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double q1 = R::plogis(v, 0, 1, 0, 0);
        double q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*(pow((q1*(1-q1) - q2*(1-q2))/(q1-q2), 2) +
                        (q1*(1-q1)*(2*q1-1) - q2*(1-q2)*(2*q2-1))/(q1-q2));
        double c2 = wt*((q1*(1-q1) - q2*(1-q2))*
                        (q1*(1-q1)*v - q2*(1-q2)*u)/pow(q1-q2, 2) +
                        (q1*(1-q1)*(1+(2*q1-1)*v) -
                        q2*(1-q2)*(1+(2*q2-1)*u))/(q1-q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((q1*(1-q1)*v - q2*(1-q2)*u)/(q1-q2), 2) +
          (q1*(1-q1)*(1+(2*q1-1)*v)*v - q2*(1-q2)*(1+(2*q2-1)*u)*u)/(q1-q2));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(pow(exp(u - exp(u))/(1 - exp(-exp(u))), 2) -
                        exp(u - exp(u))*(1 - exp(u))/(1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(pow(exp(u - exp(u))/(1 - exp(-exp(u))), 2) -
                        exp(u - exp(u))*(1 - exp(u))/(1 - exp(-exp(u))));
        double c2 = wt*(pow(exp(u - exp(u))/(1 - exp(-exp(u))), 2)*u -
                        exp(u - exp(u))*(1 + (1 - exp(u))*u)/
                          (1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(pow(R::dnorm(u, 0, 1, 0)/
                        R::pnorm(u, 0, 1, 1, 0), 2) +
                          R::dnorm(u, 0, 1, 0)*u/R::pnorm(u, 0, 1, 1, 0));
        double c2 = wt*(pow(R::dnorm(u, 0, 1, 0)/
                        R::pnorm(u, 0, 1, 1, 0), 2)*u -
                          R::dnorm(u, 0, 1, 0)*(1 - u*u)/
                            R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*q2*(1-q2);
        double c2 = wt*(q2*(1-q2)*u - q2);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(pow(R::dnorm(u, 0, 1, 0)/
                        R::pnorm(u, 0, 1, 1, 0), 2) +
                          R::dnorm(u, 0, 1, 0)*u/R::pnorm(u, 0, 1, 1, 0));
        double c2 = wt*(pow(R::dnorm(u, 0, 1, 0)/
                        R::pnorm(u, 0, 1, 1, 0), 2)*u -
                          R::dnorm(u, 0, 1, 0)*(1 - u*u)/
                            R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*q2*(1-q2);
        double c2 = wt*(q2*(1-q2)*u - q2);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*u;
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        double c2 = wt*exp(v)*(1+v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*v;
      } else if (param->dist == "lognormal") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(pow(R::dnorm(v, 0, 1, 0)/
                        R::pnorm(v, 0, 1, 0, 0), 2) -
                          R::dnorm(v, 0, 1, 0)*v/R::pnorm(v, 0, 1, 0, 0));
        double c2 = wt*(pow(R::dnorm(v, 0, 1, 0)/
                        R::pnorm(v, 0, 1, 0, 0), 2)*v +
                          R::dnorm(v, 0, 1, 0)*(1 - v*v)/
                            R::pnorm(v, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*v;
      } else if (param->dist == "loglogistic") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double q1 = R::plogis(v, 0, 1, 0, 0);
        double c1 = wt*q1*(1-q1);
        double c2 = wt*(1-q1+q1*(1-q1)*v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*v;
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(pow(R::dnorm(v, 0, 1, 0)/
                        R::pnorm(v, 0, 1, 0, 0), 2) -
                          R::dnorm(v, 0, 1, 0)*v/R::pnorm(v, 0, 1, 0, 0));
        double c2 = wt*(pow(R::dnorm(v, 0, 1, 0)/
                        R::pnorm(v, 0, 1, 0, 0), 2)*v +
                          R::dnorm(v, 0, 1, 0)*(1 - v*v)/
                            R::pnorm(v, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*v;
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double q1 = R::plogis(v, 0, 1, 0, 0);
        double c1 = wt*q1*(1-q1);
        double c2 = wt*(1-q1+q1*(1-q1)*v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += c2*v;
      }
    }
  }

  for (i=0; i<p-1; i++) {
    for (j=i+1; j<p; j++) {
      imat(i,j) = imat(j,i);
    }
  }

  return imat;
}


// score residual matrix
NumericMatrix f_ressco_1(NumericVector par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int p = static_cast<int>(par.size());
  int person, i, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = exp(par[k]);
    }
  }

  NumericMatrix resid(n, p);
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];
    NumericVector z = param->z(person, _)/sigma;
    k = param->strata[person] + nvar - 1;

    if (param->status[person] == 1) { // event
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -wt*(1 - exp(u));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -wt*(1 - exp(u));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*((1 - exp(u))*(-u) - 1);
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*u;
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(u*u - 1);
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c0 = 1 - 2*R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*c0;
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(c0*u - 1);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*u;
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(u*u - 1);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c0 = 1 - 2*R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*c0;
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(c0*u - 1);
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(exp(v - exp(v)) - exp(u - exp(u)))/
          (exp(-exp(v)) - exp(-exp(u)));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(exp(v - exp(v)) - exp(u - exp(u)))/
          (exp(-exp(v)) - exp(-exp(u)));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(exp(v - exp(v))*v - exp(u - exp(u))*u)/
          (exp(-exp(v)) - exp(-exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(R::dnorm(v, 0, 1, 0)*v -
          R::dnorm(u, 0, 1, 0)*u)/
            (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(R::dlogis(v, 0, 1, 0) - R::dlogis(u, 0, 1, 0))/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(R::dlogis(v, 0, 1, 0)*v -
          R::dlogis(u, 0, 1, 0)*u)/
            (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(R::dnorm(v, 0, 1, 0) - R::dnorm(u, 0, 1, 0))/
          (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(R::dnorm(v, 0, 1, 0)*v -
          R::dnorm(u, 0, 1, 0)*u)/
            (R::pnorm(v, 0, 1, 0, 0) - R::pnorm(u, 0, 1, 0, 0));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(R::dlogis(v, 0, 1, 0) - R::dlogis(u, 0, 1, 0))/
          (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = wt*(R::dlogis(v, 0, 1, 0)*v -
          R::dlogis(u, 0, 1, 0)*u)/
            (R::plogis(v, 0, 1, 0, 0) - R::plogis(u, 0, 1, 0, 0));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-exp(u - exp(u))/(1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-exp(u - exp(u))/(1 - exp(-exp(u))));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-R::dnorm(u, 0, 1, 0)/R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(-R::dnorm(u, 0, 1, 0)/R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(-R::plogis(u, 0, 1, 0, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*exp(v);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "lognormal") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*R::dnorm(v, 0, 1, 0)/R::pnorm(v, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "loglogistic") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*R::plogis(v, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*R::dnorm(v, 0, 1, 0)/R::pnorm(v, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*R::plogis(v, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      }
    }
  }

  return resid;
}


//' @title Parametric regression models for failure time data
//' @description Obtains the parameter estimates from parametric
//' regression models with uncensored, right censored, left censored, or
//' interval censored data.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{rep}: The replication for by-group processing.
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{time}: The follow-up time for right censored data, or
//'     the left end of each interval for interval censored data.
//'
//'   * \code{time2}: The right end of each interval for interval
//'     censored data.
//'
//'   * \code{event}: The event indicator, normally 1=event, 0=no event.
//'
//'   * \code{covariates}: The values of baseline covariates.
//'     This is the full-rank design matrix (excluding the intercept)
//'     for the regression model, assuming that factor variables
//'     have already been expanded into dummy variables.
//'     The intercept will be added automatically.
//'
//'   * \code{weight}: The weight for each observation.
//'
//'   * \code{id}: The optional subject ID to group the score residuals
//'     in computing the robust sandwich variance.
//'
//' @param rep The name of the replication variable in the input data.
//' @param stratum The name of the stratum variable in the input data.
//' @param time The name of the time variable or the left end of each
//'   interval for interval censored data in the input data.
//' @param time2 The name of the right end of each interval for
//'   interval censored data in the input data.
//' @param event The name of the event variable in the input data
//'   for right censored data.
//' @param covariates The vector of names of baseline covariates
//'   in the input data.
//' @param weight The name of the weighting variable in the input data.
//' @param id The name of the id variable in the input data.
//' @param dist The assumed distribution for time to event. Options include
//'   "exponential", "weibull", "lognormal", and "loglogistic" to be
//'   modeled on the log-scale, and "normal" and "logistic" to be modeled
//'   on the original scale.
//' @param robust Whether a robust sandwich variance estimate should be
//'   computed. The default is TRUE if there are fractional weights or
//'   there is at least 1 id with >1 event. In the presence of the id
//'   variable, the score residual will be aggregated for each id when
//'   computing the robust sandwich variance estimate.
//'
//' @details There are two ways to specify the model, one for right censored
//' data through the time and event variables, and the other for interval
//' censored data through the time and time2 variables. For the second form,
//' we follow the convention used in SAS PROC LIFEREG:
//'
//' * If lower is not missing, upper is not missing, and lower is equal
//'   to upper, then there is no censoring and the event occurred at
//'   time lower.
//'
//' * If lower is not missing, upper is not missing, and lower < upper,
//'   then the event time is censored within the interval (lower, upper).
//'
//' * If lower is missing, but upper is not missing, then upper will be
//'   used as the left censoring value.
//'
//' * If lower is not missing, but upper is missing, then lower will be
//'   used as the right censoring value.
//'
//' * If lower is not missing, upper is not missing, but lower > upper,
//'   or if both lower and upper are missing, then the observation will
//'   not be used.
//'
//' @return A list with the following components:
//'
//' * \code{sumstat}: The data frame of summary statistics of model fit
//'   with the following variables:
//'
//'     - \code{rep}: The replication.
//'
//'     - \code{n}: The number of observations.
//'
//'     - \code{nevents}: The number of events.
//'
//'     - \code{loglik0}: The log-likelihood under null.
//'
//'     - \code{loglik1}: The maximum log-likelihood.
//'
//'     - \code{scoretest}: The score test statistic.
//'
//' * \code{parest}: The data frame of parameter estimates with the
//'   following variables:
//'
//'     - \code{rep}: The replication.
//'
//'     - \code{param}: The name of the covariate for the parameter estimate.
//'
//'     - \code{beta}: The parameter estimate.
//'
//'     - \code{sebeta}: The standard error of parameter estimate.
//'
//'     - \code{z}: The Wald test statistic.
//'
//'     - \code{expbeta}: The exponentiated parameter.
//'
//'     - \code{vbeta}: The covariance matrix for parameter estimates.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' library(dplyr)
//'
//' # right censored data
//' liferegr(data = rawdata %>% mutate(treat = 1*(treatmentGroup == 1)),
//'          rep = "iterationNumber", stratum = "stratum",
//'          time = "timeUnderObservation", event = "event",
//'          covariates = "treat", dist = "weibull")
//'
//' # tobit regression for left censored data
//' liferegr(data = tobin %>% mutate(time = ifelse(durable>0, durable, NA)),
//'          time = "time", time2 = "durable",
//'          covariates = c("age", "quant"), dist = "normal")
//'
//' @export
// [[Rcpp::export]]
List liferegr(const DataFrame data,
              const std::string rep = "rep",
              const std::string stratum = "stratum",
              const std::string time = "time",
              const std::string time2 = "time2",
              const std::string event = "event",
              const StringVector& covariates = "treat",
              const std::string weight = "weight",
              const std::string id = "id",
              const std::string dist = "weibull",
              bool robust = 0) {
  std::string dist1 = dist;
  std::for_each(dist1.begin(), dist1.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  if ((dist1 == "log-logistic") || (dist1 == "llogistic")) {
    dist1 = "loglogistic";
  } else if  ((dist1 == "log-normal") || (dist1 == "lnormal")) {
    dist1 = "lognormal";
  } else if (dist1 == "gaussian") {
    dist1 = "normal";
  }

  if (!((dist1 == "exponential") || (dist1 == "weibull") ||
      (dist1 == "lognormal") || (dist1 == "loglogistic") ||
      (dist1 == "normal") || (dist1 == "logistic"))) {
    std::string str1 = "dist must be exponential, weibull, lognormal,";
    std::string str2 = "loglogistic, normal, or logistic";
    std::string errmsg = str1 + " " + str2;
    stop(errmsg);
  }

  int h, i, j, k, n = data.nrows();
  int nvar = static_cast<int>(covariates.size()) + 1;

  bool has_rep = hasVariable(data, rep);
  bool has_stratum = hasVariable(data, stratum);

  bool has_time = hasVariable(data, time);
  if (!has_time) {
    stop("data must contain the time variable");
  }

  NumericVector timen = data[time];
  for (i=0; i<n; i++) {
    if (!R_isnancpp(timen[i]) && ((dist1 == "exponential") ||
        (dist1 == "weibull") || (dist1 == "lognormal") ||
        (dist1 == "loglogistic")) && (timen[i] <= 0)) {
      std::string str1 = "time must be positive for each subject for the";
      std::string str2 = "distribution";
      std::string errmsg = str1 + " " + dist1 + " " + str2;
      stop(errmsg);
    }
  }

  bool has_time2 = hasVariable(data, time2);
  NumericVector time2n(n);
  if (has_time2) {
    time2n = data[time2];
    for (i=0; i<n; i++) {
      if (!R_isnancpp(time2n[i]) && ((dist1 == "exponential") ||
          (dist1 == "weibull") || (dist1 == "lognormal") ||
          (dist1 == "loglogistic")) && (time2n[i] <= 0)) {
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
    eventn = data[event];
    if (is_true(any((eventn != 1) & (eventn != 0)))) {
      stop("event must be 1 or 0 for each subject");
    }

    if (is_true(all(eventn == 0))) {
      stop("at least 1 event is needed to fit the parametric model");
    }
  }

  NumericMatrix zn(n,nvar);
  for (i=0; i<n; i++) zn(i,0) = 1; // intercept

  for (j=0; j<nvar-1; j++) {
    String zj = covariates[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in covariates");
    }
    NumericVector u = data[zj];
    for (i=0; i<n; i++) zn(i,j+1) = u[i];
  }

  bool has_weight = hasVariable(data, weight);

  NumericVector weightn(n, 1.0);
  if (has_weight) {
    weightn = data[weight];
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }

  bool has_id = hasVariable(data, id);

  // create the numeric rep variable
  IntegerVector repn(n);
  IntegerVector repwi;
  NumericVector repwn;
  StringVector repwc;
  if (!has_rep) {
    repn.fill(1);
  } else {
    if (TYPEOF(data[rep]) == INTSXP) {
      IntegerVector repv = data[rep];
      repwi = unique(repv);
      repwi.sort();
      repn = match(repv, repwi);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      NumericVector repv = data[rep];
      repwn = unique(repv);
      repwn.sort();
      repn = match(repv, repwn);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      StringVector repv = data[rep];
      repwc = unique(repv);
      repwc.sort();
      repn = match(repv, repwc);
    } else {
      stop("incorrect type for the rep variable in the input data");
    }
  }

  // create the numeric stratum variable
  IntegerVector stratumn(n);
  if (!has_stratum) {
    stratumn.fill(1);
  } else {
    if (TYPEOF(data[stratum]) == INTSXP) {
      IntegerVector stratumv = data[stratum];
      IntegerVector stratumwi = unique(stratumv);
      stratumwi.sort();
      stratumn = match(stratumv, stratumwi);
    } else if (TYPEOF(data[stratum]) == REALSXP) {
      NumericVector stratumv = data[stratum];
      NumericVector stratumwn = unique(stratumv);
      stratumwn.sort();
      stratumn = match(stratumv, stratumwn);
    } else if (TYPEOF(data[stratum]) == STRSXP) {
      StringVector stratumv = data[stratum];
      StringVector stratumwc = unique(stratumv);
      stratumwc.sort();
      stratumn = match(stratumv, stratumwc);
    } else {
      stop("incorrect type for the stratum variable in the input data");
    }
  }

  IntegerVector stratumn1 = unique(stratumn);
  int nstrata = static_cast<int>(stratumn1.size());
  int p = dist1 == "exponential" ? nvar : (nvar+nstrata);

  // create the numeric id variable
  IntegerVector idn(n);
  if (!has_id) {
    idn = seq(1,n);
  } else {
    if (TYPEOF(data[id]) == INTSXP) {
      IntegerVector idv = data[id];
      IntegerVector idwi = unique(idv);
      idwi.sort();
      idn = match(idv, idwi);
    } else if (TYPEOF(data[id]) == REALSXP) {
      NumericVector idv = data[id];
      NumericVector idwn = unique(idv);
      idwn.sort();
      idn = match(idv, idwn);
    } else if (TYPEOF(data[id]) == STRSXP) {
      StringVector idv = data[id];
      StringVector idwc = unique(idv);
      idwc.sort();
      idn = match(idv, idwc);
    } else {
      stop("incorrect type for the id variable in the input data");
    }
  }

  // check if there is as least one id with more than 1 event
  IntegerVector idn0 = idn[eventn==1];
  IntegerVector idw0 = unique(idn0);
  bool dup_id = idw0.size() < idn0.size();

  if (R_isnancpp(robust)) {
    if ((has_id && dup_id) || is_true(any(weightn != floor(weightn)))) {
      robust = 1;
    } else {
      robust = 0;
    }
  }

  // sort the data by rep
  IntegerVector order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return repn[i] < repn[j];
  });

  repn = repn[order];
  stratumn = stratumn[order];
  timen = timen[order];
  time2n = time2n[order];
  eventn = eventn[order];
  weightn = weightn[order];
  idn = idn[order];

  NumericMatrix z(n,nvar);
  for (i=0; i<n; i++) {
    for (j=0; j<nvar; j++) {
      z(i,j) = zn(order[i], j);
    }
  }
  zn = z;

  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (i=1; i<n; i++) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }

  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);

  // variables in the output data sets
  IntegerVector rep01 = seq(1,nreps);
  IntegerVector nobs(nreps), nevents(nreps);
  NumericMatrix loglik(nreps,2);
  NumericVector scoretest(nreps);

  IntegerVector rep0(nreps*p);
  StringVector par0(nreps*p);
  NumericVector beta0(nreps*p), sebeta0(nreps*p), rsebeta0(nreps*p);
  NumericMatrix vbeta0(nreps*p,p), rvbeta0(nreps*p,p);

  for (h=0; h<nreps; h++) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());

    IntegerVector stratum1 = stratumn[q1];
    NumericVector time1 = timen[q1];
    NumericVector time21 = time2n[q1];
    IntegerVector event1 = eventn[q1];
    NumericVector weight1 = weightn[q1];
    IntegerVector id1 = idn[q1];

    NumericMatrix z1(n1,nvar);
    for (i=0; i<n1; i++) {
      for (j=0; j<nvar; j++) {
        z1(i,j) = zn(q1[i],j);
      }
    }

    // unify right censored data with interval censored data
    NumericVector tstart(n1), tstop(n1);
    if (!has_time2) {
      tstart = time1;
      for (i=0; i<n1; i++) {
        tstop[i] = event1[i] == 1 ? tstart[i] : NA_REAL;
      }
    } else {
      tstart = time1;
      tstop = time21;
    }

    IntegerVector status(n1);
    for (i=0; i<n1; i++) {
      if (!R_isnancpp(tstart[i]) && !R_isnancpp(tstop[i]) &&
          (tstart[i] == tstop[i])) {
        status[i] = 1; // event
      } else if (!R_isnancpp(tstart[i]) && !R_isnancpp(tstop[i]) &&
        (tstart[i] < tstop[i])) {
        status[i] = 3; // interval censoring
      } else if (R_isnancpp(tstart[i]) && !R_isnancpp(tstop[i])) {
        status[i] = 2; // left censoring
      } else if (!R_isnancpp(tstart[i]) && R_isnancpp(tstop[i])) {
        status[i] = 0; // right censoring
      } else {
        status[i] = -1; // exclude the observation
      }
    }

    nobs[h] = n1;
    nevents[h] = sum(status == 1);

    // initial parameter values
    NumericVector time0(n1);
    for (i=0; i<n1; i++) {
      if (status[i] == 1) { // event
        time0[i] = tstart[i];
      } else if (status[i] == 3) { // interval censoring
        time0[i] = (tstart[i] + tstop[i])/2;
      } else if (status[i] == 2) { // left censoring
        time0[i] = tstop[i];
      } else if (status[i] == 0) { // right censoring
        time0[i] = tstart[i];
      } else {
        time0[i] = NA_REAL;
      }
    }

    // intercept only model
    LogicalVector sub = !is_na(time0);
    NumericVector y0 = time0[sub];
    if ((dist1 == "exponential") || (dist1 == "weibull") ||
        (dist1 == "lognormal") || (dist1 == "loglogistic")) {
      y0 = log(y0);
    }

    double int0 = mean(y0);
    double logsig0 = log(sd(y0));

    int pint = dist1 == "exponential" ? 1 : (nstrata+1);
    NumericVector bint0(pint);
    if (dist1 == "exponential") {
      bint0[0] = int0;
    } else {
      bint0[0] = int0;
      for (i=0; i<nstrata; i++) {
        bint0[i+1] = logsig0;
      }
    }

    NumericMatrix zi(n1,1);
    for (i=0; i<n1; i++) zi(i,0) = 1;

    // parameter estimates and standard errors for the null model
    aftparams parami = {dist1, stratum1, tstart, tstop, status, weight1, zi};
    List outint = bmini(bint0, f_nllik_1, f_nscore_1, &parami, 1e-9);
    NumericVector bint = outint["par"];
    std::vector<double> parbint(bint.begin(), bint.end());

    NumericVector b0(p);
    if (dist1 == "exponential") {
      b0[0] = bint[0];
    } else {
      b0[0] = bint[0];
      for (i=0; i<nstrata; i++) {
        b0[nvar+i] = bint[i+1];
      }
    }

    // exclude observations with missing covariates
    for (i=0; i<n1; i++) sub[i] = is_false(any(is_na(z1(i,_))));
    int nsub = sum(sub);

    if (nsub < n1) {
      stratum1 = stratum1[sub];
      tstart = tstart[sub];
      tstop = tstop[sub];
      status = status[sub];
      weight1 = weight1[sub];
      id1 = id1[sub];

      NumericMatrix z2(nsub, nvar);
      j = 0;
      for (i = 0; i < n1; i++) {
        if (sub[i]) {
          z2(j, _) = z1(i, _);
          j++;
        }
      }
      z1 = z2;
    }

    // parameter estimates and standard errors for the full model
    aftparams param = {dist1, stratum1, tstart, tstop, status, weight1, z1};
    List out = bmini(b0, f_nllik_1, f_nscore_1, &param, 1e-9);
    NumericVector b = out["par"];

    std::vector<double> parb(b.begin(), b.end());
    NumericMatrix infob = f_info_1(p, parb.data(), &param);
    NumericMatrix vb = invsympd(infob);

    NumericVector seb(p);
    for (j=0; j<p; j++) {
      seb[j] = sqrt(vb(j,j));
    }

    for (i=0; i<p; i++) {
      rep0[h*p+i] = h+1;

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
      for (j=0; j<p; j++) {
        vbeta0(h*p+i,j) = vb(i,j);
      }
    }

    // log-likelihoods and score test statistic
    std::vector<double> parb0(b0.begin(), b0.end());
    loglik(h,0) = -f_nllik_1(p, parb0.data(), &param);
    loglik(h,1) = -as<double>(out["value"]);

    NumericVector score(p);
    std::vector<double> pars(score.begin(), score.end());
    f_nscore_1(p, parb0.data(), pars.data(), &param);
    for (j=0; j<p; j++) score[j] = -pars[j];

    NumericMatrix infob0 = f_info_1(p, parb0.data(), &param);
    NumericMatrix vb0 = invsympd(infob0);
    for (i=1; i<nvar; i++) {
      for (j=1; j<nvar; j++) {
        scoretest[h] += score[i]*vb0(i,j)*score[j];
      }
    }

    // robust variance estimates
    if (robust) {
      NumericMatrix ressco = f_ressco_1(b, &param); // score residuals

      int nr; // number of rows in the score residual matrix
      if (!has_id) {
        nr = nsub;
      } else { // need to sum up score residuals by id
        IntegerVector order = seq(0, nsub-1);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return id1[i] < id1[j];
        });

        IntegerVector id2 = id1[order];
        IntegerVector idx(1,0);
        for (i=1; i<nsub; i++) {
          if (id2[i] != id2[i-1]) {
            idx.push_back(i);
          }
        }

        int nids = static_cast<int>(idx.size());
        idx.push_back(nsub);

        NumericMatrix resid(nsub,p);
        for (i=0; i<nsub; i++) {
          for (j=0; j<p; j++) {
            resid(i,j) += ressco(order[i],j);
          }
        }

        NumericMatrix ressco2(nids,p);
        for (i=0; i<nids; i++) {
          for (j=0; j<p; j++) {
            for (k=idx[i]; k<idx[i+1]; k++) {
              ressco2(i,j) += resid(k,j);
            }
          }
        }

        ressco = ressco2;  // update the score residuals
        nr = nids;
      }

      NumericMatrix D(nr,p); // DFBETA
      for (i=0; i<nr; i++) {
        for (j=0; j<p; j++) {
          for (k=0; k<p; k++) {
            D(i,j) += weight1[i]*ressco(i,k)*vb(k,j);
          }
        }
      }

      NumericMatrix rvb(p,p); // robust variance matrix for betahat
      for (j=0; j<p; j++) {
        for (k=0; k<p; k++) {
          for (i=0; i<nr; i++) {
            rvb(j,k) += D(i,j)*D(i,k);
          }
        }
      }

      NumericVector rseb(p);  // robust standard error for betahat
      for (i=0; i<p; i++) rseb[i] = sqrt(rvb(i,i));

      for (i=0; i<p; i++) {
        rsebeta0[h*p+i] = rseb[i];
        for (j=0; j<p; j++) {
          rvbeta0(h*p+i,j) = rvb(i,j);
        }
      }
    }
  }

  NumericVector expbeta0 = exp(beta0);
  NumericVector z0(nreps*p);
  if (!robust) z0 = beta0/sebeta0;
  else z0 = beta0/rsebeta0;

  DataFrame sumstat = List::create(
    _["n"] = nobs,
    _["nevents"] = nevents,
    _["loglik0"] = loglik(_,0),
    _["loglik1"] = loglik(_,1),
    _["scoretest"] = scoretest);

  DataFrame parest;
  if (!robust) {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = sebeta0,
      _["z"] = z0,
      _["expbeta"] = expbeta0,
      _["vbeta"] = vbeta0);
  } else {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = sebeta0,
      _["rsebeta"] = rsebeta0,
      _["z"] = z0,
      _["expbeta"] = expbeta0,
      _["vbeta"] = vbeta0,
      _["rvbeta"] = rvbeta0);
  }

  if (has_rep) {
    if (TYPEOF(data[rep]) == INTSXP) {
      sumstat.push_back(repwi[rep01-1], rep);
      parest.push_back(repwi[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      sumstat.push_back(repwn[rep01-1], rep);
      parest.push_back(repwn[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      sumstat.push_back(repwc[rep01-1], rep);
      parest.push_back(repwc[rep0-1], rep);
    }
  }

  List result = List::create(
    _["sumstat"] = sumstat,
    _["parest"] = parest);

  return result;
}


// define functions in likelihood inference, algorithms adapted from coxph

// negative log likelihood
double f_nllik_2(int p, double *par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, k, person;

  NumericVector beta(p);
  for (i=0; i<p; i++) beta[i] = par[i];

  double loglik = 0;        // log-likelihood value
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double zbeta;             // linear predictor
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    zbeta = 0;
    for (i=0; i<p; i++) zbeta += beta[i]*param->z(person,i);
    eta[person] = zbeta;
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      denom = 0;
    }

    dtime = param->tstop[person];
    while (person < param->nused && param->tstop[person] == dtime) {
      // walk through this set of tied times
      risk = param->weight[person]*exp(eta[person]);
      if (param->event[person] == 0) {
        denom += risk;
      } else {
        ndead++;
        deadwt += param->weight[person];
        denom2 += risk;
        loglik += param->weight[person]*eta[person];
      }

      person++;

      if (person < param->nused && param->strata[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1 < param->nused; i1++) {
      p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      risk = param->weight[p1]*exp(eta[p1]);
      denom -= risk;
    }

    // add to the main terms
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        loglik -= deadwt*log(denom);
      } else {
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          loglik -= meanwt*log(denom);
        }
      }

      // reset for the next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
    }
  }

  return -loglik;
}


// negative score vector
void f_nscore_2(int p, double *par, double *gr, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, k, person;

  NumericVector beta(p);
  for (i=0; i<p; i++) beta[i] = par[i];

  NumericVector u(p);       // score vector
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double zbeta;             // linear predictor
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;              // zbar(beta,k,t)

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    zbeta = 0;
    for (i=0; i<p; i++) zbeta += beta[i]*param->z(person,i);
    eta[person] = zbeta;
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      denom = 0;
      for (i=0; i<p; i++) a[i] = 0;
    }

    dtime = param->tstop[person];
    while (person < param->nused && param->tstop[person] == dtime) {
      // walk through this set of tied times
      risk = param->weight[person]*exp(eta[person]);
      if (param->event[person] == 0) {
        denom += risk;
        for (i=0; i<p; i++) a[i] += risk*param->z(person,i);
      } else {
        ndead++;
        deadwt += param->weight[person];
        denom2 += risk;
        for (i=0; i<p; i++) {
          a2[i] += risk*param->z(person,i);
          u[i] += param->weight[person]*param->z(person,i);
        }
      }

      person++;

      if (person < param->nused && param->strata[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1<param->nused; i1++) {
      p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      risk = param->weight[p1]*exp(eta[p1]);
      denom -= risk;
      for (i=0; i<p; i++) a[i] -= risk*param->z(p1,i);
    }

    // add to the main terms
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        for (i=0; i<p; i++) {
          a[i] += a2[i];
          xbar = a[i]/denom;
          u[i] -= deadwt*xbar;
        }
      } else {
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            xbar = a[i]/denom;
            u[i] -= meanwt*xbar;
          }
        }
      }

      // reset for the next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
      for (i=0; i<p; i++) a2[i] = 0;
    }
  }

  for (i=0; i<p; i++) gr[i] = -u[i];
}


// observed information matrix
NumericMatrix f_info_2(int p, double *par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, j, k, person;

  NumericVector beta(p);
  for (i=0; i<p; i++) beta[i] = par[i];

  NumericMatrix imat(p,p);  // information matrix
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double zbeta;             // linear predictor
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;              // zbar(beta,k,t)
  NumericMatrix cmat(p,p);  // s2(beta,k,t)
  NumericMatrix cmat2(p,p); // sum of w*exp(zbeta)*z*z' for the deaths

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    zbeta = 0;
    for (i=0; i<p; i++) zbeta += beta[i]*param->z(person,i);
    eta[person] = zbeta;
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) {  // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      denom = 0;
      for (i=0; i<p; i++) {
        a[i] = 0;
        for (j=0; j<p; j++) cmat(i,j) = 0;
      }
    }

    dtime = param->tstop[person];
    while (person < param->nused && param->tstop[person] == dtime) {
      // walk through this set of tied times
      risk = param->weight[person]*exp(eta[person]);
      if (param->event[person] == 0) {
        denom += risk;
        for (i=0; i<p; i++) {
          a[i] += risk*param->z(person,i);
          for (j=0; j<=i; j++) {
            cmat(i,j) += risk*param->z(person,i)*param->z(person,j);
          }
        }
      } else {
        ndead++;
        deadwt += param->weight[person];
        denom2 += risk;
        for (i=0; i<p; i++) {
          a2[i] += risk*param->z(person,i);
          for (j=0; j<=i; j++) {
            cmat2(i,j) += risk*param->z(person,i)*param->z(person,j);
          }
        }
      }

      person++;

      if (person < param->nused && param->strata[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1<param->nused; i1++) {
      p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      risk = param->weight[p1]*exp(eta[p1]);
      denom -= risk;
      for (i=0; i<p; i++) {
        a[i] -= risk*param->z(p1,i);
        for (j=0; j<=i; j++) {
          cmat(i,j) -= risk*param->z(p1,i)*param->z(p1,j);
        }
      }
    }

    // add to the main terms
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        for (i=0; i<p; i++) {
          a[i] += a2[i];
          xbar = a[i]/denom;
          for (j=0; j<=i; j++) {
            cmat(i,j) += cmat2(i,j);
            imat(i,j) += deadwt*(cmat(i,j) - xbar*a[j])/denom;
          }
        }
      } else {
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            xbar = a[i]/denom;
            for (j=0; j<=i; j++) {
              cmat(i,j) += cmat2(i,j)/ndead;
              imat(i,j) += meanwt*(cmat(i,j) - xbar*a[j])/denom;
            }
          }
        }
      }

      // reset for next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
      for (i=0; i<p; i++) {
        a2[i] = 0;
        for (j=0; j<=i; j++) cmat2(i,j) = 0;
      }
    }
  }

  for (i=0; i<p-1; i++) {
    for (j=i+1; j<p; j++) {
      imat(i,j) = imat(j,i);
    }
  }

  return imat;
}


// score residual matrix
NumericMatrix f_ressco_2(NumericVector beta, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, j, k, person, n = static_cast<int>(param->tstart.size());
  int p = static_cast<int>(beta.size());

  NumericMatrix resid(n,p);
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double zbeta;             // linear predictor
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;              // zbar(beta,k,t)

  double downwt, hazard, cumhaz = 0;
  NumericVector xhaz(p), mh1(p), mh2(p), mh3(p);

  NumericVector score(param->nused);
  for (person = 0; person < param->nused; person++) {
    zbeta = 0;
    for (i=0; i<p; i++) zbeta += beta[i]*param->z(person,i);
    score[person] = exp(zbeta);
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person];
      // first obs of a new stratum, finish off the prior stratum
      for (; i1 < param->nused && param->order1[i1] < person; i1++) {
        p1 = param->order1[i1];
        for (i=0; i<p; i++) {
          resid(p1,i) -= score[p1]*(param->z(p1,i)*cumhaz - xhaz[i]);
        }
      }
      denom = 0; // reset temporary variables
      cumhaz = 0;
      for (i=0; i<p; i++) {
        a[i] =0;
        xhaz[i] =0;
      }
    }

    dtime = param->tstop[person];
    while (person < param->nused && param->tstop[person] == dtime) {
      // walk through this set of tied times
      // initialize residuals to score[i]*(x[i]*cumhaz - xhaz), before
      // updating cumhaz and xhaz
      for (i=0; i<p; i++) {
        resid(person,i) = score[person]*(param->z(person,i)*cumhaz - xhaz[i]);
      }

      risk = param->weight[person]*score[person];
      if (param->event[person] == 0) {
        denom += risk;
        for (i=0; i<p; i++) a[i] += risk*param->z(person,i);
      } else {
        ndead++;
        deadwt += param->weight[person];
        denom2 += risk;
        for (i=0; i<p; i++) {
          a2[i] += risk*param->z(person,i);
        }
      }

      person++;

      if (person < param->nused && param->strata[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1<param->nused; i1++) {
      p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      risk = param->weight[p1]*score[p1];
      denom -= risk;
      for (i=0; i<p; i++) {
        // finish the residual by subtracting score[i]*(x[i]*cumhaz - xhaz)
        resid(p1,i) -= score[p1]*(param->z(p1,i)*cumhaz - xhaz[i]);
        a[i] -= risk*param->z(p1,i);
      }
    }

    // update the cumulative sums at death times
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        hazard = deadwt/denom;
        cumhaz += hazard;
        for (i=0; i<p; i++) {
          a[i] += a2[i];
          xbar = a[i]/denom;
          xhaz[i] += xbar*hazard;
          for (j=person-1; j>=person-ndead; j--) {
            resid(j,i) += param->z(j,i) - xbar;
          }
        }
      } else {
        for (i=0; i<p; i++) {
          mh1[i] = 0;
          mh2[i] = 0;
          mh3[i] = 0;
        }
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          hazard = meanwt/denom;
          cumhaz += hazard;
          downwt = (ndead-k-1.0)/ndead;
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            xbar = a[i]/denom;
            xhaz[i] += xbar*hazard;
            mh1[i]  += hazard*downwt;
            mh2[i]  += xbar*hazard*downwt;
            mh3[i]  += xbar/ndead;
          }
        }

        for (j=person-1; j>=person-ndead; j--) {
          for (i=0; i<p; i++) {
            resid(j,i) += (param->z(j,i) - mh3[i]) +
              score[j]*(param->z(j,i)*mh1[i] - mh2[i]);
          }
        }
      }

      // reset for next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
      for (i=0; i<p; i++) a2[i] = 0;
    }
  }

  // finish those remaining in the final stratum
  for (; i1<param->nused; i1++) {
    p1 = param->order1[i1];
    for (i=0; i<p; i++)
      resid(p1,i) -= score[p1]*(param->z(p1,i)*cumhaz - xhaz[i]);
  }

  return resid;
}


//' @title Proportional hazards regression model
//' @description Obtains the hazard ratio estimates from the proportional
//' hazards regression model with right censored or counting process data.
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
//'     data only. Intervals are assumed to be open on the left
//'     and closed on the right, and event indicates whether an event
//'     occurred at the right end of each interval.
//'
//'   * \code{event}: The event indicator, normally 1=event, 0=no event.
//'
//'   * \code{covariates}: The values of baseline covariates (and
//'     time-dependent covariates in each interval for counting
//'     process data). This is the full-rank design matrix for the Cox
//'     model, assuming that factor variables have already been
//'     expanded into dummy variables.
//'
//'   * \code{weight}: The weight for each observation.
//'
//'   * \code{id}: The optional subject ID for counting process data
//'     with time-dependent covariates.
//'
//' @param rep The name of the replication variable in the input data.
//' @param stratum The name of the stratum variable in the input data.
//' @param time The name of the time variable or the left end of each
//'   interval for counting process data in the input data.
//' @param time2 The name of the right end of each interval for counting
//'   process data in the input data.
//' @param event The name of the event variable in the input data.
//' @param covariates The vector of names of baseline and time-dependent
//'   covariates in the input data.
//' @param weight The name of the weighting variable in the input data.
//' @param id The name of the id variable in the input data.
//' @param ties The method for handling ties with options including
//'   "breslow" and "efron" (default).
//' @param robust Whether a robust sandwich variance estimate should be
//'   computed. The default is TRUE if there are fractional weights or
//'   there is at least 1 id with >1 event. In the presence of the id
//'   variable, the score residual will be aggregated for each id when
//'   computing the robust sandwich variance estimate.
//'
//' @return A list with the following components:
//'
//' * \code{sumstat}: The data frame of summary statistics of model fit
//'   with the following variables:
//'
//'     - \code{rep}: The replication.
//'
//'     - \code{n}: The number of observations.
//'
//'     - \code{nevents}: The number of events.
//'
//'     - \code{loglik0}: The log-likelihood under null.
//'
//'     - \code{loglik1}: The maximum log-likelihood.
//'
//'     - \code{scoretest}: The score test statistic.
//'
//' * \code{parest}: The data frame of parameter estimates with the
//'   following variables:
//'
//'     - \code{rep}: The replication.
//'
//'     - \code{param}: The name of the covariate for the parameter estimate.
//'
//'     - \code{beta}: The log hazard ratio estimate.
//'
//'     - \code{sebeta}: The standard error of log hazard ratio estimate.
//'
//'     - \code{rsebeta}: The robust standard error of log hazard ratio
//'       estimate if robust variance is requested.
//'
//'     - \code{z}: The Wald test statistic for log hazard ratio. The
//'       \code{rsebeta} will be used if robust variance is requested.
//'
//'     - \code{hazardRatio}: The hazard ratio estimate.
//'
//'     - \code{vbeta}: The covariance matrix for parameter estimates.
//'
//'     - \code{rvbeta}: The robust covariance matrix for parameter
//'       estimates if robust variance is requested.
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' library(dplyr)
//'
//' # Example 1 with right-censored data
//' phregr(data = rawdata %>% mutate(treat = 1*(treatmentGroup == 1)),
//'        rep = "iterationNumber", stratum = "stratum",
//'        time = "timeUnderObservation", event = "event",
//'        covariates = "treat")
//'
//' # Example 2 with counting process data and robust variance estimate
//' phregr(data = heart %>% mutate(rx = as.numeric(transplant) - 1),
//'        time = "start", time2 = "stop", event = "event",
//'        covariates = c("rx", "age"), id = "id", robust = 1)
//'
//' @export
// [[Rcpp::export]]
List phregr(const DataFrame data,
            const std::string rep = "rep",
            const std::string stratum = "stratum",
            const std::string time = "time",
            const std::string time2 = "time2",
            const std::string event = "event",
            const StringVector& covariates = "treat",
            const std::string weight = "weight",
            const std::string id = "id",
            const std::string ties = "efron",
            bool robust = 0) {

  int h, i, j, k, n = data.nrows(), p = static_cast<int>(covariates.size());

  bool has_rep = hasVariable(data, rep);
  bool has_stratum = hasVariable(data, stratum);

  bool has_time = hasVariable(data, time);
  if (!has_time) {
    stop("data must contain the time variable");
  }

  NumericVector timen = data[time];
  if (is_true(any(timen < 0))) {
    stop("time must be nonnegative for each subject");
  }

  bool has_time2 = hasVariable(data, time2);

  NumericVector time2n(n);
  if (has_time2) {
    time2n = data[time2];
    if (is_true(any(time2n <= timen))) {
      stop("time2 must be greater than time for each observation");
    }
  }

  bool has_event = hasVariable(data, event);
  if (!has_event) {
    stop("data must contain the event variable");
  }

  IntegerVector eventn = data[event];
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0 for each subject");
  }

  if (is_true(all(eventn == 0))) {
    stop("at least 1 event is needed to fit the Cox model");
  }

  NumericMatrix zn(n,p);
  for (j=0; j<p; j++) {
    String zj = covariates[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in covariates");
    }
    NumericVector u = data[zj];
    for (i=0; i<n; i++) zn(i,j) = u[i];
  }

  bool has_weight = hasVariable(data, weight);

  NumericVector weightn(n, 1.0);
  if (has_weight) {
    weightn = data[weight];
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }

  bool has_id = hasVariable(data, id);

  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  int method = meth == "efron" ? 1 : 0;

  // create the numeric rep variable
  IntegerVector repn(n);
  IntegerVector repwi;
  NumericVector repwn;
  StringVector repwc;
  if (!has_rep) {
    repn.fill(1);
  } else {
    if (TYPEOF(data[rep]) == INTSXP) {
      IntegerVector repv = data[rep];
      repwi = unique(repv);
      repwi.sort();
      repn = match(repv, repwi);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      NumericVector repv = data[rep];
      repwn = unique(repv);
      repwn.sort();
      repn = match(repv, repwn);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      StringVector repv = data[rep];
      repwc = unique(repv);
      repwc.sort();
      repn = match(repv, repwc);
    } else {
      stop("incorrect type for the rep variable in the input data");
    }
  }

  // create the numeric stratum variable
  IntegerVector stratumn(n);
  if (!has_stratum) {
    stratumn.fill(1);
  } else {
    if (TYPEOF(data[stratum]) == INTSXP) {
      IntegerVector stratumv = data[stratum];
      IntegerVector stratumwi = unique(stratumv);
      stratumwi.sort();
      stratumn = match(stratumv, stratumwi);
    } else if (TYPEOF(data[stratum]) == REALSXP) {
      NumericVector stratumv = data[stratum];
      NumericVector stratumwn = unique(stratumv);
      stratumwn.sort();
      stratumn = match(stratumv, stratumwn);
    } else if (TYPEOF(data[stratum]) == STRSXP) {
      StringVector stratumv = data[stratum];
      StringVector stratumwc = unique(stratumv);
      stratumwc.sort();
      stratumn = match(stratumv, stratumwc);
    } else {
      stop("incorrect type for the stratum variable in the input data");
    }
  }

  // create the numeric id variable
  IntegerVector idn(n);
  if (!has_id) {
    idn = seq(1,n);
  } else {
    if (TYPEOF(data[id]) == INTSXP) {
      IntegerVector idv = data[id];
      IntegerVector idwi = unique(idv);
      idwi.sort();
      idn = match(idv, idwi);
    } else if (TYPEOF(data[id]) == REALSXP) {
      NumericVector idv = data[id];
      NumericVector idwn = unique(idv);
      idwn.sort();
      idn = match(idv, idwn);
    } else if (TYPEOF(data[id]) == STRSXP) {
      StringVector idv = data[id];
      StringVector idwc = unique(idv);
      idwc.sort();
      idn = match(idv, idwc);
    } else {
      stop("incorrect type for the id variable in the input data");
    }
  }

  // check if there is as least one id with more than 1 event
  IntegerVector idn0 = idn[eventn==1];
  IntegerVector idw0 = unique(idn0);
  bool dup_id = idw0.size() < idn0.size();

  if (R_isnancpp(robust)) {
    if ((has_id && dup_id) || is_true(any(weightn != floor(weightn)))) {
      robust = 1;
    } else {
      robust = 0;
    }
  }

  if (robust && has_time2 && !has_id) {
    stop("id is needed for counting process data with robust variance");
  }

  // create an unique id variable to recover the order of input data set
  IntegerVector uid = seq(1,n);

  // sort the data by rep
  IntegerVector order = seq(0, n-1);
  std::sort(order.begin(), order.end(), [&](int i, int j) {
    return repn[i] < repn[j];
  });

  repn = repn[order];
  stratumn = stratumn[order];
  timen = timen[order];
  time2n = time2n[order];
  eventn = eventn[order];
  weightn = weightn[order];
  idn = idn[order];
  uid = uid[order];

  NumericMatrix z(n,p);
  for (i=0; i<n; i++) {
    for (j=0; j<p; j++) {
      z(i,j) = zn(order[i], j);
    }
  }
  zn = z;

  // identify the locations of the unique values of rep
  IntegerVector idx(1,0);
  for (i=1; i<n; i++) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }

  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);

  // variables in the output data sets
  IntegerVector rep01 = seq(1,nreps);
  IntegerVector nobs(nreps), nevents(nreps);
  NumericMatrix loglik(nreps,2);
  NumericVector scoretest(nreps);

  IntegerVector rep0(nreps*p);
  StringVector par0(nreps*p);
  NumericVector beta0(nreps*p), sebeta0(nreps*p), rsebeta0(nreps*p);
  NumericMatrix vbeta0(nreps*p,p), rvbeta0(nreps*p,p);

  for (h=0; h<nreps; h++) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());

    IntegerVector stratum1 = stratumn[q1];
    NumericVector time1 = timen[q1];
    NumericVector time21 = time2n[q1];
    IntegerVector event1 = eventn[q1];
    NumericVector weight1 = weightn[q1];
    IntegerVector id1 = idn[q1];
    IntegerVector uid1 = uid[q1];

    NumericMatrix z1(n1,p);
    for (i=0; i<n1; i++) {
      for (j=0; j<p; j++) {
        z1(i,j) = zn(q1[i],j);
      }
    }

    nobs[h] = n1;
    nevents[h] = sum(event1);

    // unify right censored data with counting process data
    NumericVector tstart(n1), tstop(n1);
    if (!has_time2) {
      tstop = time1;
    } else {
      tstart = time1;
      tstop = time21;
    }

    // ignore subjects not at risk for any event time
    double delta = max(tstop) + 1.0; // ensure no overlap between strata
    for (i=0; i<n1; i++) {
      tstart[i] = tstart[i] + stratum1[i]*delta;
      tstop[i] = tstop[i] + stratum1[i]*delta;
    }

    NumericVector etime = tstop[event1==1];
    etime = unique(etime);
    etime.sort();

    IntegerVector index1 = findInterval3(tstart, etime);
    IntegerVector index2 = findInterval3(tstop, etime);
    IntegerVector ignore1(n1);
    for (i=0; i<n1; i++) {
      if ((index1[i] == index2[i]) || is_true(any(is_na(z1(i,_))))) {
        ignore1[i] = 1;
      } else {
        ignore1[i] = 0;
      }
    }

    int nused = n1 - sum(ignore1);

    // sort by descending stopping time within each stratum
    IntegerVector order2 = seq(0, n1-1);
    std::sort(order2.begin(), order2.end(), [&](int i, int j) {
      return (ignore1[i] < ignore1[j]) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] < stratum1[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstop[i] > tstop[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstop[i] == tstop[j]) && (event1[i] < event1[j]));
    });

    stratum1 = stratum1[order2];
    tstart = tstart[order2];
    tstop = tstop[order2];
    event1 = event1[order2];
    weight1 = weight1[order2];
    id1 = id1[order2];
    uid1 = uid1[order2];
    ignore1 = ignore1[order2];

    NumericMatrix zy(n1,p);
    for (i=0; i<n1; i++) {
      for (j=0; j<p; j++) {
        zy(i,j) = z1(order2[i],j);
      }
    }
    z1 = zy;

    // sort by starting time in descending order within each stratum
    IntegerVector order1 = seq(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      return (ignore1[i] < ignore1[j]) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] < stratum1[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstart[i] > tstart[j]));
    });


    // parameter estimates and standard errors
    coxparams param = {nused, stratum1, tstart, tstop, event1,
                       weight1, z1, order1, method};

    NumericVector b0(p);
    List out = bmini(b0, f_nllik_2, f_nscore_2, &param, 1e-9);
    NumericVector b = out["par"];

    std::vector<double> parb(b.begin(), b.end());
    NumericMatrix infob = f_info_2(p, parb.data(), &param);
    NumericMatrix vb = invsympd(infob);

    NumericVector seb(p);
    for (j=0; j<p; j++) seb[j] = sqrt(vb(j,j));

    for (i=0; i<p; i++) {
      rep0[h*p+i] = h+1;
      par0[h*p+i] = covariates[i];
      beta0[h*p+i] = b[i];
      sebeta0[h*p+i] = seb[i];
      for (j=0; j<p; j++) {
        vbeta0(h*p+i,j) = vb(i,j);
      }
    }

    // log-likelihoods and score test statistic
    std::vector<double> parb0(b0.begin(), b0.end());
    loglik(h,0) = -f_nllik_2(p, parb0.data(), &param);
    loglik(h,1) = -as<double>(out["value"]);

    NumericVector score(p);
    std::vector<double> pars(score.begin(), score.end());
    f_nscore_2(p, parb0.data(), pars.data(), &param);
    for (j=0; j<p; j++) score[j] = -pars[j];

    NumericMatrix infob0 = f_info_2(p, parb0.data(), &param);
    NumericMatrix vb0 = invsympd(infob0);
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        scoretest[h] += score[i]*vb0(i,j)*score[j];
      }
    }

    // robust variance estimates
    if (robust) {
      NumericMatrix ressco = f_ressco_2(b, &param); // score residuals

      int nr; // number of rows in the score residual matrix
      if (!has_id) {
        nr = n1;
      } else { // need to sum up score residuals by id
        IntegerVector order = seq(0, n1-1);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
          return id1[i] < id1[j];
        });

        IntegerVector id2 = id1[order];
        IntegerVector idx(1,0);
        for (i=1; i<n1; i++) {
          if (id2[i] != id2[i-1]) {
            idx.push_back(i);
          }
        }

        int nids = static_cast<int>(idx.size());
        idx.push_back(n1);

        NumericMatrix resid(n1,p);
        for (i=0; i<n1; i++) {
          for (j=0; j<p; j++) {
            resid(i,j) += ressco(order[i],j);
          }
        }

        NumericMatrix ressco2(nids,p);
        for (i=0; i<nids; i++) {
          for (j=0; j<p; j++) {
            for (k=idx[i]; k<idx[i+1]; k++) {
              ressco2(i,j) += resid(k,j);
            }
          }
        }

        ressco = ressco2;  // update the score residuals
        nr = nids;
      }

      NumericMatrix D(nr,p); // DFBETA
      for (i=0; i<nr; i++) {
        for (j=0; j<p; j++) {
          for (k=0; k<p; k++) {
            D(i,j) += weight1[i]*ressco(i,k)*vb(k,j);
          }
        }
      }

      NumericMatrix rvb(p,p); // robust variance matrix for betahat
      for (j=0; j<p; j++) {
        for (k=0; k<p; k++) {
          for (i=0; i<nr; i++) {
            rvb(j,k) += D(i,j)*D(i,k);
          }
        }
      }

      NumericVector rseb(p);  // robust standard error for betahat
      for (i=0; i<p; i++) rseb[i] = sqrt(rvb(i,i));

      for (i=0; i<p; i++) {
        rsebeta0[h*p+i] = rseb[i];
        for (j=0; j<p; j++) {
          rvbeta0(h*p+i,j) = rvb(i,j);
        }
      }
    }
  }

  NumericVector hazardRatio0 = exp(beta0);
  NumericVector z0(nreps*p);
  if (!robust) z0 = beta0/sebeta0;
  else z0 = beta0/rsebeta0;

  DataFrame sumstat = List::create(
    _["n"] = nobs,
    _["nevents"] = nevents,
    _["loglik0"] = loglik(_,0),
    _["loglik1"] = loglik(_,1),
    _["scoretest"] = scoretest);

  DataFrame parest;
  if (!robust) {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = sebeta0,
      _["z"] = z0,
      _["hazardRatio"] = hazardRatio0,
      _["vbeta"] = vbeta0);
  } else {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = sebeta0,
      _["rsebeta"] = rsebeta0,
      _["z"] = z0,
      _["hazardRatio"] = hazardRatio0,
      _["vbeta"] = vbeta0,
      _["rvbeta"] = rvbeta0);
  }

  if (has_rep) {
    if (TYPEOF(data[rep]) == INTSXP) {
      sumstat.push_back(repwi[rep01-1], rep);
      parest.push_back(repwi[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == REALSXP) {
      sumstat.push_back(repwn[rep01-1], rep);
      parest.push_back(repwn[rep0-1], rep);
    } else if (TYPEOF(data[rep]) == STRSXP) {
      sumstat.push_back(repwc[rep01-1], rep);
      parest.push_back(repwc[rep0-1], rep);
    }
  }

  List result = List::create(
      _["sumstat"] = sumstat,
      _["parest"] = parest);

  return result;
}
