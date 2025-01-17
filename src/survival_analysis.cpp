#include <Rcpp.h>
#include <R_ext/Applic.h>
#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


NumericVector fsurvci(double surv, double sesurv, String ct, double z) {
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
//'   * \code{time}: The possibly right-censored survival time.
//'
//'   * \code{event}: The event indicator.
//'
//' @param rep The name(s) of the replication variable(s) in the input data.
//' @param stratum The name(s) of the stratum variable(s) in the input data.
//' @param time The name of the time variable in the input data.
//' @param event The name of the event variable in the input data.
//' @param conftype The type of the confidence interval. One of "none",
//'   "plain", "log", "log-log" (the default), or "arcsin".
//'   The arcsin option bases the intervals on asin(sqrt(survival)).
//' @param conflev The level of the two-sided confidence interval for
//'   the survival probabilities. Defaults to 0.95.
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
//' * \code{survival}: The Kaplan-Meier estimate of the survival probability.
//'
//' * \code{stderr}: The standard error of the estimated survival
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
                const std::string event = "event",
                const std::string conftype = "log-log",
                const double conflev = 0.95) {
  int h, i, j, n = data.nrows();

  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(1);
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
    stratumn.fill(1);
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
  double z = R::qnorm((1.0 + conflev)/2.0, 0, 1, 1, 0);

  // sort the data by rep
  IntegerVector order = seq(0, n-1);
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
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
    IntegerVector event1 = eventn[q1];

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
      IntegerVector event2 = event1[q2];

      int s = stratum1[q2[0]], n2 = static_cast<int>(q2.size());
      double t = 0, nrisk = n2, nevent = 0, surv = 1, vcumhaz = 0, sesurv;
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
              NumericVector ci = fsurvci(surv, sesurv, ct, z);
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
            NumericVector ci = fsurvci(surv, sesurv, ct, z);
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
          NumericVector ci = fsurvci(surv, sesurv, ct, z);
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
    result.push_back(conflev, "conflev");
    result.push_back(conftype, "conftype");
  }

  if (has_stratum) {
    for (i=0; i<p_stratum; i++) {
      String s = stratum[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector stratumwi = u_stratum[s];
        result.push_back(stratumwi[stratum0-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector stratumwn = u_stratum[s];
        result.push_back(stratumwn[stratum0-1], s);
      } else if (TYPEOF(data[s]) == STRSXP) {
        StringVector stratumwc = u_stratum[s];
        result.push_back(stratumwc[stratum0-1], s);
      }
    }
  }

  if (has_rep) {
    for (i=0; i<p_rep; i++) {
      String s = rep[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector repwi = u_rep[s];
        result.push_back(repwi[rep0-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector repwn = u_rep[s];
        result.push_back(repwn[rep0-1], s);
      } else if (TYPEOF(data[rep]) == STRSXP) {
        StringVector repwc = u_rep[s];
        result.push_back(repwc[rep0-1], s);
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
//' * \code{vsurvDiff}: The variance for survDiff.
//'
//' * \code{survDiffZ}: The Z-statistic value.
//'
//' * \code{survDiffPValue}: The one-sided p-value.
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
  int h, i, j, k, n = data.nrows();

  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(1);
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
    stratumn.fill(1);
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
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    IntegerVector treatv = data[treat];
    treatwi = unique(treatv);
    if (treatwi.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwi == 0) | (treatwi == 1)))) {
      treatwi = IntegerVector::create(1,0);
      treatn = 2 - treatv;
    } else {
      treatwi.sort();
      treatn = match(treatv, treatwi);
    }
  } else if (TYPEOF(data[treat]) == REALSXP) {
    NumericVector treatv = data[treat];
    treatwn = unique(treatv);
    if (treatwn.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwn == 0) | (treatwn == 1)))) {
      treatwn = NumericVector::create(1,0);
      treatn = 2 - as<IntegerVector>(treatv);
    } else {
      treatwn.sort();
      treatn = match(treatv, treatwn);
    }
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
  IntegerVector order = seq(0, n-1);
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
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
  NumericVector surv10(nreps), surv20(nreps), survDiff0(nreps);
  NumericVector vsurv10(nreps), vsurv20(nreps), vsurvDiff0(nreps);
  NumericVector survDiffZ0(nreps), survDiffPValue0(nreps);
  NumericVector lower0(nreps), upper0(nreps);

  double z = R::qnorm((1.0 + conflev)/2.0, 0, 1, 1, 0);

  bool noerr = 1;
  int index = 0;
  for (h=0; h<nreps; h++) {
    bool skip = 0;
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
    for (i=1; i<n1; i++) {
      if (stratum1[i] != stratum1[i-1]) {
        idx1.push_back(i);
      }
    }

    int nstrata = static_cast<int>(idx1.size());
    idx1.push_back(n1);

    // whether the milestone exceeds the largest observed time
    for (i=0; i<nstrata; i++) {
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
          for (j=0; j<p_rep; j++) {
            std::string s = as<std::string>(rep[j]);
            if (TYPEOF(data[s]) == INTSXP) {
              IntegerVector repwi = u_rep[s];
              reperr = reperr + " " + s + " = " +
                std::to_string(repwi[repn[idx[h]]-1]);
            } else if (TYPEOF(data[s]) == REALSXP) {
              NumericVector repwn = u_rep[s];
              reperr = reperr + " " + s + " = " +
                std::to_string(repwn[repn[idx[h]]-1]);
            } else {
              StringVector repwc = u_rep[s];
              reperr = reperr + " " + s + " = " +
                repwc[repn[idx[h]]-1];
            }
          }
        }


        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else {
          for (j=0; j<p_stratum; j++) {
            std::string s = as<std::string>(stratum[j]);
            if (TYPEOF(data[s]) == INTSXP) {
              IntegerVector stratumwi = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                std::to_string(stratumwi[stratum1[idx1[i]]-1]);
            } else if (TYPEOF(data[s]) == REALSXP) {
              NumericVector stratumwn = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                std::to_string(stratumwn[stratum1[idx1[i]]-1]);
            } else {
              StringVector stratumwc = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                stratumwc[stratum1[idx1[i]]-1];
            }
          }
        }


        int k = milestone > max(time21) ? 0 : 1;
        std::string treaterr;
        if ((TYPEOF(data[treat]) == LGLSXP) ||
            (TYPEOF(data[treat]) == INTSXP)) {
          treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
        } else if (TYPEOF(data[treat]) == REALSXP) {
          treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
        } else {
          treaterr = " " + treat + " = " + treatwc[k];
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
          noerr = 0;
        }

        skip = 1;
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

    DataFrame dfout = kmest(dfin, "stratum", "treat", "time", "event",
                            "none", 0.95);

    IntegerVector stratum2 = dfout["stratum"];
    IntegerVector treat2 = dfout["treat"];
    IntegerVector treatsize = dfout["size"];
    NumericVector time2 = dfout["time"];
    NumericVector survival2 = dfout["survival"];
    NumericVector stderr2 = dfout["stderr"];

    int n2 = static_cast<int>(stratum2.size());

    // identify the locations of the unique values of stratum
    IntegerVector idx2(1,0);
    for (i=1; i<n2; i++) {
      if (stratum2[i] != stratum2[i-1]) {
        idx2.push_back(i);
      }
    }

    idx2.push_back(n2);

    IntegerVector m(nstrata, 0); // number of subjects in each stratum
    for (i=0; i<nstrata; i++) {
      int j1 = idx2[i], j2 = idx2[i+1] - 1;
      if ((treat2[j1] != 1) || (treat2[j2] != 2)) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else {
          for (j=0; j<p_rep; j++) {
            std::string s = as<std::string>(rep[j]);
            if (TYPEOF(data[s]) == INTSXP) {
              IntegerVector repwi = u_rep[s];
              reperr = reperr + " " + s + " = " +
                std::to_string(repwi[repn[idx[h]]-1]);
            } else if (TYPEOF(data[s]) == REALSXP) {
              NumericVector repwn = u_rep[s];
              reperr = reperr + " " + s + " = " +
                std::to_string(repwn[repn[idx[h]]-1]);
            } else {
              StringVector repwc = u_rep[s];
              reperr = reperr + " " + s + " = " +
                repwc[repn[idx[h]]-1];
            }
          }
        }


        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else {
          for (j=0; j<p_stratum; j++) {
            std::string s = as<std::string>(stratum[j]);
            if (TYPEOF(data[s]) == INTSXP) {
              IntegerVector stratumwi = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                std::to_string(stratumwi[stratum2[j1]-1]);
            } else if (TYPEOF(data[s]) == REALSXP) {
              NumericVector stratumwn = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                std::to_string(stratumwn[stratum2[j1]-1]);
            } else {
              StringVector stratumwc = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                stratumwc[stratum2[j1]-1];
            }
          }
        }


        int k = treat2[j1] != 1 ? 0 : 1;
        std::string treaterr;
        if ((TYPEOF(data[treat]) == LGLSXP) ||
            (TYPEOF(data[treat]) == INTSXP)) {
          treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
        } else if (TYPEOF(data[treat]) == REALSXP) {
          treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
        } else {
          treaterr = " " + treat + " = " + treatwc[k];
        }

        std::string str1 = "The data set does not contain";
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

      m[i] += treatsize[j1] + treatsize[j2];
    }

    // skip the replication if there is a stratum without both treatments
    if (skip) continue;

    double M = sum(m);
    NumericVector p(nstrata);

    double surv1 = 0.0, surv2 = 0.0, vsurv1 = 0.0, vsurv2 = 0.0;
    for (i=0; i<nstrata; i++) {
      p[i] = m[i]/M; // fraction of subjects in the stratum
      IntegerVector q = Range(idx2[i], idx2[i+1]-1);
      IntegerVector treatx = treat2[q];
      NumericVector timex = time2[q];
      NumericVector survivalx = survival2[q];
      NumericVector stderrx = stderr2[q];

      NumericVector surv(2), vsurv(2);
      for (j=0; j<2; j++) {
        LogicalVector sub = (treatx == j+1);
        NumericVector time0 = timex[sub];
        NumericVector survival0 = survivalx[sub];
        NumericVector stderr0 = stderrx[sub];
        int K = sum(sub);

        // find the latest event time before milestone for each treat
        for (k = 0; k < K; k++) {
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

    rep0[index] = repn[idx[h]];
    surv10[index] = surv1;
    surv20[index] = surv2;
    vsurv10[index] = vsurv1;
    vsurv20[index] = vsurv2;
    survDiff0[index] = surv1 - surv2;
    vsurvDiff0[index] = vsurv1 + vsurv2;
    double sesurvDiff = sqrt(vsurvDiff0[index]);
    survDiffZ0[index] = (survDiff0[index] - survDiffH0)/sesurvDiff;
    survDiffPValue0[index] = 1.0 - R::pnorm(survDiffZ0[index], 0, 1, 1, 0);
    lower0[index] = survDiff0[index] - z*sesurvDiff;
    upper0[index] = survDiff0[index] + z*sesurvDiff;

    index++;
  }

  // only keep nonmissing records
  LogicalVector sub = !is_na(rep0);
  if (is_false(any(sub))) {
    stop("no replication enables valid inference");
  }

  rep0 = rep0[sub];
  surv10 = surv10[sub];
  surv20 = surv20[sub];
  survDiff0 = survDiff0[sub];
  vsurv10 = vsurv10[sub];
  vsurv20 = vsurv20[sub];
  vsurvDiff0 = vsurvDiff0[sub];
  survDiffZ0 = survDiffZ0[sub];
  survDiffPValue0 = survDiffPValue0[sub];
  lower0 = lower0[sub];
  upper0 = upper0[sub];

  DataFrame result = DataFrame::create(
    _["milestone"] = milestone,
    _["survDiffH0"] = survDiffH0,
    _["surv1"] = surv10,
    _["surv2"] = surv20,
    _["survDiff"] = survDiff0,
    _["vsurv1"] = vsurv10,
    _["vsurv2"] = vsurv20,
    _["vsurvDiff"] = vsurvDiff0,
    _["survDiffZ"] = survDiffZ0,
    _["survDiffPValue"] = survDiffPValue0,
    _["lower"] = lower0,
    _["upper"] = upper0,
    _["conflev"] = conflev);

  if (has_rep) {
    for (i=0; i<p_rep; i++) {
      String s = rep[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector repwi = u_rep[s];
        result.push_back(repwi[rep0-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector repwn = u_rep[s];
        result.push_back(repwn[rep0-1], s);
      } else if (TYPEOF(data[rep]) == STRSXP) {
        StringVector repwc = u_rep[s];
        result.push_back(repwc[rep0-1], s);
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
//'   * \code{time}: The possibly right-censored survival time.
//'
//'   * \code{event}: The event indicator.
//'
//' @param rep The name of the replication variable in the input data.
//' @param stratum The name of the stratum variable in the input data.
//' @param treat The name of the treatment variable in the input data.
//' @param time The name of the time variable in the input data.
//' @param event The name of the event variable in the input data.
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
//' * \code{logRankPValue}: The one-sided p-value.
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
                 const std::string event = "event",
                 const double rho1 = 0,
                 const double rho2 = 0) {
  int h, i, j, k, n = data.nrows();

  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(1);
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
    stratumn.fill(1);
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
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    IntegerVector treatv = data[treat];
    treatwi = unique(treatv);
    if (treatwi.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwi == 0) | (treatwi == 1)))) {
      treatwi = IntegerVector::create(1,0);
      treatn = 2 - treatv;
    } else {
      treatwi.sort();
      treatn = match(treatv, treatwi);
    }
  } else if (TYPEOF(data[treat]) == REALSXP) {
    NumericVector treatv = data[treat];
    treatwn = unique(treatv);
    if (treatwn.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwn == 0) | (treatwn == 1)))) {
      treatwn = NumericVector::create(1,0);
      treatn = 2 - as<IntegerVector>(treatv);
    } else {
      treatwn.sort();
      treatn = match(treatv, treatwn);
    }
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


  if (rho1 < 0) stop("rho1 must be non-negative");
  if (rho2 < 0) stop("rho2 must be non-negative");

  // sort the data by rep
  IntegerVector order = seq(0, n-1);
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
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
    IntegerVector event1 = eventn[q1];

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
      IntegerVector event2 = event1[q];
      int n2 = static_cast<int>(q.size());

      // running index of number of subjects at risk
      int nriskx = n2;
      int nrisk1x = sum(treat2 == 1), nrisk2x = sum(treat2 == 2);

      if ((nrisk1x == 0) || (nrisk2x == 0)) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else {
          for (j=0; j<p_rep; j++) {
            std::string s = as<std::string>(rep[j]);
            if (TYPEOF(data[s]) == INTSXP) {
              IntegerVector repwi = u_rep[s];
              reperr = reperr + " " + s + " = " +
                std::to_string(repwi[repn[idx[h]]-1]);
            } else if (TYPEOF(data[s]) == REALSXP) {
              NumericVector repwn = u_rep[s];
              reperr = reperr + " " + s + " = " +
                std::to_string(repwn[repn[idx[h]]-1]);
            } else {
              StringVector repwc = u_rep[s];
              reperr = reperr + " " + s + " = " +
                repwc[repn[idx[h]]-1];
            }
          }
        }

        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else {
          for (j=0; j<p_stratum; j++) {
            std::string s = as<std::string>(stratum[j]);
            if (TYPEOF(data[s]) == INTSXP) {
              IntegerVector stratumwi = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                std::to_string(stratumwi[stratumn[idx[h]]-1]);
            } else if (TYPEOF(data[s]) == REALSXP) {
              NumericVector stratumwn = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                std::to_string(stratumwn[stratumn[idx[h]]-1]);
            } else {
              StringVector stratumwc = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                stratumwc[stratumn[idx[h]]-1];
            }
          }
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
      double nevent = 0, nevent1 = 0, nevent2 = 0;
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
    for (i=0; i<p_rep; i++) {
      String s = rep[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector repwi = u_rep[s];
        result.push_back(repwi[rep0-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector repwn = u_rep[s];
        result.push_back(repwn[rep0-1], s);
      } else if (TYPEOF(data[rep]) == STRSXP) {
        StringVector repwc = u_rep[s];
        result.push_back(repwc[rep0-1], s);
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
                const bool biascorrection = 0) {
  int h, i, j, k, n = data.nrows();

  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(1);
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
    stratumn.fill(1);
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
  IntegerVector order = seq(0, n-1);
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
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
  NumericVector rmst0(n), stderr0(n);
  NumericVector lower0(n), upper0(n);

  double z = R::qnorm((1.0 + conflev)/2.0, 0, 1, 1, 0);

  bool noerr = 1;
  int index = 0;
  for (h=0; h<nreps; h++) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());

    IntegerVector stratum1 = stratumn[q1];
    NumericVector time1 = timen[q1];
    IntegerVector event1 = eventn[q1];

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
      IntegerVector event2 = event1[q2];
      int n2 = static_cast<int>(q2.size());

      if (milestone > max(time2)) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else {
          for (j=0; j<p_rep; j++) {
            std::string s = as<std::string>(rep[j]);
            if (TYPEOF(data[s]) == INTSXP) {
              IntegerVector repwi = u_rep[s];
              reperr = reperr + " " + s + " = " +
                std::to_string(repwi[repn[idx[h]]-1]);
            } else if (TYPEOF(data[s]) == REALSXP) {
              NumericVector repwn = u_rep[s];
              reperr = reperr + " " + s + " = " +
                std::to_string(repwn[repn[idx[h]]-1]);
            } else {
              StringVector repwc = u_rep[s];
              reperr = reperr + " " + s + " = " +
                repwc[repn[idx[h]]-1];
            }
          }
        }


        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else {
          for (j=0; j<p_stratum; j++) {
            std::string s = as<std::string>(stratum[j]);
            if (TYPEOF(data[s]) == INTSXP) {
              IntegerVector stratumwi = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                std::to_string(stratumwi[stratum1[idx1[i]]-1]);
            } else if (TYPEOF(data[s]) == REALSXP) {
              NumericVector stratumwn = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                std::to_string(stratumwn[stratum1[idx1[i]]-1]);
            } else {
              StringVector stratumwc = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                stratumwc[stratum1[idx1[i]]-1];
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
          noerr = 0;
        }

        continue;
      }

      NumericVector time0(n2, NA_REAL);
      NumericVector nrisk0(n2), nevent0(n2), surv0(n2);
      int index1 = 0;

      double t = 0, nrisk = n2, nevent = 0, surv = 1.0;
      bool cache = 0;
      for (j=0; j<n2; j++) {
        if (((j == 0) && (event2[j] == 1)) ||
            ((j >= 1) && (event2[j] == 1) && (time2[j] > time2[j-1]))) {
          // new event
          // add the info for the previous event
          if (cache) {
            surv = surv*(1.0 - nevent/nrisk);

            time0[index1] = t;
            nrisk0[index1] = nrisk;
            nevent0[index1] = nevent;
            surv0[index1] = surv;

            index1++;
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

          time0[index1] = t;
          nrisk0[index1] = nrisk;
          nevent0[index1] = nevent;
          surv0[index1] = surv;

          index1++;

          // empty the cache for the current event time
          cache = 0;
        }
      }

      // add the info for the last event
      if (cache) {
        surv = surv*(1.0 - nevent/nrisk);

        time0[index1] = t;
        nrisk0[index1] = nrisk;
        nevent0[index1] = nevent;
        surv0[index1] = surv;

        index1++;
      }

      // only keep nonmissing records
      int N;
      LogicalVector sub = !is_na(time0);
      if (is_false(any(sub))) { // no event
        N = 0;
        time0 = NumericVector::create(0.0);
        nrisk0 = NumericVector::create(n2);
        nevent0 = NumericVector::create(0.0);
        surv0 = NumericVector::create(1.0);
      } else { // at least 1 event
        time0 = time0[sub];
        nrisk0 = nrisk0[sub];
        nevent0 = nevent0[sub];
        surv0 = surv0[sub];

        // locate the latest event time before milestone
        NumericVector milestone1(1, milestone);
        N = findInterval3(milestone1, time0)[0];

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
      for (k=1; k<=N; k++) {
        rmstx[k] = rmstx[k-1] + surv0[k]*(time0[k+1] - time0[k]);
      }

      // calculate rmst and its variance
      double u = rmstx[N];
      double v = 0.0;
      for (k=1; k<=N; k++) {
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
        for (k=1; k<=N; k++) {
          m1 += nevent0[k];
        }

        if (m1 <= 1.0) {
          std::string reperr;
          if (!has_rep) {
            reperr = "";
          } else {
            for (j=0; j<p_rep; j++) {
              std::string s = as<std::string>(rep[j]);
              if (TYPEOF(data[s]) == INTSXP) {
                IntegerVector repwi = u_rep[s];
                reperr = reperr + " " + s + " = " +
                  std::to_string(repwi[repn[idx[h]]-1]);
              } else if (TYPEOF(data[s]) == REALSXP) {
                NumericVector repwn = u_rep[s];
                reperr = reperr + " " + s + " = " +
                  std::to_string(repwn[repn[idx[h]]-1]);
              } else {
                StringVector repwc = u_rep[s];
                reperr = reperr + " " + s + " = " +
                  repwc[repn[idx[h]]-1];
              }
            }
          }


          std::string stratumerr;
          if (!has_stratum) {
            stratumerr = "";
          } else {
            for (j=0; j<p_stratum; j++) {
              std::string s = as<std::string>(stratum[j]);
              if (TYPEOF(data[s]) == INTSXP) {
                IntegerVector stratumwi = u_stratum[s];
                stratumerr = stratumerr + " " + s + " = " +
                  std::to_string(stratumwi[stratum1[idx1[i]]-1]);
              } else if (TYPEOF(data[s]) == REALSXP) {
                NumericVector stratumwn = u_stratum[s];
                stratumerr = stratumerr + " " + s + " = " +
                  std::to_string(stratumwn[stratum1[idx1[i]]-1]);
              } else {
                StringVector stratumwc = u_stratum[s];
                stratumerr = stratumerr + " " + s + " = " +
                  stratumwc[stratum1[idx1[i]]-1];
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
            noerr = 0;
          }
        } else {
          v = m1/(m1 - 1.0)*v;
        }
      }

      rep0[index] = repn[idx[h]];
      stratum0[index] = stratum1[idx1[i]];
      size0[index] = n2;
      rmst0[index] = u;
      stderr0[index] = sqrt(v);
      lower0[index] = u - z*stderr0[index];
      upper0[index] = u + z*stderr0[index];

      index++;
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
  rmst0 = rmst0[sub];
  stderr0 = stderr0[sub];
  lower0 = lower0[sub];
  upper0 = upper0[sub];


  DataFrame result = DataFrame::create(
    _["size"] = size0,
    _["milestone"] = milestone,
    _["rmst"] = rmst0,
    _["stderr"] = stderr0,
    _["lower"] = lower0,
    _["upper"] = upper0,
    _["conflev"] = conflev,
    _["biascorrection"] = biascorrection);

  if (has_stratum) {
    for (i=0; i<p_stratum; i++) {
      String s = stratum[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector stratumwi = u_stratum[s];
        result.push_back(stratumwi[stratum0-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector stratumwn = u_stratum[s];
        result.push_back(stratumwn[stratum0-1], s);
      } else if (TYPEOF(data[s]) == STRSXP) {
        StringVector stratumwc = u_stratum[s];
        result.push_back(stratumwc[stratum0-1], s);
      }
    }
  }

  if (has_rep) {
    for (i=0; i<p_rep; i++) {
      String s = rep[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector repwi = u_rep[s];
        result.push_back(repwi[rep0-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector repwn = u_rep[s];
        result.push_back(repwn[rep0-1], s);
      } else if (TYPEOF(data[rep]) == STRSXP) {
        StringVector repwc = u_rep[s];
        result.push_back(repwc[rep0-1], s);
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
//' * \code{vrmstDiff}: The variance for rmstDiff.
//'
//' * \code{rmstDiffZ}: The Z-statistic value.
//'
//' * \code{rmstDiffPValue}: The one-sided p-value.
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
                 const bool biascorrection = 0) {
  int h, i, j, n = data.nrows();

  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(1);
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
    stratumn.fill(1);
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
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    IntegerVector treatv = data[treat];
    treatwi = unique(treatv);
    if (treatwi.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwi == 0) | (treatwi == 1)))) {
      treatwi = IntegerVector::create(1,0);
      treatn = 2 - treatv;
    } else {
      treatwi.sort();
      treatn = match(treatv, treatwi);
    }
  } else if (TYPEOF(data[treat]) == REALSXP) {
    NumericVector treatv = data[treat];
    treatwn = unique(treatv);
    if (treatwn.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwn == 0) | (treatwn == 1)))) {
      treatwn = NumericVector::create(1,0);
      treatn = 2 - as<IntegerVector>(treatv);
    } else {
      treatwn.sort();
      treatn = match(treatv, treatwn);
    }
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
  IntegerVector order = seq(0, n-1);
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
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
  NumericVector rmst10(nreps), rmst20(nreps), rmstDiff0(nreps);
  NumericVector vrmst10(nreps), vrmst20(nreps), vrmstDiff0(nreps);
  NumericVector rmstDiffZ0(nreps), rmstDiffPValue0(nreps);
  NumericVector lower0(nreps), upper0(nreps);

  double z = R::qnorm((1.0 + conflev)/2.0, 0, 1, 1, 0);

  bool noerr = 1;
  int index = 0;
  for (h=0; h<nreps; h++) {
    bool skip = 0;
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
    for (i=1; i<n2; i++) {
      if (stratum2[i] != stratum2[i-1]) {
        idx2.push_back(i);
      }
    }

    int nstrata = static_cast<int>(idx2.size());
    idx2.push_back(n2);

    IntegerVector m(nstrata, 0); // number of subjects in each stratum
    for (i=0; i<nstrata; i++) {
      int j1 = idx2[i], j2 = idx2[i+1] - 1;
      if ((treat2[j1] != 1) || (treat2[j2] != 2)) {
        std::string reperr;
        if (!has_rep) {
          reperr = "";
        } else {
          for (j=0; j<p_rep; j++) {
            std::string s = as<std::string>(rep[j]);
            if (TYPEOF(data[s]) == INTSXP) {
              IntegerVector repwi = u_rep[s];
              reperr = reperr + " " + s + " = " +
                std::to_string(repwi[repn[idx[h]]-1]);
            } else if (TYPEOF(data[s]) == REALSXP) {
              NumericVector repwn = u_rep[s];
              reperr = reperr + " " + s + " = " +
                std::to_string(repwn[repn[idx[h]]-1]);
            } else {
              StringVector repwc = u_rep[s];
              reperr = reperr + " " + s + " = " +
                repwc[repn[idx[h]]-1];
            }
          }
        }

        std::string stratumerr;
        if (!has_stratum) {
          stratumerr = "";
        } else {
          for (j=0; j<p_stratum; j++) {
            std::string s = as<std::string>(stratum[j]);
            if (TYPEOF(data[s]) == INTSXP) {
              IntegerVector stratumwi = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                std::to_string(stratumwi[stratumn[idx[h]]-1]);
            } else if (TYPEOF(data[s]) == REALSXP) {
              NumericVector stratumwn = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                std::to_string(stratumwn[stratumn[idx[h]]-1]);
            } else {
              StringVector stratumwc = u_stratum[s];
              stratumerr = stratumerr + " " + s + " = " +
                stratumwc[stratumn[idx[h]]-1];
            }
          }
        }

        int k = treat2[j1] != 1 ? 0 : 1;
        std::string treaterr;
        if ((TYPEOF(data[treat]) == LGLSXP) ||
            (TYPEOF(data[treat]) == INTSXP)) {
          treaterr = " " + treat + " = " + std::to_string(treatwi[k]);
        } else if (TYPEOF(data[treat]) == REALSXP) {
          treaterr = " " + treat + " = " + std::to_string(treatwn[k]);
        } else {
          treaterr = " " + treat + " = " + treatwc[k];
        }

        std::string str1 = "The data set does not contain";
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

      m[i] += treatsize[j1] + treatsize[j2];
    }

    // skip the replication if there is a stratum without both treatments
    if (skip) continue;

    double M = sum(m);
    NumericVector p(nstrata);

    double rmst1 = 0.0, rmst2 = 0.0, vrmst1 = 0.0, vrmst2 = 0.0;
    for (i=0; i<nstrata; i++) {
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

    rep0[index] = repn[idx[h]];
    rmst10[index] = rmst1;
    rmst20[index] = rmst2;
    vrmst10[index] = vrmst1;
    vrmst20[index] = vrmst2;
    rmstDiff0[index] = rmst1 - rmst2;
    vrmstDiff0[index] = vrmst1 + vrmst2;
    double sermstDiff = sqrt(vrmstDiff0[index]);
    rmstDiffZ0[index] = (rmstDiff0[index] - rmstDiffH0)/sermstDiff;
    rmstDiffPValue0[index] = 1.0 - R::pnorm(rmstDiffZ0[index], 0, 1, 1, 0);
    lower0[index] = rmstDiff0[index] - z*sermstDiff;
    upper0[index] = rmstDiff0[index] + z*sermstDiff;

    index++;
  }

  // only keep nonmissing records
  LogicalVector sub = !is_na(rep0);
  if (is_false(any(sub))) {
    stop("no replication enables valid inference");
  }

  rep0 = rep0[sub];
  rmst10 = rmst10[sub];
  rmst20 = rmst20[sub];
  rmstDiff0 = rmstDiff0[sub];
  vrmst10 = vrmst10[sub];
  vrmst20 = vrmst20[sub];
  vrmstDiff0 = vrmstDiff0[sub];
  rmstDiffZ0 = rmstDiffZ0[sub];
  rmstDiffPValue0 = rmstDiffPValue0[sub];
  lower0 = lower0[sub];
  upper0 = upper0[sub];

  DataFrame result = DataFrame::create(
    _["milestone"] = milestone,
    _["rmstDiffH0"] = rmstDiffH0,
    _["rmst1"] = rmst10,
    _["rmst2"] = rmst20,
    _["rmstDiff"] = rmstDiff0,
    _["vrmst1"] = vrmst10,
    _["vrmst2"] = vrmst20,
    _["vrmstDiff"] = vrmstDiff0,
    _["rmstDiffZ"] = rmstDiffZ0,
    _["rmstDiffPValue"] = rmstDiffPValue0,
    _["lower"] = lower0,
    _["upper"] = upper0,
    _["conflev"] = conflev,
    _["biascorrection"] = biascorrection);

  if (has_rep) {
    for (i=0; i<p_rep; i++) {
      String s = rep[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector repwi = u_rep[s];
        result.push_back(repwi[rep0-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector repwn = u_rep[s];
        result.push_back(repwn[rep0-1], s);
      } else if (TYPEOF(data[rep]) == STRSXP) {
        StringVector repwc = u_rep[s];
        result.push_back(repwc[rep0-1], s);
      }
    }
  }

  return result;
}


// define functions in likelihood inference for the AFT model
// algorithms adapted from the survreg function in the survival package

// log likelihood
double f_llik_1(int p, NumericVector par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int person, i, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
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

  return loglik;
}


// score vector
NumericVector f_score_1(int p, NumericVector par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int person, i, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
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
        double d1 = R::dnorm(v, 0, 1, 0), d2 = R::dnorm(u, 0, 1, 0);
        double p1 = R::pnorm(v, 0, 1, 0, 0), p2 = R::pnorm(u, 0, 1, 0, 0);
        double c1 = wt*(d1 - d2)/(p1 - p2);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(d1*v - d2*u)/(p1 - p2);
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double d1 = R::dlogis(v, 0, 1, 0), d2 = R::dlogis(u, 0, 1, 0);
        double p1 = R::plogis(v, 0, 1, 0, 0), p2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*(d1 - d2)/(p1 - p2);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(d1*v - d2*u)/(p1 - p2);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = R::dnorm(v, 0, 1, 0), d2 = R::dnorm(u, 0, 1, 0);
        double p1 = R::pnorm(v, 0, 1, 0, 0), p2 = R::pnorm(u, 0, 1, 0, 0);
        double c1 = wt*(d1 - d2)/(p1 - p2);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(d1*v - d2*u)/(p1 - p2);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = R::dlogis(v, 0, 1, 0), d2 = R::dlogis(u, 0, 1, 0);
        double p1 = R::plogis(v, 0, 1, 0, 0), p2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*(d1 - d2)/(p1 - p2);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(d1*v - d2*u)/(p1 - p2);
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

  return score;
}


// observed information matrix
NumericMatrix f_info_1(int p, NumericVector par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int person, i, j, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
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
        double w1 = exp(v), w2 = exp(u);
        double q1 = exp(-w1), q2 = exp(-w2);
        double d1 = w1*q1, d2 = w2*q2;
        double c1 = wt*(pow((d1 - d2)/(q1 - q2), 2) +
                        (d1*(1-w1) - d2*(1-w2))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double w1 = exp(v), w2 = exp(u);
        double q1 = exp(-w1), q2 = exp(-w2);
        double d1 = w1*q1, d2 = w2*q2;
        double c1 = wt*(pow((d1 - d2)/(q1 - q2), 2) +
                        (d1*(1-w1) - d2*(1-w2))/(q1 - q2));
        double c2 = wt*((d1 - d2)*(d1*v - d2*u)/pow(q1 - q2, 2) +
                        (d1*(1 + (1-w1)*v) - d2*(1 + (1-w2)*u))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((d1*v - d2*u)/(q1 - q2), 2) +
          (d1*(1 + (1-w1)*v)*v - d2*(1 + (1-w2)*u)*u)/(q1 - q2));
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double d1 = R::dnorm(v, 0, 1, 0), d2 = R::dnorm(u, 0, 1, 0);
        double q1 = R::pnorm(v, 0, 1, 0, 0), q2 = R::pnorm(u, 0, 1, 0, 0);
        double c1 = wt*(pow((d1 - d2)/(q1 -  q2), 2) +
                        (-d1*v + d2*u)/(q1 - q2));
        double c2 = wt*((d1 - d2)*(d1*v - d2*u)/pow(q1 - q2, 2) +
                        (d1*(1 - v*v) - d2*(1 - u*u))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((d1*v - d2*u)/(q1 - q2), 2) +
          (d1*(1 - v*v)*v - d2*(1 - u*u)*u)/(q1 - q2));
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double d1 = R::dlogis(v, 0, 1, 0), d2 = R::dlogis(u, 0, 1, 0);
        double q1 = R::plogis(v, 0, 1, 0, 0), q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*(pow((d1 - d2)/(q1 - q2), 2) +
                        (d1*(2*q1-1) - d2*(2*q2-1))/(q1 - q2));
        double c2 = wt*((d1 - d2)*(d1*v - d2*u)/pow(q1 - q2, 2) +
                        (d1*(1+(2*q1-1)*v) - d2*(1+(2*q2-1)*u))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((d1*v - d2*u)/(q1 - q2), 2) +
          (d1*(1+(2*q1-1)*v)*v - d2*(1+(2*q2-1)*u)*u)/(q1 - q2));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = R::dnorm(v, 0, 1, 0), d2 = R::dnorm(u, 0, 1, 0);
        double q1 = R::pnorm(v, 0, 1, 0, 0), q2 = R::pnorm(u, 0, 1, 0, 0);
        double c1 = wt*(pow((d1 - d2)/(q1 -  q2), 2) +
                        (-d1*v + d2*u)/(q1 - q2));
        double c2 = wt*((d1 - d2)*(d1*v - d2*u)/pow(q1 - q2, 2) +
                        (d1*(1 - v*v) - d2*(1 - u*u))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((d1*v - d2*u)/(q1 - q2), 2) +
          (d1*(1 - v*v)*v - d2*(1 - u*u)*u)/(q1 - q2));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = R::dlogis(v, 0, 1, 0), d2 = R::dlogis(u, 0, 1, 0);
        double q1 = R::plogis(v, 0, 1, 0, 0), q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*(pow((d1 - d2)/(q1 - q2), 2) +
                        (d1*(2*q1-1) - d2*(2*q2-1))/(q1 - q2));
        double c2 = wt*((d1 - d2)*(d1*v - d2*u)/pow(q1 - q2, 2) +
                        (d1*(1+(2*q1-1)*v) - d2*(1+(2*q2-1)*u))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat(k,j) += c2*z[j];
        }
        imat(k,k) += wt*(pow((d1*v - d2*u)/(q1 - q2), 2) +
          (d1*(1+(2*q1-1)*v)*v - d2*(1+(2*q2-1)*u)*u)/(q1 - q2));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double w2 = exp(u), q2 = exp(-w2), d2 = w2*q2;
        double c1 = wt*(pow(d2/(1 - q2), 2) - d2*(1-w2)/(1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat(i,j) += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double w2 = exp(u), q2 = exp(-w2), d2 = w2*q2;
        double c1 = wt*(pow(d2/(1 - q2), 2) - d2*(1-w2)/(1 - q2));
        double c2 = wt*(pow(d2/(1 - q2), 2)*u - d2*(1 + (1-w2)*u)/(1 - q2));
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
        double d2 = R::dnorm(u, 0, 1, 0), p2 = R::pnorm(u, 0, 1, 1, 0);
        double c1 = wt*(pow(d2/p2, 2) + d2*u/p2);
        double c2 = wt*(pow(d2/p2, 2)*u - d2*(1 - u*u)/p2);
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
        double d2 = R::dlogis(u, 0, 1, 0), q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*d2;
        double c2 = wt*(d2*u - q2);
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
        double d2 = R::dnorm(u, 0, 1, 0), p2 = R::pnorm(u, 0, 1, 1, 0);
        double c1 = wt*(pow(d2/p2, 2) + d2*u/p2);
        double c2 = wt*(pow(d2/p2, 2)*u - d2*(1 - u*u)/p2);
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
        double d2 = R::dlogis(u, 0, 1, 0), q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = wt*d2;
        double c2 = wt*(d2*u - q2);
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
        double d1 = R::dnorm(v, 0, 1, 0), q1 = R::pnorm(v, 0, 1, 0, 0);
        double c1 = wt*(pow(d1/q1, 2) - d1*v/q1);
        double c2 = wt*(pow(d1/q1, 2)*v + d1*(1 - v*v)/q1);
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
        double d1 = R::dlogis(v, 0, 1, 0), q1 = R::plogis(v, 0, 1, 0, 0);
        double c1 = wt*d1;
        double c2 = wt*(1-q1 + d1*v);
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
        double d1 = R::dnorm(v, 0, 1, 0), q1 = R::pnorm(v, 0, 1, 0, 0);
        double c1 = wt*(pow(d1/q1, 2) - d1*v/q1);
        double c2 = wt*(pow(d1/q1, 2)*v + d1*(1 - v*v)/q1);
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
        double d1 = R::dlogis(v, 0, 1, 0), q1 = R::plogis(v, 0, 1, 0, 0);
        double c1 = wt*d1;
        double c2 = wt*(1-q1 + d1*v);
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
NumericMatrix f_ressco_1(int p, NumericVector par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int nvar = param->z.ncol();
  int person, i, k;

  NumericVector eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
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
    double sigma = sig[person];
    NumericVector z = param->z(person, _)/sigma;
    k = param->strata[person] + nvar - 1;

    if (param->status[person] == 1) { // event
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -(1 - exp(u));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -(1 - exp(u));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = (1 - exp(u))*(-u) - 1;
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        for (i=0; i<nvar; i++) {
          resid(person,i) = u*z[i];
        }
        resid(person,k) = u*u - 1;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = 1 - 2*R::plogis(u, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u - 1;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        for (i=0; i<nvar; i++) {
          resid(person,i) = u*z[i];
        }
        resid(person,k) = u*u - 1;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = 1 - 2*R::plogis(u, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u - 1;
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double w1 = exp(v), w2 = exp(u);
        double q1 = exp(-w1), q2 = exp(-w2);
        double d1 = w1*q1, d2 = w2*q2;
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double w1 = exp(v), w2 = exp(u);
        double q1 = exp(-w1), q2 = exp(-w2);
        double d1 = w1*q1, d2 = w2*q2;
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = (d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double d1 = R::dnorm(v, 0, 1, 0), d2 = R::dnorm(u, 0, 1, 0);
        double q1 = R::pnorm(v, 0, 1, 0, 0), q2 = R::pnorm(u, 0, 1, 0, 0);
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = (d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double d1 = R::dlogis(v, 0, 1, 0), d2 = R::dlogis(u, 0, 1, 0);
        double q1 = R::plogis(v, 0, 1, 0, 0), q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = (d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = R::dnorm(v, 0, 1, 0), d2 = R::dnorm(u, 0, 1, 0);
        double q1 = R::pnorm(v, 0, 1, 0, 0), q2 = R::pnorm(u, 0, 1, 0, 0);
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = (d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = R::dlogis(v, 0, 1, 0), d2 = R::dlogis(u, 0, 1, 0);
        double q1 = R::plogis(v, 0, 1, 0, 0), q2 = R::plogis(u, 0, 1, 0, 0);
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = (d1*v - d2*u)/(q1 - q2);
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double w2 = exp(u), q2 = exp(-w2), d2 = w2*q2;
        double c1 = -d2/(1 - q2);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double w2 = exp(u), q2 = exp(-w2), d2 = w2*q2;
        double c1 = -d2/(1 - q2);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "lognormal") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = (-R::dnorm(u, 0, 1, 0)/R::pnorm(u, 0, 1, 1, 0));
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -R::plogis(u, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = -R::dnorm(u, 0, 1, 0)/R::pnorm(u, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = -R::plogis(u, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*u;
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = exp(v);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = exp(v);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "lognormal") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = R::dnorm(v, 0, 1, 0)/R::pnorm(v, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "loglogistic") {
        double v = (log(param->tstart[person]) - eta[person])/sigma;
        double c1 = R::plogis(v, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = R::dnorm(v, 0, 1, 0)/R::pnorm(v, 0, 1, 0, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = R::plogis(v, 0, 1, 1, 0);
        for (i=0; i<nvar; i++) {
          resid(person,i) = c1*z[i];
        }
        resid(person,k) = c1*v;
      }
    }
  }

  return resid;
}


// substitute information matrix guaranteed to be positive definite
NumericMatrix f_jj_1(int p, NumericVector par, void *ex) {
  aftparams *param = (aftparams *) ex;
  int n = param->z.nrow();
  int person, i, j;

  NumericMatrix resid = f_ressco_1(p, par, param);
  NumericMatrix jj(p,p);
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        jj(i,j) += wt*resid(person,i)*resid(person,j);
      }
    }
  }

  return jj;
}


// underlying optimization algorithm for lifereg
List liferegloop(int p, NumericVector par, void *ex,
                 int maxiter, double eps,
                 IntegerVector colfit, int ncolfit) {
  aftparams *param = (aftparams *) ex;
  int i, j, iter, halving = 0;
  bool fail;

  int nstrata = param->nstrata;
  int nvar = param->z.ncol();
  int nsub = param->z.nrow();
  NumericMatrix z1 = param->z;

  // standardize the design matrix
  NumericVector mu(nvar), sigma(nvar);
  NumericMatrix z2(nsub, nvar);
  mu[0] = 0;
  sigma[0] = 1;
  z2(_,0) = z1(_,0);
  for (i=1; i<nvar; i++) {
    NumericVector u = z1(_,i);
    mu[i] = mean(u);
    sigma[i] = sd(u);
    // no standardization for constant columns
    if (sigma[i] == 0) sigma[i] = 1;
    z2(_,i) = (u - mu[i])/sigma[i];
  }

  aftparams para = {param->dist, param->strata, param->tstart, param->tstop,
                    param->status, param->weight, param->offset, z2, nstrata};

  double toler = 1e-12;
  NumericVector beta(p), newbeta(p);
  double loglik, newlk;
  NumericVector u(p);
  NumericMatrix imat(p,p);
  NumericMatrix jj(p,p);
  NumericVector u1(ncolfit);
  NumericMatrix imat1(ncolfit, ncolfit);
  NumericMatrix jj1(ncolfit, ncolfit);

  // initial beta and log likelihood
  beta[0] = par[0];
  for (i=1; i<nvar; i++) {
    beta[i] = par[i]*sigma[i];
    beta[0] += par[i]*mu[i];
  }

  loglik = f_llik_1(p, beta, &para);
  u = f_score_1(p, beta, &para);
  for (i=0; i<ncolfit; i++) {
    u1[i] = u[colfit[i]];
  }

  imat = f_info_1(p, beta, &para);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      imat1(i,j) = imat(colfit[i], colfit[j]);
    }
  }

  i = cholesky2(imat1, ncolfit, toler);
  if (i < 0) {
    jj = f_jj_1(p, beta, &para);
    for (i=0; i<ncolfit; i++) {
      for (j=0; j<ncolfit; j++) {
        jj1(i,j) = jj(colfit[i], colfit[j]);
      }
    }

    i = cholesky2(jj1, ncolfit, toler);
    chsolve2(jj1, ncolfit, u1);
  } else {
    chsolve2(imat1, ncolfit, u1);
  }

  u.fill(0);
  for (i=0; i<ncolfit; i++) {
    u[colfit[i]] = u1[i];
  }

  // new beta
  for (i=0; i<p; i++) {
    newbeta[i] = beta[i] + u[i];
  }

  for (iter=0; iter<maxiter; iter++) {
    // new log likelihood
    newlk = f_llik_1(p, newbeta, &para);

    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk) == 1;

    if (!fail && halving == 0 && fabs(1 - (loglik/newlk)) < eps) {
      break;
    }

    if (fail || newlk < loglik) { // adjust step size if likelihood decreases
      halving++;
      for (i=0; i<p; i++) {
        newbeta[i] = (beta[i] + newbeta[i]) /2;
      }

      // special handling of sigmas
      if (halving == 1 && param->dist != "exponential") {
        for (i=0; i<nstrata; i++) {
          if (beta[nvar+i] - newbeta[nvar+i] > 1.1) {
            newbeta[nvar+i] = beta[nvar+i] - 1.1;
          }
        }
      }
    } else { // update beta normally
      halving = 0;

      for (i=0; i<p; i++) {
        beta[i] = newbeta[i];
      }
      loglik = newlk;

      u = f_score_1(p, beta, &para);
      for (i=0; i<ncolfit; i++) {
        u1[i] = u[colfit[i]];
      }

      imat = f_info_1(p, beta, &para);
      for (i=0; i<ncolfit; i++) {
        for (j=0; j<ncolfit; j++) {
          imat1(i,j) = imat(colfit[i], colfit[j]);
        }
      }

      i = cholesky2(imat1, ncolfit, toler);
      if (i < 0) {
        jj = f_jj_1(p, beta, &para);
        for (i=0; i<ncolfit; i++) {
          for (j=0; j<ncolfit; j++) {
            jj1(i,j) = jj(colfit[i], colfit[j]);
          }
        }

        i = cholesky2(jj1, ncolfit, toler);
        chsolve2(jj1, ncolfit, u1);
      } else {
        chsolve2(imat1, ncolfit, u1);
      }

      u.fill(0);
      for (i=0; i<ncolfit; i++) {
        u[colfit[i]] = u1[i];
      }

      for (i=0; i<p; i++) {
        newbeta[i] = beta[i] + u[i];
      }
    }
  }

  if (iter == maxiter) fail = 1;

  // parameter estimates on the original scale of the design matrix
  for (i=1; i<nvar; i++) {
    newbeta[i] = newbeta[i]/sigma[i];
    newbeta[0] = newbeta[0] - newbeta[i]*mu[i];
  }

  imat = f_info_1(p, newbeta, param);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      imat1(i,j) = imat(colfit[i], colfit[j]);
    }
  }

  NumericMatrix var1 = invsympd(imat1, ncolfit, toler);
  NumericMatrix var(p,p);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      var(colfit[i], colfit[j]) = var1(i,j);
    }
  }

  return List::create(
    Named("coef") = newbeta,
    Named("iter") = iter,
    Named("var") = var,
    Named("loglik") = newlk,
    Named("fail") = fail);
}


// confidence limit of profile likelihood method
double liferegplloop(int p, NumericVector par, void *ex,
                     int maxiter, double eps,
                     int k, int which, double l0) {
  aftparams *param = (aftparams *) ex;

  int i, j, iter;
  bool fail = 0;

  NumericMatrix z1 = param->z;

  double toler = 1e-12;
  NumericVector beta(p), newbeta(p);
  double loglik, newlk;
  NumericVector u(p);
  NumericVector delta(p);
  NumericMatrix imat(p,p);
  NumericMatrix jj(p,p);
  NumericMatrix v(p,p);

  // initial beta and log likelihood
  for (i=0; i<p; i++) {
    beta[i] = par[i];
  }

  loglik = f_llik_1(p, beta, param);
  u = f_score_1(p, beta, param);

  imat = f_info_1(p, beta, param);
  jj = clone(imat);

  i = cholesky2(jj, p, toler);
  if (i < 0) {
    jj = f_jj_1(p, beta, param);
    v = invsympd(jj, p, toler);
  } else {
    v = invsympd(imat, p, toler);
  }
  v = -1.0*v;

  // Lagrange multiplier method as used in SAS PROC LOGISTIC
  double w = 0;
  for (i=0; i<p; i++) {
    for (j=0; j<p; j++) {
      w += u[i]*v(i,j)*u[j];
    }
  }

  double underroot = 2*(l0 - loglik + 0.5*w)/v(k,k);
  double lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
  u[k] += lambda;

  delta.fill(0.0);
  for (i=0; i<p; i++) {
    for (j=0; j<p; j++) {
      delta[i] -= v(i,j)*u[j];
    }
  }

  // update beta
  for (i=0; i<p; i++) {
    newbeta[i] = beta[i] + delta[i];
  }

  for (iter=0; iter<maxiter; iter++) {
    newlk = f_llik_1(p, newbeta, param);

    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk) == 1;

    if (!fail && fabs((newlk - l0)/l0) < eps) {
      break;
    }

    for (i=0; i<p; i++) {
      beta[i] = newbeta[i];
    }
    loglik = newlk;

    u = f_score_1(p, beta, param);

    imat = f_info_1(p, beta, param);
    jj = clone(imat);

    i = cholesky2(jj, p, toler);
    if (i < 0) {
      jj = f_jj_1(p, beta, param);
      v = invsympd(jj, p, toler);
    } else {
      v = invsympd(imat, p, toler);
    }
    v = -1.0*v;

    // Lagrange multiplier method as used in SAS PROC LOGISTIC
    w = 0;
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        w += u[i]*v(i,j)*u[j];
      }
    }

    underroot = 2*(l0 - loglik + 0.5*w)/v(k,k);
    lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
    u[k] += lambda;

    delta.fill(0.0);
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        delta[i] -= v(i,j)*u[j];
      }
    }
    // update beta
    for (i=0; i<p; i++) {
      newbeta[i] = beta[i] + delta[i];
    }
  }

  if (iter == maxiter) fail = 1;

  if (fail) {
    stop("The algorithm in liferegplloop did not converge");
  }

  return newbeta[k];
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
                const bool robust = 0,
                const bool plci = 0,
                const double alpha = 0.05) {

  int h, i, j, k, n = data.nrows();
  int nvar = static_cast<int>(covariates.size()) + 1;
  if (nvar == 2 && (covariates[0] == "" || covariates[0] == "none")) {
    nvar = 1;
  }

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


  bool has_rep;
  IntegerVector repn(n);
  DataFrame u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    repn.fill(1);
  } else {
    List out = bygroup(data, rep);
    has_rep = 1;
    repn = out["index"];
    u_rep = DataFrame(out["lookup"]);
  }


  IntegerVector stratumn(n);
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    stratumn.fill(1);
  } else {
    List out = bygroup(data, stratum);
    stratumn = out["index"];
  }


  IntegerVector stratumn1 = unique(stratumn);
  int nstrata = static_cast<int>(stratumn1.size());
  int p = dist1 == "exponential" ? nvar : (nvar+nstrata);

  if (dist1 == "exponential" && nstrata > 1) {
    stop("Stratification is not valid with the exponential distribution");
  }


  bool has_time = hasVariable(data, time);
  if (!has_time) stop("data must contain the time variable");

  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  for (i=0; i<n; i++) {
    if (!std::isnan(timen[i]) && ((dist1 == "exponential") ||
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
    NumericVector time2nz = data[time2];
    time2n = clone(time2nz);
    for (i=0; i<n; i++) {
      if (!std::isnan(time2n[i]) && ((dist1 == "exponential") ||
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
    IntegerVector eventnz = data[event];
    eventn = clone(eventnz);
    if (is_true(any((eventn != 1) & (eventn != 0)))) {
      stop("event must be 1 or 0 for each subject");
    }

    if (is_true(all(eventn == 0))) {
      stop("at least 1 event is needed to fit the parametric model");
    }
  }

  NumericMatrix zn(n,nvar);
  for (i=0; i<n; i++) {
    zn(i,0) = 1; // intercept
  }

  for (j=0; j<nvar-1; j++) {
    String zj = covariates[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in covariates");
    }
    NumericVector u = data[zj];
    for (i=0; i<n; i++) {
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

  bool has_id = hasVariable(data, id);

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


  // sort the data by rep
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

  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (i=0; i<n; i++) {
    if ((repn[i] == NA_INTEGER) || (stratumn[i] == NA_INTEGER) ||
        (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
        (std::isnan(weightn[i])) || (std::isnan(offsetn[i])) ||
        (idn[i] == NA_INTEGER)) {
      sub[i] = 0;
    }
    for (j=0; j<nvar-1; j++) {
      if (std::isnan(zn(i,j+1))) sub[i] = 0;
    }
  }

  order = which(sub);
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
  NumericVector loglik0(nreps), loglik1(nreps);
  IntegerVector niter(nreps);

  IntegerVector rep0(nreps*p);
  StringVector par0(nreps*p);
  NumericVector beta0(nreps*p), sebeta0(nreps*p), rsebeta0(nreps*p);
  NumericMatrix vbeta0(nreps*p,p), rvbeta0(nreps*p,p);
  NumericVector lb0(nreps*p), ub0(nreps*p), prob0(nreps*p);
  StringVector clparm0(nreps*p);

  for (h=0; h<nreps; h++) {
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
      for (i=0; i<n1; i++) {
        tstop[i] = event1[i] == 1 ? tstart[i] : NA_REAL;
      }
    } else {
      tstart = time1;
      tstop = time21;
    }

    IntegerVector status(n1);
    for (i=0; i<n1; i++) {
      if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
          (tstart[i] == tstop[i])) {
        status[i] = 1; // event
      } else if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
        (tstart[i] < tstop[i])) {
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
    for (i=0; i<n2; i++) {
      if (status[i] == 0 || status[i] == 1) { // right censoring or event
        time0[i] = tstart[i];
      } else if (status[i] == 2) { // left censoring
        time0[i] = tstop[i];
      } else if (status[i] == 3) { // interval censoring
        time0[i] = (tstart[i] + tstop[i])/2;
      }
    }

    NumericVector y0 = clone(time0);
    if ((dist1 == "exponential") || (dist1 == "weibull") ||
        (dist1 == "lognormal") || (dist1 == "loglogistic")) {
      y0 = log(y0);
    }

    double int0 = mean(y0);
    double logsig0 = log(sd(y0));

    NumericVector bint0(p);
    int ncolfit0 = dist1 == "exponential" ? 1 : nstrata + 1;
    IntegerVector colfit0(ncolfit0);
    if (dist1 == "exponential") {
      bint0[0] = int0;
      ncolfit0 = 1;
      colfit0[0] = 0;
    } else {
      bint0[0] = int0;
      for (i=0; i<nstrata; i++) {
        bint0[nvar+i] = logsig0;
      }

      colfit0[0] = 0;
      for (i=0; i<nstrata; i++) {
        colfit0[i+1] = nvar+i;
      }
    }

    // parameter estimates and standard errors for the null model
    aftparams param = {dist1, stratum1, tstart, tstop, status, weight1,
                          offset1, z1, nstrata};

    List outint = liferegloop(p, bint0, &param, 100, 1.0e-9,
                              colfit0, ncolfit0);
    NumericVector bint = outint["coef"];
    NumericMatrix vbint = as<NumericMatrix>(outint["var"]);

    NumericVector b(p);
    NumericMatrix vb(p,p);
    List out;

    if (nvar > 1) {
      // parameter estimates and standard errors for the full model
      IntegerVector colfit = seq(0,p-1);
      out = liferegloop(p, bint, &param, 100, 1.0e-9, colfit, p);

      bool fail = out["fail"];
      if (fail) stop("The algorithm in liferegr did not converge");

      b = out["coef"];
      vb = as<NumericMatrix>(out["var"]);
    } else {
      b = bint;
      vb = vbint;
      out = outint;
    }

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

    niter[h] = out["iter"];

    // robust variance estimates
    NumericVector rseb(p);  // robust standard error for betahat
    if (robust) {
      NumericMatrix ressco = f_ressco_1(p, b, &param);

      int nr; // number of rows in the score residual matrix
      if (!has_id) {
        for (i=0; i<n2; i++) {
          for (j=0; j<p; j++) {
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
        for (i=1; i<n2; i++) {
          if (id2[i] != id2[i-1]) {
            idx.push_back(i);
          }
        }

        int nids = static_cast<int>(idx.size());
        idx.push_back(n2);

        NumericMatrix resid(n2,p);
        for (i=0; i<n2; i++) {
          for (j=0; j<p; j++) {
            resid(i,j) += ressco(order[i],j);
          }
        }

        NumericVector weight2 = weight1[order];

        NumericMatrix ressco2(nids,p);
        for (i=0; i<nids; i++) {
          for (j=0; j<p; j++) {
            for (k=idx[i]; k<idx[i+1]; k++) {
              ressco2(i,j) += weight2[k]*resid(k,j);
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
            D(i,j) += ressco(i,k)*vb(k,j);
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

      for (i=0; i<p; i++) {
        rseb[i] = sqrt(rvb(i,i));
      }

      for (i=0; i<p; i++) {
        rsebeta0[h*p+i] = rseb[i];
        for (j=0; j<p; j++) {
          rvbeta0(h*p+i,j) = rvb(i,j);
        }
      }
    }

    // profile likelihood confidence interval for regression coefficients
    NumericVector lb(p), ub(p), prob(p);
    StringVector clparm(p);

    double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
    if (plci) {
      double lmax = f_llik_1(p, b, &param);
      double l0 = lmax - 0.5*R::qchisq(1-alpha, 1, 1, 0);

      for (k=0; k<p; k++) {
        lb[k] = liferegplloop(p, b, &param, 30, 1.0e-8, k, -1, l0);
        ub[k] = liferegplloop(p, b, &param, 30, 1.0e-8, k, 1, l0);

        IntegerVector colfit1(p-1);
        for (i=0; i<k; i++) {
          colfit1[i] = i;
        }
        for (i=k+1; i<p; i++) {
          colfit1[i-1] = i;
        }

        NumericVector b0(p);
        List out0 = liferegloop(p, b0, &param, 30, 1.0e-8, colfit1, p-1);
        double lmax0 = out0["loglik"];
        prob[k] = R::pchisq(-2*(lmax0 - lmax), 1, 0, 0);
        clparm[k] = "PL";
      }
    } else {
      for (k=0; k<p; k++) {
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

    for (i=0; i<p; i++) {
      lb0[h*p+i] = lb[i];
      ub0[h*p+i] = ub[i];
      prob0[h*p+i] = prob[i];
      clparm0[h*p+i] = clparm[i];
    }

    // log-likelihoods
    loglik0[h] = outint["loglik"];
    loglik1[h] = out["loglik"];
  }

  NumericVector expbeta0 = exp(beta0);
  NumericVector z0(nreps*p);
  if (!robust) z0 = beta0/sebeta0;
  else z0 = beta0/rsebeta0;

  DataFrame sumstat = List::create(
    _["n"] = nobs,
    _["nevents"] = nevents,
    _["loglik0"] = loglik0,
    _["loglik1"] = loglik1,
    _["niter"] = niter,
    _["dist"] = dist1,
    _["p"] = p,
    _["nvar"] = nvar-1,
    _["robust"] = robust);

  DataFrame parest;
  if (!robust) {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = sebeta0,
      _["z"] = z0,
      _["expbeta"] = expbeta0,
      _["vbeta"] = vbeta0,
      _["lower"] = lb0,
      _["upper"] = ub0,
      _["p"] = prob0,
      _["method"] = clparm0);
  } else {
    parest = DataFrame::create(
      _["param"] = par0,
      _["beta"] = beta0,
      _["sebeta"] = rsebeta0,
      _["z"] = z0,
      _["expbeta"] = expbeta0,
      _["vbeta"] = rvbeta0,
      _["lower"] = lb0,
      _["upper"] = ub0,
      _["p"] = prob0,
      _["method"] = clparm0,
      _["sebeta_naive"] = sebeta0,
      _["vbeta_naive"] = vbeta0);
  }

  if (has_rep) {
    for (i=0; i<p_rep; i++) {
      String s = rep[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector repwi = u_rep[s];
        sumstat.push_back(repwi[rep01-1], s);
        parest.push_back(repwi[rep0-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector repwn = u_rep[s];
        sumstat.push_back(repwn[rep01-1], s);
        parest.push_back(repwn[rep0-1], s);
      } else if (TYPEOF(data[rep]) == STRSXP) {
        StringVector repwc = u_rep[s];
        sumstat.push_back(repwc[rep01-1], s);
        parest.push_back(repwc[rep0-1], s);
      }
    }
  }

  List result = List::create(
    _["sumstat"] = sumstat,
    _["parest"] = parest);

  return result;
}


// define functions in likelihood inference for the Cox model
// algorithms adapted from the coxph function in the survival package

// log likelihood
double f_llik_2(int p, NumericVector par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, k, person;

  double loglik = 0;        // log-likelihood value
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
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

  return loglik;
}


// score vector
NumericVector f_score_2(int p, NumericVector par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, k, person;

  NumericVector u(p);       // score vector
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;              // zbar(beta,k,t)

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      denom = 0;
      for (i=0; i<p; i++) {
        a[i] = 0;
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
        }
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
      for (i=0; i<p; i++) {
        a[i] -= risk*param->z(p1,i);
      }
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
      for (i=0; i<p; i++) {
        a2[i] = 0;
      }
    }
  }

  return u;
}


// observed information matrix
NumericMatrix f_info_2(int p, NumericVector par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, j, k, person;

  NumericMatrix imat(p,p);  // information matrix
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
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
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
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
        for (j=0; j<=i; j++) {
          cmat(i,j) = 0;
        }
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
        for (j=0; j<=i; j++) {
          cmat2(i,j) = 0;
        }
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


// Firth's penalized log likelihood
double f_pen_llik_2(int p, NumericVector par, void *ex) {
  coxparams *param = (coxparams *) ex;
  double loglik = f_llik_2(p, par, param);
  NumericMatrix imat = f_info_2(p, par, param);

  // obtain the determinant of imat
  double toler = 1e-12;
  int i = cholesky2(imat, p, toler);

  double v = 0;
  for (i=0; i<p; i++) {
    v += log(imat(i,i));
  }

  return loglik + 0.5*v;
}


// score vector
NumericVector f_pen_score_2(int p, NumericVector par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int h, i, j, k, l, person;

  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;              // zbar(beta,k,t)
  NumericMatrix cmat(p,p);  // s2(beta,k,t)
  NumericMatrix cmat2(p,p); // sum of w*exp(zbeta)*z*z' for the deaths
  NumericMatrix dmat(p,p*p);  // q2(beta,k,t)
  NumericMatrix dmat2(p,p*p); // sum of w*exp(zbeta)*z*z*z' for the deaths

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector u(p);       // score vector
  NumericMatrix imat(p,p);  // information matrix
  NumericMatrix dimat(p,p*p);

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  double t1, t2, t3;
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      denom = 0;
      for (i=0; i<p; i++) {
        a[i] = 0;
        for (j=0; j<=i; j++) {
          cmat(i,j) = 0;
          for (k=0; k<=j; k++) {
            dmat(i,k*p+j) = 0;
          }
        }
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
            for (k=0; k<=j; k++) {
              dmat(i,k*p+j) += risk*param->z(person,i)*param->z(person,j)*
                param->z(person,k);
            }
          }
        }
      } else {
        ndead++;
        deadwt += param->weight[person];
        denom2 += risk;
        for (i=0; i<p; i++) {
          a2[i] += risk*param->z(person,i);
          u[i] += param->weight[person]*param->z(person,i);
          for (j=0; j<=i; j++) {
            cmat2(i,j) += risk*param->z(person,i)*param->z(person,j);
            for (k=0; k<=j; k++) {
              dmat2(i,k*p+j) += risk*param->z(person,i)*param->z(person,j)*
                param->z(person,k);
            }
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
          for (k=0; k<=j; k++) {
            dmat(i,k*p+j) -= risk*param->z(p1,i)*param->z(p1,j)*
              param->z(p1,k);
          }
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
          u[i] -= deadwt*xbar;
          for (j=0; j<=i; j++) {
            cmat(i,j) += cmat2(i,j);
            imat(i,j) += deadwt*(cmat(i,j) - xbar*a[j])/denom;
            for (k=0; k<=j; k++) {
              dmat(i,k*p+j) += dmat2(i,k*p+j);
              t1 = dmat(i,k*p+j) - cmat(i,j)*a[k]/denom;
              t2 = (cmat(i,k) - a[i]*a[k]/denom)*a[j]/denom;
              t3 = a[i]/denom*(cmat(j,k) - a[j]*a[k]/denom);
              dimat(i,k*p+j) += deadwt*(t1 - t2 - t3)/denom;
            }
          }
        }
      } else {
        meanwt = deadwt/ndead;
        for (l=0; l<ndead; l++) {
          denom += denom2/ndead;
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            xbar = a[i]/denom;
            u[i] -= meanwt*xbar;
            for (j=0; j<=i; j++) {
              cmat(i,j) += cmat2(i,j)/ndead;
              imat(i,j) += meanwt*(cmat(i,j) - xbar*a[j])/denom;
              for (k=0; k<=j; k++) {
                dmat(i,k*p+j) += dmat2(i,k*p+j)/ndead;
                t1 = dmat(i,k*p+j) - cmat(i,j)*a[k]/denom;
                t2 = (cmat(i,k) - a[i]*a[k]/denom)*a[j]/denom;
                t3 = a[i]/denom*(cmat(j,k) - a[j]*a[k]/denom);
                dimat(i,k*p+j) += meanwt*(t1 - t2 - t3)/denom;
              }
            }
          }
        }
      }

      // reset for the next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
      for (i=0; i<p; i++) {
        a2[i] = 0;
        for (j=0; j<=i; j++) {
          cmat2(i,j) = 0;
          for (k=0; k<=j; k++) {
            dmat2(i,k*p+j) = 0;
          }
        }
      }
    }
  }

  // fill the symmetric elements of the imat matrix
  for (i=0; i<p-1; i++) {
    for (j=i+1; j<p; j++) {
      imat(i,j) = imat(j,i);
    }
  }

  // fill the symmetric elements of the dimat array
  for (i=0; i<p-1; i++) {
    for (j=i+1; j<p; j++) {
      for (k=0; k<=i; k++) {
        dimat(i,k*p+j) = dimat(j,k*p+i);
      }
    }
  }

  for (j=0; j<p-1; j++) {
    for (i=j; i<p; i++) {
      for (k=j+1; k<p; k++) {
        dimat(i,k*p+j) = dimat(i,j*p+k);
      }
    }
  }

  for (i=0; i<p-1; i++) {
    for (j=i+1; j<p; j++) {
      for (k=i+1; k<p; k++) {
        dimat(i,k*p+j) = dimat(k,i*p+j);
      }
    }
  }

  // obtain the penalized score vector
  double toler = 1e-12;
  i = cholesky2(imat, p, toler);

  double temp;
  NumericMatrix y(p,p);
  NumericVector g(p);
  for (k=0; k<p; k++) {

    // partial derivative of the information matrix w.r.t. beta[k]
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        y(i,j) = dimat(i,k*p+j);
      }
    }

    // solve(imat, y)
    for (h=0; h<p; h++) {
      for (i=0; i<p; i++) {
        temp = y(i,h);
        for (j=0; j<i; j++)
          temp -= y(j,h)*imat(j,i);
        y(i,h) = temp;
      }

      for (i=p-1; i>=0; i--) {
        if (imat(i,i) == 0) y(i,h) = 0;
        else {
          temp = y(i,h)/imat(i,i);
          for (j=i+1; j<p; j++)
            temp -= y(j,h)*imat(i,j);
          y(i,h) = temp;
        }
      }
    }

    // trace
    for (i=0; i<p; i++) {
      g[k] += y(i,i);
    }

    g[k] = u[k] + 0.5*g[k];
  }

  return g;
}


// score residual matrix
NumericMatrix f_ressco_2(int p, NumericVector par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, j, k, person, n = static_cast<int>(param->tstart.size());

  NumericMatrix resid(n,p);
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;              // zbar(beta,k,t)

  double downwt, hazard, cumhaz = 0;
  NumericVector xhaz(p), mh1(p), mh2(p), mh3(p);

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  NumericVector score = exp(eta);

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
        a[i] = 0;
        xhaz[i] = 0;
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
        for (i=0; i<p; i++) {
          a[i] += risk*param->z(person,i);
        }
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
      for (i=0; i<p; i++) {
        a2[i] = 0;
      }
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


// substitute information matrix guaranteed to be positive definite
NumericMatrix f_jj_2(int p, NumericVector par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int n = param->z.nrow();
  int person, i, j;

  NumericMatrix resid = f_ressco_2(p, par, param);
  NumericMatrix jj(p,p);
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        jj(i,j) += wt*resid(person,i)*resid(person,j);
      }
    }
  }

  return jj;
}


// underlying optimization algorithm for phreg
List phregloop(int p, NumericVector par, void *ex,
               int maxiter, double eps, bool firth,
               IntegerVector colfit, int ncolfit) {
  coxparams *param = (coxparams *) ex;
  int i, j, iter, halving = 0;
  bool fail;

  NumericMatrix z1 = param->z;

  double toler = 1e-12;
  NumericVector beta(p), newbeta(p);
  double loglik, newlk = 0;
  NumericVector u(p);
  NumericMatrix imat(p,p);
  NumericVector u1(ncolfit);
  NumericMatrix imat1(ncolfit, ncolfit);

  // initial beta and log likelihood
  for (i=0; i<p; i++) {
    beta[i] = par[i];
  }

  if (!firth) {
    loglik = f_llik_2(p, beta, param);
    u = f_score_2(p, beta, param);
  } else{
    loglik = f_pen_llik_2(p, beta, param);
    u = f_pen_score_2(p, beta, param);
  }
  for (i=0; i<ncolfit; i++) {
    u1[i] = u[colfit[i]];
  }

  imat = f_info_2(p, beta, param);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      imat1(i,j) = imat(colfit[i], colfit[j]);
    }
  }

  i = cholesky2(imat1, ncolfit, toler);

  chsolve2(imat1, ncolfit, u1);

  u.fill(0);
  for (i=0; i<ncolfit; i++) {
    u[colfit[i]] = u1[i];
  }

  for (i=0; i<p; i++) {
    newbeta[i] = beta[i] + u[i];
  }

  for (iter=0; iter<maxiter; iter++) {
    // new log likelihood
    if (!firth) newlk = f_llik_2(p, newbeta, param);
    else newlk = f_pen_llik_2(p, newbeta, param);

    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk) == 1;

    if (!fail && halving == 0 && fabs(1 - (loglik/newlk)) < eps) {
      break;
    }

    if (fail || newlk < loglik) { // adjust step size if likelihood decreases
      halving++;
      for (i=0; i<p; i++) {
        newbeta[i] = (beta[i] + newbeta[i])/2;
      }
    } else { // update beta normally
      halving = 0;

      for (i=0; i<p; i++) {
        beta[i] = newbeta[i];
      }
      loglik = newlk;

      if (!firth) u = f_score_2(p, beta, param);
      else u = f_pen_score_2(p, beta, param);
      for (i=0; i<ncolfit; i++) {
        u1[i] = u[colfit[i]];
      }

      imat = f_info_2(p, beta, param);
      for (i=0; i<ncolfit; i++) {
        for (j=0; j<ncolfit; j++) {
          imat1(i,j) = imat(colfit[i], colfit[j]);
        }
      }

      i = cholesky2(imat1, ncolfit, toler);

      chsolve2(imat1, ncolfit, u1);
      u.fill(0);
      for (i=0; i<ncolfit; i++) {
        u[colfit[i]] = u1[i];
      }

      for (i=0; i<p; i++) {
        newbeta[i] = beta[i] + u[i];
      }
    }
  }

  if (iter == maxiter) fail = 1;

  imat = f_info_2(p, newbeta, param);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      imat1(i,j) = imat(colfit[i], colfit[j]);
    }
  }

  NumericMatrix var1 = invsympd(imat1, ncolfit, toler);
  NumericMatrix var(p,p);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      var(colfit[i], colfit[j]) = var1(i,j);
    }
  }

  return List::create(
    Named("coef") = newbeta,
    Named("iter") = iter,
    Named("var") = var,
    Named("loglik") = newlk,
    Named("fail") = fail);
}


// confidence limit of profile likelihood method
double phregplloop(int p, NumericVector par, void *ex,
                   int maxiter, double eps, bool firth,
                   int k, int which, double l0) {
  coxparams *param = (coxparams *) ex;

  int i, j, iter;
  bool fail = 0;

  NumericMatrix z1 = param->z;

  NumericVector beta(p), newbeta(p);
  double loglik, newlk;
  NumericVector u(p);
  NumericVector delta(p);
  NumericMatrix imat(p,p);
  NumericMatrix v(p,p);

  // initial beta and log likelihood
  for (i=0; i<p; i++) {
    beta[i] = par[i];
  }

  if (!firth) {
    loglik = f_llik_2(p, beta, param);
    u = f_score_2(p, beta, param);
  } else {
    loglik = f_pen_llik_2(p, beta, param);
    u = f_pen_score_2(p, beta, param);
  }

  imat = f_info_2(p, beta, param);

  // Lagrange multiplier method as used in SAS PROC LOGISTIC
  double toler = 1e-12;
  v = invsympd(imat, p, toler);
  v = -1.0*v;

  double w = 0;
  for (i=0; i<p; i++) {
    for (j=0; j<p; j++) {
      w += u[i]*v(i,j)*u[j];
    }
  }
  double underroot = 2*(l0 - loglik + 0.5*w)/v(k,k);
  double lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
  u[k] += lambda;

  delta.fill(0.0);
  for (i=0; i<p; i++) {
    for (j=0; j<p; j++) {
      delta[i] -= v(i,j)*u[j];
    }
  }

  // update beta
  for (i=0; i<p; i++) {
    newbeta[i] = beta[i] + delta[i];
  }

  for (iter=0; iter<maxiter; iter++) {
    if (!firth) newlk = f_llik_2(p, newbeta, param);
    else newlk = f_pen_llik_2(p, newbeta, param);

    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk) == 1;

    if (!fail && fabs((newlk - l0)/l0) < eps) {
      break;
    }

    for (i=0; i<p; i++) {
      beta[i] = newbeta[i];
    }
    loglik = newlk;

    if (!firth) u = f_score_2(p, beta, param);
    else u = f_pen_score_2(p, beta, param);

    imat = f_info_2(p, beta, param);

    // Lagrange multiplier method as used in SAS PROC LOGISTIC
    v = invsympd(imat, p, toler);
    v = -1.0*v;

    w = 0;
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        w += u[i]*v(i,j)*u[j];
      }
    }
    underroot = 2*(l0 - loglik + 0.5*w)/v(k,k);
    lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
    u[k] += lambda;

    delta.fill(0.0);
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        delta[i] -= v(i,j)*u[j];
      }
    }
    // update beta
    for (i=0; i<p; i++) {
      newbeta[i] = beta[i] + delta[i];
    }
  }

  if (iter == maxiter) fail = 1;

  if (fail) {
    stop("The algorithm in phregplloop did not converge");
  }

  return newbeta[k];
}


// baseline hazard estimates
List f_basehaz(int p, NumericVector par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, k, person;

  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths
  int natrisk = 0;

  // unique times
  IntegerVector idx = Range(0, param->nused-1);
  NumericVector time0 = param->tstop[idx];
  NumericVector time1 = unique(time0);
  int J = static_cast<int>(time1.size());

  IntegerVector stratum(J);
  NumericVector time(J), nrisk(J), nevent(J), haz(J), varhaz(J);
  NumericMatrix gradhaz(J,p);

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  int j = J; // index the unique time in ascending order
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      natrisk = 0;
      denom = 0;
      for (i=0; i<p; i++) {
        a[i] = 0;
      }
    }

    dtime = param->tstop[person];
    bool first = 1;
    while (person < param->nused && param->tstop[person] == dtime) {
      if (first) { // first incidence at this time
        j--;
        stratum[j] = param->strata[person];
        time[j] = dtime;
        first = 0;
      }

      // walk through this set of tied times
      risk = param->weight[person]*exp(eta[person]);
      natrisk++;
      if (param->event[person] == 0) {
        denom += risk;
        for (i=0; i<p; i++) {
          a[i] += risk*param->z(person,i);
        }
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
      risk = param->weight[p1]*exp(eta[p1]);
      natrisk--;
      denom -= risk;
      for (i=0; i<p; i++) {
        a[i] -= risk*param->z(p1,i);
      }
    }

    // add to the main terms
    nrisk[j] = natrisk;
    nevent[j] = ndead;
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        haz[j] = deadwt/denom;
        varhaz[j] = deadwt/denom/denom;
        for (i=0; i<p; i++) {
          a[i] += a2[i];
          gradhaz(j,i) = deadwt*a[i]/denom/denom;
        }
      } else {
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          haz[j] += meanwt/denom;
          varhaz[j] += meanwt/denom/denom;
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            gradhaz(j,i) += meanwt*a[i]/denom/denom;
          }
        }
      }

      // reset for the next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
      for (i=0; i<p; i++) {
        a2[i] = 0;
      }
    }
  }

  List result = List::create(
    _["stratum"] = stratum,
    _["time"] = time,
    _["nrisk"] = nrisk,
    _["nevent"] = nevent,
    _["haz"] = haz,
    _["varhaz"] = varhaz
  );

  if (p > 0) result.push_back(gradhaz, "gradhaz");

  return result;
}


// martingale residuals
NumericVector f_resmart(int p, NumericVector par, void *ex) {
  int i, person;
  coxparams *param = (coxparams *) ex;

  double dtime;             // distinct time
  double ndead = 0;         // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  int n = static_cast<int>(param->tstop.size());
  NumericVector resid(n);
  for (person = 0; person < param->nused; person++) {
    resid[person] = param->event[person];
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  int j, j0 = 0, j1, j2;
  double hazard;

  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      j0 = person;  // first person in the stratum
      denom = 0;
    }

    dtime = param->tstop[person];

    j1 = person;   // first person in the stratum with the tied time
    while (person < param->nused && param->tstop[person] == dtime) {
      // walk through this set of tied times
      risk = param->weight[person]*exp(eta[person]);
      if (param->event[person] == 0) {
        denom += risk;
      } else {
        ndead++;
        deadwt += param->weight[person];
        denom2 += risk;
      }

      person++;

      if (person < param->nused && param->strata[person] != istrata) break;
    }

    j2 = person-1; // last person in the stratum with the tied time

    // remove subjects no longer at risk
    for (; i1<param->nused; i1++) {
      p1 = param->order1[i1];
      if (param->tstart[p1] < dtime || param->strata[p1] != istrata) break;
      risk = param->weight[p1]*exp(eta[p1]);
      denom -= risk;
    }

    // add to the main terms
    if (ndead > 0) {
      denom += denom2;

      for (j=j0; j<=j2; j++) {
        if (param->tstart[j] < dtime) {
          if (param->method == 0 || ndead == 1) {
            hazard = deadwt/denom;
          } else {
            meanwt = deadwt/ndead;
            hazard = 0;

            if (j<j1 || param->event[j]==0) {
              for (i=0; i<ndead; i++) {
                hazard += meanwt/(denom - i/ndead*denom2);
              }
            } else {
              for (i=0; i<ndead; i++) {
                hazard += (1 - i/ndead)*meanwt/(denom - i/ndead*denom2);
              }
            }
          }
          resid[j] -= hazard*exp(eta[j]);
        }
      }

      // reset for the next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
    }
  }

  return resid;
}


// schoenfeld residuals
List f_ressch(int p, NumericVector par, void *ex) {
  coxparams *param = (coxparams *) ex;
  int i, k, person;

  int nevent = sum(param->event);
  NumericMatrix resid(nevent, p);
  IntegerVector index(nevent);

  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  NumericVector a(p);       // s1(beta,k,t)
  NumericVector a2(p);      // sum of w*exp(zbeta)*z for the deaths

  NumericVector eta(param->nused);
  for (person = 0; person < param->nused; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  int istrata = param->strata[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  int j = nevent;
  NumericVector xbar(p);
  for (person = 0; person < param->nused; ) {
    if (param->strata[person] != istrata) { // hit a new stratum
      istrata = param->strata[person]; // reset temporary variables
      i1 = person;
      denom = 0;
      for (i=0; i<p; i++) {
        a[i] = 0;
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
        }
      } else {
        j--;
        resid(j,_) = param->z(person,_);
        index[j] = person;

        ndead++;
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
      risk = param->weight[p1]*exp(eta[p1]);
      denom -= risk;
      for (i=0; i<p; i++) {
        a[i] -= risk*param->z(p1,i);
      }
    }

    // add to the main terms
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        for (i=0; i<p; i++) {
          a[i] += a2[i];
          xbar[i] = a[i]/denom;
        }
      } else {
        xbar.fill(0);
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            xbar[i] += a[i]/denom;
          }
        }
        xbar = xbar/ndead;
      }

      for (k=0; k<ndead; k++) {
        for (i=0; i<p; i++) {
          resid(j+k,i) -= xbar[i];
        }
      }

      // reset for the next death time
      ndead = 0;
      denom2 = 0;
      for (i=0; i<p; i++) {
        a2[i] = 0;
      }
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
              const bool robust = 0,
              const bool est_basehaz = 1,
              const bool est_resid = 1,
              const bool firth = 0,
              const bool plci = 0,
              const double alpha = 0.05) {

  int h, i, j, k, n = data.nrows();
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
    repn.fill(1);
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
    stratumn.fill(1);
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

  if (is_true(all(eventn == 0))) {
    stop("at least 1 event is needed to fit the Cox model");
  }

  NumericMatrix zn(n,p);
  if (p > 0) {
    for (j=0; j<p; j++) {
      String zj = covariates[j];
      if (!hasVariable(data, zj)) {
        stop("data must contain the variables in covariates");
      }
      NumericVector u = data[zj];
      for (i=0; i<n; i++) {
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


  bool has_id = hasVariable(data, id);

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

  if (robust && has_time2 && !has_id) {
    stop("id is needed for counting process data with robust variance");
  }


  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  int method = meth == "efron" ? 1 : 0;


  // sort the data by rep
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
  if (p > 0) zn = subset_matrix_by_row(zn, order);

  // exclude observations with missing values
  LogicalVector sub(n,1);
  for (i=0; i<n; i++) {
    if ((repn[i] == NA_INTEGER) || (stratumn[i] == NA_INTEGER) ||
        (std::isnan(timen[i])) || (eventn[i] == NA_INTEGER) ||
        (std::isnan(weightn[i])) || (std::isnan(offsetn[i])) ||
        (idn[i] == NA_INTEGER)) {
      sub[i] = 0;
    }
    for (j=0; j<p; j++) {
      if (std::isnan(zn(i,j))) sub[i] = 0;
    }
  }

  order = which(sub);
  repn = repn[order];
  stratumn = stratumn[order];
  timen = timen[order];
  time2n = time2n[order];
  eventn = eventn[order];
  weightn = weightn[order];
  offsetn = offsetn[order];
  idn = idn[order];
  if (p > 0) zn = subset_matrix_by_row(zn, order);
  n = sum(sub);

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
  NumericMatrix regloglik(nreps,2);
  NumericVector scoretest(nreps);
  IntegerVector niter(nreps);

  IntegerVector rep0(nreps*p);
  StringVector par0(nreps*p);
  NumericVector beta0(nreps*p), sebeta0(nreps*p), rsebeta0(nreps*p);
  NumericMatrix vbeta0(nreps*p,p), rvbeta0(nreps*p,p);
  NumericVector lb0(nreps*p), ub0(nreps*p), prob0(nreps*p);
  StringVector clparm0(nreps*p);

  // baseline hazards
  IntegerVector drep(n, NA_INTEGER), dstratum(n);
  NumericVector dtime(n), dnrisk(n), dnevent(n), dhaz(n), dvarhaz(n);
  NumericMatrix dgradhaz(n,p);
  int n0 = 0;
  int bign0 = 0;
  double toler = 1e-12;

  DataFrame basehaz2;
  NumericVector resmart(n);

  for (h=0; h<nreps; h++) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());

    IntegerVector stratum1 = stratumn[q1];
    NumericVector time1 = timen[q1];
    NumericVector time21 = time2n[q1];
    IntegerVector event1 = eventn[q1];
    NumericVector weight1 = weightn[q1];
    NumericVector offset1 = offsetn[q1];
    IntegerVector id1 = idn[q1];

    NumericMatrix z1(n1,p);
    if (p > 0) z1 = subset_matrix_by_row(zn, q1);

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
      if (index1[i] == index2[i]) {
        ignore1[i] = 1;
      } else {
        ignore1[i] = 0;
      }
    }

    int nused = n1 - sum(ignore1);

    // sort by stopping time in descending order within each stratum
    IntegerVector order2 = seq(0, n1-1);
    std::sort(order2.begin(), order2.end(), [&](int i, int j) {
      return (ignore1[i] < ignore1[j]) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] < stratum1[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstop[i] > tstop[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstop[i] == tstop[j]) && (event1[i] < event1[j]));
    });

    IntegerVector stratum1a = stratum1[order2];
    NumericVector tstarta = tstart[order2];
    NumericVector tstopa = tstop[order2];
    IntegerVector event1a = event1[order2];
    NumericVector weight1a = weight1[order2];
    NumericVector offset1a = offset1[order2];
    IntegerVector id1a = id1[order2];
    IntegerVector ignore1a = ignore1[order2];
    NumericMatrix z1a;
    if (p > 0) z1a = subset_matrix_by_row(z1, order2);

    // sort by starting time in descending order within each stratum
    IntegerVector order1 = seq(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      return (ignore1a[i] < ignore1a[j]) ||
        ((ignore1a[i] == ignore1a[j]) && (stratum1a[i] < stratum1a[j])) ||
        ((ignore1a[i] == ignore1a[j]) && (stratum1a[i] == stratum1a[j]) &&
        (tstarta[i] > tstarta[j]));
    });

    coxparams param = {nused, stratum1a, tstarta, tstopa, event1a,
                       weight1a, offset1a, z1a, order1, method};

    NumericVector bint(p);


    // prepare the data for estimating the baseline hazards at all time points
    List basehaz1;

    // sort by stopping time in descending order within each stratum
    IntegerVector order2x = seq(0, n1-1);
    std::sort(order2x.begin(), order2x.end(), [&](int i, int j) {
      return (stratum1[i] > stratum1[j]) ||
        ((stratum1[i] == stratum1[j]) && (tstop[i] > tstop[j]));
    });

    IntegerVector stratum1x = stratum1[order2x];
    NumericVector tstartx = tstart[order2x];
    NumericVector tstopx = tstop[order2x];
    IntegerVector event1x = event1[order2x];
    NumericVector weight1x = weight1[order2x];
    NumericVector offset1x = offset1[order2x];
    IntegerVector id1x = id1[order2x];
    IntegerVector ignore1x = ignore1[order2x];
    NumericMatrix z1x(n1,p);
    if (p > 0) z1x = subset_matrix_by_row(z1, order2x);

    // sort by starting time in descending order within each stratum
    IntegerVector order1x = seq(0, n1-1);
    std::sort(order1x.begin(), order1x.end(), [&](int i, int j) {
      return (stratum1x[i] > stratum1x[j]) ||
        ((stratum1x[i] == stratum1x[j]) && (tstartx[i] > tstartx[j]));
    });

    coxparams paramx = {n1, stratum1x, tstartx, tstopx, event1x,
                        weight1x, offset1x, z1x, order1x, method};

    NumericVector b(p);
    NumericMatrix vb(p,p);
    if (p > 0) {
      IntegerVector colfit = seq(0,p-1);
      List out = phregloop(p, bint, &param, 30, 1.0e-9, firth, colfit, p);

      bool fail = out["fail"];
      if (fail) stop("The algorithm in phregr did not converge");

      b = out["coef"];
      vb = as<NumericMatrix>(out["var"]);

      NumericVector seb(p);
      for (j=0; j<p; j++) {
        seb[j] = sqrt(vb(j,j));
      }

      for (i=0; i<p; i++) {
        rep0[h*p+i] = h+1;
        par0[h*p+i] = covariates[i];
        beta0[h*p+i] = b[i];
        sebeta0[h*p+i] = seb[i];
        for (j=0; j<p; j++) {
          vbeta0(h*p+i,j) = vb(i,j);
        }
      }

      // score statistic
      NumericVector scoreint = f_score_2(p, bint, &param);
      NumericMatrix infobint = f_info_2(p, bint, &param);
      NumericMatrix vbint = invsympd(infobint, p, toler);
      for (i=0; i<p; i++) {
        for (j=0; j<p; j++) {
          scoretest[h] += scoreint[i]*vbint(i,j)*scoreint[j];
        }
      }

      niter[h] = out["iter"];

      // robust variance estimates
      NumericVector rseb(p);  // robust standard error for betahat
      if (robust) {
        NumericMatrix ressco = f_ressco_2(p, b, &param);

        int nr; // number of rows in the score residual matrix
        if (!has_id) {
          for (i=0; i<n1; i++) {
            for (j=0; j<p; j++) {
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
              resid(i,j) = ressco(order[i],j);
            }
          }

          NumericVector weight2 = weight1a[order];

          NumericMatrix ressco2(nids,p);
          for (i=0; i<nids; i++) {
            for (j=0; j<p; j++) {
              for (k=idx[i]; k<idx[i+1]; k++) {
                ressco2(i,j) += weight2[k]*resid(k,j);
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
              D(i,j) += ressco(i,k)*vb(k,j);
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

        for (i=0; i<p; i++) {
          rseb[i] = sqrt(rvb(i,i));
        }

        for (i=0; i<p; i++) {
          rsebeta0[h*p+i] = rseb[i];
          for (j=0; j<p; j++) {
            rvbeta0(h*p+i,j) = rvb(i,j);
          }
        }
      }

      // profile likelihood confidence interval for regression coefficients
      NumericVector lb(p), ub(p), prob(p);
      StringVector clparm(p);

      double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
      if (plci) {
        double lmax;
        if (firth) {
          lmax = f_pen_llik_2(p, b, &param);
        } else {
          lmax = f_llik_2(p, b, &param);
        }
        double l0 = lmax - 0.5*R::qchisq(1-alpha, 1, 1, 0);

        for (k=0; k<p; k++) {
          lb[k] = phregplloop(p, b, &param, 30, 1.0e-8, firth, k, -1, l0);
          ub[k] = phregplloop(p, b, &param, 30, 1.0e-8, firth, k, 1, l0);

          IntegerVector colfit1(p-1);
          for (i=0; i<k; i++) {
            colfit1[i] = i;
          }
          for (i=k+1; i<p; i++) {
            colfit1[i-1] = i;
          }

          NumericVector b0(p);
          List out0 = phregloop(p, b0, &param, 30, 1.0e-8, firth,
                                colfit1, p-1);
          double lmax0 = out0["loglik"];
          prob[k] = R::pchisq(-2*(lmax0 - lmax), 1, 0, 0);
          clparm[k] = "PL";
        }
      } else {
        for (k=0; k<p; k++) {
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

      for (i=0; i<p; i++) {
        lb0[h*p+i] = lb[i];
        ub0[h*p+i] = ub[i];
        prob0[h*p+i] = prob[i];
        clparm0[h*p+i] = clparm[i];
      }
    }

    // log-likelihoods
    if (firth) {
      loglik(h,0) = f_pen_llik_2(p, bint, &param);
      loglik(h,1) = f_pen_llik_2(p, b, &param);

      regloglik(h,0) = f_llik_2(p, bint, &param);
      regloglik(h,1) = f_llik_2(p, b, &param);
    } else {
      loglik(h,0) = f_llik_2(p, bint, &param);
      loglik(h,1) = f_llik_2(p, b, &param);
    }


    // estimate baseline hazard
    if (est_basehaz) {
      basehaz1 = f_basehaz(p, b, &paramx);

      IntegerVector dstratum1 = basehaz1["stratum"];
      NumericVector dtime1 = basehaz1["time"];
      NumericVector dnrisk1 = basehaz1["nrisk"];
      NumericVector dnevent1 = basehaz1["nevent"];
      NumericVector dhaz1 = basehaz1["haz"];
      NumericVector dvarhaz1 = basehaz1["varhaz"];
      int J = static_cast<int>(dstratum1.size());

      // recover original event times
      for (j=0; j<J; j++) {
        dtime1[j] = dtime1[j] - dstratum1[j]*delta;
      }

      // add to output data frame
      for (j=0; j<J; j++) {
        drep[n0+j] = h+1;
        dstratum[n0+j] = dstratum1[j];
        dtime[n0+j] = dtime1[j];
        dnrisk[n0+j] = dnrisk1[j];
        dnevent[n0+j] = dnevent1[j];
        dhaz[n0+j] = dhaz1[j];
        dvarhaz[n0+j] = dvarhaz1[j];

        if (p > 0) {
          NumericMatrix dgradhaz1 = basehaz1["gradhaz"];
          for (i=0; i<p; i++) {
            dgradhaz(n0+j,i) = dgradhaz1(j,i);
          }
        }
      }

      n0 += J;
    }

    // martingale residuals
    if (est_resid) {
      NumericVector resid = f_resmart(p, b, &param);

      for (i=0; i<n1; i++) {
        resmart[bign0 + order2[i]] = resid[i];
      }

      bign0 += n1;
    }
  }


  if (est_basehaz) {
    IntegerVector sub = which(!is_na(drep));
    drep = drep[sub];
    dstratum = dstratum[sub];
    dtime = dtime[sub];
    dnrisk = dnrisk[sub];
    dnevent = dnevent[sub];
    dhaz = dhaz[sub];
    dvarhaz = dvarhaz[sub];
    if (p > 0) dgradhaz = subset_matrix_by_row(dgradhaz, sub);
  }

  List result;
  if (p > 0) {
    NumericVector expbeta0 = exp(beta0);
    NumericVector z0(nreps*p);
    if (!robust) z0 = beta0/sebeta0;
    else z0 = beta0/rsebeta0;

    DataFrame sumstat = List::create(
      _["n"] = nobs,
      _["nevents"] = nevents,
      _["loglik0"] = loglik(_,0),
      _["loglik1"] = loglik(_,1),
      _["scoretest"] = scoretest,
      _["niter"] = niter,
      _["ties"] = meth,
      _["p"] = p,
      _["robust"] = robust,
      _["firth"] = firth);

    if (firth) {
      sumstat.push_back(regloglik(_,0), "loglik0_unpenalized");
      sumstat.push_back(regloglik(_,1), "loglik1_unpenalized");
    }

    DataFrame parest;
    if (!robust) {
      parest = DataFrame::create(
        _["param"] = par0,
        _["beta"] = beta0,
        _["sebeta"] = sebeta0,
        _["z"] = z0,
        _["expbeta"] = expbeta0,
        _["vbeta"] = vbeta0,
        _["lower"] = lb0,
        _["upper"] = ub0,
        _["p"] = prob0,
        _["method"] = clparm0);
    } else {
      parest = DataFrame::create(
        _["param"] = par0,
        _["beta"] = beta0,
        _["sebeta"] = rsebeta0,
        _["z"] = z0,
        _["expbeta"] = expbeta0,
        _["vbeta"] = rvbeta0,
        _["lower"] = lb0,
        _["upper"] = ub0,
        _["p"] = prob0,
        _["method"] = clparm0,
        _["sebeta_naive"] = sebeta0,
        _["vbeta_naive"] = vbeta0);
    }

    if (has_rep) {
      for (i=0; i<p_rep; i++) {
        String s = rep[i];
        if (TYPEOF(data[s]) == INTSXP) {
          IntegerVector repwi = u_rep[s];
          sumstat.push_back(repwi[rep01-1], s);
          parest.push_back(repwi[rep0-1], s);
        } else if (TYPEOF(data[s]) == REALSXP) {
          NumericVector repwn = u_rep[s];
          sumstat.push_back(repwn[rep01-1], s);
          parest.push_back(repwn[rep0-1], s);
        } else if (TYPEOF(data[rep]) == STRSXP) {
          StringVector repwc = u_rep[s];
          sumstat.push_back(repwc[rep01-1], s);
          parest.push_back(repwc[rep0-1], s);
        }
      }
    }

    result = List::create(
      _["sumstat"] = sumstat,
      _["parest"] = parest);

    if (est_basehaz) {
      DataFrame basehaz = DataFrame::create(
        _["time"] = dtime,
        _["nrisk"] = dnrisk,
        _["nevent"] = dnevent,
        _["haz"] = dhaz,
        _["varhaz"] = dvarhaz,
        _["gradhaz"] = dgradhaz
      );

      if (has_stratum) {
        for (i=0; i<p_stratum; i++) {
          String s = stratum[i];
          if (TYPEOF(data[s]) == INTSXP) {
            IntegerVector stratumwi = u_stratum[s];
            basehaz.push_back(stratumwi[dstratum-1], s);
          } else if (TYPEOF(data[s]) == REALSXP) {
            NumericVector stratumwn = u_stratum[s];
            basehaz.push_back(stratumwn[dstratum-1], s);
          } else if (TYPEOF(data[s]) == STRSXP) {
            StringVector stratumwc = u_stratum[s];
            basehaz.push_back(stratumwc[dstratum-1], s);
          }
        }
      }

      if (has_rep) {
        for (i=0; i<p_rep; i++) {
          String s = rep[i];
          if (TYPEOF(data[s]) == INTSXP) {
            IntegerVector repwi = u_rep[s];
            basehaz.push_back(repwi[drep-1], s);
          } else if (TYPEOF(data[s]) == REALSXP) {
            NumericVector repwn = u_rep[s];
            basehaz.push_back(repwn[drep-1], s);
          } else if (TYPEOF(data[rep]) == STRSXP) {
            StringVector repwc = u_rep[s];
            basehaz.push_back(repwc[drep-1], s);
          }
        }
      }

      result.push_back(basehaz, "basehaz");
    }

    if (est_resid) {
      result.push_back(resmart, "residuals");
    }
  } else {
    DataFrame sumstat = List::create(
      _["n"] = nobs,
      _["nevents"] = nevents,
      _["loglik0"] = loglik(_,0),
      _["loglik1"] = loglik(_,0),
      _["scoretest"] = 0,
      _["niter"] = niter,
      _["ties"] = meth,
      _["p"] = p,
      _["robust"] = robust,
      _["firth"] = firth);

    if (firth) {
      sumstat.push_back(regloglik(_,0), "loglik0_unpenalized");
      sumstat.push_back(regloglik(_,1), "loglik1_unpenalized");
    }

    if (has_rep) {
      for (i=0; i<p_rep; i++) {
        String s = rep[i];
        if (TYPEOF(data[s]) == INTSXP) {
          IntegerVector repwi = u_rep[s];
          sumstat.push_back(repwi[rep01-1], s);
        } else if (TYPEOF(data[s]) == REALSXP) {
          NumericVector repwn = u_rep[s];
          sumstat.push_back(repwn[rep01-1], s);
        } else if (TYPEOF(data[rep]) == STRSXP) {
          StringVector repwc = u_rep[s];
          sumstat.push_back(repwc[rep01-1], s);
        }
      }
    }

    result = List::create(
      _["sumstat"] = sumstat);

    if (est_basehaz) {
      DataFrame basehaz = DataFrame::create(
        _["time"] = dtime,
        _["nrisk"] = dnrisk,
        _["nevent"] = dnevent,
        _["haz"] = dhaz,
        _["varhaz"] = dvarhaz
      );

      if (has_stratum) {
        for (i=0; i<p_stratum; i++) {
          String s = stratum[i];
          if (TYPEOF(data[s]) == INTSXP) {
            IntegerVector stratumwi = u_stratum[s];
            basehaz.push_back(stratumwi[dstratum-1], s);
          } else if (TYPEOF(data[s]) == REALSXP) {
            NumericVector stratumwn = u_stratum[s];
            basehaz.push_back(stratumwn[dstratum-1], s);
          } else if (TYPEOF(data[s]) == STRSXP) {
            StringVector stratumwc = u_stratum[s];
            basehaz.push_back(stratumwc[dstratum-1], s);
          }
        }
      }

      if (has_rep) {
        for (i=0; i<p_rep; i++) {
          String s = rep[i];
          if (TYPEOF(data[s]) == INTSXP) {
            IntegerVector repwi = u_rep[s];
            basehaz.push_back(repwi[drep-1], s);
          } else if (TYPEOF(data[s]) == REALSXP) {
            NumericVector repwn = u_rep[s];
            basehaz.push_back(repwn[drep-1], s);
          } else if (TYPEOF(data[rep]) == STRSXP) {
            StringVector repwc = u_rep[s];
            basehaz.push_back(repwc[drep-1], s);
          }
        }
      }

      result.push_back(basehaz, "basehaz");
    }

    if (est_resid) {
      result.push_back(resmart, "residuals");
    }
  }

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
                           const bool sefit = 1,
                           const String conftype = "log-log",
                           const double conflev = 0.95) {

  int h, i, j, k, n0 = basehaz.nrows(), n = newdata.nrows();
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

  double zcrit = R::qnorm((1.0 + conflev)/2.0, 0, 1, 1, 0);

  IntegerVector stratumn0(n0);
  DataFrame u_stratum0;
  IntegerVector nlevels;
  List lookups;
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    stratumn0.fill(1);
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
  for (j=0; j<p; j++) {
    String zj = covariates[j];
    if (!hasVariable(newdata, zj)) {
      stop("newdata must contain the variables in covariates");
    }
    NumericVector u = newdata[zj];
    for (i=0; i<n; i++) {
      zn(i,j) = u[i];
    }
  }

  bool has_stratum;
  IntegerVector stratumn(n);
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(1);
  } else {
    has_stratum = 1;
    // match stratum in newdata to stratum in basehaz
    int orep = u_stratum0.nrow();
    for (i=0; i<p_stratum; i++) {
      orep = orep/nlevels[i];
      String s = stratum[i];
      IntegerVector idx;
      if (TYPEOF(newdata[s]) == LGLSXP || TYPEOF(newdata[s]) == INTSXP) {
        IntegerVector v = newdata[s];
        IntegerVector w = lookups[i];
        idx = match(v,w) - 1;
      } else if (TYPEOF(newdata[s]) == REALSXP) {
        NumericVector v = newdata[s];
        NumericVector w = lookups[i];
        idx = match(v,w) - 1;
      } else if (TYPEOF(newdata[s]) == STRSXP) {
        StringVector v = newdata[s];
        StringVector w = lookups[i];
        idx = match(v,w) - 1;
      }

      stratumn = stratumn + idx*orep;
    }

    stratumn = stratumn + 1;
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
  NumericVector haz0 = basehaz["haz"];
  NumericVector vhaz0 = basehaz["varhaz"];
  NumericMatrix ghaz0(n0,p);
  for (j=0; j<p; j++) {
    NumericVector u = basehaz[j+5]; // gradhaz starts at col 5 in basehaz
    ghaz0(_,j) = u;
  }


  // whether the input data has the counting process style of input
  bool has_id = hasVariable(newdata, id);

  // create the numeric id variable
  IntegerVector idn(n);
  IntegerVector idwi;
  NumericVector idwn;
  StringVector idwc;
  if (!has_id) {
    idn = seq(1,n);
  } else {
    if (TYPEOF(newdata[id]) == INTSXP) {
      IntegerVector idv = newdata[id];
      idwi = unique(idv);
      idwi.sort();
      idn = match(idv, idwi);
    } else if (TYPEOF(newdata[id]) == REALSXP) {
      NumericVector idv = newdata[id];
      idwn = unique(idv);
      idwn.sort();
      idn = match(idv, idwn);
    } else if (TYPEOF(newdata[id]) == STRSXP) {
      StringVector idv = newdata[id];
      idwc = unique(idv);
      idwc.sort();
      idn = match(idv, idwc);
    } else {
      stop("incorrect type for the id variable in newdata");
    }
  }

  NumericVector risk(n);
  NumericVector zbeta = clone(offsetn);
  for (j=0; j<p; j++) {
    zbeta = zbeta + beta[j]*zn(_,j);
  }
  risk = exp(zbeta);

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
    return ((idn[i] < idn[j]) ||
            ((idn[i] == idn[j]) && (tstopn[i] < tstopn[j])));
  });

  idn = idn[order];
  stratumn = stratumn[order];
  tstopn = tstopn[order];
  risk = risk[order];
  zn = subset_matrix_by_row(zn, order);

  // count number of observations for each id
  IntegerVector idx(1,0);
  for (i=1; i<n; i++) {
    if (idn[i] != idn[i-1]) {
      idx.push_back(i);
    }
  }

  int nids = static_cast<int>(idx.size());
  idx.push_back(n);

  int m = 0;
  for (h=0; h<nids; h++) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());
    IntegerVector stratum1 = stratumn[q1];
    NumericVector tstop1 = tstopn[q1];

    // match the stratum in basehaz
    IntegerVector idx1 = which(stratumn0 == stratum1[0]);
    NumericVector time1 = time0[idx1];

    // left-open and right-closed intervals containing the event time
    IntegerVector idx2 = rev(n1 - findInterval3(rev(-time1), rev(-tstop1)));
    int m1 = max(which(idx2 < n1));

    if (m1 != NA_INTEGER) {
      m1 = m1 + 1;
      m += m1;
    }
  }

  NumericVector time(m);
  NumericVector nrisk(m), nevent(m);
  NumericVector cumhaz(m), vcumhaz(m), secumhaz(m);
  IntegerVector strata(m);
  NumericMatrix z(m,p);
  IntegerVector ids(m);

  // process by id
  int l = 0;
  for (h=0; h<nids; h++) {
    IntegerVector q1 = Range(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());

    IntegerVector id1 = idn[q1];
    IntegerVector stratum1 = stratumn[q1];
    NumericVector tstop1 = tstopn[q1];
    NumericVector risk1 = risk[q1];
    NumericMatrix z1 = subset_matrix_by_row(zn, q1);

    // match the stratum in basehaz
    IntegerVector idx1 = which(stratumn0 == stratum1[0]);
    NumericVector time1 = time0[idx1];
    NumericVector nrisk1 = nrisk0[idx1];
    NumericVector nevent1 = nevent0[idx1];

    // left-open and right-closed intervals containing the event time
    IntegerVector idx2 = rev(n1 - findInterval3(rev(-time1), rev(-tstop1)));
    int m1 = max(which(idx2 < n1));

    if (m1 != NA_INTEGER) {
      m1 = m1 + 1;
      IntegerVector idx3 = idx1[Range(0, m1-1)];
      NumericVector haz1 = haz0[idx3];

      // cumulative hazards
      for (i=0; i<m1; i++) {
        time[l+i] = time1[i];
        nrisk[l+i] = nrisk1[i];
        nevent[l+i] = nevent1[i];

        k = idx2[i];
        ids[l+i] = id1[k];
        strata[l+i] = stratum1[k];
        for (j=0; j<p; j++) {
          z(l+i,j) = z1(k,j);
        }

        if (i==0) {
          cumhaz[l+i] = haz1[i]*risk1[k];
        } else {
          cumhaz[l+i] = cumhaz[l+i-1] + haz1[i]*risk1[k];
        }
      }

      if (sefit) {
        NumericVector vhaz1 = vhaz0[idx3];
        NumericMatrix ghaz1(m1,p);
        for (j=0; j<p; j++) {
          for (i=0; i<m1; i++) {
            ghaz1(i,j) = ghaz0(idx3[i],j);
          }
        }

        NumericMatrix a(m1,p);
        for (j=0; j<p; j++) {
          for (i=0; i<m1; i++) {
            k = idx2[i];
            if (i==0) {
              a(i,j) = (haz1[i]*z1(k,j) - ghaz1(i,j))*risk1[k];
            } else {
              a(i,j) = a(i-1,j) + (haz1[i]*z1(k,j) - ghaz1(i,j))*risk1[k];
            }
          }
        }

        // calculate the first component of variance
        for (i=0; i<m1; i++) {
          k = idx2[i];
          if (i==0) {
            vcumhaz[l+i] = vhaz1[i]*risk1[k]*risk1[k];
          } else {
            vcumhaz[l+i] = vcumhaz[l+i-1] + vhaz1[i]*risk1[k]*risk1[k];
          }
        }

        // add the second component of variance
        for (i=0; i<m1; i++) {
          for (j=0; j<p; j++) {
            for (k=0; k<p; k++) {
              vcumhaz[l+i] += a(i,j)*vbeta(j,k)*a(i,k);
            }
          }
          secumhaz[l+i] = sqrt(vcumhaz[l+i]);
        }
      }

      l += m1;
    }
  }

  NumericVector surv = exp(-cumhaz);

  DataFrame result = DataFrame::create(
    _["time"] = time,
    _["nrisk"] = nrisk,
    _["nevent"] = nevent,
    _["cumhaz"] = cumhaz,
    _["surv"] = surv);

  if (sefit) {
    NumericVector sesurv = surv*secumhaz;

    NumericVector lower(m), upper(m);
    for (i=0; i<m; i++) {
      NumericVector ci = fsurvci(surv[i], sesurv[i], ct, zcrit);
      lower[i] = ci[0];
      upper[i] = ci[1];
    }

    result.push_back(sesurv, "sesurv");
    result.push_back(lower, "lower");
    result.push_back(upper, "upper");
    result.push_back(conflev, "conflev");
    result.push_back(conftype, "conftype");
  }

  for (j=0; j<p; j++) {
    NumericVector u = z(_,j);
    String zj = covariates[j];
    result.push_back(u, zj);
  }

  if (has_stratum) {
    for (i=0; i<p_stratum; i++) {
      String s = stratum[i];
      if (TYPEOF(newdata[s]) == INTSXP) {
        IntegerVector stratumwi = u_stratum0[s];
        result.push_back(stratumwi[strata-1], s);
      } else if (TYPEOF(newdata[s]) == REALSXP) {
        NumericVector stratumwn = u_stratum0[s];
        result.push_back(stratumwn[strata-1], s);
      } else if (TYPEOF(newdata[s]) == STRSXP) {
        StringVector stratumwc = u_stratum0[s];
        result.push_back(stratumwc[strata-1], s);
      }
    }
  }

  if (has_id) {
    if (TYPEOF(newdata[id]) == INTSXP) {
      result.push_front(idwi[ids-1], id);
    } else if (TYPEOF(newdata[id]) == REALSXP) {
      result.push_front(idwn[ids-1], id);
    } else if (TYPEOF(newdata[id]) == STRSXP) {
      result.push_front(idwc[ids-1], id);
    }
  }

  return result;
}


// [[Rcpp::export]]
List residuals_phregcpp(const int p,
                        const NumericVector& beta,
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
                        const std::string type = "schoenfeld") {

  if (p == 0) {
    std::string errmsg = "The mode must include >=1 covariate to yield " +
      type + " residuals";
    stop(errmsg);
  }

  int i, j, n = data.nrows();

  bool has_stratum;
  IntegerVector stratumn(n);
  DataFrame u_stratum;
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(1);
  } else {
    has_stratum = 1;
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

  if (is_true(all(eventn == 0))) {
    stop("at least 1 event is needed to fit the Cox model");
  }

  NumericMatrix zn(n,p);
  if (p > 0) {
    for (j=0; j<p; j++) {
      String zj = covariates[j];
      if (!hasVariable(data, zj)) {
        stop("data must contain the variables in covariates");
      }
      NumericVector u = data[zj];
      for (i=0; i<n; i++) {
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

  bool has_id = hasVariable(data, id);

  IntegerVector idn(n);
  IntegerVector idwi(n);
  NumericVector idwn(n);
  StringVector idwc(n);
  if (!has_id) {
    idn = seq(1,n);
  } else {
    if (TYPEOF(data[id]) == INTSXP) {
      IntegerVector idv = data[id];
      idwi = unique(idv);
      idwi.sort();
      idn = match(idv, idwi);
    } else if (TYPEOF(data[id]) == REALSXP) {
      NumericVector idv = data[id];
      idwn = unique(idv);
      idwn.sort();
      idn = match(idv, idwn);
    } else if (TYPEOF(data[id]) == STRSXP) {
      StringVector idv = data[id];
      idwc = unique(idv);
      idwc.sort();
      idn = match(idv, idwc);
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
  NumericVector tstart(n), tstop(n);
  if (!has_time2) {
    tstop = timen;
  } else {
    tstart = timen;
    tstop = time2n;
  }

  // ignore subjects not at risk for any event time
  double delta = max(tstop) + 1.0; // ensure no overlap between strata
  for (i=0; i<n; i++) {
    tstart[i] = tstart[i] + stratumn[i]*delta;
    tstop[i] = tstop[i] + stratumn[i]*delta;
  }

  NumericVector etime = tstop[eventn==1];
  etime = unique(etime);
  etime.sort();

  IntegerVector index1 = findInterval3(tstart, etime);
  IntegerVector index2 = findInterval3(tstop, etime);
  IntegerVector ignore(n);
  for (i=0; i<n; i++) {
    if (index1[i] == index2[i]) {
      ignore[i] = 1;
    } else {
      ignore[i] = 0;
    }
  }
  int nused = n - sum(ignore);

  List result;
  NumericMatrix resid(n,p);
  if (type == "score") {
    // sort by stopping time in descending order within each stratum
    IntegerVector order2 = seq(0, n-1);
    std::sort(order2.begin(), order2.end(), [&](int i, int j) {
      return (ignore[i] < ignore[j]) ||
        ((ignore[i] == ignore[j]) && (stratumn[i] < stratumn[j])) ||
        ((ignore[i] == ignore[j]) && (stratumn[i] == stratumn[j]) &&
        (tstop[i] > tstop[j])) ||
        ((ignore[i] == ignore[j]) && (stratumn[i] == stratumn[j]) &&
        (tstop[i] == tstop[j]) && (eventn[i] < eventn[j]));
    });

    IntegerVector stratum1 = stratumn[order2];
    NumericVector tstart1 = tstart[order2];
    NumericVector tstop1 = tstop[order2];
    IntegerVector event1 = eventn[order2];
    NumericVector weight1 = weightn[order2];
    NumericVector offset1 = offsetn[order2];
    IntegerVector id1 = idn[order2];
    IntegerVector ignore1 = ignore[order2];
    NumericMatrix z1 = subset_matrix_by_row(zn, order2);

    // sort by starting time in descending order within each stratum
    IntegerVector order1 = seq(0, n-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      return (ignore1[i] < ignore1[j]) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] < stratum1[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstart1[i] > tstart1[j]));
    });

    coxparams param = {nused, stratum1, tstart1, tstop1, event1,
                       weight1, offset1, z1, order1, method};

    NumericMatrix ressco = f_ressco_2(p, beta, &param);

    for (i=0; i<n; i++) {
      resid(order2[i],_) = ressco(i,_);
    }

    result = List::create(
      Named("type") = type,
      Named("resid") = resid);

    if (has_id) {
      if (TYPEOF(data[id]) == INTSXP) {
        result.push_back(idwi[idn-1], id);
      } else if (TYPEOF(data[id]) == REALSXP) {
        result.push_back(idwn[idn-1], id);
      } else if (TYPEOF(data[id]) == STRSXP) {
        result.push_back(idwc[idn-1], id);
      }
    }
    result.push_back(idn, "obs");
  } else if (type == "schoenfeld") {
    // sort by stopping time in descending order within each stratum
    IntegerVector order2 = seq(0, n-1);
    std::sort(order2.begin(), order2.end(), [&](int i, int j) {
      return (ignore[i] < ignore[j]) ||
        ((ignore[i] == ignore[j]) && (stratumn[i] > stratumn[j])) ||
        ((ignore[i] == ignore[j]) && (stratumn[i] == stratumn[j]) &&
        (tstop[i] > tstop[j])) ||
        ((ignore[i] == ignore[j]) && (stratumn[i] == stratumn[j]) &&
        (tstop[i] == tstop[j]) && (eventn[i] < eventn[j])) ||
        ((ignore[i] == ignore[j]) && (stratumn[i] == stratumn[j]) &&
        (tstop[i] == tstop[j]) && (eventn[i] == eventn[j]) &&
        (idn[i] > idn[j]));
    });

    IntegerVector stratum1 = stratumn[order2];
    NumericVector tstart1 = tstart[order2];
    NumericVector tstop1 = tstop[order2];
    IntegerVector event1 = eventn[order2];
    NumericVector weight1 = weightn[order2];
    NumericVector offset1 = offsetn[order2];
    IntegerVector id1 = idn[order2];
    IntegerVector ignore1 = ignore[order2];
    NumericMatrix z1 = subset_matrix_by_row(zn, order2);

    // sort by starting time in descending order within each stratum
    IntegerVector order1 = seq(0, n-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      return (ignore1[i] < ignore1[j]) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] > stratum1[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstart1[i] > tstart1[j]));
    });

    coxparams param = {nused, stratum1, tstart1, tstop1, event1,
                       weight1, offset1, z1, order1, method};

    List out = f_ressch(p, beta, &param);

    resid = as<NumericMatrix>(out["resid"]);

    IntegerVector index = out["index"];
    IntegerVector stratum2 = stratum1[index];
    NumericVector time2x = tstop1[index];

    // recover original event times
    int J = static_cast<int>(stratum2.size());
    for (j=0; j<J; j++) {
      time2x[j] = time2x[j] - stratum2[j]*delta;
    }

    IntegerVector id2 = id1[index];

    result = List::create(
      Named("type") = type,
      Named("resid") = resid,
      Named("time") = time2x);

    if (has_id) {
      if (TYPEOF(data[id]) == INTSXP) {
        result.push_back(idwi[id2-1], id);
      } else if (TYPEOF(data[id]) == REALSXP) {
        result.push_back(idwn[id2-1], id);
      } else if (TYPEOF(data[id]) == STRSXP) {
        result.push_back(idwc[id2-1], id);
      }
    }
    result.push_back(id2, "obs");

    if (has_stratum) {
      for (i=0; i<p_stratum; i++) {
        String s = stratum[i];
        if (TYPEOF(data[s]) == INTSXP) {
          IntegerVector stratumwi = u_stratum[s];
          result.push_back(stratumwi[stratum2-1], s);
        } else if (TYPEOF(data[s]) == REALSXP) {
          NumericVector stratumwn = u_stratum[s];
          result.push_back(stratumwn[stratum2-1], s);
        } else if (TYPEOF(data[s]) == STRSXP) {
          StringVector stratumwc = u_stratum[s];
          result.push_back(stratumwc[stratum2-1], s);
        }
      }
    }
    result.push_back(stratum2, "stratumn");

  }

  return result;
}
