#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;


DataFrame untreated(
    const double psi,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const bool recensor,
    const bool autoswitch) {

  NumericVector u = time*((1 - rx) + rx*exp(psi));
  NumericVector t_star = clone(u);
  IntegerVector d_star = clone(event);

  if (recensor) {
    NumericVector c_star = pmin(censor_time, censor_time*exp(psi));

    if (autoswitch) {
      NumericVector rx1 = rx[treat == 1];
      NumericVector rx0 = rx[treat == 0];
      if (is_true(all(rx1 == 1))) c_star[treat == 1] = R_PosInf;
      if (is_true(all(rx0 == 0))) c_star[treat == 0] = R_PosInf;
    }

    t_star = pmin(u, c_star);
    d_star[c_star < u] = 0;
  }

  DataFrame result = DataFrame::create(
    Named("t_star") = t_star,
    Named("d_star") = d_star
  );

  return result;
}


double est_eqn(
    const double psi,
    const IntegerVector& stratum,
    const NumericVector& time,
    const IntegerVector& event,
    const IntegerVector& treat,
    const NumericVector& rx,
    const NumericVector& censor_time,
    const double treat_modifier,
    const bool recensor,
    const bool autoswitch,
    double target = 0) {

  DataFrame Sstar = untreated(psi*treat_modifier, time, event, treat, rx,
                              censor_time, recensor, autoswitch);

  NumericVector t_star = Sstar["t_star"];
  IntegerVector d_star = Sstar["d_star"];

  DataFrame data = DataFrame::create(
    Named("stratum") = stratum,
    Named("treat") = treat,
    Named("time") = t_star,
    Named("event") = d_star);

  DataFrame df = lrtest(data, "none", "stratum", "treat", "time", "event",
                        0, 0);

  double result = as<double>(df["logRankZ"]) - target;
  return result;
}


//' @title Rank preserving structured failure time model (RPSFTM) for
//' treatment switching
//' @description Obtains the causal parameter estimate of the RPSFTM from
//' the log-rank test and the hazard ratio estimate from the Cox model.
//'
//' @param data The input data frame that contains the following variables:
//'
//'   * \code{stratum}: The stratum.
//'
//'   * \code{time}: The survival time for right censored data.
//'
//'   * \code{event}: The event indicator, 1=event, 0=no event.
//'
//'   * \code{treat}: The randomized treatment indicator, 1=treatment,
//'     0=control.
//'
//'   * \code{rx}: The proportion of time on active treatment.
//'
//'   * \code{censor_time}: The administrative censoring time. It should
//'     be provided for all subjects including those who had events.
//'
//'   * \code{base_cov}: The values of baseline covariates.
//'     This is the full-rank design matrix (excluding treat)
//'     for the Cox model, assuming that factor variables
//'     have already been expanded into dummy variables.
//'
//' @param stratum The name of the stratum variable in the input data.
//' @param time The name of the time variable in the input data.
//' @param event The name of the event variable in the input data.
//' @param treat The name of the treatment variable in the input data.
//' @param rx The name of the rx variable in the input data.
//' @param censor_time The name of the censor_time variable in the input data.
//' @param base_cov The vector of names of baseline covariates (excluding
//'   treat) in the input data.
//' @param low_psi The lower limit of the causal parameter of RPSFTM.
//' @param hi_psi The upper limit of the causal parameter of RPSFTM.
//' @param n_eval_z The number of points between low_psi and hi_psi at which
//'   to evaluate the log-rank Z-statistics.
//' @param alpha The significance level to calculate confidence intervals.
//' @param treat_modifier The optional sensitivity parameter for the
//'   constant treatment effect assumption.
//' @param recensor Whether to apply recensoring to counter-factual
//'   survival times. Defaults to \code{TRUE}.
//' @param autoswitch Whether to exclude recensoring for treatment arms
//'   with no switching. Defaults to \code{TRUE}.
//' @param gridsearch Whether to use grid search to estimate the causal
//'   parameter psi. Defaults to \code{FALSE}, in which case, a root
//'   finding algorithm will be used.
//' @param boot Whether to use bootstrap to obtain the confidence
//'   interval for hazard ratio. Defaults to \code{FALSE}, in which case,
//'   the confidence interval will be constructed to match the log-rank
//'   test p-value.
//' @param n_boot The number of bootstrap samples.
//'
//' @details We use the following steps to obtain the hazard ratio estimate
//' and confidence interval had there been no treatment switching:
//'
//' * use RPSFTM to estimate the causal parameter psi based on the log-rank
//'   test for counter-factual untreated survival times for both arms:
//'   \eqn{U = T_{off} + T_{on} e^{\psi}}.
//'
//' * Fit the Cox proportional hazards model to the observed survival times
//'   on the treatment arm and the counter-factual untreated survival times
//'   on the control arm to obtain the hazard ratio estimate.
//'
//' * Use either the log-rank test p-value for the treatment policy strategy
//'   or bootstrap to construct the confidence interval for hazard ratio.
//'
//' @return A list with the following components:
//'
//' * \code{psi}: The estimated causal parameter for RPSFTM.
//'
//' * \code{psi_CI}: The confidence interval for psi.
//'
//' * \code{psi_type}: The type of psi estimate, either "grid search" or
//'   "root finding".
//'
//' * \code{Sstar}: A data frame containing the counter-factual untreated
//'   survival times and the event indicators.
//'
//' * \code{kmstar}: A data frame containing the Kaplan-Meier estimates
//'   based on the counter-factual untreated survival times by treatment arm.
//'
//' * \code{eval_z}: A data frame containing the log-rank test Z-statistics
//'   evaluated at a sequence of psi values. Used to plot and to check
//'   if the range of psi values to search for the solution and
//'   limits of confidence interval of psi need be modified.
//'
//' * \code{pvalue}: The p-value of the log-rank test based on
//'   the treatment policy strategy.
//'
//' * \code{hr}: The estimated hazard ratio from the Cox model.
//'
//' * \code{hr_CI}: The confidence interval for hazard ratio.
//'
//' * \code{hr_CI_type}: The type of confidence interval for hazard ratio,
//'   either "log-rank p-value" or "bootstrap quantile".
//'
//' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
//'
//' @examples
//'
//' library(dplyr)
//'
//' data <- immdef %>% mutate(rx = 1-xoyrs/progyrs)
//'
//' fit <- rpsft(data, time = "progyrs", event = "prog", treat = "imm",
//'              rx = "rx", censor_time = "censyrs", boot = 0)
//'
//' c(fit$hr, fit$hr_CI)
//'
//' @export
// [[Rcpp::export]]
List rpsft(
    const DataFrame data,
    const std::string stratum = "stratum",
    const std::string time = "time",
    const std::string event = "event",
    const std::string treat = "treat",
    const std::string rx = "rx",
    const std::string censor_time = "censor_time",
    const StringVector& base_cov = "none",
    const double low_psi = -1,
    const double hi_psi = 1,
    const int n_eval_z = 100,
    const double alpha = 0.05,
    const double treat_modifier = 1,
    const bool recensor = 1,
    const bool autoswitch = 1,
    const bool gridsearch = 0,
    const bool boot = 0,
    const int n_boot = 1000) {

  int i, j, k, l, n = data.nrow();

  int p = static_cast<int>(base_cov.size());
  if (p == 1 && base_cov[0] == "none") p = 0;


  bool has_stratum = hasVariable(data, stratum);
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);
  bool has_treat = hasVariable(data, treat);
  bool has_rx = hasVariable(data, rx);
  bool has_censor_time = hasVariable(data, censor_time);


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


  if (!has_time) {
    stop("data must contain the time variable");
  }

  if (TYPEOF(data[time]) != INTSXP && TYPEOF(data[time]) != REALSXP) {
    stop("time must take numeric values");
  }

  NumericVector timen = data[time];
  if (is_true(any(timen <= 0))) {
    stop("time must be positive");
  }


  if (!has_event) {
    stop("data must contain the event variable");
  }

  if (TYPEOF(data[event]) != INTSXP) {
    stop("event must take integer values");
  }

  IntegerVector eventn = data[event];
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0");
  }

  if (is_true(all(eventn == 0))) {
    stop("at least 1 event is needed");
  }


  if (!has_treat) {
    stop("data must contain the treat variable");
  }

  // create the numeric treat variable
  IntegerVector treatn(n);
  if (TYPEOF(data[treat]) == LGLSXP) {
    IntegerVector treatv = data[treat];
    IntegerVector treatwi = unique(treatv);
    if (treatwi.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    treatwi = IntegerVector::create(1,0);
    treatn = 2 - treatv;
  } else if (TYPEOF(data[treat]) == INTSXP) {
    IntegerVector treatv = data[treat];
    IntegerVector treatwi = unique(treatv);
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
    NumericVector treatwn = unique(treatv);
    if (treatwn.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    treatwn.sort();
    treatn = match(treatv, treatwn);
  } else if (TYPEOF(data[treat]) == STRSXP) {
    StringVector treatv = data[treat];
    StringVector treatwc = unique(treatv);
    if (treatwc.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    treatwc.sort();
    treatn = match(treatv, treatwc);
  } else {
    stop("incorrect type for the treat variable in the input data");
  }

  treatn = 2 - treatn; // use the 1/0 treatment coding


  if (!has_rx) {
    stop("data must contain the rx variable");
  }

  if (TYPEOF(data[rx]) != INTSXP && TYPEOF(data[rx]) != REALSXP) {
    stop("rx must take numeric values");
  }

  NumericVector rxn = data[rx];
  if (is_true(any((rxn < 0) | (rxn > 1)))) {
    stop("rx must take values between 0 and 1");
  }


  if (!has_censor_time) {
    stop("data must contain the censor_time variable");
  }

  if (TYPEOF(data[censor_time]) != INTSXP &&
      TYPEOF(data[censor_time]) != REALSXP) {
    stop("censor_time must take numeric values");
  }

  NumericVector censor_timen = data[censor_time];
  if (is_true(any(censor_timen < timen))) {
    stop("censor_time must be greater than or equal to time");
  }


  if (low_psi >= hi_psi) {
    stop("low_psi must be less than hi_psi");
  }


  if (n_eval_z < 2) {
    stop("n_eval_z must be greater than or equal to 2");
  }


  if (alpha <= 0 || alpha >= 0.5) {
    stop("alpha must lie between 0 and 0.5");
  }


  if (treat_modifier <= 0) {
    stop("treat_modifier must be positive");
  }


  if (n_boot < 100) {
    stop("n_boot must be greater than or equal to 100");
  }


  DataFrame lr = lrtest(data, "none", stratum, treat, time, event, 0, 0);
  double logRankPValue = as<double>(lr["logRankPValue"]);


  // evaluate the log-rank test statistic at each psi
  double step_psi = (hi_psi - low_psi)/(n_eval_z - 1);
  NumericVector psi(n_eval_z), Z(n_eval_z);
  for (i=0; i<n_eval_z; i++) {
    psi[i] = low_psi + i*step_psi;
    Z[i] = est_eqn(psi[i], stratumn, timen, eventn, treatn, rxn,
                   censor_timen, treat_modifier, recensor,
                   autoswitch, 0);
  }

  DataFrame eval_z = DataFrame::create(
    Named("psi") = psi,
    Named("Z") = Z);


  // obtain the estimate and confidence interval of psi
  double psihat, psilower, psiupper;
  String psi_type;
  if (gridsearch) {
    auto f = [psi, Z](double target)->double{
      NumericVector Z1 = Z - target;
      NumericVector Zsq = Z1*Z1;
      return psi[which_min(Zsq)];
    };

    psihat = f(0);
    psilower = f(R::qnorm(1-alpha/2, 0, 1, 1, 0));
    psiupper = f(R::qnorm(alpha/2, 0, 1, 1, 0));
    psi_type = "grid search";
  } else {
    double target = 0;
    auto f = [stratumn, timen, eventn, treatn, rxn, censor_timen,
              treat_modifier, recensor, autoswitch,
              &target](double x)->double{
                return est_eqn(x, stratumn, timen, eventn, treatn, rxn,
                               censor_timen, treat_modifier, recensor,
                               autoswitch, target);
              };

    psihat = brent(f, low_psi, hi_psi, 1.0e-6);
    target = R::qnorm(1-alpha/2, 0, 1, 1, 0);
    psilower = brent(f, low_psi, psihat, 1.0e-6);
    target = R::qnorm(alpha/2, 0, 1, 1, 0);
    psiupper = brent(f, psihat, hi_psi, 1.0e-6);
    psi_type = "root finding";
  }


  // construct the counter-factual untreated survival times
  double psi1 = psihat*treat_modifier;
  DataFrame Sstar = untreated(psi1, timen, eventn, treatn, rxn,
                              censor_timen, recensor, autoswitch);

  NumericVector t_star = Sstar["t_star"];
  IntegerVector d_star = Sstar["d_star"];


  // obtain the Kaplan-Meier estimates
  DataFrame kmdata = DataFrame::create(
    Named("treat") = treatn,
    Named("time") = t_star,
    Named("event") = d_star);

  DataFrame kmstar = kmest(kmdata, "none", "treat", "time", "event",
                           "log-log", 1-alpha);


  // run Cox model to obtain the hazard ratio estimate
  NumericVector rx1 = rxn[treatn == 1];
  bool switch1 = is_false(all(rx1 == 1)); // treated arm has switcher(s)
  double t_tilde, c_tilde;

  NumericVector timewn(n);
  IntegerVector eventwn(n);
  for (i=0; i<n; i++) {
    if (treatn[i] == 1) {
      if (switch1) {
        // counter-factual survival time on treatment
        t_tilde = timen[i]*((1 - rxn[i])*exp(-psi1) + rxn[i]);
        if (recensor) {
          c_tilde = std::min(censor_timen[i], censor_timen[i]*exp(-psi1));
          timewn[i] = std::min(t_tilde, c_tilde);
          eventwn[i] = c_tilde < t_tilde ? 0 : eventn[i];
        } else {
          timewn[i] = t_tilde;
          eventwn[i] = eventn[i];
        }
      } else {
        timewn[i] = timen[i];
        eventwn[i] = eventn[i];
      }
    } else {
      timewn[i] = t_star[i];
      eventwn[i] = d_star[i];
    }
  }

  DataFrame phdata = DataFrame::create(
    Named("stratum") = stratumn,
    Named("time") = timewn,
    Named("event") = eventwn,
    Named("treat") = treatn);

  StringVector covariates(p+1);
  NumericMatrix zn(n,p+1);
  covariates[0] = "treat";
  zn(_,0) = treatn;
  for (j=0; j<p; j++) {
    String zj = base_cov[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in base_cov");
    }
    if (zj == treat) {
      stop("treat should be excluded from base_cov");
    }
    NumericVector u = data[zj];
    phdata.push_back(u, zj);
    covariates[j+1] = zj;
    zn(_,j+1) = u;
  }

  List fit = phregr(phdata, "none", "stratum", "time", "none",
                    "event", covariates, "none", "none", "efron", 0);

  DataFrame parest = DataFrame(fit["parest"]);
  StringVector param = parest["param"];
  NumericVector beta = parest["beta"];
  double loghrhat = 0;
  for (j=0; j<=p; j++) {
    if (param[j] == "treat") {
      loghrhat = beta[j];
      break;
    }
  }

  double hrhat = exp(loghrhat);


  // construct the confidence interval for HR
  double hrlower, hrupper;
  String hr_CI_type;
  if (!boot) { // use log-rank p-value to construct CI for HR if no boot
    double zcox = R::qnorm(logRankPValue, 0, 1, 1, 0);
    double seloghr = loghrhat/zcox;
    double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
    hrlower = exp(loghrhat - zcrit*seloghr);
    hrupper = exp(loghrhat + zcrit*seloghr);
    hr_CI_type = "log-rank p-value";
  } else { // bootstrap the entire process to construct CI for HR
    NumericVector hrhats(n_boot);
    IntegerVector stratumb(n), eventb(n), treatb(n);
    NumericVector timeb(n), rxb(n), censor_timeb(n);
    NumericMatrix zb(n,p+1);

    for (k=0; k<n_boot; k++) {

      // sample the data with replacement
      for (i=0; i<n; i++) {
        double u = R::runif(0,1);
        j = static_cast<int>(floor(u*n));

        stratumb[i] = stratumn[j];
        timeb[i] = timen[j];
        eventb[i] = eventn[j];
        treatb[i] = treatn[j];
        rxb[i] = rxn[j];
        censor_timeb[i] = censor_timen[j];

        for (l=0; l<=p; l++) {
          zb(i,l) = zn(j,l);
        }
      }


      // obtain psihat
      double psihat;
      if (gridsearch) {
        NumericVector Z(n_eval_z);
        for (i=0; i<n_eval_z; i++) {
          Z[i] = est_eqn(psi[i], stratumb, timeb, eventb, treatb, rxb,
                         censor_timeb, treat_modifier, recensor,
                         autoswitch, 0);
        }

        NumericVector Zsq = Z*Z;
        psihat = psi[which_min(Zsq)];
      } else {
        auto f = [stratumb, timeb, eventb, treatb, rxb, censor_timeb,
                  treat_modifier, recensor, autoswitch](double x)->double{
                    return est_eqn(x, stratumb, timeb, eventb, treatb, rxb,
                                   censor_timeb, treat_modifier, recensor,
                                   autoswitch, 0);
                  };

        psihat = brent(f, low_psi, hi_psi, 1.0e-6);
      }


      // run Cox model to obtain the hazard ratio estimate
      double psi1 = psihat*treat_modifier;
      DataFrame Sstar = untreated(psi1, timeb, eventb, treatb, rxb,
                                  censor_timeb, recensor, autoswitch);

      NumericVector t_star = Sstar["t_star"];
      IntegerVector d_star = Sstar["d_star"];

      NumericVector rx1 = rxb[treatb == 1];
      bool switch1 = is_false(all(rx1 == 1));

      NumericVector timewb(n);
      IntegerVector eventwb(n);
      for (i=0; i<n; i++) {
        if (treatb[i] == 1) {
          if (switch1) {
            t_tilde = timeb[i]*((1 - rxb[i])*exp(-psi1) + rxb[i]);
            if (recensor) {
              c_tilde = std::min(censor_timeb[i], censor_timeb[i]*exp(-psi1));
              timewb[i] = std::min(t_tilde, c_tilde);
              eventwb[i] = c_tilde < t_tilde ? 0 : eventb[i];
            } else {
              timewb[i] = t_tilde;
              eventwb[i] = eventb[i];
            }
          } else {
            timewb[i] = timeb[i];
            eventwb[i] = eventb[i];
          }
        } else {
          timewb[i] = t_star[i];
          eventwb[i] = d_star[i];
        }
      }

      DataFrame phdata = DataFrame::create(
        Named("stratum") = stratumb,
        Named("time") = timewb,
        Named("event") = eventwb,
        Named("treat") = treatb);

      for (j=0; j<p; j++) {
        String zj = base_cov[j];
        NumericVector u(n);
        for (i=0; i<n; i++) u[i] = zb(i,j+1);
        phdata.push_back(u, zj);
      }

      List fit = phregr(phdata, "none", "stratum", "time", "none",
                        "event", covariates, "none", "none", "efron", 0);

      DataFrame parest = DataFrame(fit["parest"]);
      StringVector param = parest["param"];
      NumericVector beta = parest["beta"];
      double loghrhat = 0;
      for (j=0; j<=p; j++) {
        if (param[j] == "treat") {
          loghrhat = beta[j];
          break;
        }
      }

      hrhats[k] = exp(loghrhat);
    }


    // obtain bootstrap confidence interval for HR
    hrlower = quantilecpp(hrhats, alpha/2);
    hrupper = quantilecpp(hrhats, 1-alpha/2);
    hr_CI_type = "bootstrap quantile";
  }


  List result = List::create(
    Named("psi") = psihat,
    Named("psi_CI") = NumericVector::create(psilower, psiupper),
    Named("psi_type") = psi_type,
    Named("Sstar") = Sstar,
    Named("kmstar") = kmstar,
    Named("eval_z") = eval_z,
    Named("pvalue") = logRankPValue,
    Named("hr") = hrhat,
    Named("hr_CI") = NumericVector::create(hrlower, hrupper),
    Named("hr_CI_type") = hr_CI_type);

  return result;
}
