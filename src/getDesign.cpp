#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;

//' @title Get group sequential design
//' @description Obtains the drift parameter and stopping boundaries for a
//' generic group sequential design assuming a constant treatment effect.
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
//'
//' @return A list of S3 class \code{design}.
//'
//' @examples
//'
//' getDesign(beta = 0.2,
//'           kMax = 2,
//'           informationRates = c(0.5,1),
//'           alpha = 0.025,
//'           typeAlphaSpending = "sfOF",
//'           typeBetaSpending = "sfP")
//'
//' @export
// [[Rcpp::export]]
List getDesign(const double beta = 0.2,
               const int kMax = NA_INTEGER,
               const NumericVector& informationRates = NA_REAL,
               const LogicalVector& efficacyStopping = NA_LOGICAL,
               const LogicalVector& futilityStopping = NA_LOGICAL,
               const NumericVector& criticalValues = NA_REAL,
               const double alpha = 0.025,
               const String typeAlphaSpending = "sfOF",
               const double parameterAlphaSpending = NA_REAL,
               const NumericVector& userAlphaSpending = NA_REAL,
               const NumericVector& futilityBounds = NA_REAL,
               const String typeBetaSpending = "none",
               const double parameterBetaSpending = NA_REAL,
               const NumericVector& userBetaSpending = NA_REAL) {
  
  
  NumericVector informationRates1 = clone(informationRates);
  LogicalVector efficacyStopping1 = clone(efficacyStopping);
  LogicalVector futilityStopping1 = clone(futilityStopping);
  NumericVector criticalValues1 = clone(criticalValues);
  NumericVector futilityBounds1 = clone(futilityBounds);
  
  String asf = typeAlphaSpending;
  double asfpar = parameterAlphaSpending;
  
  String bsf = typeBetaSpending;
  double bsfpar = parameterBetaSpending;
  
  
  
  if (R_isnancpp(beta)) {
    stop("beta must be provided");
  }
  
  if (beta >= 1-alpha || beta < 0.0001) {
    stop("beta must lie in [0.0001, 1-alpha)");
  }
  
  
  if (kMax < 1) {
    stop("kMax must be a positive integer");
  }
  
  if (is_false(any(is_na(informationRates)))) {
    if (informationRates.size() != kMax) {
      stop("Invalid length for informationRates");
    } else if (informationRates[0] <= 0) {
      stop("Elements of informationRates must be positive");
    } else if (kMax > 1 && is_true(any(diff(informationRates) <= 0))) {
      stop("Elements of informationRates must be increasing");
    } else if (informationRates[kMax-1] != 1) {
      stop("informationRates must end with 1");
    }
  } else {
    IntegerVector tem = seq_len(kMax);
    informationRates1 = as<NumericVector>(tem)/(kMax+0.0);
  }
  
  if (is_false(any(is_na(efficacyStopping)))) {
    if (efficacyStopping.size() != kMax) {
      stop("Invalid length for efficacyStopping");
    } else if (efficacyStopping[kMax-1] != 1) {
      stop("efficacyStopping must end with 1");
    }
  } else {
    efficacyStopping1 = rep(1, kMax);
  }
  
  if (is_false(any(is_na(futilityStopping)))) {
    if (futilityStopping.size() != kMax) {
      stop("Invalid length for futilityStopping");
    } else if (futilityStopping[kMax-1] != 1) {
      stop("futilityStopping must end with 1");
    }
  } else {
    futilityStopping1 = rep(1, kMax);
  }
  
  if (is_false(any(is_na(criticalValues)))) {
    if (criticalValues.size() != kMax) {
      stop("Invalid length for criticalValues");
    }
  }
  
  if (!R_isnancpp(alpha)) {
    if (alpha < 0.00001 || alpha >= 0.5) {
      stop("alpha must lie in [0.00001, 0.5)");
    }
  }
  
  if (is_true(any(is_na(criticalValues))) && !(asf=="OF" || asf=="P" ||
      asf=="WT" || asf=="sfOF" || asf=="sfP" ||
      asf=="sfKD" || asf=="sfHSD" || asf=="user" || asf=="none")) {
    stop("Invalid type for alpha spending");
  }
  
  if ((asf=="WT" || asf=="sfKD" || asf=="sfHSD") && R_isnancpp(asfpar)) {
    stop("Missing parameter for the alpha spending function");
  }
  
  if (asf=="sfKD" && asfpar <= 0) {
    stop ("parameterAlphaSpending must be positive for sfKD");
  }
  
  if (asf=="user") {
    if (is_true(any(is_na(userAlphaSpending)))) {
      stop("userAlphaSpending must be specified");
    } else if (userAlphaSpending.size() < kMax) {
      stop("Insufficient length of userAlphaSpending");
    } else if (userAlphaSpending[0] < 0) {
      stop("Elements of userAlphaSpending must be nonnegnative");
    } else if (kMax > 1 && is_true(any(diff(userAlphaSpending) < 0))) {
      stop("Elements of userAlphaSpending must be nondecreasing");
    } else if (userAlphaSpending[kMax-1] != alpha) {
      stop("userAlphaSpending must end with specified alpha");
    }
  }
  
  if (is_false(any(is_na(futilityBounds)))) {
    if (!(futilityBounds.size() == kMax-1 ||
        futilityBounds.size() == kMax)) {
      stop("Invalid length for futilityBounds");
    }
  }
  
  if (is_false(any(is_na(criticalValues))) &&
      is_false(any(is_na(futilityBounds)))) {
    for (int i=0; i<kMax-1; i++) {
      if (futilityBounds[i] > criticalValues[i]) {
        stop("futilityBounds must lie below criticalValues");
      }
    }
    
    if (futilityBounds.size() == kMax &&
        futilityBounds[kMax-1] != criticalValues[kMax-1]) {
      stop("futilityBounds and criticalValues must meet at final analysis");
    }
  }
  
  if (is_true(any(is_na(futilityBounds))) && !(bsf=="sfOF" || bsf=="sfP" ||
      bsf=="sfKD" || bsf=="sfHSD" || bsf=="user" || bsf=="none")) {
    stop("Invalid type for beta spending");
  }
  
  if ((bsf=="sfKD" || bsf=="sfHSD") && R_isnancpp(bsfpar)) {
    stop("Missing parameter for the beta spending function");
  }
  
  if (bsf=="sfKD" && bsfpar <= 0) {
    stop ("parameterBetaSpending must be positive for sfKD");
  }
  
  if (bsf=="user") {
    if (is_true(any(is_na(userBetaSpending)))) {
      stop("userBetaSpending must be specified");
    } else if (userBetaSpending.size() < kMax) {
      stop("Insufficient length of userBetaSpending");
    } else if (userBetaSpending[0] < 0) {
      stop("Elements of userBetaSpending must be nonnegnative");
    } else if (kMax > 1 && is_true(any(diff(userBetaSpending) < 0))) {
      stop("Elements of userBetaSpending must be nondecreasing");
    } else if (userBetaSpending[kMax-1] != beta) {
      stop("userBetaSpending must end with specified beta");
    }
  }
  
  
  
  bool missingCriticalValues = is_true(any(is_na(criticalValues)));
  bool missingFutilityBounds = is_true(any(is_na(futilityBounds)));
  
  if (missingCriticalValues) {
    criticalValues1 = rep(6.0, kMax);
    
    NumericVector theta(kMax); // mean values under H0, initialized to zero
    NumericVector t = clone(informationRates1); // information time
    
    
    if (asf == "none") {
      for (int i=0; i<kMax-1; i++) {
        criticalValues1[i] = 6.0;
      }
      criticalValues1[kMax-1] = R::qnorm(1-alpha, 0, 1, 1, 0);
    } else if (asf == "OF" || asf == "P" || asf == "WT") {
      double Delta;
      if (asf == "OF") {
        Delta = 0;
      } else if (asf == "P") {
        Delta = 0.5;
      } else {
        Delta = asfpar;
      }
      
      auto f = [kMax, alpha, Delta, theta, t, 
                efficacyStopping1] (double aval)->double {
                  NumericVector u(kMax), l(kMax);
                  for (int i=0; i<kMax; i++) {
                    u[i] = aval*pow(t[i], Delta-0.5);
                    if (!efficacyStopping1[i]) u[i] = 6.0;
                    l[i] = -6.0;
                  }
                  
                  List probs = exitprob(u, l, theta, t);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - alpha;
                };
      
      double cwt = brent(f, 0, 10, 1e-6);
      for (int i=0; i<kMax; i++) {
        criticalValues1[i] = cwt*pow(t[i], Delta-0.5);
        if (!efficacyStopping1[i]) criticalValues1[i] = 6.0;
      }
    } else if (asf == "sfOF" || asf == "sfP" || asf == "sfKD" ||
      asf == "sfHSD" || asf == "user") {
      
      // stage 1
      double cumAlphaSpent;
      if (asf == "user") {
        cumAlphaSpent = userAlphaSpending[0];
      } else {
        cumAlphaSpent = errorSpent(t[0], alpha, asf, asfpar);
      }
      
      if (!efficacyStopping1[0]) {
        criticalValues1[0] = 6.0;
      } else {
        criticalValues1[0] = R::qnorm(1 - cumAlphaSpent, 0, 1, 1, 0);
      }
      
      
      // lambda expression for finding the critical Values at stage k
      int k=0;
      auto f = [&k, &cumAlphaSpent, &criticalValues1,
                theta, t](double aval)->double {
                  NumericVector u(k+1), l(k+1);
                  for (int i=0; i<k; i++) {
                    u[i] = criticalValues1[i];
                    l[i] = -6.0;
                  }
                  u[k] = aval;
                  l[k] = -6.0;
                  
                  IntegerVector idx = Range(0,k);
                  List probs = exitprob(u, l, theta[idx], t[idx]);
                  double cpu = sum(NumericVector(probs[0]));
                  return cpu - cumAlphaSpent;
                };
      
      // subsequent stages
      for (k=1; k<kMax; k++) {
        if (asf == "user") {
          cumAlphaSpent = userAlphaSpending[k];
        } else {
          cumAlphaSpent = errorSpent(t[k], alpha, asf, asfpar);
        }
        
        if (!efficacyStopping1[k]) {
          criticalValues1[k] = 6.0;
        } else {
          if (f(6) > 0) { // no alpha spent at current visit
            criticalValues1[k] = 6.0;
          } else {
            criticalValues1[k] = brent(f, 0, 6, 1e-6);
          }
        }
      }
    } else {
      stop("Invalid type of alpha spending");
    }
    
    
  }
  
  
  
  auto f = [beta, kMax, informationRates1, futilityStopping1,
            criticalValues1, &futilityBounds1, 
            bsf, bsfpar, userBetaSpending, 
            missingFutilityBounds](double aval)->double {
              
              if (kMax > 1) {
                if (missingFutilityBounds && bsf=="none") {
                  futilityBounds1 = rep(-6.0, kMax);
                  futilityBounds1[kMax-1] = criticalValues1[kMax-1];
                } else if (!missingFutilityBounds && 
                  futilityBounds1.size() == kMax-1) {
                  futilityBounds1.push_back(criticalValues1[kMax-1]);
                } else if (!missingFutilityBounds && 
                  futilityBounds1.size() < kMax-1) {
                  stop("Insufficient length of futilityBounds");
                }
              } else {
                if (missingFutilityBounds) {
                  futilityBounds1 = criticalValues1[kMax-1];
                }
              }
              
              NumericVector theta = rep(aval, kMax);
              NumericVector t = clone(informationRates1);
              
              // compute the stagewise exit probabilities for efficacy and 
              // futility
              
              if (!missingFutilityBounds || bsf=="none" || kMax==1) {
                List probs = exitprob(criticalValues1, futilityBounds1, 
                                      theta, t);
                NumericVector pu = NumericVector(probs[0]);
                double overallReject = sum(pu);
                return overallReject - (1-beta);
              } else {
                // initialize futilityBounds to be updated
                futilityBounds1 = NumericVector(kMax);
                double epsilon;
                
                // first stage
                int k = 0;
                double cumBetaSpent;
                if (bsf=="user") {
                  cumBetaSpent = userBetaSpending[0];
                } else {
                  cumBetaSpent = errorSpent(t[0], beta, bsf, bsfpar);
                }
                
                if (!futilityStopping1[0]) {
                  futilityBounds1[0] = -6.0;
                } else {
                  epsilon = R::pnorm(criticalValues1[0] -
                    theta[0]*sqrt(t[0]), 0, 1, 1, 0) - cumBetaSpent;
                  if (epsilon < 0) return -1.0;
                  futilityBounds1[0] = R::qnorm(cumBetaSpent, 0, 1, 1, 0) +
                    theta[0]*sqrt(t[0]);
                }
                
                
                // lambda expression for finding the futility bound at stage k
                auto g = [&k, &cumBetaSpent, criticalValues1, &futilityBounds1,
                          theta, t](double aval)->double {
                            NumericVector u(k+1);
                            NumericVector l(k+1);
                            for (int i=0; i<k; i++) {
                              u[i] = criticalValues1[i];
                              l[i] = futilityBounds1[i];
                            }
                            u[k] = 6.0;
                            l[k] = aval;
                            
                            IntegerVector idx = Range(0,k);
                            List probs = exitprob(u, l, theta[idx], t[idx]);
                            double cpl = sum(NumericVector(probs[1]));
                            return cpl - cumBetaSpent;
                          };
                
                
                for (k=1; k<kMax; k++) {
                  if (bsf == "user") {
                    cumBetaSpent = userBetaSpending[k];
                  } else {
                    cumBetaSpent = errorSpent(t[k], beta, bsf, bsfpar);
                  }
                  
                  if (!futilityStopping1[k]) {
                    futilityBounds1[k] = -6.0;
                  } else {
                    epsilon = g(criticalValues1[k]);
                    
                    if (g(-6.0) > 0) { // no beta spent at the current visit
                      futilityBounds1[k] = -6.0;
                    } else if (epsilon > 0) {
                      futilityBounds1[k] = brent(g, -6.0, criticalValues1[k], 
                                                 1e-6);
                    } else if (k < kMax-1) {
                      return -1.0;
                    }
                    
                  }
                }
                
                return epsilon;
                
              }
            };
  
  double drift = brent(f, 0, 6, 0.0001);
  double driftf = R::qnorm(1-alpha, 0, 1, 1, 0) + R::qnorm(1-beta, 0, 1, 1, 0);
  double inflationFactor = pow(drift/driftf, 2);
  
  futilityBounds1[kMax-1] = criticalValues1[kMax-1];
  
  NumericVector theta = rep(drift, kMax);
  NumericVector t = clone(informationRates1);
  List probs = exitprob(criticalValues1, futilityBounds1, theta, t);
  
  // output the results
  NumericVector efficacyP(kMax);
  NumericVector futilityP(kMax);
  for (int i=0; i<kMax; i++) {
    efficacyP[i] = 1 - R::pnorm(criticalValues1[i], 0, 1, 1, 0);
    futilityP[i] = 1 - R::pnorm(futilityBounds1[i], 0, 1, 1, 0);
  }
  
  
  // stagewise total exit probabilities
  NumericVector pu(kMax), pl(kMax), ptotal(kMax);
  pu = NumericVector(probs[0]);
  pl = NumericVector(probs[1]);
  ptotal = pu + pl;
  
  double overallReject = sum(pu);
  NumericVector cpu = cumsum(pu);
  NumericVector cpl = cumsum(pl);
  
  NumericVector futilityBounds0 = rep(-6.0, kMax);
  futilityBounds0[kMax-1] = criticalValues1[kMax-1];
  NumericVector theta0(kMax);
  List probs0 = exitprob(criticalValues1, futilityBounds0, theta0, t);
  NumericVector cumAlphaSpent = cumsum(NumericVector(probs0[0]));
  
  IntegerVector stageNumber = seq_len(kMax);
  
  for (int i=0; i<kMax; i++) {
    if (criticalValues1[i] == 6) {
      efficacyStopping1[i] = 0;
    }
    
    if (futilityBounds1[i] == -6) {
      futilityStopping1[i] = 0;
    }
  }
  
  
  DataFrame byStageResults = DataFrame::create(
    _["informationRates"] = informationRates1,
    _["efficacyBounds"] = criticalValues1,
    _["futilityBounds"] = futilityBounds1,
    _["rejectPerStage"] = pu,
    _["futilityPerStage"] = pl,
    _["cumulativeRejection"] = cpu,
    _["cumulativeFutility"] = cpl,
    _["cumulativeAlphaSpent"] = cumAlphaSpent,
    _["efficacyP"] = efficacyP,
    _["futilityP"] = futilityP,
    _["efficacyStopping"] = efficacyStopping1,
    _["futilityStopping"] = futilityStopping1);
  
  DataFrame overallResults = DataFrame::create(
    _["overallReject"] = overallReject,
    _["alpha"] = (cumAlphaSpent[kMax-1]),
    _["kMax"] = kMax,
    _["drift"] = drift,
    _["inflationFactor"] = inflationFactor);
  
  List settings = List::create(
    _["typeAlphaSpending"] = typeAlphaSpending,
    _["parameterAlphaSpending"] = parameterAlphaSpending,
    _["userAlphaSpending"] = userAlphaSpending,
    _["typeBetaSpending"] = typeBetaSpending,
    _["parameterBetaSpending"] = parameterBetaSpending,
    _["userBetaSpending"] = userBetaSpending);
  
  List result = List::create(
    _["byStageResults"] = byStageResults,
    _["overallResults"] = overallResults,
    _["settings"] = settings);
  
  result.attr("class") = "design";
  
  return result;
  
}

