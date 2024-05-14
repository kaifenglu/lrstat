#include <Rcpp.h>
#include "utilities.h"

using namespace Rcpp;

//' @title Analysis of Simon's Bayesian basket trials
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
List simonBayesAnalysis(
    const int nstrata = NA_INTEGER,
    const IntegerVector& r = NA_INTEGER,
    const IntegerVector& n = NA_INTEGER,
    const double lambda = NA_REAL,
    const double gamma = NA_REAL,
    const double phi = NA_REAL,
    const double plo = NA_REAL) {

  if (nstrata == NA_INTEGER) {
    stop("nstrata must be provided");
  }

  if (nstrata <= 0) {
    stop("nstrata must be a positive integer");
  }

  if (is_true(any(is_na(r)))) {
    stop("r must be provided");
  }

  if (is_true(any(is_na(n)))) {
    stop("n must be provided");
  }

  if (r.size() != nstrata) {
    stop("Invalid length for r");
  }

  if (n.size() != nstrata) {
    stop("Invalid length for n");
  }

  if (is_true(any(r < 0))) {
    stop("r must be nonnegative");
  }

  if (is_true(any(r > n))) {
    stop("r must be less than or equal to n");
  }

  if (R_isnancpp(lambda)) {
    stop("lambda must be provided");
  }

  if (R_isnancpp(gamma)) {
    stop("gamma must be provided");
  }

  if (R_isnancpp(phi)) {
    stop("phi must be provided");
  }

  if (R_isnancpp(plo)) {
    stop("plo must be provided");
  }

  if (lambda <= 0 || lambda >= 1) {
    stop("lambda must lie between 0 and 1");
  }

  if (gamma <= 0 || gamma >= 1) {
    stop("gamma must lie between 0 and 1");
  }

  if (phi <= 0 || phi >= 1) {
    stop("phi must lie between 0 and 1");
  }

  if (plo <= 0 || plo >= 1) {
    stop("plo must lie between 0 and 1");
  }

  if (plo >= phi) {
    stop("plo must be less than phi");
  }

  NumericVector r1 = NumericVector(r);
  NumericVector n1 = NumericVector(n);

  int ncases = (int)std::pow(2, nstrata);
  NumericMatrix incid(ncases, nstrata);
  NumericVector prior(ncases), like(ncases);
  for (int i=0; i<ncases; i++) {
    int number = ncases - i - 1;
    NumericVector cc(nstrata);
    for (int j=0; j<nstrata; j++) {
      cc[j] = (number/(int)std::pow(2, nstrata-1-j)) % 2;
    }

    if (is_true(all(cc==1)) || is_true(all(cc==0))) {
      prior[i] = lambda*(gamma*(cc[0]==1) + (1-gamma)*(cc[0]==0)) +
        (1-lambda)*(pow(gamma,nstrata)*(cc[0]==1) +
        pow(1-gamma,nstrata)*(cc[0]==0));
    } else {
      double y = sum(cc);
      prior[i] = (1-lambda)*pow(gamma,y)*pow(1-gamma,nstrata-y);
    }

    NumericVector x = phi*cc + plo*(1-cc);
    like[i] = exp(sum(r1*log(x) + (n1-r1)*log(1-x)));
    incid(i,_) = cc;
  }
  NumericVector post = prior*like;
  post = post/sum(post);

  NumericVector q_prior(nstrata), q_post(nstrata);
  for (int j=0; j<nstrata; j++) {
    q_prior[j] = sum(incid(_,j)*prior);
    q_post[j] = sum(incid(_,j)*post);
  }

  return List::create(
    _["case"] = incid,
    _["prior_case"] = prior,
    _["prior_stratum"] = q_prior,
    _["post_case"] = post,
    _["post_stratum"] = q_post);
}



//' @title Simulation of Simon's Bayesian basket trials
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
List simonBayesSim(
    const NumericVector& p = NA_REAL,
    const NumericVector& accrualTime = 0,
    const NumericVector& accrualIntensity = NA_REAL,
    const NumericVector& stratumFraction = 1,
    const double lambda = NA_REAL,
    const double gamma = NA_REAL,
    const double phi = NA_REAL,
    const double plo = NA_REAL,
    const double T = NA_REAL,
    const int maxSubjects = NA_INTEGER,
    const IntegerVector& plannedSubjects = NA_INTEGER,
    const int maxNumberOfIterations = 1000,
    const int maxNumberOfRawDatasets = 1,
    const int seed = NA_INTEGER) {

  int nstrata = stratumFraction.size();

  if (is_true(any(is_na(p)))) {
    stop("p must be provided");
  }

  if (is_true(any((p <= 0) | (p >= 1)))) {
    stop("p must lie between 0 and 1");
  }

  if (accrualTime[0] != 0) {
    stop("accrualTime must start with 0");
  }

  if (accrualTime.size() > 1 && is_true(any(diff(accrualTime) <= 0))) {
    stop("accrualTime should be increasing");
  }

  if (is_true(any(is_na(accrualIntensity)))) {
    stop("accrualIntensity must be provided");
  }

  if (accrualTime.size() != accrualIntensity.size()) {
    stop("accrualTime must have the same length as accrualIntensity");
  }

  if (is_true(any(accrualIntensity < 0))) {
    stop("accrualIntensity must be non-negative");
  }

  if (is_true(any(stratumFraction <= 0))) {
    stop("stratumFraction must be positive");
  }

  if (fabs(sum(stratumFraction) - 1.0) > 1.0e-8) {
    stop("stratumFraction must sum to 1");
  }

  if (R_isnancpp(lambda)) {
    stop("lambda must be provided");
  }

  if (R_isnancpp(gamma)) {
    stop("gamma must be provided");
  }

  if (R_isnancpp(phi)) {
    stop("phi must be provided");
  }

  if (R_isnancpp(plo)) {
    stop("plo must be provided");
  }

  if (R_isnancpp(T)) {
    stop("T must be provided");
  }

  if (lambda <= 0 || lambda >= 1) {
    stop("lambda must lie between 0 and 1");
  }

  if (gamma <= 0 || gamma >= 1) {
    stop("gamma must lie between 0 and 1");
  }

  if (phi <= 0 || phi >= 1) {
    stop("phi must lie between 0 and 1");
  }

  if (plo <= 0 || plo >= 1) {
    stop("plo must lie between 0 and 1");
  }

  if (plo >= phi) {
    stop("plo must be less than phi");
  }

  if (T <= 0 || T >= 1) {
    stop("T must lie between 0 and 1");
  }

  if (maxSubjects == NA_INTEGER) {
    stop("maxSubjects must be provided");
  }

  if (maxSubjects < 1) {
    stop("maxSubjects must be a positive integer");
  }

  if (is_true(any(is_na(plannedSubjects)))) {
    stop("plannedSubjects must be provided");
  }

  if (plannedSubjects[0] <= 0) {
    stop("Elements of plannedSubjects must be positive");
  }

  if (plannedSubjects.size() > 1 &&
      is_true(any(diff(plannedSubjects) <= 0))) {
    stop("Elements of plannedSubjects must be increasing");
  }

  if (maxNumberOfIterations < 1) {
    stop("maxNumberOfIterations must be a positive integer");
  }

  if (maxNumberOfRawDatasets < 0) {
    stop("maxNumberOfRawDatasets must be a non-negative integer");
  }


  int iter, i, j, k, l, stage;
  double enrollt, u;
  int kMax = plannedSubjects.size();
  int nreps = maxNumberOfIterations;

  LogicalVector act = (p == phi);
  int nactive = sum(act);

  NumericVector cumStratumFraction = cumsum(stratumFraction);
  NumericVector post_stratum(nstrata);
  IntegerVector n(nstrata), r(nstrata);
  LogicalVector open(nstrata), pos(nstrata), neg(nstrata);
  List bayes;

  NumericVector arrivalTime(maxSubjects);
  IntegerVector stratum(maxSubjects), y(maxSubjects);
  IntegerVector iterationNumber(nreps), N(nreps), nact(nreps), nopen(nreps);
  IntegerVector tpos(nreps), fneg(nreps), fpos(nreps), tneg(nreps);
  IntegerVector numberOfStrata(nreps, nstrata);

  // cache for the patient-level raw data to extract
  int nrow1 = kMax*maxNumberOfRawDatasets*maxSubjects;
  IntegerVector iterationNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector stageNumberx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector subjectIdx = IntegerVector(nrow1, NA_INTEGER);
  NumericVector arrivalTimex = NumericVector(nrow1, NA_REAL);
  IntegerVector stratumx = IntegerVector(nrow1, NA_INTEGER);
  IntegerVector yx = IntegerVector(nrow1, NA_REAL);

  // cache for the summary data to extract
  int nrow2 = kMax*maxNumberOfIterations*nstrata;
  IntegerVector iterationNumbery = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector stageNumbery = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector stratumy = IntegerVector(nrow2, NA_INTEGER);
  LogicalVector activey = LogicalVector(nrow2, NA_LOGICAL);
  IntegerVector ny = IntegerVector(nrow2, NA_INTEGER);
  IntegerVector ry = IntegerVector(nrow2, NA_INTEGER);
  NumericVector posty = NumericVector(nrow2, NA_REAL);
  LogicalVector openy = LogicalVector(nrow2, NA_LOGICAL);
  LogicalVector posy = LogicalVector(nrow2, NA_LOGICAL);
  LogicalVector negy = LogicalVector(nrow2, NA_LOGICAL);

  int index1 = 0, index2 = 0;

  // set up random seed
  if (seed != NA_INTEGER) {
    set_seed(seed);
  }

  for (iter=0; iter<nreps; iter++) {
    // initialize the contents in each stratum
    n.fill(0);
    r.fill(0);
    open.fill(1);
    pos.fill(0);
    neg.fill(0);

    k = 0;      // index of the number of subjects included in analysis
    stage = 0;
    enrollt = 0;
    for (i=0; i<100000; i++) {
      // generate accrual time
      u = R::runif(0,1);
      enrollt = qtpwexpcpp1(u, accrualTime, accrualIntensity, enrollt, 1, 0);

      // generate stratum information
      u = R::runif(0,1);
      for (j=0; j<nstrata; j++) {
        if (cumStratumFraction[j] > u) {
          break;
        }
      }

      // if the stratum is open, generate the response for the subject
      if (open[j]) {
        arrivalTime[k] = enrollt;
        stratum[k] = j+1;
        u = R::runif(0,1);
        y[k] = u < p[j] ? 1 : 0;

        // update the number of subjects and responders in the stratum
        n[j] = n[j] + 1;
        r[j] = r[j] + y[k];

        k++;

        // interim analysis
        if (is_true(any(plannedSubjects == k))) {

          // output raw data
          if (iter < maxNumberOfRawDatasets) {
            for (int idx=0; idx<k; idx++) {
              iterationNumberx[index1] = iter+1;
              stageNumberx[index1] = stage+1;
              subjectIdx[index1] = idx+1;
              arrivalTimex[index1] = arrivalTime[idx];
              stratumx[index1] = stratum[idx];
              yx[index1] = y[idx];

              index1++;
            }
          }

          // calculate the posterior probabilities
          bayes = simonBayesAnalysis(nstrata, r, n, lambda, gamma, phi, plo);
          post_stratum = bayes["post_stratum"];

          // whether to close the stratum due to positive or negative results
          for (l=0; l<nstrata; l++) {
            if (open[l]) {
              if (post_stratum[l] > T) {
                pos[l] = 1;
                open[l] = 0;
              } else if (post_stratum[l] < 1 - T) {
                neg[l] = 1;
                open[l] = 0;
              }
            }

            // output strata-level data
            iterationNumbery[index2] = iter+1;
            stageNumbery[index2] = stage+1;
            stratumy[index2] = l+1;
            activey[index2] = act[l];
            ny[index2] = n[l];
            ry[index2] = r[l];
            posty[index2] = post_stratum[l];
            openy[index2] = open[l];
            posy[index2] = pos[l];
            negy[index2] = neg[l];

            index2++;
          }

          stage++;
        }


        // stop the trial if all strata are closed or max subjects reached
        if (is_true(all(open == 0)) || (k == maxSubjects)) {
          iterationNumber[iter] = iter+1;
          N[iter] = k;
          nact[iter] = nactive;
          tpos[iter] = sum(pos & act);
          fneg[iter] = sum(neg & act);
          fpos[iter] = sum(pos & !act);
          tneg[iter] = sum(neg & !act);
          nopen[iter] = sum(open);

          break;
        }
      }
    }
  }


  DataFrame rawdata;
  if (maxNumberOfRawDatasets > 0) {
    LogicalVector sub1 = !is_na(iterationNumberx);
    iterationNumberx = iterationNumberx[sub1];
    stageNumberx = stageNumberx[sub1];
    subjectIdx = subjectIdx[sub1];
    arrivalTimex = arrivalTimex[sub1];
    stratumx = stratumx[sub1];
    yx = yx[sub1];

    rawdata = DataFrame::create(
      _["iterationNumber"] = iterationNumberx,
      _["stageNumber"] = stageNumberx,
      _["subjectId"] = subjectIdx,
      _["arrivalTime"] = arrivalTimex,
      _["stratum"] = stratumx,
      _["y"] = yx);
  }


  // simulation summary data set
  LogicalVector sub2 = !is_na(iterationNumbery);
  iterationNumbery = iterationNumbery[sub2];
  stageNumbery = stageNumbery[sub2];
  stratumy = stratumy[sub2];
  activey = activey[sub2];
  ny = ny[sub2];
  ry = ry[sub2];
  posty = posty[sub2];
  openy = openy[sub2];
  posy = posy[sub2];
  negy = negy[sub2];

  DataFrame sumdata1 = DataFrame::create(
    _["iterationNumber"] = iterationNumbery,
    _["stageNumber"] = stageNumbery,
    _["stratum"] = stratumy,
    _["active"] = activey,
    _["n"] = ny,
    _["r"] = ry,
    _["posterior"] = posty,
    _["open"] = openy,
    _["positive"] = posy,
    _["negative"] = negy);


  DataFrame sumdata2 = DataFrame::create(
    _["iterationNumber"] = iterationNumber,
    _["numberOfStrata"] = numberOfStrata,
    _["n_active_strata"] = nact,
    _["true_positive"] = tpos,
    _["false_negative"] = fneg,
    _["false_positive"] = fpos,
    _["true_negative"] = tneg,
    _["n_indet_strata"] = nopen,
    _["numberOfSubjects"] = N);


  double mn_nact = mean(NumericVector(nact));
  double mn_tpos = mean(NumericVector(tpos));
  double mn_fneg = mean(NumericVector(fneg));
  double mn_fpos = mean(NumericVector(fpos));
  double mn_tneg = mean(NumericVector(tneg));
  double mn_inde = mean(NumericVector(nopen));
  double mn_N = mean(NumericVector(N));

  DataFrame overview = DataFrame::create(
    _["numberOfStrata"] = nstrata,
    _["n_active_strata"] = mn_nact,
    _["true_positive"] = mn_tpos,
    _["false_negative"] = mn_fneg,
    _["false_positive"] = mn_fpos,
    _["true_negative"] = mn_tneg,
    _["n_indet_strata"] = mn_inde,
    _["numberOfSubjects"] = mn_N);

  List result;

  if (maxNumberOfRawDatasets > 0) {
    result = List::create(
      _["rawdata"] = rawdata,
      _["sumdata1"] = sumdata1,
      _["sumdata2"] = sumdata2,
      _["overview"] = overview);
  } else {
    result = List::create(
      _["sumdata1"] = sumdata1,
      _["sumdata2"] = sumdata2,
      _["overview"] = overview);
  }

  return result;
}


