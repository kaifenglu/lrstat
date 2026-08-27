# lrstat 0.3.4

* updated `lrsim_bmTrtSel` TSSSD rank-method calibration and boundary
  calculations: added the known-correlation real rank-based method
  (`tsssd.k.rank`) and its conditional-error counterpart. The existing
  effective-multiplicity approximations are named `tsssd.k.modrank` and
  `tsssd.uk.modrank`, with corresponding `*.modrank.ce` methods.
* optimized selection-boundary calculations by making use of no efficacy
  stopping in stage 1.
* renamed the `fPCStagewise` input parameters `stg1_inthyp_nr` and
  `stg2_elemhyp` to `stg1_inthyp_nr_idx` and `stg2_elemhyp_idx`, respectively,
  to make their index semantics explicit and align them with the output names
  and related adaptive-testing interfaces
* standardized `lrsim_bmTrtSel` testing method identifiers to lowercase,
  including inputs, returned `methods` and `byMethod` names, and
  `sumdataTTE` rejection indicator columns
* added `sumdataBIN` and `sumdataTTE` summary datasets to `lrsim_bmTrtSel`,
  along with optional subject-level `rawdataBIN` and `rawdataTTE` datasets
  controlled by `maxNumberOfRawDatasets`

# lrstat 0.3.3

* standardized C++ pointer/reference declaration style across src with a
  project-level `.clang-format` (80-column, 2-space indent, `type *name` and
  `type &name` alignment) and aligned headers/sources to that convention
* excluded `.clang-format` from package build artifacts via `.Rbuildignore`

* corrected conditional Type I error calculations in `adaptDesign`,
  `adaptDesign_multiarm`, `adaptDesign_seamless`, `getCP`, `getCP_multiarm`, and
  `getCP_seamless` to use normalized post-adaptation information rates
  consistently, including when the Muller and Schafer method is selected
* added `finthyp` to construct the intersection-hypothesis indicator matrix for
  a given number of elementary hypotheses, matching the row order used by
  `fwgtmat` and `fDefaultWgtmat`, for use when migrating custom weight matrices
  to the new `fadjpbon`/`fadjpsim`/`fadjpdun` interfaces
* reduced per-iteration buffer allocation in `fadjpsim`'s Simes p-value
  adjustment by reusing preallocated vectors instead of allocating new ones for
  each family block
* avoided recomputing the transform gradient for each probability in
  `survQuantile`, computing it once per event time instead
* simplified graphical and adaptive multiple-testing p-value adjustment
  interfaces by making p-values the first argument and automatically supplying
  equal intersection weights, a single family, and default within-family
  correlations when these inputs are omitted
* added public R wrappers with `nthreads` controls for the multi-arm and
  seamless design, conditional power, and confidence interval functions
* simplified one-active-arm `getBound_multiarm` and `getBound_seamless`
  calculations by delegating to `getBound`
* added rank-adjusted and conditional-error-updated TSSSD methods to
  `lrsim_bmTrtSel` (`TSSSD.k`, `TSSSD.uk`, `TSSSD.k.rank`, `TSSSD.uk.rank`,
  `TSSSD.k.ce`, `TSSSD.uk.ce`, `TSSSD.k.rank.ce`, `TSSSD.uk.rank.ce`) and
  improved TSSSD boundary computation efficiency, including cached single-arm
  boundaries and precomputed CE nominal boundaries
* replaced the lpSolve dependency in getDesignAgreement with a native C++17
  simplex solver exposed through Rcpp, preserving the documented design results
  and adding regression tests for the agreement constraints
* added lrsim_bmTrtSel and print.lrsim_bmTrtSel for simulation and summary
  output of phase 2/3 seamless design selecting one dose at the end of phase 2
  based on the posterior benefit-risk tradeoff of correlated binary efficacy and
  toxicity endpoints, and comparing the closed testing procedures with inverse
  normal combination, the conditional error rate method, the two-stage seamless
  design boundaries, and the unadjusted tests for the confirmatory analysis of
  the long-term endpoint
* added unit tests for the adaptive two-stage multiple testing functions
* allowed `alpha1 = 0` in `fCERStageBound` when no stage 1 efficacy stopping is
  planned
* set the stage 2 boundary to infinity in `fCERNewBound` when the conditional
  error rate exceeds 1
* used struct output to reduce overhead of returning multiple values for
  exitprobcpp, exitprob_seamless_cpp and exitprob_multiarm_cpp
* allowed the retesting of an active hypothesis if its weight has increased due
  to the rejection of another hypothesis in fseqbon
* added getDesign2 for futility boundary under the null hypothesis for group
  sequential designs
* renamed mams to multiarm for more informative description
* added rdsim_seamless for simulation of phase 2/3 seamless design for risk
  difference
* updated lrsim_seamless and rdsim_seamless to use rank-based phase-2 selection
  naming and outputs, including rankp0 in the wrapper inputs, overview output,
  and print methods
* updated seamless_design.cpp related functions to use rank-based treatment
  selection with rankp0 throughout seamless design, conditional power,
  confidence interval, and print outputs
* added rdsim_multiarm for simulation of multi-arm multi-stage design for risk
  difference
* added lrsim_mcpmod and print.lrsim_mcpmod for simulation and summary output of
  MCPMod design using the log-rank test
* updated mvnormr to support positive semidefinite covariance matrix via minimal
  diagonal stabilization during factorization
* used seed = 314159 as the default for pmvnormr and qmvnormr to ensure the same
  random numbers are generated and results are reproducible across different
  runs
* added the selection to Stage 2 probability in lrsim_seamless
* fixed liferegr to handle empty optional variable inputs (stratum, covariates,
  weight, offset, id) without triggering a formula parse error from an empty
  right-hand side
* updated fadjpbon, fadjpsim, and fadjpdun to add the local p-values and the
  incidence matrix for the intersection hypotheses
* updated fadjpdun for calculating the local p-values for the intersection
  hypotheses
* updated ftrunc, fstdmix, and fmodmix to return a list with padj and pinter
  components instead of a bare matrix/vector of adjusted p-values

# lrstat 0.3.2

* removed the pivoting step for pmvnormr and qmvnormr to take advantage of
  independent increment structure of score process
* removed the fast parameter for pmvnormr and qmvnormr so that the default is to
  use fast approximation
* updated the point estimate for getRCI
* renamed tsssd to seamless in the function names for two-stage seamless
  sequential design for treatment/dose selection
* updated the critical values when typeAlphaSpending = "none" in
  getDesign_seamless
* added integratedTrial to adaptDesign
* updated the backward image confidence interval method to account for the case
  when astar > conditional_power
* updated point estimate for repeated confidence interval
* added functions for multi-arm multi-stage (MAMS) design and analysis,
  including exitprob_multiarm, getBound_multiarm, getDesign_multiarm,
  adaptDesign_multiarm, getCI_multiarm, getADCI_multiarm
* ensured that if typeAlphaSpending is "OF", "P", "WT", or "none", then
  informationRates, efficacyStopping, and spendingTime must be of full length
  kMax, and informationRates and spendingTime must end with 1 for getBound,
  getBound_seamless, and getBound_multiarm
* ensured that if typeAlphaSpending is "OF", "P", "WT", or "none", then
  informationRates, efficacyStopping, and spendingTime must be of full length
  kMax, and informationRates and spendingTime must end with 1 for getCI, getRCI
* ensured that if typeAlphaSpendingNew is "OF", "P", "WT", or "none", then
  informationRatesNew, efficacyStoppingNew, and spendingTimeNew must be of full
  length kNew, and informationRatesNew and spendingTimeNew must end with 1 for
  getADCI, getADRCI
* moved up the position of parameter MullerSchafer in getADCI, getADRCI
* added getCP_seamless for conditional power calculation for two-stage seamless
  sequential design for treatment/dose selection
* added getCP_multiarm for conditional power calculation for multi-arm
  multi-stage design
* added adaptDesign_seamless for power and sample size calculation of adaptive
  two-stage seamless sequential design for treatment/dose selection
* added getCI_seamless for confidence interval calculation for two-stage
  seamless sequential design for treatment/dose selection
* added getADCI_seamless for confidence interval calculation using the backward
  image method for adaptive two-stage seamless sequential design for
  treatment/dose selection
* moved the f_pvalue function from confidence_interval.cpp to generic_design.cpp
* updated infoRatesNew when informationRatesNew is not missing in
  adaptDesign_seamless_cpp and adaptDesign_multiarm_cpp
* added default value of 1 for parameter r of adaptDesign_seamless
* updated pmvnormmccpp to use s and T instead of Ivec
* added byLevelBounds in the output of getDesign_multiarm
* updated the condition for calculating cpu0 in adaptDesign_multiarm_cpp
* updated the description for adaptDesign_multiarm output
* updated the default value of r and rNew to 1 in adaptDesign_multiarm and
  getADCI_multiarm
* updated the stage‑1 covariance used for bound calculation in
  getBound_seamless to account for corr_known when using alpha‑spending
  functions
* used memoization to avoid recalculating p0 and p1 in the lambda h of
  preject_by_arm (exitprob_seamless)
* added stopStage to output data sets and added bestArm and reject to sumdata2
  output data set for lrsim_seamless
* used NULL as the default value for lambdas and gammas for lrsim_seamless
* added lrsim_multiarm for simulation of multi-arm multi-stage design
* added futilityBounds, futilityCP, futilityTheta, futilityBoundsInt,
  futilityCPInt, futilityThetaInt to adaptDesign, adaptDesign_seamless, and
  adaptDesign_multiarm
* removed typeBetaSpending and parameterBetaSpending from adaptDesign,
  adaptDesign_seamless, and adaptDesign_multiarm

# lrstat 0.3.1

* calculated eventsPerStage, dropoutPerStage, subjectsPerStage, and
  analysisTimePerstage in the lrsim output based on trials having the stage
  (i.e., did not stop before the stage).
* added mvnormr for multivariate normal probability and quantile calculations
* added tsssd_design: analytical calculation for two-stage seamless sequential
  designs for treatment/dose selection.
* added lrstim_tsssd: simulation function for two-stage seamless sequential
  designs for treatment/dose selection.
* printed power with 4 instead of 3 decimal places in the output
* used std::size_t instead of size_t in all header files
* added pnorm_fast, qnorm_acklam, expand_stratified_to_slice, and invchol
  utility functions

# lrstat 0.3.0

* replaced accrualDuration with n for lrsim, lrsim2e, lrsim3a, lrsim2e3a,
  lrsimsub, and binary\_tte\_sim to allow the simulation to run with a fixed
  number of subjects
* added time2 and weight parameters to kmest
* added time2, weight, and weight\_readj parameters to lrtest
* added assess\_phregr and zph\_phregr functions
* updated findInterval3 to add rightmost\_closed, all\_inside, and left\_open
  parameters
* added dtpwexp for density function of truncated piecewise exponential
  distribution
* updated hazard\_pd to match probabilities at quantiles in addition to cut
  points
* added corr\_pfs\_os to obtain correlation between PFS and OS given correlation
  between PD and OS
* added hazard\_pd2 to obtain the hazards for PD to match distribution for PFS
  and correlation between PFS and OS
* renamed rho to rho\_pd\_os for added clarity for lrsim2e and lrsim2e3a
* updated hazard\_sub to match probabilities at quantiles in addition to cut
  points
* moved logparams to logistic\_regression.cpp and removed logparams, f\_der\_0,
  f\_ressco\_0, logisregloop, and logisregplloop from logistic\_regression.h
* updated lrsim2e and lrsim2e3a to account for additional quantile-based cut
  points for pd
* updated lrsimsub to account for additional quantile-based cut points for
  biomarker negative sub population
* moved aftparams and coxparams to survival\_analysis.cpp and removed aftparams,
  f\_der\_1, f\_ressco\_1, f\_jj\_1, liferegloop, liferegplloop, f\_ld\_1,
  coxparams, f\_der\_2, phregloop, phregplloop, f\_basehaz, f\_ressco\_2,
  f\_resmart, f\_ressch, and f\_der\_i\_2 from survival\_analysis.h
* added ptpwexpcpp1 and updated ptpwexpcpp accordingly
* added mtpwexp for the mean and variance of a truncated piecewise exponential
  distribution
* added quad2d for numerical integration in two dimensions
* removed hd, pd, ad, house, row\_house, col\_house, givens, row\_rot, col\_rot,
  house\_bidiag, zero\_diagonal, and svd\_step from utilities.h
* renamed the parameters in lrsim2e and lrsim2e3a, replacing e1 by pfs and e2 by
  os
* ran the simulations in serial mode to pass CRAN checks

# lrstat 0.2.15

* added na.action = na.pass for model frame construction involving covariates
  for all methods and updated lrstat-package.R to add na.omit and na.pass to
  reflect this change
* added pbvnorm for the distribution function for standard bivariate normal
* added hazard\_pd for the hazard function for progressive disease given the
  piecewise hazard functions for PFS and OS and the correlation between PD and
  OS
* added hazard\_sub for the hazard function for the biomarker negative
  subpopulation given the hazard function of the ITT population, the hazard
  function and the prevalence of the biomarker positive subpopulation
* added the init parameter and the fail flag to the output of logisregr,
  liferegr, and phregr
* added the nullVariance parameter to getDesignOneProportion
* added details section to the documentation of getDesignPairedPropMcNemar,
  getDesignRiskDiff, getDesignRiskRatio, getDesignRiskRatioFM,
  getDesignOddsRatio, getDesignOneMultinom, getDesignTwoMultinom,
  getDesignTwoOrdinal, getDesignOrderedBinom, getDesignUnorderedBinom,
  getDesignUnorderedMultinom, pwexploglik, pwexpcuts, and lrschoenfeld
* added references to the documentation of ClopperPearsonCI, BOINTable,
  mTPI2Table
* added the tol parameter to pwexpcuts
* added the lrsimsub function for log-rank test simulation for enrichment design
* updated lrsim2e and lrsim2e3a to use the hazard\_pd function when generating
  PFS time
* updated kmstat.cpp and rmstat.cpp to report an error message if the
  information at the first interim is less than the information at the milestone
  time
* updated lrstat.cpp to avoid the redundant calculation of the final analysis
  time in the lrpower, lrsamplesize, lrpowerequiv, and lrsamplesizeequiv
  functions
* updated getDesignMeanDiffCarryover and getDesignMeanDiffCarryoverEquiv to fix
  the calculation of p and v
* updated getDesignMeanDiff, getDesignMeanDiffXO, and getDesignSlopeDiff to add
  the parameter value at the null hypothesis to the boundaries for the parameter
  of interest for noninferiority trials
* updated getDesignPairedMeanDiffEquiv, getDesignMeanDiffEquiv,
  getDesignMeanDiffXOEquiv, and getDesignSlopeDiffMMRM to simplify the
  calculation of attainedAlpha for equivalence trials
* updated getDesignWilcoxon, getDesignMeanDiffMMRM, getDesignOneSlope,
  getDesignSlopeDiff to add details to the documentation
* updated getDesignMeanDiffMMRM to add numberOfCompleters to the output and a
  reference to the documentation
* updated getDesignMeanDiffCarryover and getDesignMeanDiffCarryoverEquiv to add
  half\_width and nu to the output
* updated getDesignTwoWayANOVA to ensure the sample sizes are multiples of 4
* updated getDesignSlopeDiffMMRM to add more detailed description for the
  parameters w and N
* updated getDesignSlopeDiffMMRM to add the fixedFollowup parameter to allow
  fixed follow-up designs and modify the function for computing information and
  add number of completers accordingly
* updated getDesignSlopeDiffMMRM to add more content to the details section of
  the documentation
* updated the definition of the Z statistic in the details section of the
  documentation for getDesignEquiv
* updated print.designMeanDiffMMRM to add numberOfCompleters and analysisTime to
  the output
* updated print.designMeanDiffCarryover and print.designMeanDiffCarryoverEquiv
  to add half\_width to the output
* updated print.designSlopeDiffMMRM to add numberOfCompleters to the output for
  fixed follow-up designs
* updated the documentation of the output of lrschoenfeld
* replaced the survreg initial value method with the OLS initial value method
  for liferegr

# lrstat 0.2.14

* updated survival\_analysis to ignore intervals not at risk within each stratum
  without creating non overlapping times across strata
* updated documentation for the survsplit utility function
* removed bc from logisregr
* added residuals\_liferegr for residuals from parametric regression models for
  failure time data
* added maxiter and eps to logisregr, liferegr, and phregr
* added initial values for liferegcpp
* moved survQuantile from getDesignSurvival to survival\_analysis.cpp
* moved the detailed implementation from residuals\_phregr to
  residuals\_phregcpp
* added detail description of the spending functions for errorSpent

# lrstat 0.2.13

* checked whether the rounded up n is different from the original n0 before
  updating accrual intensity or accrual duration
* changed accrualIntensity to accrualIntensity1 in lrsamplesizeequiv when
  obtaining the timing of interim analysis
* updated kmpower1s and rmpower1s to remove drift and inflation factor from the
  output while adding number of events, number of dropouts, and number of
  milestone subjects to the output
* simplified the calculation of underlying survival probability for each
  treatment group in kmstat1
* updated kmpowersamplesize, kmsamplesize1s, kmsamplesizeequiv,
  rmpowersamplesize, rmsamplesize1s, rmsamplesizeequiv to remove the impossible
  case when followupTime is missing for fixed follow-up
* updated kmpowersamplesize, kmsamplesize1s, kmsamplesizeequiv,
  rmpowersamplesize, rmsamplesize1s, rmsamplesizeequiv to add a check for
  solving accrualIntensity when study duration is shorter than or equal to
  milestone
* updated kmpowersamplesize, kmsamplesize1s, kmsamplesizeequiv,
  rmpowersamplesize, rmsamplesize1s, rmsamplesizeequiv to add a case to decrease
  accrual duration for fixed follow-up under H0
* updated getDesignMeanDiffCarryover to allow specification of the treatment
  comparison of interest and whether to account for carryover effects in power
  calculation
* added getDesignMeanDiffCarryoverEquiv to test for equivalence in mean
  difference for the direct treatment effect of interest in crossover designs
* added the keep\_censor parameter to the kmest function to specify whether to
  retain censoring time points in the output data frame
* added interval to survsplit output
* used the same stopping criteria for SAS PROC LOGISTIC profile likelihood
* added the float\_to\_fraction utility function
* renamed runShinyApp to runShinyApp\_lrstat

# lrstat 0.2.12

* allowed time to take zero values in kmest, lrtest, rmest, and rmdiff
* added the bc parameter for bias correction for noncanonical parameterization
  and weighted logistic regression
* added loop to look for appropriate bracketing interval for the brent algorithm
  in lrstat, kmstat, rmstat, and nbstat
* updated f\_info nbstat to calculate the information directly
* updated the mean exposure calculations in nbstat for starting values for reml
* added svdcpp for singular value decomposition of a rectangular matrix
* added number of events, number of dropouts and number of milestone subjects to
  kmstat and kmpower
* added number of events, number of dropouts and number of milestone subjects to
  rmstat and rmpower
* removed variance ratios for equivalence trials assuming Wald test statistics
  with variance of parameter estimator evaluated under the alternative
  hypothesis
* removed the interval argument from getDurationFromNevents
* added lrschoenfeld for sample size calculation using the Schoenfeld formula
  under proportional hazards

# lrstat 0.2.11

* used numerical integration for lrstat for cases including unequal dropouts
  between treatment groups, hence removed the numSubintervals argument from
  lrstat1, lrstat, lrpower, lrsamplesize, lrpowerequiv, and lrsamplesizeequiv
* updated shiny app accordingly
* added noncanonical parametrizations for logisregr

# lrstat 0.2.10

* used newton-raphson instead of the vmmin algorithm for liferegr and phregr
* removed rpsft with all treatment switching functions moved to trtswitch
  package
* added the survfit\_phregr function to obtain the survivor curve estimates
* added the residuals\_phregr function to obtain the martingale, deviance,
  score, and schoenfeld residuals
* added the shilong dataset
* replaced std::ceil(x) with std::ceil(x - 1.0e-12) to handle superficial
  precision in binary presentation of a double number
* allowed rep and stratum parameters to represent more than one character
  variables in the source data
* added the logisregr function for logistic regression of a binary endpoint.

# lrstat 0.2.9

* fixed the runtime error: nan is outside the range of representable values of
  type 'int'
* added the binary\_tte\_sim function to simulate two endpoints with one being a
  binary endpoint and the other being a time-to-event endpoint
* added the rpsft function to estimate hazard ratio using rank-preserving
  structured failure time model

# lrstat 0.2.8

* updated std::sort handling for Rcpp character vectors
* used static cast to explicitly cast "long"" of the member function size() to
  "int" and cast "double" of the floor, ceil, round functions to "int"

# lrstat 0.2.7

* added the aml, heart, and tobin datasets
* in basket.cpp, replace y\[k] = R::rbinom(1, p\[j]) with u = R::runif(0,1);
  y\[k] = u < p\[j] ? 1 : 0;
* in misc.cpp, replace const in x = NA\_REAL with const int x = NA\_INTEGER for
  x = n, n1, y1, n2, y2
* added delete\[] for iwork and work pointers in the quad function
* added simonBayesAnalysis and simonBayesSim for Simon's Bayesian analysis of
  basket discovery trials
* updated the algorithm for finding the backward image (J, zJ) in getADCI
* updated the initial values for finding the confidence interval in getADCI
* added bmini for optimization and invsympd for inverting symmetric and positive
  defined matrices to utilities.cpp
* added phregr to estimate hazard ratios from the Cox model for right-censored
  or counting process data with robust variance
* added liferegr for parametric regression of failure time data with uncensored,
  right censored, left censored, and interval censored data

# lrstat 0.2.6

* replaced "R\_isnancpp(x)" with "x == NA\_INTEGER" for integer parameter x
* replaced the check for equality of spendingTime and informationRates in
  settings in prints.R

# lrstat 0.2.5

* updated survQuantile in getDesignMeans.R to handle new censoring and
  inestimable quantiles
* updated simon2stage to increase the search space of sample size
* added zstatRiskDiff, zstatRiskRatio, zstatOddsRatio, zstatRateDiff,
  zstatRateRatio for the score test of two-sample crude rates and
  exposure-adjusted incidence rates
* added mini for minimization of a univariate function over a finite interval
* updated powerRiskDiffExact, samplesizeRiskDiffExact, powerRiskRatioExact,
  samplesizeRiskRatioExact, powerRiskDiffExactEquiv,
  samplesizeRiskDiffExactEquiv, powerRiskRatioExactEquiv, and
  samplesizeRiskRatioExactEquiv to divide the interval of the response rate in
  the control group into 100 subintervals and obtain the minimum value within
  each subinterval in order to obtain the global minimum. This replaces the
  random draw of 500 values of the control group response rate for global
  minimum. In addition, the attained alpha calculation is updated to use the
  global maximum
* added riskDiffExactPValue, riskDiffExactCI, riskRatioExactPValue, and
  riskRatioExactCI for exact unconditional test and confidence limits based on
  the Miettinen \& Nurminen score statistics
* added rawdata to be used as an input data frame to test survival analysis
  functions
* added kmest and kmdiff for Kaplan-Meier estimates of survival curves and
  milestone survival difference
* added lrtest for log-rank test using the Fleming-Harrington family of weights
* added rmest and rmdiff for estimates of restricted mean survival times and
  restricted mean survival time difference
* added lower.tail and log.p arguments to ptpwexp and qtpwexp functions
* added the fquantile function to obtain the quantiles of a survival function
* added the pwexploglik function to obtain the profile log-likelihood function
  for the change points in the piecewise exponential approximation to a survival
  function
* added the pwexpcuts function to obtain a piecewise exponential distribution
  that approximates a survival distribution
* updated getDesignMeanDiffMMRM and getDesignSlopeDiffMMRM to add accrual and
  dropout information

# lrstat 0.2.4

* updated the check for pconfigs of getDesignLogistic to avoid an error on macOS
* replaced which\_min with which\_max in lrsim2e and lrsim2e3a
* replaced qromb with quad for numerical integration
* used numerical integration for lrstat
* renamed kmest to kmstat and use numerical integration
* added kmpower and kmsamplesize for power and sample size calculation for
  difference in milestone survival probability
* added kmpower1s and kmsamplesize1s for power and sample size calculation for
  one-sample milestone survival probability
* added rmst, covrmst, and rmstat for restricted mean survival time analysis
* added rmpower and rmsamplesize for power and sample size calculation for
  difference in restricted mean survival time
* added rmpower1s and rmsamplesize1s for power and sample size calculation for
  one-sample restricted mean survival time
* updated rmstat, rmpower, rmsamplesize, rmpower1s, and rmsamplesize1s to
  account for stratification
* added lrpowerequiv and lrsamplesizeequiv for power and sample size calculation
  for equivalence in hazard ratio
* added kmpowerequiv and kmsamplesizeequiv for power and sample size calculation
  for equivalence in milestone survival probability difference
* added rmpowerequiv and rmsamplesizeequiv for power and sample size calculation
  for equivalence in restricted mean survival time difference

# lrstat 0.2.3

* issued a warning message when unequal spacing is used with O'Brien-Fleming,
  Pocock, or Wang-Tsiatis boundaries in getBound
* renamed maxInformation to information in the overallResults data frame of the
  getDesign function output
* added simon2stage for Simon's two-stage design
* added nbstat to calculate the number of events and information for negative
  binomial rate ratio
* added nbpower to calculate the power for negative binomial rate ratio test
* added nbpower1s to calculate the power for one-sample negative binomial rate
* added nbpowerequiv to calculate the power for equivalence in negative binomial
  rate ratio
* added nbsamplesize to calculate the sample size for negative binomial rate
  ratio
* added nbsamplesize1s to calculate the sample size for one-sample negative
  binomial rate
* added nbsamplesizeequiv to calculate the sample size for equivalence in
  negative binomial rate ratio
* added runShinyApp to run a Shiny app for power and sample size calculation for
  log-rank tests
* added getDesignOneMean for group sequential design for one-sample mean
* added getDesignPairedMeanDiff for group sequential design for paired mean
  difference
* added getDesignPairedMeanRatio for group sequential design for paired mean
  ratio
* added getDesignMeanDiff for group sequential design for two-sample mean
  difference
* added getDesignMeanRatio for group sequential design for two-sample mean ratio
* added getDesignMeanDiffXO for group sequential design for mean difference in
  2x2 crossover
* added getDesignMeanRatioXO for group sequential design for mean ratio in 2x2
  crossover
* added getDesignPairedMeanDiffEquiv for group sequential design for equivalence
  in paired mean difference
* added getDesignPairedMeanRatioEquiv for group sequential design for
  equivalence in paired mean ratio
* added getDesignMeanDiffEquiv for group sequential design for equivalence in
  two-sample mean difference
* added getDesignMeanRatioEquiv for group sequential design for equivalence in
  two-sample mean ratio
* added getDesignMeanDiffXOEquiv for group sequential design for equivalence in
  mean difference in 2x2 crossover
* added getDesignMeanRatioXOEquiv for group sequential design for equivalence in
  mean ratio in 2x2 crossover
* added getDesignWilcoxon for group sequential design for two-sample Wilcoxon
  test
* added getDesignMeanDiffMMRM for two-sample mean difference at the last time
  point from the MMRM model
* added getDesignMeanDiffCarryover for direct treatment effects in crossover
  trials accounting for the carryover effects
* added getDesignANOVA for one-way analysis of variance
* added getDesignANOVAContrast for one-way analysis of variance contrast
* added getDesignRepeatedANOVA for one-way repeated analysis of variance
* added getDesignRepeatedANOVAContrast for one-way repeated analysis of variance
  contrast
* added getDesignTwoWayANOVA for two-way analysis of variance
* added getDesignOneSlope for group sequential design for one-sample slope
* added getDesignSlopeDiff for group sequential design for two-sample slope
  difference
* added getDesignSlopeDiffMMRM for two-sample slope difference from the MMRM
  model
* added getDesignOneProportion for group sequential design for one-sample
  proportion
* added getDesignPairedPropMcNemar for group sequential design for McNemar's
  test for paired proportions
* added getDesignRiskDiff for group sequential design for two-sample risk
  difference
* added getDesignRiskDiffExact for exact unconditional test for risk difference
* added getDesignRiskRatio for group sequential design for two-sample risk ratio
* added getDesignRiskRatioFM for the Farrington-Manning score test for risk
  ratio
* added getDesignRiskRatioExact for exact unconditional test for risk ratio
* added getDesignOddsRatio for group sequential design for two-sample odds ratio
* added getDesignRiskDiffEquiv for group sequential design for equivalence in
  two-sample risk difference
* added getDesignRiskDiffExactEquiv for exact unconditional test for equivalence
  in risk difference
* added getDesignRiskRatioEquiv for group sequential design for equivalence in
  two-sample risk ratio
* added getDesignRiskRatioExactEquiv for exact unconditional test for
  equivalence in risk ratio
* added getDesignOddsRatioEquiv for group sequential design for equivalence in
  two-sample odds ratio
* added getDesignFisherExact for Fisher's exact conditional test for two
  proportions
* added ClopperPearsonCI for Clopper-Pearson confidence interval for one-sample
  proportion
* added survQuantile for Brookmeyer-Crowley confidence interval of quantiles of
  right-censored time-to-event data
* added mTPI2Table for mTPI-2 decision table
* added BOINTable for BOIN decision table
* added mnRiskDiffCI for the Miettinen-Nurminen score confidence interval for
  two-sample risk difference
* added mnRiskRatioCI for the Miettinen-Nurminen score confidence interval for
  two-sample risk ratio
* added mnOddsRatioCI for the Miettinen-Nurminen score confidence interval for
  two-sample odds ratio
* added mnRateDiffCI for the Miettinen-Nurminen score confidence interval for
  two-sample rate difference
* added mnRateRatioCI for the Miettinen-Nurminen score confidence interval for
  two-sample rate ratio
* added getDesignOneMultinom for one-sample multinomial response
* added getDesignTwoMultinom for difference in two-sample multinomial response
* added getDesignTwoOrdinal for Wilcoxon test for two-sample ordinal response
* added getDesignOrderedBinom for Cochran-Armitage trend test for ordered
  multi-sample binomial response
* added getDesignUnorderedBinom for unordered multi-sample binomial response
* added getDesignUnorderedMultinom for unordered multi-sample multinomial
  response
* added getDesignLogistic for logistic regression
* added getDesignAgreement for Cohen's kappa agreement coefficient
* added getDesignOneRateExact for exact test of one-sample Poisson rate
* added ptpwexp for distribution function of truncated piecewise exponential
  distribution
* added rtpwexp for random number generation of truncated piecewise exponential
  distribution
* added hedgesg for Hedges' g effect size estimate and confidence interval
* added getDesignEquiv for a generic group sequential equivalence design
* added remlRiskDiff for REML estimates of individual proportions with specified
  risk difference
* added remlRiskRatio for REML estimates of individual proportions with
  specified risk ratio
* added remlOddsRatio for REML estimates of individual proportions with
  specified odds ratio
* added remlRateDiff for REML estimates of individual rates with specified rate
  difference
* added remlRateRatio for REML estimates of individual rates with specified rate
  ratio

# lrstat 0.2.2

* added the intnorm utility function to integrate a function with respect to a
  normal density
* added predictive power calculation to adaptDesign
* added the ftrunc function to calculate the adjusted p-values for truncated
  Holm, Hochberg, or Hommel procedures
* reuse the efficacy and futility stopping boundaries calculated under H1 for H0
  in lrsamplesize
* added capabilities to calculate Haybittle \& Peto boundaries in getDesign,
  lrpower, and lrsamplesize
* used informationRates as event fractions for conventional log-rank test and
  information fractions for weighted log-rank tests in lrpower and lrsamplesize
* match the number of events under H0 with the number of events under H1 for
  conventional log-rank test and match the information under H0 with the
  information under H1 for weighted log-rank tests
* removed getCriticalValues and getCumAlphaSpent function in lrstat.cpp
* adjusted test-f\_lrpower and test-f\_lrsamplesize to reflect changes to the
  definition of informationRates
* renamed informationTime to informationRates in lrsim for consistency

# lrstat 0.2.1

* used Markdown for Roxygen documentation
* renamed getAccrualDuration to getAccrualDurationFromN
* replaced predictEventOnly with predictTarget for the lrstat function
* added number of subjects reaching the maximum follow-up for fixed follow-up
  design for the lrstat function
* added efficacyStopping to the getBound function to improve coding efficiency
* added the getPower utility function to improve coding efficiency
* applied only equal spacing of looks for typeAlphaSpending of "OF", "P", or
  "WT" in the getBound function
* replaced the drift parameter with Imax and theta parameters in the getDesign
  function
* calculated alpha when critical values are not missing for the getDesign and
  lrpower functions
* added expected information under H0 to the getDesign function output
* added rejectPerStageH0, futilityPerStageH0, cumulativeRejectionH0,
  cumulativeFutilityH0, and attainedAlpha to the output of the getDesign
  function
* added the getCI function for parameter estimation after termination of a group
  sequential trial
* added the getRCI function to calculate repeated confidence intervals of a
  group sequential trial
* added the adaptDesign function for sample size re-estimation and conditional
  power calculation
* added the getADCI function for parameter estimation using the backward image
  method after termination of an adaptive group sequential trial
* added the getADRCI function to calculate repeated confidence intervals for an
  adaptive group sequential trial
* added the getCP function to calculate the conditional power when the parameter
  value may vary over time

# lrstat 0.2.0

* added fadjpdun to calculate the adjusted p-values for Dunnett-based graphical
  approaches.

# lrstat 0.1.15

* added fstp2seq for stepwise gatekeeping procedures with or without retesting
  for multiplicity problems involving two sequences of hypotheses.
* added fstdmix to obtain adjusted p-values for standard mixture gatekeeping
  procedures
* added fmodmix to obtain adjusted p-values for modified mixture gatekeeping
  procedures

# lrstat 0.1.14

* added the getAccrualDuration function to obtain the accrual duration to enroll
  the target number of subjects.
* added the getDurationFromNevents function to obtain a range of accrual
  duration to reach the target number of events.
* added the getNeventsFromHazardRatio function to obtain the required number of
  events given the hazard ratios under the null and alternative hypotheses for a
  group sequential design.
* allowed studyDuration < accrualDuration + followupTime for fixed follow-up in
  lrpower
* updated the handling of rounding for fixed follow-up design in lrsamplesize
* updated the handling of null hypothesis for fixed follow-up design

# lrstat 0.1.13

* added a rounding argument to lrsamplesize to round up the total sample size
  and events at each stage.
* added by treatment counts of events, counts, and subjects to lrpower output.
* added results under H0 to lrsamplesize output.

# lrstat 0.1.12

* used tolower to make typeAlphaSpending and typeBetaSpending into case
  insensitive inputs.

# lrstat 0.1.11

* added the Kaplan-Meier estimate of milestone survival, Greenwood variance
  estimate, difference in milestone survival, and Z test statistic for survival
  difference.

# lrstat 0.1.10

* added drift parameter to the getDesign function to compute power given the
  drift parameter.
* updated the repeatedPValue function to respect the range of repeated p-values
  and to allow matrix input of raw p-values.
* removed repeatedPValueFlag from the fseqbon function.
* removed numSubintervals from the caltime function.
* updated the description of selected functions, parameters, and output.

# lrstat 0.1.9

* added fwgtmat and fadjpsim to calculate the adjusted p-values for Simes-based
  graphical approaches.
* updated the print method for design, lrpower, and lrsim.

# lrstat 0.1.8

* added spendingTime to getDesign, lrpower, and lrsamplesize to allow the error
  spending time to be different from the information time.
* rewrote lrsamplesize to simplify and accelerate the computation for
  typeOfComputation == "Schoenfeld".
* added getBound to obtain the efficacy stopping boundaries for a group
  sequential design allowing the error spending time to be different from the
  information time.
* added fadjpbon to obtain the adjusted p-values for graphical approaches using
  weighted Bonferroni tests for fixed design.
* added updateGraph to update the weights and transition matrix after removing a
  hypothesis from the set of indices of yet to be rejected null hypotheses.
* added repeatedPValue to Obtain the repeated p-values for a group sequential
  design based on a given alpha spending function.
* added fseqbon to obtain the test results for group sequential trials using
  graphical approaches based on weighted Bonferroni tests with the option to
  provide repeated p-values for each hypothesis over time.
* added lrsim3a to perform simulation for three-arm group sequential trials
  based on weighted log-rank test. The looks are driven by the total number of
  events in Arm A and Arm C combined.
* added lrsim2e to perform simulation for two-endpoint two-arm group sequential
  trials based on weighted log-rank test. The first few looks are driven by the
  total number of PFS events in two arms combined, and the subsequent looks are
  driven by the total number of OS events in two arms combined.
* added lrsim2e3a to perform simulation for two-endpoint three-arm group
  sequential trials based on weighted log-rank test. The first few looks are
  driven by the total number of PFS events in Arm A and Arm C combined, and the
  subsequent looks are driven by the total number of OS events in Arm A and Arm
  C combined.

# lrstat 0.1.7

* added getDesign for creating a generic group sequential design with constant
  treatment effect over time.

# lrstat 0.1.6

* added capability for performing noninferiority tests in lrpower, lrsamplesize,
  and lrsim.
* added capability for simulating analyses based on calendar times in lrsim.
* adjusted the critical value at the final look if the observed total number of
  events is less than the planned total number of events in lrsim.
* retained summary statistics for all stages even after crossing the efficacy
  and futility boundaries in lrsim.
* added number of dropouts to lrpower/lrsamplesize and lrsim output.
* added Schoenfeld method for proportional hazards and conventional log-rank
  test in lrpower and lrsamplesize.

# lrstat 0.1.5

* replaced Inf with 6 and -Inf with -6 for test statistic stopping boundaries to
  avoid potential memory issue.

# lrstat 0.1.4

New features

* added capability for lrstat to calculate hazard ratios from weighted Cox
  regression model.
* added capability for lrsamplesize to calculate absolute accrual rate from
  relative accrual rate given power, accrual duration, and follow-up duration.

Bug fixes

* used specified informationRates to calculate Wang-Tsiatis boundaries.
* used hazard ratios from weighted Cox regression model to determine crossing
  boundaries on the hazard ratio scale for lrpower.
* replaced stratum-specific output with overall results for lrstat.
* removed hazard ratio estimate from weighted log-rank test from lrsim output.

# lrstat 0.1.3

* added more statistics to lrpower output.

# lrstat 0.1.2

New features

* added capability for lrpower and lrsamplesize to use error spending functions.
* added more statistics to lrstat, lrpower and lrsim output.
* allowed user to specify numSubintervals to control approximation.

Bug fixes

* added parameter checking for lrpower, lrsamplesize, and lrsim.
* added test files.
* added print\_lrpower.R to print lrpower objects.
* used informationTime instead of informationRates in lrsim to differentiate
  information based on weighted log-rank tests score statistic variance from
  information based on number of events.
* renamed sumstat to overview in lrsim output.

# lrstat 0.1.1

* fixed hyperlinks.

# lrstat 0.1.0

* Initial release.
