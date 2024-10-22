# lrstat 0.2.10

- use newton-raphson instead of the vmmin algorithm for liferegr and phregr
- removed rpsft with all treatment switching functions moved to trtswitch package
- add the survfit_phregr function to obtain the survivor curve estimates
- add the residuals_phregr function to obtain the martingale, deviance, score, and schoenfeld residuals
- add the shilong dataset
- replace std::ceil(x) with std::ceil(x - 1.0e-12) to handle superficial precision in binary presentation of a double number
- allow rep and stratum parameters to represent more than one character variables in the source data
- add the logisregr function for logistic regression of a binary endpoint.

# lrstat 0.2.9

- fix the runtime error: nan is outside the range of representable values of type 'int'
- add the binary_tte_sim function to simulate two endpoints with one being 
  a binary endpoint and the other being a time-to-event endpoint
- add the rpsft function to estimate hazard ratio using rank-preserving structured failure time model

# lrstat 0.2.8

- update std::sort handling for Rcpp character vectors
- use static cast to explicitly cast "long"" of the member function size() 
  to "int" and cast "double" of the floor, ceil, round functions to "int"

# lrstat 0.2.7

- add the aml, heart, and tobin datasets
- in basket.cpp, replace y[k] = R::rbinom(1, p[j]) with u = R::runif(0,1); y[k] = u < p[j] ? 1 : 0;
- in misc.cpp, replace const in x = NA_REAL with const int x = NA_INTEGER for x = n, n1, y1, n2, y2
- add delete[] for iwork and work pointers in the quad function
- add simonBayesAnalysis and simonBayesSim for Simon's Bayesian analysis of basket discovery trials
- update the algorithm for finding the backward image (J, zJ) in getADCI
- update the initial values for finding the confidence interval in getADCI
- add bmini for optimization and invsympd for inverting symmetric and positive defined matrices to utilities.cpp
- add phregr to estimate hazard ratios from the Cox model for right-censored or counting process data with robust variance
- add liferegr for parametric regression of failure time data with uncensored, right censored, left censored, and interval censored data

# lrstat 0.2.6

- replace "R_isnancpp(x)" with "x == NA_INTEGER" for integer parameter x
- replace "!all.equal(x$settings$spendingTime, s$informationRates)" with "!all(x$settings$spendingTime == s$informationRates)" in prints.R

# lrstat 0.2.5

- update survQuantile in getDesignMeans.R to handle new censoring and inestimable quantiles
- update simon2stage to increase the search space of sample size
- add zstatRiskDiff, zstatRiskRatio, zstatOddsRatio, zstatRateDiff, zstatRateRatio for the score test of two-sample crude rates and exposure-adjusted incidence rates
- add mini for minimization of a univariate function over a finite interval
- update powerRiskDiffExact, samplesizeRiskDiffExact, powerRiskRatioExact, samplesizeRiskRatioExact, powerRiskDiffExactEquiv, samplesizeRiskDiffExactEquiv, powerRiskRatioExactEquiv, and samplesizeRiskRatioExactEquiv to divide the interval of the response rate in the control group into 100 subintervals and obtain the minimum value within each subinterval in order to obtain the global minimum. This replaces the random draw of 500 values of the control group response rate for global minimum. In addition, the attained alpha calculation is updated to use the global maximum
- add riskDiffExactPValue, riskDiffExactCI, riskRatioExactPValue, and riskRatioExactCI for exact unconditional test and confidence limits based on the Miettinen & Nurminen score statistics
- add rawdata to be used as an input data frame to test survival analysis functions
- add kmest and kmdiff for Kaplan-Meier estimates of survival curves and milestone survival difference
- add lrtest for log-rank test using the Fleming-Harrington family of weights
- add rmest and rmdiff for estimates of restricted mean survival times and restricted mean survival time difference
- add lower.tail and log.p arguments to ptpwexp and qtpwexp functions
- add the fquantile function to obtain the quantiles of a survival function
- add the pwexploglik function to obtain the profile log-likelihood function for the change points in the piecewise exponential approximation to a survival function
- add the pwexpcuts function to obtain a piecewise exponential distribution that approximates a survival distribution
- update getDesignMeanDiffMMRM and getDesignSlopeDiffMMRM to add accrual and dropout information

# lrstat 0.2.4

- update the check for pconfigs of getDesignLogistic to avoid an error on macOS
- replace which_min with which_max in lrsim2e and lrsim2e3a
- replace qromb with quad for numerical integration
- use numerical integration for lrstat
- rename kmest to kmstat and use numerical integration
- add kmpower and kmsamplesize for power and sample size calculation for difference in milestone survival probability
- add kmpower1s and kmsamplesize1s for power and sample size calculation for one-sample milestone survival probability
- add rmst, covrmst, and rmstat for restricted mean survival time analysis
- add rmpower and rmsamplesize for power and sample size calculation for difference in restricted mean survival time
- add rmpower1s and rmsamplesize1s for power and sample size calculation for one-sample restricted mean survival time
- update rmstat, rmpower, rmsamplesize, rmpower1s, and rmsamplesize1s to account for stratification
- add lrpowerequiv and lrsamplesizeequiv for power and sample size calculation for equivalence in hazard ratio
- add kmpowerequiv and kmsamplesizeequiv for power and sample size calculation for equivalence in milestone survival probability difference
- add rmpowerequiv and rmsamplesizeequiv for power and sample size calculation for equivalence in restricted mean survival time difference


# lrstat 0.2.3

- issue a warning message when unequal spacing is used with O'Brien-Fleming, Pocock, or Wang-Tsiatis boundaries in getBound
- rename maxInformation to information in the overallResults data frame of the getDesign function output
- add simon2stage for Simon's two-stage design
- add nbstat to calculate the number of events and information for negative binomial rate ratio 
- add nbpower to calculate the power for negative binomial rate ratio test
- add nbpower1s to calculate the power for one-sample negative binomial rate
- add nbpowerequiv to calculate the power for equivalence in negative binomial rate ratio
- add nbsamplesize to calculate the sample size for negative binomial rate ratio
- add nbsamplesize1s to calculate the sample size for one-sample negative binomial rate
- add nbsamplesizeequiv to calculate the sample size for equivalence in negative binomial rate ratio
- add runShinyApp to run a Shiny app for power and sample size calculation for log-rank tests
- add getDesignOneMean for group sequential design for one-sample mean
- add getDesignPairedMeanDiff for group sequential design for paired mean difference
- add getDesignPairedMeanRatio for group sequential design for paired mean ratio
- add getDesignMeanDiff for group sequential design for two-sample mean difference
- add getDesignMeanRatio for group sequential design for two-sample mean ratio
- add getDesignMeanDiffXO for group sequential design for mean difference in 2x2 crossover
- add getDesignMeanRatioXO for group sequential design for mean ratio in 2x2 crossover
- add getDesignPairedMeanDiffEquiv for group sequential design for equivalence in paired mean difference
- add getDesignPairedMeanRatioEquiv for group sequential design for equivalence in paired mean ratio
- add getDesignMeanDiffEquiv for group sequential design for equivalence in two-sample mean difference
- add getDesignMeanRatioEquiv for group sequential design for equivalence in two-sample mean ratio
- add getDesignMeanDiffXOEquiv for group sequential design for equivalence in mean difference in 2x2 crossover
- add getDesignMeanRatioXOEquiv for group sequential design for equivalence in mean ratio in 2x2 crossover
- add getDesignWilcoxon for group sequential design for two-sample Wilcoxon test
- add getDesignMeanDiffMMRM for two-sample mean difference at the last time point from the MMRM model
- add getDesignMeanDiffCarryover for direct treatment effects in crossover trials accounting for the carryover effects
- add getDesignANOVA for one-way analysis of variance
- add getDesignANOVAContrast for one-way analysis of variance contrast
- add getDesignRepeatedANOVA for one-way repeated analysis of variance
- add getDesignRepeatedANOVAContrast for one-way repeated analysis of variance contrast
- add getDesignTwoWayANOVA for two-way analysis of variance
- add getDesignOneSlope for group sequential design for one-sample slope
- add getDesignSlopeDiff for group sequential design for two-sample slope difference
- add getDesignSlopeDiffMMRM for two-sample slope difference from the MMRM model
- add getDesignOneProportion for group sequential design for one-sample proportion
- add getDesignPairedPropMcNemar for group sequential design for McNemar's test for paired proportions
- add getDesignRiskDiff for group sequential design for two-sample risk difference
- add getDesignRiskDiffExact for exact unconditional test for risk difference
- add getDesignRiskRatio for group sequential design for two-sample risk ratio
- add getDesignRiskRatioFM for the Farrington-Manning score test for risk ratio
- add getDesignRiskRatioExact for exact unconditional test for risk ratio
- add getDesignOddsRatio for group sequential design for two-sample odds ratio
- add getDesignRiskDiffEquiv for group sequential design for equivalence in two-sample risk difference
- add getDesignRiskDiffExactEquiv for exact unconditional test for equivalence in risk difference
- add getDesignRiskRatioEquiv for group sequential design for equivalence in two-sample risk ratio
- add getDesignRiskRatioExactEquiv for exact unconditional test for equivalence in risk ratio
- add getDesignOddsRatioEquiv for group sequential design for equivalence in two-sample odds ratio
- add getDesignFisherExact for Fisher's exact conditional test for two proportions
- add ClopperPearsonCI for Clopper-Pearson confidence interval for one-sample proportion
- add survQuantile for Brookmeyer-Crowley confidence interval of quantiles of right-censored time-to-event data
- add mTPI2Table for mTPI-2 decision table
- add BOINTable for BOIN decision table
- add mnRiskDiffCI for the Miettinen-Nurminen score confidence interval for two-sample risk difference
- add mnRiskRatioCI for the Miettinen-Nurminen score confidence interval for two-sample risk ratio
- add mnOddsRatioCI for the Miettinen-Nurminen score confidence interval for two-sample odds ratio
- add mnRateDiffCI for the Miettinen-Nurminen score confidence interval for two-sample rate difference
- add mnRateRatioCI for the Miettinen-Nurminen score confidence interval for two-sample rate ratio
- add getDesignOneMultinom for one-sample multinomial response
- add getDesignTwoMultinom for difference in two-sample multinomial response
- add getDesignTwoOrdinal for Wilcoxon test for two-sample ordinal response
- add getDesignOrderedBinom for Cochran-Armitage trend test for ordered multi-sample binomial response
- add getDesignUnorderedBinom for unordered multi-sample binomial response
- add getDesignUnorderedMultinom for unordered multi-sample multinomial response
- add getDesignLogistic for logistic regression
- add getDesignAgreement for Cohen's kappa agreement coefficient
- add getDesignOneRateExact for exact test of one-sample Poisson rate
- add ptpwexp for distribution function of truncated piecewise exponential distribution
- add rtpwexp for random number generation of truncated piecewise exponential distribution
- add hedgesg for Hedges' g effect size estimate and confidence interval
- add getDesignEquiv for a generic group sequential equivalence design
- add remlRiskDiff for REML estimates of individual proportions with specified risk difference
- add remlRiskRatio for REML estimates of individual proportions with specified risk ratio
- add remlOddsRatio for	REML estimates of individual proportions with specified odds ratio
- add remlRateDiff for REML estimates of individual rates with specified rate difference
- add remlRateRatio for REML estimates of individual rates with specified rate ratio

# lrstat 0.2.2

- add the intnorm utility function to integrate a function with respect to a normal density
- add predictive power calculation to adaptDesign
- add the ftrunc function to calculate the adjusted p-values for truncated Holm, Hochberg, or Hommel procedures
- reuse the efficacy and futility stopping boundaries calculated under H1 for H0 in lrsamplesize
- add capabilities to calculate Haybittle & Peto boundaries in getDesign, lrpower, and lrsamplesize
- use informationRates as event fractions for conventional log-rank test and information fractions for weighted log-rank tests in lrpower and lrsamplesize
- match the number of events under H0 with the number of events under H1 for conventional log-rank test and match the information under H0 with the information under H1 for weighted log-rank tests
- remove getCriticalValues and getCumAlphaSpent function in lrstat.cpp
- adjust test-f_lrpower and test-f_lrsamplesize to reflect changes to the definition of informationRates
- rename informationTime to informationRates in lrsim for consistency


# lrstat 0.2.1

- use markdown for Roxygen documentations
- rename getAccrualDuration to getAccrualDurationFromN
- replace predictEventOnly with predictTarget for the lrstat function
- add number of subjects reaching the maximum follow-up for fixed follow-up design for the lrstat function
- add efficacyStopping to the getBound function to improve coding efficiency
- add the getPower utility function to improve coding efficiency
- apply only equal spacing of looks for typeAlphaSpending of "OF", "P", or "WT" in the getBound function
- replace the drift parameter with Imax and theta parameters in the getDesign function
- calculate alpha when critical values are not missing for the getDesign and lrpower functions
- add expected information under H0 to the getDesign function output
- add rejectPerStageH0, futilityPerStageH0, cumulativeRejectionH0, cumulativeFutilityH0, and attainedAlpha to the output of the getDesign function
- add the getCI function for parameter estimation after termination of a group sequential trial
- add the getRCI function to calculate repeated confidence intervals of a group sequential trial
- add the adaptDesign function for sample size re-estimation and conditional power calculation
- add the getADCI function for parameter estimation using the backward image method after termination of an adaptive group sequential trial
- add the getADRCI function to calculate repeated confidence intervals for an adpaptive group sequential trial
- add the getCP function to calculate the conditional power when the parameter value may vary over time

# lrstat 0.2.0

- add fadjpdun to calculate the adjusted p-values for Dunnett-based graphical approaches. 

# lrstat 0.1.15

- add fstp2seq for stepwise gatekeeping procedures with or without retesting for multiplicity problems involving two sequences of hypotheses.
- add fstdmix to obtain adjusted p-values for standard mixture gatekeeping procedures
- add fmodmix to obtain adjusted p-values for modified mixture gatekeeping procedures

# lrstat 0.1.14

- add the getAccrualDuration function to obtain the accrual duration to enroll the target number of subjects.
- add the getDurationFromNevents function to obtain a range of accrual duration to reach the target number of events.
- add the getNeventsFromHazardRatio function to obtain the required number of events given the hazard ratios under the null and alternative hypotheses for a group sequential design.
- allow studyDuration < accrualDuration + followupTime for fixed follow-up in lrpower
- update the handling of rounding for fixed follow-up design in lrsamplesize
- update the handling of null hypothesis for fixed follow-up design

# lrstat 0.1.13

- add a rounding argument to lrsamplesize to round up the total sample size and events at each stage.
- add by treatment counts of events, counts, and subjects to lrpower output.
- add results under H0 to lrsamplesize output.

# lrstat 0.1.12

- use tolower to make typeAlphaSpending and typeBetaSpending into case insensitive inputs.

# lrstat 0.1.11

- Add Kaplan-Meier estimate of milestone survival, Greendwood variance estimate, difference in milestone survival, and Z test statistic for survival difference.


# lrstat 0.1.10

- Add drift parameter to the getDesign function to compute power given the drift parameter.
- Update the repeatedPValue function to respect the range of repeated p-values and to allow matrix input of raw p-values.
- Remove repeatedPValueFlag from the fseqbon function.
- Remove numSubintervals from the caltime function.
- Update the description of selected functions, parameters, and output.

# lrstat 0.1.9

- Add fwgtmat and fadjpsim to calculate the adjusted p-values for Simes-based graphical approaches.
- update the print method for design, lrpower, and lrsim.

# lrstat 0.1.8

- Add spendingTime to getDesign, lrpower, and lrsamplesize to allow the error spending time to be different from the information time.
- Rewrite lrsamplesize to simplify and accelerate the computation for typeOfComputation == "Schoenfeld".
- Add getBound to obtain the efficacy stopping boundaries for a group sequential design allowing the error spending time to be different from the information time.
- Add fadjpbon to obtain the adjusted p-values for graphical approaches using weighted Bonferroni tests for fixed design.
- Add updateGraph to update the weights and transition matrix after removing a hypothesis from the set of indices of yet to be rejected null hypotheses. 
- Add repeatedPValue to Obtain the repeated p-values for a group sequential design based on a given alpha spending function.
- Add fseqbon to obtain the test results for group sequential trials using graphical approaches based on weighted Bonferroni tests with the option to provide repeated p-values for each hypothesis over time.
- Add lrsim3a to perform simulation for three-arm group sequential trials based on weighted log-rank test. The looks are driven by the total number of events in Arm A and Arm C combined.
- Add lrsim2e to perform simulation for two-endpoint two-arm group sequential trials based on weighted log-rank test. The first few looks are driven by the total number of PFS events in two arms combined, and the subsequent looks are driven by the total number of OS events in two arms combined.
- Add lrsim2e3a to perform simulation for two-endpoint three-arm group sequential trials based on weighted log-rank test. The first few looks are driven by the total number of PFS events in Arm A and Arm C combined, and the subsequent looks are driven by the total number of OS events in Arm A and Arm C combined.


# lrstat 0.1.7

- Add getDesign for creating a generic group sequential design with constant treatment effect over time.

# lrstat 0.1.6

- Add capability for performing noninferiority tests in lrpower, lrsamplesize, and lrsim.
- Add capability for simulating analyses based on calendar times in lrsim.
- Adjust the critical value at the final look if the observed total number of events is less than the planned total number of events in lrsim.
- Retain summary statistics for all stages even after crossing the efficacy and futility boundaries in lrsim.
- Add number of dropouts to lrpower/lrsamplesize and lrsim output.
- Add Schoenfeld method for proportional hazards and conventional log-rank test in lrpower and lrsamplesize.

# lrstat 0.1.5

- Replace Inf with 6 and -Inf with -6 for test statistic stopping boundaries to avoid potential memory issue.

# lrstat 0.1.4

New features

- Add capability for lrstat to calculate hazard ratios from weighted Cox regression model.
- Add capability for lrsamplesize to calculate absolute accrual rate from  relative accrual rate given power, accrual duration, and follow-up duration.

Bug fixes

- Use specified informationRates to calculate Wang-Tsiatis boundaries.
- Use hazard ratios from weighted Cox regression model to determine crossing boundaries on the hazard ratio scale for lrpower.
- Replace stratum-specific output with overall results for lrstat.
- Remove hazard ratio estimate from weighted log-rank test from lrsim output.


# lrstat 0.1.3

- Add more statistics to lrpower output.


# lrstat 0.1.2

New features

- Add capability for lrpower and lrsamplesize to use error spending functions.
- Add more statistics to lrstat, lrpower and lrsim output.
- Allow user to specify numSubintervals to control approximation.

Bug fixes

- Add parameter checking for lrpower, lrsamplesize, and lrsim.
- Add test files.
- Add print_lrpower.R to print lrpower objects.
- Use informationTime instead of informationRates in lrsim to differentiate information based on weighted log-rank tests score statistic variance from information based on number of events.
- Rename sumstat to overview in lrsim output.


# lrstat 0.1.1


- Fix hyperlinks.


# lrstat 0.1.0

- Initial release.
