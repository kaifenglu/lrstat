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
