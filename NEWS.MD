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
