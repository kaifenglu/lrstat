# lrstat 0.1.6

New features

- Add capability for performing noninferiority tests in lrpower, lrsamplesize, and lrsim.
- Add capability for simulating analyses based on calendar times in lrsim.
- Adjust the critical value at the final look if the observed total number of events is less than the planned total number of events in lrsim.
- Retain summary statistics for all stages even after crossing the efficacy and futility boundaries in lrsim.
- Add number of dropouts 

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


- Fix hyperlinks


# lrstat 0.1.0

- Initial release
