# Conditional Power Allowing for Varying Parameter Values

Obtains the conditional power for specified incremental information
given the interim results, parameter values, and data-dependent changes
in the error spending function, as well as the number and spacing of
interim looks.

## Usage

``` r
getCP(
  INew = NA_real_,
  L = NA_integer_,
  zL = NA_real_,
  theta = NA_real_,
  IMax = NA_real_,
  kMax = NA_integer_,
  informationRates = NA_real_,
  efficacyStopping = NA_integer_,
  futilityStopping = NA_integer_,
  criticalValues = NA_real_,
  alpha = 0.025,
  typeAlphaSpending = "sfOF",
  parameterAlphaSpending = NA_real_,
  userAlphaSpending = NA_real_,
  futilityBounds = NA_real_,
  typeBetaSpending = "none",
  parameterBetaSpending = NA_real_,
  spendingTime = NA_real_,
  MullerSchafer = 0L,
  kNew = NA_integer_,
  informationRatesNew = NA_real_,
  efficacyStoppingNew = NA_integer_,
  futilityStoppingNew = NA_integer_,
  typeAlphaSpendingNew = "sfOF",
  parameterAlphaSpendingNew = NA_real_,
  typeBetaSpendingNew = "none",
  parameterBetaSpendingNew = NA_real_,
  spendingTimeNew = NA_real_,
  varianceRatio = 1
)
```

## Arguments

- INew:

  The maximum information of the secondary trial.

- L:

  The interim adaptation look of the primary trial.

- zL:

  The z-test statistic at the interim adaptation look of the primary
  trial.

- theta:

  A scalar or a vector of parameter values of length `kMax + kMax - L`
  if `MullerSchafer = FALSE` or length `kMax + kNew` if
  `MullerSchafer = TRUE`.

- IMax:

  The maximum information of the primary trial.

- kMax:

  The maximum number of stages of the primary trial.

- informationRates:

  The information rates of the primary trial.

- efficacyStopping:

  Indicators of whether efficacy stopping is allowed at each stage of
  the primary trial. Defaults to true if left unspecified.

- futilityStopping:

  Indicators of whether futility stopping is allowed at each stage of
  the primary trial. Defaults to true if left unspecified.

- criticalValues:

  The upper boundaries on the z-test statistic scale for efficacy
  stopping for the primary trial.

- alpha:

  The significance level of the primary trial. Defaults to 0.025.

- typeAlphaSpending:

  The type of alpha spending for the primary trial. One of the
  following: "OF" for O'Brien-Fleming boundaries, "P" for Pocock
  boundaries, "WT" for Wang & Tsiatis boundaries, "sfOF" for
  O'Brien-Fleming type spending function, "sfP" for Pocock type spending
  function, "sfKD" for Kim & DeMets spending function, "sfHSD" for
  Hwang, Shi & DeCani spending function, "user" for user defined
  spending, and "none" for no early efficacy stopping. Defaults to
  "sfOF".

- parameterAlphaSpending:

  The parameter value of alpha spending for the primary trial.
  Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".

- userAlphaSpending:

  The user defined alpha spending for the primary trial. Cumulative
  alpha spent up to each stage.

- futilityBounds:

  The lower boundaries on the z-test statistic scale for futility
  stopping for the primary trial. Defaults to `rep(-6, kMax-1)` if left
  unspecified.

- typeBetaSpending:

  The type of beta spending for the primary trial. One of the following:
  "sfOF" for O'Brien-Fleming type spending function, "sfP" for Pocock
  type spending function, "sfKD" for Kim & DeMets spending function,
  "sfHSD" for Hwang, Shi & DeCani spending function, and "none" for no
  early futility stopping. Defaults to "none".

- parameterBetaSpending:

  The parameter value of beta spending for the primary trial.
  Corresponds to rho for "sfKD", and gamma for "sfHSD".

- spendingTime:

  The error spending time of the primary trial. Defaults to missing, in
  which case, it is the same as `informationRates`.

- MullerSchafer:

  Whether to use the Muller and Schafer (2001) method for trial
  adaptation.

- kNew:

  The number of looks of the secondary trial.

- informationRatesNew:

  The spacing of looks of the secondary trial.

- efficacyStoppingNew:

  The indicators of whether efficacy stopping is allowed at each look of
  the secondary trial. Defaults to true if left unspecified.

- futilityStoppingNew:

  The indicators of whether futility stopping is allowed at each look of
  the secondary trial. Defaults to true if left unspecified.

- typeAlphaSpendingNew:

  The type of alpha spending for the secondary trial. One of the
  following: "OF" for O'Brien-Fleming boundaries, "P" for Pocock
  boundaries, "WT" for Wang & Tsiatis boundaries, "sfOF" for
  O'Brien-Fleming type spending function, "sfP" for Pocock type spending
  function, "sfKD" for Kim & DeMets spending function, "sfHSD" for
  Hwang, Shi & DeCani spending function, and "none" for no early
  efficacy stopping. Defaults to "sfOF".

- parameterAlphaSpendingNew:

  The parameter value of alpha spending for the secondary trial.
  Corresponds to Delta for "WT", rho for "sfKD", and gamma for "sfHSD".

- typeBetaSpendingNew:

  The type of beta spending for the secondary trial. One of the
  following: "sfOF" for O'Brien-Fleming type spending function, "sfP"
  for Pocock type spending function, "sfKD" for Kim & DeMets spending
  function, "sfHSD" for Hwang, Shi & DeCani spending function, and
  "none" for no early futility stopping. Defaults to "none".

- parameterBetaSpendingNew:

  The parameter value of beta spending for the secondary trial.
  Corresponds to rho for "sfKD", and gamma for "sfHSD".

- spendingTimeNew:

  The error spending time of the secondary trial. Defaults to missing,
  in which case, it is the same as `informationRatesNew`.

- varianceRatio:

  The ratio of the variance under H0 to the variance under H1.

## Value

The conditional power given the interim results, parameter values, and
data-dependent design changes.

## References

Cyrus R. Mehta and Stuart J. Pocock. Adaptive increase in sample size
when interim results are promising: A practical guide with examples.
Stat Med. 2011;30:3267–3284.

## See also

[`getDesign`](https://kaifenglu.github.io/lrstat/reference/getDesign.md)

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
# Conditional power calculation with delayed treatment effect

# Two interim analyses have occurred with 179 and 266 events,
# respectively. The observed hazard ratio at the second interim
# look is 0.81.

trialsdt = as.Date("2020-03-04")                       # trial start date
iadt = c(as.Date("2022-02-01"), as.Date("2022-11-01")) # interim dates
mo1 = as.numeric(iadt - trialsdt + 1)/30.4375          # interim months

# Assume a piecewise Poisson enrollment process with a 8-month ramp-up
# and 521 patients were enrolled after 17.94 months
N = 521                   # total number of patients
Ta = 17.94                # enrollment duration
Ta1 = 8                   # assumed end of enrollment ramp-up
enrate = N / (Ta - Ta1/2) # enrollment rate after ramp-up

# Assume a median survival of 16.7 months for the control group, a
# 5-month delay in treatment effect, and a hazard ratio of 0.7 after
# the delay
lam1 = log(2)/16.7  # control group hazard of exponential distribution
t1 = 5              # months of delay in treatment effect
hr = 0.7            # hazard ratio after delay
lam2 = hr*lam1      # treatment group hazard after delay

# Assume an annual dropout rate of 5%
gam = -log(1-0.05)/12  # hazard for dropout

# The original target number of events was 298 and the new target is 335
mo2 <- caltime(
  nevents = c(298, 335),
  allocationRatioPlanned = 1,
  accrualTime = seq(0, Ta1),
  accrualIntensity = enrate*seq(1, Ta1+1)/(Ta1+1),
  piecewiseSurvivalTime = c(0, t1),
  lambda1 = c(lam1, lam2),
  lambda2 = c(lam1, lam1),
  gamma1 = gam,
  gamma2 = gam,
  accrualDuration = Ta,
  followupTime = 1000)

# expected number of events and average hazard ratios
(lr1 <- lrstat(
  time = c(mo1, mo2),
  accrualTime = seq(0, Ta1),
  accrualIntensity = enrate*seq(1, Ta1+1)/(Ta1+1),
  piecewiseSurvivalTime = c(0, t1),
  lambda1 = c(lam1, lam2),
  lambda2 = c(lam1, lam1),
  gamma1 = gam,
  gamma2 = gam,
  accrualDuration = Ta,
  followupTime = 1000,
  predictTarget = 3))
#>       time subjects  nevents  nevents1  nevents2 ndropouts ndropouts1
#> 1 22.99795      521 184.4200  85.65447  98.76552  20.64012   10.46884
#> 2 31.96715      521 267.7201 122.68031 145.03981  30.85289   15.91609
#> 3 36.23725      521 298.0000 136.77727 161.22273  34.59342   17.99003
#> 4 42.68854      521 335.0000 154.60965 180.39035  39.19088   20.61354
#>   ndropouts2 nfmax nfmax1 nfmax2     uscore   vscore  logRankZ hazardRatioH0
#> 1   10.17128     0      0      0  -7.773866 46.07799 -1.145223             1
#> 2   14.93680     0      0      0 -15.171213 66.79076 -1.856360             1
#> 3   16.60338     0      0      0 -17.860967 74.26549 -2.072581             1
#> 4   18.57735     0      0      0 -21.139606 83.32376 -2.315861             1
#>          HR     vlogHR    zlogHR
#> 1 0.8445734 0.02180878 -1.143865
#> 2 0.7966399 0.01506393 -1.852382
#> 3 0.7861833 0.01353711 -2.067617
#> 4 0.7760959 0.01204385 -2.309722


hr2 = 0.81                    # observed hazard ratio at interim 2
z2 = (-log(hr2))*sqrt(266/4)  # corresponding z-test statistic value

# expected mean of -log(HR) at the original looks and the new final look
theta = -log(lr1$HR[c(1,2,3,4)])

# conditional power with sample size increase
getCP(INew = (335 - 266)/4,
      L = 2, zL = z2, theta = theta,
      IMax = 298/4, kMax = 3,
      informationRates = c(179, 266, 298)/298,
      alpha = 0.025, typeAlphaSpending = "sfOF")
#> [1] 0.5550156
```
