# Error Spending

Obtains the error spent at given spending times for the specified error
spending function.

## Usage

``` r
errorSpent(t, error = 0.025, sf = "sfOF", sfpar = NA_real_)
```

## Arguments

- t:

  A vector of spending times, typically equal to information fractions.

- error:

  The total error to spend.

- sf:

  The spending function. One of the following: "sfOF" for
  O'Brien-Fleming type spending function, "sfP" for Pocock type spending
  function, "sfKD" for Kim & DeMets spending function, and "sfHSD" for
  Hwang, Shi & DeCani spending function. Defaults to "sfOF".

- sfpar:

  The parameter for the spending function. Corresponds to rho for "sfKD"
  and gamma for "sfHSD".

## Value

A vector of errors spent up to the interim look.

## Details

This function implements a variety of error spending functions commonly
used in group sequential designs, assuming one-sided hypothesis testing.

**O'Brien-Fleming-Type Spending Function**

This spending function allocates very little alpha early on and more
alpha later in the trial. It is defined as: \$\$ \alpha(t) = 2 -
2\Phi\left(\frac{z\_{\alpha/2}}{\sqrt{t}}\right), \$\$ where \\\Phi\\ is
the standard normal cumulative distribution function, \\z\_{\alpha/2}\\
is the critical value from the standard normal distribution, and \\t \in
\[0, 1\]\\ denotes the information fraction.

**Pocock-Type Spending Function**

This function spends alpha more evenly throughout the study: \$\$
\alpha(t) = \alpha \log(1 + (e - 1)t), \$\$ where \\e\\ is Euler's
number (approximately 2.718).

**Kim and DeMets Power-Type Spending Function**

This family of spending functions is defined as: \$\$ \alpha(t) = \alpha
t^{\rho}, \quad \rho \> 0. \$\$

- When \\\rho = 1\\, the function mimics Pocock-type boundaries.

- When \\\rho = 3\\, it approximates O’Brien-Fleming-type boundaries.

**Hwang, Shih, and DeCani Spending Function**

This flexible family of functions is given by: \$\$ \alpha(t) =
\begin{cases} \alpha \frac{1 - e^{-\gamma t}}{1 - e^{-\gamma}}, &
\text{if } \gamma \ne 0 \\ \alpha t, & \text{if } \gamma = 0.
\end{cases} \$\$

- When \\\gamma = -4\\, the spending function resembles O’Brien-Fleming
  boundaries.

- When \\\gamma = 1\\, it resembles Pocock boundaries.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>

## Examples

``` r
errorSpent(t = 0.5, error = 0.025, sf = "sfOF")
#> [1] 0.001525323

errorSpent(t = c(0.5, 0.75, 1), error = 0.025, sf = "sfHSD", sfpar = -4)
#> [1] 0.002980073 0.008902144 0.025000000
```
