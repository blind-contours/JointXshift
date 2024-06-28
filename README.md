
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`JointXshift`

<!-- badges: start -->

[![R-CMD-check](https://github.com/blind-contours/JointXshift/workflows/R-CMD-check/badge.svg)](https://github.com/blind-contours/JointXshift/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/blind-contours/JointXshift/master.svg)](https://codecov.io/github/blind-contours/JointXshift?branch=master)
[![CRAN](https://www.r-pkg.org/badges/version/JointXshift)](https://www.r-pkg.org/pkg/JointXshift)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/JointXshift)](https://CRAN.R-project.org/package=JointXshift)
[![CRAN total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/JointXshift)](https://CRAN.R-project.org/package=JointXshift)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4070042.svg)](https://doi.org/10.5281/zenodo.4070042) -->
<!-- [![DOI](https://joss.theoj.org/papers/10.21105/joss.02447/status.svg)](https://doi.org/10.21105/joss.02447) -->
<!-- badges: end -->

> Joint Shift Estimation using Stochastic Interventions **Authors:**
> [David McCoy](https://davidmccoy.org)

------------------------------------------------------------------------

## What’s `JointXshift`?

The `JointXshift` is a simple package that compares the expected outcome
change under a joint stochastic shift intervention (the change in
outcome compared to the observed outcome if all exposures were shifted
by a fixed amount), to the additive expected outcome change, the
additive impact of each exposure shifted individually. The idea is that
if the joint shift is positively larger, this indicates synergistic
relationships in the mixed exposure, if it is negative this indicates
more antagonistic relationships exist. Of course, if there is no
difference between the joint and additive this could indicate no
interaction or interactions that are both synergistic and antagonistic
and cancel each other out.

------------------------------------------------------------------------

## Installation

*Note:* Because the `JointXshift` package (currently) depends on `sl3`
that allows ensemble machine learning to be used for nuisance parameter
estimation and `sl3` is not on CRAN the `JointXshift` package is not
available on CRAN and must be downloaded here.

There are many depedencies for `JointXshift` so it’s easier to break up
installation of the various packages to ensure proper installation.

First install the basis estimators used in the data-adaptive variable
discovery of the exposure and covariate space:

`JointXshift` uses the `sl3` package to build ensemble machine learners
for each nuisance parameter.

``` r
remotes::install_github("tlverse/sl3@devel")
```

Make sure `sl3` installs correctly then install `JointXshift`

``` r
remotes::install_github("blind-contours/JointXshift@main")
```

------------------------------------------------------------------------

## Example

To illustrate how `JointXshift` may be used to ascertain the effect of a
mixed exposure, consider the following example:

``` r
library(JointXshift)
library(devtools)
#> Loading required package: usethis
library(kableExtra)
library(sl3)

seed <- 429153
set.seed(seed)
```

We will directly use synthetic data from the NIEHS used to test new
mixture methods. This data has built in strong positive and negative
marginal effects and certain interactions. Found here:
<https://github.com/niehs-prime/2015-NIEHS-MIxtures-Workshop>

``` r
data("NIEHS_data_1", package = "JointXshift")
```

``` r
NIEHS_data_1$W <- rnorm(nrow(NIEHS_data_1), mean = 0, sd = 0.1)
w <- NIEHS_data_1[, c("W", "Z")]
a <- NIEHS_data_1[, c("X1", "X2", "X3", "X4", "X5", "X6", "X7")]
y <- NIEHS_data_1$Y

deltas <- list(
  "X1" = 1, "X2" = 1, "X3" = 1,
  "X4" = 1, "X5" = 1, "X6" = 1, "X7" = 1
)
head(NIEHS_data_1) %>%
  kbl(caption = "NIEHS Data") %>%
  kable_classic(full_width = F, html_font = "Cambria")
```

<table class=" lightable-classic" style="font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
NIEHS Data
</caption>
<thead>
<tr>
<th style="text-align:right;">
obs
</th>
<th style="text-align:right;">
Y
</th>
<th style="text-align:right;">
X1
</th>
<th style="text-align:right;">
X2
</th>
<th style="text-align:right;">
X3
</th>
<th style="text-align:right;">
X4
</th>
<th style="text-align:right;">
X5
</th>
<th style="text-align:right;">
X6
</th>
<th style="text-align:right;">
X7
</th>
<th style="text-align:right;">
Z
</th>
<th style="text-align:right;">
W
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7.534686
</td>
<td style="text-align:right;">
0.4157066
</td>
<td style="text-align:right;">
0.5308077
</td>
<td style="text-align:right;">
0.2223965
</td>
<td style="text-align:right;">
1.1592634
</td>
<td style="text-align:right;">
2.4577556
</td>
<td style="text-align:right;">
0.9438601
</td>
<td style="text-align:right;">
1.8714406
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.1335790
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
19.611934
</td>
<td style="text-align:right;">
0.5293572
</td>
<td style="text-align:right;">
0.9339570
</td>
<td style="text-align:right;">
1.1210595
</td>
<td style="text-align:right;">
1.3350074
</td>
<td style="text-align:right;">
0.3096883
</td>
<td style="text-align:right;">
0.5190970
</td>
<td style="text-align:right;">
0.2418065
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0585291
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
12.664050
</td>
<td style="text-align:right;">
0.4849759
</td>
<td style="text-align:right;">
0.7210988
</td>
<td style="text-align:right;">
0.4629027
</td>
<td style="text-align:right;">
1.0334138
</td>
<td style="text-align:right;">
0.9492810
</td>
<td style="text-align:right;">
0.3664090
</td>
<td style="text-align:right;">
0.3502445
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.1342057
</td>
</tr>
<tr>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
15.600288
</td>
<td style="text-align:right;">
0.8275456
</td>
<td style="text-align:right;">
1.0457137
</td>
<td style="text-align:right;">
0.9699040
</td>
<td style="text-align:right;">
0.9045099
</td>
<td style="text-align:right;">
0.9107914
</td>
<td style="text-align:right;">
0.4299847
</td>
<td style="text-align:right;">
1.0007901
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0734320
</td>
</tr>
<tr>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
18.606498
</td>
<td style="text-align:right;">
0.5190363
</td>
<td style="text-align:right;">
0.7802400
</td>
<td style="text-align:right;">
0.6142188
</td>
<td style="text-align:right;">
0.3729743
</td>
<td style="text-align:right;">
0.5038126
</td>
<td style="text-align:right;">
0.3575472
</td>
<td style="text-align:right;">
0.5906156
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
-0.0148427
</td>
</tr>
<tr>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
18.525890
</td>
<td style="text-align:right;">
0.4009491
</td>
<td style="text-align:right;">
0.8639886
</td>
<td style="text-align:right;">
0.5501847
</td>
<td style="text-align:right;">
0.9011016
</td>
<td style="text-align:right;">
1.2907615
</td>
<td style="text-align:right;">
0.7990418
</td>
<td style="text-align:right;">
1.5097039
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.1749775
</td>
</tr>
</tbody>
</table>

We will apply our method to see if there is largely more synergistic or
antagonistic relationships in the data.

``` r

ptm <- proc.time()
sim_results <- JointXshift(
  w = w,
  a = a,
  y = y,
  delta = deltas,
  n_folds = 3,
  num_cores = 6,
  outcome_type = "continuous",
  seed = seed)
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found
#> Error in get(family, mode = "function", envir = parent.frame()) : 
#>   object 'multinomial' of mode 'function' was not found

proc.time() - ptm
#>    user  system elapsed 
#>  99.657   5.496 141.089

## marginal effects
additive_results <- sim_results$`Additive Effects`
joint_results <- sim_results$`Joint Effects`
joint_vs_additive <- sim_results$`Joint vs Additive Effects`
```

The expected outcome change under additive assumptions is:

``` r
additive_results
#>        psi  psi_var  se_ests CI_lower CI_upper        p_vals
#> 1 54.74664 4.605863 2.146127  50.5403   58.953 1.174941e-305
```

The expected outcome change under a joint shift is:

``` r
joint_results
#>        psi   psi_var   se_ests CI_lower CI_upper       p_vals
#> 1 7.017974 0.1021843 0.3196628   6.3914   7.6445 2.230457e-35
```

The difference is:

``` r
joint_vs_additive
#>         psi  psi_var  se_ests CI_lower CI_upper        p_vals
#> 1 -47.72867 3.373185 1.836623 -51.3284  -44.129 1.046561e-271
```

This indicates that, compared to additive, the expected outcome change
is -47.7, indicating possible more antagonistic relationships in the
mixed exposure compared to synergistic.

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/blind-contours/JointXshift/issues).
Further details on filing issues are provided in our [contribution
guidelines](https://github.com/blind-contours/%20JointXshift/main/contributing.md).

------------------------------------------------------------------------

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/blind-contours/JointXshift/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

------------------------------------------------------------------------

## Citation

After using the `JointXshift` R package, please cite the following:

------------------------------------------------------------------------

## Related

- [R/`tmle3shift`](https://github.com/tlverse/tmle3shift) - An R package
  providing an independent implementation of the same core routines for
  the TML estimation procedure and statistical methodology as is made
  available here, through reliance on a unified interface for Targeted
  Learning provided by the [`tmle3`](https://github.com/tlverse/tmle3)
  engine of the [`tlverse` ecosystem](https://github.com/tlverse).

- [R/`medshift`](https://github.com/nhejazi/medshift) - An R package
  providing facilities to estimate the causal effect of stochastic
  treatment regimes in the mediation setting, including classical (IPW)
  and augmented double robust (one-step) estimators. This is an
  implementation of the methodology explored by Dı́az and Hejazi (2020).

- [R/`haldensify`](https://github.com/nhejazi/haldensify) - A minimal
  package for estimating the conditional density treatment mechanism
  component of this parameter based on using the [highly adaptive
  lasso](https://github.com/tlverse/hal9001) (Coyle et al. 2022; Hejazi,
  Coyle, and van der Laan 2020) in combination with a pooled hazard
  regression. This package implements a variant of the approach
  advocated by Dı́az and van der Laan (2011).

------------------------------------------------------------------------

## Funding

The development of this software was supported in part through NIH grant
P42ES004705 from NIEHS

------------------------------------------------------------------------

## License

© 2020-2022 [David B. McCoy](https://davidmccoy.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    Copyright (c) 2020-2022 David B. McCoy
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

------------------------------------------------------------------------

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-coyle-hal9001-rpkg" class="csl-entry">

Coyle, Jeremy R, Nima S Hejazi, Rachael V Phillips, Lars W van der Laan,
and Mark J van der Laan. 2022. “<span class="nocase">hal9001</span>: The
Scalable Highly Adaptive Lasso.”
<https://doi.org/10.5281/zenodo.3558313>.

</div>

<div id="ref-diaz2020causal" class="csl-entry">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 82 (3): 661–83.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-diaz2011super" class="csl-entry">

Dı́az, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *The International Journal of Biostatistics* 7 (1): 1–20.

</div>

<div id="ref-hejazi2020hal9001-joss" class="csl-entry">

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020.
“<span class="nocase">hal9001</span>: Scalable Highly Adaptive Lasso
Regression in R.” *Journal of Open Source Software* 5 (53): 2526.
<https://doi.org/10.21105/joss.02526>.

</div>

</div>
