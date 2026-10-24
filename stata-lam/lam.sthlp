{smcl}
{* *! version 0.1.0 23oct2026}{...}
{vieweralsosee "lam_plot" "help lam_plot"}{...}
{viewerjumpto "Syntax" "lam##syntax"}{...}
{viewerjumpto "Description" "lam##description"}{...}
{viewerjumpto "Options" "lam##options"}{...}
{viewerjumpto "Postestimation" "lam##postestimation"}{...}
{viewerjumpto "Stored results" "lam##results"}{...}
{viewerjumpto "Examples" "lam##examples"}{...}
{title:Title}

{p 4 4 2}
{bf:lam} {hline 2} Latent attitude trends from pooled ordinal survey series

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:lam} {it:depvar} {ifin} {weight}{cmd:,}
{opt survey(varname)} {opt time(varname)}
[{opt baseperiod(#)} {opt level(#)} {opt iterate(#)} {opt nolog}]

{p 4 4 2}
{it:depvar} must be numeric and ordered so that larger values indicate the same
attitude direction in every survey. Survey identifiers may be numeric or string.
The time variable must be numeric. Only {cmd:fweight}s are accepted.

{marker description}{...}
{title:Description}

{p 4 4 2}
{cmd:lam} implements the constant-cutpoint latent attitude method of Wang, Xie,
and Huang (2023). Each survey has its own ordered cutpoints, held fixed across
periods. All surveys observed in one period share that period's latent mean.
The latent residual distribution is standard normal.

{p 8 8 2}
Pr(Y = k | q,t) = Phi(tau[q,k] - mu[t]) - Phi(tau[q,k-1] - mu[t])

{p 4 4 2}
The first and last category boundaries are negative and positive infinity.
The mean in one observed period is fixed at zero to identify location. Other
means measure differences from that reference period in latent standard
deviation units. They are not approval proportions.

{p 4 4 2}
Surveys must be connected through overlapping periods. Indirect overlap is
sufficient. For example, A can overlap B and B can overlap C even if A and C
never overlap. A common period shared by every survey is unnecessary.
Disconnected groups are rejected because a single reference period cannot
identify their relative locations.

{p 4 4 2}
The method assumes consistent measurement within each series and representative
samples of the same target population in overlapping periods. Connectivity
checks cannot establish these substantive assumptions.

{p 4 4 2}
Missing outcomes, identifiers, periods, and weights are excluded. Zero-frequency
rows are omitted from the estimation sample. Frequency weights must be
nonnegative integer counts. The estimator combines identical survey-period-category
cells internally without changing the likelihood. Frequency-weighted counts
therefore give the same estimates and observed-information variance as expanded
respondent records.

{p 4 4 2}
Each survey must have at least two observed response categories in the estimation
sample. Nonconsecutive numeric response codes are permitted. Categories are
ranked separately within each survey. An entirely absent category is outside
the fitted support and has no separately estimated boundaries. Value labels do
not determine category order.

{p 4 4 2}
This version implements the paper's main model. It does not support covariates,
sampling weights, {cmd:svy}, cluster-robust variance, time-varying thresholds,
or interpolation of unobserved periods. Source code targets Stata 16 or later.
This project was runtime-tested with StataNow 19.5.

{marker options}{...}
{title:Options}

{phang}
{opt survey(varname)} specifies the survey-series identifier. It is required.
Different instruments or instrument versions should be assigned distinct series
identifiers when their measurement differs.

{phang}
{opt time(varname)} specifies the numeric period variable. It is required.
Periods need not be consecutive integers or calendar years.

{phang}
{opt baseperiod(#)} sets the mean of the specified observed period to zero.
The default is the smallest observed period in the estimation sample.
Changing this reference shifts every mean and cutpoint without changing
fitted probabilities or the likelihood.

{phang}
{opt level(#)} sets the confidence level for the returned trend intervals.
The default is Stata's current confidence level.

{phang}
{opt iterate(#)} sets the maximum optimizer iterations. The default is 300.

{phang}
{opt nolog} suppresses the iteration log.

{marker postestimation}{...}
{title:Postestimation}

{p 8 16 2}
{cmd:predict} [{it:type}] {it:newvar} {ifin}
[{cmd:,} {opt xb} | {opt stdp} | {opt pr} {opt outcome(#)}]

{phang}
{opt xb} returns the fitted latent population mean for each record's period.
This is the default. It is not a respondent-specific posterior latent score.
Only periods present in the fitted model can be predicted.

{phang}
{opt stdp} returns the standard error of the fitted period mean. The reference
period has a zero standard error because its mean is fixed by normalization.

{phang}
{opt pr} with {opt outcome(#)} returns the fitted response probability for the
specified original numeric response code, using each record's survey and period.
The response variable need not remain in the data. Unknown surveys, periods,
or category codes yield missing where required by the prediction.

{p 4 4 2}
{help lam_plot} displays the fitted trend and its confidence interval.
Use {cmd:graph export "filename.png", replace} to save a PNG image.

{p 4 4 2}
{cmd:estimates store}, {cmd:estimates save}, and {cmd:lincom} use the posted
coefficient vector and covariance matrix. Free means have names such as
{cmd:mu:t2}, where the number identifies the sorted period index, not the
original period value. The reference-period coefficient is omitted.
Cutpoints have names such as {cmd:cut1:c1}, using survey and cutpoint indices.

{p 4 4 2}
Default variance estimates and confidence intervals use the observed information
matrix. A respondent bootstrap must resample
respondents under a suitable sampling design. Resampling aggregated count rows
does not perform a respondent bootstrap. Bootstrap analyses must also check
that resamples retain the required survey, period, and category support.

{p 4 4 2}
After {cmd:bootstrap} or a manual covariance replacement, predictions and plots
use the current {cmd:e(V)}. The cached {cmd:e(trend)} and {cmd:e(cutpoints)}
matrices are not refreshed by an external prefix. {cmd:lam_plot} returns its
current table in {cmd:r(trend)}. With bootstrap standard errors the plotted
interval is a normal approximation, not a percentile or bias-corrected interval.

{p 4 4 2}
The estimator checks for separated response support and rejects a model without
a finite maximum-likelihood estimate. Connected survey overlap alone does not
guarantee a finite estimate.

{marker results}{...}
{title:Stored results}

{p 4 4 2}
{cmd:lam} stores the following in {cmd:e()}.

{p2colset 5 25 27 2}{...}
{p2col:{cmd:e(N)}}Total positive frequency count{p_end}
{p2col:{cmd:e(N_rows)}}Number of contributing input rows{p_end}
{p2col:{cmd:e(N_cells)}}Number of positive survey-period-category cells{p_end}
{p2col:{cmd:e(n_surveys)}}Number of survey series{p_end}
{p2col:{cmd:e(n_periods)}}Number of observed periods{p_end}
{p2col:{cmd:e(k)}}Number of free parameters{p_end}
{p2col:{cmd:e(ll)}}Log likelihood, without multinomial combinatorial constants{p_end}
{p2col:{cmd:e(baseperiod)}}Reference-period value{p_end}
{p2col:{cmd:e(level)}}Confidence level{p_end}
{p2col:{cmd:e(converged)}}Optimizer convergence indicator{p_end}
{p2col:{cmd:e(iterations)}}Optimizer iterations{p_end}
{p2col:{cmd:e(b)}}Free means and natural-scale cutpoints{p_end}
{p2col:{cmd:e(V)}}Observed-information covariance on the posted scale{p_end}
{p2col:{cmd:e(trend)}}Columns {cmd:time estimate se ll ul}{p_end}
{p2col:{cmd:e(cutpoints)}}Columns {cmd:survey_id category threshold se}{p_end}
{p2col:{cmd:e(categories)}}Columns {cmd:survey_id category rank}{p_end}
{p2col:{cmd:e(periods)}}Sorted observed period values{p_end}
{p2col:{cmd:e(surveyvalues)}}Exact numeric survey identifiers, missing for string surveys{p_end}
{p2col:{cmd:e(survey1)}}Original value of the first survey identifier{p_end}
{p2col:{cmd:e(survey2)}}Original value of the second survey identifier, and so forth{p_end}
{p2col:{cmd:e(timevar)}}Time variable name{p_end}
{p2col:{cmd:e(surveyvar)}}Survey variable name{p_end}
{p2col:{cmd:e(cmd)}}{cmd:lam}{p_end}
{p2col:{cmd:e(sample)}}Indicator for contributing input rows{p_end}
{p2colreset}{...}

{p 4 4 2}
The {cmd:category} column in {cmd:e(cutpoints)} is the original code immediately
below the threshold. Survey indices are consecutive internal identifiers.
The reference-period row in {cmd:e(trend)} has estimate, standard error, and
interval endpoints equal to zero by construction.

{marker examples}{...}
{title:Examples}

{p 4 4 2}
Using the original example variables

{p 8 8 2}{cmd:. lam k, survey(q) time(s)}{p_end}
{p 8 8 2}{cmd:. matrix list e(trend)}{p_end}
{p 8 8 2}{cmd:. predict double attitude, xb}{p_end}
{p 8 8 2}{cmd:. predict double probability, pr outcome(2)}{p_end}
{p 8 8 2}{cmd:. lam_plot}{p_end}

{p 4 4 2}
Using counts and a specified reference year

{p 8 8 2}{cmd:. lam category [fw=count], survey(source) time(year) baseperiod(2000) nolog}{p_end}

{p 4 4 2}
Run {cmd:do examples/synthetic.do} from the source package directory for a
standalone reproducible example. The empirical example
{cmd:do examples/china.do} checks for {cmd:../lam-main/lam-china.dta}, then
{cmd:../lam-china.dta}, supporting both the local workspace and the GitHub
repository layout. The empirical data are not redistributed by this package.

{title:Reference}

{p 4 4 2}
Wang, D., Y. Xie, and J. Huang. 2023. Trend Analysis with Pooled Data from
Different Survey Series. The Latent Attitude Method. Sociological Methodology.
{browse "https://doi.org/10.1177/00811750231193641":doi.org/10.1177/00811750231193641}.
