# Method and legacy-code audit

This audit compares the supplied article, Wang, Xie, and Huang's 2023 *Trend Analysis with Pooled Data from Different Survey Series. The Latent Attitude Method*, with the supplied `lam-ml.do`, `lam-sim.do`, and data files. The original publication is identified by DOI `10.1177/00811750231193641`. The wording above replaces punctuation in the original title to respect the project writing preference.

## Statistical target

The published core model is an ordered probit model with a common population mean for each observed time and separate, time-constant cutpoints for each survey series. A series represents a consistent question wording and response scale, rather than necessarily an entire survey organization. Higher response values must represent more of the same latent construct in every series.

For series q, time t, and ordered category k, the model is

```text
Y* = mu[t] + epsilon
epsilon ~ Normal(0, 1)
P(Y = k | q,t) = Phi(tau[q,k] - mu[t]) - Phi(tau[q,k-1] - mu[t])
tau[q,0] = -infinity
tau[q,K[q]] = +infinity
log L = sum(q,t,k) count[q,t,k] * log P(Y = k | q,t)
mu[reference time] = 0
```

The likelihood omits multinomial combinatorial constants, as does the original individual-record Stata implementation. This convention must remain consistent when comparing log likelihoods and information criteria. Grouping identical survey-time-response records and using their frequency counts gives exactly the same objective and information matrix as expanding the records.

The latent variance is fixed to one, and fixing a reference mean to zero removes the common location indeterminacy. The first observed time is a natural default reference. Changing the reference must shift every mean and every cutpoint by the same amount and leave fitted probabilities and the likelihood unchanged. It is a change of coordinates, not evidence that the baseline attitude is known without uncertainty.

## Identification and assumptions

The survey overlap graph must be connected. Two series are adjacent when they are observed at a common time. A path through other surveys is sufficient, and neither a single common time for all surveys nor direct overlap of every pair is necessary. A disconnected graph has separate unidentified locations and cannot produce a single comparable trend without additional assumptions. Connectivity should be checked on the actual estimation sample after exclusions and zero-frequency removal.

The model assumes that surveys at the same time sample the same latent population distribution and that cutpoints remain constant within each series. The graph check establishes structural connectivity but does not prove these substantive assumptions, guarantee adequate overlap precision, or rule out separation and nonfinite estimates. At least two observed response categories are needed per series. A category absent throughout a series cannot have an independently estimated positive-width interval. Any mapping from observed numeric values to ordered categories must be explicit in stored results and documentation.

Observed times may be irregular. The model estimates a separate mean at each observed time and does not interpolate omitted years. Each fitted probability is conditional on its survey's response scale. Means are relative locations on a latent standard-deviation scale and are not favorability percentages.

## Safe numerical implementation

Cutpoints must be strictly increasing. A first unrestricted cutpoint followed by exponentiated positive gaps enforces this throughout optimization. If internal parameters use log gaps, reported cutpoints and their covariance must be transformed back with the full Jacobian.

Log probabilities should use stable Gaussian tail calculations. Directly subtracting two Gaussian CDF values near one, or computing one minus a CDF, can lose precision. The implementation should reject invalid probability intervals, check convergence, and assess observed-information rank before reporting ordinary standard errors. Nonconvergence must remain visible to callers and bootstrap workflows.

The implementation also checks for separation before optimization. It builds a directed inequality graph with one node for every time mean and every finite cutpoint. Each occupied response interval contributes an edge from its lower cutpoint to its mean and an edge from its mean to its upper cutpoint, omitting infinite endpoints. Both the graph and its transpose must be reachable from the baseline node.

The reason follows directly from the likelihood. Along an unbounded parameter direction, keeping a positive-count category probability away from zero requires the lower-end slope to be no greater than the mean slope and the upper-end slope to be no smaller. If all nodes are strongly connected, those inequalities force a common slope, which is zero after anchoring the baseline. Otherwise, increasing levels on the graph's ordered components produce a direction that expands at least one occupied interval without shrinking any other. Because each mapped interior category appears somewhere in its series, the observed-bound inequalities already imply adjacent-cutpoint ordering. The test therefore catches separation involving several means and cutpoints, in addition to simple periods containing only extreme responses. A positive-definite numerical Hessian alone does not provide this protection.

Model-based observed-information standard errors are a valid first-release choice under the independent-response likelihood. They do not account automatically for complex survey designs, clustering, sampling weights, or uncertainty in reconstructed frequencies. The paper's plotted bootstrap intervals and model-based intervals are different inferential procedures and must be labeled accordingly. An ordinary row bootstrap on aggregated count cells is not a respondent bootstrap.

## What the supplied files establish

The empirical file has 108,036 expanded records, 13 survey series, 40 observed years from 1974 through 2019, and 50 finite cutpoints. A reference-year normalization leaves 39 free means plus 50 cutpoints, or 89 parameters, agreeing with the core model's parameter count in Table A1. The published log likelihood is rounded to -155,574.5 and should be used as a coarse benchmark. Direct numerical agreement with a runnable version of the supplied core likelihood is a stronger implementation check.

The supplied simulation CSV has 100,000 records in total, with 10,000 per year across ten years. The article describes 100,000 per year. The first true latent mean is approximately -0.05790020. Simulation truth comparisons on a zero-baseline scale must subtract this value from both the true means and the true cutpoints.

## Legacy issues to fix or document

- Both scripts contain machine-specific Windows paths. The empirical script installs an unrelated plotting dependency and changes the user's permanent graph scheme.
- The likelihood evaluator and initialization are manually tied to 13 series, selected variable names, and specific response counts. They cannot safely serve as a general command.
- Direct CDF subtraction and unrestricted cutoff parameters can create numerical failures or crossed cutpoints during optimization.
- The empirical plot names the lower and upper confidence-limit columns in reverse order. It extracts model-based intervals from `r(table)`, although the article describes bootstrap intervals.
- The second alternative model calls `ml model lf ml` instead of the previously defined likelihood evaluator. The four threshold specifications selected for its pre/post split also do not match a literal interpretation of changing every survey's cutpoints.
- Allowing every cutoff a free linear time slope while estimating free time means introduces a common slope alias. Adding a common multiple of time minus reference time to all means and all time-varying cutpoints leaves every probability unchanged. Baseline anchoring alone does not fix this additional indeterminacy.
- The simulation uses malformed temporary-file macro references and refers to shortened response variable names before renaming the imported CSV fields. Its initial sample draw has no fixed seed.
- Survey B's code samples 15 percent, while the article states 10 percent. Several comments describe different year spans from the executed loops.
- The simulation recenters true means but compares estimated normalized cutpoints to unnormalized true cutpoints. The same location transformation must be applied to both.
- The supplied example data contain reconstructed individual records from marginal distributions. Their appearance as individual rows does not supply unobserved respondent covariates, design variables, or clustering information.

## First-release scope

A defensible first release should implement the article's constant-cutpoint model for arbitrary connected survey series, provide an explicit reference time, accept individual records and integer frequency counts, preserve the user's dataset, and expose means, cutpoints, standard errors, likelihood, and convergence status through standard Stata results. A reproducible example, a help file, a small simulation, and checks against the original core likelihood are more useful than carrying forward exploratory specifications with unresolved identification.

Time-varying cutpoints, covariates, alternative links, complex survey inference, and formal overlap-equivalence tests should be treated as separate extensions that require their own identification conditions and validation. These boundaries concern the initial software implementation and do not narrow the broader research possibilities discussed in the article.

## Independent implementation review

An independent finite-difference check at an off-optimum seven-parameter, two-survey configuration found a maximum absolute gradient discrepancy of approximately 1.48e-8 and a Hessian discrepancy of approximately 2.08e-8. This checks the analytic normal-interval derivatives, the ordered log-gap transformation, and its second-derivative correction. The covariance transformation uses the full cutpoint Jacobian.

The review supplied regression cases for complete separation, one-sided separation in a nonbaseline period, and a joint mean-cutpoint divergence involving two overlapping surveys. The one-sided case initially returned a finite estimate with an enormous standard error despite an infinite likelihood maximum, motivating the directed-graph check above. Additional checks covered very large numeric survey identifiers, Unicode labels, leading and trailing spaces, embedded quotes, and literal Stata macro delimiters in labels. Survey strings are now stored directly through Mata in `e(survey#)` to bypass repeated Stata macro evaluation. Numeric survey values are preserved as doubles in `e(surveyvalues)`, and probability prediction uses that exact numeric mapping.

After those corrections, fresh Stata runs rejected all supplied separation cases with return code 430 and passed the string-mapping, prediction-probability, large-identifier, and derivative checks. These targeted checks supplement the package's main regression suite and do not establish behavior under unsupported survey designs or alternative model specifications.
