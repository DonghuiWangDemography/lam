{smcl}
{* *! version 0.1.0 23oct2026}{...}
{vieweralsosee "lam" "help lam"}{...}
{title:Title}

{p 4 4 2}
{bf:lam_plot} {hline 2} Plot a fitted latent attitude trend

{title:Syntax}

{p 8 16 2}
{cmd:lam_plot} [{cmd:,} {opt level(#)} {it:twoway_options}]

{title:Description}

{p 4 4 2}
{cmd:lam_plot} plots the period means and confidence interval from the current
{help lam} estimation results. The default confidence level is the level used
by the fitted model. {opt level(#)} changes the interval displayed in the plot.

{p 4 4 2}
The horizontal axis uses the original period values. The vertical axis measures
latent mean differences from the reference period, in latent standard deviation
units. A line joins observed periods for display. No unobserved period is
estimated or interpolated by the model.

{p 4 4 2}
The reference-period mean and confidence limits are zero by normalization.
Intervals for other periods use the current {cmd:e(V)}. The default fit uses
observed-information standard errors. After the {cmd:bootstrap} prefix,
intervals use a normal approximation with bootstrap standard errors, rather
than percentile or bias-corrected bootstrap limits.

{p 4 4 2}
The command leaves the dataset and estimation results available. It displays a
Stata graph and does not export a file automatically.

{p 4 4 2}
The plotted table is returned as {cmd:r(trend)}, with columns {cmd:time},
{cmd:estimate}, {cmd:se}, {cmd:ll}, and {cmd:ul}. The original fit's cached
{cmd:e(trend)} is unchanged.

{title:Options}

{phang}
{opt level(#)} specifies the confidence level shown in the graph.

{phang}
{it:twoway_options} customize the title, axis labels, graph name, and other
standard graph settings. See {help twoway_options}.

{title:Examples}

{p 8 8 2}{cmd:. lam response, survey(source) time(year)}{p_end}
{p 8 8 2}{cmd:. lam_plot, title("Latent attitude trend") xtitle("Year")}{p_end}
{p 8 8 2}{cmd:. graph export "trend.png", replace width(1800)}{p_end}
{p 8 8 2}{cmd:. lam_plot, level(90)}{p_end}
