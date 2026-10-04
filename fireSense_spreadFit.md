---
title: "fireSense_spreadFit Manual"
subtitle: "v.1.1.2"
date: "Last updated: 2026-10-01"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
bibliography: citations/references_fireSense_spreadFit.bib
link-citations: true
always_allow_html: true
---

# fireSense_spreadFit Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:fireSense-spreadFit) *fireSense_spreadFit*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Eliot McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut, cre], Tati Micheletti <tati.micheletti@gmail.com> [aut], Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut], Jean Marchal <jean.d.marchal@gmail.com> [aut], Alex M. Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Fit statistical models that can be used to parameterize the fire spread component of simulation models (e.g., fireSense [@Marchal:2017a; @Marchal:2017b; @Marchal:2019]).
This module implement a Pattern Oriented Modelling (POM) approach to derive spread probabilities from final fire sizes. <!-- TODO add citation for POM -->
Spread probabilities can vary between pixels, and thus reflect local heterogeneity in environmental conditions.

The fit is a differential evolution search (`DEoptim`, run by `fireSenseUtils::runDEoptim()` on a cluster built by the `clusters` package).
Each candidate parameter set is scored by simulating the historical fires and comparing simulated with observed fire sizes (`fireSenseUtils::.objfunSpreadFit()`).
The 5 best parameter sets, and the covariate ranges used to rescale the covariates, are written as one row per polygon to a shared "fit ledger" on Google Drive (`spreadFitFilename` in `spreadFitGoogleDriveFolder`), keyed by `.ELFind`.
If the ledger already holds a row for the polygon, the module does nothing unless `refitExisting = TRUE`.
By default (`stopIfNoPreRunFit = TRUE`) the module stops rather than start a fit; set it to `FALSE` to fit.

#### `refitExisting`

`refitExisting = TRUE` forces a fit for this polygon even when the ledger already holds a row for it.
Use it when the fit's INPUTS have changed -- new land cover, new vegetation parameters, a new objective
function -- so the stored row is stale and the polygon must be fitted again.

`refitExisting` **overrides** `stopIfNoPreRunFit`: with `refitExisting = TRUE`, `init` schedules the fit
rather than stopping, whatever `stopIfNoPreRunFit` is set to. Setting `stopIfNoPreRunFit = FALSE` is only
needed when the polygon has no ledger row.

`refitExisting` is intended **for developers** who have access to **at least 40 cores**: it triggers a
full DEoptim run (see `cores` and `nCoresNeeded`), which is not practical on a small machine.

### Module inputs and parameters

The covariate tables, fire buffers, fire points and formula are made by `fireSense_dataPrepFit`.
`fireSense_spreadFormula` must be supplied; `.ELFind` defaults to `.runName`.

Table \@ref(tab:moduleInputs-fireSense-spreadFit) shows the full list of module inputs.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-fireSense-spreadFit)(\#tab:moduleInputs-fireSense-spreadFit)List of (ref:fireSense-spreadFit) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
   <th style="text-align:left;"> sourceURL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> .runName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> Some descriptive, short name for this fitting, e.g., ELF14.1 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .ELFind </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> Identifier of the polygon being fit, e.g. '6.1.1'. This becomes the `polygonID` of the row this module writes to the shared cloud fit ledger (`spreadFitFilename` in `spreadFitGoogleDriveFolder`), which `fireSense_dataPrepFit` matches against the polygon ids carried by `rasterToMatchELF`. It must therefore be the polygon's identity, not a run label: `.runName` encodes the whole scenario (climate period, GCM, SSP, rep) in some projects, and keying the ledger on it writes rows no other run can find and trips dataPrepFit's id match. Defaults to `.runName` for backwards compatibility. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireBufferedListDT </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> list of data.tables with fire id, pixelID, and buffer status </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_annualSpreadFitCovariates </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> table of climate and/or veg covariates, burn status, polyID, and pixelID </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_nonAnnualSpreadFitCovariates </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> table of veg covariates, burn status, polyID, and pixelID </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitAdditionalColNames </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> Names of the list-columns of the ledger row. Reset to `fireSenseUtils::spreadFitAdditionalColNamesTxt` if different. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_spreadFormula </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> a formula that contains the annual and non-annual covariates e.g. `~ 0 + MDC + class2 + class3 + youngAge`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> parsKnown </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> Optional vector of known parameters, e.g., from a previous `DEoptim` run. If this is supplied, then 'mode' will be automatically converted to 'debug' </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> template raster for study area </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFirePoints </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> list of `sf` points, one element per year, of fire ignition locations </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyArea </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Polygon being fit; its geometry and crs go in the ledger row. Defaults to NWT. </td>
   <td style="text-align:left;"> https://drive.google.com/open?id=1LUxoY2-pgkCmmNH5goagBp3IMpj6YrdU </td>
  </tr>
</tbody>
</table>

Summary of user-visible parameters (Table \@ref(tab:moduleParams-fireSense-spreadFit))


<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-fireSense-spreadFit)(\#tab:moduleParams-fireSense-spreadFit)List of (ref:fireSense-spreadFit) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramClass </th>
   <th style="text-align:left;"> default </th>
   <th style="text-align:left;"> min </th>
   <th style="text-align:left;"> max </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> .plots </td>
   <td style="text-align:left;"> characte.... </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Plot types passed to `Plots()`, e.g. 'png' or 'screen'; NULL or NA for none. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 25 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> DEoptim generations between DEoptim progress figures; the final figures are always drawn. Passed to `fireSenseUtils::runDEoptim()` as `plotEvery`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotSize </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> 1600, 2000 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> List specifying height and width of plotting device (in pixels) used to plot DEoptim histograms when `visualizeDEoptim` is TRUE. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .runInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> when to start this module? By default, the start time of the simulation. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .studyAreaName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Human-readable name for the study area used. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical,.... </td>
   <td style="text-align:left;"> init </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cores </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Passed to `cores` in `fireSenseUtils::runDEoptim()`: a number of local cores, or a character vector of machine names, one element per core wanted on that machine. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DEoptimTests </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> adTest, .... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Currently either `'SNLL_FS'` or `'adTest'` or a length 2 character vector of both. Passed to `tests` in `fireSenseUtils::.objfunSpreadFit()`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> doObjFunAssertions </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Passed to `fireSenseUtils::.objfunSpreadFit()`; TRUE runs diagnostics but is slower; FALSE for operational runs </td>
  </tr>
  <tr>
   <td style="text-align:left;"> initialpop </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> A numeric matrix of dimensions `NCOL = length(lower)` and `NROW = NP`. This will be passed into DEoptim through `control$initialpop = P(sim)$initialpop` if it is not NULL </td>
  </tr>
  <tr>
   <td style="text-align:left;"> iterDEoptim </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 5000 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> integer defining the maximum number of iterations allowed (DEoptim optimizer). A ceiling: clusters (&gt;= 0.0.46) stops the fit earlier, once the population's median value has stopped improving. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> iterThresh </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 96 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Number of random parameter sets tried when calibrating `SNLL_FS_thresh`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> thresholdMargin </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> When calibrating `SNLL_FS_thresh`, the threshold is this multiple of the best usable trial's first-block average annual SNLL (trials run with no early stop; a trial that saturates spreadProb is not usable). Must be &gt;= 1. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> libPathDEoptim </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> /home/ru.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Absolute path specifying R package directory location to use when running DEotpim. NOTE: this path must be read/write accessible on ALL machines used for fitting (identified in cores). Therefore, it's best use a directory in your user's `~` directory. If the directory does not exist at this path, will attempt to create it. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lower </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> see `?DEoptim`. Lower limits for the logistic function parameters (lower bound, upper bound, slope, asymmetry) and the statistical model parameters (named in the order they appear in the formula). Do not include `hillSlope1`: it is fixed at 1, not fitted (see `estimateSpreadParams()`); supplying it is an error. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> maxFireSpread </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.28 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. Maximum fire spread average to be passed to the `.objFun`. This puts an upper limit on `spreadProb` during optimization. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> link </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> logistic3p </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The spread link. 'logistic3p', or 'logistic3pUpper': the same curve with Stukel's upper tail, one more parameter `upperTail1` that changes only how the curve approaches its ceiling (`fireSenseUtils::logistic3pUpper()`). Its default bounds are `upperTailBounds`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mode </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> fit </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Options: debug, fit, visualize, validate. Can use multiples. 'debug' runs the objective function with visuals instead of DEoptim; 'fit' runs DEoptim; 'visualize' adds the `debug` and `plot` events after the fit; 'validate' adds `crossValidate`, two more fits, each on half the years, predicting the other half (`sim$spreadFitHeldOut`). Validation never writes the ledger, but does write `sim$spreadFitHeldOut` to `outputPath(sim)`, since a batch run typically stops after `crossValidate` and the simList is discarded. See `heldOutFold` to run a single fold as its own job instead of both folds together. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> heldOutFold </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA (default): unchanged behaviour, governed by `mode`. `1` or `2`: run ONLY that cross-validation fold, as its own job. `init` then schedules `spreadFitPrepare`, `estimateThreshold` and `crossValidate` -- never `run`, so the full fit and the ledger write never happen, and the ledger (`stopIfNoPreRunFit`/`refitExisting`) is not consulted. `crossValidate` fits on the OTHER fold's years and scores this fold's held-out years (`cvFolds()` in `R/fitSpread.R`), and writes `spreadFitHeldOut_&lt;.runName&gt;_fold&lt;heldOutFold&gt;.rds` instead of `spreadFitHeldOut_&lt;.runName&gt;.rds`. A run script stops after `crossValidate`: `events = list(.stopAfter = list(fireSense_spreadFit = "crossValidate"))`. Any other value is an error. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> profileReps </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> After the fit, each covariate coefficient in turn is set to 0 and to 5 values across the final population, the others held at the best member, and each point is evaluated this many times (`fireSenseUtils::profileCoefficients()`). About `6 nCoefficients profileReps` evaluations, on the fit's workers. It decides which coefficients are identified in isolation (`sim$spreadFitIdentifiability`). 0 skips it. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mutuallyExclusiveCols </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> c("class.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> a named list of mutually exclusive covariates - see `fireSenseUtils::makeMutuallyExclusive` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nCoresNeeded </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> How many workers to request for the DEoptim cluster. This IS the population size: `clusters::clusterSetup()` sets NP to the workers it builds. `NULL` leaves `fireSenseUtils::runDEoptim()`'s default of 10 per estimated parameter. A generation costs the slowest of NP evaluations and that barely falls as NP falls, so a smaller NP buys throughput by allowing more fits at once rather than by shortening generations (measured 2026-09-16). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> simulateMembers </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> After the fit, this many best members simulate the observed fires, for `sim$spreadFitSizes` and `sim$spreadFitLinkSaturation`; also the members each `crossValidate` fold predicts with. 0 skips it after the fit. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sizeLik </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> t </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Likelihood of fire size in the objective, 'kde' or 't', passed to `fireSenseUtils::runDEoptim()`. 't' with `weighted = FALSE` predicted held-out years best in the 2026-09-21 cross-validation. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> escapeSizeHa </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 50 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Size (ha) a fire must reach to count as escaped. The spread model is fitted to escaped fires only, and each simulated fire first grows to this size with its own spread probabilities (its burning cells stay active until it gets there), then spreads normally. `NULL` or `NA` gives the old fit (any fire over 1 pixel). Passed to `fireSenseUtils::runDEoptim()`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sizeLikDf </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Degrees of freedom of the 't' size likelihood. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> weighted </td>
   <td style="text-align:left;"> logical&amp;#124;.... </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Weight of each fire in the size likelihood: FALSE (none), TRUE (log size) or 'sqrt'. Passed to `fireSenseUtils::runDEoptim()`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> adWeight </td>
   <td style="text-align:left;"> characte.... </td>
   <td style="text-align:left;"> auto </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Weight of the Anderson-Darling term against the size likelihood; 'auto' is `fireSenseUtils::adWeightAuto()`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> yearAreaWeight </td>
   <td style="text-align:left;"> numeric&amp;#124;.... </td>
   <td style="text-align:left;"> auto </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Weight of the annual-area term in the objective: each fit year's observed area burned is scored against that year's simulated totals (one per replicate) with the size likelihood. 0 leaves it out; 'auto' (default) is (number of fitted fires) / (number of fit years), so the year view and the per-fire view weigh the same. Passed to `fireSenseUtils::runDEoptim()`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> areaDistWeight </td>
   <td style="text-align:left;"> numeric&amp;#124;.... </td>
   <td style="text-align:left;"> auto </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Weight of the area-weighted size-distribution term: simulated and observed fires compared by the share of area burned that fires up to each size make up (`fireSenseUtils::areaWeightedCvM()`). 0 leaves it out; 'auto' (default) uses the Anderson-Darling term's weight (`fireSenseUtils::adWeightAuto()`). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> penaliseRunaways </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> A simulated fire that burns any pixel of the outer edge of its own buffer is scored as a runaway (at least that big), not as a fire of the size it reached: in the size likelihood it has no density at the observed size, and in the Anderson-Darling and annual-area terms its size is the landscape's pixel count. Fires are not capped at a size; spread is bounded by the buffers. FALSE scores the simulated size. Passed to `fireSenseUtils::runDEoptim()`; the threshold calibration uses the same setting. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> runawayEdgeFrac </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> A simulated fire is a runaway when it burns at least `max(runawayEdgeMin, ceiling(runawayEdgeFrac n))` of the `n` pixels of the edge ring of its own buffer (never more than `n`); one touched pixel is luck. Passed to `fireSenseUtils::runDEoptim()`. Not part of the DEoptim cache key. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> runawayEdgeMin </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The least number of edge-ring pixels that makes a fire a runaway; see `runawayEdgeFrac`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jumpTries </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 20 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> With `escapeSizeHa`: how many attempts a simulated fire that is still below the escape size, with no burnable neighbour left, may make to jump to burnable land nearby (`SpaDES.tools::spreadCpp()`). Default 20; 0 is off. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jumpMeanDist </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Mean jump distance (pixels) for `jumpTries`; distances are exponential, truncated to 1.5-20 pixels. No effect while `jumpTries` is 0. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> objFunCoresInternal </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Integer defining the number of cores to pass to `mcmapply(mc.cores = ...)` This will fork this many to do the years loop internally. This would be in addition to `cores` and is effecively a multiplier. The computer needs to have `cores objFunCoresInternal` threads or it will stall. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> objfunFireReps </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 50 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> integer defining the number of replicates the objective function will attempt each fire. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .rep </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> An optional integer indicating which replicate run this represents. This is used to identify unique runs of `runDEoptim`, from a Cache perspective. For example, if this module is run twice with all the same data, Cache will think that the second run should recover the cache result, unless this `.rep` is modified </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .c </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> the `c` argument passed to DEoptim.control. `iterStep` is hard-coded to 1, so DEoptim's adaptation restarts every generation, and `c` has no effect either way. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DEoptimControl </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Further `DEoptim.control()` settings, e.g. `list(CR = 0.7, F = 0.6)`, passed through `fireSenseUtils::runDEoptim()` to DEoptim. Names must be `DEoptim.control()` arguments. `strategy`, `trace`, `initialpop` and `.c` have their own parameters; `NP` is the number of workers the cluster gets. The default `p = 0.1` is for `strategy = 6`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rescaleAll </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> rescale covariates for `DEOptim` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitGoogleDriveFolder </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> https://.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Google Drive folder url holding the shared fit ledger (`spreadFitFilename`). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitFilename </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> latest </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> File name of the shared fit ledger: an `sf` object with one row of fitted parameters per polygon. `"latest"` (the default) writes to the file named for this fit's fire years and model, `fireSenseUtils::spreadFitFilenameFor()`, e.g. `fireSenseParams_1985-2024_linearFuel.rds`; readers then find it with `fireSenseUtils::latestSpreadFits()`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> strategy </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Passed to `DEoptim.control`. 6 (DE/current-to-p-best/1) with `p = 0.1` did best in a settings study (see NEWS). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SNLL_FS_thresh </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Threshold multiplier used in objective function SNLL fire size test. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> refitExisting </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> FOR DEVELOPERS ONLY: a re-fit is a full DEoptim run and is only practical with access to at least 40 cores. Fit this polygon even when the ledger already holds parameters for it. A ledger row normally means the fit is done, and the run event skips it. Set this when the fit's INPUTS have changed -- new land cover, new vegetation parameters, a new objective -- so the stored row is stale and the polygon must be fitted again. When TRUE it OVERRIDES `stopIfNoPreRunFit`: `init` schedules the fit instead of stopping. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stopIfNoPreRunFit </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> If TRUE, `init` stops with an error when this polygon would have to be fitted, instead of fitting it. Ignored when `refitExisting` is TRUE. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> trace </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> non-negative integer. If &gt; 0, tracing information on the progress of the optimization are printed every `trace` iteration. Default is 1, i.e. every iteration. Setting to 0 turns off tracing. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> upper </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> see `?DEoptim`. Upper limits for the logistic function parameters (lower bound, upper bound, slope, asymmetry) and the statistical model parameters (named in the order they appear in the formula). Do not include `hillSlope1`: it is fixed at 1, not fitted (see `estimateSpreadParams()`); supplying it is an error. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> useCache_DE </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> should `DEoptim` use `Cache`? to do multiple independent runs, use FALSE </td>
  </tr>
  <tr>
   <td style="text-align:left;"> verbose </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. With increasing number, more verbosity. Level 1 is normal reproducible (e.g., Cache), level 2 includes objective function e.g., print median of spreadProb during calculations </td>
  </tr>
  <tr>
   <td style="text-align:left;"> visualizeDEoptim </td>
   <td style="text-align:left;"> Path </td>
   <td style="text-align:left;"> /tmp/Rtm.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Directory where `runDEoptim` saves parameter plots every `.plotInterval` generations. Reset to `figurePath(sim)` unless its last folder is the module name. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> covFixedRange </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> c(0, 100.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Named list of `c(min, max)`: covariates rescaled with this FIXED range and not with the range of this polygon's data. `CMDsm = c(0, 100)` makes the covariate CMDsm / 100 in every polygon. With the data's range, 1 meant a CMDsm of 104 in one polygon and 297 in another, so the coefficient could not be compared across polygons, and a polygon that never gets dry stretched its small range over [0, 1]. Names not among the covariates are ignored. `fireSense_spreadPredict` rescales with the stored `covMinMax_spread`, so it follows. CMD, CMDsp and cumMDC (also mm) are the other candidates of fireSense_dataPrepFit's `spread = 'auto'`, so an ELF that picks one of them gets the same fixed scale. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> yearSpreadSDBounds </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0, 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Bounds of `yearSpreadSD`, the sd of a per-year random effect on logit spread probability (`fireSenseUtils::.objfunSpreadFit()`), when `lower`/`upper` are not supplied. A seasonal departure: each year draws one eps, so all of a year's fires burn hotter or cooler together, which widens the simulated fire-size distribution. `NA` turns it off. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> upperTailBounds </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> -1, 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Bounds of `upperTail1` when `link` is 'logistic3pUpper' and `lower`/`upper` are not supplied. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> upperAndLowerVal </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 50 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Bound given to each covariate coefficient (`upper` = this, `lower` = minus this) when `upper` or `lower` is not supplied. Bounds should be wide enough that they do not influence the fitted value; only the sign of drought-index and `youngAge` terms is constrained (see `estimateSpreadParams()`). A held-out experiment (7 ELFs x 2 folds) found estimates up to 25.7 and youngAge medians down to -23.0 with the previous default of 9. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> upperAndLowerValFuel </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 100 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> As `upperAndLowerVal`, for the fuel biomass covariates. They are biomass / 1e4, so their coefficients are larger than those of covariates rescaled to [0, 1]. The same held-out experiment found fuel estimates up to 54.5 (29 of 82 above 25) with the previous default of 60. </td>
  </tr>
</tbody>
</table>

### Events

- `init`: schedules `spreadFitPrepare` and, if the polygon needs a fit, `estimateThreshold` then `run` (or `debug` when `mode` includes `"debug"`; `debug` and `plot` after `run` when it includes `"visualize"`);
- `spreadFitPrepare`: sets default `lower`/`upper`, computes `covMinMax_spread` and `lociList`, converts covariates to integers (x 1000);
- `estimateThreshold`: uses `SNLL_FS_thresh`, or calibrates it from `iterThresh` random parameter sets (cached, with a seed derived from `.ELFind`);
- `run`: runs DEoptim (cached per generation) and writes the polygon's row to the ledger;
- `debug`: evaluates the objective function without DEoptim;
- `plot`: histograms of the final population;

### Plotting

With `.plots` set, histograms of the annual and non-annual covariates.
During the fit, `runDEoptim()` saves parameter histograms and trace plots to `visualizeDEoptim` after each generation.

### Saving

The fit's row in the cloud ledger, and the DEoptim cache entries. Nothing else is saved.

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-fireSense-spreadFit)).

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-fireSense-spreadFit)(\#tab:moduleOutputs-fireSense-spreadFit)List of (ref:fireSense-spreadFit) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> covMinMax_spread </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> `data.table` of covariates min and max </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DE </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> list of `DEoptim` objects, one per generation, ordered by best objective value </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaWithSpreadParams </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Rows of the shared fit ledger that intersect `studyArea`, including the row this fit writes: `studyArea` geometry, `polygonID`, and list-columns named by `spreadFitAdditionalColNames` (the 5 best parameter sets are in `params`). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fsSpreadFit_hists </td>
   <td style="text-align:left;"> ggplot </td>
   <td style="text-align:left;"> histograms of each parameter used in `DEoptim` fitting. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lociList </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> per-year `data.table`s of fire start cells and sizes, from `fireSenseUtils::makeLociList()` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitConvergence </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> The objective across the fit's generations (`fireSenseUtils::fitConvergence()`). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitRescore </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> The final population, one row per member, with the mean and sd of its replicated re-scores (`reMean`, `reSD`). The ledger's parameter sets are the best of these. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitIdentifiability </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> One row per covariate coefficient: how tightly the population pins it (`fireSenseUtils::coefIdentifiability()`) and, with `profileReps &gt; 0`, whether dropping it worsens the fit; `identified` = identified in isolation (`fireSenseUtils::identifiedInIsolation()`). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitProfile </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> The one-at-a-time profile around the best member (`fireSenseUtils::profileCoefficients()`). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitSizes </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Observed against simulated fire sizes of the fitted years (`fireSenseUtils::scoreFireSizes()`): bias, error, quantiles. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitLinkSaturation </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Per member, the share of pixel-years at the spread-probability ceiling and the quantiles of spread probability (`fireSenseUtils::linkSaturation()`). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitHeldOut </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> mode 'validate', or `heldOutFold` in `1:2`: `sims`, the held-out years simulated from the fit to the other years (column `fold`), and `score`, from `fireSenseUtils::scoreFireSizes()`. With `heldOutFold`, `sims` holds only that fold, and the list also has `fit` (the fold's fitted parameters as a one-row ledger `sf` object, the same columns the `run` event writes, with all `simulateMembers` members in `params`), `heldOutFold`, `fitYears`, `heldOutYears`, `formula` (`fireSense_spreadFormula`) and `link`, so `fireSense_spreadPredict` can predict with the fold's fit. Also written to `file.path(outputPath(sim), currentModule(sim), "spreadFitHeldOut_&lt;.runName&gt;.rds")` (mode 'validate') or `"...spreadFitHeldOut_&lt;.runName&gt;_fold&lt;heldOutFold&gt;.rds"` (`heldOutFold`). </td>
  </tr>
</tbody>
</table>

### Links to other modules

Inputs come from `fireSense_dataPrepFit`. `fireSense_spreadPredict` uses `studyAreaWithSpreadParams` (the ledger rows) to predict spread probabilities.

### Usage


``` r
## in a project that also runs fireSense_dataPrepFit
params <- list(fireSense_spreadFit = list(
  stopIfNoPreRunFit = FALSE,          # allow a fit to start
  cores = rep(c("hostA", "hostB"), each = 20), # host name repeated once per worker; or a number for localhost
  nCoresNeeded = 40,                  # = DEoptim population size (NP)
  iterDEoptim = 500,
  mode = "fit"
))
objects <- list(.ELFind = "6.1.1")    # polygon id used as the ledger key
```

### Getting help

- <https://github.com/PredictiveEcology/fireSense_spreadFit/issues>

## References

<!-- autogenerated from bibligraphy -->
