# fireSense_spreadFit (development version)

- The parameter `rep` is renamed `.rep`, a SpaDES-aware parameter like `.studyAreaName`: `SpaDES.project::setupProject()` sets `.globals$.rep` from the experiment's `.rep` (PredictiveEcology/SpaDES.project#190), as `fireSense_spreadPredict` uses it (PredictiveEcology/fireSense_spreadPredict#24). The DEoptim cache names use its value, as before, so cached generations are still found; the event cache key changes once.

# fireSense_spreadFit 1.1.2

- New parameters `runawayEdgeFrac` (`fireSenseUtils::fireSenseRunawayEdgeFrac`, 0.01) and `runawayEdgeMin` (`fireSenseUtils::fireSenseRunawayEdgeMin`, 3L), passed to `fireSenseUtils::runDEoptim()`: a simulated fire is a runaway only when it burns at least `max(runawayEdgeMin, ceiling(runawayEdgeFrac * ring size))` pixels of its buffer's edge ring, not one. Both are also passed to the threshold calibration, and are in `omitArgs` of the DEoptim and calibration `Cache()` calls: the rule only changes how quickly DEoptim moves away from an unlucky draw, so fits cached under the 1-cell rule stay valid. Needs fireSenseUtils >= 0.2.3.9077.

- Fixed: threshold calibration trials were uniform draws over the wide bounds, and most saturate the logistic, which the objective refuses ("Too burny a landscape", "Not spread out enough"); a fold with no usable trial (ELF 14.3 fold 1: 0 of 96) got no threshold and no early stop. Only about 5% of such draws pass those gates for 14.3 fold 1 (103 of 2000). `runSpreadWithoutDEoptim()` now rejection-samples: it keeps the draws that pass `fireSenseUtils::spreadProbGates()` on the first block of years (no fire simulation, the same gate code as the objective) until `iterThresh` are kept or `maxDrawsPerTrial` (new argument, default 100, a constant not a module parameter) x `iterThresh` have been drawn, logs "accepted A of D draws", and runs the accepted ones as before; with none accepted the threshold is `NA` (`Inf`, no early stop). The calibration's cache key includes `drawActivePars()`, `spreadProbGates()` and `spreadProbGateTest()`. This changes the DEoptim cache key. Needs a fireSenseUtils with `spreadProbGates()` (the PR on fireSenseUtils); the `reqdPkgs` floor is to be set once that merges.
- Fixed: the SNLL threshold was calibrated by luck. `runSpreadWithoutDEoptim()` paired each of the `iterThresh` random parameter sets with an independent random threshold and kept the smallest threshold of the trials that passed it, so the result had no relation to the trials' own SNLLs (ELF 6.2.1 fold 1: non-saturated trials scored 201, 212, 226, but the paired thresholds were 17, 209, 144 and none passed; 14.3 fold 1: 533 and 541 against a largest possible threshold of 381). The trials now run with no early stop, and the threshold is `thresholdMargin` (new parameter, default 2) times the best trial's first-block average annual SNLL, the unit the objective compares `thresh` with. Trials that saturate spreadProb ("Too burny", "Not spread out enough") are not usable; with none usable the threshold is `Inf`. `pickThreshold()` is rewritten to this rule; each trial's value and bailed flag are read from the `firstBlockSNLL` and `bailed` fields of `fireSenseUtils::.objfunSpreadFit(returnTerms = TRUE)` (needs the fireSenseUtils with #122), not parsed from its printed output; and the calibration's cache key now includes its code. This changes the DEoptim cache key.
- Fixed: when every threshold-calibration trial failed, `pickThreshold()` returned `NA` and the fit used it: every DEoptim cluster node died with "missing value where TRUE/FALSE needed" (ELF 6.2.1, heldOutFold 1). `estimateSNLLThresholdPostLargeFires()` now turns `NA` (also a cached one) into `Inf`, i.e. no early stop, with a message; `fitSpread()` stops if `thresh` is `NA`.
- `penaliseCapHits` is renamed `penaliseRunaways` (default `TRUE`): fireSenseUtils no longer caps a simulated fire at a size, and a fire that burns any pixel of the outer edge of its own buffer is the runaway that is scored as at least that big. Passed to `fireSenseUtils::runDEoptim()` and to the threshold calibration. Needs the fireSenseUtils release with PredictiveEcology/fireSenseUtils#122; the `reqdPkgs` floor is to be set once that merges. The threshold calibration must use this objective, because the first-block SNLL has a different scale without the cap. This changes the DEoptim cache key.
- The held-out object of `heldOutFold` (`spreadFitHeldOut_<run>_fold<k>.rds`) now also holds `fit`, the fold's parameters as a one-row ledger `sf` object (all `simulateMembers` members, built by the same `spreadFitLedgerRow()` as the `run` event's row), plus `formula` and `link`, so `fireSense_spreadPredict` can predict with the fold's fit.
- Fixed: a held-out fold's fit (`crossValidate`, `heldOutFold`) used the SNLL threshold calibrated on all years. The threshold bounds the SNLL of the two largest fire years, and a fold's two largest years are not the full data's, so for ELFs 4.3 (fold 2) and 5.2.1 (fold 1) no parameter set ever passed it: every evaluation returned the fail value 1e6 for 5000 generations. Each fold now calibrates its own threshold on the years it fits (`estimateSNLLThresholdPostLargeFires(sim, covs)`), unless `SNLL_FS_thresh` is set. This changes the DEoptim cache key of fold fits.
- `fitSpread()` now stops when every member of the final DEoptim population has the fail value. Before, such a fit went on to score or save parameters that were random draws.

# fireSense_spreadFit 1.1.1

- reqdPkgs now lists `dplyr`, `sf`, `withr` and `reproducible`, which the module calls (`Cache`, `CacheGeo`, `asPath`) but did not list. Version 1.1.1.

# fireSense_spreadFit 1.1.0

- Renamed from `fireSense_SpreadFit` to `fireSense_spreadFit` (module naming convention `<model>_<camelCaseComponent>`); projects must rename the module and its `params` key. Version 1.1.0.
- `fitSpread()` no longer calls `termsInDEoptim()`, which printed "logit1, logit2, ..." for
  maxAsymptote, inflectionPoint1 and the other non-formula parameters. `fireSenseUtils::runDEoptim()`
  now prints the real names from `names(lower)`. Requires `fireSenseUtils@development (>= 0.2.3.9066)`.
  Version 1.0.6.9024.
- New parameter `.plotInterval` (default 25): DEoptim generations between the DEoptim progress
  figures, passed to `fireSenseUtils::runDEoptim()` as `plotEvery`; the final figures are always drawn.
  Drawing them after every generation took 8.3 s of each 53 s generation (16% of a fit's wall time)
  with every worker idle. It is left out of the `runDEoptim()` cache key, so changing it does not refit.
  Needs fireSenseUtils >= 0.2.3.9064 and clusters >= 0.0.52 (now from its `development` branch).
- After a fit (`run`) and for each held-out fold (`crossValidate`, `heldOutFold`), two figures
  compare the fit with its data, through `Plots()` under `figurePath(sim)`, when `.plots` asks for
  them (default `NULL`: none, and nothing extra is simulated). `spreadFitObservedVsSimulated_<run>`
  shows, per covariate, the share of pixel-years that burned against the share that burned in
  simulations from the best parameter set (`fireSenseUtils::plotSpreadFitValidation()`); the
  fold's figure uses only its held-out years. `spreadFitResponseCurves_<run>` shows the fitted
  response curves, titled as the model's response (`fireSenseUtils::plotSpreadFitResponse()`).
  Both come from one simulation of `objfunFireReps` replicates (18 s on ELF 5.3.2). Requires
  `fireSenseUtils@development (>= 0.2.3.9063)`. Version 1.0.6.9022.
- Fixed: the held-out simulation (`spreadObjFunArgs()`, R/fitSpread.R) left out `escapeSizeHa`,
  `jumpTries` and `jumpMeanDist`, so held-out years were simulated without the escape rule the fit
  used (the default `escapeSizeHa` is 50 ha). The in-sample simulations already had it.
- New parameter `heldOutFold` (default `NA`, unchanged behaviour). Set to `1` or `2` to run that
  cross-validation fold as its own job: `init` schedules only `spreadFitPrepare`,
  `estimateThreshold` and `crossValidate`, never the full fit or the ledger write, and
  `crossValidate` fits on the other fold's years and scores this fold's held-out years, writing
  `spreadFitHeldOut_<.runName>_fold<heldOutFold>.rds`. A run script stops after `crossValidate`
  (`events = list(.stopAfter = list(fireSense_SpreadFit = "crossValidate"))`), same as mode
  "validate". Lets the two folds of a held-out experiment run as separate jobs instead of one job
  doing both. Version 1.0.6.9021.
- `hillSlope1` (the spread link's slope) is fixed at 1, not fitted by `DEoptim`.
  `estimateSpreadParams()` (fireSense_SpreadFit.R:886-921 pre-fix) put it in the default `upper`/
  `lower` bounds with `[0.2, 2]`, but with the link's linear predictor `x = covariates %*% beta`,
  `hillSlope1` enters only as `hillSlope1 * x`, so scaling every covariate coefficient by `k` and
  dividing `hillSlope1` by `k` leaves every prediction unchanged: it was never identifiable, and let
  every coefficient drift along that 10x ridge. `estimateSpreadParams()` no longer emits
  `hillSlope1`; `fireSenseUtils::.objfunSpreadFit()` (>= 0.2.3.9049) reinserts `hillSlope1 = 1`
  before evaluating the fit, and the `run` event's ledger row gets it back too
  (`addHillSlope1ToLedger()`), so an old ledger row keeps predicting with its own fitted
  `hillSlope1` and a new one predicts with 1. A supplied `upper`/`lower` naming `hillSlope1` is now
  an error. Requires `fireSenseUtils@development (>= 0.2.3.9049)`. Version 1.0.6.9020.
- `spreadFitPrep()` (fireSense_SpreadFit.R:519-526 pre-fix) appended every non-annual covariate
  name to youngAge's own `mutuallyExclusiveCols` entry, including `youngAge` itself when it is a
  non-annual column. `fireSenseUtils::makeMutuallyExclusive()` then zeroed `youngAge` on young
  pixels instead of leaving it at 1, and once zeroed, later columns (e.g. `nfLCC_*`) were left
  un-zeroed too. `youngAge` is now excluded from its own pattern list. Requires
  `fireSenseUtils@development (>= 0.2.3.9048)`, which fixes the same root cause inside
  `makeMutuallyExclusive()`. Version 1.0.6.9019.
- The `iterStep` parameter (fireSense_SpreadFit.R:63, default 25L) is removed; `iterStep` is now
  hard-coded to 1 in `fitSpread()`. `iterStep` is supposed to always be 1: with more than one
  generation per DEoptim call, the `run` event's `vapply(sim$DE, function(D) D$member$bestvalit, ...)`
  and its `numIterations <- length(sim$DE)` both assume one generation per block, so a fit with
  `iterStep = 25` crashed after converging ("values must be length 1, but FUN(X[[1]]) result is length
  25"). Projects used to set `iterStep = 1` themselves; that setting was dropped somewhere along the
  way. Version 1.0.6.9018.
- The `crossValidate` event (fireSense_SpreadFit.R:476-477 pre-fix) put mode "validate"'s result in
  `sim$spreadFitHeldOut` but never wrote it to disk. Batch runs stop after `crossValidate`
  (`events = list(.stopAfter = list(fireSense_SpreadFit = "crossValidate"))`), so the simList is
  discarded and the held-out validation was lost. `crossValidate` now also writes
  `sim$spreadFitHeldOut` to `file.path(outputPath(sim), currentModule(sim),
  "spreadFitHeldOut_<.runName>.rds")`. Version 1.0.6.9017.
- `estimateSpreadParams()` (fireSense_SpreadFit.R:886-889 pre-fix) set the sign of a covariate's
  DEoptim bound by whether its name appeared in the annual covariates table, so a non-drought
  annual covariate (e.g. `PPT_sm`) was wrongly floored at 0 like a drought index, and the default
  bounds (`upperAndLowerVal = 9`, `upperAndLowerValFuel = 60`) were narrow enough to bind: a
  held-out experiment (7 ELFs x 2 folds) found climate estimates up to 25.7, youngAge top-10
  medians down to -23.0, fuel estimates up to 54.5 (29 of 82 above 25), and non-forest classes
  reaching +-9. Sign is now decided by term name: drought-index terms (`CMD` or `MDC` anywhere in
  the name) get a lower bound of 0, `youngAge` gets an upper bound of 0, and every other term,
  including any other annual covariate, is symmetric. `upperAndLowerVal` defaults to 50 and
  `upperAndLowerValFuel` to 100, wide enough that they constrain sign, not magnitude. Version
  1.0.6.9016.
- `runSpreadWithoutDEoptim()` drew its threshold-calibration parameter sets unnamed, so
  `fireSenseUtils:::.objfunSpreadFit` (which tells a trailing `yearSpreadSD` bound apart from a
  logistic parameter only by name) miscounted the logistic parameters and every trial errored;
  `mod$thresh` came back `NA`. Drawn (and unnamed user-supplied) parameter sets are now named with
  `names(lower)`. Version 1.0.6.9015.
- `histOfCovariates()` plotted a hard-coded `CMDsm` column regardless of which annual climate
  covariate the ELF actually used, so any ELF with a different column (e.g. `CMD`, `CMD_sp`,
  `cumMDC`-derived columns) failed inside `spreadFitPrepare` with "object 'CMDsm' not found" as
  soon as the plot was drawn. It now plots every annual covariate column, faceted by covariate and
  year, and no longer uses the deprecated `aes_string()`. Version 1.0.6.9014.
- New parameters for fireSenseUtils >= 0.2.3.9045's objective options, ON by default: `yearAreaWeight = "auto"`
  (annual area burned scored against each year's simulated totals), `areaDistWeight = "auto"` (area-weighted
  size distribution), and `jumpTries = 20`/`jumpMeanDist = 3` (a fire stuck below the escape size may jump to
  nearby burnable land). They reach the fit and both threshold calibrations. This changes the objective, and
  so the cache key, of every fit; set the weights and `jumpTries` to 0 for the previous objective.
- The threshold calibration now uses the fit's objective settings: `weighted`, `sizeLik`, `sizeLikDf`,
  `adWeight` and `link` were not passed, so it ran with `weighted = TRUE` and the "kde" likelihood whatever
  the fit used. Calibrated thresholds change, and so does the `estimateThreshold` cache key.
  `weighted = "sqrt"` no longer errors in the rough threshold estimate. Version 1.0.6.9013.
- New parameter `escapeSizeHa` (default 50): the spread model is fitted to escaped fires only, fires that
  reached that size, and each simulated fire burns that area first before spreading normally. Before, any fire
  over 1 pixel counted, and many simulated fires never left their first pixel. It reaches the fit
  (`runDEoptim()`) and the threshold calibration (`runSpreadWithoutDEoptim()`), so both evaluate the same
  objective. `NULL`/`NA` gives the old fit. Needs fireSenseUtils >= 0.2.3.9044. Version 1.0.6.9012.

## The fit ledger

- `spreadFitFilename` now defaults to `"latest"`. A fit is written to the file named for its fire years and model,
  `fireSenseUtils::spreadFitFilenameFor()` (e.g. `fireSenseParams_1985-2024_linearFuel.rds`; the years are
  fireSense_dataPrepFit's `fireYears`, else those of the annual covariates), and readers find each polygon's most
  recent fit with `fireSenseUtils::latestSpreadFits()`. A named file is used as before.
- New parameter `.studyAreaName` (default `NA`), the name PredictiveEcology modules use for the study area. This module does not use it yet.

## DEoptim defaults

- New defaults, so a project need not set them: `strategy = 6` with `DEoptimControl = list(p = 0.1)`,
  `.c = 0`, `iterDEoptim = 5000`, `objfunFireReps = 50` and `DEoptimTests = c("adTest", "SNLL_FS")`.
- The strategy comes from a settings study on ELF 13.1 (NP 110, 3 seeds per setting, 150 generations, 2026-09-15).
  Strategy 6 with p = 0.1 had the lowest median best value (3810, against 3829-3869 for strategies 1, 2, 3 and 6
  with p = 0.2) and stayed best when each run's best members were re-scored 10 times. This is provisional: one
  ELF, and fits far from converged. The study ran with `c = 0.1` in 5-generation DEoptim calls; `c` cannot take
  effect now (see `.c`), so the default is 0.
- `iterDEoptim` is a ceiling: clusters >= 0.0.46 (now required) stops a fit once the population's median value
  has stopped improving.
- `objfunFireReps = 50` and both tests are what production fits have used; no other combination was tested.
- These change the cache key of any fit that relied on the old defaults.

## DEoptim crossover adaptation

- Requires clusters >= 0.0.42 (was 0.0.41). With `iterStep` > 1, that version runs DEoptim with `c = 0`, because
  DEoptim's F adaptation turns every trial vector into NaN once a call's first generation has no successful trial.
  Without it, this module's defaults (`iterStep = 25`, `.c = 0.5`) could crash a fit on every worker, so projects
  had to set `iterStep = 1`.

## Per-year random effect

- New parameter `yearSpreadSDBounds` (default `c(0, 1)`): the default bounds get `yearSpreadSD` last, and
  `fireSenseUtils::runDEoptim()` (>= 0.2.3.9041) fits it as the sd of a per-year random effect on logit spread
  probability, a seasonal departure: each year draws one eps, so all of a year's fires burn hotter or cooler
  together, which widens the simulated fire-size distribution. `NA` turns it off. If only one bound is supplied,
  the other includes `yearSpreadSD` only if the supplied one does. (Briefly `fireSpreadSDBounds`, per fire.) This
  changes every fit's cache key.
- `covFixedRange` also fixes the scale of CMD, CMDsp and cumMDC (all / 100), the other climate candidates of
  fireSense_dataPrepFit's `spread = "auto"`, so their coefficients compare across ELFs as CMDsm's do.

## Objective and link

- New parameters `sizeLik` (default "t"), `sizeLikDf`, `weighted` (default FALSE) and `adWeight` reach the objective in the fit and in the re-score. Fits used "kde" with a log(size) weight before, only because the module could not ask for anything else; "t" without a weight predicted held-out years best in the 2026-09-21 cross-validation. This changes every fit's cache key.
- New parameter `link`: "logistic3pUpper" adds `upperTail1` (bounds `upperTailBounds`, default c(-1, 1)), which changes only how the spread probability approaches its ceiling (`fireSenseUtils::logistic3pUpper()`). The default stays "logistic3p".

## Diagnostics after every fit

- New event `postFitDiagnostics`, scheduled after `run`. It makes `spreadFitRescore`, `spreadFitIdentifiability` (which covariates are identified in isolation), `spreadFitProfile`, `spreadFitSizes` (observed against uncapped simulated fire sizes), `spreadFitLinkSaturation` and `spreadFitConvergence`. The costly parts run on the fit's workers inside `runDEoptim()`: `profileReps` (default 10) and `simulateMembers` (default 10).
- `mode = "validate"` adds `crossValidate`: two fits, each on every other year, predicting the years it did not see (`spreadFitHeldOut`). It never writes the ledger.

## Climate covariate

- New parameter `covFixedRange` (default `list(CMDsm = c(0, 100))`): covariates rescaled with a fixed range, not the range of the polygon's data. CMDsm is now CMDsm / 100 everywhere. With the data's range, 1 meant a CMDsm of 104 in ELF 5.3.2 and 297 in ELF 13.1, so the coefficient meant something different in each polygon, and one that never gets dry stretched its small range over [0, 1]. Pooled over six ELFs on the absolute scale, fire size is flat below a CMDsm of about 125 and about twice as large above 150; no single polygon's own scale shows that. `covFixedRange = list()` restores the old behaviour. An NA in the data still reaches the NA check.

## Fuel covariates

- Fuel biomass is fitted on the linear scale, divided by a fixed 1e4. It arrives from `fireSense_dataPrepFit` logged (`fireSenseUtils::logMinB()`); `spreadFitPrepare` undoes that with `fireSenseUtils::fuelLogToLinear()` on a copy, and the supplied covariates are not changed. On the log scale the treed pixels of ELF 5.3.2 fell in 16% of the covariate range and 45% of pixels sat on the floor, so the fuel coefficients estimated little more than treed against treeless. Needs fireSenseUtils >= 0.2.3.9029.
- `covMinMax_spread` gives every fuel column `c(0, 1e4)` (`fireSenseUtils::fuelLinearRange`) and no longer the data's shared range. `fireSense_SpreadPredict` recognises a linear fit by that range, so parameters fitted earlier, on the log scale, still predict as they did.
- New parameter `upperAndLowerValFuel` (default 60): the default bound of the fuel coefficients. With 9, as for every other covariate, the fuel coefficient sat on its bound.
- In 36 model-selection fits over six ELFs, one column per fuel class, as here, was best or within 0.5% of the best in every ELF; collapsing fuel to a total cost 9.6% in the most mixed ELF (three replicate runs each, p = 0.007).

# fireSense_SpreadFit 1.0.6

First release from `development` since `master` was last updated (2023-09-06). Full history: https://github.com/PredictiveEcology/fireSense_SpreadFit/compare/cfc4cd7...v1.0.6

## Breaking changes

- Removed input `dataFireSense_SpreadFit` (RasterLayer, RasterStack).
- Removed input `firePoints` (SpatialPointsDataFrame).
- Removed input `firePolys` (list).
- Removed input `flammableRTM` (RasterLayer).
- Removed input `polyCentroids` (list).
- Input `rasterToMatch` is now `SpatRaster` (was `RasterLayer`).
- Input `studyArea` is now `sf` (was `SpatialPolygonDataFrame`).
- Removed output `covMinMax` (data.table).
- Removed parameters: `.plot`, `debugMode`, `fireYears`, `formula`, `minBufferSize`, `parallelMachinesIP`, `useCentroids`.

## New features

- New inputs: `.ELFind`, `.runName`, `fireBufferedListDT`, `fireSense_annualSpreadFitCovariates`, `fireSense_nonAnnualSpreadFitCovariates`, `fireSense_spreadFormula`, `parsKnown`, `spreadFirePoints`, `spreadFitAdditionalColNames`.
- New outputs: `covMinMax_spread`, `fsSpreadFit_hists`, `lociList`, `studyAreaWithSpreadParams`.
- New parameters: `.c`, `.plotSize`, `.plots`, `DEoptimTests`, `SNLL_FS_thresh`, `doObjFunAssertions`, `iterThresh`, `libPathDEoptim`, `mode`, `mutuallyExclusiveCols`, `rep`, `spreadFitFilename`, `spreadFitGoogleDriveFolder`, `stopIfNoPreRunFit`, `upperAndLowerVal`, `useCache_DE`.

## Dependencies

- No longer depends on `fastdigest`, `rgeos`.
- Now depends on `clusters`, `fpCompare`, `ggplot2`, `munsell`, `scales`, `terra`, `tidyr`.

## Testing

- testthat suite and CI (`testthat-module`), including a snapshot of the module's inputs, outputs and parameters in `tests/testthat/test-metadata.R`.
