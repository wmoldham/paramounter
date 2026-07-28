# paramounter 0.2.0

Acts on a field report from an untargeted Orbitrap lipidomics study (438
injections), which found that the parameters measured on high-dynamic-range data were
far more permissive than the instrument warranted. `legacy = TRUE` still reproduces the
original method exactly.

## The distinction behind these changes

The report separated the parameters by what they control, and tested both halves against
xcms with everything else held constant.

* A **tolerance** — `ppm`, `mz_diff`, `peakwidth` — decides how a peak is built. Setting
  it from the largest value measured cannot be undone downstream: a `ppm` of 46 rather
  than 10 found 1.97× the chromatographic peaks and produced ~6% *fewer* features, with
  QC RSD unchanged and per-feature m/z spans nearly doubled.
* A **threshold** — `noise`, `snthresh`, `min_peak_height`, `width_scans` — decides how
  many candidates pass, and downstream filtering removes what does not reproduce.
  Tightening these was measured to be harmful: features with QC RSD below 30% fell from
  5,880 to 3,944.

**Only the tolerances changed.** The thresholds are left on their extremes, on the
report's evidence rather than on principle.

## Behaviour changes

* **Tolerances are drawn at a quantile rather than a maximum.** The new
  `pm_config(tolerance_quantile = )`, default `0.95`, sets `ppm`, `mz_diff` and the upper
  peak-width bound. It is taken from the untrimmed distribution rather than applied on top
  of `trim`, because the two compose: 0.97 then 0.95 gives 0.9215, not 0.95.

* **The wide-peak branch no longer applies by default.** It halved the upper peak-width
  bound and dropped the lower bound to *zero* when the widest peak exceeded 35 seconds and
  a height-to-width ratio exceeded 515. Both parts caused trouble. A lower bound of zero
  asks `CentWave` for peaks of no width, which it cannot distinguish from a spike. And 515
  compares an intensity-per-second quantity against a constant calibrated on Bruker Q-TOF
  counts, so on instruments with larger intensities it is always true. Once the upper bound
  is a quantile it is already robust to a few unusually wide peaks, which is all the branch
  was defending against.

* On the shipped demo data the default now gives `ppm` 20 rather than 30 and `peakwidth`
  `c(5, 41)` rather than `c(0, 28.5)`. Under `legacy = TRUE` both are unchanged, and 11 of
  the 12 published values still reproduce exactly.

## New

* `check_estimates()` reports, for each measured quantity, the parameter that will be
  used, the median of its distribution, and the ratio between them — labelled `tolerance`
  or `threshold` so you can tell a problem from the intended design. Worth running before
  committing to a long peak-picking run.

## Fixes

* `to_xcms(params, sample_groups = )` no longer requires one entry per *measured* file.
  Measuring a few representative injections and processing a much larger experiment is the
  intended workflow, and the check made the argument unusable for it. Nothing downstream
  needed it: `PeakDensityParam` accepts a vector of any length.

# paramounter 0.1.0

First release. The package installs from GitHub only.

The measurement and translation layers are complete and verified against the original
implementation. The version is `0.1.0` rather than `1.0.0` because the interface is new
and may still move; the numbers it produces are the part that is settled.

## Behaviour change

* **The corrections to the original method are now applied by default.**
  `pm_config()`'s `legacy` argument defaults to `FALSE` rather than `TRUE`, and the
  `legacy_scan_count`, `legacy_isolation` and `legacy_first_match` arguments of
  `zoi_features()`, `flag_isolated_zoi()` and `match_zoi_across_files()` follow it.

  This changes results. Peak width in scans is one lower, a peak counts as isolated
  only when both its neighbours are distant, cross-file matching takes the nearest
  candidate rather than the first in the window, and the xcms grouping bandwidth comes
  from the measured drift instead of being fixed at 5. On the shipped demo data the
  visible effect on `to_xcms()` is `prefilter` 3 to 2 and `bw` 5 to 7.3.

  Pass `pm_config(legacy = TRUE)` for the previous behaviour. Reproducing published
  values needs it, and the README's reproduction section sets it explicitly.

## Packaging

* `R CMD check` is clean: 0 errors, 0 warnings, 0 notes. The fixes were all to
  documentation and packaging metadata, not to behaviour — the class-valued property
  defaults of `universal_parameters` are now deferred with `quote()`, which is what
  produced a valid `\usage` section for it.

## Performance

* The measurement loop is about twice as fast: the five shipped demo files went from
  44 s to 21 s end to end, roughly 8.9 s to 4.2 s per file. The measured values are
  bit-for-bit identical, not merely close.

  The gains came from grouping peaks by bin without a factor round-trip, filling the
  per-zone results by index instead of growing them with `c()` and `rbind()`, checking
  for non-finite values with a single `range()` pass, and skipping the per-bin
  `data.frame()` that the loop immediately took apart again.

## Measuring parameters

* `paramounter()` measures universal LC-MS parameters directly from centroided data
  files and returns a `universal_parameters` object. It takes a list of files, a
  `pm_config()` of settings, and a reader, and does a single pass over the data — no
  search over a parameter space against an objective function.

* The quantities measured are mass tolerance (ppm and absolute), noise level, peak width
  (in seconds and in scans), signal-to-noise ratio, peak height, and — given two or more
  files — instrument mass and retention-time shift.

* The public interface is deliberately small: `paramounter()`, the `pm_config()` and
  `universal_parameters` classes, `read_ms_data()`, and the three `to_*()` translators.
  Every step of the measurement chain is internal. They are documented and separately
  tested, but they are implementation, not interface — most take intermediate structures
  that only another internal function produces, and they do not compose into the pipeline
  without reimplementing the two orchestrators around them.

* `pm_config()` collects every tunable setting in one validated object, defaulting to the
  published Paramounter values. `universal_parameters` carries the measured
  distributions, with `summary`, `print` and `plot` methods.

* `read_ms_data()` reads centroided files through `Spectra`, so the package sits on the
  modern xcms/`Spectra` stack rather than legacy MSnbase/xcmsSet.

## Converting to software settings

* `to_xcms()` returns ready-to-use `CentWaveParam`, `PeakDensityParam` and
  `ObiwarpParam` objects.

* `to_msdial()` and `to_mzmine()` return parameter-value tables, optionally written to
  CSV. These are values to enter into those tools yourself; the package does not generate
  importable configuration files, which would tie it to particular tool versions.

* All three translators derive from a single `point_estimates()` step, so a change to a
  conversion rule follows through to every one of them.

## Reproducing the original

* This is a port of the method of Guo and Huan (2022,
  [doi:10.1021/acs.analchem.1c04758](https://doi.org/10.1021/acs.analchem.1c04758)).
  Original authors Jian Guo and Tao Huan are credited as copyright holders.

* Each computational function was verified against a verbatim port of the corresponding
  original script across thousands of synthetic inputs. Those comparisons live in the
  test suite and are the correctness contract.

* On the shipped demo data, `paramounter()` followed by `to_xcms()` reproduces 11 of the
  12 values published in Table S-7 of the paper's supporting information. The twelfth,
  the grouping bin width, differs for a diagnosed reason: a single mis-paired feature
  under the original's first-match rule. See the README.

* `pm_config(legacy = TRUE)` reproduces the original end-to-end, including its known
  quirks. It is the setting to use when comparing against published values.
