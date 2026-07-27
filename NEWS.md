# paramounter 0.0.0.9000

First working version. The package is not yet released; it installs from GitHub only.
(R's `NEWS.md` parser only recognises a literal version in the heading, so this says
`0.0.0.9000` rather than the usual "development version".)

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

* Every step of the measurement chain is exported, documented and separately testable, so
  the pipeline can be composed by hand as well as run through `paramounter()`:
  `bin_peaks()`, `assemble_bin_traces()`, `smooth_intensity()`, `estimate_noise()`,
  `find_zoi()`, `flag_isolated_zoi()`, `collect_zoi_masses()`, `zoi_features()`,
  `match_zoi_across_files()`, and `estimate_instrument_shift()`.

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
