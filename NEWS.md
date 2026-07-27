# paramounter (development version)

First working version. The package is not yet released; it installs from GitHub only.

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
