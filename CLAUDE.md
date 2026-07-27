# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A port of HuanLab/Paramounter (Guo et al. 2022, *Anal. Chem.* 94:4260) into a clean,
tested, modern R package. It measures universal LC-MS parameters (mass tolerance, noise,
peak width, S/N, peak height, instrument mass/RT shift) directly from centroided MS data
files, then translates them into settings for xcms (primary), MS-DIAL, and MZmine.
Original authors Jian Guo and Tao Huan are credited as copyright holders in
`DESCRIPTION`.

Distribution is GitHub-only. Object system is S7. The MS stack is modern
xcms/`Spectra`, not legacy MSnbase/xcmsSet.

### Reference material

The verbatim original scripts are the ground truth for reproduction and live in
`reference/` — gitignored and `.Rbuildignore`d, so they are local-only:

- `Paramounter_part1 (V2).R` — ppm-cutoff determination
- `Paramounter_part2 (V2).R` — the measurement loop, aggregation, and point estimates
- `XCMS.R` — a hardcoded demo of the xcms workflow
- `Paramounter User Manual Version 3.0.pdf`

Refetch them from `https://raw.githubusercontent.com/HuanLab-backup/Paramounter/main/`.
The paper PDF is in Zotero storage
(`~/Zotero/storage/I3PYUUUS/2022_guo_paramounter.pdf`); its supporting information
carries the published parameter tables.

### Demo data

`inst/extdata/UrineOriginal1-5.mzXML` (committed, 61 MB) are the upstream demo files
from the same repository's `Demo Data.zip`: five technical replicates of a urine sample
on a **Bruker maXis impact Q-TOF, ESI RP(+), centroided, 7473 scans per file**. Reach
them with `system.file("extdata", package = "paramounter")`.

They are the reproduction target. Table S-7 of the SI publishes the Paramounter → xcms
values for exactly this dataset (its "Urine in Bruker Q-TOF RP(+) DDA mode" row):

| ppm | peakwidth | mzdiff | snthresh | integrate | prefilter | noise | bw | minfrac | mzwid | minsamp | max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 0, 28.5 | -0.01 | 3 | 1 | 3, 298 | 298 | 5 | 0.5 | 0.006 | 1 | 100 |

### Reproduction status (2026-07-26)

`paramounter(files, pm_config(ppm_cutoff = 30))` then `to_xcms()` reproduces **11 of 12**
published values exactly: `peakwidth 0, 28.5`, `snthresh 3`, `prefilter 3, 298`,
`noise 298`, `bw 5`, `mzdiff -0.01`, `integrate 1`, `minfrac 0.5`, `minsamp 1`,
`max 100`, and `ppm 30`. Note this now needs `legacy = TRUE` passed explicitly, since
the default flipped to the corrected behaviour.

Runtime is ~4.2 s per file, ~21 s for all five, measured end to end through
`paramounter()`. (It was ~8.9 s per file before the measurement loop was optimised; an
earlier note here claiming ~19 s per file was stale.)

The one gap is grouping `mzwid`/`binSize`: 0.00708 vs the published 0.006.

**Diagnosed.** A verbatim port of part 2's matching block (lines 259–332) run on the same
clean-ZOI input returns results *bit-identical* to
`match_zoi_across_files(legacy_first_match = TRUE)` — 95 features, `mass_shift`
differences exactly 0. The port is faithful; there is no bug to fix. The ppm cutoff is
not the lever either: sweeping it over 29–30.5 leaves the clean set unchanged at 95
features and `mzwid` at 0.00708.

The cause is the first-match rule itself. Exactly **one feature of 95** is mis-paired:
anchor m/z 481.12942 matches file 5's 481.118316 (11 mDa away) when 481.128420 (1 mDa)
sat in the same ±0.015 Da window. That inflated value changes which values survive the
97% trim, lifting the trimmed maximum from 0.00579 to 0.00708. Correcting the matching
rule *alone* gives 0.00579 → 0.006, exactly the published value.

Note this is **not** what `legacy = FALSE` does: that flag also switches isolation to
both-sides, which rebuilds the clean set and lands back at 0.00708 (with `prefilter` 2
and `bw` 7.3). So no supported setting reproduces all 12 values at once — 11 of 12 under
`legacy = TRUE` is the honest ceiling, and why the twelfth differs is now understood.
Why the authors' own run avoided the mis-pairing is not recoverable: it depends on their
within-window zone ordering, which their published output does not record.

**Why `ppm_cutoff` must be passed explicitly to reproduce.** In the original this is a
*human-entered* number: part 1 plots the ppm distribution with a dashed line at the
`ppm_quantile` position and prints "find the cutoff line ... and run part 2 using the ppm
cutoff"; part 2 then hardcodes it (`ppmCut <- 40`). The published 30 is that dashline
rounded down by eye. Our automatic path computes the dashline faithfully — here ~32.4 —
so `ppm_cutoff = NULL` yields `ppm 33`, not 30. Neither is wrong; the original just has a
human in the loop that the port cannot reproduce.

## Commands

Everything runs through `Rscript -e` from the package root.

```r
devtools::test()                                   # full suite (~13k assertions, a few minutes)
devtools::test(filter = "zoi_features")            # one test file (matches tests/testthat/test-<filter>.R)
devtools::load_all()                               # attach for interactive poking
devtools::document()                               # regenerate NAMESPACE, man/, and DESCRIPTION Collate
devtools::check()                                  # full R CMD check
```

`man/` and `NAMESPACE` are roxygen2 output — never hand-edit them; edit the roxygen
block above the function and re-run `document()`.

`README.md` is knitted from `README.Rmd` — never hand-edit it either. Regenerating it
needs pandoc, which is not installed standalone on this machine; point R at the copy
bundled with quarto first, or `build_readme()` fails with "pandoc version 1.12.3 or
higher is required":

```sh
export RSTUDIO_PANDOC=/Applications/quarto/bin/tools/aarch64
Rscript -e 'devtools::build_readme()'
```

`knitr::knit()` alone runs without pandoc but is **not** a substitute: pandoc rewraps
prose, pads table cells, and emits ` ``` r ` fences, so a knitr-only pass rewrites most
of the file. The README's two live chunks run the real pipeline, so a knit takes a
couple of minutes.

Neither `faahKO` nor `msdata` (both Suggests) is installed locally, so a clean run
reports one skip at `test-read_ms_data.R:22`. Suggested packages are used behind
`skip_if_not_installed()` guards. A green baseline is currently
**13,335 passing, 0 failures, 1 skip**. Bioconductor dependencies
(S7, Spectra, xcms) install via BiocManager; when bioconductor.org returns 504s, the
reliable mirror is `https://ftp.gwdg.de/pub/misc/bioconductor`.

## The verification convention

**This is how work gets done here.** Every non-trivial computational function was
verified by running the new implementation against a **verbatim port of the original
code** across thousands of synthetic inputs, *before* being accepted. When adding or
changing computational code, keep this up: write a verbatim port of the relevant
original, generate synthetic inputs, compare, and only then finalize.

Reproduction is **tolerance-based, not bit-identical** — the new code should land within
per-parameter tolerances. **Intentional divergence is allowed** where the original is
wrong or suboptimal, and is documented each time (in the roxygen block, and behind
`legacy = FALSE` when it changes numbers).

This covers pure logic. The two Bioconductor-facing pieces — reading files via `Spectra`
and constructing xcms S4 objects — cannot be verified this way and are checked by
actually running them.

In practice the verbatim ports live at the top of the corresponding test file (typically
`orig_*` functions); `tests/testthat/test-measure_file.R` is the fullest example. That
comparison, not the assertion count, is the correctness contract.

## Architecture

### The pipeline

`paramounter(files, config, reader)` (`R/paramounter.R`) is the only entry point and is
thin: it validates, calls `reader` once per file, calls `measure_file()` per file, and
hands the results to `aggregate_files()`.

- **`measure_file()`** (`R/measure_file.R`, not exported) is the per-file chain and the
  best single file to read to understand the method. For each file it runs
  `bin_peaks()` → per bin `assemble_bin_traces()` → `smooth_intensity()` →
  `estimate_noise()` → `find_zoi()` → `flag_isolated_zoi()` → per zone
  `collect_zoi_masses()` → `zoi_features()`. Returns per-bin noise values and a data
  frame with one row per measured zone.
- **`aggregate_files()`** (`R/aggregate_files.R`, not exported) pools across files:
  derives the ppm cutoff, filters zones, trims each distribution's outlier tail, and —
  with two or more files — runs `match_zoi_across_files()` → `estimate_instrument_shift()`
  to get the mass- and RT-shift distributions. Returns a `universal_parameters`.

Each step function is separately exported, documented, and unit-tested, so they can be
composed by hand as well as by the orchestrator.

Deliberate division of labour: `measure_file()` does **not** apply the ppm cutoff — every
zone clearing the mass-count bounds is returned so its ppm value can enter the pooled
distribution the cutoff is derived from. The cutoff is applied once, in
`aggregate_files()`.

### Two S7 classes

`S7::new_class()`, registered in `.onLoad` via `S7::methods_register()` (`R/zzz.R`).

- **`pm_config`** (`R/pm_config.R`) — every tunable setting in one validated object,
  defaulting to the published Paramounter values. Threaded through the whole pipeline;
  each step function also takes its own bare arguments so it stays testable without a
  config. Properties are declared `class_any` with a `validator` that wraps a `check_*`
  call in `as_message()`.
- **`universal_parameters`** (`R/universal_parameters.R`) — the result: the
  `distributions` list (one finite numeric vector per quantity in `PM_QUANTITIES`,
  empty when unmeasured), plus `files` and `config`. `summary` is a computed getter, so
  it always tracks `distributions`. `print` and `plot` methods live in the same file.

### Translation layer

`to_xcms()`, `to_msdial()`, `to_mzmine()` all read from `point_estimates()`
(`R/point_estimates.R`), the single place where distributions collapse to point values
(max ppm, min noise, the peak-width bounds with the wide-peak adjustment, …). Change a
conversion rule there and all three translators follow.

Each translator splits into a `*_values()` / `*_table()` internal that is pure numerics
(testable without xcms installed) and a thin exported wrapper that builds the
software-specific objects. Missing distributions produce a named error listing which
quantities are absent.

### Validation

`R/validate.R` holds internal `check_*` helpers that `stop()` on failure
(`check_positive_number`, `check_count`, `check_flag`, `check_scalar`, …). All are
`@noRd`. `as_message()` converts a `check_*` call into the `NULL`-or-string form S7
validators expect — that's the bridge between the two validation styles.

Functions validate their arguments up front and error rather than silently coercing.

## Key design decisions

### The `legacy` flag

There is **one master flag**, `pm_config@legacy` (default `FALSE`), and setting it `TRUE`
reproduces the original end-to-end rather than switching a single behaviour. It drives
four things, each reached through its own named parameter at the call site:

| Effect | Where | `legacy = TRUE` | `legacy = FALSE` |
|---|---|---|---|
| `width_scans` off-by-one | `zoi_features(legacy_scan_count=)` | `+2` | `+1` |
| Clean-ZOI isolation rule | `flag_isolated_zoi(legacy_isolation=)` | forward-only | both-sides |
| Cross-file matching | `match_zoi_across_files(legacy_first_match=)` | first match | nearest match |
| xcms grouping `bw` | `to_xcms()` | hardcoded `5` | data-driven `max(rt_shift)` |

`legacy = FALSE` gives the corrected/improved behaviour throughout. When touching any of
these paths, keep the legacy branch intact — it is what the verbatim-port comparisons
are checked against. New intentional divergences belong behind this flag, not in the
shared path.

Two asymmetries in that table worth knowing:

- The first three effects have their own **exported sub-parameter**, so each can be
  toggled independently of `pm_config`. The fourth does not — `xcms_values()` reads
  `params@config@legacy` directly, so `bw` can only be switched through the config.
- Because the sub-parameters carry their own defaults, they must be kept **in step with**
  `pm_config@legacy`'s default. A mismatch would mean `zoi_features()` called directly
  behaves differently from the same function reached through `paramounter()`.

### Other decisions

- **ppm cutoff**: `ppm_cutoff = NULL` means auto — computed from the pooled ppm
  distribution at `ppm_quantile` (default `0.95`), reproducing part 1 of the original.
  A number sets it manually. Note this is the one place the original has a human in the
  loop, so reproducing published values needs the number passed explicitly — see
  **Reproduction status**.
- **Point estimates live in the translation layer**, not in `universal_parameters`,
  which stays software-agnostic. The ceiling/floor rounding, the peak-width
  `+4`/`+5`/`+7` conditional with its 5%-trimmed-mean height:width ratio, the
  `max(3, ·)` S/N floor, and unit conversions are all in `point_estimates()`, shared by
  the three translators.
- **Resolved original inconsistencies**: the demo `XCMS.R` used `integrate = 2` and
  `snthresh = 1.065`, but part 2's computed output uses `integrate = 1` and floors S/N
  at 3. We follow part 2 — the tool's actual output — not the demo script.
- **Export policy**: internal orchestration (`measure_file`, `aggregate_files`, the
  `*_values`/`*_table` helpers, `point_estimates`) is `@noRd`; the building-block step
  functions and user-facing functions are `@export`.

## S7 gotcha (already handled — don't regress)

Files with top-level **construction** (not just definition) create load-order
dependencies. `universal_parameters` has
`config = new_property(pm_config, default = pm_config())`, which constructs at load time
and invokes validators. So `pm_config.R` carries `@include validate.R` and
`universal_parameters.R` carries `@include pm_config.R` to force the `Collate` order.

Default function *arguments* like `config = pm_config()` are lazy and therefore safe —
only top-level construction needs `@include`.

## Style

- snake_case throughout; roxygen2 >= 7.3.0 (required for S7); testthat 3e;
  `test-<file>.R` mirrors `R/<file>.R`.
- Base-R-leaning implementation — `stats`, `graphics`, `S7`, plus `Spectra`/`xcms` at
  the I/O and translation edges. No tidyverse.
- Multi-argument **calls** get one argument per line (trivial constructors stay inline);
  multi-line function **definitions** get one argument per line.
- `stop()` uses `call. = FALSE` with backtick-quoted argument names. Functions validate
  up front and error rather than silently coercing.
- Roxygen blocks are prose-heavy: they explain *why* a step exists, how it relates to its
  neighbours, and where behaviour matches or deviates from the original. Match that
  depth on new functions.
- Every R file opens with a `# <filename>.R` comment.

## Roadmap

The measurement and translation layers are complete and verified, `DESCRIPTION` and the
README are written, and `README.Rmd` knits the demo-data workflow end to end. What
remains:

1. **pkgdown site + CI** (GitHub Actions), if wanted. Note the constraint before
   starting: `README.Rmd` has two live chunks that run the real pipeline, so any job
   that knits it pulls the 61 MB committed data and the full Bioconductor stack.

**Decided against — do not revisit without new information:**

- **Importable-file templates** for MS-DIAL / MZmine. `to_msdial()` and `to_mzmine()`
  return parameter-value tables (plus optional CSV) and that is the intended output;
  users enter the values themselves. Producing drop-in configs would need the full
  format template for a specific named version of each tool, which is an acquisition
  problem rather than a coding one, and the payoff does not justify pinning the package
  to particular tool versions.
- **A vignette.** The README already walks the demo-data workflow step by step with live
  output; a vignette would duplicate it and double the maintenance surface. If long-form
  docs are wanted, write an article on something the README does *not* cover — the
  `legacy` flag's effects, or the `mzwid` diagnosis as a methods note.
