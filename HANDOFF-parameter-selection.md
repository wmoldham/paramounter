# Handoff — revisiting how point estimates are drawn from the distributions

## What prompted this

Applying `paramounter` 0.1.0 to an untargeted Orbitrap lipidomics dataset (438 Thermo
`.raw` injections, plasma, 40 min gradient, centroided mzML) produced a `CentWaveParam`
that was far more permissive than the instrument warrants. The measured *distributions*
look right; it is the step from distribution to point estimate that does not transfer to
this data.

This document records the evidence and proposes changes. It does not assume the current
behaviour is wrong in general — the argument is that it is right for tight distributions
and degrades for heavy-tailed ones, and that the package can detect and handle both.

## The observation

Every estimate in `point_estimates()` (`R/point_estimates.R`) is a distribution extreme:

| estimate | rule |
|---|---|
| `max_ppm` | `ceiling(max(ppm))` |
| `max_mz_diff` | `ceiling(max(mz_diff) * 100) / 100` |
| `min_noise` | `floor(min(noise))` |
| `min_peak_scan` | `floor(min(width_scans))` |
| `min_peak_height` | `floor(min(height))` |
| `min_sn` | `max(3, min(sn))` |
| `peakwidth` | `c(ceiling(min(width_seconds)) + 4, ceiling(max(width_seconds)) + 5)` |
| `max_mass_shift`, `max_rt_shift` | `max(...)` |

That is a coherent design: be inclusive at detection so no real feature is lost, and let
downstream QC remove the noise. It works when the measured distribution is tight, because
then the extreme sits close to the bulk. It stops working when the distribution has a
heavy tail, because the extreme is then set by the *worst-measured* zones in the file
rather than by the instrument's capability.

## Evidence

Measured on 5 pooled-QC injections, negative mode (~600,000 zones):

| quantity | min | q25 | median | q75 | max | estimate used |
|---|---|---|---|---|---|---|
| `ppm` | 0 | — | 2.03 | 3.39 | 45.4 | **46** |
| `noise` | 339 | 3,737 | 7,344 | 14,833 | 72,713 | **339** |
| `sn` | 3.7 | 11 | 20.2 | 36 | 772,426 | **3.69** |
| `height` | 3,188 | 6,887 | 10,429 | 17,878 | 1.1e8 | **339** (prefilter) |

The extreme is 22× the median for `ppm` and about 1/20th the median for `noise`.

### The tail is low-abundance zones, not instrument behaviour

Per-zone measurements from one QC file (192,409 zones), binned by peak height:

| height decile | median height | median ppm | 95th pct ppm | % over 20 ppm |
|---|---|---|---|---|
| 1 (lowest) | 3,708 | 2.88 | **135.2** | 24.0 |
| 5 | 8,494 | 2.48 | 53.6 | 8.6 |
| 10 (highest) | 56,442 | 1.63 | **7.3** | 2.8 |

The *median* ppm is essentially flat across the whole intensity range (2.9 → 1.6). Only
the *tail* moves, and it collapses by a factor of 18 from the lowest to the highest
decile. This is the expected behaviour of a trapping mass analyser: at low ion counts the
m/z estimate is dominated by counting statistics, so `sd(masses)` within a zone reflects
measurement noise rather than the instrument's mass accuracy. Spearman ρ between
`log10(height)` and `log10(ppm)` is only −0.145 precisely because the effect lives in the
tail, not the centre.

Zones flagged `isolated` show the same thing from another angle: 2.7% exceed 20 ppm versus
10.3% for the rest.

### Downstream consequence, measured

Six pooled QCs through xcms 4.8 with everything held constant except `ppm`:

| | ppm = 10 | ppm = 46 (as estimated) |
|---|---|---|
| chromatographic peaks | 135,279 | 266,066 (**1.97×**) |
| features after grouping | **15,138** | 14,216 |
| features complete across all 6 | 86.6% | 86.9% |
| QC RSD (p25 / median / p75) | 0.22 / 0.36 / 0.49 | 0.22 / 0.36 / 0.49 |
| per-feature m/z span, p99 | **4.6 ppm** | 8.2 ppm |
| runtime | 7.9 min | 8.9 min |

The estimated value detects **twice the chromatographic peaks and yields ~6% fewer
features**, with RSD quartiles identical to two decimals. The additional detections do not
survive grouping; they cost compute and widen features without adding information.

## The distinction that matters: tolerances vs thresholds

An obvious follow-up hypothesis was that the same problem afflicts the other
extreme-based estimates — `noise` at 339 against a median of 7,344, `snthresh` at 3.69
against a median of 20. **Testing it showed the opposite, and the result is the most
useful thing in this document.**

Same 6 QCs, `ppm` fixed at 10, varying only the intensity and S/N cutoffs:

| arm | peaks | features | **features with QC RSD < 30%** | median RSD | complete |
|---|---|---|---|---|---|
| measured (noise 339, sn 3.69) | 135,279 | 15,138 | **5,880** | 0.36 | 86.6% |
| noise at q25 (3,736) | 124,453 | 14,667 | 5,778 | 0.36 | 86.5% |
| snthresh 10 | 72,402 | 8,784 | 3,944 | 0.33 | 89.6% |
| both, prefilter 3/6,886 | 65,781 | 8,474 | 3,834 | 0.33 | 89.7% |

Tightening `snthresh` improves every ratio — median RSD, the RSD<30% *fraction*,
completeness — while eliminating a third of the features that actually reproduce across
QC injections. The apparent quality gain is an artefact of discarding marginal features;
the count of genuinely reproducible features falls from 5,880 to 3,944. `noise` is close
to inert: an 11-fold increase costs 8% of peaks and 3% of features.

The permissive thresholds are therefore **correct**, and for the reason the design
intends: over-detection is cheap and reversible, because downstream QC and blank filtering
remove what does not reproduce, whereas a feature never detected cannot be recovered.

The difference between the two cases is what the parameter controls:

- **Tolerance parameters** (`ppm`, `mz_diff`, and by the same structure `peakwidth`)
  govern *how a peak is constructed*. Setting them from the extreme degrades the peaks
  themselves — wider m/z spans, doubled peak counts, fewer features. The damage is not
  reversible downstream.
- **Threshold parameters** (`noise`, `snthresh`, `min_peak_height`) govern *how many
  candidates pass*. Setting them from the extreme only adds candidates, which downstream
  filtering removes. The permissiveness is reversible and therefore harmless.

**Any revision should apply to the tolerance parameters only.** Left alone, the threshold
estimates are doing exactly what they should.

## Proposed changes

Three options, roughly in increasing order of invasiveness. They are not exclusive. All
of them should be scoped to the tolerance parameters, per the section above.

### 1. Report tail heaviness (non-breaking, do first)

Add a column to the `@summary` data frame — or a small `check_estimates()` helper —
reporting the ratio of the estimate to the distribution's central value, e.g.
`max/median` for max-like quantities and `median/min` for min-like ones. Optionally warn
when the ratio exceeds some threshold.

This changes no behaviour and would have surfaced the problem immediately: a user seeing
`ppm: max/median = 22.4` knows the estimate is not representative before they run xcms for
an hour. Even if nothing else changes, this is worth having.

### 2. Configurable quantile estimators

Replace the hard `min`/`max` with a configurable position, defaulting to current
behaviour so nothing changes for existing users:

```r
pm_config(estimator_quantile = 1)      # current: max for max-like, min for min-like
pm_config(estimator_quantile = 0.95)   # q95 for max-like, q05 for min-like
```

`point_estimates()` gains the quantile as an argument and applies it symmetrically. This
works on the marginal distributions alone, so it needs no architectural change, and
`legacy = TRUE` would pin it to 1 for reproduction of the published values.

The weakness is that a quantile is a blunt instrument: it discards a fixed fraction
regardless of whether the tail is real. On this data q95 gives ppm ≈ 9.1, which is about
right, but that is partly luck.

### 3. Condition the estimate on intensity (the principled fix)

The diagnostic says the tail is *caused by* low ion counts. The targeted fix is therefore
to estimate mass tolerance from zones with enough signal to measure it, e.g. zones above
the median height, rather than from all zones:

```r
pm_config(estimate_min_height_quantile = 0.25)
```

**This needs an architectural change.** `point_estimates()` currently receives only the
marginal distributions (`@distributions`), which have lost the per-zone linkage — you
cannot condition `ppm` on `height` from marginals. `aggregate_files()` already has the
joined table (`kept`, the per-zone data frame in `R/aggregate_files.R`), so the change is
to retain it, e.g. as a `@zones` property on `universal_parameters`, and have
`point_estimates()` take that instead.

Retaining the per-zone table has independent value: it makes exactly the diagnostic above
(ppm by intensity decile) a one-liner for any user, and it would let `plot()` show
conditional structure rather than only marginals. The cost is object size — roughly
200,000 rows × 9 columns per file — which may warrant storing a subsample or making
retention optional.

My suggestion is to implement 1 immediately, 2 as the general-purpose knob, and treat 3 as
the real fix if the per-zone table is worth keeping for other reasons.

## Compatibility notes

- `legacy = TRUE` currently means "reproduce the published method". Whatever is added
  should be pinned under that flag so the Table S-7 reproduction stays exact.
- Tests in `tests/testthat/` that assert specific estimate values will need updating if
  defaults change. Keeping the default at current behaviour avoids this entirely.
- `NEWS.md` should distinguish "new option, default unchanged" from any default change.

## The wide-peak branch produces a zero lower bound on `peakwidth`

Running the two polarities separately surfaced a second, independent issue.

```r
peakwidth <- if (hi > 35 && ratio > 515) {
  c(0, (ceiling(hi) + 7) / 2)
} else {
  c(ceiling(lo) + 4, ceiling(hi) + 5)
}
```

Positive mode triggered the wide-peak branch and received `peakwidth = c(0, 22.5)`.

**Problem 1 — the lower bound of 0 is degenerate.** `CentWave` converts `peakwidth` to
wavelet scales as roughly `round(peakwidth / mean(diff(scantime)) / 2)`. At this scan rate
(40 min, ~6,200 MS1 scans) that gives `scalerange = c(0, 29)` for positive mode against
`c(6, 22)` for negative. A scale of zero is not a meaningful wavelet width, and a peak
narrower than two or three scans cannot be distinguished from a spike. Negative mode,
which took the normal branch, got a sensible `c(5, 17)`.

**Problem 2 — the `ratio > 515` guard does not discriminate on high-dynamic-range data.**
The guard is presumably there to make the branch selective, but the trimmed
height-to-width ratio on this data is 32,279 (pos) and 14,656 (neg) — 63× and 28× the
threshold. *Both* polarities satisfy it comfortably, so the branch is governed entirely by
`hi > 35`, i.e. by whether a single unusually wide peak was measured anywhere in the file.
Positive mode crossed it at `max(width_seconds) = 37.7` while negative mode did not at
11.5. A threshold calibrated against QTOF peak heights will effectively always fire on
Orbitrap data, where heights run 10^4–10^9.

Measured widths, for reference:

| | min | q25 | median | q95 | max | normal branch would give |
|---|---|---|---|---|---|---|
| pos | 0.267 | 1.01 | 2.02 | 18.9 | 37.7 | `c(5, 43)` |
| neg | 0.267 | 0.34 | 0.66 | 5.03 | 11.5 | `c(5, 17)` |

### Suggested fix

The wide-peak adjustment is really about the *upper* bound — not searching for peaks far
wider than the bulk of the data. The lower bound has no reason to change with it:

```r
peakwidth <- c(
  ceiling(lo) + 4,
  if (hi > 35 && ratio > 515) (ceiling(hi) + 7) / 2 else ceiling(hi) + 5
)
```

For positive mode here that gives `c(5, 22.5)` — the wide-peak halving of the upper bound
is retained, and the lower bound stays data-driven and non-degenerate. A hard floor
(`max(3, ...)`) would also work as a safety net, though on this data the data-driven bound
is 5 and would bind first.

The `ratio > 515` threshold deserves separate attention: as an absolute cutoff on a
quantity with units of intensity-per-second, it is not scale-free and cannot transfer
between instrument classes. If the intent is "peaks are unusually wide relative to their
height", a dimensionless formulation — or a comparison against the distribution's own
spread rather than a constant — would travel better.

## Two smaller issues found while using the package

**`to_xcms(params, sample_groups)` cannot be used as documented for a real dataset.**
It validates `length(sample_groups) == length(params@files)`, but `paramounter()` is
intended to run on a handful of representative files while `PeakDensityParam` needs one
group per file in the *full* experiment (here 5 measured files versus 216 to be grouped).
As it stands the argument is only usable when you measure every file you intend to
process. Options: drop the length check and document `sample_groups` as belonging to the
target experiment; or remove the argument and document setting
`param@sampleGroups` afterwards. Note `xcms` does not export a `sampleGroups<-` method, so
the documented workaround has to be direct slot assignment.

**`to_xcms()` emits `Warning: stack imbalance in '<-', 2 then 4`.** Reproducible on this
data, appears benign, likely from S4 slot assignment inside the constructor chain. Worth a
look since it will alarm users.

## Open questions

1. Was the extreme-based rule a deliberate inclusiveness choice, or inherited from the
   original implementation without revisiting? That determines whether option 2 is a fix
   or a philosophy change.
2. What do the distributions look like on the QTOF data the method was developed against?
   If the tails are tight there, the extreme rule is fine for that regime and this is
   properly a dynamic-range issue to document rather than a defect.
3. Is `min_sn = max(3, min(sn))` intended to floor at 3 because low S/N is expected, or as
   a guard against degenerate input? On this data it returns 3.69 against a median of 20 —
   and the pilot above says that permissiveness is correct, so this is a question about
   intent rather than a problem to fix.
4. Should `peakwidth`, which uses both extremes of `width_seconds`, follow the tolerance
   treatment? It was not limiting on this data, but structurally it constrains peak
   construction rather than candidate count, which puts it on the tolerance side.
5. Does the tolerance/threshold split hold on the QTOF development data? If it does, it is
   worth stating in the package documentation as the rule for when the extreme is safe —
   it is a more useful framing than "the extremes are inclusive by design", because it
   tells a user which parameters to scrutinise when results look wrong.
