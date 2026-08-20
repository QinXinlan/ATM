# EEG + fNIRS Avalanche Transition Matrix (ATM) Analysis

This repository contains Avalanche Transition Matrices (ATMs) and branching-parameter
(σ) results computed from simultaneous EEG + fNIRS recordings during a motor
imagery (MI) / motor execution (ME), left/right hand cued-movement task. This
README explains **what signal each ATM is computed from, what modelling choices
were made, why, and what the alternatives were**, so that the results can be
interpreted (and reproduced) without needing to read the full pipeline source.

> [!WARNING]
> Before using any specific subject/session, check §14. A subset of
> sessions have known data-quality issues — device failures, noisy
> recordings, incomplete sessions, or subject-specific physical/cognitive
> confounds — that were not excluded from this repository. Depending on the
> issue, these can invalidate the RT-window definition (EMG failure),
> depress signal-to-noise in a given band/modality (loose optodes, room
> noise, motion), bias σ toward artefactual transitions (motion artefacts,
> incomplete task engagement), or make a condition's data non-representative
> of neurotypical performance (e.g. a hemiparetic hand in the left/right
> task). See §14 for the full per-subject/session list before drawing
> conclusions from an affected session.

---

## 1. Pipeline overview

```
raw EEG (.csv) + raw fNIRS (.snirf)
        │
        ├─ clock alignment on one shared anchor event (fNIRS device clock ↔ PC clock)
        │
        ├─ EEG preprocessing: bad-channel detection (ABCD/BRP+LCSS), blink
        │  correction, re-referencing, notch filter
        │
        ├─ fNIRS preprocessing: optical density → Beer-Lambert → motion
        │  correction (peak_power) → optional short-channel regression → bandpass
        │
        ├─ EEG source localisation (dSPM, fsaverage, aparc)  ─┐
        ├─ fNIRS→parcel mapping (nearest_vertex, aparc — §8   ┘  → common atlas
        │  lists other channel→parcel methods, e.g. DOT)
        │
        └─ Avalanche Transition Matrix (ATM) + branching parameter σ
           computed independently per: modality × space × frequency band ×
           time window × condition
```

Every ATM in this repository is therefore tagged by **five axes**:

| Axis | Values |
|---|---|
| **Modality** | EEG, fNIRS |
| **Space** | `source_space` (EEG, dSPM parcels) / `sensor_space` (raw channels) · `parcel_space` (fNIRS, aparc parcels) / `sensor_space` (raw optode channels) |
| **Band** | EEG: `delta, theta, alpha, beta, gamma, broadband`. fNIRS: `broadband, vlf, slow_hrf, task_hrf, lf_mayer, resp` |
| **Window** | EEG: `full_window`, `rt_window` (+ `common_parcels_fnirs` variant). fNIRS: `full_window`, `hrf_window` |
| **Condition** | `MOTOR_IMAGERY=LEFT/RIGHT`, `MOTOR_EXECUTION=LEFT/RIGHT`, plus pooled/contrast rows (see §8) |

---

## 2. What signal actually feeds each ATM

| Modality × Space | Input signal to the ATM | How it was produced |
|---|---|---|
| **EEG `source_space`** | Per-parcel dSPM source time course (68 `aparc` cortical labels, `mean_flip` aggregation across vertices) | Band-filtered sensor data (IIR Butterworth, order 4) → forward model (3-layer BEM, `fsaverage`) → dSPM inverse → `mne.extract_label_time_course` (`dSPM` inverse method; §7 lists alternatives such as `MNE`, `sLORETA`, `eLORETA`, and other parcellation/source-space choices) |
| **EEG `sensor_space`** | Band-filtered raw EEG channel voltages (no source modelling) | Same IIR band filter applied directly to the re-referenced/cleaned sensor data |
| **fNIRS `parcel_space`** | Per-parcel HbO (and, by default, HbR — see §6) time course | Each optode-pair midpoint geometrically assigned to the nearest `fsaverage` pial-surface vertex (≤ 30 mm) → parcel = that vertex's `aparc` label → channels within a parcel averaged (`nearest_vertex` method; §8 lists alternatives such as distance-weighted, sensitivity/PMDF-weighted, and DOT) |
| **fNIRS `sensor_space`** | Band-filtered HbO channel concentrations (no parcellation) | Same IIR band filter applied directly to Beer-Lambert-converted, motion-corrected channel data |

All four variants use the **same ATM algorithm** (§3) — only the definition of
"node" changes (anatomical parcel vs. raw channel).

---

## 3. ATM algorithm (identical for every variant)

1. **Per-node thresholding.** For each node (parcel or channel) independently,
   over the analysed window/trial set:
   `threshold = mean(signal) + 2.0 × std(signal)` (`ATM_KWARGS: threshold_method='std', threshold_value=2.0`).
2. **Binarisation.** `active(t) = 1` if `signal(t) > threshold`, else `0`.
   Samples that are missing (e.g. a trial truncated by a device stopping
   early) are **NaN-padded upstream and always treated as inactive**, never
   as "below threshold" evidence.
3. **Avalanche segmentation.** An avalanche is any maximal contiguous stretch
   of time where **at least one node is active**. Avalanches are detected on
   the *concatenation of all trials* for a given condition (see caveat in
   §10.1 — this is important).
4. **Transition counting.** Within each avalanche, every pair `(node i active
   at t, node j active at t+1)` increments `ATM[i, j]` by 1.
5. **Row normalisation.** Each row of the raw transition-count matrix is
   divided by its sum, turning `ATM[i, j]` into `P(node j active at t+1 |
   node i active at t)` — a **transition-probability matrix**, not a raw
   co-occurrence count. Group-level ATM comparisons therefore compare
   *propagation structure*, not *how much* activity there was.

This is a **binary, single-lag (Markov, t→t+1), amplitude-agnostic** model of
propagation — see §10 for the modelling consequences.

---

## 4. Branching parameter σ

For every ATM, two independent estimators of the branching/criticality
parameter are computed and reported side by side:

- **σ_ratio** (primary/most sensitive): mean, over all frames where ≥1 node
  is active, of `Σ active(t+1) / Σ active(t)`. Computed directly from the
  same binarised activity that produced the ATM.
- **σ_eigenvalue**: dominant real eigenvalue of the row-normalised ATM. This
  is the theoretically "clean" estimator (a row-stochastic matrix has a
  fixed point at eigenvalue 1) but requires no time-series data — it can be
  computed for group-average ATMs where no `binary_data` exists.

**Interpretation:** `σ = 1` → critical, `σ > 1` → supercritical (explosive
propagation), `σ < 1` → subcritical (damped propagation). Prior work (Javed
et al., *J Neurosci* 2025) reports that Alzheimer's disease dynamics shift
toward supercriticality — this motivates using σ as a candidate biomarker
here.

**Confidence intervals** on σ_ratio are computed with
`ATM_CI_METHOD = 'bootstrap'` (non-parametric percentile bootstrap, 500
resamples) by default. The bootstrap seed is **deterministically derived**
per `(subject, session, window, band, condition)` from a fixed base seed
(`ATM_BOOTSTRAP_SEED_BASE = 20250715`), so re-running the pipeline on
unchanged data reproduces identical CIs. Other CI estimators are implemented
but not used by default: `wald_laplace` (normal approximation),
`wilson`/`agresti_coull` (binomial-proportion score intervals, ratios
clamped to [0,1]), `jeffreys` (Beta credible interval).

---

## 5. Time windows

| Modality | Window | Bounds | Rationale |
|---|---|---|---|
| EEG | `full_window` | MI: `0 → 5.0 s` (fixed) · ME: `0 → representative RT` (or 2.0 s fixed if no EMG) | Diagnostic / within-condition only. **Not a valid MI-vs-ME contrast** — the two conditions use different-length, different-reference windows. This is flagged explicitly in the output filenames/titles. |
| EEG | `rt_window` | `0 → representative RT` for **both** MI and ME | Valid MI-vs-ME contrast: both conditions share an identical window definition. |
| EEG | `common_parcels_fnirs` | Same full/RT windows, restricted to the subset of `aparc` parcels also covered by this session's fNIRS optode array | Session-matched, equal-node-count EEG/fNIRS comparison — **not** directly comparable to the full-parcellation EEG σ values (fewer nodes changes avalanche/A(t) statistics; see §10.5). |
| fNIRS | `full_window` | `3 – 12 s` post-trigger | Covers the full canonical HRF (peak ≈ 6 s). |
| fNIRS | `hrf_window` | `2 – 8 s` post-trigger | Isolates the rising flank of the HRF (initial cortical engagement), timed to be closer to the EEG RT window. |

`representative RT` = median (default; configurable to mean) of the
per-trial offline-detected EMG onset latency (slope-based burst detector on
a two-level rolling-SD ("global SD") stack of the AUX EMG channels; see
`EMG_DETECTION_METHOD='slope'`). It is a **single session-level scalar**,
not matched per trial — this is a modelling simplification (§10.6).

---

## 6. Frequency bands

| Modality | Bands | Notes |
|---|---|---|
| EEG | `delta 1–4`, `theta 4–8`, `alpha 8–13`, `beta 12.5–30`, `gamma 30–50`, `broadband 1–50` Hz | 4th-order Butterworth IIR. |
| fNIRS | `broadband 0.01–0.5`, `vlf 0.01–0.04`, `slow_hrf 0.04–0.1`, `task_hrf 0.01–0.1`, `lf_mayer 0.07–0.15` (Mayer/vasomotor waves), `resp 0.15–0.4` (respiration) Hz | 4th-order Butterworth IIR, applied to the parcel/channel time courses. `lf_mayer`/`resp` isolate physiological noise bands rather than task signal — useful as negative controls. |


---

## 7. EEG source localisation — choices made and alternatives

| Setting | Value used | Alternatives available in the pipeline |
|---|---|---|
| Anatomy | `fsaverage` template (no individual MRI) | Any subject with a FreeSurfer reconstruction |
| Inverse method | `dSPM` | `MNE`, `sLORETA`, `eLORETA` |
| Parcellation | `aparc` (Desikan-Killiany, 68 cortical labels, "unknown"/corpus callosum removed) | `aparc.a2009s` (Destrieux, ~150 labels), `HCPMMP1` (Glasser, 360 labels) |
| Source-space spacing | `oct6` (~1.6 mm) | `oct4`, `oct5`, `ico4`, `ico5` |
| BEM | 3-layer, conductivity `(0.3, 0.006, 0.3)`, `ico4` | 1-layer (MEG-style) |
| Orientation | loose `0.2`, depth `0.8` | fixed orientation |
| Label aggregation | `mean_flip` (sign-flip-corrected mean across vertices in a label) | `mean`, `pca_flip` |
| SNR / regularisation | single-trial SNR = 1.0 (λ² = 1, more regularised); condition-average SNR = 3.0 (λ² = 1/9) | any SNR |

Using a template brain means anatomical parcel boundaries are only as
accurate as the standard `fsaverage` surface warp for this population — no
individual cortical folding is captured.

---

## 8. fNIRS parcellation — choices made and alternatives

| Setting | Value used | Alternatives available in the pipeline / literature |
|---|---|---|
| Channel→parcel method | `nearest_vertex` — each optode-pair channel projected to nearest `fsaverage` pial-surface vertex (rejected if > 30 mm away), inherits that vertex's `aparc` label | Distance-weighted / spatial-kernel projection · Sensitivity-profile / PMDF-weighted (fOLD-style) · Diffuse Optical Tomography (DOT) — see comparison table below |
| Aggregation | `mean` across all channels assigned to the same parcel (min. 1 channel/parcel) | `median` |
| Chromophores parcellated | `both` — HbO and HbR, independently | HbO only |
| ATM node set (`FNIRS_ATM_USE_HBR`) | `True` (default) — dual-chromophore nodes `[parcel₁_HbO…parcel_N_HbO, parcel₁_HbR…parcel_N_HbR]`, giving a `(2N × 2N)` matrix capturing HbO↔HbO, HbO↔HbR, HbR↔HbO, HbR↔HbR transitions | `False` — HbO-only nodes, `(N × N)` matrix |

Because fNIRS has much coarser native spatial resolution (~3 cm) than the
68-parcel cortical atlas, many parcels are covered by very few — sometimes
just one — optode channel; parcel-level fNIRS "resolution" should not be
over-interpreted as anatomically precise.

### fNIRS channel→parcel methods: what's used vs. what's possible

`nearest_vertex` is the simplest of several established ways to go from
fNIRS channel space to an anatomical parcel representation. They differ in
how much of the actual light-propagation physics they model, and therefore
in cost and required inputs. Only `nearest_vertex` is implemented in this
pipeline; the rest are listed here as the alternatives that exist in the
wider fNIRS literature/toolchain, in order of increasing physical realism:

| Method | What it does | Key parameters it would require |
|---|---|---|
| **`nearest_vertex`** (used here) | Hard-assigns each channel's MNI midpoint to the closest `aparc` vertex/parcel; channels in a parcel are averaged. | `max_distance_mm` (rejection radius), `min_channels_per_parcel`, `aggregation` (`mean`/`median`) |
| **Distance-weighted / spatial-kernel projection** | Spreads each channel's signal across nearby parcels with a distance-decay kernel instead of a hard cutoff, so a channel can contribute partially to more than one parcel. | Kernel type (`gaussian`, `inverse_distance`), `fwhm_mm`/`sigma_mm`, `cutoff_radius_mm`, weight-normalisation scheme |
| **Sensitivity-profile / PMDF-weighted (fOLD-style)** | Uses a precomputed Photon Measurement Density Function per channel (Monte Carlo photon migration or diffusion approximation on a head model) to weight each channel's contribution by how sensitive it actually is to absorption changes there, rather than by geometric distance alone. | Head model (`Colin27`/`ICBM152`/subject MRI), tissue optical properties (μa, μs′ per layer), simulation engine (`MCX`, `tMCimg`, `MMC`) or precomputed fOLD tables, photon count, wavelength(s), sensitivity threshold |
| **Diffuse Optical Tomography (DOT)** — full image reconstruction | Solves an inverse problem to reconstruct absorption-change (→ HbO/HbR) images at every voxel/vertex of a head model, then aggregates reconstructed voxels into parcels afterward. | Forward sensitivity/Jacobian matrix (Monte Carlo or diffusion-equation solver, e.g. NIRFAST/Toast++/Redbird), head mesh, regularisation method + parameter λ, dual-wavelength combination, voxel/mesh resolution |

`nearest_vertex` and the distance-weighted variant are purely geometric,
data-independent operators — they don't depend on the recorded signal at
all. The sensitivity-weighted and DOT methods are *not* purely geometric:
they depend on a modelled or simulated light-propagation forward operator,
and DOT in particular is a regularised inverse problem in the same family
as EEG source localisation (§7). See §9 for what that difference means for
the order trial-averaging is done in.

---

## 9. Averaging order: trial-average-then-map, or map-then-trial-average?

EEG source localisation and fNIRS parcellation resolve the "average
trials, or map to anatomy, first?" question in opposite ways — and this is
intentional, not an inconsistency. The table below summarises the order
used by every method mentioned in §7/§8; the reasoning for each follows.

| Modality / method | Order used | Because |
|---|---|---|
| EEG source localisation — any of `MNE`, `dSPM`, `sLORETA`, `eLORETA` (§7) | **Average, then localise** | Regularised inverse problem: optimal regularisation depends on the SNR of the data being inverted |
| fNIRS geometric parcellation — `nearest_vertex` (used here), distance-weighted (§8) | **Order doesn't matter** — parcellate-then-average ≡ average-then-parcellate | Fixed, data-independent linear spatial average; no regularisation term |
| fNIRS sensitivity/PMDF-weighted (§8) | Same as geometric methods, **provided** the sensitivity profiles are precomputed from a template/simulation rather than re-fit per recording | Still a fixed, data-independent linear weighting in that case |
| fNIRS DOT — not implemented here (§8) | Would need **average, then reconstruct**, like EEG | Regularised inverse problem; same SNR-dependence as EEG source localisation |

### Why EEG must average before it localises

EEG source localisation is a regularised, ill-posed inverse problem. All
four candidate methods solve a weighted minimum-norm problem under the
same regularisation parameter, λ² = 1/SNR², and MNE-Python exposes this
explicitly: `apply_inverse` (condition-average evoked data, SNR ≈ 3) and
`apply_inverse_epochs` (single-trial data, SNR ≈ 1, scaled via `nave`)
accept the identical `method` argument but need a different `lambda2`
depending on whether the data being inverted has already been averaged.
Hincapié et al. (2016, *Comput. Intell. Neurosci.*,
doi:10.1155/2016/3979547) show that minimum-norm-based estimates — the
family `MNE`/`dSPM`/`sLORETA`/`eLORETA` all belong to — require different
regularisation depending on whether they are computed from averaged or
single-trial data; the optimal λ² is a function of the SNR of the data
actually being inverted, not a fixed constant. Averaging first and then
inverting the evoked response at the condition-average SNR is therefore
the correct order for the primary/ATM source estimate, whichever of the
four methods is selected; the single-trial SNR (§7) is used separately,
for genuine single-trial inversion (the per-trial CSVs in §13).

### Why the fNIRS geometric methods don't care about order

fNIRS parcellation, by contrast — for the geometric methods (§8's
`nearest_vertex`, used here, and the distance-weighted variant) — has no
regularisation term: it is a fixed linear spatial average of the channels
whose projected position falls inside a parcel (Tsuzuki & Dan, 2014,
*NeuroImage*, doi:10.1016/j.neuroimage.2013.07.025; see also `mne-nirs`'s
region-of-interest averaging). Being a data-independent linear operation,
it commutes with trial-averaging: parcellate-then-average and
average-then-parcellate are numerically identical. The per-trial fNIRS
parcel time courses (§13) are therefore parcellated per-trial first only
to retain genuine single-trial parcel time courses for downstream
per-trial ATM/GLM/hub-figure use, and averaged afterward as a byproduct —
the fNIRS analogue of EEG's separate single-trial pass.

### The exception: DOT would behave like EEG, not like the geometric methods

This commutativity argument does not extend to DOT. Full diffuse optical
tomography reconstruction is, like `MNE`/`dSPM`/`sLORETA`/`eLORETA` source
localisation, a regularised inverse problem — its Tikhonov (or equivalent)
regularisation parameter λ is a function of the SNR of the data actually
being inverted, exactly the property that makes EEG's average-before-
inverting order matter in the first place. Were DOT ever adopted in this
pipeline, it would need the **same** treatment §7 gives EEG: reconstruct
once at the single-trial SNR for genuine per-trial CSVs, and separately at
the (higher) condition-average SNR for the ATM/GLM primary estimate —
parcellate-then-average and average-then-parcellate would **not** be
interchangeable for that method. The sensitivity/PMDF-weighted method sits
in between: it is still a fixed, data-independent linear weighting (so it
commutes like `nearest_vertex` does) as long as the sensitivity profiles
themselves are precomputed from a template/simulation rather than re-fit
per recording.

---

## 10. Modelling hypotheses and how they shape the ATM

This section makes explicit every assumption baked into the numbers, since
each one measurably changes σ and the ATM structure.

1. **Trials are concatenated end-to-end before avalanche detection, with no
   gap inserted at trial boundaries.** For a condition with N trials in a
   given window, the analysed signal is one continuous array formed by
   concatenating the N per-trial windows in trial order. An avalanche can
   therefore span the boundary between trial *k* and trial *k+1* if activity
   is still above threshold at the last sample of trial *k* and the first
   sample of trial *k+1*. This is a deliberate simplification for
   statistical power (more usable frames per condition) but means avalanche
   *durations/sizes* — and therefore σ_ratio — are not strictly
   trial-isolated.
2. **Binary, amplitude-agnostic activation.** Once a node crosses threshold
   it contributes identically to every avalanche/transition regardless of
   how far above threshold it is. Two very different-magnitude bursts are
   indistinguishable to the ATM.
3. **Single-lag (t → t+1) Markov transition model.** Only immediately
   consecutive time samples are counted; no longer-range or lagged
   coupling is captured (this is standard for ATM/avalanche analyses but
   worth stating).
4. **Per-node threshold is fixed within each window/trial-set** (mean + 2 SD
   computed once over the whole analysed signal for that node), not
   adaptive over time — so slow drift within a condition could bias which
   samples cross threshold.
5. **Common-parcels-with-fNIRS EEG variant uses fewer nodes** than the
   full-parcellation EEG ATM. Because A(t) (the number of simultaneously
   active nodes) and avalanche statistics are node-count-dependent, σ from
   this variant is **not** on the same scale as σ from the full-parcellation
   variant, even for the same session/condition — it exists purely to give
   an EEG/fNIRS comparison at matched spatial resolution.
6. **`representative_rt` is a single scalar per session** (median/mean EMG
   onset latency across all execution trials), used to define the RT window
   for every trial equally. It is not trial-matched — a genuinely fast or
   slow trial is windowed identically to the session average.
7. **fsaverage / mean_flip / dSPM assumptions** (§7): anatomical precision
   is bounded by the template warp; `mean_flip` assumes broadly coherent
   dipole orientation within a parcel (opposing orientations within one
   parcel can partially cancel).
8. **Row-normalisation discards magnitude information.** ATM contrasts (and
   the eigenvalue estimator of σ) compare the *shape* of propagation
   (transition probabilities), not the raw amount of avalanche activity;
   raw avalanche count/size statistics are recorded in each ATM's metadata
   but are not what group_contrasts/connectivity figures visualise.
9. **Group-average and contrast σ rows are unweighted means** across the
   contributing per-condition σ_ratio values — sessions/conditions are not
   weighted by trial count when pooled (see `SIGMA_COMBINATIONS` in §11).
10. **Full-window MI-vs-ME is explicitly not a valid contrast** (§5) — any
    comparison across these two conditions should use the RT-window
    variant instead.

---

## 11. Group-level σ rows exported per session

Every session's `sigma_ratio_*.csv` contains these 12 standard rows:

| Row | Type | Definition |
|---|---|---|
| `MOTOR_EXECUTION-LEFT/RIGHT`, `MOTOR_IMAGERY-LEFT/RIGHT` | per-condition | direct σ for that single condition |
| `MOTOR_EXECUTION`, `MOTOR_IMAGERY`, `LEFT`, `RIGHT` | group average | unweighted mean of σ_ratio across the matching conditions |
| `LEFT_vs_RIGHT`, `IMAGERY_vs_EXECUTION` | cross-group delta | `mean(group A) − mean(group B)`; CI = conservative widest-gap interval (`ci_lo_A − ci_hi_B`, `ci_hi_A − ci_lo_B`) |
| `ME_LEFT_vs_ME_RIGHT`, `MI_LEFT_vs_MI_RIGHT` | within-condition laterality delta | same delta logic, restricted to one task (execution or imagery) |

`sigma_eig` is left `0`/undefined for all delta rows — the dominant
eigenvalue of a *difference* matrix is not a meaningful branching-process
quantity.

---

## 12. Population / group labelling used in figures

Clinical grouping is taken directly from the `Group` column of
`Participants*.csv`, matched to each session by subject ID and handling
longitudinally re-tested subjects via the `Sessions` column. It is **not** a
simple MoCA cutoff applied post-hoc; it's a pre-computed clinical
categorisation that combines MoCA and MMSE (education-adjusted):

| Group | Criterion |
|---|---|
| `HC` | MoCA **and** MMSE both indicate normal cognition |
| `MCI` | MoCA **and** MMSE both indicate impairment (concordant) |
| `MCI_Discordant` | MoCA and MMSE **disagree** — one indicates impairment, the other doesn't (most commonly MoCA below threshold with MMSE still ≥ 24, a known MMSE ceiling-effect pattern) |
| `AD` | MoCA and MMSE both indicate impairment, with MoCA in the severe range — clinically consistent with dementia |
| `NA` | Both MoCA and MMSE missing — left ungrouped |

Normal-vs-impaired thresholds are **not a single fixed number** — the
per-subject `Reason` field shows education-adjusted MoCA cutoffs (standard
≥ 23; +1 point credit for ≤ 12 years education per Carson et al.; a
Youden-optimal ≥ 24 cutoff for high-education subjects per Pugh et al.) and
an analogous MMSE split (standard ≥ 24; ≥ 27 for 16+ years education per
O'Bryant et al.). When only one of MoCA/MMSE is available, that score alone
determines the group (e.g. MMSE-only classification for illiterate subjects
with no MoCA).

---

## 13. Repository layout (per session)

```
per_session/S<id>_Sess<id>/
  atm_eeg/<method>/{full_window,rt_window,common_parcels_fnirs}/<band>/
      csv/            per-trial + average ATM matrices (CSV)
      figures/        per-trial + average ATM heatmaps (PNG)
      hub_csv/, hub_figures/   same, restricted to top-N "hub" parcels
                                (RE-COMPUTED from only those parcels, not a
                                 sub-index of the full ATM)
      group_contrasts/         per-condition + pooled-group ATM heatmaps
      connectivity/            MI-vs-ME / LEFT-vs-RIGHT ATM contrasts +
                                network-circle figures
      sigma_ratio_eeg_<window>_<band>.csv
  atm_fnirs/{full_window,hrf_window}/<band>/   (same structure)
  glm/, source_localised/, signal_quality/, EMG/, ...
group/
  branching_parameters.csv       one row per session, every σ combination
  emg_latency_summary.csv
  participant_comparison_{moca,mmse}.png
  group_overview.png, longitudinal_trajectory.png
```

## 14. Data-quality caveats 

The remarks below are copied/summarised from the Remark column
of Participants*.csv, filled in only when something notable happened
during acquisition for that subject/session. An empty Remark field means
nothing unusual was logged — not that the session was independently
verified clean.

### EMG / input-device issues
- **S03, Sess01** — EMG not working; keyboard used as input instead. Any
  `representative_rt` / RT-window results (§5) for this session are **not**
  EMG-derived and are not comparable to other sessions on that basis.

### fNIRS signal / montage issues
- **S22, Sess01** — fNIRS optode D06 ("hat") fell down during the session.
- **S25** — fNIRS optode D13 may have shifted ("might have fell a bit").
- **S31, Sess01** — fNIRS: no signal. EEG: very noisy. Session was not
  completed to the end.
- **S48, Sess02** — 6th lantern trial: subject moved their whole body to
  move the fNIRS cable box (likely a motion artefact at that trial).

### Noise issues
- **S16, Sess05** — Phone rang during the session.
- **S32** — Subject scratched their head several times during the
  sessions (possible movement artefact).
- **S46, Sess01** — Very noisy (other people playing cards in the room
  during recording). 

### Task comprehension / session completion issues
- **S26** — Could not understand the task instructions.
- **S30** — Could not understand the task instructions (audio cues).
- **S34** — Sess01 not completed; Sess02 completed but with verbal cues
  given to the subject (deviation from the standard protocol).
- **S35** — Sess01 not completed; Sess02 completed but with verbal cues
  given to the subject (deviation from the standard protocol).
- **S43, Sess02** — Not completed; many errors, including the lantern/arrow
  code being wrong on trials 5, 13, 16, and 19.

### Physical / motor confounds (relevant to left/right hand task)
- **S47** — Gave up on several MMSE items (read, write, draw) and the MoCA
  draw item; recent surgery on the right arm meant she could not hold a
  pen. Also unfamiliar with the local area (came to Hangzhou only for the
  surgery), so location-related cognitive-test questions may
  under-represent her baseline cognition.
- **S48** — History of stroke (1999, 2010) in the right hemisphere; left
  hand does not work well. Also has high blood pressure. **Left-hand
  motor execution/imagery data for this subject should be interpreted with
  this impairment in mind**, not treated as a neurotypical left-hand
  condition.


None of these subjects/sessions have been excluded from the repository —
the remarks are preserved so each downstream analysis can decide per-case
whether to exclude, flag, or keep them.
