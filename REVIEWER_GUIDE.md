# Independent replication reviewer guide

## Study to be checked

Original: Gervais et al., *Molecular Systems Biology* 2026, DOI https://doi.org/10.1038/s44320-026-00226-5. Data: Zenodo DOI https://doi.org/10.5281/zenodo.18059846.

This is a **partial computational** replication of the archived R analysis aggregations for Figures 6B/C/D. Upstream microscopy processing and the RealTrace biological inference model were **not** independently re-fitted from image intensities.

## Execute from the original ZIP

1. Download the author ZIP (842,127,932 bytes) separately from the official Zenodo DOI, keeping it outside this repository.
2. On a clean Python 3.13 machine, install packages with `pip install -r requirements-full.txt` and save `python -m pip freeze` and OS details.
3. Run `python -m unittest discover -s tests -v` and preserve full stdout/stderr.
4. Run `python src/reviewer_run.py --author-zip /path/to/Gervais2026_data.zip --out /new/unoccupied/directory`. Existing output paths must be rejected.
5. Verify `ANALYSIS_RUN_RECEIPT.json`, the measured SHA-256 of original members, the number of 21 ZIP member CRCs and the output counts in `EXPECTED_RESULTS.json`.
6. Run `python src/replot_all_figures.py --computed-run /new/unoccupied/directory --out /new/figure_directory` to recreate the six manuscript figures **from the source-run outputs**, not the frozen references. The command refuses an absent or failed run verdict. Separately, `--references references` is a visual smoke test only. Compare the regenerated PNGs with the manuscript figures and record exact differences. Check `FIGURE_SOURCE_SHA256_RECEIPT.json` in the figure output: it records the fresh run-verdict SHA, five derived input SHA-256 values, six image SHA-256 values and the certified original-source locks. Across operating systems or Matplotlib versions, visual images may differ bytewise while the underlying numeric results agree; use numerical and visual checks as the substantive replication criterion.

## Interpretation checks

- Figure 6B: 1,153 software cycles: 375 Survive, 533 No Growth, 245 Unknown; source-published plotting code yields `S/(S+N)` for the fraction and `sqrt(p(1-p)/(S+N+U))` for error whiskers. The journal's figure legend calls the latter 's.d.'; report the semantic distinction without inventing experiment-level uncertainty.
- Figure 6C: 616 bolA group-time summaries of author-inferred corrected concentrations, not 616 independent cultures.
- Figure 6D: eight reporter-by-treatment summaries, 32,400 valid corrected observations in the pooled record-level analysis; the published original caption says two independent experiments and two concentrations. Do not recast time points as biological replicates.
- Cross-panel software keys: 1,145 ID overlaps, 1,139 with both reporters; software ID equality does not certify physical barcode identity.
- Finite missing-outcome sensitivity: all 245 Unknown outcomes hypothetically completed to either class. Overall group ranking is invariant only to that binary completion with original **cycle-equal weights**.
- A second descriptive estimand gives equal weights to 37 source-ID position-like partitions. It permits a contrast below 0 by 0.0718623 percentage points for Unperturbed-minus-CAM-after. The positions are not certified randomized cultures or independent experimental units.
- The two September 2022 CAM date×dose strata do not guarantee the pooled ranking under all binary outcome completions; both June 2023 CAM strata do.
- Retrospective AUC sensitivity bounds cover 0.5 in all four reporter/condition combinations; no predictive out-of-sample validation or causal mechanism is claimed.

A substantive independent review should identify the actual reviewer, the code version, original ZIP input hash, environment, exact tests and reproducibility exceptions, and distinguish computation from biological inference.
