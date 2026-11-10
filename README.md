# Independent numerical replication of E. coli Figure 6B/C/D

This repository contains an independent implementation for a partial computational replication of the source article [Gervais et al. (2026)](https://doi.org/10.1038/s44320-026-00226-5). The article in `paper.md` and the original software in `src/` form a coherent review candidate. 

## Scientific boundary

This implementation independently reproduces **downstream numerical summaries** of the authors' cell-tracking and RealTrace-corrected data. It does not reimplement fluorescence inference from raw image pixels, prove physical identity across segmentation, assess independently randomized cultures, or identify a molecular memory carrier. The source of the 842 MB archive is [Zenodo 10.5281/zenodo.18059846](https://doi.org/10.5281/zenodo.18059846). Original author files are **not bundled** in this repository.

## Run

Use a fresh Python 3.13 virtual environment and the original file `Gervais2026_data.zip` downloaded separately from the public Zenodo source.

```bash
python -m pip install -r requirements-full.txt
python -m unittest discover -s tests -v
python src/reviewer_run.py --author-zip /path/to/Gervais2026_data.zip --out /new/fig6_analysis
```

Or, from authenticated extracted source members for a lower-disk-memory run:

```bash
python src/reviewer_run.py \
  --stress-feather /path/to/myframes_stress_exp.feather \
  --rmd /path/to/Gervais_et_al_analysis.Rmd \
  --realtrace-csv /path/to/realtrace_selected.csv.gz \
  --realtrace-receipt /path/to/realtrace_source_receipt.json \
  --out /new/fig6_analysis
```

The script rejects SHA-mismatched source inputs and refuses to overwrite any existing output directory. The official archive is expected to have SHA256 `0697803959ffa64b8b2bc849d716bc841bcb5afe8a42e92b53f397f8c7d7494b`.

## Expected outputs and testing

See `EXPECTED_RESULTS.json` for verified group-count, timepoint and ID-join invariants and the separate cloud execution modes. `references/` contains **small source-derived expected summaries**, not unpublished experimental measurements. `src/partition_weighted_bounds.py` computes source-ID partition bounds; `src/date_dose_audit.py` audits row weighting. `src/replot_all_figures.py` independently recreates **all six** figures from either freshly computed results or the shipped reference summaries (visual-only mode).

```
python src/reviewer_figures.py --groups references/fig6b_group_counts.csv --strata references/date_dose_binary_completion_bounds.csv --out /new/diagnostic_figures
```

The remaining four figures are created from the verified small summary tables in the same manner. This is only a visual smoke test and does not establish source-to-figure replay. To draw all six figures **from the actual newly computed source-run outputs**, first run the full source program, then:

```bash
python src/replot_all_figures.py --computed-run /new/fig6_analysis --out /new/all_figures
```

See `REVIEWER_GUIDE.md` for checks at the original source-data, observation-population and scientific-interpretation levels.

## Provenance and rights

The independently authored **code** is MIT licensed (`LICENSE`). The source authors' data, images, scripts, and any third-party publications retain their original rights; their ZIP and R Markdown are not assigned this code license. Paper/figure licensing is not represented by the code LICENSE. No author-source large file, secret or original author R Markdown is included.

An earlier independent complete-ZIP cold run in Google Colab verified 842,127,932 source bytes and all 21 ZIP-member CRCs for the five-stage core software. This revised extended seven-stage release was separately run using byte-authenticated extracted original Feather/Rmd and Colab-selected RealTrace records. A fresh complete-archive run of this revised version has **not** been completed. Independent external human peer review is not represented by any internal software test.

### Distinguishing fresh source replay from reference-only rendering

The shell entrypoint `sh run_source_and_figures.sh <original_zip> <new_output>` now
passes **the just-computed output directory** to the six-figure renderer. It refuses
missing or failed seven-stage receipts. For quick plot styling checks *only*, use
`python src/replot_all_figures.py --references references --out /new/smoke_figures`;
this reference-only command does not demonstrate complete source-to-figure reproduction.
The two modes are explicit in the JSON stdout receipt.


## Software citations and scientific scope

`paper.md` and `paper.bib` provide the companion article source; the publisher-ready PDF is separately prepared. `CITATION.cff` identifies the independently authored code and the source work; `AI_DISCLOSURE.md` records assistance. The original study, code provenance, and data rights are attributed in `DATA_AND_LICENSE.md`. This public source repository is distinct from a Zenodo software archive and from a ReScience C journal publication. No software-archive DOI is claimed.

The GitHub Actions offline job runs tests and reference-only rendering as a *smoke check*; no hosted runner is claimed to have already executed the original 842 MB ZIP.

## Exact tipping-point sensitivity

`src/tipping_point.py` calculates the minimum equal-partition contrast, its exact
rational representation, the smallest single-label departure from the two-arm
adversarial completion and an *illustrative* shared missing-outcome propensity.
It reads `references/source_id_partition_counts.csv` as a small precomputed
summary; the same source table is rebuilt by the full analysis runner. Run:

```sh
python src/tipping_point.py --partitions references/source_id_partition_counts.csv --out /new/tipping_audit
```

The reversal under unrestricted binary completion requires 130 jointly
opposing assignments across the two compared groups. This is a mathematical
identifiability stress test, not a behavioral model for missing cells.

## Verification status (as of 2026-11-10)

The complete upstream ZIP was independently checked in Google Colab for all 21 CRC members and the five core analysis stages. The expanded seven-stage implementation was separately replayed from authenticated original source members and six figures were regenerated. No new biological experiments, independent-culture variance estimates, or original-image reprocessing are claimed. Run the included 63 offline unit tests before assessing the replicated article scope.


## Companion preprint (distinct from the code and journal publication)

The author reports an independent preprint archived on Zenodo as version v1,
10 November 2026: **https://doi.org/10.5281/zenodo.24115006**. The
all-versions DOI is **https://doi.org/10.5281/zenodo.24115005**. These
identifiers describe the author's *preprint*; they do **not** identify this
source code release, the Gervais original data, or a ReScience C article.

The publicly accessible source code for this replication is at
https://github.com/dairongzhen3-creator/ecoli-stress-rescience-replication. The author may respond to journal review comments through GitHub.
The companion preprint and this manuscript may differ in publication date
or file bytes; no byte-identical match is claimed without an external comparison.
