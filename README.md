# Biological findings that survive reannotation — ICDE 2027 artifact

**Paper:** *Biological Findings That Survive Reannotation: Stable Query Processing Across Taxonomic Reference Systems*

This is a **public-release candidate** assembled from four frozen research archives. It contains the final IEEE paper in LaTeX/PDF, primary sample-level results, the two-anchored budgeted certificate method, hierarchical and incremental benchmarks, an independently computed second-data-set validation, figures, tests and an integrity manifest.

> **Status:** The files are prepared locally. They have **not** been uploaded to GitHub by preparing this archive. The authors must approve the code/data redistribution terms before making a repository public. There is no standalone software license grant in this candidate.

## Two datasets — different scientific claims

| Cohort | Fixed observations | Reference assignments | Validated claims |
| --- | --- | --- | --- |
| QIIME 2 **Moving Pictures** | 770 shared features, 34 samples, 153,807 reads | SILVA 138 / Greengenes2 2024.09 | Additive within-site contrasts, sign- and top-k certificates, budgeted and block-conditional stress tests. Two subjects only. |
| Independent **gut-to-soil** | 41,328 paired ASV IDs, **21,227,600 pooled reads** | SILVA 144 / GTDB r232 | Pooled rank robustness: certified Pseudomonadota top-1 and a two-ASV counterexample to an agreed informative-genus top-1. NOT subject-level or phenotype contrast validation. |

**Neither dataset validates the biological probability of counterfactual mixed-reference assignment worlds.** A family switch is an explicit model constraint, not a learned biological law. The synthetic benchmark tests the same operator under matched computational semantics; it is not a production F-IVM/DBSP comparison.

## Quick start (does not download any data)

Requires Python >= 3.11. Fast checks use only the standard library:

```bash
python verify_artifact.py
python verify_results.py
python external/check_frozen_results.py
python external/reproduce_independent_validation.py /tmp/not_present --skip-computation --out build/external_smoke
```

The last command runs the exact finite-world enumeration tests without reading the external QZV. It is **not** a second-data-set rerun.

## Full primary-data regeneration (Moving Pictures)

`work/inputs/` contains three small QIIME 2 artifacts derived from public tutorial data and the team's paired classification outputs. A rights/attribution review is required before publicly sharing their binaries (see `DATA_AVAILABILITY.md`).

```bash
python -m pip install -r requirements.txt
python work/stage25/run_stage25.py
python correlated_certificates.py
python hierarchical_certificates.py
python matched_semantics_benchmark.py
python make_figures.py
python verify_results.py
```

The budgeted scripts write into versioned `results/` and may replace the frozen files; run in a fresh clone for comparison. Stage 25 writes its output into `work/stage25/`, while the original baseline is in `frozen/primary/`. The full regeneration is more computationally expensive than the fast checks.

## Independent external-data regeneration (gut-to-soil)

Retrieve the *exact* externally published paired-classifier QIIME 2 `.qzv` via its verified source URL. The downloader checks SHA-256; the original binary is **not included** here.

```bash
python external/download_external.py
python external/reproduce_independent_validation.py external/input/asv-seqs-gtdb-silva-fixed.qzv --out build/external_rerun
python external/check_frozen_results.py
```

`external/reproduce_independent_validation.py` enumerates all relevant top-1, top-3 and top-5 pair margins for phylum, family and genus; `external/results/` holds the frozen original summary and margin curves, including the two switched ASV identifiers. Public source: [QIIME 2 tutorial](https://forum.qiime2.org/t/using-multiple-taxonomy-classifiers-for-comparison-or-in-filtering-workflows/34349). Original dataset: [Zenodo 10.5281/zenodo.15390940](https://doi.org/10.5281/zenodo.15390940).

## Directory map

- `main.tex`, `figures/`, `paper.pdf`: checked 13-page IEEE manuscript and figures; last page references only.
- `correlated_certificates.py`, `hierarchical_certificates.py`, `matched_semantics_benchmark.py`: the new query method, tree DP and incremental benchmark.
- `work/inputs/`, `work/stage25/`, `work/results/`: real paired primary inputs, Stage 25 script, and frozen reference tables.
- `frozen/primary/`, `results/`: original primary certificate outputs and budgeted/hierarchical results.
- `external/`: public provenance, re-run script, downloader, label comparison, figures and external validation outputs.
- `tests/`: earlier stage 18--20 test sources, interval/ranking audit scripts.
- `reports/`: prior independent quality-control PDFs and text reports.
- `verify_artifact.py`, `verify_results.py`, `SHA256SUMS.txt`: hash integrity, frozen-result checks, automated GitHub Actions workflow.

## Exact status of reproducibility

- **Pass:** frozen artifact SHA-256 verification, primary result consistency checks, and external synthetic finite-world tests.
- **Full rerun:** primary 34-sample experiments can be re-executed from the included QIIME inputs; requires installed dependencies and runtime.
- **External original-source rerun:** requires fetching the externally hosted QZV; hash check avoids accidental source drift.
- **Not yet established:** independent biological calibration of the mixed-world constraints, subject-level contrasts on the external dataset, native database performance comparisons.

## Submission and rights

ICDE 2027 Research is single blind. The CMT supplemental-material URL must be a **real, anonymously accessible public GitHub release or repository**. Do not use an invented or private URL. No observer-tracking links are required. See `RIGHTS_AND_RELEASE.md` and `CMT_SUPPLEMENTARY_FIELD.md` before publication.
