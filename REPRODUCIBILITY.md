# Reproduction checklist

## Test matrix

| Experiment | Source | Exact run | Frozen results | Status |
| --- | --- | --- | --- | --- |
| Primary additive contrasts | Moving Pictures `work/inputs/` | `python work/stage25/run_stage25.py` | `frozen/primary/` | Inputs & code included. |
| Budgeted and block certificates | Same primary inputs | `python correlated_certificates.py` | `results/` | Inputs & code included. |
| Hierarchical clade-event DP | Same primary inputs | `python hierarchical_certificates.py` | `results/` | Includes exhaustive tests. |
| Incremental exact-operator benchmark | Synthetic aggregate updates | `python matched_semantics_benchmark.py` | `results/` | Measures Python implementations only. |
| External pooled rankings | Published gut-to-soil QZV | `python external/download_external.py` then `python external/reproduce_independent_validation.py external/input/asv-seqs-gtdb-silva-fixed.qzv --out build/external_rerun` | `external/results/` | Full external binary is retrieved from publisher source, SHA-256 checked. |
| Fast release QC | no outside data | `python verify_results.py` | all frozen results | Runs on Python standard library. |

All reproducible changes are restricted to binary per-feature or block assignments from two observed taxonomic classifiers. A negative witness is a mathematical constructed world, not an observed third taxonomy release.
