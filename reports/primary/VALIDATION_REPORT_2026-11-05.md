# Independent validation of stable query certificates in biological reannotation

**Date:** 5 November 2026  
**Manuscript:** Biological Findings That Survive Reannotation: Stable Query Processing Across Taxonomic Reference Systems  
**Test scope:** newly executed tests from available, public QIIME 2 Moving Pictures input artifacts.  
**Status:** MAIN CERTIFICATES REPRODUCED; MODEL-CONDITIONAL COUNTERFACTUAL WITNESSES VERIFIED; SAME-SEMANTICS SYNTHETIC BENCHMARK RE-RUN. EXTERNAL COHORT VALIDATION NOT ESTABLISHED.

## 1. Audited inputs

- Public QIIME 2 Moving Pictures input: **770 shared sequence features, 34 samples, 153,807 reads**, two subjects.
- SILVA 138 and Greengenes2 2024.09 classification of the **same 770 features** (local qza files in the package).
- Nine query contexts: three anatomical-site contrasts × phylum/family/genus levels.
- This is not an additional independently sampled human study.

## 2. Exact repeat of the original stage 25 experiment

The original `run_stage25.py` completed successfully. Five frozen outputs were reproduced **byte-for-byte**: `contrast_certificates.csv`, `contrast_summary.csv`, `topk_sufficient_certificates.csv`, `agreement_failure_examples.csv`, `validation.json`.

- **2,097** additive contrast queries.
- **3** cases: observed same-sign agreement under both complete references but no certificate under independent per-feature switches.
- **27** top-k settings, of which **11** show observed set agreement and **10** satisfy the interval certificate.
- Important: all reference-switch worlds apart from the two complete references are modeling scenarios, not observed biological outcomes.

## 3. New exhaustive correctness tests

- **1,440** exhaustive finite-world comparisons, including **720** random disjoint-block monotonicity checks.
- Maximum absolute difference between brute-force and algebraic certificate boundaries: **1.78e-15**.
- All checks PASS, under the declared independent / disjoint-block choice semantics.

## 4. New exact pairwise top-k analysis and explicit witness

We replaced independent score-range separation with the **exact pairwise worst-case margin** under the same feature-independent model: 

`min_(world) [Q_a - Q_b] = sum_feature min(c^S_(f,a)-c^S_(f,b), c^G_(f,a)-c^G_(f,b))`.

- **3,860** pairwise checks over 11 observed-agreement settings.
- **10/11** strictly certified via exact pairwise bounds; the original sufficient test already certified 10/11, so this extension provides **no additional positive certificates on the current data**.
- The remaining observed-agreement setting has an **explicit constructed ranking reversal**: `left palm - right palm`, phylum, top-1: Pseudomonadota loses to Firmicutes in an independently recombined assignment world. The exact pairwise margin is **-0.00521511** (normalized relative-abundance contrast).
- This is a **counterfactual mathematical witness**; it does not claim that a future taxonomy update actually generates that mixed world.

## 5. Three independently constructed sign-reversal witnesses

For each earlier non-certified contrast we constructed the exact extremizing assignment by choosing, feature by feature, the assignment that minimizes a positive contrast or maximizes a negative contrast. All three witnesses **strictly reverse the observed sign** under the declared independent model.

| Contrast (left palm minus right palm) | SILVA | GG2 | Constructed mixed world | Number of GG2 choices |
|---|---:|---:|---:|---:|
| Aerococcaceae (family) | -0.000071 | -0.002187 | +0.000193 | 1 / 770 |
| Corynebacterium (genus) | +0.002459 | +0.006047 | -0.003976 | 4 / 770 |
| Prevotella (genus) | +0.002765 | +0.016298 | -0.000117 | 5 / 770 |

**Interpretation:** such a witness proves that two complete-reference agreements do not imply stability in a larger hypothetical assignment space. It is not an empirical biological sign reversal and not a test of statistical significance.

## 6. Same-semantics incremental maintenance benchmarks (synthetic trees)

Original stage19 and stage20 code was re-executed on this host with NumPy; both completed without failure. Stage19: 80 exhaustive small-tree tests PASS (maximum error 3.55e-15), and incremental vs full tree-DP outputs agree. Stage20: 18 runs × 120 updates covering two tree sizes, three query frequencies and three repeats; eager, lazy and full recomputation return identical certificates (reported maximum absolute error 0).

| Tree nodes | Read frequency | Median eager (total ms) | Median lazy (total ms) | Median full rebuild (total ms) |
|---|---|---:|---:|---:|
| 255 | Every 1 update(s) | 31.75 | 34.07 | 546.65 |
| 255 | Every 5 update(s) | 26.21 | 19.74 | 107.92 |
| 255 | Every 20 update(s) | 27.81 | 11.45 | 27.01 |
| 1023 | Every 1 update(s) | 39.11 | 43.48 | 2410.75 |
| 1023 | Every 5 update(s) | 37.91 | 28.30 | 452.96 |
| 1023 | Every 20 update(s) | 34.08 | 19.93 | 115.71 |

These benchmarks operate on **synthetic dependency trees**, not biologically validated taxon-transition constraints. The tested methods share semantics, but are written in Python and NumPy; this is **not** a comparison against compiled native incremental-view maintenance or a mature database engine. Timings should not be compared across different machines or interpreted as general speed-ups.

## 7. Biological constraint audit

Stage18 input audit was re-executed on 770 real paired feature labels. Literal read-weighted agreement: phylum **58.27%**, family **78.93%**, genus **28.13%**. Among SILVA groups having two or more features, mixed changed/unchanged label groups appear in **9/20 phyla** and **44/78 families**. These observations **do not validate** group-wise independent reference-switch probabilities or universal hard constraints.

## 8. Separate 62-sample rodent dataset: failed comparability gate

Checked the public `tallon-research/Dual-Reference-Taxonomic-and-Diversity-Profiling-of-Microbial-16S-rRNA-Amplicon-Data-using-QIIME-2` repository via GitHub connector on 5 Nov 2026. The provided `results/silva/taxonomy.tsv` file contains 1,335 IDs, whereas `results/gg2/taxonomy.tsv` has 507 IDs; their **literal feature-ID intersection is empty (0)**. These exported rows are not in an established shared ASV identifier domain. Paired genus-abundance summaries also do not guarantee equal read totals. Therefore **do not** present this cohort as a second independent, exact test of our feature-switch certificates without reconstruction from a common feature table and stable identifier mapping.

## 9. Overall conclusion and remaining ICDE risk

- **PASS:** code-level reproducibility of the stage25 real-data experiment; exact independent/block semantics; exact coupled pairwise top-k; explicit conditional witnesses; synthetic same-semantics maintenance.
- **NOT PROVEN:** biologically plausible independent switching, causal or statistically significant reversals, external-cohort replication, equivalence of taxon synonyms across references, parity with optimized production-grade database IVM.
- **Scientific decision:** the new exact pairwise witness strengthens the correctness story but does not alone remove the prior high ICDE Research Track novelty risk.

Figures `Fig_NEW_1_contrast_witnesses.png` and `Fig_NEW_2_topk_coupled_certification.png` report actual computed outcomes. Both captions must state that derived mixed-world results are conditional on an explicit counterfactual model.
