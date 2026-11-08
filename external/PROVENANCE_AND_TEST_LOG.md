# Evidence and reproducibility audit — 2026-11-06

## Public source directly retrieved and evaluated

- Article/tutorial: https://forum.qiime2.org/t/using-multiple-taxonomy-classifiers-for-comparison-or-in-filtering-workflows/34349
- Data DOI linked there: https://doi.org/10.5281/zenodo.15390940
- Paired QIIME 2 visualization (hosted by tutorial author):
  `https://www.dropbox.com/scl/fi/dycygxvmvjzrn0qdgiqvs/asv-seqs-gtdb-silva-fixed.qzv?rlkey=yemtsq1s5qprs4e808878z73s&dl=1`
- Retrieved 2026-11-06 into a **network-enabled analysis workspace**, not into the offline local container.
- SHA-256: `e8aa584032374072c3b8607a870fadd179bd9c6e5198d37f605c71fcddb03cfb`
- Original file size: 21,065,310 bytes.
- Embedded QIIME 2 `tabulate-seqs` JSON was inspected directly: 41,328 records, with aligned feature IDs, SILVA 144 and GTDB r232 labels, and per-feature frequencies. No ASV feature crosswalk was inferred from literal taxon labels.

## Actual computations

The script executed in that environment parsed the original QZV and checked unique ASV IDs, integer count integrity and common 21,227,600 total reads, then computed:

- read-weighted literal agreement under phylum/family/genus;
- top-1/top-3/top-5 reference agreements;
- exact worst-case shared-world pairwise margins of top-k versus **all** out-of-set competitors;
- comparison of per-ASV edit units and SILVA-family edit units;
- a secondary genus-label-filter sensitivity analysis;
- constructive witness for radius 2 after informative-genus filtering;
- exhaustive brute-force correctness of the optimized budgeted-margin primitive: **1,320 cases**, error max `8.881784197001252e-16`.

The optimized primitive was also checked separately using the local package with the same RNG seed and 1,320 cases. Derived numeric results were exported and frozen into `frozen_independent_validation.json`; the original 21 MB QZV was not exported or redistributed.

## Key results recorded in actual network-workspace stdout

- 41,328 ASVs; 21,227,600 pooled reads; 314,497 summed ASV-by-sample observation occurrences; max 825 samples for a single ASV.
- Phylum: both annotated 39,214 ASVs / 20,630,628 reads; same nonempty literal label 31,250 ASVs / 17,377,393 reads.
- Family: both annotated 32,754 ASVs / 18,483,547 reads; same nonempty literal label 16,042 ASVs / 9,528,607 reads.
- Genus: both annotated 21,286 ASVs / 11,283,622 reads; same nonempty literal label 6,243 ASVs / 5,669,841 reads.
- Phylum Pseudomonadota: top-1 under both complete references, certified against 115 competitors in all binary per-ASV worlds, worst margin 2,891,234 pooled reads.
- Genus label sensitivity: first named genus in both full references is Rhodococcus after excluding non-informative labels; exact critical feature-edit radius **2**, challenger Rhodanobacter.
- GTDB anchor, first switch `f80c0a86f49e5740576bda5302b82d72` (71,255 reads), second `a5a401d92657b7ac6e4d96863ad2e1d7` (6,167 reads); worst pair margin at B=1: +3,874; at B=2: -2,293 pooled reads.
- The analogous SILVA-family-block partition gives positive genus top-1 margin even under unrestricted block switches; this is only a model-conditional scenario.

## Distinction between verification levels

- The **original paired dataset** was processed and numerically analyzed remotely in the network-enabled runtime, with actual output verified as shown above.
- The **packaged local script** was syntax-checked and its **1,320 exact enumeration tests passed**, but the full public 21 MB QZV was not mounted in the local container. Full independent rerun from ZIP requires retrieving the source QZV.
- No per-sample feature table or phenotype contrast was used in this external experiment; do not call it biological clinical/phenotypic validation.
- The QIIME visualization and its source sequence may have reuse restrictions; the original file is not redistributed in this package.
