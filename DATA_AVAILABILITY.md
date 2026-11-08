# Data availability and provenance

## Primary paired sample-level data

The feature count matrix and metadata are taken from the QIIME 2 *Moving Pictures* educational dataset: https://docs.qiime2.org/2024.10/tutorials/moving-pictures/ . The two referenced taxonomic annotations were generated from the same 770 identified representative sequence features using SILVA 138 and Greengenes2 2024.09. QIIME artifacts stored in `work/inputs/` are included in this **publication candidate**, but the authors must confirm the distribution rights of the original tutorial materials and classifications before making the repository public.

The derived results are exact for the stated additive query definitions and binary/blocked switching models, not evidence that individual-classifier or family switches are biologically likely.

## Independent paired-annotation experiment

External source: https://forum.qiime2.org/t/using-multiple-taxonomy-classifiers-for-comparison-or-in-filtering-workflows/34349

External dataset DOI: https://doi.org/10.5281/zenodo.15390940

External paired QIIME 2 visualization (SILVA 144 vs GTDB r232):
`https://www.dropbox.com/scl/fi/dycygxvmvjzrn0qdgiqvs/asv-seqs-gtdb-silva-fixed.qzv?rlkey=yemtsq1s5qprs4e808878z73s&dl=1`

Expected SHA-256: `e8aa584032374072c3b8607a870fadd179bd9c6e5198d37f605c71fcddb03cfb`.

The external binary is not redistributed in this archive; download the original artifact and check its hash. Its publicly exposed data are per-feature **pooled frequencies**: no sample-by-feature table was used to evaluate contrasts or cross-person uncertainty. No clinical, causal or subject-level conclusion is claimed.

## Public-release checklist

Before uploading this candidate, verify the right to redistribute every original/derived `.qza`, `.qzv`, and classifier output, the contributors' code rights, citation requirements, and any data licenses. If some artifacts cannot be redistributed, replace them with lawful source links and checksummed download scripts, then reassess whether the paper's declared submission category still satisfies reproducibility requirements.
