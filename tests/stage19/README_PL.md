# Biological ICDE 2027 — Stage 19: Same-semantics tree certificates

## Findings
The exact binary-state, edge-constrained tree DP and its leaf-to-root update scheme are compared against full recomputation of **the same DP** on **the same sequence of updates**. 80 small synthetic tests check min and max against exhaustive enumeration. 180 updates (60 for each tree size) are used for benchmarking 48 additive queries.

These are **synthetic dependency trees**, not observed or independently validated biological reannotation dependencies. No DuckDB/DBToaster/F-IVM/DBSP native incremental benchmark is claimed.

## Mathematical statement
Let X be the finite Cartesian product of per-node state spaces, F a nonempty subset defined by edge compatibility, and f any fixed scalar-valued query. Define L_ind=min_{x in X} f(x), U_ind=max_{x in X} f(x), L_tree=min_{x in F}f(x), U_tree=max_{x in F}f(x). Then L_ind <= L_tree <= U_tree <= U_ind. The first inequality is strict precisely if no minimizer of f over X belongs to F; the final inequality is strict precisely if no maximizer over X belongs to F. The theorem is not specific to additive f; the implemented efficient DP assumes node-additive contributions and tree-factorized constraints.

For each node v and state s, lower value is c(v,s)+ sum over child u of min_{t compatible(s,t)} lower(u,t), upper analogously. Empty feasible sets use +/- infinity and require infeasibility detection at root. In tested model states are globally feasible.

## Results
Run `python run_stage19.py` with NumPy and `hierarchical_certificates.py` from stage 17 on the adjacent parent folder (included in this ZIP under `Biological_ICDE_stage17`). The benchmark uses wall-clock `perf_counter_ns`; results are hardware- and runtime-dependent, not directly comparable to Windows DuckDB benchmark.

## Important novelty caveat
Path-to-root maintenance and tree DP are standard algorithmic strategies. Published F-IVM (Kara et al., VLDB Journal 2024, DOI 10.1007/s00778-023-00817-w) supports factorized incremental maintenance; DBSP (Budiu et al., VLDB Journal 2025, DOI 10.1007/s00778-025-00922-y) is another strong incremental approach. This prototype is NOT a head-to-head implementation of either. A paper claiming general IVM novelty is unjustified.

## Proposed next empirical milestone
Construct admissibility constraints from independently verified taxonomic crosswalks or multiple real reference releases; test nontrivial certification gains on matched read tables. Implement a strong incremental baseline at the same semantics (at minimum memoized subtrees with dependency tracking, and a factorized relational implementation) and benchmark varying edit batch sizes, query selectivity, depth, and constraint density. Add complete memory and construction-time costs.
