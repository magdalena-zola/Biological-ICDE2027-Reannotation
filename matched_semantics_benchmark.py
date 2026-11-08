#!/usr/bin/env python3
"""Exact same-operator pair certificate maintenance, full vs affected-pair replay.
Synthetic aggregate payload updates, not a DBSP/F-IVM/DuckDB native benchmark.
"""
import json,time
from pathlib import Path
import numpy as np,pandas as pd
import correlated_certificates as c
OUT=Path(__file__).resolve().parent/'results'

def run_case(G=206,T=64,K=4,B=4,steps=120,repetitions=3):
  results=[]
  pairs=[(i,j) for i in range(K) for j in range(K,T)]
  includes={t:np.array([idx for idx,(i,j) in enumerate(pairs) if i==t or j==t],dtype=int) for t in range(T)}
  for seed in range(repetitions):
    rng=np.random.default_rng(26362+seed)
    A=rng.normal(0,.1,(G,T));D=rng.normal(0,.04,(G,T));C=A+D
    ops=[(int(rng.integers(G)),int(rng.integers(T)),float(rng.normal(0,.01)),float(rng.normal(0,.01))) for _ in range(steps)]
    def calc_pair(i,j,X,Y):return c.pair_margin_budget(X[:,i]-X[:,j],Y[:,i]-Y[:,j],B)
    def run_full():
      X=A.copy();Y=C.copy();last=np.zeros(len(pairs))
      start=time.perf_counter()
      for g,t,da,db in ops:
        X[g,t]+=da;Y[g,t]+=db
        for p,(i,j) in enumerate(pairs):last[p]=calc_pair(i,j,X,Y)
      return time.perf_counter()-start,last
    def run_indexed():
      X=A.copy();Y=C.copy()
      last=np.array([calc_pair(i,j,X,Y) for i,j in pairs])
      build=time.perf_counter()
      start=time.perf_counter()
      for g,t,da,db in ops:
        X[g,t]+=da;Y[g,t]+=db
        for p in includes[t]:
          i,j=pairs[int(p)];last[p]=calc_pair(i,j,X,Y)
      return time.perf_counter()-start,last
    tf,full=run_full(); ti,idx=run_indexed()
    err=float(np.max(np.abs(full-idx)));assert err<1e-10
    results.append({'groups':G,'taxa':T,'top_k':K,'budget':B,'updates':steps,'repetition':seed,
                    'comparisons':len(pairs),'baseline_full_ms':tf*1000,'indexed_ms':ti*1000,
                    'speedup':tf/ti,'max_abs_error':err,'execution_note':'Python; same exact budgeted pair-minimum operator; single taxon-block update'})
  return results

if __name__=='__main__':
  rows=run_case()
  df=pd.DataFrame(rows);df.to_csv(OUT/'same_operator_incremental_benchmark.csv',index=False)
  summary={
      'repetitions':len(df),'updates_per_repetition':120,
      'pairwise_comparisons_per_full_rebuild':int(df.comparisons.iloc[0]),
      'median_total_full_ms':float(df.baseline_full_ms.median()),
      'median_total_indexed_ms':float(df.indexed_ms.median()),
      'median_speedup_per_run':float(df.speedup.median()),
      'max_error':float(df.max_abs_error.max()),
      'limitations':'Synthetic aggregate payloads; no native DB engine or actual taxonomy reload; same Python pair operator only.'}
  (OUT/'same_operator_incremental_summary.json').write_text(json.dumps(summary,indent=2))
  print(json.dumps(summary,indent=2))
