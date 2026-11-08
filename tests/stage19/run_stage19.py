"""Stage 19: conditional exactness, equal-semantics tree-DP benchmark.
Synthetic constraint trees; NOT evidence for biological coupling. Python + numpy.
"""
import csv,json,time,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from hierarchical_certificates import TreeCertificate

rng=np.random.default_rng(20261031)
allowed=np.array([[True,True],[False,True]])
# Since any root assignment is permitted, all-0 and all-1 satisfy these constraints;
# nevertheless mixed optima need not be feasible.

def independent(c):
    return c.min(axis=1).sum(axis=0),c.max(axis=1).sum(axis=0)

# counterexample: root=0 child=1 is allowed, root=1 child=0 is prohibited.
# Min per-node contributions favors root 1 (cost 0) and child 0 (cost 0) but cannot co-occur.
counter=np.array([[[1.],[0.]],[[0.],[1.]]])
obj=TreeCertificate([-1,0],allowed,counter)
lo,hi=obj.certificate();ilo,ihi=independent(counter)
assert np.isclose(ilo[0],0) and np.isclose(lo[0],1)
assert np.isclose(ihi[0],2) and np.isclose(hi[0],2)
# Exhaustive equivalence and bound monotonicity
max_error=0.;strict=0;count=0
for n in (4,6,8,10):
    for run in range(20):
        par=[-1]+[(i-1)//2 for i in range(1,n)]
        c=rng.normal(size=(n,2,8))
        obj=TreeCertificate(par,allowed,c)
        lower,upper=obj.certificate();a,b,num=obj.brute();il,iu=independent(c)
        err=max(float(np.max(abs(lower-a))),float(np.max(abs(upper-b))))
        assert err < 1e-10 and np.all(lower >= il-1e-10) and np.all(upper <= iu+1e-10)
        strict+=int(np.any(lower > il+1e-9) or np.any(upper < iu-1e-9))
        max_error=max(max_error,err);count+=1

# Identical event stream, identical feasibility model, identical contribution tensor.
bench=[]
for n in (63,255,1023):
    q=48;par=[-1]+[(v-1)//2 for v in range(1,n)]
    c=rng.normal(size=(n,2,q))*0.01
    events=[(int(rng.integers(0,n)),rng.normal(size=(2,q))*0.01) for i in range(60)]
    incremental=TreeCertificate(par,allowed,c)
    full=TreeCertificate(par,allowed,c)
    inc_times=[];full_times=[];touched=[]
    for v,new in events:
        t0=time.perf_counter_ns();touched.append(incremental.edit(v,new)); inc_times.append((time.perf_counter_ns()-t0)/1e6)
        t0=time.perf_counter_ns();full.c[v]=new
        for x in reversed(range(n)):full._calc(x)
        full_times.append((time.perf_counter_ns()-t0)/1e6)
        il,iu=incremental.certificate();fl,fu=full.certificate()
        assert np.allclose(il,fl,atol=1e-10) and np.allclose(iu,fu,atol=1e-10)
    bench.append({'nodes':n,'queries':q,'events':len(events),'incremental_median_ms':float(np.median(inc_times)),'full_dp_median_ms':float(np.median(full_times)),'full_over_incremental':float(np.median(full_times)/np.median(inc_times)),'median_recomputed_nodes':float(np.median(touched)),'full_recomputed_nodes':n})
with (ROOT/'same_semantics_benchmark.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(bench[0]));w.writeheader();w.writerows(bench)
report={'validation_cases':count,'cases_with_strictly_tighter_bounds':strict,'max_abs_error':max_error,'small_counterexample':{'independent_lower':float(ilo[0]),'tree_lower':float(lo[0]),'independent_upper':float(ihi[0]),'tree_upper':float(hi[0])},'model':'synthetic binary pairwise state compatibility; fixed additive query contribution; not biological dependency evidence','baseline':'Full tree DP, same constraints and query vectors, not database SQL','benchmarks':bench}
(ROOT/'validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
