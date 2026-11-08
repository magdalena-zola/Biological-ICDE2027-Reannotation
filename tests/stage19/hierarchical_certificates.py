"""Exact min/max additive certificates on a tree of coupled annotation states.
Demo uses synthetic constraints, NOT inferred biological dependencies.
Dependencies: Python 3.9+, numpy. Reproduce: python hierarchical_certificates.py
"""
import json, time, itertools, csv
from pathlib import Path
import numpy as np

class TreeCertificate:
    """Each node v has a state s in [0,K); allowed[parent_state,child_state] is Boolean.
    Query score is sum_v contribution[v,state,query]. Constraints only on edges.
    """
    def __init__(self, parents, allowed, contribution):
        self.parent = list(parents)
        self.n = len(parents)
        self.children = [[] for _ in range(self.n)]
        for v,p in enumerate(parents):
            if p >= 0: self.children[p].append(v)
        assert self.parent[0] == -1 and all(0 <= p < v for v,p in enumerate(parents) if v)
        self.allowed = np.asarray(allowed, bool)
        self.c = np.asarray(contribution, dtype=float).copy()
        self.K = self.c.shape[1]; self.Q = self.c.shape[2]
        assert self.c.shape[0] == self.n and self.allowed.shape == (self.K,self.K)
        self.lo = np.zeros_like(self.c); self.hi = np.zeros_like(self.c)
        for v in reversed(range(self.n)): self._calc(v)
    def _calc(self,v):
        low = self.c[v].copy(); high = low.copy()
        for ch in self.children[v]:
            for s in range(self.K):
                opts = self.allowed[s]
                if not opts.any():
                    low[s] = np.inf; high[s] = -np.inf
                else:
                    low[s] += self.lo[ch,opts].min(axis=0)
                    high[s] += self.hi[ch,opts].max(axis=0)
        self.lo[v]=low; self.hi[v]=high
    def certificate(self):
        return np.min(self.lo[0],axis=0),np.max(self.hi[0],axis=0)
    def edit(self,v,new_contribution):
        self.c[v] = new_contribution
        path = []
        while v >= 0:
            path.append(v);self._calc(v);v=self.parent[v]
        return len(path)
    def brute(self):
        mn=np.full(self.Q,np.inf);mx=np.full(self.Q,-np.inf);valid=0
        for states in itertools.product(range(self.K),repeat=self.n):
            if any(not self.allowed[states[self.parent[v]],states[v]] for v in range(1,self.n)):continue
            val = sum(self.c[v,states[v]] for v in range(self.n))
            mn=np.minimum(mn,val);mx=np.maximum(mx,val);valid+=1
        return mn,mx,valid

def main():
    rng=np.random.default_rng(20261031)
    out=Path(__file__).parent
    tests=0;maxerr=0.;paths=[];bench=[]
    # Rich but not vacuous constraint: parent=0 permits {0,1}, parent=1 only {1}
    allow=np.array([[1,1],[0,1]],dtype=bool)
    for n in (4,6,8,10):
        for trial in range(15):
            par=[-1]+[(i-1)//2 for i in range(1,n)]
            c=rng.normal(size=(n,2,7))
            t=TreeCertificate(par,allow,c)
            a,b,valid=t.brute();x,y=t.certificate()
            maxerr=max(maxerr,float(np.max(abs(x-a))),float(np.max(abs(y-b))))
            assert np.allclose(x,a,atol=1e-10) and np.allclose(y,b,atol=1e-10)
            v=int(rng.integers(n));new=rng.normal(size=(2,7));paths.append(t.edit(v,new))
            a,b,valid=t.brute();x,y=t.certificate()
            maxerr=max(maxerr,float(np.max(abs(x-a))),float(np.max(abs(y-b))))
            assert np.allclose(x,a,atol=1e-10) and np.allclose(y,b,atol=1e-10)
            tests+=2
    # Fairer tree baseline: same DP, full rebuild vs updated leaf-to-root path.
    for n in (63,255,1023):
        q=48;par=[-1]+[(i-1)//2 for i in range(1,n)]
        base=rng.normal(size=(n,2,q));edits=[(int(rng.integers(n)),rng.normal(size=(2,q))) for _ in range(25)]
        t=TreeCertificate(par,allow,base)
        full=TreeCertificate(par,allow,base)
        ti=[];tf=[];depths=[]
        for v,new in edits:
            st=time.perf_counter_ns();depths.append(t.edit(v,new));ti.append((time.perf_counter_ns()-st)/1e6)
            st=time.perf_counter_ns();full.c[v]=new
            for node in reversed(range(n)):full._calc(node)
            tf.append((time.perf_counter_ns()-st)/1e6)
            assert np.allclose(t.certificate()[0],full.certificate()[0],atol=1e-10)
            assert np.allclose(t.certificate()[1],full.certificate()[1],atol=1e-10)
        bench.append({'nodes':n,'queries':q,'edits':len(edits),'median_path_ms':round(float(np.median(ti)),5),'median_full_ms':round(float(np.median(tf)),5),'ratio_full_over_path':round(float(np.median(tf)/np.median(ti)),3),'median_path_nodes':float(np.median(depths))})
    with (out/'benchmark_stage17.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(bench[0]));w.writeheader();w.writerows(bench)
    report={'test_cases':tests,'max_abs_error':maxerr,'model':'synthetic parent-child compatibility, 2 states; additive queries, fixed denominators','validity':'exact conditional on modeled tree structure, not empirically validated taxonomy constraints','benchmarks':bench}
    (out/'validation_stage17.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
