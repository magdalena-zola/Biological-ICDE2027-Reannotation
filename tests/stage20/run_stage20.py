"""Comparison of eager vs lazy cached tree DP vs full rebuild, identical semantics.
All dependencies and contributions are synthetic; no empirical taxonomic edge claims.
"""
import csv, json, time
from pathlib import Path
import numpy as np
from hierarchical_certificates import TreeCertificate
ROOT=Path(__file__).resolve().parent
ALLOWED=np.array([[True,True],[False,True]])
class LazyTree(TreeCertificate):
    def __init__(self, parents, allowed, contribution):
        super().__init__(parents,allowed,contribution)
        self.dirty=np.zeros(self.n,dtype=bool)
    def edit_lazy(self,v,x):
        self.c[v]=x
        while v>=0:
            self.dirty[v]=True
            v=self.parent[v]
    def query(self):
        # Only dirty nodes, bottom-up; their children already up-to-date.
        for v in reversed(np.flatnonzero(self.dirty)):
            self._calc(int(v))
        self.dirty[:]=False
        return self.certificate()

def experiment(n,q,query_every,seed):
    rng=np.random.default_rng(seed)
    par=[-1]+[(v-1)//2 for v in range(1,n)]
    c=rng.normal(size=(n,2,q))*0.01
    edits=[(int(rng.integers(n)),rng.normal(size=(2,q))*0.01) for _ in range(120)]
    eager=TreeCertificate(par,ALLOWED,c)
    lazy=LazyTree(par,ALLOWED,c)
    full=TreeCertificate(par,ALLOWED,c)
    times={'eager':0,'lazy':0,'full':0}; maxerr=0; checks=0
    for i,(v,new) in enumerate(edits):
        t=time.perf_counter_ns();eager.edit(v,new)
        if i%query_every==query_every-1 or i==len(edits)-1: ea=eager.certificate()
        times['eager']+=time.perf_counter_ns()-t
        t=time.perf_counter_ns();lazy.edit_lazy(v,new)
        if i%query_every==query_every-1 or i==len(edits)-1: la=lazy.query()
        times['lazy']+=time.perf_counter_ns()-t
        t=time.perf_counter_ns();full.c[v]=new
        if i%query_every==query_every-1 or i==len(edits)-1:
            for j in reversed(range(n)):full._calc(j)
            fu=full.certificate()
        times['full']+=time.perf_counter_ns()-t
        if i%query_every==query_every-1 or i==len(edits)-1:
            checks+=1
            for result in (ea,la):
                maxerr=max(maxerr,*[float(np.max(abs(x-y))) for x,y in zip(result,fu)])
            assert maxerr<1e-10,(n,q,query_every,maxerr)
    return {'nodes':n,'queries':q,'query_every_updates':query_every,'events':len(edits),'certificate_reads':checks,'eager_total_ms':times['eager']/1e6,'lazy_total_ms':times['lazy']/1e6,'full_total_ms':times['full']/1e6,'lazy_vs_eager':times['eager']/times['lazy'],'max_abs_error':maxerr}

if __name__=='__main__':
    rows=[]
    for n in (255,1023):
        for query_every in (1,5,20):
            for repeat in range(3):
                r=experiment(n,24,query_every,20261031+repeat)
                r['repeat']=repeat
                rows.append(r)
                print('PASS',n,query_every,repeat,'lazy/eager',round(r['lazy_vs_eager'],3),'error',r['max_abs_error'],flush=True)
    with (ROOT/'same_semantics_eager_lazy_full.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    s=[]
    for n in (255,1023):
        for freq in (1,5,20):
            a=[r for r in rows if r['nodes']==n and r['query_every_updates']==freq]
            s.append({'nodes':n,'query_every_updates':freq,'median_eager_ms':float(np.median([x['eager_total_ms'] for x in a])),'median_lazy_ms':float(np.median([x['lazy_total_ms'] for x in a])),'median_full_ms':float(np.median([x['full_total_ms'] for x in a])),'max_abs_error':max(x['max_abs_error'] for x in a)})
    with (ROOT/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(s[0]));w.writeheader();w.writerows(s)
    (ROOT/'validation.json').write_text(json.dumps({'runs':len(rows),'all_passed':True,'synthetic_hierarchies':True,'same_state_semantics':True,'summary':s},indent=2))
