#!/usr/bin/env python3
"""Exact tree-coherent edit budget for additive two-reference pair differences.
This is standard min-plus tree dynamic programming specialized to taxonomy worlds;
its value is a query contract, not a claim to invent generic tree DP.
"""
from pathlib import Path
import csv, json, time, itertools
import numpy as np
import correlated_certificates as core

OUT=Path(__file__).resolve().parent/'results'


def tree_min_pair(d0,d1,children,leaves,budget):
    """Rooted tree switches: selecting any node toggles reference (S/G) for
    all its leaf descendants; at most budget node toggles across the tree.
    Incoming ref b and outgoing ref b XOR toggle. Return exact min.
    """
    budget=int(budget); n=len(children)
    INF=np.inf;DP=[None]*n
    def conv(a,b):
        c=np.full(budget+1,INF)
        for i in range(budget+1):
            if not np.isfinite(a[i]):continue
            c[i:]=np.minimum(c[i:], a[i]+b[:budget+1-i])
        return c
    for v in reversed(range(n)):
        arr=np.full((2,budget+1),INF)
        if len(children[v])==0:
            f=leaves[v]
            for incoming in (0,1):
                arr[incoming,0]=d0[f] if incoming==0 else d1[f]
                if budget>=1:arr[incoming,1]=d1[f] if incoming==0 else d0[f]
        else:
            for incoming in (0,1):
                for event in (0,1):
                    if event>budget:continue
                    parity=incoming^event
                    acc=np.full(budget+1,INF);acc[0]=0
                    for child in children[v]:
                        acc=conv(acc,DP[child][parity])
                    arr[incoming,event:]=np.minimum(arr[incoming,event:],acc[:budget+1-event])
        DP[v]=arr
    return min(np.min(DP[0][0]),np.min(DP[0][1]))


def build_tax_tree(feat):
    # Nodes constructed preorder, children indexed larger than parent.
    # Every node is a permissible binary-toggle operation. Group names
    # depend only on SILVA's taxonomy, which is an explicit modeling choice.
    import io,zipfile,pandas as pd
    from pathlib import Path
    f=core.SRC/'Biological_ICDE_TAXONOMY_RETURN.zip'
    with zipfile.ZipFile(f) as zp:
        with zipfile.ZipFile(io.BytesIO(zp.read('taxonomy-silva-138.qza'))) as z:
            fn=next(n for n in z.namelist() if n.endswith('/data/taxonomy.tsv'))
            tax=pd.read_csv(io.BytesIO(z.read(fn)),sep='\t').set_index('Feature ID')
    keys=['root']
    parents=[-1]
    leaf={}
    index={():0}
    for i,feature in enumerate(feat):
        parts=[x.strip().split('__',1)[-1].strip() for x in tax.loc[feature,'Taxon'].split(';')]
        par=0;prefix=[]
        for ix in (1,4,5):
            lbl=parts[ix] if ix<len(parts) and parts[ix] else 'undefined:'+feature
            prefix.append((ix,lbl))
            tup=tuple(prefix)
            if tup not in index:
                index[tup]=len(keys);keys.append(f'rank{ix}:{lbl}');parents.append(par)
            par=index[tup]
        k=len(keys);keys.append('feature:'+feature);parents.append(par);leaf[k]=i
    child=[[] for _ in keys]
    for j,p in enumerate(parents):
        if p>=0:child[p].append(j)
    assert len(leaf)==len(feat)
    return child,leaf,keys


def random_trees_check():
    rng=np.random.default_rng(5062)
    done=0;err=0
    for N in range(2,8):
        for repetition in range(60):
            # A rooted random tree with N-1 leaves (some internal degree 1 nodes possible).
            child=[[] for _ in range(N)]
            for v in range(1,N):child[int(rng.integers(0,v))].append(v)
            leaves={v:i for i,v in enumerate(x for x in range(N) if not child[x])}
            d0=rng.normal(size=len(leaves));d1=rng.normal(size=len(leaves))
            for budget in range(0,min(4,N+1)):
                dp=tree_min_pair(d0,d1,child,leaves,budget)
                vals=[]
                for anchor in (0,1):
                 for selections in itertools.product([0,1],repeat=N):
                    if sum(selections)>budget:continue
                    s=0.
                    for node,f in leaves.items():
                        parity=anchor;v=node
                        # retrieve parent by traversal for the small random tree
                        while v>=0:
                            parity^=selections[v]
                            if v==0:break
                            v=next(j for j,c in enumerate(child) if v in c)
                        s+=d1[f] if parity else d0[f]
                    vals.append(s)
                truth=min(vals)
                err=max(err,abs(dp-truth)); assert abs(dp-truth)<1e-10
                done+=1
    return {'tree_DP_brute_cases':done,'max_abs_error':float(err)}


def run_real():
    feat,labels,groups,W=core.input_real()
    tree,leaves,keys=build_tax_tree(feat)
    results=[]
    for ci,(a,b) in enumerate(core.CONTRASTS):
        cname=f'{a} - {b}'
        if cname!='left palm - right palm':continue
        w=W[ci]
        for rank,taxa in labels.items():
            la,lb=taxa; names=sorted((set(la)|set(lb))-{''})
            if rank=='phylum':
                pairs=[('Pseudomonadota','Firmicutes')]
            elif rank=='genus':
                pairs=[('Corynebacterium','__ZERO__'),('Prevotella','__ZERO__')]
            elif rank=='family':pairs=[('Aerococcaceae','__ZERO__')]
            for q,r in pairs:
                f0=w*(la==q);f1=w*(lb==q)
                if r!='__ZERO__': f0=f0-w*(la==r);f1=f1-w*(lb==r)
                # choose direction where the observed sign agreement is positive,
                # for sign tests q may be robust-negative.
                orient=1 if (sum(f0)+sum(f1))>0 else -1
                if r=='__ZERO__':f0=f0*orient;f1=f1*orient
                for budget in (0,1,2,3,4):
                    start=time.perf_counter(); val=tree_min_pair(f0,f1,tree,leaves,budget)
                    ms=(time.perf_counter()-start)*1000
                    results.append({'query':q,'competitor':r,'rank':rank,'clade_event_budget':budget,
                                    'worst_margin':val,'strict_certified':bool(val>core.TOL),
                                    'milliseconds_python':ms, 'taxonomy_nodes':len(tree),
                                    'orientation':orient})
    import pandas as pd
    df=pd.DataFrame(results);df.to_csv(OUT/'hierarchical_clade_budget_checks.csv',index=False)
    return {'taxonomy_tree_nodes':len(tree),'hierarchical_real_queries':int(df[['query','competitor']].drop_duplicates().shape[0]),
            'real_DP_evaluations':len(df),'max_budget':4}

if __name__=='__main__':
    v=random_trees_check();r=run_real();v.update(r)
    (OUT/'hierarchy_validation.json').write_text(json.dumps(v,indent=2))
    print(json.dumps(v,indent=2))
