"""Independent post-stage25 tests: exact pairwise top-k stability and witnesses.

Requires only the publicly available QIIME2 Moving Pictures artifacts already in
our local work directory. The admissible mixed worlds are hypothetical.
"""
from pathlib import Path
import io, zipfile, json, tempfile, itertools
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csr_matrix
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parent/'results'
OUT.mkdir(parents=True,exist_ok=True)
EPS=1e-12

def get_tax(fn):
    with zipfile.ZipFile(ROOT/fn) as z:
        name=next(f for f in z.namelist() if f.endswith('/data/taxonomy.tsv'))
        return pd.read_csv(io.BytesIO(z.read(name)),sep='\t').set_index('Feature ID')
for f in ['taxonomy-silva-138.qza', 'taxonomy-greengenes2-2024.09.qza']:
    if not (ROOT/f).exists():
        raise FileNotFoundError(f'Copy public taxonomy input {f} into {ROOT}')
S,G=get_tax('taxonomy-silva-138.qza'),get_tax('taxonomy-greengenes2-2024.09.qza')
with zipfile.ZipFile(ROOT/'table.qza') as z:
    name=next(n for n in z.namelist() if n.endswith('/data/feature-table.biom'))
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)/'x.biom';p.write_bytes(z.read(name))
        with h5py.File(p) as h:
            features=[v.decode() for v in h['observation/ids'][:]]
            samples=[v.decode() for v in h['sample/ids'][:]]
            a=h['observation/matrix']
            counts=csr_matrix((a['data'][:],a['indices'][:],a['indptr'][:]),shape=(len(features),len(samples))).toarray()
with zipfile.ZipFile(ROOT/'sample-metadata-viz.qzv') as z:
    name=next(n for n in z.namelist() if n.endswith('/data/metadata.tsv'))
    meta=pd.read_csv(io.BytesIO(z.read(name)),sep='\t',comment='#').set_index('sample-id').loc[samples]
assert set(S.index)==set(G.index)==set(features)
assert len(features)==770 and len(samples)==34 and counts.sum()==153807
normalized=counts/counts.sum(axis=0,keepdims=True)
contrasts=[('gut','tongue'),('left palm','right palm'),('gut','left palm')]
wts=[]
for a,b in contrasts:
    w=np.zeros(len(samples))
    for site,sign in [(a,1),(b,-1)]:
        subjects=[s for s in sorted(meta.subject.unique()) if ((meta['body-site']==site)&(meta.subject==s)).any()]
        for subj in subjects:
            ix=np.flatnonzero(((meta['body-site']==site)&(meta.subject==subj)).to_numpy())
            w[ix]+=sign/(len(subjects)*len(ix))
    wts.append(normalized@w)

def label(x,idx):
    parts=str(x).split(';')
    return parts[idx].strip().split('__',1)[-1].strip() if idx<len(parts) else ''

rank_idx={'phylum':1,'family':4,'genus':5}
labels={}
for rk,ix in rank_idx.items():
    x=np.array([label(v,ix) for v in S.loc[features,'Taxon']]); y=np.array([label(v,ix) for v in G.loc[features,'Taxon']])
    if rk=='phylum':
        x=np.where(x=='Proteobacteria','Pseudomonadota',x)
        y=np.where(y=='Proteobacteria','Pseudomonadota',y)
    labels[rk]=(x,y)

original=pd.read_csv(ROOT/'Biological_ICDE_stage25/topk_sufficient_certificates.csv').set_index(['contrast','rank','k'])
new=[]; counterexamples=[]; worst=0; pair_checks=0
for (a,b),w in zip(contrasts,wts):
  cn=f'{a} - {b}'
  for rk,(x,y) in labels.items():
    taxa=sorted((set(x)|set(y))-{''})
    # same scores and deterministic ordering as stage25
    zero=np.array([w[x==t].sum() for t in taxa]);one=np.array([w[y==t].sum() for t in taxa]);
    lo=np.array([np.minimum(np.where(x==t,w,0),np.where(y==t,w,0)).sum() for t in taxa]);
    hi=np.array([np.maximum(np.where(x==t,w,0),np.where(y==t,w,0)).sum() for t in taxa]);
    df=pd.DataFrame({'taxon':taxa,'silva_effect':zero,'gg2_effect':one,'L':lo,'U':hi})
    for k in (1,3,5):
      si=set(df.nlargest(k,'silva_effect').taxon)
      sj=set(df.nlargest(k,'gg2_effect').taxon)
      agree=si==sj
      old=original.loc[cn,rk,k]
      assert bool(old['observed_topk_same'])==agree
      if not agree:
        new.append(dict(contrast=cn,rank=rk,k=k,observed_agreement=False,interval_certified=False,exact_pairwise_certified=False,min_exact_pairwise_margin=np.nan,negative_witness=False,witness_margin=np.nan,witness_new_topk='',worst_pair=''))
        continue
      chosen=sorted(si); unchosen=sorted(set(taxa)-si)
      old_interval=bool(old['feature_topk_certified_sufficient'])
      min_int=min(lo[taxa.index(t)] for t in chosen)-max(hi[taxa.index(t)] for t in unchosen)
      assert old_interval==(min_int>EPS)
      pair_min=np.inf;worst_pair=None; witness_top=''; witness_margin=np.nan; negative=False
      for t in chosen:
        tx=(x==t);ty=(y==t)
        for u in unchosen:
          cS=w*(tx.astype(float)-(x==u).astype(float))
          cG=w*(ty.astype(float)-(y==u).astype(float))
          z=np.minimum(cS,cG).sum()
          pair_checks+=1
          if z<pair_min:
            pair_min=float(z);worst_pair=(t,u)
            # Per-feature switches independently minimize a single pair margin.
            choices=cS>cG # True means choose G
            choose=np.where(choices,y,x)
            actual_scores={tax:float(w[choose==tax].sum()) for tax in taxa}
            rank_order=sorted(taxa,key=lambda t:(-actual_scores[t],t))
            actual_top=set(rank_order[:k]); witness_top='|'.join(sorted(actual_top))
            witness_margin=actual_scores[t]-actual_scores[u]
            negative=(witness_margin < -EPS and actual_top!=si)
      certified=pair_min>EPS
      assert not old_interval or certified
      assert abs(witness_margin-pair_min)<1e-10
      new.append(dict(contrast=cn,rank=rk,k=k,observed_agreement=True,interval_certified=old_interval,exact_pairwise_certified=certified,min_exact_pairwise_margin=pair_min,negative_witness=negative,witness_margin=witness_margin,witness_new_topk=witness_top,worst_pair=' vs '.join(worst_pair)))
      if not certified:
        counterexamples.append({'contrast':cn,'rank':rk,'k':k,'worst_pair':worst_pair,'exact_margin':pair_min,'witness_negative_and_swaps':negative,'witness_topk':witness_top})
new=pd.DataFrame(new)
new.to_csv(OUT/'exact_pairwise_topk.csv',index=False)

# Independent finite exhaustive verification of mathematically exact additive bounds.
rng=np.random.default_rng(20261105)
tests=0; max_err=0.0; monotonic=0
for n in range(2,11):
    for rep in range(80):
        c=rng.normal(0,1,(n,2))
        all_states=np.array(list(itertools.product([0,1],repeat=n)),dtype=int)
        values=c[np.arange(n)[:,None], all_states.T].sum(axis=0)
        l0=c.min(axis=1).sum();u0=c.max(axis=1).sum()
        e=max(abs(l0-values.min()),abs(u0-values.max()))
        assert e<1e-10
        max_err=max(max_err,float(e));tests+=1
        # Random disjoint 2-4 feature blocks, all blocks choose one shared reference.
        z=rng.integers(0,max(2,n//2),n)
        buckets=list(dict.fromkeys(z))
        allowed=np.array([all((state[z==b]==state[z==b][0]).all() for b in buckets) for state in all_states])
        bvals=values[allowed]
        l=sum(min(c[z==b,0].sum(),c[z==b,1].sum()) for b in buckets)
        u=sum(max(c[z==b,0].sum(),c[z==b,1].sum()) for b in buckets)
        e=max(abs(l-bvals.min()),abs(u-bvals.max()))
        assert e<1e-10
        assert l>=l0-1e-10 and u<=u0+1e-10
        max_err=max(max_err,float(e));tests+=1;monotonic+=1

report={
 'execution_date':'2026-11-05',
 'real_input':{'features':len(features),'samples':len(samples),'reads':int(counts.sum()),'contrast_rank_settings':9},
 'exact_pairwise_topk':{'scenarios':len(new),'observed_agreement':int(new.observed_agreement.sum()),'interval_certified':int(new.interval_certified.sum()),'pairwise_certified':int(new.exact_pairwise_certified.sum()),'additional_pairwise_certifications':int((new.exact_pairwise_certified&~new.interval_certified).sum()),'observed_agreement_cases_with_pairwise_failure':int((new.observed_agreement&~new.exact_pairwise_certified).sum()),'negative_topk_witnesses':int(new.negative_witness.sum()),'pairwise_comparisons':pair_checks,'counterexamples':counterexamples},
 'exhaustive_finite_worlds':{'test_instances':tests,'max_abs_error':max_err,'block_monotonicity_checks':monotonic},
 'interpretation':'Pairwise minima are exact for the feature-independent additive world model. They do NOT establish biological likelihood of mixed worlds. Strict top-k certification assumes positive margins and no cross-boundary ties.'
}
(OUT/'additional_tests.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf8')
print(json.dumps(report,indent=2,ensure_ascii=False));print(new[['contrast','rank','k','observed_agreement','interval_certified','exact_pairwise_certified','min_exact_pairwise_margin','negative_witness']].to_string(index=False))
