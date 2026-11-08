"""Stage 25: exact additive contrast certificates and sufficient top-k certificates.
Same 770 features, 34 samples, two observed reference assignments.
For each feature, both reference labels switch together for all its reads and all samples.
Independent per-feature worlds are a sensitivity scenario, not empirically observed worlds.
"""
import io,zipfile,tempfile,json
from pathlib import Path
import pandas as pd,numpy as np,h5py
from scipy.sparse import csr_matrix
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent / 'inputs'
OUT=Path(__file__).resolve().parent
with zipfile.ZipFile(ROOT/'Biological_ICDE_TAXONOMY_RETURN.zip') as zz:
 tax=[]
 for f in ['taxonomy-silva-138.qza','taxonomy-greengenes2-2024.09.qza']:
  with zipfile.ZipFile(io.BytesIO(zz.read(f))) as z:
   f2=next(x for x in z.namelist() if x.endswith('/data/taxonomy.tsv'))
   tax.append(pd.read_csv(io.BytesIO(z.read(f2)),sep='\t').set_index('Feature ID'))
with zipfile.ZipFile(ROOT/'table.qza') as z:
 f=next(x for x in z.namelist() if x.endswith('/data/feature-table.biom'))
 with tempfile.TemporaryDirectory() as td:
  p=Path(td)/'table.h5';p.write_bytes(z.read(f))
  with h5py.File(p) as h:
   features=[x.decode() for x in h['observation/ids'][:]]
   samples=[x.decode() for x in h['sample/ids'][:]]
   m=h['observation/matrix'];counts=csr_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=(len(features),len(samples))).toarray()
with zipfile.ZipFile(ROOT/'sample-metadata-viz.qzv') as z:
 f=next(x for x in z.namelist() if x.endswith('/data/metadata.tsv'))
 meta=pd.read_csv(io.BytesIO(z.read(f)),sep='\t',comment='#').set_index('sample-id').loc[samples]
assert len(features)==770 and len(samples)==34 and set(features)==set(tax[0].index)==set(tax[1].index)
X=counts/counts.sum(axis=0,keepdims=True)
contrasts=[('gut','tongue'),('left palm','right palm'),('gut','left palm')]
W=[]
for a,b in contrasts:
 w=np.zeros(len(samples))
 for site,sign in [(a,1),(b,-1)]:
  subjects=[s for s in sorted(meta.subject.unique()) if ((meta['body-site']==site)&(meta.subject==s)).any()]
  for subj in subjects:
   ix=np.flatnonzero(((meta['body-site']==site)&(meta.subject==subj)).to_numpy());w[ix]+=sign/(len(subjects)*len(ix))
 W.append(X@w)
W=np.array(W)
def label(v,k):
 ps=str(v).split(';');return ps[k].strip().split('__',1)[-1].strip() if k<len(ps) else ''
indices={'phylum':1,'family':4,'genus':5}
labels={}
for rank,k in indices.items():
 a=np.array([label(v,k) for v in tax[0].loc[features,'Taxon']]);b=np.array([label(v,k) for v in tax[1].loc[features,'Taxon']])
 if rank=='phylum':
  a=np.where(a=='Proteobacteria','Pseudomonadota',a)
  b=np.where(b=='Proteobacteria','Pseudomonadota',b)
 labels[rank]=(a,b)
rows=[];rankrows=[];toprows=[]
for ci,(aa,bb) in enumerate(contrasts):
 w=W[ci]; cn=f'{aa} - {bb}'
 for rank,(la,lb) in labels.items():
  names=sorted((set(la)|set(lb))-{''});group=[]
  for name in names:
   c0=np.where(la==name,w,0.);c1=np.where(lb==name,w,0.)
   v0=float(c0.sum());v1=float(c1.sum());low=float(np.minimum(c0,c1).sum());high=float(np.maximum(c0,c1).sum())
   obs_robust=(v0>1e-12 and v1>1e-12) or (v0< -1e-12 and v1< -1e-12)
   cert=(low>1e-12 or high< -1e-12)
   group.append({'contrast':cn,'rank':rank,'taxon':name,'silva_effect':v0,'gg2_effect':v1,'feature_L':low,'feature_U':high,'observed_same_nonzero_sign':obs_robust,'feature_sign_certified':cert,'observed_agreement_but_uncertified':bool(obs_robust and not cert),'width':high-low})
  rows.extend(group)
  g=pd.DataFrame(group)
  rankrows.append({'contrast':cn,'rank':rank,'queries':len(g),'observed_same_nonzero_sign':int(g.observed_same_nonzero_sign.sum()),'feature_sign_certified':int(g.feature_sign_certified.sum()),'agreement_but_uncertified':int(g.observed_agreement_but_uncertified.sum())})
  # Top-k ranking of taxon effect (not abundance) separately for each observed reference, and conservative sufficient certificate:
  # for every chosen q, L_q > max U_notchosen. Such pairwise interval test is valid, but may be conservative because comparisons share worlds.
  for k in (1,3,5):
   if len(g)<=k:continue
   s0=set(g.nlargest(k,'silva_effect').taxon);s1=set(g.nlargest(k,'gg2_effect').taxon)
   same=s0==s1
   if same:
    chosen=g[g.taxon.isin(s0)];others=g[~g.taxon.isin(s0)]
    margin=float(chosen.feature_L.min()-others.feature_U.max())
    certified=margin>1e-12
   else:margin=None;certified=False
   toprows.append({'contrast':cn,'rank':rank,'k':k,'observed_topk_same':same,'feature_topk_certified_sufficient':certified,'sufficient_margin':margin,'silva_topk':'|'.join(sorted(s0)),'gg2_topk':'|'.join(sorted(s1))})
q=pd.DataFrame(rows);summary=pd.DataFrame(rankrows);top=pd.DataFrame(toprows)
assert ((q.feature_L<=q[['silva_effect','gg2_effect']].min(axis=1)+1e-10)&(q.feature_U>=q[['silva_effect','gg2_effect']].max(axis=1)-1e-10)).all()
assert not ((q.feature_sign_certified)&(~q.observed_same_nonzero_sign)).any()
q.to_csv(OUT/'contrast_certificates.csv',index=False)
summary.to_csv(OUT/'contrast_summary.csv',index=False)
top.to_csv(OUT/'topk_sufficient_certificates.csv',index=False)
ex=q[q.observed_agreement_but_uncertified].copy();ex['abs_observed_effect_min']=ex[['silva_effect','gg2_effect']].abs().min(axis=1);ex=ex.sort_values('abs_observed_effect_min',ascending=False)
ex.to_csv(OUT/'agreement_failure_examples.csv',index=False)
fig,ax=plt.subplots(figsize=(10,5));idx=np.arange(len(summary));ax.bar(idx-.2,summary.observed_same_nonzero_sign,width=.4,label='Same sign, 2 observed references');ax.bar(idx+.2,summary.feature_sign_certified,width=.4,label='Certified across per-feature switches');ax.set_xticks(idx,[f'{r.rank}\n{r.contrast}' for r in summary.itertuples()],rotation=25,ha='right',fontsize=8);ax.set_ylabel('Number of queries');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'fig25_contrast_certification.png',dpi=180);plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4));tp=top.groupby('k')[['observed_topk_same','feature_topk_certified_sufficient']].sum();ix=np.arange(len(tp));ax.bar(ix-.2,tp.observed_topk_same,.4,label='Same top-k in 2 observed references');ax.bar(ix+.2,tp.feature_topk_certified_sufficient,.4,label='Certified (sufficient condition)');ax.set_xticks(ix,[f'Top-{k}' for k in tp.index]);ax.set_ylabel('Number of contrast × rank settings');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(OUT/'fig25_topk.png',dpi=180);plt.close(fig)
manifest={'features':len(features),'samples':len(samples),'reads':int(counts.sum()),'total_contrast_queries':len(q),'observed_agreement_but_uncertified':len(ex),'topk_scenarios':len(top),'topk_observed_same':int(top.observed_topk_same.sum()),'topk_certified_sufficient':int(top.feature_topk_certified_sufficient.sum()),'caveats':['Per-feature alternative worlds are model scenarios, not independently observed biological reality','No biological synonym matching apart from explicit Proteobacteria/Pseudomonadota canonicalization','Top-k interval separation sufficient but not necessary; failures do not prove top-k can change','Subject-balanced within-site contrast; not a significance test','Taxon label comparisons literal; unclassified placeholders may occur']}
(OUT/'validation.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
print(json.dumps(manifest,indent=2));print(summary.to_string(index=False));print(top.groupby('k')[['observed_topk_same','feature_topk_certified_sufficient']].sum().to_string())
