#!/usr/bin/env python3
"""Two-anchored, budgeted, group-coherent query certificates.
Recomputes results using public QIIME2 Moving Pictures paired reference inputs.
Do not interpret anchored clade switching as an empirically calibrated biological process.
"""
from __future__ import annotations
import csv, io, itertools, json, math, os, statistics, tempfile, time, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csr_matrix

ROOT=Path(__file__).resolve().parent
SRC=ROOT/'work'/'inputs'
OUT=ROOT/'results'
OUT.mkdir(exist_ok=True,parents=True)
CONTRASTS=[('gut','tongue'),('left palm','right palm'),('gut','left palm')]
RANKS={'phylum':1,'family':4,'genus':5}
BUDGETS=[0,1,2,4,8,16,32,64,128,'all']
TOL=1e-11


def input_real():
    with zipfile.ZipFile(SRC/'Biological_ICDE_TAXONOMY_RETURN.zip') as z:
        tax=[]
        for f in ['taxonomy-silva-138.qza','taxonomy-greengenes2-2024.09.qza']:
            with zipfile.ZipFile(io.BytesIO(z.read(f))) as inner:
                name=next(x for x in inner.namelist() if x.endswith('/data/taxonomy.tsv'))
                tax.append(pd.read_csv(io.BytesIO(inner.read(name)),sep='\t').set_index('Feature ID'))
    with zipfile.ZipFile(SRC/'table.qza') as z:
        name=next(x for x in z.namelist() if x.endswith('/data/feature-table.biom'))
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'table.h5';p.write_bytes(z.read(name))
            with h5py.File(p) as h:
                feat=[v.decode() for v in h['observation/ids'][:]]
                sam=[v.decode() for v in h['sample/ids'][:]]
                m=h['observation/matrix']
                counts=csr_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=(len(feat),len(sam))).toarray()
    with zipfile.ZipFile(SRC/'sample-metadata-viz.qzv') as z:
        f=next(x for x in z.namelist() if x.endswith('/data/metadata.tsv'))
        meta=pd.read_csv(io.BytesIO(z.read(f)),sep='\t',comment='#').set_index('sample-id').loc[sam]
    assert len(feat)==770 and len(sam)==34 and int(counts.sum())==153807
    assert set(feat)==set(tax[0].index)==set(tax[1].index)
    X=counts/counts.sum(axis=0,keepdims=True)
    W=[]
    for aa,bb in CONTRASTS:
        w=np.zeros(len(sam))
        for site,sgn in [(aa,1),(bb,-1)]:
            subjects=[s for s in sorted(meta.subject.unique()) if ((meta['body-site']==site)&(meta.subject==s)).any()]
            for s in subjects:
                ix=np.flatnonzero(((meta['body-site']==site)&(meta.subject==s)).to_numpy())
                w[ix]+=sgn/(len(subjects)*len(ix))
        W.append(X@w)
    def label(v,k):
        parts=str(v).split(';')
        return parts[k].strip().split('__',1)[-1].strip() if len(parts)>k else ''
    labels={}
    for rank,k in RANKS.items():
        la=np.array([label(v,k) for v in tax[0].loc[feat,'Taxon']]); lb=np.array([label(v,k) for v in tax[1].loc[feat,'Taxon']])
        if rank=='phylum':
            la=np.where(la=='Proteobacteria','Pseudomonadota',la)
            lb=np.where(lb=='Proteobacteria','Pseudomonadota',lb)
        labels[rank]=(la,lb)
    silva_family=np.array([label(v,4) for v in tax[0].loc[feat,'Taxon']]);
    silva_phylum=np.array([label(v,1) for v in tax[0].loc[feat,'Taxon']]);
    # Missing labels are not merged into a giant spurious "unknown family" block.
    groups={'feature':[f'feature:{x}' for x in feat],
            'silva_family':[f'family:{s}' if s else f'unknown_feature:{i}' for i,s in enumerate(silva_family)],
            'silva_phylum':[f'phylum:{s}' if s else f'unknown_feature:{i}' for i,s in enumerate(silva_phylum)]}
    return feat,labels,groups,W


def grouped(A,B,items):
    labels,inv=np.unique(items,return_inverse=True)
    F,T=A.shape
    x=np.zeros((len(labels),T));y=x.copy()
    np.add.at(x,inv,A);np.add.at(y,inv,B)
    return x,y,labels


def pair_margin_budget(d0,d1,budget,witness=False):
    """Return exact worst margin over symmetric union of B-edit balls at 2 anchors.
    d0,d1 = blockwise contribution to q_i-q_j in S and GG2.
    Choice assignment in witness = 0/S, 1/G, one per block.
    """
    d0=np.asarray(d0,dtype=float);d1=np.asarray(d1,dtype=float)
    G=len(d0); b=G if budget=='all' else min(int(budget),G)
    delta=d1-d0
    neg=np.argsort(delta,kind='stable')[:b]
    takeS=neg[delta[neg]<0]
    vS=float(np.sum(d0)+delta[takeS].sum())
    pos=np.argsort(delta,kind='stable')[::-1][:b]
    takeG=pos[delta[pos]>0]
    vG=float(np.sum(d1)-delta[takeG].sum())
    if vS<=vG:
        val=vS; z=np.zeros(G,dtype=np.uint8); z[takeS]=1; anchor=0
    else:
        val=vG; z=np.ones(G,dtype=np.uint8); z[takeG]=0; anchor=1
    if witness:
        actual=float(np.where(z,d1,d0).sum())
        assert abs(actual-val)<1e-9
        return val,z,anchor
    return val


def fast_min_margin(d0,d1,brange):
    """Exact shared-world pair margin vector for multiple budgets using sorted deltas."""
    delta=d1-d0
    G=len(delta)
    sorted_delta=np.sort(delta)
    cneg=np.concatenate([[0.0],np.cumsum(np.minimum(sorted_delta,0.0))])
    sorted_pos=sorted_delta[::-1]
    cpos=np.concatenate([[0.0],np.cumsum(np.maximum(sorted_pos,0.0))])
    s0=float(d0.sum());s1=float(d1.sum())
    return {str(b):min(s0+cneg[min(int(b),G)] if b!='all' else s0+cneg[G],
                        s1-cpos[min(int(b),G)] if b!='all' else s1-cpos[G]) for b in brange}


def brute_tests():
    rng=np.random.default_rng(20261106)
    tests=0;maxerr=0.;monotone=0;witnesses=0
    for G in range(1,10):
        for trial in range(55):
            d0=rng.normal(size=G);d1=rng.normal(size=G)
            budgets=[0,1,2,3,G,'all']
            prev=1e100
            for b in budgets:
                k=G if b=='all' else min(b,G)
                v,z,anchor=pair_margin_budget(d0,d1,b,witness=True)
                vals=[]
                for bits in itertools.product((0,1),repeat=G):
                    n=sum(bits)
                    if min(n,G-n)<=k:vals.append(sum(d1[j] if bits[j] else d0[j] for j in range(G)))
                exact=min(vals)
                maxerr=max(maxerr,abs(v-exact));assert abs(v-exact)<3e-12
                assert v<=prev+3e-12
                prev=v;tests+=1;monotone+=1
                witnesses+=int(abs(np.where(z,d1,d0).sum()-exact)<1e-11)
    return {'exhaustive_cases':tests,'max_abs_error':maxerr,'budget_monotonic_checks':monotone,'witness_checks':witnesses}


def calc():
    feat,labels,groups,W=input_real()
    record=[];summary=[];rad_rows=[];witnesses=[];signw=[]
    orig_top=pd.read_csv(ROOT/'work'/'results'/'topk_sufficient_certificates.csv')
    original_summary=pd.read_csv(ROOT/'work'/'results'/'contrast_summary.csv')
    npairs=0
    for ci,(aa,bb) in enumerate(CONTRASTS):
        cn=f'{aa} - {bb}'
        w=W[ci]
        for rank,(la,lb) in labels.items():
            names=sorted((set(la)|set(lb))-{''})
            names=np.array(names)
            ca=(la[:,None]==names[None,:])*w[:,None]
            cb=(lb[:,None]==names[None,:])*w[:,None]
            assert len(names) in (34,197,468)
            orig=original_summary[(original_summary.contrast==cn)&(original_summary['rank']==rank)].iloc[0]
            a0=ca.sum(axis=0);a1=cb.sum(axis=0)
            same_sign=(a0>TOL)&(a1>TOL)|(a0< -TOL)&(a1< -TOL)
            assert int(same_sign.sum())==int(orig.observed_same_nonzero_sign)
            for partition,blocks in groups.items():
                A,B,blabels=grouped(ca,cb,blocks)
                G=len(blabels)
                # Sign radius for each query: first B under which a zero/opposite sign is possible
                for t in np.flatnonzero(same_sign):
                    sgn=1. if a0[t]>0 else -1.
                    res=fast_min_margin(A[:,t]*sgn,B[:,t]*sgn,range(G+1))
                    # W_B = union of at most B switches from either observed full reference.
                    radius=next((b for b in range(G+1) if res[str(b)]<=TOL),None)
                    rad_rows.append({'contrast':cn,'rank':rank,'partition':partition,'taxon':names[t],
                              'blocks':G,'observed_sign':int(sgn),'critical_edit_radius':radius if radius is not None else -1,
                              'critical_edit_fraction':(radius/G if radius is not None else np.nan),
                              'margin_b0':res['0'],'margin_all':res[str(G)]})
                    if radius is not None:
                        value,z,anchor=pair_margin_budget(A[:,t]*sgn,B[:,t]*sgn,radius,witness=True)
                        signw.append({'contrast':cn,'rank':rank,'partition':partition,'taxon':names[t],
                            'critical_edit_radius':radius,'oriented_margin_at_radius':value,
                            'anchor':'SILVA' if anchor==0 else 'Greengenes2',
                            'number_switched_blocks':int(np.sum(z!=anchor)),
                            'switched_block_labels':'|'.join(blabels[np.flatnonzero(z!=anchor)])})
                for k in (1,3,5):
                    s0=set(np.argsort(-a0,kind='stable')[:k]);s1=set(np.argsort(-a1,kind='stable')[:k]);same=s0==s1
                    # Exact labels may tie; verify no relevant ties at input top-k boundary.
                    assert not(np.isclose(np.sort(a0)[-k],np.sort(a0)[-k-1],atol=1e-13) and np.sort(a0)[-k]!=0)
                    assert not(np.isclose(np.sort(a1)[-k],np.sort(a1)[-k-1],atol=1e-13) and np.sort(a1)[-k]!=0)
                    if same:
                        pairs=[(i,j) for i in sorted(s0) for j in range(len(names)) if j not in s0]
                        margins={str(b):float('inf') for b in BUDGETS}
                        radius=None;critical_pair=None
                        for i,j in pairs:
                            d0=A[:,i]-A[:,j];d1=B[:,i]-B[:,j]
                            mm=fast_min_margin(d0,d1,range(G+1))
                            cand_r=next((b for b in range(G+1) if mm[str(b)]<=TOL),None)
                            if cand_r is not None and (radius is None or cand_r<radius):
                                radius=cand_r;critical_pair=(i,j)
                            for b in BUDGETS:
                                v=mm[str(min(int(b),G))] if b!='all' else mm[str(G)]
                                margins[str(b)]=min(margins[str(b)],v)
                            npairs+=1
                        for budget in BUDGETS:
                            m=margins[str(budget)]
                            record.append({'contrast':cn,'rank':rank,'partition':partition,'k':k,'blocks':G,'budget':budget,
                                          'two_observed_agree':True,'certified_strict':bool(m>TOL),'worst_margin':m,
                                          'critical_edit_radius':radius if radius is not None else -1})
                        if radius is not None and critical_pair is not None:
                            i,j=critical_pair
                            m,z,anchor=pair_margin_budget(A[:,i]-A[:,j],B[:,i]-B[:,j],radius,witness=True)
                            witnesses.append({'contrast':cn,'rank':rank,'partition':partition,'k':k,
                                'critical_edit_radius':radius,'chosen':names[i],'challenger':names[j],
                                'margin_at_radius':m,'anchor':'SILVA' if anchor==0 else 'Greengenes2',
                                'switch_count':int(np.sum(z!=anchor)),
                                'switched_block_labels':'|'.join(blabels[np.flatnonzero(z!=anchor)])})
                    else:
                        for budget in BUDGETS:
                            record.append({'contrast':cn,'rank':rank,'partition':partition,'k':k,'blocks':G,'budget':budget,
                                        'two_observed_agree':False,'certified_strict':False,'worst_margin':np.nan,
                                        'critical_edit_radius':0})
                summary.append({'contrast':cn,'rank':rank,'partition':partition,'blocks':G,
                                'observed_same_sign':int(same_sign.sum()),'taxa':len(names)})
    df=pd.DataFrame(record);df.to_csv(OUT/'budgeted_topk_results.csv',index=False)
    rr=pd.DataFrame(rad_rows);rr.to_csv(OUT/'sign_robustness_radii.csv',index=False)
    ww=pd.DataFrame(witnesses);ww.to_csv(OUT/'budgeted_topk_witnesses.csv',index=False)
    pd.DataFrame(signw).to_csv(OUT/'sign_counterexample_assignments.csv',index=False)
    pd.DataFrame(summary).to_csv(OUT/'input_summary.csv',index=False)
    assert df[df.budget==0].query('two_observed_agree').certified_strict.all()
    obs=int(df[(df.budget==0)&(df.partition=='feature')].two_observed_agree.sum())
    assert obs==11
    real_summary=[]
    for partition in groups:
        tmp=df[(df.partition==partition)&(df.two_observed_agree)]
        for b in BUDGETS:
            sub=tmp[tmp.budget.astype(str)==str(b)]
            real_summary.append({'partition':partition,'budget':str(b),'observed_agreement':11,'certified_topk':int(sub.certified_strict.sum()),
                                 'min_margin':float(sub.worst_margin.min()),'block_count':int(sub.blocks.iloc[0]) if partition=='feature' else None})
    pd.DataFrame(real_summary).to_csv(OUT/'certificate_frontiers.csv',index=False)
    integrity={'features':len(feat),'samples':34,'reads':153807,'rank_contrast_contexts':9,'topk_scenarios':27,
               'observed_topk_agreements':obs,'total_tested_pair_comparisons':npairs,
               'partition_block_count':{p:len(set(g)) for p,g in groups.items()},
               'observed_agreement_sign':int(rr[rr.partition=='feature'].shape[0]),
               'sign_query_radii_computed':int(len(rr)),
               'topk_witnesses_computed':len(ww)}
    (OUT/'dataset_integrity.json').write_text(json.dumps(integrity,indent=2))
    print('INTEGRITY',json.dumps(integrity),flush=True)
    print(pd.DataFrame(real_summary).to_string(index=False),flush=True)
    return df,rr,ww,integrity

if __name__=='__main__':
    start=time.perf_counter()
    v=brute_tests();print('BRUTE',v,flush=True)
    a,b,c,info=calc()
    v['elapsed_sec_total']=time.perf_counter()-start
    (OUT/'validation.json').write_text(json.dumps(v,indent=2))
    print('DONE',v,flush=True)
