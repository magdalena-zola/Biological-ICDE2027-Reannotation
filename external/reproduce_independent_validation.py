#!/usr/bin/env python3
"""Independent pooled-ASV ranking validation, using the gut-to-soil paired QIIME 2 export.

Input: asv-seqs-gtdb-silva-fixed.qzv, publicly linked by Greg Caporaso,
QIIME 2 Forum, 9 October 2026. The *same* ASV IDs/frequencies are classified
by SILVA 144 and GTDB r232; this is not a sample-level cohort-contrast test.

Pure Python 3.11+, no external packages required. Exact calculations for the
symmetric two-anchor block-edit worlds; independently brute-force tested.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import itertools
import json
import math
import random
import re
import sys
import zipfile
from pathlib import Path

SOURCE_URL = ('https://www.dropbox.com/scl/fi/dycygxvmvjzrn0qdgiqvs/'
              'asv-seqs-gtdb-silva-fixed.qzv?rlkey=yemtsq1s5qprs4e808878z73s&dl=1')
SOURCE_SHA256 = 'e8aa584032374072c3b8607a870fadd179bd9c6e5198d37f605c71fcddb03cfb'
RANKS = ('p','f','g')
BUDGETS = (0,1,2,4,8,16,32,64,128,256,512,1024,'all')


def tax(text: str) -> dict[str,str]:
    out={}
    for a in text.split(';'):
        a=a.strip()
        if '__' in a:
            key,value=a.split('__',1)
            if key and value and value.lower() not in ('uncultured','unclassified','unassigned','unknown','metagenome'):
                out[key]=value
    return out


def informative(label: str) -> bool:
    return bool(label) and not bool(re.search(r'incertae|uncult|unknown|unclass|unident|--other|metagenom|^other$|^norank',label,re.I))


def read_paired_qzv(path: Path):
    raw=path.read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=SOURCE_SHA256:
        raise ValueError(f'Unexpected QZV SHA-256: {digest}. Expected {SOURCE_SHA256}; refusing uncontrolled input')
    with zipfile.ZipFile(path) as z:
        candidates=[n for n in z.namelist() if n.endswith('/data/index.html')]
        assert len(candidates)==1
        html=z.read(candidates[0]).decode('utf-8')
    match=re.search(r'<script id="data" type="application/json">(.*?)</script>',html,re.S)
    if not match:
        raise ValueError('Paired-tabulation JSON not found in QIIME 2 visualization')
    parsed=json.loads(match.group(1))
    columns=[c[0] for c in parsed['columns']]
    expected=['Feature ID','Taxon: gtdb_r232','Taxon: silva_144','Frequency','No. of Samples Observed In']
    assert all(c in columns for c in expected),columns
    ix={k:columns.index(k) for k in expected}
    features=[];ids=set();total=0;sample_hits=0
    for row in parsed['data']:
        fid=row[ix['Feature ID']]
        if fid in ids:raise ValueError('Duplicate feature ID: '+fid)
        ids.add(fid)
        count=int(row[ix['Frequency']]);assert count==row[ix['Frequency']] and count>=0
        total+=count;sample_hits+=int(row[ix['No. of Samples Observed In']])
        features.append((fid,count,tax(row[ix['Taxon: silva_144']]),tax(row[ix['Taxon: gtdb_r232']])))
    assert len(features)==41328 and total==21227600,(len(features),total)
    return features,dict(features=len(features),reads=total,feature_sample_observations=sample_hits,
                         qzv_sha256=digest,qzv_bytes=len(raw))


def pair_operator(s_i,s_j,g_i,g_j):
    """Return exact worst pair margin and critical edit radius under a union of edit balls.

    The units of margins are *read counts*. For block b: d0_b=contrib_S(i)-contrib_S(j),
    d1_b=contrib_G(i)-contrib_G(j). B edits from either anchor, globally coherent per block.
    """
    d0=sum(s_i.values())-sum(s_j.values())
    d1=sum(g_i.values())-sum(g_j.values())
    diff=collections.defaultdict(int)
    for k,v in s_i.items():diff[k]-=v
    for k,v in s_j.items():diff[k]+=v
    for k,v in g_i.items():diff[k]+=v
    for k,v in g_j.items():diff[k]-=v
    neg=sorted(v for v in diff.values() if v<0)
    pos=sorted((v for v in diff.values() if v>0),reverse=True)
    a=[0];b=[0]
    for v in neg:a.append(a[-1]+v)
    for v in pos:b.append(b[-1]+v)
    def margin(k):
        k=max(len(neg),len(pos)) if k=='all' else int(k)
        return min(d0+a[min(k,len(neg))],d1-b[min(k,len(pos))])
    if margin('all')>0:
        radius=None
    else:
        lo,hi=0,max(len(neg),len(pos))
        while lo<hi:
            mid=(lo+hi)//2
            if margin(mid)<=0:hi=mid
            else:lo=mid+1
        radius=lo
    return margin,radius


def exact_brute_tests():
    rng=random.Random(20261106);n=0;max_error=0
    for blocks in range(1,9):
        for trial in range(30):
            a=[rng.uniform(-1,1) for _ in range(blocks)]
            b=[rng.uniform(-1,1) for _ in range(blocks)]
            m,_=pair_operator(dict(enumerate(a)),{},dict(enumerate(b)),{})
            for k in range(blocks+1):
                brute=min(sum(b[i] if mask>>i&1 else a[i] for i in range(blocks))
                          for mask in range(1<<blocks)
                          if min(mask.bit_count(),blocks-mask.bit_count())<=k)
                error=abs(m(k)-brute)
                max_error=max(max_error,error)
                assert error<1e-9,(blocks,trial,k,error)
                n+=1
    return n,max_error


def supports_for(features,rank,partition,strict_label_filter=False):
    support=[collections.defaultdict(lambda:collections.defaultdict(int)) for _ in range(2)]
    counts=[collections.defaultdict(int),collections.defaultdict(int)]
    blocks=set()
    for fid,c,s,g in features:
        block=fid if partition=='feature' else ('f:'+s['f'] if s.get('f') else 'asv:'+fid)
        blocks.add(block)
        for ref,anno in ((0,s),(1,g)):
            name=anno.get(rank,'')
            if name and (not strict_label_filter or informative(name)):
                support[ref][name][block]+=c
                counts[ref][name]+=c
    return support,counts,len(blocks)


def run_rank(features, rank, partition, strict_label_filter):
    sup,counts,nblocks=supports_for(features,rank,partition,strict_label_filter)
    names=set(counts[0])|set(counts[1])
    orders=[sorted(counts[r],key=lambda name:(-counts[r][name],name)) for r in (0,1)]
    output=dict(rank=rank,partition=partition,strict_label_filter=strict_label_filter,
                blocks=nblocks,labels=len(names),top5_silva=[(x,counts[0][x]) for x in orders[0][:5]],
                top5_gtdb=[(x,counts[1][x]) for x in orders[1][:5]],cases={})
    rows=[]
    for k in (1,3,5):
        matched=set(orders[0][:k])==set(orders[1][:k])
        case=dict(observed_agreement=matched,critical_radius=None,critical_pair=None,examined_pairs=0,certified={})
        if matched:
            chosen=orders[0][:k]
            minima={str(b):math.inf for b in BUDGETS}
            critical=None
            witness=None
            outsiders=names-set(chosen)
            for a in chosen:
                for challenger in outsiders:
                    margin,rho=pair_operator(sup[0].get(a,{}),sup[0].get(challenger,{}),
                                             sup[1].get(a,{}),sup[1].get(challenger,{}))
                    case['examined_pairs']+=1
                    if rho is not None and (critical is None or rho<critical):
                        critical=rho;witness=(a,challenger)
                    for b in BUDGETS:
                        minima[str(b)]=min(minima[str(b)],margin(b))
            case['critical_radius']=critical
            case['critical_pair']=witness
            case['min_margin_counts']=minima
            case['certified']={str(b):(m>0) for b,m in minima.items()}
            for budget,m in minima.items():
                rows.append(dict(rank=rank,partition=partition,strict_label_filter=int(strict_label_filter),
                                 top_k=k,budget=budget,observed_agreement=1,certified=int(m>0),
                                 worst_margin_read_counts=m,critical_radius=critical or (-1 if critical is None else 0),
                                 challenger=','.join(witness) if witness else ''))
        else:
            for budget in BUDGETS:
                rows.append(dict(rank=rank,partition=partition,strict_label_filter=int(strict_label_filter),
                                 top_k=k,budget=budget,observed_agreement=0,certified=0,
                                 worst_margin_read_counts='',critical_radius='',challenger=''))
        output['cases'][str(k)]=case
    return output,rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('qzv',type=Path,help='public QIIME 2 paired taxonomy .qzv file')
    p.add_argument('--out',type=Path,default=Path('results'))
    p.add_argument('--skip-computation',action='store_true',help='only verify arithmetic on synthetic inputs')
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    n,err=exact_brute_tests()
    summary=dict(algorithm_tests=n,max_absolute_error=err,source_url=SOURCE_URL)
    if not args.skip_computation:
        feats,meta=read_paired_qzv(args.qzv)
        summary.update(meta)
        summary['label_coverage']={}
        for rank in RANKS:
            both=equal=miss=both_read=equal_read=0
            for _,c,s,g in feats:
                a=s.get(rank,'');b=g.get(rank,'')
                if a and b:
                    both+=1;both_read+=c
                    if a==b:equal+=1;equal_read+=c
                else:miss+=1
            summary['label_coverage'][rank]=dict(both_features=both,both_reads=both_read,
                                                  literal_agreement_features=equal,
                                                  literal_agreement_reads=equal_read,
                                                  missing_at_least_one_annotation_features=miss)
        summary['evaluations']=[];rows=[]
        for rank in RANKS:
            for part in ('feature','silva_family'):
                for filter_ in ([False,True] if rank=='g' else [False]):
                    s,r=run_rank(feats,rank,part,filter_)
                    summary['evaluations'].append(s);rows.extend(r)
                    print(rank,part,'informative_only=',filter_,'top1=',s['cases']['1']['observed_agreement'],
                          'radius=',s['cases']['1']['critical_radius'],file=sys.stderr,flush=True)
        keys=list(rows[0]);csvpath=args.out/'margin_by_budget.csv'
        with csvpath.open('w',newline='',encoding='utf8') as f:
            w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    (args.out/'independent_validation.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('evaluations',)},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
