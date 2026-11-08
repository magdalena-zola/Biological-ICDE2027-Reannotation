from pathlib import Path
import pandas as pd, numpy as np, json, matplotlib.pyplot as plt, zipfile
import argparse
ap=argparse.ArgumentParser()
ap.add_argument('--input', default='feature_pairs_audit.csv')
ap.add_argument('--out', default='stage18_results')
args=ap.parse_args()
src=Path(args.input)
out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(src).fillna('[missing]'); n=len(d)
levels=['phylum','family','genus']
for l in levels:
    for p in ('silva','gg2'): d[f'{p}_{l}']=d[f'{p}_{l}'].astype(str).str.strip()
# Conservative synonym at phylum, no unverified equivalence.
d['silva_phylum_norm']=d.silva_phylum.replace({'Proteobacteria':'Pseudomonadota'})
d['gg2_phylum_norm']=d.gg2_phylum.replace({'Proteobacteria':'Pseudomonadota'})
rows=[]
for level in levels:
    a=d[f'silva_{level}_norm'] if level=='phylum' else d[f'silva_{level}']
    b=d[f'gg2_{level}_norm'] if level=='phylum' else d[f'gg2_{level}']
    eq=a.eq(b) & a.ne('[unassigned]') & b.ne('[unassigned]'); rows.append({'rank':level,'n_features':n,'exact_match_features':int(eq.sum()),'exact_match_read_fraction':float(d.loc[eq,'reads'].sum()/d.reads.sum()),'nonmatching_reads':float(d.loc[~eq,'reads'].sum())})
pd.DataFrame(rows).to_csv(out/'literal_agreement_audit.csv',index=False)
# Non-bijective mapping audit (same SILVA label maps to many GG2 labels and vice versa).
maprows=[]
for level in levels:
    for left,right in [('silva','gg2'),('gg2','silva')]:
        x=d.groupby(f'{left}_{level}').agg(num_targets=(f'{right}_{level}','nunique'),features=('feature_id','size'),reads=('reads','sum')).reset_index()
        for _,r in x.iterrows(): maprows.append({'rank':level,'source':left,'source_label':r[f'{left}_{level}'],'targets':int(r.num_targets),'features':int(r.features),'reads':float(r.reads)})
mapdf=pd.DataFrame(maprows);mapdf.to_csv(out/'mapping_ambiguity.csv',index=False)
# Cross-reference assignment conflicts for SILVA families in same phylum, and SILVA genera in same family.
groups=[]
for level in ['phylum','family']:
    k=f'silva_{level}'
    for label,part in d.groupby(k):
        if len(part)<2: continue
        out_rank='family' if level=='phylum' else 'genus'
        changes=part[f'silva_{out_rank}'].ne(part[f'gg2_{out_rank}'])
        groups.append({'group_level':level,'silva_group':label,'features':len(part),'reads':float(part.reads.sum()),'changed_fraction':float(changes.mean()),'changed_read_fraction':float(part.loc[changes,'reads'].sum()/part.reads.sum()),'gg2_distinct_target_'+out_rank:int(part[f'gg2_{out_rank}'].nunique())})
pd.DataFrame(groups).to_csv(out/'observed_group_heterogeneity.csv',index=False)
# observed mixed-change groups demonstrate joint all-or-nothing choice not compelled by data
summary=[]
for level in ['phylum','family']:
    x=pd.DataFrame(groups);x=x[x.group_level==level];mixed=x[(x.changed_fraction>0)&(x.changed_fraction<1)]
    summary.append({'level':level,'groups_2plus_features':len(x),'mixed_changed_unchanged_groups':len(mixed),'mixed_group_fraction':len(mixed)/len(x) if len(x) else None,'reads_in_mixed_groups':float(mixed.reads.sum())})
pd.DataFrame(summary).to_csv(out/'group_heterogeneity_summary.csv',index=False)
# actual difference concentration across groups; counts different annotations not independence observations
family=d.assign(changed=d.silva_genus.ne(d.gg2_genus)).groupby('silva_family').apply(lambda g:pd.Series({'features':len(g),'reads':g.reads.sum(),'changed_reads':g.loc[g.changed,'reads'].sum()}),include_groups=False).reset_index()
family=family.sort_values('changed_reads',ascending=False);family.to_csv(out/'family_disagreement_concentration.csv',index=False)
fig,ax=plt.subplots(figsize=(8.2,4.5));a=pd.DataFrame(rows);ax.bar(a['rank'],a['exact_match_read_fraction']*100);ax.set_ylabel('Read-weighted identical labels (%)');ax.set_ylim(0,100);ax.set_xlabel('Taxonomic rank');fig.tight_layout();fig.savefig(out/'fig18_read_agreement.png',dpi=210);plt.close(fig)
fig,ax=plt.subplots(figsize=(10,5));top=family.head(12).sort_values('changed_reads');ax.barh(top.silva_family.astype(str),top.changed_reads);ax.set_xlabel('Reads assigned to different genus labels');fig.tight_layout();fig.savefig(out/'fig18_disagreement_concentration.png',dpi=210);plt.close(fig)
readme=f'''# Stage 18: Empirical taxonomy-constraint audit\n\nInput: stage12 paired SILVA/Greengenes2 mapping of {n} features; total reads {int(d.reads.sum())}.\n\nA paired pair of annotations gives two observed complete assignments. It does NOT identify a stochastic or hard constraint relationship between local assignment variables. SILVA phylum/family trees are a taxonomy hierarchy, not an experimentally verified dependency tree of reannotation switches.\n\nWe audit literal match (one explicitly declared nomenclature rule only), one-to-many labels and heterogeneity of switching within SILVA groups. Avoid asserting that mixed changed/unchanged labels conclusively falsify group switching: identical labels can persist under either global classifier.\n\nResults files: literal_agreement_audit.csv; mapping_ambiguity.csv; observed_group_heterogeneity.csv; group_heterogeneity_summary.csv; family_disagreement_concentration.csv; two figures.\n\nRecommended next evidence: external taxonomy crosswalk, multiple database releases or perturbation replicates, and a complete native incremental-baseline implementation. Stage 17 hierarchy results are synthetic and should not be presented as biologically validated.\n'''
(out/'README.md').write_text(readme)
print('rows',n,'reads',d.reads.sum());print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(summary).to_string(index=False))
