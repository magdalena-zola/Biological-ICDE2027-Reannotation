from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent/'results'
D=Path(__file__).resolve().parent/'figures'
D.mkdir(exist_ok=True)
A=pd.read_csv(R/'certificate_frontiers.csv')
# 1. certification frontier
fig,ax=plt.subplots(figsize=(7.5,3.4))
for group,label,sty in [('feature','Independent features','o-'),('silva_family','SILVA-family blocks','s--'),('silva_phylum','SILVA-phylum blocks','^:')]:
 a=A[A.partition==group].copy();a['x']=a.budget.replace({'all':770}).astype(int)
 if group=='silva_family':a.loc[a.budget=='all','x']=206
 if group=='silva_phylum':a.loc[a.budget=='all','x']=40
 a=a.sort_values('x')
 ax.plot(a.x,a.certified_topk,sty,markersize=4,label=label)
ax.set_xscale('symlog',linthresh=2)
ax.set_ylim(9.6,11.35);ax.set_yticks([10,11]);ax.set_xlabel('Maximum coherent switches from either observed reference')
ax.set_ylabel('Certified observed top-k agreements (of 11)')
ax.legend(loc='lower left',fontsize=8,ncol=1)
ax.grid(alpha=.22)
fig.tight_layout();fig.savefig(D/'fig_budgeted_frontier.png',dpi=215);fig.savefig(D/'fig_budgeted_frontier.pdf');plt.close(fig)

RADIUS=pd.read_csv(R/'sign_robustness_radii.csv');B=pd.read_csv(R/'budgeted_topk_witnesses.csv');f=RADIUS[(RADIUS.partition=='feature')&(RADIUS.critical_edit_radius>=0)].sort_values('critical_edit_radius')
rows=list(zip(f.taxon.tolist(),f.critical_edit_radius.tolist()))+[('Top-1 Pseudomonadota',44)]
fig,ax=plt.subplots(figsize=(7.0,2.9)); labs=[x[0] for x in rows]; rad=[x[1] for x in rows]
rect=ax.barh(range(len(rows)),rad);ax.set_yticks(range(len(rows)),labs);ax.invert_yaxis();ax.set_xlim(0,50)
ax.set_xlabel('Minimum number of feature switches yielding a nonpositive margin')
for b,v in zip(rect,rad):ax.text(v+.7,b.get_y()+b.get_height()/2,str(v),va='center',fontsize=9)
ax.grid(axis='x',alpha=.22);fig.tight_layout();fig.savefig(D/'fig_robustness_radii.png',dpi=215);fig.savefig(D/'fig_robustness_radii.pdf');plt.close(fig)

perf=pd.read_csv(R/'same_operator_incremental_benchmark.csv')
fig,ax=plt.subplots(figsize=(5.8,3.0)); vals=[perf.baseline_full_ms.median(),perf.indexed_ms.median()]
bar=ax.bar(['Full pair recertification','Affected pairs only'],vals)
for p,v in zip(bar,vals):ax.text(p.get_x()+p.get_width()/2,v+18,f'{v:.1f} ms',ha='center',fontsize=9)
ax.set_ylim(0,1100);ax.set_ylabel('Median total runtime (ms) over 120 updates');ax.tick_params(axis='x',labelsize=8)
ax.grid(axis='y',alpha=.2);fig.tight_layout();fig.savefig(D/'fig_same_semantics.png',dpi=215);fig.savefig(D/'fig_same_semantics.pdf');plt.close(fig)
print('figures:',list(D.glob('*.png')))
