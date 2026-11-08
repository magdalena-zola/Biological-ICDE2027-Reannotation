from pathlib import Path
import json,matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
P=Path(__file__).resolve().parent
D=json.loads((P/'frozen_independent_validation.json').read_text())
N=D['read_count']
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220,'font.family':'DejaVu Sans'})
# Figure 1: cross-reference annotation agreement on the same ASVs, measured over the same read mass.
ranks=['Phylum','Family','Genus'];key=['phylum','family','genus']
same=[100*D['coverage'][k]['same_literal_reads']/N for k in key]
diff=[100*(D['coverage'][k]['both_annotated_reads']-D['coverage'][k]['same_literal_reads'])/N for k in key]
missing=[100*(N-D['coverage'][k]['both_annotated_reads'])/N for k in key]
fig,ax=plt.subplots(figsize=(8.3,3.25));yy=range(3)
ax.barh(yy,same,label='Same literal label',color='#236F7B',height=.56)
ax.barh(yy,diff,left=same,label='Different literal label',color='#C18D4A',height=.56)
ax.barh(yy,missing,left=[a+b for a,b in zip(same,diff)],label='Missing in ≥1 reference',color='#BBC6CB',height=.56)
ax.set_yticks(list(yy),ranks);ax.invert_yaxis();ax.set_xlim(0,100);ax.set_xlabel('Share of pooled read counts (%)')
for y,v in enumerate(same):ax.text(v/2,y,f'{v:.1f}%',va='center',ha='center',color='white',fontweight='bold',fontsize=10)
ax.legend(loc='lower center',bbox_to_anchor=(.5,1.02),ncol=3,frameon=False,fontsize=9)
ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
fig.tight_layout();fig.savefig(P/'Fig_external_1_taxonomic_agreement.png',bbox_inches='tight');fig.savefig(P/'Fig_external_1_taxonomic_agreement.pdf',bbox_inches='tight');plt.close(fig)
# Figure 2: robust top-1 phylum margin over feature/family edits.
fig,ax=plt.subplots(figsize=(8.3,3.65))
budgets=['0','1','2','4','8','16','32','64','128','256','512','all']; xx=list(range(len(budgets)))
for k,col,name in [('worst_margin_read_counts_by_budget_feature','#236F7B','Independent ASV switches'),('worst_margin_read_counts_by_budget_silva_family','#B9783E','SILVA-family block switches')]:
 v=[D['phylum_top1'][k][b]/1e6 for b in budgets]
 ax.plot(xx,v,linewidth=2.1,marker='o',markersize=4.6,label=name,color=col)
ax.set_xticks(xx,budgets);ax.set_xlabel('Maximum switches from either observed reference ("all" = unrestricted)');ax.set_ylabel('Guaranteed minimum lead (million reads)');ax.set_ylim(2.84,3.025)
ax.grid(alpha=.18);ax.legend(loc='upper right',frameon=False,fontsize=9)
fig.tight_layout();fig.savefig(P/'Fig_external_2_phylum_robustness.png',bbox_inches='tight');fig.savefig(P/'Fig_external_2_phylum_robustness.pdf',bbox_inches='tight');plt.close(fig)
# Figure 3: informative genus sensitivity, real counterfactual witness.
fig,ax=plt.subplots(figsize=(8.3,3.6));budgets=['0','1','2','4','8','16','32','all'];xx=list(range(len(budgets)))
for k,col,name in [('pair_min_margin_read_counts_by_budget_feature','#236F7B','Independent ASV switches'),('pair_min_margin_read_counts_by_budget_family','#B9783E','SILVA-family block switches')]:
 v=[D['genus_sensitivity'][k][b]/1000 for b in budgets]
 ax.plot(xx,v,linewidth=2.1,marker='o',markersize=4.6,label=name,color=col)
ax.axhline(0,color='#8B2F38',linestyle='--',linewidth=1)
ax.scatter([2],[-2.293],s=95,facecolors='none',edgecolors='#8B2F38',linewidth=2.2,zorder=4)
ax.annotate('First certified counterexample (B = 2)',(2,-2.293),xytext=(2.7,21),fontsize=9,color='#7B343B',arrowprops=dict(arrowstyle='->',color='#7B343B'))
ax.set_xticks(xx,budgets);ax.set_xlabel('Maximum switches from either observed reference');ax.set_ylabel('Worst Rhodococcus − Rhodanobacter margin (thousand reads)');ax.grid(alpha=.18);ax.legend(loc='center right',bbox_to_anchor=(0.99,0.62),frameon=False,fontsize=9)
fig.tight_layout();fig.savefig(P/'Fig_external_3_genus_witness.png',bbox_inches='tight');fig.savefig(P/'Fig_external_3_genus_witness.pdf',bbox_inches='tight');plt.close(fig)
print('Figures created')
