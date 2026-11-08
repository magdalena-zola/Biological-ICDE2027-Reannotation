from pathlib import Path
import json,math
p=Path(__file__).resolve().parent
r=json.loads((p/'results'/'frozen_independent_validation.json').read_text())
N=r['read_count'];assert N==21227600 and r['feature_count']==41328
for k,v in r['coverage'].items():
 assert 0<=v['same_literal_reads']<=v['both_annotated_reads']<=N
 assert 0<=v['same_literal_asv']<=v['both_annotated_asv']<=41328
f=r['phylum_top1']['worst_margin_read_counts_by_budget_feature']
assert all(x>0 for x in f.values()) and f['all']==2891234 and r['phylum_top1']['n_competitors_checked']==115
g=r['genus_sensitivity']
assert g['feature_level_critical_radius_all_competitors']==2
assert g['pair_min_margin_read_counts_by_budget_feature']['1']>0
assert g['pair_min_margin_read_counts_by_budget_feature']['2']<0
assert g['pair_min_margin_read_counts_by_budget_family']['all']>0
assert sum(x['reads'] for x in g['witness']['switched_features'])==77422
print('PASS: provenance metadata, aggregate data constraints, top-1 stability, genus two-ASV witness')
