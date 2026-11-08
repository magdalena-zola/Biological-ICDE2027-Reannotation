#!/usr/bin/env python3
"""Frozen scientific-output consistency tests; does not re-run either raw-data pipeline."""
import csv,json,math,pathlib
P=pathlib.Path(__file__).resolve().parent

def j(name):return json.loads((P/name).read_text())
def csvrows(name):
 with (P/name).open(newline='',encoding='utf8') as f:return list(csv.DictReader(f))
def close(x,y,atol=1e-10):assert abs(x-y)<=atol,(x,y)
primary=j('frozen/primary/validation.json')
assert primary['features']==770 and primary['samples']==34 and primary['reads']==153807
assert primary['total_contrast_queries']==2097 and primary['topk_scenarios']==27
assert primary['observed_agreement_but_uncertified']==3
assert primary['topk_observed_same']==11 and primary['topk_certified_sufficient']==10
b=j('results/dataset_integrity.json')
assert b['features']==770 and b['samples']==34 and b['reads']==153807
assert b['observed_topk_agreements']==11
r=csvrows('results/certificate_frontiers.csv')
assert len(r)==30 and all(int(x['observed_agreement'])==11 for x in r)
by={(x['partition'],str(x['budget'])):int(x['certified_topk']) for x in r}
assert by[('feature','all')]==10
assert by[('silva_family','all')]==11 and by[('silva_phylum','all')]==11
rad=csvrows('results/sign_robustness_radii.csv')
interesting={(x['taxon'],x['partition']):int(x['critical_edit_radius']) for x in rad if x['contrast']=='left palm - right palm'}
assert interesting['Aerococcaceae','feature']==1
assert interesting['Corynebacterium','feature']==1
assert interesting['Prevotella','feature']==4
rank=csvrows('results/budgeted_topk_witnesses.csv')
assert any(x['chosen']=='Pseudomonadota' and int(x['critical_edit_radius'])==44 and x['partition']=='feature' for x in rank)
correct=j('results/validation.json')
assert correct['exhaustive_cases']==2970
hier=j('results/hierarchy_validation.json')
assert hier['tree_DP_brute_cases']==1380
bench=j('results/same_operator_incremental_summary.json')
assert bench['repetitions']==3 and bench['updates_per_repetition']==120
assert bench['max_error']<1e-10
ext=j('external/results/frozen_independent_validation.json')
assert ext['feature_count']==41328 and ext['read_count']==21227600
assert ext['exhaustive_unit_test_worlds']==1320
assert ext['phylum_top1']['worst_margin_read_counts_by_budget_feature']['all']==2891234
assert ext['genus_sensitivity']['feature_level_critical_radius_all_competitors']==2
assert ext['genus_sensitivity']['pair_min_margin_read_counts_by_budget_feature']['2']==-2293
print('PASS frozen primary: 2097 queries, 27 top-k, 3 sign witnesses, 44-edit top1')
print('PASS budgeted/hierarchy: 2970 + 1380 brute cases (frozen metadata), 3 runs × 120 benchmark updates')
print('PASS external pooled ranks: 41328 ASVs / 21227600 reads, 2891234 margin, genus radius 2')
print('NOTE: this test validates frozen outputs; full raw data reconstruction is separate.')
