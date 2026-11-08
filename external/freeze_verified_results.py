from pathlib import Path
import json,csv,hashlib
P=Path(__file__).resolve().parent
counts=21227600
result={
 'study':'Gut-to-soil QIIME 2 dual-classifier sequence tabulation',
 'source_publication':'QIIME 2 Forum, Greg Caporaso, 2026-10-09',
 'source_url':'https://forum.qiime2.org/t/using-multiple-taxonomy-classifiers-for-comparison-or-in-filtering-workflows/34349',
 'original_experiment_doi':'https://doi.org/10.5281/zenodo.15390940',
 'reference_pair':['SILVA 144','GTDB r232'],
 'data_regime':'Single pooled feature-frequency vector; not per-sample group contrasts',
 'input_file':'asv-seqs-gtdb-silva-fixed.qzv',
 'input_sha256':'e8aa584032374072c3b8607a870fadd179bd9c6e5198d37f605c71fcddb03cfb',
 'input_file_bytes':21065310,
 'feature_count':41328,
 'read_count':counts,
 'feature_sample_observations_sum':314497,
 'max_samples_observed_in_per_feature':825,
 'validation_date':'2026-11-06',
 'execution_location':'Public-source network workspace; summary exported to local package; original binary is not redistributed',
 'exhaustive_unit_test_worlds':1320,
 'exhaustive_unit_test_max_abs_error':8.881784197001252e-16,
 'coverage':{
 'phylum':{'both_annotated_asv':39214,'both_annotated_reads':20630628,'same_literal_asv':31250,'same_literal_reads':17377393},
 'family':{'both_annotated_asv':32754,'both_annotated_reads':18483547,'same_literal_asv':16042,'same_literal_reads':9528607},
 'genus':{'both_annotated_asv':21286,'both_annotated_reads':11283622,'same_literal_asv':6243,'same_literal_reads':5669841}},
 'observed_rankings':{
 'phylum':{'silva_top5':['Pseudomonadota','Actinomycetota','Bacteroidota','Clostridiota','Bacillota'],'gtdb_top5':['Pseudomonadota','Actinomycetota','Bacillota','Bacteroidota','Verrucomicrobiota'],'same_top1':True,'same_top3':False,'same_top5':False},
 'family':{'silva_top5':['Alcaligenaceae','Flavobacteriaceae','Sphingobacteriaceae','Rhizobiaceae--other','Nocardiaceae--other'],'gtdb_top5':['Mycobacteriaceae','Chitinophagaceae','Burkholderiaceae_B','Rhizobiaceae','Sphingobacteriaceae'],'same_top1':False,'same_top3':False,'same_top5':False},
 'genus_unfiltered':{'silva_top5':['Incertae_Sedis','Rhodococcus','Acinetobacter','Devosia','Pseudomonas'],'gtdb_top5':['Rhodococcus','Romboutsia','Rhodanobacter','Caryophanon','Devosia'],'same_top1':False,'same_top3':False,'same_top5':False},
 'genus_informative_only':{'silva_top5':['Rhodococcus','Acinetobacter','Devosia','Pseudomonas','Albibacterium'],'gtdb_top5':['Rhodococcus','Romboutsia','Rhodanobacter','Caryophanon','Devosia'],'same_top1':True,'same_top3':False,'same_top5':False}},
 'phylum_top1':{
 'taxon':'Pseudomonadota','silva_reads':7806350,'gtdb_reads':7041539,
 'certified_for_all_feature_level_binary_worlds':True,'critical_radius':None,
 'worst_margin_read_counts_by_budget_feature':{'0':2984494,'1':2955664,'2':2945571,'4':2931474,'8':2922520,'16':2911954,'32':2903762,'64':2897820,'128':2893816,'256':2891794,'512':2891234,'1024':2891234,'all':2891234},
 'worst_margin_read_counts_by_budget_silva_family':{'0':2984494,'1':2926204,'2':2912690,'4':2905773,'8':2900430,'16':2896311,'32':2894192,'64':2893482,'128':2893257,'256':2893195,'512':2893195,'1024':2893195,'all':2893195},
 'n_competitors_checked':115},
 'genus_sensitivity':{
 'filter':'Exclude labels matching incertae/uncult/unknown/unclass/unident/--other/metagenom/^other$/^norank',
 'taxon':'Rhodococcus', 'silva_reads':651711,'gtdb_reads':523231,
 'best_challenger':'Rhodanobacter',
 'feature_level_critical_radius_all_competitors':2,
 'family_grouped_critical_radius_all_competitors':None,
 'pair_min_margin_read_counts_by_budget_feature':{'0':75129,'1':3874,'2':-2293,'4':-3362,'8':-4152,'16':-4353,'32':-4366,'64':-4366,'128':-4366,'all':-4366},
 'pair_min_margin_read_counts_by_budget_family':{'0':75129,'1':75118,'2':75100,'4':75096,'8':75096,'16':75096,'32':75096,'64':75096,'128':75096,'all':75096},
 'witness':{'anchor':'GTDB r232','switched_features':[
 {'id':'f80c0a86f49e5740576bda5302b82d72','reads':71255,'silva_genus':'Rhodanobacter','gtdb_genus':'unassigned'},
 {'id':'a5a401d92657b7ac6e4d96863ad2e1d7','reads':6167,'silva_genus':'Rhodanobacter','gtdb_genus':'unassigned'}],
 'margin_after_two_switches':-2293}},
 'limitations':['Reference versions differ from the 2024.09/138 pair in the primary paper',
 'Pooled counts cannot evaluate phenotype or site contrasts or subject-level biological uncertainty',
 'Literal matching does not validate synonym/lineage-level biological equivalence',
 'Independent ASV or SILVA-family switches are scenario constraints, not estimated evolutionary dynamics',
 'No raw data or original sequences are redistributed in this derivative release']
}
(P/'frozen_independent_validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
with (P/'taxonomic_agreement.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=['rank','asv_count','both_annotated_asv','both_annotated_reads','same_literal_asv','same_literal_reads','same_pct_total_reads','different_pct_total_reads','missing_pct_total_reads']);w.writeheader()
 for rank,v in result['coverage'].items():
  w.writerow(dict(rank=rank,asv_count=41328,**v,same_pct_total_reads=100*v['same_literal_reads']/counts,different_pct_total_reads=100*(v['both_annotated_reads']-v['same_literal_reads'])/counts,missing_pct_total_reads=100*(counts-v['both_annotated_reads'])/counts))
with (P/'certified_margin_frontiers.csv').open('w',newline='',encoding='utf8') as f:
 w=csv.DictWriter(f,fieldnames=['query','partition','budget','worst_margin_read_counts','worst_margin_percent_of_total_reads']);w.writeheader()
 for query,parties in [('phylum_Pseudomonadota_top1',('phylum_top1','worst_margin_read_counts_by_budget_feature','worst_margin_read_counts_by_budget_silva_family')),
                        ('genus_Rhodococcus_minus_Rhodanobacter',('genus_sensitivity','pair_min_margin_read_counts_by_budget_feature','pair_min_margin_read_counts_by_budget_family'))]:
  rec=result[parties[0]]
  for p,k in [('feature',parties[1]),('silva_family',parties[2])]:
   for b,v in rec[k].items():w.writerow(dict(query=query,partition=p,budget=b,worst_margin_read_counts=v,worst_margin_percent_of_total_reads=100*v/counts))
print('Wrote frozen results, curves and agreement CSV')
