from pathlib import Path
import runpy,json
import numpy as np, pandas as pd
root=Path(__file__).resolve().parent
v=runpy.run_path(str(root/'new_validation.py'))
labels=v['labels'];wts=v['wts']; contrasts=v['contrasts'];S=v['S'];G=v['G'];features=v['features']
q=pd.read_csv(root/'Biological_ICDE_stage25/contrast_certificates.csv')
bad=q[q.observed_agreement_but_uncertified]
rows=[]
for r in bad.itertuples():
    idx=[f'{a} - {b}' for a,b in contrasts].index(r.contrast)
    w=wts[idx]
    x,y=labels[r.rank]
    ca=np.where(x==r.taxon,w,0);cb=np.where(y==r.taxon,w,0)
    pos=r.silva_effect>0 and r.gg2_effect>0
    # If both observed positive, choose minimum over independent switches.
    # If both negative, choose maximum.
    switch=(ca>cb) if pos else (cb>ca)
    chosen=np.where(switch,cb,ca)
    witness=chosen.sum()
    assert np.isclose(witness,r.feature_L if pos else r.feature_U,atol=1e-13)
    rows.append({'contrast':r.contrast,'rank':r.rank,'taxon':r.taxon,'observed_silva':r.silva_effect,'observed_gg2':r.gg2_effect,'mixed_world_witness':witness,'witness_opposite_sign':bool(witness*r.silva_effect < -1e-12),'reference_choices_gg2':int(switch.sum()),'reference_choices_silva':int((~switch).sum()),'status':'counterfactual feature-independent mixed world, NOT biological observation'})
result=pd.DataFrame(rows)
result.to_csv(root.parent/'results'/'sign_reversal_witnesses.csv',index=False)
print(result.to_string(index=False))
