from pathlib import Path
import ast,json,hashlib,datetime,fcntl
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[1];V3=ROOT.parent;PR=ROOT/'private'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 with (PR/'.analysis.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  gate=json.loads((ROOT/'results/INDEPENDENT_RETRIEVAL_QC_GATE.json').read_text());assert gate['independently_integrity_verified']==1672 and all(g['passes'] for g in gate['gates'])
  source=V3.parent/'POTASSIUM_SERIAL_AUDIT_V1/src/analyze.py';tree=ast.parse(source.read_text());fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='effect');scope={'np':np};exec(compile(ast.Module(body=[fn],type_ignores=[]),str(source),'exec'),scope);effect=scope['effect'];(ROOT/'src/frozen_effect.py').write_text('import numpy as np\n\n'+ast.unparse(fn)+'\n')
  inputs={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'reports/PHASE2_ANALYSIS_PLAN.md',ROOT/'results/INDEPENDENT_RETRIEVAL_QC_GATE.json',PR/'pairs_verified.pkl',PR/'qc_verified.json']};inputs['frozen_effect_source_sha256']=sha(source)
  receipt=ROOT/'results/ANALYSIS_COMPLETED.json'
  if receipt.exists():assert json.loads(receipt.read_text())['inputs']==inputs;print('Identical completed fixed analysis retained; no duplicate effect execution.');return
  if (ROOT/'results/EFFECTS.json').exists():raise RuntimeError('Partial effect output exists; investigate instead of silently repeat')
  p=pd.read_pickle(PR/'pairs_verified.pkl');p=p[p.both_qc].copy();qmap={q['study_id']:q for q in json.loads((PR/'qc_verified.json').read_text())};deltas=[]
  for r in p.itertuples():
   L=[]
   for study in [r.study_id_a,r.study_id_b]:
    if r.group=='STABLE_CONTROL':q=qmap[study];v=np.log(q['original_post_rms']/q['original_qrs_rms'])
    else:
     with np.load(PR/'templates'/f'{study}.npz',allow_pickle=False) as z:v=np.log((float(np.sqrt(np.mean(z['post']**2)))+1e-12)/(float(np.sqrt(np.mean(z['qrs']**2)))+1e-12))
    L.append(float(v))
   delta=L[1]-L[0]
   if r.group!='STABLE_CONTROL' and bool(r.abnormal_first):delta=-delta
   deltas.append(delta)
  p['signed_change']=deltas;p['abs_change']=np.abs(deltas);p.to_pickle(PR/'analysis_pairs.pkl');C=p[p.group!='STABLE_CONTROL'];N=p[p.group=='STABLE_CONTROL'];assert not C.duplicated(['subject_id','group']).any() and not N.subject_id.duplicated().any() and not set(C.subject_id)&set(N.subject_id)
  rows=[];descriptive=[]
  for group in ['HYPERKALEMIA','HYPOKALEMIA']:
   c=C[C.group==group]
   masks=[('primary_30min',c,N,'Primary frozen cohort'),('both_labs_15min',c[c.offset_min_a.abs().le(15)&c.offset_min_b.abs().le(15)],N[N.offset_min_a.abs().le(15)&N.offset_min_b.abs().le(15)],'Subset only;no new matches'),('window_60min_frozen_matches',c,N,'Identical by construction;not a newly linked60min cohort'),('both_collections_at_or_before',c[c.offset_min_a.le(0)&c.offset_min_b.le(0)],N[N.offset_min_a.le(0)&N.offset_min_b.le(0)],'Mixed directions omitted'),('both_collections_at_or_after',c[c.offset_min_a.ge(0)&c.offset_min_b.ge(0)],N[N.offset_min_a.ge(0)&N.offset_min_b.ge(0)],'Mixed directions omitted;zero-gap pairs may overlap'),('exclude_hemolysis_mentions',c[~(c.hemolysis_mentioned_a.fillna(False)|c.hemolysis_mentioned_b.fillna(False))],N[~(N.hemolysis_mentioned_a.fillna(False)|N.hemolysis_mentioned_b.fillna(False))],'Exclusion in both arms;text is unadjudicated'),('exclude_source_time_disagreement',c,N,'Already required for all primary/QC pairs;identical by construction')]
   for name,cf,nf,note in masks:
    z=effect(cf,nf);ok=z['status']=='estimated';rows.append(dict(group=group,analysis=name,eligible_cases=len(cf),eligible_controls=len(nf),common_cases=z['cases'] if ok else 0,common_controls=z['controls'] if ok else 0,strata=z.get('strata',0),effect=z.get('effect'),CI95=z.get('CI95'),status=z['status'],reason='' if ok else 'No shared cell with>=3patients per arm',note=note,cells=z.get('cells',[])))
  for group,g in p.groupby('group'):
   descriptive.append(dict(group=group,patients=g.subject_id.nunique(),pairs=len(g),signed_change_median=float(g.signed_change.median()),signed_ratio_percent=100*float(np.expm1(g.signed_change.median())),absolute_change_median=float(g.abs_change.median()),age_median=float(g.age.median()),male=int((g.sex=='M').sum()),female=int((g.sex=='F').sum()),interval_median=float(g.interval_hours.median()),interval_q25=float(g.interval_hours.quantile(.25)),interval_q75=float(g.interval_hours.quantile(.75)),potassium_earlier_median=float(g.valuenum_a.median()),potassium_later_median=float(g.valuenum_b.median()),hemolysis_mention_pairs=int((g.hemolysis_mentioned_a.fillna(False)|g.hemolysis_mentioned_b.fillna(False)).sum()),both_labs_before=int((g.offset_min_a.le(0)&g.offset_min_b.le(0)).sum()),both_labs_after=int((g.offset_min_a.ge(0)&g.offset_min_b.ge(0)).sum()),both_labs_exact_zero=int((g.offset_min_a.eq(0)&g.offset_min_b.eq(0)).sum())))
  (ROOT/'results/EFFECTS.json').write_text(json.dumps(rows,indent=2));pd.DataFrame([{k:v for k,v in r.items() if k!='cells'} for r in rows]).to_csv(ROOT/'results/COMPLETE_SENSITIVITY_TABLE.csv',index=False);(ROOT/'results/DESCRIPTIVE_RESULTS.json').write_text(json.dumps(descriptive,indent=2));pd.DataFrame(descriptive).to_csv(ROOT/'results/CHARACTERISTICS.csv',index=False)
  receipt.write_text(json.dumps(dict(inputs=inputs,completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),bootstrap_replicates=2000,seed=20261119,analysis_rows=len(rows),one_fixed_execution=True,selection_or_tuning_after_effect_results=False),indent=2));print('One fixed primary/sensitivity execution completed;all14rows saved.')
if __name__=='__main__':main()
