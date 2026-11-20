from pathlib import Path
import json,hashlib,shutil,datetime
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[1];V3=ROOT.parent/'POTASSIUM_SERIAL_V3';SOURCE=V3/'phase2_post_retrieval_v1';PR=ROOT/'private'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def describe(x):
 x=np.asarray(x,dtype=float);return dict(n=len(x),mean=float(x.mean()),SD=float(x.std(ddof=1)),median=float(np.median(x)),q25=float(np.quantile(x,.25)),q75=float(np.quantile(x,.75)),minimum=float(x.min()),maximum=float(x.max()))
def main():
 preserve=PR/'PRESERVATION.json'
 if not preserve.exists():preserve.write_text(json.dumps({str(p.relative_to(V3)):sha(p) for p in V3.rglob('*') if p.is_file() and not p.is_symlink()},indent=2))
 p=pd.read_pickle(SOURCE/'private/analysis_pairs.pkl');qc={r['study_id']:r for r in json.loads((SOURCE/'private/qc_verified.json').read_text())};effects=json.loads((SOURCE/'results/EFFECTS.json').read_text());assert len(effects)==14
 counts=p.groupby('group').size().to_dict();assert counts=={'HYPERKALEMIA':122,'HYPOKALEMIA':200,'STABLE_CONTROL':650};assert not p.duplicated(['subject_id','group']).any();assert not set(p.loc[p.group=='STABLE_CONTROL','subject_id'])&set(p.loc[p.group!='STABLE_CONTROL','subject_id'])
 K=[];D=[];errors=[]
 for group,g in p.groupby('group'):
  if group=='STABLE_CONTROL':
   assert g.valuenum_a.between(3.5,5).all() and g.valuenum_b.between(3.5,5).all();roles=[('earlier_normal',g.valuenum_a),('later_normal',g.valuenum_b)]
  else:
   first=g.abnormal_first.astype(bool);a=np.where(first,g.valuenum_a,g.valuenum_b);n=np.where(first,g.valuenum_b,g.valuenum_a);assert np.all((n>=3.5)&(n<=5));assert np.all(a>=5.5) if group=='HYPERKALEMIA' else np.all(a<3.5);roles=[('abnormal',a),('paired_normal',n)]
  for role,x in roles:K.append(dict(group=group,role=role,**describe(x)))
 for r in p.itertuples():
  native=[];reduced64=[];original64=[]
  for s in [r.study_id_a,r.study_id_b]:
   with np.load(SOURCE/'private/templates'/f'{s}.npz',allow_pickle=False) as z:
    q32=float(np.sqrt(np.mean(z['qrs']**2)));p32=float(np.sqrt(np.mean(z['post']**2)));q64=float(np.sqrt(np.mean(z['qrs'].astype('float64')**2)));p64=float(np.sqrt(np.mean(z['post'].astype('float64')**2)))
   native.append(float(np.log((p32+1e-12)/(q32+1e-12))));reduced64.append(float(np.log(p64/q64)));original64.append(float(np.log(qc[s]['original_post_rms']/qc[s]['original_qrs_rms'])))
  sign=-1 if r.group!='STABLE_CONTROL' and bool(r.abnormal_first) else 1
  v32=sign*(native[1]-native[0]);v64=sign*(original64[1]-original64[0]);vr=sign*(reduced64[1]-reduced64[0]);legacy=v64 if r.group=='STABLE_CONTROL' else v32
  errors.append(abs(legacy-r.signed_change));assert np.isfinite([v32,v64,vr,legacy]).all()
  D.append(dict(group=r.group,subject_id=r.subject_id,signed_change_legacy=legacy,absolute_change_legacy=abs(legacy),uniform_original64_signed=v64,pair_native32_vs_original64=abs(v32-v64),pair_native32_vs_reduction64=abs(v32-vr),absolute_outcome_discrepancy_to_uniform64=abs(abs(legacy)-abs(v64)),max_ECG_logscore_native32_vs_original64=float(np.max(np.abs(np.array(native)-np.array(original64))))))
 d=pd.DataFrame(D);d.to_pickle(PR/'precision_per_pair.pkl');summary=[]
 for group,g in d.groupby('group'):
  z=dict(group=group,pairs=len(g))
  for col in ['pair_native32_vs_original64','pair_native32_vs_reduction64','absolute_outcome_discrepancy_to_uniform64','max_ECG_logscore_native32_vs_original64']:z[col+'_max']=float(g[col].max());z[col+'_median']=float(g[col].median())
  summary.append(z)
 by={r['group']:r for r in summary};bounds=[]
 for group in ['HYPERKALEMIA','HYPOKALEMIA']:
  bound=by[group]['absolute_outcome_discrepancy_to_uniform64_max']+by['STABLE_CONTROL']['absolute_outcome_discrepancy_to_uniform64_max'];r=next(r for r in effects if r['group']==group and r['analysis']=='primary_30min');values=[r['effect']]+r['CI95'];margins=[min((abs(v)*10000)%1,1-((abs(v)*10000)%1))/10000 for v in values];zero_margin=min(abs(v) for v in r['CI95']);bounds.append(dict(group=group,max_possible_fixed_weight_contrast_shift=bound,nearest_four_decimal_rounding_boundary=min(abs((abs(v)*10000)%1-.5)/10000 for v in values),CI_zero_margin=zero_margin))
 material=max(errors)>1e-12 or any(r['max_possible_fixed_weight_contrast_shift']>1e-5 or r['max_possible_fixed_weight_contrast_shift']>=r['nearest_four_decimal_rounding_boundary'] or r['max_possible_fixed_weight_contrast_shift']>=r['CI_zero_margin'] for r in bounds)
 out=dict(fixed_counts=counts,unique_patients=int(p.subject_id.nunique()),case_group_overlap=len(set(p.loc[p.group=='HYPERKALEMIA','subject_id'])&set(p.loc[p.group=='HYPOKALEMIA','subject_id'])),saved_outcome_reconstruction_max_error=max(errors),numerical_materiality_threshold_logunits=1e-5,contrast_shift_bounds=bounds,group_discrepancies=summary,material_problem=material,no_effect_function_or_bootstrap_called=True,results_modified=False,interpretation='Individual numerical-representation check,not a new effect analysis',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
 (ROOT/'results/FIXED_PAIR_PRECISION_CHECK.json').write_text(json.dumps(out,indent=2));pd.DataFrame(K).to_csv(ROOT/'results/POTASSIUM_DISTRIBUTIONS.csv',index=False);pd.DataFrame(summary).to_csv(ROOT/'results/PRECISION_DISCREPANCIES.csv',index=False)
 if material:raise RuntimeError('Material precision/consistency problem: stop manuscript finishing and report')
 for n in ['EFFECTS.json','COMPLETE_SENSITIVITY_TABLE.csv','ANALYSIS_COMPLETED.json','INDEPENDENT_RETRIEVAL_QC_GATE.json','QC_SUMMARY.json','POST_RETRIEVAL_ATTRITION.csv','SOURCE_COHORT_FLOW.csv','ECG_LAB_TIMING.csv','CARE_SETTING.csv','DESCRIPTIVE_RESULTS.json','VERIFIED_STRATA.csv','SOURCE_INVENTORY.json']:shutil.copy2(SOURCE/'results'/n,ROOT/'results'/n)
 (ROOT/'results/SAVED_RESULTS_LOCK.json').write_text(json.dumps({n:sha(SOURCE/'results'/n) for n in ['EFFECTS.json','COMPLETE_SENSITIVITY_TABLE.csv','ANALYSIS_COMPLETED.json']},indent=2));print('Fixed counts and numerical guard passed;all14saved rows copied unchanged. No effect analysis rerun.')
if __name__=='__main__':main()
