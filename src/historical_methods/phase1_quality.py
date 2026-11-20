from pathlib import Path
import ast,json,hashlib,collections
from concurrent.futures import ProcessPoolExecutor
import numpy as np,pandas as pd,wfdb
from wfdb import processing
from scipy.signal import butter,sosfiltfilt
ROOT=Path(__file__).resolve().parents[1];PR=ROOT/'private';OLD=Path.home()/'Desktop/POTASSIUM_ECG_RESEARCH';EX=None
LEADS=['I','II','V1','V2','V3','V4','V5','V6'];EXPECTED=set(LEADS+['III','aVR','aVL','aVF'])
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def worker(item):
 global EX
 study,paths,expected_time=item;q=dict(study_id=int(study),passes=False,flags=[])
 try:
  if len(paths)>1 and len({(sha(p+'.hea'),sha(p+'.dat')) for p in paths})>1:q['flags'].append('local_copy_disagreement')
  base=paths[0];r=wfdb.rdrecord(base);x=r.p_signal;digital=wfdb.rdrecord(base,physical=False)
  header_time=pd.Timestamp.combine(r.base_date,r.base_time) if r.base_date and r.base_time else pd.NaT
  if pd.isna(header_time) or header_time!=pd.Timestamp(expected_time):q['flags'].append('header_index_time_disagreement')
  if x.shape!=(5000,12) or r.fs!=500:q['flags'].append('shape_or_sampling')
  if set(r.sig_name)!=EXPECTED:q['flags'].append('lead_set')
  if any(u!='mV' for u in r.units) or any(g<=0 for g in r.adc_gain):q['flags'].append('unit_or_calibration')
  if not all(int(np.sum(digital.d_signal[:,i],dtype=np.int64))%65536==int(c)%65536 for i,c in enumerate(r.checksum)):q['flags'].append('checksum')
  if np.mean(~np.isfinite(x))>.01:q['flags'].append('nonfinite_gt1percent')
  for lead in LEADS:
   if lead not in r.sig_name:continue
   y=x[:,r.sig_name.index(lead)];y=y[np.isfinite(y)]
   if len(y)<100 or np.std(y)<.01:q['flags'].append('flat_or_insufficient_lead')
   if len(y) and np.percentile(np.abs(y),99.5)>10:q['flags'].append('large_amplitude')
  if not q['flags']:
   if EX is None:
    tree=ast.parse((OLD/'analyze_potassium_geometry_v2.py').read_text());nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['repair_small_gaps','geometry','axis_angle','extract']];scope=dict(np=np,wfdb=wfdb,processing=processing,butter=butter,sosfiltfilt=sosfiltfilt,ROOT=Path('/'),LEADS=LEADS);exec(compile(ast.Module(body=nodes,type_ignores=[]),'unchanged_v2_functions','exec'),scope);EX=scope['extract']
   f,a,b=EX(base);q.update(passes=True,rr_ms=f['median_rr_ms'],beats_used=f['beats_used']);dest=PR/'templates';dest.mkdir(exist_ok=True);np.savez_compressed(dest/(study+'.npz'),qrs=a,post=b)
 except Exception as ex:q['flags'].append('read_or_extraction_'+type(ex).__name__)
 return q

def stratum(h,r):
 if not np.isfinite(r) or not 2<=h<=48:return None
 return ('02-12' if h<12 else '12-24' if h<24 else '24-48','RR<=100' if abs(r)<=100 else 'RR100-200' if abs(r)<=200 else 'RR>200')
def support(p,label):
 records=[]
 for setting in ['ALL']+sorted(p.care_setting.unique().tolist()):
  d=p if setting=='ALL' else p[p.care_setting==setting];n=d[d.group=='STABLE_CONTROL']
  for group in ['HYPERKALEMIA','HYPOKALEMIA']:
   c=d[d.group==group];keys=sorted(set(c.stratum.dropna())|set(n.stratum.dropna()))
   for key in keys:
    a=int((c.stratum==key).sum());b=int((n.stratum==key).sum());records.append(dict(scope=label,care_setting=setting,group=group,stratum=key,cases=a,controls=b,shared_3perarm=a>=3 and b>=3))
 return records

def main():
 p=pd.read_pickle(PR/'selected_pairs.pkl');locations=json.loads((PR/'waveform_locations.json').read_text());wanted=set(p.study_id_a.astype(str))|set(p.study_id_b.astype(str));times={str(int(row['study_id_'+suffix])):str(row['ecg_time_'+suffix]) for row in p.to_dict('records') for suffix in ['a','b']};items=[(s,locations[s],times[s]) for s in sorted(wanted) if s in locations];qc=[]
 with ProcessPoolExecutor(max_workers=2) as pool:
  for i,q in enumerate(pool.map(worker,items),1):
   qc.append(q)
   if i%100==0:print('Available ECGs QC checked:',i,'/',len(items),flush=True)
 (PR/'waveform_qc.json').write_text(json.dumps(qc,indent=2));by={q['study_id']:q for q in qc};p['both_present']=p.study_id_a.astype(str).isin(locations)&p.study_id_b.astype(str).isin(locations);p['both_qc']=p.study_id_a.map(lambda s:by.get(s,{}).get('passes',False))&p.study_id_b.map(lambda s:by.get(s,{}).get('passes',False));p['delta_rr_ms']=p.study_id_b.map(lambda s:by.get(s,{}).get('rr_ms',np.nan))-p.study_id_a.map(lambda s:by.get(s,{}).get('rr_ms',np.nan));p['stratum']=[str(stratum(h,r)) if stratum(h,r) is not None else None for h,r in zip(p.interval_hours,p.delta_rr_ms)];p.to_pickle(PR/'quality_pairs.pkl');clean=p[p.both_qc].copy();records=support(clean,'actual_waveform_RR_QC_passing')
 meta=p.copy();a=meta.rr_interval_a.where(meta.rr_interval_a.between(250,3000));b=meta.rr_interval_b.where(meta.rr_interval_b.between(250,3000));meta['stratum']=[str(stratum(h,r)) if stratum(h,r) is not None else None for h,r in zip(meta.interval_hours,b-a)];records+=support(meta,'machine_RR_metadata_provisional');pd.DataFrame(records).to_csv(ROOT/'results/COMMON_SUPPORT.csv',index=False)
 desc=[]
 for stage,d in [('metadata',p),('both_waveform_QC',clean)]:
  for g,x in d.groupby('group'):
   desc.append(dict(stage=stage,group=g,pairs=len(x),patients=x.subject_id.nunique(),male=int((x.sex=='M').sum()),age_median=float(x.age.median()),interval_hours_median=float(x.interval_hours.median()),interval_hours_q25=float(x.interval_hours.quantile(.25)),interval_hours_q75=float(x.interval_hours.quantile(.75)),k_a_median=float(x.valuenum_a.median()),k_b_median=float(x.valuenum_b.median()),absolute_lab_gap_a_median=float(x.offset_min_a.abs().median()),absolute_lab_gap_b_median=float(x.offset_min_b.abs().median()),both_present=int(x.both_present.sum()),both_QC=int(x.both_qc.sum()),hemolysis_mention_pairs=int((x.hemolysis_mentioned_a.fillna(False)|x.hemolysis_mentioned_b.fillna(False)).sum()),RR_available=int(x.delta_rr_ms.notna().sum())))
 pd.DataFrame(desc).to_csv(ROOT/'results/CHARACTERISTICS.csv',index=False);pd.crosstab(p.group,p.care_setting).to_csv(ROOT/'results/CARE_SETTING_METADATA.csv');pd.crosstab(clean.group,clean.care_setting).to_csv(ROOT/'results/CARE_SETTING_QC.csv')
 gates=[]
 for group in ['HYPERKALEMIA','HYPOKALEMIA']:
  cells=[r for r in records if r['scope']=='actual_waveform_RR_QC_passing' and r['care_setting']=='ALL' and r['group']==group and r['shared_3perarm']];a=sum(r['cases'] for r in cells);b=sum(r['controls'] for r in cells);go=a>=100 and b>=100;gates.append(dict(group=group,common_cases=a,common_controls=b,shared_strata=len(cells),passes_gate=go,reason='Operational support gate met; not a power guarantee' if go else 'Insufficient QC-passing comparable patients: require100distinct cases and100distinct controls in shared strata with>=3perarm'))
 out=dict(waveforms_with_local_pairs=len(items),waveforms_QC_pass=sum(q['passes'] for q in qc),failure_flags=dict(collections.Counter(f for q in qc for f in set(q['flags']))),pair_QC_counts=clean.groupby('group').size().to_dict(),gates=gates,no_effect_computed=True)
 (ROOT/'results/QUALITY_AND_GATE.json').write_text(json.dumps(out,indent=2));print('Feasibility gate assessed; no effects computed.')
if __name__=='__main__':main()
