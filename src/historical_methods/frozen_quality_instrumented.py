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
   f,a,b=EX(base);q.update(passes=True,rr_ms=f['median_rr_ms'],beats_used=f['beats_used'],original_qrs_rms=f['qrs_rms'],original_post_rms=f['post_rms']);dest=PR/'templates';dest.mkdir(exist_ok=True);np.savez_compressed(dest/(study+'.npz'),qrs=a,post=b)
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

