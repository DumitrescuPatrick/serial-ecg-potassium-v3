from pathlib import Path
import json,hashlib,os,collections,time
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];PR=ROOT/'private';HOME=Path.home()
from inventory import SOURCES
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def save(n,x):(ROOT/'results'/n).write_text(json.dumps(x,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)))
def main():
 assert (ROOT/'reports/PHASE1_PROTOCOL.md').exists();flow=[]
 def count(stage,d):flow.append(dict(stage=stage,ecgs=len(d),patients=d.subject_id.nunique()))
 idx=pd.read_csv(SOURCES['index']);idx['ecg_time']=pd.to_datetime(idx.ecg_time,errors='coerce');count('full_provider_index',idx);assert not idx.study_id.duplicated().any()
 m=pd.read_csv(SOURCES['machine'],usecols=['subject_id','study_id','ecg_time','rr_interval']);m['ecg_time']=pd.to_datetime(m.ecg_time,errors='coerce');assert not m.study_id.duplicated().any();e=idx.merge(m.rename(columns={'ecg_time':'machine_time'}),on=['subject_id','study_id'],validate='one_to_one');e=e[e.ecg_time.notna()&(e.ecg_time==e.machine_time)].copy();count('index_machine_time_agreement',e)
 p=pd.read_csv(SOURCES['patients'],usecols=['subject_id','gender','anchor_age','anchor_year']);e=e.merge(p,on='subject_id',validate='many_to_one');e['age']=e.anchor_age+e.ecg_time.dt.year-e.anchor_year;e=e[e.age>=18];count('adult',e)
 a=pd.read_csv(SOURCES['admissions'],usecols=['subject_id','hadm_id','admittime','dischtime','edregtime']);assert not a.hadm_id.duplicated().any()
 for c in ['admittime','dischtime','edregtime']:a[c]=pd.to_datetime(a[c],errors='coerce')
 a['start']=a[['admittime','edregtime']].min(axis=1);a=a[a.start.notna()&a.dischtime.notna()&(a.start<=a.dischtime)]
 j=e[['subject_id','study_id','ecg_time']].merge(a,on='subject_id');j=j[(j.ecg_time>=j.start)&(j.ecg_time<=j.dischtime)];n=j.groupby('study_id').size();j=j[j.study_id.isin(n[n==1].index)];e=e.merge(j[['study_id','hadm_id','admittime','start','dischtime']],on='study_id',validate='one_to_one');count('unique_corroborated_hospital_encounter',e)
 ed=pd.read_csv(SOURCES['ed'],usecols=['subject_id','hadm_id','stay_id','intime','outtime']);ed['intime']=pd.to_datetime(ed.intime,errors='coerce');ed['outtime']=pd.to_datetime(ed.outtime,errors='coerce');ed=ed[ed.hadm_id.notna()&(ed.outtime>=ed.intime)].copy();ed.hadm_id=ed.hadm_id.astype('int64')
 ej=e[['subject_id','hadm_id','study_id','ecg_time']].merge(ed,on=['subject_id','hadm_id']);ej=ej[(ej.ecg_time>=ej.intime)&(ej.ecg_time<=ej.outtime)];ed_ids=set(ej.study_id);e['care_phase']=np.where(e.study_id.isin(ed_ids),'corroborated_ED','inpatient');e=e[(e.ecg_time>=e.admittime)|e.study_id.isin(ed_ids)].copy();count('pre_admission_ECGs_require_matching_ED_stay',e)
 dictionary=pd.read_csv(SOURCES['dictionary']);concept=dictionary[dictionary.itemid==50971];assert len(concept)==1 and concept.iloc[0].label=='Potassium' and concept.iloc[0].fluid=='Blood' and concept.iloc[0].category=='Chemistry'
 cache=PR/'raw_item50971.pkl'
 if cache.exists():raise RuntimeError('New study already has a lab cache; inspect saved state rather than overwrite')
 chunks=[];rows=0;cols=['subject_id','labevent_id','hadm_id','specimen_id','itemid','charttime','storetime','value','valuenum','valueuom','comments']
 for ch in pd.read_csv(SOURCES['labs'],usecols=cols,chunksize=500000,low_memory=False):
  rows+=len(ch);chunks.append(ch[ch.itemid==50971])
  if rows%25000000==0:print('Raw rows scanned:',rows,flush=True)
 lab=pd.concat(chunks,ignore_index=True);lab['charttime']=pd.to_datetime(lab.charttime,errors='coerce');lab['storetime']=pd.to_datetime(lab.storetime,errors='coerce');lab.to_pickle(cache)
 initial=len(lab);unit=lab.valueuom.fillna('').str.lower().str.replace(' ','',regex=False);lab['hemolysis_mentioned']=lab.comments.fillna('').str.contains('hemoly|haemoly',case=False,regex=True);invalidity=lab.comments.fillna('').str.contains('error|cancel|invalid|unable|disregard',case=False,regex=True)
 lab['unsupported_unit']=~unit.isin(['meq/l','mmol/l']);lab['missing_lab_admission']=lab.hadm_id.isna();lab['missing_specimen']=lab.specimen_id.isna()
 valid=lab.valuenum.between(1.5,10)&~lab.unsupported_unit&~lab.value.fillna('').str.contains('[<>]',regex=True)&~invalidity&lab.charttime.notna()&lab.storetime.notna()&(lab.storetime>=lab.charttime)
 save('LAB_FLAGS.json',dict(item50971_rows=initial,hemolysis_mentions=int(lab.hemolysis_mentioned.sum()),unsupported_or_missing_units=int(lab.unsupported_unit.sum()),missing_admission_ids=int(lab.missing_lab_admission.sum()),missing_specimen_ids=int(lab.missing_specimen.sum()),valid_numeric_time_unit_rows=int(valid.sum()),raw_rows_scanned=rows,hemolysis_mentions_retained_primary=True))
 lab=lab[valid].copy();lab['same_time_count']=lab.groupby(['subject_id','charttime']).valuenum.transform('nunique');lab=lab.sort_values(['charttime','labevent_id']).drop_duplicates(['subject_id','charttime']);lab=lab.rename(columns={'hadm_id':'lab_hadm_id'});lab=lab.sort_values('charttime');e=e.sort_values('ecg_time')
 lc=['subject_id','charttime','labevent_id','specimen_id','valuenum','same_time_count','lab_hadm_id','hemolysis_mentioned'];back=pd.merge_asof(e,lab[lc],left_on='ecg_time',right_on='charttime',by='subject_id',direction='backward',tolerance=pd.Timedelta('30min'));forward=pd.merge_asof(e,lab[lc],left_on='ecg_time',right_on='charttime',by='subject_id',direction='forward',tolerance=pd.Timedelta('30min'))
 db=(back.ecg_time-back.charttime).dt.total_seconds().abs();df=(forward.charttime-forward.ecg_time).dt.total_seconds().abs();choosef=df.notna()&(db.isna()|(df<db));chosen=back.copy()
 for c in lc[1:]:chosen.loc[choosef,c]=forward.loc[choosef,c]
 chosen['tie_conflict']=db.notna()&df.notna()&(db==df)&(back.valuenum!=forward.valuenum);chosen=chosen[chosen.charttime.notna()].copy();count('nearest_valid_lab_within30min_before_admission_filter',chosen)
 chosen.to_pickle(PR/'nearest_lab_candidates.pkl');save('MATCHED_LAB_FLAGS.json',dict(nearest_candidates=len(chosen),missing_lab_admission=int(chosen.lab_hadm_id.isna().sum()),lab_admission_disagrees=int((chosen.lab_hadm_id.notna()&(chosen.lab_hadm_id!=chosen.hadm_id)).sum()),ambiguous_same_time=int((chosen.same_time_count>1).sum()),ambiguous_equal_distance=int(chosen.tie_conflict.sum()),hemolysis_mentions=int(chosen.hemolysis_mentioned.fillna(False).sum())))
 chosen=chosen[chosen.lab_hadm_id.notna()&(chosen.lab_hadm_id==chosen.hadm_id)&(chosen.charttime>=chosen.start)&(chosen.charttime<=chosen.dischtime)].copy();count('recorded_lab_admission_equal_to_corroborated_ECG_encounter',chosen)
 chosen=chosen[~chosen.tie_conflict&(chosen.same_time_count==1)&chosen.specimen_id.notna()].copy();count('unambiguous_nonmissing_specimen',chosen);chosen['offset_min']=(chosen.charttime-chosen.ecg_time).dt.total_seconds()/60;chosen.to_pickle(PR/'eligible_ecgs.pkl')
 best={};candidate_pairs=collections.Counter()
 for (sid,hadm),g in chosen.groupby(['subject_id','hadm_id'],sort=False):
  rr=g.sort_values(['ecg_time','study_id']).to_dict('records')
  for i,x in enumerate(rr):
   for y in rr[i+1:]:
    h=(y['ecg_time']-x['ecg_time']).total_seconds()/3600
    if h>48:break
    if h<2:continue
    group=None;abnormal_first=None
    nx=3.5<=x['valuenum']<=5;ny=3.5<=y['valuenum']<=5
    if nx and ny and abs(x['valuenum']-y['valuenum'])<=.30+1e-9:group='STABLE_CONTROL'
    elif (x['valuenum']>=5.5 and ny) or (y['valuenum']>=5.5 and nx):group='HYPERKALEMIA';abnormal_first=x['valuenum']>=5.5
    elif (x['valuenum']<3.5 and ny) or (y['valuenum']<3.5 and nx):group='HYPOKALEMIA';abnormal_first=x['valuenum']<3.5
    if group is None:continue
    candidate_pairs[group]+=1;key=(sid,group);rank=(h,x['ecg_time'],x['study_id'],y['study_id']);care=x['care_phase'] if x['care_phase']==y['care_phase'] else 'mixed_phase'
    out=dict(subject_id=sid,group=group,abnormal_first=abnormal_first,hadm_id=hadm,interval_hours=h,care_setting=care,age=x['age'],sex=x['gender'])
    for suffix,r in [('a',x),('b',y)]:
     for c in ['study_id','path','ecg_time','valuenum','labevent_id','charttime','offset_min','rr_interval','hemolysis_mentioned','lab_hadm_id']:out[c+'_'+suffix]=r[c]
    if key not in best or rank<best[key][0]:best[key]=(rank,out)
 pairs=pd.DataFrame([v[1] for v in best.values()]);caseids=set(pairs.loc[pairs.group!='STABLE_CONTROL','subject_id']);excluded=int(((pairs.group=='STABLE_CONTROL')&pairs.subject_id.isin(caseids)).sum());pairs=pairs[~((pairs.group=='STABLE_CONTROL')&pairs.subject_id.isin(caseids))].copy();pairs.to_pickle(PR/'selected_pairs.pkl')
 available=json.loads((PR/'waveform_locations.json').read_text());summary=[]
 for group,g in pairs.groupby('group'):
  aa=g.study_id_a.astype(str).isin(available);bb=g.study_id_b.astype(str).isin(available);summary.append(dict(group=group,metadata_pairs=len(g),patients=g.subject_id.nunique(),both_waveforms_present=int((aa&bb).sum()),one_waveform_present=int((aa^bb).sum()),neither_present=int((~aa&~bb).sum())))
 save('COHORT_FLOW.json',flow);save('WAVEFORM_AVAILABILITY.json',summary);save('PAIR_SELECTION.json',dict(candidate_pairs=dict(candidate_pairs),selected_pairs=pairs.groupby('group').size().to_dict(),controls_excluded_for_case_patient_overlap=excluded,hyper_hypo_patient_overlap=len(set(pairs.loc[pairs.group=='HYPERKALEMIA','subject_id'])&set(pairs.loc[pairs.group=='HYPOKALEMIA','subject_id'])),lab_item50971_rows=initial,valid_lab_time_representatives_before_admission_filter=len(lab),raw_rows_scanned=rows,no_icd_used=True,window_minutes=30));print('Metadata cohort and availability complete; no effects calculated.')
if __name__=='__main__':main()
