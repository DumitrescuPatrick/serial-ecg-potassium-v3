from pathlib import Path
import json,hashlib,os,time,csv,gzip
ROOT=Path(__file__).resolve().parents[1];HOME=Path.home()
SOURCES=dict(index=HOME/'MIMIC-DATABASES/ECG/record_list.csv',machine=HOME/'Downloads/machine_measurements.csv',labs=HOME/'Downloads/labevents.csv.gz',dictionary=HOME/'Downloads/d_labitems.csv.gz',patients=HOME/'Downloads/patients.csv.gz',admissions=HOME/'Downloads/admissions.csv.gz',ed=HOME/'Downloads/edstays.csv.gz',extractor=HOME/'Desktop/POTASSIUM_ECG_RESEARCH/analyze_potassium_geometry_v2.py')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def main():
 rows=[]
 for k,p in SOURCES.items():
  r=dict(role=k,path=str(p),exists=p.is_file(),readable=os.access(p,os.R_OK))
  if p.is_file():
   r.update(bytes=p.stat().st_size,sha256=sha(p))
   if p.suffix in ['.csv','.gz']:
    with (gzip.open(p,'rt',newline='') if p.suffix=='.gz' else p.open(newline='')) as f:cols=next(csv.reader(f));assert len(cols)==len(set(cols));r['columns']=cols
  rows.append(r)
 (ROOT/'private/SOURCE_INVENTORY_PRIVATE.json').write_text(json.dumps(rows,indent=2));public=[{k:v for k,v in r.items() if k!='path'} for r in rows];(ROOT/'results/SOURCE_INVENTORY.json').write_text(json.dumps(public,indent=2))
 roots=[HOME/'Desktop/POTASSIUM_ECG_RESEARCH/waveforms_mimic_iv_ecg_v1',HOME/'MIMIC-DATABASES/ECG/matched_echo_waveforms',HOME/'MIMIC-DATABASES/PE-ECG-Research'];files=[];locations={}
 for base in roots:
  if not base.exists():continue
  for h in sorted(base.rglob('*.hea')):
   d=h.with_suffix('.dat');entry=dict(study_key=h.stem,header_path=str(h),header_sha256=sha(h),header_bytes=h.stat().st_size,data_present=d.is_file())
   if d.is_file():entry.update(data_path=str(d),data_sha256=sha(d),data_bytes=d.stat().st_size);locations.setdefault(h.stem,[]).append(str(h.with_suffix('')))
   files.append(entry)
 (ROOT/'private/WAVEFORM_FILE_HASHES.json').write_text(json.dumps(files,indent=2));(ROOT/'private/waveform_locations.json').write_text(json.dumps(locations));(ROOT/'results/WAVEFORM_SOURCE_INVENTORY.json').write_text(json.dumps(dict(discovered_headers=len(files),headers_with_data=sum(r['data_present'] for r in files),unique_paired_keys=len(locations),total_header_bytes=sum(r['header_bytes'] for r in files),total_data_bytes=sum(r.get('data_bytes',0) for r in files),all_discovered_header_and_data_files_sha256_hashed=True,full_provider_waveform_coverage=False,versions='ECGv1.0 indicated by paths; exact provider hash authenticity not established',private_manifest=True),indent=2))
 prior={}
 for name in ['POTASSIUM_SERIAL_AUDIT_V1','POTASSIUM_SERIAL_AUDIT_V2','POTASSIUM_BENCHMARK_V1','POTASSIUM_SERIAL_V3_FEASIBILITY']:
  p=ROOT.parent/name/'FILE_MANIFEST.json'
  if p.exists():prior[str(p)]=sha(p)
 (ROOT/'private/PRIOR_MANIFEST_HASHES.json').write_text(json.dumps(prior,indent=2));print('Source and local waveform hash inventory complete.')
if __name__=='__main__':main()
