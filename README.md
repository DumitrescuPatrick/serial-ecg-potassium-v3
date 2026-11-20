# Serial ECG–Potassium Study: Public Code

Code and aggregate results for the exploratory same-source serial ECG–potassium V3 study. This is separate from the prediction benchmark. The original study and manuscript remain unchanged.

## Contents

- `src/reproduce_figures.py`: regenerate the complete 14-row forest plot and cohort-flow figure from saved aggregates; no statistical refitting.
- `src/historical_methods/`: preserved cohort, waveform-extraction, QC and effect source references. These require authorized MIMIC sources, private intermediate files and the original project layout. They are not a turnkey public raw-data pipeline.
- `src/check_fixed_pairs.py`: finishing check used with authorized local V3 private inputs; no effect estimation. It expects the original sibling POTASSIUM_SERIAL_V3 project.
- `results/`: aggregate results, fixed-pair potassium distributions and numerical precision summaries; no individual records.
- `docs/`: dated protocols and reference verification.

## Run without patient data

Use Python 3.12 and an isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python src/reproduce_figures.py
```

Figures are written to `figures/`. The commands do not rerun effect estimation or bootstrap intervals. The reported contrasts remain 0.0522 (95% CI −0.0210 to 0.1321) and 0.0196 (−0.0182 to 0.0582) log units; all 14 intervals include zero. Fixed groups comprise 122 hyperkalemia, 200 hypokalemia and 650 control pairs, with 967 distinct patients. Five patients occur in both case groups.

## Interpretation and access

This is exploratory work informed by earlier analyses, not independent confirmation. Clock synchronization, treatment timing, source hospital/ED release versions and independent clinical waveform interpretation remain unresolved. No clinical utility, equivalence or novel biomarker claim is supported. See protocols for the redundant 60-minute envelope and source-time exclusion rows.

No patient-level records, identifiers, waveforms, waveform paths, private manifests, individual predictions or credentials are distributed. MIMIC data require authorized access and must not be added to this repository. Source-level reproduction requires the original protected workflow and inputs, not just this public archive. AI-assisted preparation requires complete author verification and disclosure before publication. Authors must approve code ownership and a license; no license grant is invented here.

Upload the contents of this folder to your GitHub repository. Do not upload the original study folders or their private directories.
