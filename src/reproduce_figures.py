from pathlib import Path
import json,hashlib,ast
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
(R/"figures").mkdir(exist_ok=True)
import figure_flow
R=Path(__file__).resolve().parents[1]
for n,h in json.loads((R/'results/SAVED_RESULTS_LOCK.json').read_text()).items():
 assert hashlib.sha256((R/'results'/n).read_bytes()).hexdigest()==h,n
E=pd.read_csv(R/'results/COMPLETE_SENSITIVITY_TABLE.csv');assert len(E)==14
assert not json.loads((R/'results/FIXED_PAIR_PRECISION_CHECK.json').read_text())['material_problem']
labels=['Primary ±30 min','Both gaps ≤15 min','60-min frozen envelope*','Both labs at/before ECG','Both labs at/after ECG','Exclude hemolysis mentions','Exclude source-time conflicts*']
f,axs=plt.subplots(1,2,figsize=(11,5.3),sharey=True)
for ax,(group,g) in zip(axs,E.groupby('group',sort=False)):
 for i,r in enumerate(g.itertuples()):
  lo,hi=ast.literal_eval(r.CI95);ax.errorbar(r.effect,i,xerr=[[r.effect-lo],[hi-r.effect]],fmt='o',color='#114B5F' if i==0 else '#555555',capsize=3,ms=5)
 ax.axvline(0,color='black',lw=.8,ls='--');ax.set_title(group.capitalize());ax.set_xlabel('Difference in absolute change (log units)');ax.set_yticks(range(7),labels);ax.set_xlim(-.13,.28);ax.grid(axis='x',alpha=.2);ax.tick_params(labelleft=True,labelsize=9)
axs[0].invert_yaxis();f.tight_layout(w_pad=2)
for ext in ['png','pdf']:f.savefig(R/f'figures/Figure_1_Complete_Effects.{ext}',dpi=300,bbox_inches='tight')
plt.close(f)
print('Aggregate figures regenerated; no effects or bootstrap rerun.')
