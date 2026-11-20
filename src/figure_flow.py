from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
R=Path(__file__).resolve().parents[1]

f,axs=plt.subplots(1,2,figsize=(10,6.3),gridspec_kw={'width_ratios':[1.15,1]})
for ax in axs:ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
left=['Full ECG index\n800,035 ECGs; 161,352 patients','Index–machine time agreement\n782,185 ECGs','Adult eligibility\n782,079 ECGs','Corroborated encounter / ED check\n399,847 ECGs','Nearest valid lab within 30 min\n74,662 ECGs','Recorded matching admission and\nunambiguous specimen: 33,440 ECGs']
for i,t in enumerate(left):
 y=.92-i*.155;axs[0].add_patch(FancyBboxPatch((.04,y-.105),.92,.105,boxstyle='round,pad=.009',facecolor='#edf2f7',edgecolor='#6b7280'));axs[0].text(.5,y-.052,t,ha='center',va='center',fontsize=9)
 if i<5:axs[0].annotate('',xy=(.5,y-.145),xytext=(.5,y-.112),arrowprops={'arrowstyle':'->','color':'#555'})
axs[0].set_title('A. Source linkage',fontsize=12)
right=['Fixed metadata pairs\nHyperkalemia 126; hypokalemia 211\nStable controls 682','All selected records locally available\n2,037 unique ECGs\n1,990 pass individual QC','Pair-level exclusions\nSignal QC: 2 / 4 / 8\nExtraction: 2 / 7 / 24','Final shared-stratum pairs\nHyperkalemia 122; hypokalemia 200\nStable controls 650\n967 distinct patients']
for i,t in enumerate(right):
 y=.92-i*.23;axs[1].add_patch(FancyBboxPatch((.04,y-.16),.92,.16,boxstyle='round,pad=.009',facecolor='#edf2f7',edgecolor='#6b7280'));axs[1].text(.5,y-.08,t,ha='center',va='center',fontsize=9)
 if i<3:axs[1].annotate('',xy=(.5,y-.21),xytext=(.5,y-.17),arrowprops={'arrowstyle':'->','color':'#555'})
axs[1].set_title('B. Fixed-pair quality control',fontsize=12)
f.tight_layout();
for ext in ['png','pdf']:f.savefig(R/f'figures/Figure_2_Cohort_Flow.{ext}',dpi=300,bbox_inches='tight')
plt.close(f)
print('Presentation polished without changing analysis numbers.')
