import sys
from pathlib import Path
import pandas as pd, numpy as np
for c in [Path("/mnt/c/Users/DELL/neo_project"), Path.home()/"neo_project", Path(".")]:
    if (c/"data/output/Phase4/predictions/mhc1_all_predictions.tsv").exists(): R=c; break
else: sys.exit("not found")
m=pd.read_csv(R/"data/output/Phase4/predictions/mhc1_all_predictions.tsv",sep="\t",low_memory=False)
fc=m.Fold_Change
print("Does 24.6% (FC>=10) or 78.3% (FC>1) appear under ANY denominator?\n")
maxfc=m.groupby('mt_peptide').Fold_Change.max()
winfc=m.groupby('peptide_id').Fold_Change.max()
rows=[
 ("all 429,464 predictions",                fc>=10,  len(m)),
 ("all 429,464 predictions, FC>1",          fc>1,    len(m)),
 ("33,405 unique peptides (any allele)",    maxfc>=10, len(maxfc)),
 ("33,405 unique peptides, FC>1",           maxfc>1,   len(maxfc)),
 ("53,683 windows (any allele)",            winfc>=10, len(winfc)),
 ("53,683 windows, FC>1",                   winfc>1,   len(winfc)),
]
for lab,mask,n in rows:
    print(f"  {lab:<42s} {mask.sum():>8,} / {n:>8,} = {mask.sum()/n*100:6.2f}%")
print("\n  restricted to binders:")
for thr in (5000,500,50):
    s=m[m.MT_IC50_nM<thr]
    if len(s):
        mx=s.groupby('mt_peptide').Fold_Change.max()
        print(f"    IC50<{thr:>5,}: predictions FC>=10 {(s.Fold_Change>=10).mean()*100:6.2f}% | "
              f"unique peptides {(mx>=10).mean()*100:6.2f}%  (n={len(mx):,})")
print(f"\n  median FC = {fc.median():.2f} | mean = {fc.mean():.2f}")
print(f"  FC>=2: {(maxfc>=2).mean()*100:.1f}%  FC>=5: {(maxfc>=5).mean()*100:.1f}%  "
      f"FC>=10: {(maxfc>=10).mean()*100:.1f}%  FC>=20: {(maxfc>=20).mean()*100:.1f}%")
