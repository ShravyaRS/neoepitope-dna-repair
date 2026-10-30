#!/usr/bin/env python3
"""Find which filter gives 73 MHC-I peptides, and where 1,419 vs 1,417 differs."""
import sys
from pathlib import Path
import pandas as pd, numpy as np

def find_root():
    for c in [Path("/mnt/c/Users/DELL/neo_project"), Path.home()/"neo_project", Path(".")]:
        if (c/"data/output/Phase4/predictions/mhc1_all_predictions.tsv").exists(): return c
    sys.exit("neo_project not found")
R = find_root()
m1 = pd.read_csv(R/"data/output/Phase4/predictions/mhc1_all_predictions.tsv", sep="\t", low_memory=False)
ic, fc, wt = m1.MT_IC50_nM.values, m1.Fold_Change.values, m1.WT_IC50_nM.values

print("="*72); print("A. WHICH FILTER GIVES 73 UNIQUE MHC-I PEPTIDES?"); print("="*72)
tests = [
 ("IC50<50  & FC>=10 & WT>5000  (as the caption states)", (ic<50)&(fc>=10)&(wt>5000)),
 ("IC50<50  & FC>=10",                                     (ic<50)&(fc>=10)),
 ("IC50<50  & FC>=10 & WT>500",                            (ic<50)&(fc>=10)&(wt>500)),
 ("IC50<50  & FC>=5",                                      (ic<50)&(fc>=5)),
 ("IC50<100 & FC>=10 & WT>5000",                           (ic<100)&(fc>=10)&(wt>5000)),
 ("IC50<100 & FC>=10",                                     (ic<100)&(fc>=10)),
 ("IC50<150 & FC>=10 & WT>5000",                           (ic<150)&(fc>=10)&(wt>5000)),
 ("IC50<500 & FC>=10 & WT>5000",                           (ic<500)&(fc>=10)&(wt>5000)),
 ("IC50<50  & FC>=20 & WT>5000",                           (ic<50)&(fc>=20)&(wt>5000)),
 ("IC50<50  & FC>=10 & WT>1000",                           (ic<50)&(fc>=10)&(wt>1000)),
]
for label, mask in tests:
    n_u = m1[mask].mt_peptide.nunique(); n_w = m1[mask].peptide_id.nunique()
    flag = "   <<< 73" if n_u == 73 else ("   <<< 73 windows" if n_w == 73 else "")
    print(f"  {label:<48s} unique {n_u:>5,} | windows {n_w:>5,} | rows {int(mask.sum()):>6,}{flag}")

print("\n" + "="*72); print("B. 1,419 vs 1,417"); print("="*72)
print(f"  distinct mutation_id in class I predictions : {m1.mutation_id.nunique():,}")
for name, p in [("peptide_table", "data/output/Phase3/peptides/peptide_table.tsv"),
                ("mutation_summary", "data/output/Phase3/analysis_data/mutation_summary.tsv")]:
    f = R/p
    if f.exists():
        t = pd.read_csv(f, sep="\t", low_memory=False)
        col = 'mutation_id' if 'mutation_id' in t.columns else t.columns[0]
        print(f"  distinct {col} in {name:<18s}: {t[col].nunique():,}  (rows {len(t):,})")
        if name == "peptide_table":
            missing = set(t[col].unique()) - set(m1.mutation_id.unique())
            print(f"  in peptide_table but NOT in predictions: {len(missing)} -> {sorted(missing)[:10]}")
            if missing:
                sub = t[t[col].isin(missing)]
                cols = [c for c in ['mutation_id','gene','mt_peptide','length'] if c in sub.columns]
                print(sub[cols].head(12).to_string(index=False))
