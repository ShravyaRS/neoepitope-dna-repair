#!/usr/bin/env python3
"""
=============================================================================
EXACT composite score — transcribed verbatim from PHASE12_WOW_FACTOR.py
=============================================================================
  python3 score_exact.py
Then send me  score_exact_report.txt

Verified by hand against two published values before writing this:
  HYNYMCNSF  -> 91.176 (published 91.175980)
  NWRPILTII  -> 85.305 (published 85.304847)
=============================================================================
"""
import sys
from pathlib import Path
import pandas as pd, numpy as np

OUT = open("score_exact_report.txt", "w", encoding="utf-8")
def say(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.write(s + "\n")
def hdr(t):
    say("\n" + "=" * 76); say(t); say("=" * 76)

def find_root():
    for c in [Path("/mnt/c/Users/DELL/neo_project"), Path.home()/"neo_project", Path(".")]:
        if (c/"data/output/Phase4/predictions/mhc1_all_predictions.tsv").exists():
            return c
    sys.exit("neo_project not found")
R = find_root()

df = pd.read_csv(R/"data/output/Phase4/predictions/mhc1_all_predictions.tsv",
                 sep="\t", low_memory=False)
say(f"root: {R}")
say(f"loaded {len(df):,} predictions")

# ---------------------------------------------------- EXACT Phase12 function
GENE_IMPORTANCE = {'TP53':1.0,'BRCA1':0.95,'BRCA2':0.95,'ATM':0.9,
                   'MSH2':0.85,'MLH1':0.85,'POLE':0.8,'ATR':0.8}

def compute(df, use_ic=True, use_fc=True, use_share=True, use_gene=True):
    s = np.zeros(len(df))
    ic = df['MT_IC50_nM']
    ic50_norm = 1 - (ic - ic.min()) / (ic.max() - ic.min())
    if use_ic:    s = s + ic50_norm * 40
    fold_norm = np.clip(df['Fold_Change'].fillna(1.0), 0, 100) / 100
    if use_fc:    s = s + fold_norm * 20
    counts = df['mt_peptide'].value_counts()
    sharing = df['mt_peptide'].map(counts)
    share_norm = np.clip(sharing, 1, 100) / 100
    if use_share: s = s + share_norm * 20
    gene_score = df['gene'].map(GENE_IMPORTANCE).fillna(0.7)
    if use_gene:  s = s + gene_score * 20
    return s, ic50_norm, fold_norm, share_norm, gene_score, sharing

df['score'], ic50_norm, fold_norm, share_norm, gene_score, sharing = compute(df)
df['epitope_sharing'] = sharing

hdr("A. VERIFICATION AGAINST THE PUBLISHED TOP-30 FILE")
ref = R/"data/output/Phase12_Advanced/top_30_immunogenic_epitopes.csv"
if ref.exists():
    t = pd.read_csv(ref)
    say(f"\n  published file has {len(t)} rows\n")
    worst = 0
    for _, r in t.iterrows():
        m = df[(df.mt_peptide == r.mt_peptide) &
               (np.isclose(df.MT_IC50_nM, r.MT_IC50_nM))]
        if len(m):
            v = float(m.score.iloc[0]); d = abs(v - r.immunogenicity_score)
            worst = max(worst, d)
            say(f"    {r.mt_peptide:<12s} published {r.immunogenicity_score:9.5f} "
                f"| rebuilt {v:9.5f} | diff {d:.6f}")
    say(f"\n  worst difference: {worst:.6f}")
    say("  -> function confirmed exactly" if worst < 0.01
        else "  -> MISMATCH, tell me and I will re-derive")

hdr("B. SCORE DISTRIBUTION — what the manuscript reports vs reality")
for label, s in [("per prediction row (n=429,464)", df.score),
                 ("per unique peptide, best row",
                  df.loc[df.groupby('mt_peptide').score.idxmax()].score)]:
    say(f"\n  {label}:  n = {len(s):,}")
    say(f"    mean {s.mean():.2f} | SD {s.std():.2f} | median {s.median():.2f}")
    say(f"    min {s.min():.2f} | max {s.max():.2f}")
    for cut in (75, 80, 85, 90):
        say(f"    n > {cut}: {int((s>cut).sum()):>8,}  ({(s>cut).mean()*100:6.2f}%)")
say("\n  MANUSCRIPT CLAIMS: mean 62.4, SD 18.7, score > 85 -> n = 2,147")

hdr("C. WHY THE RANKING BEHAVES AS IT DOES (component degeneracy)")
say(f"\n  IC50 range used for min-max: {df.MT_IC50_nM.min():.2f} to {df.MT_IC50_nM.max():.2f} nM")
say(f"  s(IC50) term, out of 40:")
for q in (0.01, 0.25, 0.5, 0.75, 0.99):
    say(f"    {int(q*100):>3d}th pct of predictions: {np.quantile(ic50_norm*40, q):6.3f}")
for thr in (50, 500, 5000):
    sub = ic50_norm[df.MT_IC50_nM < thr] * 40
    if len(sub):
        say(f"    peptides with IC50 < {thr:>5,} nM: term ranges "
            f"{sub.min():.3f} to {sub.max():.3f}")
say("""
  If that range is narrow, the affinity term barely discriminates among
  binders and the ranking is driven by the other three terms.""")
say(f"\n  s(FC): {int((df.Fold_Change.fillna(1)>=100).sum()):,} of {len(df):,} "
    f"predictions ({(df.Fold_Change.fillna(1)>=100).mean()*100:.2f}%) hit the 100 cap")
say(f"  s(share): {int((sharing>=100).sum()):,} predictions "
    f"({(sharing>=100).mean()*100:.2f}%) hit the 100 cap")
say(f"  gene weights present: {sorted(df.gene.map(GENE_IMPORTANCE).fillna(0.7).unique())}")
say(f"  genes at the 1.0 maximum: {sorted(df[df.gene.map(GENE_IMPORTANCE).fillna(0.7)==1.0].gene.unique())}")

hdr("D. TOP 14 — reproducing the published table")
df['rank'] = df.score.rank(ascending=False)
pub = (df.nsmallest(100, 'rank')[['mt_peptide','gene','mutation_id','MT_IC50_nM',
        'Fold_Change','epitope_sharing','score']]
       .drop_duplicates('mt_peptide').head(30))
say(f"\n  replicating the exact Phase12 recipe "
    f"(nsmallest(100,'rank') -> drop_duplicates -> head(30)):")
say(f"  yields {len(pub)} distinct peptides")
say(pub.to_string(index=False))
pub.to_csv("exact_top30_replicated.csv", index=False)

say("\n  TOP 20 unique peptides by best score (the defensible ranking):")
best = (df.loc[df.groupby('mt_peptide').score.idxmax()]
          .nlargest(20, 'score')[['mt_peptide','gene','HLA','MT_IC50_nM',
                                  'Fold_Change','epitope_sharing','score']])
say(best.to_string(index=False))
best.to_csv("exact_top20_unique.csv", index=False)
for n in (14, 100, 1000):
    t = df.loc[df.groupby('mt_peptide').score.idxmax()].nlargest(n, 'score')
    say(f"\n  top {n}: {t.gene.nunique()} gene(s) -> {dict(t.gene.value_counts().head(8))}")

hdr("E. REVIEWER 2 COMMENT 8 — ranking with s(gene) and s(share) removed")
variants = {
    "full 40/20/20/20 (as published)": dict(),
    "without s(gene)":                 dict(use_gene=False),
    "without s(share)":                dict(use_share=False),
    "without s(gene) and s(share)":    dict(use_gene=False, use_share=False),
    "s(IC50) only":                    dict(use_fc=False, use_share=False, use_gene=False),
}
base_top = None; rows = []
for label, kw in variants.items():
    s, *_ = compute(df, **kw)
    d = df.assign(v=s)
    top14 = (d.loc[d.groupby('mt_peptide').v.idxmax()]
               .nlargest(14, 'v')[['mt_peptide','gene','v']])
    peps = top14.mt_peptide.tolist()
    if base_top is None: base_top = peps
    ov = len(set(peps) & set(base_top))
    say(f"\n  --- {label}")
    say(f"      overlap with published top 14: {ov}/14   "
        f"genes: {dict(top14.gene.value_counts())}")
    say(top14.to_string(index=False))
    rows.append({"variant": label, "overlap_with_published_top14": f"{ov}/14",
                 "n_TP53_in_top14": int((top14.gene=='TP53').sum()),
                 "n_genes_in_top14": int(top14.gene.nunique())})
pd.DataFrame(rows).to_csv("weight_sensitivity_exact.csv", index=False)
say("\n  -> weight_sensitivity_exact.csv  (this becomes Supplementary Table 8)")

OUT.close()
print("\nWrote score_exact_report.txt + 3 CSVs")
