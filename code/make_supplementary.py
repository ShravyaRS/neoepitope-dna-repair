#!/usr/bin/env python3
"""
=============================================================================
BUILD SUPPLEMENTARY TABLES S1-S7  --  CBAC-D-26-04047
=============================================================================
  python3 make_supplementary.py

Writes ./supplementary/Supplementary_Table_S1.csv ... S7.csv
S8 and S9 already exist (weight_sensitivity_exact.csv, zenodo_deposit/SCHEMA.csv)
and are copied in automatically if found.

Takes a few minutes; the betweenness step is the slow one.
=============================================================================
"""
import sys, shutil
from pathlib import Path
import numpy as np, pandas as pd

try:
    import networkx as nx
except ImportError:
    sys.exit("pip install networkx")
from scipy.stats import mannwhitneyu

OUT = Path("./supplementary"); OUT.mkdir(exist_ok=True)

def find_root():
    for c in [Path("/mnt/c/Users/DELL/neo_project"), Path.home()/"neo_project", Path(".")]:
        if (c/"data/output/Phase4/predictions/mhc1_all_predictions.tsv").exists():
            return c
    sys.exit("neo_project not found")
R = find_root(); print("root:", R, "\n")

GENE_W = {'TP53':1.0,'BRCA1':0.95,'BRCA2':0.95,'ATM':0.9,
          'MSH2':0.85,'MLH1':0.85,'POLE':0.8,'ATR':0.8}

print("loading class I predictions ...")
m1 = pd.read_csv(R/"data/output/Phase4/predictions/mhc1_all_predictions.tsv",
                 sep="\t", low_memory=False)
print(f"  {len(m1):,} rows")

# ---- composite score, exactly as PHASE12_WOW_FACTOR.py ------------------
ic = m1['MT_IC50_nM']
m1['s_IC50']  = (1 - (ic - ic.min())/(ic.max() - ic.min())) * 40
m1['s_FC']    = np.clip(m1['Fold_Change'].fillna(1.0), 0, 100)/100 * 20
rowcount      = m1['mt_peptide'].value_counts()
m1['pred_rows'] = m1['mt_peptide'].map(rowcount)
m1['s_share'] = np.clip(m1['pred_rows'], 1, 100)/100 * 20
m1['s_gene']  = m1['gene'].map(GENE_W).fillna(0.7) * 20
m1['composite_score'] = m1.s_IC50 + m1.s_FC + m1.s_share + m1.s_gene

# ---- bipartite graph: unique peptide sequence <-> parental call ---------
print("building bipartite graph ...")
edges = m1[['mt_peptide','mutation_id']].drop_duplicates()
G = nx.Graph()
G.add_edges_from(zip("P_"+edges.mt_peptide, "M_"+edges.mutation_id.astype(str)))
deg = {p[2:]: d for p, d in G.degree() if p.startswith("P_")}
m1['degree'] = m1.mt_peptide.map(deg)
print(f"  {len(deg):,} unique peptides, {edges.mutation_id.nunique():,} parental calls")

best = m1.loc[m1.groupby('mt_peptide').composite_score.idxmax()].copy()
print(f"  {len(best):,} unique peptides after best-row selection")

# ========================= S1 ===========================================
s1 = best[['mt_peptide','gene','mutation_id','HLA','MT_IC50_nM','WT_IC50_nM',
           'Fold_Change','degree','composite_score']].copy()
s1.columns = ['peptide','gene','mutation_id','HLA','MT_IC50_nM','WT_IC50_nM',
              'fold_change','sharing_frequency_mutations','composite_score']
s1 = s1.sort_values('composite_score', ascending=False)
s1.to_csv(OUT/"Supplementary_Table_S1.csv", index=False)
print(f"S1  {len(s1):,} peptides")

# ========================= S2 ===========================================
src = R/"data/output/Phase11_True_Comparative/multi_mutation_epitopes_detailed.csv"
if src.exists():
    shutil.copy(src, OUT/"Supplementary_Table_S2.csv")
    print(f"S2  {len(pd.read_csv(src)):,} multi-mutation peptides (copied)")
else:
    multi = best[best.degree > 1]
    multi.to_csv(OUT/"Supplementary_Table_S2.csv", index=False)
    print(f"S2  {len(multi):,} multi-mutation peptides (rebuilt)")

# ========================= S3 ===========================================
g = m1.groupby('mutation_id')
s3 = pd.DataFrame({
    'gene':                g.gene.first(),
    'n_peptide_windows':   g.peptide_id.nunique(),
    'n_unique_peptides':   g.mt_peptide.nunique(),
    'n_HLA_alleles':       g.HLA.nunique(),
    'mean_IC50_nM':        g.MT_IC50_nM.mean().round(1),
    'median_IC50_nM':      g.MT_IC50_nM.median().round(1),
    'min_IC50_nM':         g.MT_IC50_nM.min().round(1),
    'max_IC50_nM':         g.MT_IC50_nM.max().round(1),
    'n_strong_binders_lt50nM':  g.MT_IC50_nM.apply(lambda s: int((s<50).sum())),
    'pct_strong_binders':  (g.MT_IC50_nM.apply(lambda s: (s<50).mean())*100).round(3),
    'mean_fold_change':    g.Fold_Change.mean().round(2),
}).reset_index()
s3.to_csv(OUT/"Supplementary_Table_S3.csv", index=False)
print(f"S3  {len(s3):,} mutations")

# ========================= S4 ===========================================
hubs = sorted([p for p, d in deg.items() if d > 5])
print(f"S4  {len(hubs):,} hubs -- computing network metrics (slow step) ...")
comp = {}
for i, c in enumerate(nx.connected_components(G)):
    for n in c: comp[n] = i
btw = nx.betweenness_centrality(G, k=min(1000, G.number_of_nodes()), seed=42)
# clustering on the hub-induced peptide projection
mut2hub = {}
for _, r in edges.iterrows():
    if r.mt_peptide in set(hubs):
        mut2hub.setdefault(r.mutation_id, []).append(r.mt_peptide)
P = nx.Graph(); P.add_nodes_from(hubs)
for peps in mut2hub.values():
    for a in range(len(peps)):
        for b in range(a+1, len(peps)):
            P.add_edge(peps[a], peps[b])
clus = nx.clustering(P)
bmap = best.set_index('mt_peptide')
s4 = pd.DataFrame({
    'peptide': hubs,
    'gene':    [bmap.gene.get(p) for p in hubs],
    'degree':  [deg[p] for p in hubs],
    'betweenness_centrality_sampled_k1000': [round(btw.get("P_"+p, 0.0), 8) for p in hubs],
    'clustering_coefficient_hub_projection': [round(clus.get(p, 0.0), 6) for p in hubs],
    'connected_component_id': [comp.get("P_"+p) for p in hubs],
    'best_MT_IC50_nM': [bmap.MT_IC50_nM.get(p) for p in hubs],
    'composite_score': [round(bmap.composite_score.get(p), 4) for p in hubs],
}).sort_values('degree', ascending=False)
s4.to_csv(OUT/"Supplementary_Table_S4.csv", index=False)
print(f"S4  {len(s4):,} hub peptides")

# ========================= S5 ===========================================
s5 = best[['mt_peptide','gene','mutation_id','HLA','MT_IC50_nM','Fold_Change',
           'degree','s_IC50','s_FC','s_share','s_gene','composite_score']].copy()
s5.columns = ['peptide','gene','mutation_id','HLA','MT_IC50_nM','fold_change',
              'sharing_frequency_mutations','s_IC50_of_40','s_FC_of_20',
              's_share_of_20','s_gene_of_20','composite_score']
s5 = s5.sort_values('composite_score', ascending=False).round(4)
s5.to_csv(OUT/"Supplementary_Table_S5.csv", index=False)
print(f"S5  {len(s5):,} peptides with components")

# ========================= S6 ===========================================
f1 = m1[(m1.MT_IC50_nM < 50) & (m1.Fold_Change >= 10)]
f1b = f1.loc[f1.groupby('mt_peptide').MT_IC50_nM.idxmin()]
f1b = f1b[['mt_peptide','gene','mutation_id','HLA','MT_IC50_nM','WT_IC50_nM',
           'Fold_Change','degree','composite_score']].copy()
f1b['clears_WT_above_5000nM'] = f1b.WT_IC50_nM > 5000
f1b.insert(0, 'class', 'MHC-I')
rows = [f1b]
m2p = R/"data/output/Phase4/mhc2/mhc2_all_predictions.tsv"
if m2p.exists():
    m2 = pd.read_csv(m2p, sep="\t", low_memory=False)
    f2 = m2[(m2.MT_IC50_nM < 100) & (m2.Fold_Change >= 10)]
    f2b = f2.loc[f2.groupby('mt_peptide').MT_IC50_nM.idxmin()]
    keep = [c for c in ['mt_peptide','gene','mutation_id','HLA','MT_IC50_nM',
                        'WT_IC50_nM','Fold_Change'] if c in f2b.columns]
    f2b = f2b[keep].copy(); f2b.insert(0, 'class', 'MHC-II')
    rows.append(f2b)
    print(f"S6  MHC-I {f1b.mt_peptide.nunique()} (of which {int(f1b.clears_WT_above_5000nM.sum())} clear WT>5000) | MHC-II {f2b.mt_peptide.nunique()}")
else:
    print(f"S6  MHC-I {f1b.mt_peptide.nunique()} | MHC-II file not found")
pd.concat(rows, ignore_index=True).to_csv(OUT/"Supplementary_Table_S6.csv", index=False)

# ========================= S7 ===========================================
print("S7  running the five comparisons ...")
bestic = m1.groupby('mt_peptide').MT_IC50_nM.min()
medic  = m1.groupby('mt_peptide').MT_IC50_nM.median()
gene_of = m1.groupby('mt_peptide').gene.first()
dser   = pd.Series(deg)
ishub  = dser > 5

def cliffs(a, b, nboot=2000, seed=0):
    a = np.asarray(a); b = np.asarray(b)
    def d(x, y):
        x = np.sort(x)
        gt = np.searchsorted(x, y, 'left').sum()
        lt = (len(x)*len(y)) - np.searchsorted(x, y, 'right').sum()
        return (lt - gt) / (len(x)*len(y))
    est = d(a, b)
    rng = np.random.default_rng(seed); boots = []
    for _ in range(nboot):
        boots.append(d(rng.choice(a, min(len(a), 2000)),
                       rng.choice(b, min(len(b), 2000))))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return est, lo, hi

COMPS = [
 ("Unstratified, best IC50 across alleles", bestic, bestic.index),
 ("Unstratified, median IC50 across alleles", medic, medic.index),
 ("Restricted to peptides with IC50 < 5,000 nM", bestic, bestic[bestic < 5000].index),
 ("Restricted to peptides with IC50 < 500 nM",   bestic, bestic[bestic < 500].index),
 ("Restricted to TP53 peptides",                 bestic, gene_of[gene_of == 'TP53'].index),
]
res = []
for label, series, idx in COMPS:
    idx = [p for p in idx if p in dser.index]
    h  = series.loc[[p for p in idx if ishub.get(p, False)]].dropna().values
    nh = series.loc[[p for p in idx if not ishub.get(p, False)]].dropna().values
    if len(h) < 3 or len(nh) < 3:
        res.append({'comparison': label, 'n_hub': len(h), 'n_nonhub': len(nh),
                    'note': 'too few observations'}); continue
    u, p = mannwhitneyu(h, nh, alternative='two-sided')
    dl, lo, hi = cliffs(h, nh)
    res.append({
        'comparison': label, 'n_hub': len(h), 'n_nonhub': len(nh),
        'median_hub_nM': round(float(np.median(h)), 1),
        'median_nonhub_nM': round(float(np.median(nh)), 1),
        'pct_difference': round((np.median(h)-np.median(nh))/np.median(nh)*100, 2),
        'cliffs_delta': round(dl, 4), 'ci_low': round(lo, 4), 'ci_high': round(hi, 4),
        'mannwhitney_U': float(u), 'p_raw': float(p)})
s7 = pd.DataFrame(res)
if 'p_raw' in s7.columns:                       # Benjamini-Hochberg
    ok = s7.p_raw.notna()
    pv = s7.loc[ok, 'p_raw'].values
    order = np.argsort(pv); m = len(pv)
    adj = np.empty(m); prev = 1.0
    for rank in range(m-1, -1, -1):
        i = order[rank]
        prev = min(prev, pv[i]*m/(rank+1)); adj[i] = prev
    s7.loc[ok, 'p_adjusted_BH'] = adj
    s7['significant_at_FDR_5pct'] = s7.p_adjusted_BH < 0.05
s7.to_csv(OUT/"Supplementary_Table_S7.csv", index=False)
print(s7.to_string(index=False))

# ---- bring in S8 / S9 if they are lying around --------------------------
for cand, dest in [
    (Path.home()/"weight_sensitivity_exact.csv", "Supplementary_Table_S8.csv"),
    (Path("/mnt/c/Users/DELL/Downloads/weight_sensitivity_exact.csv"), "Supplementary_Table_S8.csv"),
    (Path("/mnt/c/Users/DELL/Downloads/zenodo_deposit/SCHEMA.csv"), "Supplementary_Table_S9.csv"),
]:
    if cand.exists() and not (OUT/dest).exists():
        shutil.copy(cand, OUT/dest); print(f"copied {dest}")

print("\n" + "="*70)
for f in sorted(OUT.glob("*.csv")):
    print(f"  {f.name:<34s} {f.stat().st_size/1024:9.1f} KB")
print("="*70)
