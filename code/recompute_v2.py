#!/usr/bin/env python3
"""
=============================================================================
RECOMPUTE v2 — CBAC-D-26-04047
Written against the ACTUAL files in C:\\Users\\DELL\\neo_project
=============================================================================
Nothing to configure. Run it:

    python3 recompute_v2.py

It finds neo_project itself, reads your real tables, and recomputes every
number the two reviewers disputed. Output goes to ./recompute_outputs/ plus
a summary.json that I can read directly.

Needs: pandas numpy scipy   (powerlaw optional, for section E)
=============================================================================
"""
from __future__ import annotations
import sys, os, json, time
from pathlib import Path

try:
    import pandas as pd, numpy as np
    from scipy import stats
except ImportError as e:
    sys.exit(f"Missing dependency: {e}\n  pip install pandas numpy scipy")

OUT = Path("./recompute_outputs"); OUT.mkdir(exist_ok=True)
R = {}          # results collected into summary.json
NOTES = []

def hdr(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)

def save(df, name):
    df.to_csv(OUT / name, index=False)
    print(f"  -> {OUT/name}  ({len(df):,} rows)")

def find_root():
    cands = [
        Path("/mnt/c/Users/DELL/neo_project"),
        Path.home() / "neo_project",
        Path("C:/Users/DELL/neo_project"),
        Path("./neo_project"), Path("."),
    ]
    for c in cands:
        if (c / "data/output/Phase4/predictions/mhc1_all_predictions.tsv").exists():
            return c
    # last resort: search
    for base in [Path("/mnt/c/Users"), Path.home()]:
        if not base.exists(): continue
        for p in base.rglob("mhc1_all_predictions.tsv"):
            return p.parents[3]
    sys.exit("Could not locate neo_project. Run this script from inside it.")

ROOT = find_root()
print(f"neo_project root: {ROOT}")

P = {
    "mhc1":   ROOT / "data/output/Phase4/predictions/mhc1_all_predictions.tsv",
    "mhc2":   ROOT / "data/output/Phase4/mhc2/mhc2_all_predictions.tsv",
    "pepmap": ROOT / "data/output/Phase3/peptides/peptide_table.tsv",
    "mutsum": ROOT / "data/output/Phase3/analysis_data/mutation_summary.tsv",
    "muts":   ROOT / "data/Phase2_Task2/dna_repair_missense_mutations.tsv",
    "net":    ROOT / "data/output/Phase12_Advanced/epitope_network_analysis.csv",
    "score":  ROOT / "data/output/Phase12_Advanced/top_30_immunogenic_epitopes.csv",
    "genelvl":ROOT / "data/output/Phase11_True_Comparative/gene_level_patterns.csv",
    "combos": ROOT / "data/output/Phase11_True_Comparative/mutation_combinations_analysis.csv",
}
for k, v in P.items():
    print(f"  {k:8s} {'OK ' if v.exists() else 'MISSING'}  {v}")

# =============================================================================
hdr("LOADING MASTER PREDICTION TABLE")
t0 = time.time()
m1 = pd.read_csv(P["mhc1"], sep="\t", low_memory=False)
print(f"  mhc1_all_predictions.tsv: {len(m1):,} rows in {time.time()-t0:.1f}s")
print(f"  columns: {list(m1.columns)}")

PEP, MUT, GENE = "mt_peptide", "mutation_id", "gene"
IC, WT, FC = "MT_IC50_nM", "WT_IC50_nM", "Fold_Change"
HLA, LEN = "HLA", "peptide_length"

# =============================================================================
hdr("A. FILTER FLOWCHART AND THE ALLELE COUNT  [R2 figures / minor 1]")
alleles = sorted(m1[HLA].dropna().unique())
R["n_alleles"] = len(alleles); R["alleles"] = alleles
print(f"\n  HLA CLASS I ALLELES USED: {len(alleles)}")
for a in alleles: print(f"     {a}")
if len(alleles) != 7:
    NOTES.append(f"MANUSCRIPT SAYS SEVEN ALLELES. The data has {len(alleles)}.")
    print(f"\n  *** The manuscript says seven class-I alleles. This is {len(alleles)}. ***")
    print("  *** Correct the Methods, the Abstract and the Figure 1 legend.      ***")

n_pid   = m1["peptide_id"].nunique()
n_upep  = m1[PEP].nunique()
n_mut   = m1[MUT].nunique()
R.update(n_predictions=len(m1), n_peptide_windows=n_pid,
         n_unique_peptides=n_upep, n_mutations_in_predictions=n_mut)
print(f"\n  predictions (rows)        : {len(m1):,}")
print(f"  peptide windows (peptide_id): {n_pid:,}")
print(f"  unique peptide SEQUENCES  : {n_upep:,}")
print(f"  mutations                 : {n_mut:,}")
print(f"  lengths                   : {sorted(m1[LEN].dropna().unique())}")
print(f"\n  CHECK: {n_pid:,} windows x {len(alleles)} alleles = {n_pid*len(alleles):,}"
      f"   (actual rows {len(m1):,})")

for key, label in [("muts", "variants in the 29 genes (Phase2_Task2)"),
                   ("mutsum", "mutations with peptides generated (Phase3)")]:
    if P[key].exists():
        n = sum(1 for _ in open(P[key], errors="replace")) - 1
        print(f"  {label}: {n:,}")
        R[f"n_{key}"] = n

# dedup rule: is mutation_id a genomic variant or a sample-level call?
if P["mutsum"].exists():
    ms = pd.read_csv(P["mutsum"], sep="\t", low_memory=False)
    print(f"\n  mutation_summary.tsv columns: {list(ms.columns)}")
    if {"chromosome", "genomic_position", "sample_id"} <= set(ms.columns):
        n_calls = ms["mutation_id"].nunique()
        site_cols = ["gene", "chromosome", "genomic_position"]
        for extra in ("ref_aa", "alt_aa"):
            if extra in ms.columns: site_cols.append(extra)
        n_sites = ms[site_cols].drop_duplicates().shape[0]
        n_samp  = ms["sample_id"].nunique()
        R.update(n_sample_level_calls=n_calls, n_distinct_variants=n_sites,
                 n_samples=n_samp)
        print(f"\n  *** DEDUPLICATION RULE  [R2 major 6] ***")
        print(f"      sample-level calls (mutation_id) : {n_calls:,}")
        print(f"      distinct genomic variants        : {n_sites:,}")
        print(f"      tumour samples                   : {n_samp:,}")
        if n_calls > n_sites:
            print(f"\n      mutation_id is a SAMPLE-LEVEL CALL, not a distinct variant.")
            print(f"      So every 'number of mutations' in the paper — the 1,419, the")
            print(f"      240 for NWRPILTII, the 71 in Table 2 — is a RECURRENCE COUNT.")
            print(f"      That must be stated in Methods and the table footnotes.")
            NOTES.append("mutation_id = sample-level call; all mutation counts are "
                         "recurrence counts and must be relabelled.")
        ms.groupby(site_cols).size().rename("n_samples").reset_index() \
          .sort_values("n_samples", ascending=False).to_csv(
              OUT / "variant_recurrence.csv", index=False)
        print(f"  -> {OUT/'variant_recurrence.csv'}")

# =============================================================================
hdr("B. BINDER CLASSES — one table, both denominators  [R2 major 4]")
d = m1.dropna(subset=[IC])
rows = []
for label, mask in [("< 500 nM (strong)", d[IC] < 500),
                    ("500-5,000 nM (weak)", (d[IC] >= 500) & (d[IC] <= 5000)),
                    ("> 5,000 nM (non-binder)", d[IC] > 5000)]:
    s = d[mask]
    rows.append({"class": label,
                 "n_predictions": len(s),
                 "pct_of_predictions": round(len(s)/len(d)*100, 3),
                 "n_unique_peptides": s[PEP].nunique(),
                 "pct_of_unique_peptides": round(s[PEP].nunique()/n_upep*100, 3)})
bc = pd.DataFrame(rows)
print(f"\n  denominators: {len(d):,} predictions | {n_upep:,} unique peptides\n")
print(bc.to_string(index=False))
save(bc, "binder_classes.csv")
R["binder_classes"] = rows
print(f"\n  The manuscript prints 7.6% (n=3,272) and 35.7% (n=15,315).")
print(f"  Replace both with a row from this table and name the unit.")

# =============================================================================
hdr("C. DEGREE AND HUB COUNTS — one definition  [R2 major 2]")
deg = m1.groupby(PEP)[MUT].nunique().rename("degree").reset_index()
tot = len(deg)
bins = {"1": int((deg.degree == 1).sum()),
        "2": int((deg.degree == 2).sum()),
        "3-5": int(deg.degree.between(3, 5).sum()),
        "6-10": int(deg.degree.between(6, 10).sum()),
        ">10": int((deg.degree > 10).sum())}
ge3 = bins["3-5"] + bins["6-10"] + bins[">10"]
ge6 = bins["6-10"] + bins[">10"]
print(f"\n  degree = number of distinct parental mutation_ids per peptide sequence")
print(f"  unique peptides: {tot:,}   (bins sum {sum(bins.values()):,})")
for k, v in bins.items(): print(f"    {k:>5s} mutations: {v:>8,}  ({v/tot*100:5.2f}%)")
print(f"\n    degree >= 3 : {ge3:>8,}  ({ge3/tot*100:5.2f}%)")
print(f"    degree >= 6 : {ge6:>8,}  ({ge6/tot*100:5.2f}%)   <- literal 'degree > 5'")
print(f"    degree > 10 : {bins['>10']:>8,}  ({bins['>10']/tot*100:5.2f}%)  <- super-hub")
print(f"    max degree  : {int(deg.degree.max()):,}")
R.update(degree_bins=bins, n_degree_ge3=ge3, n_degree_ge6=ge6,
         n_superhub=bins[">10"], max_degree=int(deg.degree.max()),
         n_unique_peptides_network=tot)
save(deg, "degree_per_peptide.csv")
print(f"\n  Manuscript claims 836 hubs 'linked to more than five mutations' = 2.5%.")
print(f"  Pick ONE: degree>5 gives {ge6:,}; degree>=3 gives {ge3:,}.")

if P["net"].exists():
    net = pd.read_csv(P["net"])
    print(f"\n  Cross-check against epitope_network_analysis.csv ({len(net):,} rows):")
    print(f"    columns: {list(net.columns)}")
    if "Hub_Status" in net.columns:
        vc = net["Hub_Status"].value_counts().to_dict()
        print(f"    Hub_Status counts: {vc}")
        R["hub_status_counts"] = {str(k): int(v) for k, v in vc.items()}
    if "Connectivity" in net.columns:
        cn = net["Connectivity"]
        print(f"    Connectivity: max {cn.max()}, >5 = {(cn>5).sum():,}, "
              f">=3 = {(cn>=3).sum():,}, >10 = {(cn>10).sum():,}")
        R["net_connectivity_gt5"] = int((cn > 5).sum())
        R["net_connectivity_ge3"] = int((cn >= 3).sum())

# =============================================================================
hdr("D. HUB vs NON-HUB IC50  [R1 stats 2 & 4, R2 major 3]")
def cliffs(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0: return np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return 2 * (u / (len(a) * len(b))) - 1

plevel = m1.dropna(subset=[IC]).groupby(PEP).agg(
    best_ic50=(IC, "min"), median_ic50=(IC, "median"),
    mean_ic50=(IC, "mean"), gene=(GENE, lambda s: s.mode().iat[0])).reset_index()
mm = plevel.merge(deg, on=PEP)

res = []
def compare(sub, cut, col, label):
    hub = sub[sub.degree > cut][col].dropna()
    non = sub[sub.degree <= cut][col].dropna()
    if len(hub) < 3 or len(non) < 3:
        print(f"    {label}: too few (hub {len(hub)}, non-hub {len(non)})"); return
    u, p = stats.mannwhitneyu(hub, non, alternative="two-sided")
    dd = cliffs(hub, non)
    rng = np.random.default_rng(0)
    h, n_ = hub.to_numpy(), non.to_numpy()
    boot = [np.median(rng.choice(h, len(h), True)) - np.median(rng.choice(n_, len(n_), True))
            for _ in range(2000)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    bd = [cliffs(rng.choice(h, min(len(h), 3000), True),
                 rng.choice(n_, min(len(n_), 3000), True)) for _ in range(300)]
    dlo, dhi = np.nanpercentile(bd, [2.5, 97.5])
    mh, mn = float(np.median(h)), float(np.median(n_))
    print(f"\n    {label}")
    print(f"      n hub {len(hub):,} | n non-hub {len(non):,}")
    print(f"      median hub {mh:,.1f} nM | non-hub {mn:,.1f} nM | {(mh-mn)/mn*100:+.1f}%")
    print(f"      U {u:,.0f}  p {p:.3g}")
    print(f"      Cliff's delta {dd:+.3f} (95% CI {dlo:+.3f}, {dhi:+.3f})")
    print(f"      median diff 95% CI {lo:,.1f} to {hi:,.1f} nM")
    res.append(dict(stratum=label, n_hub=len(hub), n_nonhub=len(non),
                    median_hub_nM=round(mh,1), median_nonhub_nM=round(mn,1),
                    rel_diff_pct=round((mh-mn)/mn*100,1), U=float(u), p_value=float(p),
                    cliffs_delta=round(float(dd),3),
                    cliffs_ci_low=round(float(dlo),3), cliffs_ci_high=round(float(dhi),3),
                    med_diff_ci_low=round(float(lo),1), med_diff_ci_high=round(float(hi),1)))

print("\n  [D1] all peptides, hub = degree > 5")
compare(mm, 5, "best_ic50", "all peptides / best IC50 across alleles")
compare(mm, 5, "median_ic50", "all peptides / median IC50 across alleles")
print("\n  [D2] restricted to binders — the reviewer's central objection")
compare(mm[mm.best_ic50 < 5000], 5, "best_ic50", "IC50 < 5,000 nM only")
compare(mm[mm.best_ic50 < 500], 5, "best_ic50", "IC50 < 500 nM only")
print("\n  [D3] gene-stratified — separates a TP53 sequence effect from a hub effect")
for g, sub in mm.groupby("gene"):
    if len(sub) >= 50: compare(sub, 5, "best_ic50", f"gene = {g}")
print("\n  [D4] excluding TP53")
compare(mm[mm.gene != "TP53"], 5, "best_ic50", "all genes except TP53")

if res:
    t = pd.DataFrame(res)
    p = t.p_value.to_numpy(); o = np.argsort(p)
    bh = np.minimum.accumulate(((p[o]*len(p))/(np.arange(len(p))+1))[::-1])[::-1].clip(max=1)
    adj = np.empty_like(bh); adj[o] = bh
    t["p_adjusted_BH"] = np.round(adj, 6)
    t["n_comparisons"] = len(p)
    save(t, "hub_vs_nonhub_stats.csv")
    R["hub_vs_nonhub"] = res
    R["n_BH_comparisons"] = int(len(p))
    print(f"\n  [D5] {len(p)} comparisons in the BH family  [R1 stats 3]")
    print(f"       hub_vs_nonhub_stats.csv becomes Supplementary Table 7.")
    print(f"       Section 2.7 must state this number.")

# =============================================================================
hdr("E. DEGREE DISTRIBUTION — Clauset-Shalizi-Newman  [R1 stats 1]")
try:
    import powerlaw
    dv = deg.degree.to_numpy(); dv = dv[dv > 0]
    fit = powerlaw.Fit(dv, discrete=True, verbose=False)
    print(f"\n  n {len(dv):,} | alpha {fit.power_law.alpha:.3f} | "
          f"sigma {fit.power_law.sigma:.3f} | x_min {fit.power_law.xmin:.0f}")
    print(f"  n >= x_min: {int((dv >= fit.power_law.xmin).sum()):,}")
    R["powerlaw"] = {"alpha": float(fit.power_law.alpha),
                     "sigma": float(fit.power_law.sigma),
                     "xmin": float(fit.power_law.xmin)}
    comp = []
    print("\n  Likelihood ratios (R > 0 favours power law):")
    for alt in ["lognormal", "exponential", "truncated_power_law",
                "stretched_exponential"]:
        try:
            Rv, pv = fit.distribution_compare("power_law", alt, normalized_ratio=True)
            print(f"    vs {alt:<24s} R {Rv:+8.3f}  p {pv:.4g}")
            comp.append({"alternative": alt, "R": round(Rv,3), "p": float(pv)})
        except Exception as ex:
            print(f"    vs {alt:<24s} failed ({ex})")
    if comp:
        save(pd.DataFrame(comp), "powerlaw_comparisons.csv")
        R["powerlaw_comparisons"] = comp
    print("\n  Only write 'scale-free' if R > 0 and p < 0.05 against BOTH")
    print("  lognormal and exponential. Otherwise keep 'heavy-tailed'.")
except ImportError:
    print("\n  SKIPPED — pip install powerlaw, then re-run just this section.")
    NOTES.append("powerlaw not installed; 'scale-free' claim unverified.")

# =============================================================================
hdr("F. TABLE 1 — recomputed, including the four suspect cells  [R2 minor 4]")
multi = set(deg[deg.degree > 1][PEP])
g = m1.groupby(GENE).agg(unique_epitopes=(PEP, "nunique"),
                         mutations=(MUT, "nunique"),
                         predictions=(IC, "size")).reset_index()
g["epitopes_per_mutation"] = (g.unique_epitopes / g.mutations).round(2)
mmc = (m1[m1[PEP].isin(multi)].groupby(GENE)[PEP].nunique()
         .rename("multi_mutation_epitopes").reset_index())
g = g.merge(mmc, on=GENE, how="left").fillna({"multi_mutation_epitopes": 0})
g["multi_mutation_epitopes"] = g.multi_mutation_epitopes.astype(int)
g["pct_of_all_peptides"] = (g.unique_epitopes / n_upep * 100).round(2)
g["pct_of_all_mutations"] = (g.mutations / n_mut * 100).round(2)
g = g.sort_values("unique_epitopes", ascending=False)
print(); print(g.to_string(index=False))
save(g, "table1_recomputed.csv")
R["table1"] = g.to_dict("records")
print("\n  Manuscript has ATM and POLE both 228, BRCA2 and ATR both 38.")
for gn in ("ATM", "POLE", "BRCA2", "ATR"):
    r = g[g[GENE] == gn]
    if len(r):
        print(f"    {gn:<6s} multi-mutation epitopes = "
              f"{int(r.multi_mutation_epitopes.iat[0]):,}")

# =============================================================================
hdr("G. TABLE 2 — each super-hub's OWN mean IC50  [R2 major 3]")
T2 = ["MTEVVRHCPHH", "HMTEVVRHCPH", "QHMTEVVRHCP", "SQHMTEVVRH", "QSQHMTEVVRH"]
sub = m1[m1[PEP].isin(T2)]
if sub.empty:
    top = deg.nlargest(8, "degree")[PEP].tolist()
    print(f"  Table 2 sequences not found; using highest-degree peptides instead.")
    sub = m1[m1[PEP].isin(top)]
if not sub.empty:
    t2 = sub.groupby(PEP).agg(
        gene=(GENE, lambda s: ",".join(sorted(set(s.astype(str))))),
        n_mutations=(MUT, "nunique"), n_predictions=(IC, "size"),
        mean_ic50=(IC, "mean"), median_ic50=(IC, "median"),
        min_ic50=(IC, "min"), max_ic50=(IC, "max"),
        mean_fc=(FC, "mean")).round(2).reset_index()
    print(); print(t2.to_string(index=False))
    save(t2, "table2_recomputed.csv")
    R["table2"] = t2.to_dict("records")
    nb = (t2.mean_ic50 > 5000).sum()
    print(f"\n  {nb} of {len(t2)} exceed the 5,000 nM non-binder cutoff.")
    print("  The submitted table gave all five an identical 35,315 nM, which is the")
    print("  mutation-COMBINATION mean from mutation_combinations_analysis.csv,")
    print("  not a per-peptide value. Use the column above.")

# =============================================================================
hdr("H. SHARE COUNTS — what 240 and 71 actually mean  [R2 major 6]")
sc = deg.rename(columns={"degree": "n_distinct_mutation_ids"}).copy()
if P["mutsum"].exists() and "sample_id" in ms.columns:
    m2s = m1[[PEP, MUT]].drop_duplicates().merge(
        ms[["mutation_id", "sample_id"]], left_on=MUT, right_on="mutation_id", how="left")
    ns = m2s.groupby(PEP)["sample_id"].nunique().rename("n_distinct_samples").reset_index()
    sc = sc.merge(ns, on=PEP, how="left")
    site_cols = [c for c in ["gene","chromosome","genomic_position"] if c in ms.columns]
    m2v = m1[[PEP, MUT]].drop_duplicates().merge(
        ms[["mutation_id"] + site_cols], left_on=MUT, right_on="mutation_id", how="left")
    nv = (m2v.groupby(PEP)[site_cols].nunique().max(axis=1)
            .rename("n_distinct_genomic_sites").reset_index())
    sc = sc.merge(nv, on=PEP, how="left")
print(); print(sc.nlargest(15, "n_distinct_mutation_ids").to_string(index=False))
save(sc, "share_counts.csv")
for q in ("NWRPILTII", "GMNWRPILTI"):
    row = sc[sc[PEP] == q]
    if len(row):
        print(f"\n  {q}: {row.to_dict('records')[0]}")
        R[f"share_{q}"] = row.to_dict("records")[0]

# =============================================================================
hdr("I. STRICT FILTER — does 382 MHC-I + 157 MHC-II reproduce?  [S2.6]")
st1 = m1[(m1[IC] < 50) & (m1[FC] >= 10) & (m1[WT] > 5000)]
print(f"\n  MHC-I  IC50<50, FC>=10, WT>5000:")
print(f"    {len(st1):,} predictions | {st1[PEP].nunique():,} unique peptides "
      f"| {st1['peptide_id'].nunique():,} windows")
R["strict_mhc1"] = {"predictions": len(st1), "unique_peptides": int(st1[PEP].nunique()),
                    "windows": int(st1["peptide_id"].nunique())}
print(f"    manuscript says 382")

if P["mhc2"].exists():
    m2 = pd.read_csv(P["mhc2"], sep="\t", low_memory=False)
    print(f"\n  MHC-II table: {len(m2):,} rows | columns {list(m2.columns)}")
    a2 = sorted(m2[HLA].dropna().unique())
    print(f"  CLASS II ALLELES ({len(a2)}):")
    for a in a2: print(f"     {a}")
    L2 = sorted(m2[PEP].astype(str).str.len().unique())
    print(f"  peptide lengths: {L2}")
    st2 = m2[(m2[IC] < 50) & (m2[FC] >= 10) & (m2[WT] > 5000)]
    print(f"\n  MHC-II strict: {len(st2):,} predictions | "
          f"{st2[PEP].nunique():,} unique peptides")
    print(f"    manuscript says 157")
    R["mhc2"] = {"n_predictions": len(m2), "n_alleles": len(a2), "alleles": a2,
                 "lengths": [int(x) for x in L2],
                 "strict_predictions": len(st2),
                 "strict_unique_peptides": int(st2[PEP].nunique()),
                 "n_unique_peptides": int(m2[PEP].nunique()),
                 "n_mutations": int(m2[MUT].nunique())}
    print(f"\n  *** This is the material for the missing Methods 2.3b. ***")
    bands2 = []
    for label, mask in [("< 500 nM", m2[IC] < 500),
                        ("500-5,000 nM", (m2[IC] >= 500) & (m2[IC] <= 5000)),
                        ("> 5,000 nM", m2[IC] > 5000)]:
        s = m2[mask]
        bands2.append({"class": label, "n_predictions": len(s),
                       "n_unique_peptides": s[PEP].nunique()})
    save(pd.DataFrame(bands2), "mhc2_binder_classes.csv")
    save(st2, "mhc2_strict_candidates.csv")

# =============================================================================
hdr("J. COMPOSITE SCORE — empirical tail  [R2 minor 5]")
if P["score"].exists():
    sc30 = pd.read_csv(P["score"])
    print(f"\n  top_30_immunogenic_epitopes.csv: {len(sc30)} rows")
    print(f"  columns: {list(sc30.columns)}")
    print(sc30.to_string(index=False))
    R["score_file_rows"] = len(sc30)
    print(f"\n  The manuscript reports 14 top peptides, mean 62.4, SD 18.7,")
    print(f"  'approximately normal', and 2,147 peptides above 85.")
    print(f"  This file has {len(sc30)} rows, so the full 33,405-peptide score")
    print(f"  table is NOT saved anywhere. Two options:")
    print(f"    (a) find the scoring script and re-run it over all peptides, or")
    print(f"    (b) drop the 'n = 2,147 above 85' claim and the normality claim,")
    print(f"        and present the ranked list as a top-N without a tier cut.")
    NOTES.append("Full composite-score table for all 33,405 peptides not found; "
                 "only the top-N file exists. The >85 tier count cannot be verified.")

# =============================================================================
hdr("K. PER-GENE YIELD — give 'above the cohort mean' a test  [R1 stats 6]")
per = m1.groupby([GENE, MUT])[PEP].nunique().rename("peptides").reset_index()
rows = []
for gn, sub in per.groupby(GENE):
    if len(sub) < 5: continue
    rest = per[per[GENE] != gn].peptides
    u, p = stats.mannwhitneyu(sub.peptides, rest, alternative="two-sided")
    rng = np.random.default_rng(1)
    boot = [np.mean(rng.choice(sub.peptides.to_numpy(), len(sub), True)) for _ in range(2000)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    rows.append({"gene": gn, "n_mutations": len(sub),
                 "mean_peptides_per_mutation": round(sub.peptides.mean(), 2),
                 "ci_low": round(lo, 2), "ci_high": round(hi, 2),
                 "cohort_mean": round(per.peptides.mean(), 2),
                 "p_value": float(p),
                 "cliffs_delta": round(float(cliffs(sub.peptides, rest)), 3)})
t = pd.DataFrame(rows).sort_values("mean_peptides_per_mutation", ascending=False)
p = t.p_value.to_numpy(); o = np.argsort(p)
bh = np.minimum.accumulate(((p[o]*len(p))/(np.arange(len(p))+1))[::-1])[::-1].clip(max=1)
adj = np.empty_like(bh); adj[o] = bh
t["p_adjusted_BH"] = np.round(adj, 6)
print(); print(t.to_string(index=False))
save(t, "per_gene_yield_tests.csv")
R["per_gene_yield"] = t.to_dict("records")

# =============================================================================
hdr("L. EXPRESSION / FPKM — what was actually applied  [R1 repro 5]")
exp_dir = ROOT / "data/raw/Task_1B"
if exp_dir.exists():
    for f in sorted(exp_dir.glob("*")):
        try:
            n = sum(1 for _ in open(f, errors="replace"))
            first = open(f, errors="replace").readline().strip()[:200]
            print(f"\n  {f.name}  [{n} lines]")
            print(f"    {first}")
        except Exception: pass
    print(f"""
  The expression tables here cover far fewer genes than 29 and use anonymised
  Sample_N columns rather than TCGA barcodes. The Methods claim that FPKM > 1
  was applied per mutation cannot be supported by these files. Describe what
  was actually done: which genes had expression data, at what level (gene mean
  across a sample subset, not per-patient), and say plainly that expression was
  used as a gene-level inclusion screen rather than a per-mutation filter.""")
    NOTES.append("FPKM filter: expression tables cover ~12 genes / 51 anonymised "
                 "samples, so Methods cannot claim a per-patient FPKM>1 filter "
                 "across 29 genes.")

# =============================================================================
hdr("M. MISSING PREDICTIONS  [R1 stats 8]")
print(f"\n  rows with no MT_IC50 : {m1[IC].isna().sum():,}")
print(f"  rows with no WT_IC50 : {m1[WT].isna().sum():,}")
print(f"  rows with no FC      : {m1[FC].isna().sum():,}")
exp_rows = n_pid * len(alleles)
print(f"  expected rows if complete: {exp_rows:,}  actual: {len(m1):,}  "
      f"difference {exp_rows - len(m1):,}")
R["missing"] = {"mt_ic50_na": int(m1[IC].isna().sum()),
                "wt_ic50_na": int(m1[WT].isna().sum()),
                "fc_na": int(m1[FC].isna().sum()),
                "expected_rows": int(exp_rows), "actual_rows": int(len(m1))}
print("""
  Note also: this table carries ONE MT_IC50_nM per peptide-HLA pair, with no
  per-tool columns. So the 'median across three tools' in Section 2.3 cannot be
  audited from this file, and Reviewer 2's point that pVACseq wraps NetMHCpan
  and MHCflurry still has to be answered from the pipeline script.""")

# =============================================================================
hdr("SUMMARY")
R["notes"] = NOTES
(OUT / "summary.json").write_text(json.dumps(R, indent=2, default=str))
print(f"\n  -> {OUT/'summary.json'}")
if NOTES:
    print("\n  THINGS THAT NEED A DECISION OR A TEXT CHANGE:")
    for n in NOTES: print(f"    - {n}")
print(f"\n  All outputs in {OUT.resolve()}")
print("""
  Send me summary.json first — it is small and carries every headline number.
  Then the CSVs if I ask for them.""")
