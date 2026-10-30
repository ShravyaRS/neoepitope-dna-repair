#!/usr/bin/env python3
"""
Supplementary Figure 1 -- fold-change distribution.
  python3 make_suppfig.py        ->  ./figures_v2/Supplementary_Figure_1.png
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
    "figure.dpi": 300, "savefig.dpi": 300,
})
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#8C8C8C"

for c in [Path("/mnt/c/Users/DELL/neo_project"), Path.home()/"neo_project", Path(".")]:
    if (c/"data/output/Phase4/predictions/mhc1_all_predictions.tsv").exists(): R=c; break
else: sys.exit("neo_project not found")
OUT = Path("./figures_v2"); OUT.mkdir(exist_ok=True)

m = pd.read_csv(R/"data/output/Phase4/predictions/mhc1_all_predictions.tsv", sep="\t", low_memory=False)
fc_pep = m.groupby('mt_peptide').Fold_Change.max()
n = len(fc_pep)
p_gt1  = (fc_pep > 1).mean()*100
p_ge10 = (fc_pep >= 10).mean()*100
med    = fc_pep.median()
sub500 = m[m.MT_IC50_nM < 500].groupby('mt_peptide').Fold_Change.max()
p500   = (sub500 >= 10).mean()*100

fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))

# (A) distribution on log axis
v = fc_pep[fc_pep > 0]
bins = np.logspace(np.log10(max(v.min(), 1e-3)), np.log10(v.max()), 60)
ax[0].hist(v, bins=bins, color=GREY, edgecolor="white", linewidth=0.4)
ax[0].hist(v[v >= 10], bins=bins, color=ORANGE, edgecolor="white", linewidth=0.4,
           label=f"FC $\\geq$ 10  ({p_ge10:.1f}%)")
ax[0].axvline(1,  color=BLUE,  lw=1.4, ls="--")
ax[0].axvline(10, color=ORANGE, lw=1.4, ls="--")
ax[0].set_xscale("log")
ax[0].set_xlabel("Fold change  (IC50 wild-type / IC50 mutant)")
ax[0].set_ylabel(f"Unique peptides  (n = {n:,})")
ax[0].set_title("A. Fold-change distribution", loc="left", fontweight="bold", fontsize=11)
ax[0].annotate(f"median = {med:.2f}", xy=(med, ax[0].get_ylim()[1]*0.88),
               xytext=(med*2.2, ax[0].get_ylim()[1]*0.88), color=BLUE, fontsize=9,
               arrowprops=dict(arrowstyle="->", color=BLUE, lw=1))
ax[0].legend(frameon=False, fontsize=9, loc="upper right")

# (B) cumulative fractions at thresholds
th  = [1, 2, 5, 10, 20, 50, 100]
val = [(fc_pep >= t).mean()*100 for t in th]
bars = ax[1].bar([str(t) for t in th], val,
                 color=[BLUE if t < 10 else ORANGE for t in th],
                 edgecolor="white", linewidth=0.8)
for b, vv in zip(bars, val):
    ax[1].text(b.get_x()+b.get_width()/2, vv + max(val)*0.02,
               f"{vv:.1f}%", ha="center", fontsize=8.5)
ax[1].axhline(p500, color=GREEN, lw=1.4, ls=":")
ax[1].text(len(th)-0.4, p500 + max(val)*0.03,
           f"FC $\\geq$ 10 among peptides binding < 500 nM: {p500:.1f}%",
           ha="right", fontsize=8.5, color=GREEN)
ax[1].set_xlabel("Fold-change threshold")
ax[1].set_ylabel("Unique peptides at or above threshold (%)")
ax[1].set_title("B. Peptides reaching each threshold", loc="left", fontweight="bold", fontsize=11)
ax[1].set_ylim(0, max(val)*1.18)

fig.tight_layout()
fig.savefig(OUT/"Supplementary_Figure_1.png", bbox_inches="tight")
print(f"wrote {OUT/'Supplementary_Figure_1.png'}")
print(f"  n = {n:,} | FC>1 {p_gt1:.1f}% | FC>=10 {p_ge10:.1f}% | median {med:.2f} | <500nM FC>=10 {p500:.1f}%")
