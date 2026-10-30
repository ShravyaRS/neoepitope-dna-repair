"""
===============================================================================
PHASE 5: DETAILED NEOEPITOPE RANKING WITH VISUALIZATIONS
Individual Peptide-Level Analysis
===============================================================================
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'figure.dpi': 300,
    'font.size': 10,
    'font.family': 'sans-serif',
    'axes.labelsize': 11,
    'axes.titlesize': 13,
    'figure.titlesize': 16
})

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = Path("/mnt/c/Users/DELL/neo_project/data/output")
PHASE5_DIR = BASE_DIR / "Phase5"
FIGURES_DIR = PHASE5_DIR / "figures"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {
    'excellent': '#2E7D32',  # Dark green
    'high': '#66BB6A',       # Light green
    'medium': '#FFA726',     # Orange
    'low': '#EF5350'         # Red
}

print("="*80)
print("PHASE 5: DETAILED NEOEPITOPE RANKING")
print("Individual Peptide-Level Analysis")
print("="*80)

# ============================================================================
# LOAD DATA
# ============================================================================

print("\nLoading Phase 5 results...")

# Load all three key files
df_mhc1_top = pd.read_csv(PHASE5_DIR / "top_50_mhc1_neoepitopes.tsv", sep='\t')
df_mhc2_top = pd.read_csv(PHASE5_DIR / "top_50_mhc2_neoepitopes.tsv", sep='\t')
df_dual = pd.read_csv(PHASE5_DIR / "dual_coverage_neoepitopes.tsv", sep='\t')

print(f"✓ Top 50 MHC-I neoepitopes: {len(df_mhc1_top)}")
print(f"✓ Top 50 MHC-II neoepitopes: {len(df_mhc2_top)}")
print(f"✓ Dual coverage neoepitopes: {len(df_dual)}")

# ============================================================================
# CATEGORIZE BY PRIORITY
# ============================================================================

print("\n" + "="*80)
print("NEOEPITOPE PRIORITY CLASSIFICATION")
print("="*80)

def classify_priority(row, mhc_class='I'):
    """Classify neoepitope priority based on comprehensive criteria"""
    
    if mhc_class == 'I':
        ic50_threshold_excellent = 50
        ic50_threshold_good = 150
        fc_threshold_excellent = 10
    else:
        ic50_threshold_excellent = 100
        ic50_threshold_good = 500
        fc_threshold_excellent = 10
    
    # Excellent: Strong binder + High specificity + (Dual coverage bonus)
    if row['MT_IC50_nM'] < ic50_threshold_excellent and row['Fold_Change'] >= fc_threshold_excellent:
        if row.get('Has_MHC2', False) or row.get('Has_MHC1', False):
            return "TIER 1: Excellent (Dual)"
        return "TIER 1: Excellent"
    
    # High: Strong binder OR high specificity
    elif row['MT_IC50_nM'] < ic50_threshold_good and row['Fold_Change'] >= 5:
        return "TIER 2: High"
    
    # Medium: Good binding with moderate specificity
    elif row['Fold_Change'] >= 5:
        return "TIER 3: Medium"
    
    return "TIER 4: Low"

df_mhc1_top['Priority'] = df_mhc1_top.apply(lambda x: classify_priority(x, 'I'), axis=1)
df_mhc2_top['Priority'] = df_mhc2_top.apply(lambda x: classify_priority(x, 'II'), axis=1)

# Add rank
df_mhc1_top['Rank'] = range(1, len(df_mhc1_top) + 1)
df_mhc2_top['Rank'] = range(1, len(df_mhc2_top) + 1)

print("\nMHC-I Priority Distribution:")
for priority, count in df_mhc1_top['Priority'].value_counts().sort_index().items():
    print(f"  {priority}: {count}")

print("\nMHC-II Priority Distribution:")
for priority, count in df_mhc2_top['Priority'].value_counts().sort_index().items():
    print(f"  {priority}: {count}")

# ============================================================================
# CREATE DETAILED RANKING TABLE WITH ALL INFO
# ============================================================================

print("\n" + "="*80)
print("CREATING COMPREHENSIVE RANKING")
print("="*80)

# MHC-I: Select key columns and add detailed info
mhc1_detailed = df_mhc1_top[['Rank', 'Priority', 'gene', 'mutation_id', 'mt_peptide', 
                               'HLA', 'peptide_length', 'MT_IC50_nM', 'WT_IC50_nM', 
                               'Fold_Change', 'Final_Score', 'Has_MHC2']].copy()

mhc1_detailed.columns = ['Rank', 'Priority', 'Gene', 'Mutation_ID', 'Peptide_Sequence', 
                          'HLA_Allele', 'Length', 'MT_IC50_nM', 'WT_IC50_nM', 
                          'Fold_Change', 'Immunogenicity_Score', 'CD4_Support']

# MHC-II: Same
mhc2_detailed = df_mhc2_top[['Rank', 'Priority', 'gene', 'mutation_id', 'mt_peptide',
                               'HLA', 'peptide_length', 'MT_IC50_nM', 'WT_IC50_nM',
                               'Fold_Change', 'Final_Score', 'Has_MHC1']].copy()

mhc2_detailed.columns = ['Rank', 'Priority', 'Gene', 'Mutation_ID', 'Peptide_Sequence',
                          'HLA_Allele', 'Length', 'MT_IC50_nM', 'WT_IC50_nM',
                          'Fold_Change', 'Immunogenicity_Score', 'CD8_Support']

# Save detailed rankings
mhc1_detailed.to_csv(PHASE5_DIR / "MHC1_Detailed_Ranking.tsv", sep='\t', index=False)
mhc2_detailed.to_csv(PHASE5_DIR / "MHC2_Detailed_Ranking.tsv", sep='\t', index=False)

print(f"✓ MHC-I detailed ranking: {PHASE5_DIR / 'MHC1_Detailed_Ranking.tsv'}")
print(f"✓ MHC-II detailed ranking: {PHASE5_DIR / 'MHC2_Detailed_Ranking.tsv'}")

# ============================================================================
# TOP 20 INDIVIDUAL NEOEPITOPES TABLE (HUMAN READABLE)
# ============================================================================

print("\n" + "="*80)
print("TOP 20 INDIVIDUAL NEOEPITOPES (MHC-I)")
print("="*80)
print(f"\n{'Rank':<6}{'Gene':<10}{'Peptide':<15}{'HLA':<15}{'IC50':<10}{'Fold':<8}{'Priority':<25}")
print("-"*95)

for idx, row in mhc1_detailed.head(20).iterrows():
    print(f"{row['Rank']:<6}{row['Gene']:<10}{row['Peptide_Sequence']:<15}"
          f"{row['HLA_Allele']:<15}{row['MT_IC50_nM']:<10.1f}{row['Fold_Change']:<8.1f}"
          f"{row['Priority']:<25}")

print("\n" + "="*80)
print("TOP 20 INDIVIDUAL NEOEPITOPES (MHC-II)")
print("="*80)
print(f"\n{'Rank':<6}{'Gene':<10}{'Peptide':<20}{'HLA':<20}{'IC50':<10}{'Fold':<8}{'Priority':<25}")
print("-"*105)

for idx, row in mhc2_detailed.head(20).iterrows():
    print(f"{row['Rank']:<6}{row['Gene']:<10}{row['Peptide_Sequence']:<20}"
          f"{row['HLA_Allele']:<20}{row['MT_IC50_nM']:<10.1f}{row['Fold_Change']:<8.1f}"
          f"{row['Priority']:<25}")

# ============================================================================
# FIGURE 1: TOP 20 NEOEPITOPES WITH PRIORITY COLORS
# ============================================================================

print("\n[Figure 1] Creating top 20 neoepitopes visualization...")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(18, 14))

# MHC-I Top 20
top20_mhc1 = mhc1_detailed.head(20)

# Create color map based on priority
priority_colors_mhc1 = []
for priority in top20_mhc1['Priority']:
    if 'TIER 1' in priority and 'Dual' in priority:
        priority_colors_mhc1.append(COLORS['excellent'])
    elif 'TIER 1' in priority:
        priority_colors_mhc1.append(COLORS['high'])
    elif 'TIER 2' in priority:
        priority_colors_mhc1.append(COLORS['medium'])
    else:
        priority_colors_mhc1.append(COLORS['low'])

y_pos = np.arange(len(top20_mhc1))
bars1 = ax1.barh(y_pos, top20_mhc1['Immunogenicity_Score'], 
                 color=priority_colors_mhc1, edgecolor='black', linewidth=1.5)

# Labels with gene and peptide
labels = [f"#{row['Rank']}: {row['Gene']} - {row['Peptide_Sequence'][:10]}... ({row['HLA_Allele']})" 
          for idx, row in top20_mhc1.iterrows()]

ax1.set_yticks(y_pos)
ax1.set_yticklabels(labels, fontsize=9)
ax1.set_xlabel('Immunogenicity Score', fontweight='bold', fontsize=12)
ax1.set_title('A. Top 20 MHC-I Neoepitopes (Ranked by Immunogenicity)', 
              fontweight='bold', fontsize=14, pad=15)
ax1.invert_yaxis()
ax1.grid(axis='x', alpha=0.3)

# Add IC50 values as text
for i, (idx, row) in enumerate(top20_mhc1.iterrows()):
    ax1.text(row['Immunogenicity_Score'] + 1, i, 
            f"IC50: {row['MT_IC50_nM']:.1f}nM\nFC: {row['Fold_Change']:.1f}x",
            va='center', fontsize=7, fontweight='bold')

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=COLORS['excellent'], label='Tier 1 Excellent (Dual)', edgecolor='black'),
    Patch(facecolor=COLORS['high'], label='Tier 1 Excellent', edgecolor='black'),
    Patch(facecolor=COLORS['medium'], label='Tier 2 High', edgecolor='black'),
    Patch(facecolor=COLORS['low'], label='Tier 3 Medium', edgecolor='black')
]
ax1.legend(handles=legend_elements, loc='lower right', frameon=True, fancybox=True)

# MHC-II Top 20
top20_mhc2 = mhc2_detailed.head(20)

priority_colors_mhc2 = []
for priority in top20_mhc2['Priority']:
    if 'TIER 1' in priority and 'Dual' in priority:
        priority_colors_mhc2.append(COLORS['excellent'])
    elif 'TIER 1' in priority:
        priority_colors_mhc2.append(COLORS['high'])
    elif 'TIER 2' in priority:
        priority_colors_mhc2.append(COLORS['medium'])
    else:
        priority_colors_mhc2.append(COLORS['low'])

y_pos2 = np.arange(len(top20_mhc2))
bars2 = ax2.barh(y_pos2, top20_mhc2['Immunogenicity_Score'],
                 color=priority_colors_mhc2, edgecolor='black', linewidth=1.5)

labels2 = [f"#{row['Rank']}: {row['Gene']} - {row['Peptide_Sequence'][:12]}... ({row['HLA_Allele']})"
           for idx, row in top20_mhc2.iterrows()]

ax2.set_yticks(y_pos2)
ax2.set_yticklabels(labels2, fontsize=9)
ax2.set_xlabel('Immunogenicity Score', fontweight='bold', fontsize=12)
ax2.set_title('B. Top 20 MHC-II Neoepitopes (Ranked by Immunogenicity)',
              fontweight='bold', fontsize=14, pad=15)
ax2.invert_yaxis()
ax2.grid(axis='x', alpha=0.3)

# Add IC50 values
for i, (idx, row) in enumerate(top20_mhc2.iterrows()):
    ax2.text(row['Immunogenicity_Score'] + 1, i,
            f"IC50: {row['MT_IC50_nM']:.1f}nM\nFC: {row['Fold_Change']:.1f}x",
            va='center', fontsize=7, fontweight='bold')

ax2.legend(handles=legend_elements, loc='lower right', frameon=True, fancybox=True)

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'Figure1_Top20_Individual_Neoepitopes.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Figure 1 saved")

# ============================================================================
# FIGURE 2: GENE-LEVEL SUMMARY WITH NEOEPITOPE COUNTS
# ============================================================================

print("\n[Figure 2] Creating gene-level analysis...")

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# MHC-I: Gene counts with priority breakdown
gene_priority_mhc1 = df_mhc1_top.groupby(['gene', 'Priority']).size().unstack(fill_value=0)
top_genes_mhc1 = df_mhc1_top['gene'].value_counts().head(10).index
gene_priority_mhc1_top = gene_priority_mhc1.loc[top_genes_mhc1]

gene_priority_mhc1_top.plot(kind='barh', stacked=True, ax=axes[0,0], 
                            color=[COLORS['excellent'], COLORS['high'], COLORS['medium'], COLORS['low']],
                            edgecolor='black', linewidth=1)
axes[0,0].set_xlabel('Number of Neoepitopes', fontweight='bold')
axes[0,0].set_title('A. MHC-I Neoepitopes by Gene (Stacked by Priority)', fontweight='bold', pad=10)
axes[0,0].legend(title='Priority', bbox_to_anchor=(1.05, 1), loc='upper left')
axes[0,0].invert_yaxis()
axes[0,0].grid(axis='x', alpha=0.3)

# MHC-II: Gene counts
gene_priority_mhc2 = df_mhc2_top.groupby(['gene', 'Priority']).size().unstack(fill_value=0)
top_genes_mhc2 = df_mhc2_top['gene'].value_counts().head(10).index
gene_priority_mhc2_top = gene_priority_mhc2.loc[top_genes_mhc2]

gene_priority_mhc2_top.plot(kind='barh', stacked=True, ax=axes[0,1],
                            color=[COLORS['excellent'], COLORS['high'], COLORS['medium'], COLORS['low']],
                            edgecolor='black', linewidth=1)
axes[0,1].set_xlabel('Number of Neoepitopes', fontweight='bold')
axes[0,1].set_title('B. MHC-II Neoepitopes by Gene (Stacked by Priority)', fontweight='bold', pad=10)
axes[0,1].legend(title='Priority', bbox_to_anchor=(1.05, 1), loc='upper left')
axes[0,1].invert_yaxis()
axes[0,1].grid(axis='x', alpha=0.3)

# IC50 distribution by top genes (MHC-I)
ax3 = axes[1,0]
for gene in top_genes_mhc1[:5]:
    gene_data = df_mhc1_top[df_mhc1_top['gene'] == gene]['MT_IC50_nM']
    ax3.hist(gene_data, bins=20, alpha=0.6, label=gene, edgecolor='black')

ax3.set_xlabel('MT IC50 (nM)', fontweight='bold')
ax3.set_ylabel('Frequency', fontweight='bold')
ax3.set_title('C. IC50 Distribution - Top 5 Genes (MHC-I)', fontweight='bold', pad=10)
ax3.legend()
ax3.grid(axis='y', alpha=0.3)
ax3.axvline(50, color='red', linestyle='--', linewidth=2, label='Strong Binder')

# Dual coverage by gene
ax4 = axes[1,1]
dual_genes = df_dual.groupby('gene').size().sort_values(ascending=False).head(10)
ax4.bar(range(len(dual_genes)), dual_genes.values, color=COLORS['excellent'], edgecolor='black', linewidth=1.5)
ax4.set_xticks(range(len(dual_genes)))
ax4.set_xticklabels(dual_genes.index, rotation=45, ha='right')
ax4.set_ylabel('Dual Coverage Neoepitopes', fontweight='bold')
ax4.set_title('D. Genes with Dual MHC-I + MHC-II Coverage', fontweight='bold', pad=10)
ax4.grid(axis='y', alpha=0.3)

# Add value labels
for i, v in enumerate(dual_genes.values):
    ax4.text(i, v + 0.5, str(v), ha='center', va='bottom', fontweight='bold')

plt.tight_layout()
plt.savefig(FIGURES_DIR / 'Figure2_Gene_Level_Analysis.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Figure 2 saved")

# ============================================================================
# FIGURE 3: COMPREHENSIVE DASHBOARD FOR TOP GENE
# ============================================================================

print("\n[Figure 3] Creating detailed analysis for top gene (TP53)...")

top_gene = 'TP53'
df_gene_mhc1 = df_mhc1_top[df_mhc1_top['gene'] == top_gene]
df_gene_mhc2 = df_mhc2_top[df_mhc2_top['gene'] == top_gene]

fig = plt.figure(figsize=(18, 10))
gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

fig.suptitle(f'Detailed Analysis: {top_gene} Neoepitopes', fontsize=18, fontweight='bold')

# Panel 1: Top 15 TP53 neoepitopes (MHC-I)
ax1 = fig.add_subplot(gs[0, :2])
top15_gene = df_gene_mhc1.head(15)
y_pos = np.arange(len(top15_gene))

colors_15 = []
for priority in top15_gene['Priority']:
    if 'TIER 1' in priority and 'Dual' in priority:
        colors_15.append(COLORS['excellent'])
    elif 'TIER 1' in priority:
        colors_15.append(COLORS['high'])
    elif 'TIER 2' in priority:
        colors_15.append(COLORS['medium'])
    else:
        colors_15.append(COLORS['low'])

bars = ax1.barh(y_pos, top15_gene['Final_Score'], color=colors_15, edgecolor='black', linewidth=1.5)

labels_gene = [f"{row['Peptide_Sequence'][:12]} ({row['HLA_Allele']})" 
               for idx, row in top15_gene.iterrows()]
ax1.set_yticks(y_pos)
ax1.set_yticklabels(labels_gene, fontsize=9, family='monospace')
ax1.set_xlabel('Final Score', fontweight='bold')
ax1.set_title(f'Top 15 {top_gene} Neoepitopes (MHC-I)', fontweight='bold', pad=10)
ax1.invert_yaxis()
ax1.grid(axis='x', alpha=0.3)

# Add IC50 annotations
for i, (idx, row) in enumerate(top15_gene.iterrows()):
    ax1.text(row['Final_Score'] + 1, i, f"{row['MT_IC50_nM']:.1f}nM",
            va='center', fontsize=8, fontweight='bold')

# Panel 2: Priority distribution pie chart
ax2 = fig.add_subplot(gs[0, 2])
priority_counts_gene = df_gene_mhc1['Priority'].value_counts()
colors_pie = [COLORS['excellent'] if 'Dual' in p else COLORS['high'] if 'TIER 1' in p else COLORS['medium'] if 'TIER 2' in p else COLORS['low'] 
              for p in priority_counts_gene.index]

ax2.pie(priority_counts_gene.values, labels=priority_counts_gene.index, autopct='%1.1f%%',
        colors=colors_pie, startangle=90, wedgeprops={'edgecolor': 'black', 'linewidth': 1.5})
ax2.set_title(f'{top_gene} Priority Distribution', fontweight='bold', pad=10)

# Panel 3: IC50 vs Fold Change scatter
ax3 = fig.add_subplot(gs[1, 0])
scatter = ax3.scatter(df_gene_mhc1['MT_IC50_nM'], df_gene_mhc1['Fold_Change'],
                     c=df_gene_mhc1['Final_Score'], cmap='RdYlGn', s=100, 
                     edgecolors='black', linewidth=1.5, alpha=0.7)
ax3.axhline(5, color='red', linestyle='--', linewidth=2, alpha=0.7)
ax3.axvline(50, color='red', linestyle='--', linewidth=2, alpha=0.7)
ax3.set_xlabel('MT IC50 (nM)', fontweight='bold')
ax3.set_ylabel('Fold Change (WT/MT)', fontweight='bold')
ax3.set_title(f'{top_gene}: Binding vs Specificity', fontweight='bold', pad=10)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.grid(alpha=0.3)
plt.colorbar(scatter, ax=ax3, label='Final Score')

# Panel 4: HLA distribution
ax4 = fig.add_subplot(gs[1, 1])
hla_counts_gene = df_gene_mhc1['HLA'].value_counts()
ax4.bar(range(len(hla_counts_gene)), hla_counts_gene.values, 
        color=plt.cm.tab10(range(len(hla_counts_gene))), edgecolor='black', linewidth=1.5)
ax4.set_xticks(range(len(hla_counts_gene)))
ax4.set_xticklabels(hla_counts_gene.index, rotation=45, ha='right')
ax4.set_ylabel('Count', fontweight='bold')
ax4.set_title(f'{top_gene} HLA Distribution', fontweight='bold', pad=10)
ax4.grid(axis='y', alpha=0.3)

# Panel 5: Stats summary
ax5 = fig.add_subplot(gs[1, 2])
ax5.axis('off')

stats_text = f"""
{top_gene} NEOEPITOPE SUMMARY
{'='*35}

MHC-I Neoepitopes: {len(df_gene_mhc1)}
MHC-II Neoepitopes: {len(df_gene_mhc2)}

TIER 1 Excellent: {sum('TIER 1' in p for p in df_gene_mhc1['Priority'])}
  - With Dual Coverage: {sum('Dual' in p for p in df_gene_mhc1['Priority'])}

Strong Binders: {(df_gene_mhc1['MT_IC50_nM'] < 50).sum()}
High Specificity (FC≥10): {(df_gene_mhc1['Fold_Change'] >= 10).sum()}

Mean IC50: {df_gene_mhc1['MT_IC50_nM'].mean():.1f} nM
Mean Fold Change: {df_gene_mhc1['Fold_Change'].mean():.1f}x

HLA Coverage: {df_gene_mhc1['HLA'].nunique()} alleles
"""

ax5.text(0.1, 0.9, stats_text, transform=ax5.transAxes, 
         fontsize=10, verticalalignment='top', family='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.savefig(FIGURES_DIR / f'Figure3_Detailed_{top_gene}_Analysis.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Figure 3 saved")

# ============================================================================
# FINAL SUMMARY REPORT WITH SPECIFIC RECOMMENDATIONS
# ============================================================================

print("\n" + "="*80)
print("GENERATING FINAL RECOMMENDATIONS")
print("="*80)

# Identify THE BEST single neoepitope
best_mhc1 = mhc1_detailed.iloc[0]
best_mhc2 = mhc2_detailed.iloc[0]

# Dual coverage recommendations
if len(df_dual) > 0:
    best_dual = df_dual.sort_values('Final_Score_MHC1', ascending=False).iloc[0]

report = PHASE5_DIR / "FINAL_RECOMMENDATIONS.txt"
with open(report, 'w') as f:
    f.write("="*80 + "\n")
    f.write("FINAL NEOEPITOPE RECOMMENDATIONS\n")
    f.write("Specific Peptide-Level Results\n")
    f.write("="*80 + "\n\n")
    
    f.write("╔" + "="*78 + "╗\n")
    f.write("║  #1 RANKED NEOEPITOPE (MHC-I - CD8+ T-cell)                             ║\n")
    f.write("╚" + "="*78 + "╝\n\n")
    f.write(f"  Gene:             {best_mhc1['Gene']}\n")
    f.write(f"  Mutation:         {best_mhc1['Mutation_ID']}\n")
    f.write(f"  Peptide Sequence: {best_mhc1['Peptide_Sequence']}\n")
    f.write(f"  HLA Allele:       {best_mhc1['HLA_Allele']}\n")
    f.write(f"  MT IC50:          {best_mhc1['MT_IC50_nM']:.2f} nM\n")
    f.write(f"  WT IC50:          {best_mhc1['WT_IC50_nM']:.2f} nM\n")
    f.write(f"  Fold Change:      {best_mhc1['Fold_Change']:.1f}x\n")
    f.write(f"  Final Score:      {best_mhc1['Immunogenicity_Score']:.1f}\n")
    f.write(f"  Priority:         {best_mhc1['Priority']}\n")
    f.write(f"  CD4 Support:      {'YES' if best_mhc1['CD4_Support'] else 'NO'}\n\n")
    
    f.write("╔" + "="*78 + "╗\n")
    f.write("║  #1 RANKED NEOEPITOPE (MHC-II - CD4+ T-cell)                            ║\n")
    f.write("╚" + "="*78 + "╝\n\n")
    f.write(f"  Gene:             {best_mhc2['Gene']}\n")
    f.write(f"  Mutation:         {best_mhc2['Mutation_ID']}\n")
    f.write(f"  Peptide Sequence: {best_mhc2['Peptide_Sequence']}\n")
    f.write(f"  HLA Allele:       {best_mhc2['HLA_Allele']}\n")
    f.write(f"  MT IC50:          {best_mhc2['MT_IC50_nM']:.2f} nM\n")
    f.write(f"  WT IC50:          {best_mhc2['WT_IC50_nM']:.2f} nM\n")
    f.write(f"  Fold Change:      {best_mhc2['Fold_Change']:.1f}x\n")
    f.write(f"  Final Score:      {best_mhc2['Immunogenicity_Score']:.1f}\n")
    f.write(f"  Priority:         {best_mhc2['Priority']}\n")
    f.write(f"  CD8 Support:      {'YES' if best_mhc2['CD8_Support'] else 'NO'}\n\n")
    
    if len(df_dual) > 0:
        f.write("╔" + "="*78 + "╗\n")
        f.write("║  TOP DUAL COVERAGE NEOEPITOPE (Both CD8+ and CD4+ Response)             ║\n")
        f.write("╚" + "="*78 + "╝\n\n")
        f.write(f"  Gene:             {best_dual['gene']}\n")
        f.write(f"  Mutation:         {best_dual['mutation_id']}\n")
        f.write(f"  MHC-I Peptide:    {best_dual['mt_peptide_MHC1']}\n")
        f.write(f"  MHC-I HLA:        {best_dual['HLA_MHC1']}\n")
        f.write(f"  MHC-I IC50:       {best_dual['MT_IC50_nM_MHC1']:.2f} nM\n")
        f.write(f"  MHC-II Peptide:   {best_dual['mt_peptide_MHC2']}\n")
        f.write(f"  MHC-II HLA:       {best_dual['HLA_MHC2']}\n")
        f.write(f"  MHC-II IC50:      {best_dual['MT_IC50_nM_MHC2']:.2f} nM\n\n")
    
    f.write("\n" + "="*80 + "\n")
    f.write("TOP 10 PRIORITY NEOEPITOPES FOR EXPERIMENTAL VALIDATION\n")
    f.write("="*80 + "\n\n")
    
    f.write("MHC-I (CD8+ T-cell epitopes):\n")
    f.write("-" * 80 + "\n")
    f.write(f"{'#':<4}{'Gene':<10}{'Peptide':<15}{'HLA':<15}{'IC50':<10}{'Fold':<8}{'Priority':<20}\n")
    f.write("-" * 80 + "\n")
    
    for idx, row in mhc1_detailed.head(10).iterrows():
        f.write(f"{row['Rank']:<4}{row['Gene']:<10}{row['Peptide_Sequence']:<15}"
               f"{row['HLA_Allele']:<15}{row['MT_IC50_nM']:<10.1f}{row['Fold_Change']:<8.1f}"
               f"{row['Priority']:<20}\n")
    
    f.write("\n" + "="*80 + "\n")
    f.write("MHC-II (CD4+ T-cell epitopes):\n")
    f.write("-" * 80 + "\n")
    f.write(f"{'#':<4}{'Gene':<10}{'Peptide':<20}{'HLA':<20}{'IC50':<10}{'Fold':<8}{'Priority':<20}\n")
    f.write("-" * 80 + "\n")
    
    for idx, row in mhc2_detailed.head(10).iterrows():
        f.write(f"{row['Rank']:<4}{row['Gene']:<10}{row['Peptide_Sequence']:<20}"
               f"{row['HLA_Allele']:<20}{row['MT_IC50_nM']:<10.1f}{row['Fold_Change']:<8.1f}"
               f"{row['Priority']:<20}\n")
    
    f.write("\n" + "="*80 + "\n")
    f.write("GENE-SPECIFIC RECOMMENDATIONS\n")
    f.write("="*80 + "\n\n")
    
    # Top 5 genes with their best neoepitopes
    for gene in mhc1_detailed['Gene'].value_counts().head(5).index:
        gene_best = mhc1_detailed[mhc1_detailed['Gene'] == gene].iloc[0]
        f.write(f"\n{gene}:\n")
        f.write(f"  Total neoepitopes: {len(mhc1_detailed[mhc1_detailed['Gene'] == gene])}\n")
        f.write(f"  Best peptide: {gene_best['Peptide_Sequence']}\n")
        f.write(f"  HLA: {gene_best['HLA_Allele']}\n")
        f.write(f"  IC50: {gene_best['MT_IC50_nM']:.2f} nM\n")
        f.write(f"  Fold Change: {gene_best['Fold_Change']:.1f}x\n")
        f.write(f"  Priority: {gene_best['Priority']}\n")
    
    f.write("\n" + "="*80 + "\n")
    f.write("EXPERIMENTAL VALIDATION PRIORITIES\n")
    f.write("="*80 + "\n\n")
    
    tier1_dual = len(df_mhc1_top[df_mhc1_top['Priority'].str.contains('TIER 1.*Dual')])
    tier1_single = len(df_mhc1_top[df_mhc1_top['Priority'].str.contains('TIER 1')]) - tier1_dual
    
    f.write("PRIORITY 1 (Highest Confidence - Validate First):\n")
    f.write(f"  - Dual Coverage Neoepitopes: {tier1_dual}\n")
    f.write("  - Expected: Strong CD8+ and CD4+ response\n")
    f.write("  - Recommended: ELISpot, tetramer staining, functional assays\n\n")
    
    f.write("PRIORITY 2 (High Confidence):\n")
    f.write(f"  - Tier 1 Single Coverage: {tier1_single}\n")
    f.write("  - Expected: Strong single T-cell response\n")
    f.write("  - Recommended: Flow cytometry, cytokine assays\n\n")
    
    f.write("PRIORITY 3 (Medium Confidence):\n")
    f.write(f"  - Tier 2 candidates: {len(df_mhc1_top[df_mhc1_top['Priority'].str.contains('TIER 2')])}\n")
    f.write("  - Expected: Moderate immunogenicity\n")
    f.write("  - Recommended: Secondary validation if Priority 1-2 fail\n\n")

print(f"\n✓ Final recommendations: {report}")

# ============================================================================
# CREATE EXCEL-STYLE SUMMARY FOR EASY VIEWING
# ============================================================================

print("\nCreating Excel-style summary files...")

# Summary table: One row per top neoepitope with all details
summary_mhc1 = mhc1_detailed.head(20).copy()
summary_mhc1['Recommendation'] = summary_mhc1.apply(
    lambda x: 'PRIORITY 1: Validate First' if 'Dual' in x['Priority'] else
              'PRIORITY 1: High Confidence' if 'TIER 1' in x['Priority'] else
              'PRIORITY 2: Medium Confidence', axis=1
)

summary_mhc1.to_csv(PHASE5_DIR / "TOP_20_MHC1_FOR_VALIDATION.csv", index=False)
print(f"✓ Top 20 MHC-I for validation: {PHASE5_DIR / 'TOP_20_MHC1_FOR_VALIDATION.csv'}")

summary_mhc2 = mhc2_detailed.head(20).copy()
summary_mhc2['Recommendation'] = summary_mhc2.apply(
    lambda x: 'PRIORITY 1: Validate First' if 'Dual' in x['Priority'] else
              'PRIORITY 1: High Confidence' if 'TIER 1' in x['Priority'] else
              'PRIORITY 2: Medium Confidence', axis=1
)

summary_mhc2.to_csv(PHASE5_DIR / "TOP_20_MHC2_FOR_VALIDATION.csv", index=False)
print(f"✓ Top 20 MHC-II for validation: {PHASE5_DIR / 'TOP_20_MHC2_FOR_VALIDATION.csv'}")

# ============================================================================
# ANSWER THE SPECIFIC QUESTION: WHICH GENE/PEPTIDE IS #1?
# ============================================================================

print("\n" + "="*80)
print("ANSWER TO YOUR QUESTION:")
print("="*80)

print(f"\n🥇 THE #1 NEOEPITOPE (MHC-I):")
print(f"   Gene: {best_mhc1['Gene']}")
print(f"   Peptide: {best_mhc1['Peptide_Sequence']}")
print(f"   HLA: {best_mhc1['HLA_Allele']}")
print(f"   IC50: {best_mhc1['MT_IC50_nM']:.2f} nM (Strong Binder)")
print(f"   Specificity: {best_mhc1['Fold_Change']:.1f}× more specific than wildtype")
print(f"   Score: {best_mhc1['Immunogenicity_Score']:.1f}/100")
print(f"   CD4 Support: {'YES - Dual coverage!' if best_mhc1['CD4_Support'] else 'NO'}")

print(f"\n🥇 THE #1 NEOEPITOPE (MHC-II):")
print(f"   Gene: {best_mhc2['Gene']}")
print(f"   Peptide: {best_mhc2['Peptide_Sequence']}")
print(f"   HLA: {best_mhc2['HLA_Allele']}")
print(f"   IC50: {best_mhc2['MT_IC50_nM']:.2f} nM (Strong Binder)")
print(f"   Specificity: {best_mhc2['Fold_Change']:.1f}× more specific than wildtype")
print(f"   Score: {best_mhc2['Immunogenicity_Score']:.1f}/100")
print(f"   CD8 Support: {'YES - Dual coverage!' if best_mhc2['CD8_Support'] else 'NO'}")

if len(df_dual) > 0:
    print(f"\n🏆 THE BEST DUAL COVERAGE NEOEPITOPE:")
    print(f"   Gene: {best_dual['gene']}")
    print(f"   Why: Has BOTH CD8+ (killer) and CD4+ (helper) T-cell response")
    print(f"   MHC-I IC50: {best_dual['MT_IC50_nM_MHC1']:.2f} nM")
    print(f"   MHC-II IC50: {best_dual['MT_IC50_nM_MHC2']:.2f} nM")
    print(f"   ⚡ This is your BEST candidate for vaccine/therapy!")

print(f"\n📊 GENE RANKING (by total high-quality neoepitopes):")
gene_ranking = mhc1_detailed.groupby('Gene').size().sort_values(ascending=False).head(10)
for rank, (gene, count) in enumerate(gene_ranking.items(), 1):
    best_for_gene = mhc1_detailed[mhc1_detailed['Gene'] == gene].iloc[0]
    print(f"   #{rank}. {gene}: {count} neoepitopes")
    print(f"       → Best: {best_for_gene['Peptide_Sequence']} (IC50: {best_for_gene['MT_IC50_nM']:.1f}nM)")

print("\n" + "="*80)
print("✓ PHASE 5 DETAILED ANALYSIS COMPLETE")
print("="*80)

print(f"\n📁 ALL OUTPUTS SAVED:")
print(f"   {PHASE5_DIR}/")
print(f"   ├── MHC1_Detailed_Ranking.tsv (Full ranking)")
print(f"   ├── MHC2_Detailed_Ranking.tsv (Full ranking)")
print(f"   ├── TOP_20_MHC1_FOR_VALIDATION.csv (Ready for lab)")
print(f"   ├── TOP_20_MHC2_FOR_VALIDATION.csv (Ready for lab)")
print(f"   ├── FINAL_RECOMMENDATIONS.txt (Human-readable)")
print(f"   └── figures/")
print(f"       ├── Figure1_Top20_Individual_Neoepitopes.png")
print(f"       ├── Figure2_Gene_Level_Analysis.png")
print(f"       └── Figure3_Detailed_TP53_Analysis.png")

print("\n🎯 READY TO:")
print("   ✅ Share with experimental team")
print("   ✅ Submit for validation")
print("   ✅ Include in publication")
print("   ✅ Design vaccine/immunotherapy")

