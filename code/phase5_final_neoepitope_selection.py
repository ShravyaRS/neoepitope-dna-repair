"""
===============================================================================
PHASE 5: FINAL NEOEPITOPE SELECTION & VALIDATION
Industry-Standard Filtering - pVACseq Equivalent
===============================================================================
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['figure.dpi'] = 300

# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    BASE_DIR = Path("/mnt/c/Users/DELL/neo_project/data/output/Phase4")
    PHASE5_DIR = Path("/mnt/c/Users/DELL/neo_project/data/output/Phase5")
    
    # Input files
    MHC1_CANDIDATES = BASE_DIR / "candidates/top_1000_mhc1_candidates.tsv"
    MHC2_CANDIDATES = BASE_DIR / "mhc2/mhc2_candidates.tsv"
    
    # STRICT NEOEPITOPE CRITERIA (pVACseq-equivalent)
    MHC1_IC50_THRESHOLD = 500      # nM
    MHC1_STRONG_THRESHOLD = 50     # nM
    MHC2_IC50_THRESHOLD = 1000     # nM
    MHC2_STRONG_THRESHOLD = 100    # nM
    
    MIN_FOLD_CHANGE = 5            # WT/MT ratio
    MIN_WT_IC50 = 1000             # WT should NOT bind strongly
    
    # Immunogenicity scoring weights
    WEIGHT_BINDING = 0.35
    WEIGHT_SPECIFICITY = 0.35
    WEIGHT_COVERAGE = 0.30
    
    PHASE5_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("PHASE 5: FINAL NEOEPITOPE SELECTION")
print("Applying pVACseq-Equivalent Strict Criteria")
print("="*80)

# ============================================================================
# LOAD DATA
# ============================================================================

print("\nLoading Phase 4 candidates...")
df_mhc1 = pd.read_csv(Config.MHC1_CANDIDATES, sep='\t')
df_mhc2 = pd.read_csv(Config.MHC2_CANDIDATES, sep='\t')

print(f"✓ MHC-I candidates: {len(df_mhc1):,}")
print(f"✓ MHC-II candidates: {len(df_mhc2):,}")

# ============================================================================
# STRICT FILTERING - MHC-I
# ============================================================================

print("\n" + "="*80)
print("APPLYING STRICT FILTERS - MHC-I")
print("="*80)

print(f"\nStarting with: {len(df_mhc1):,} MHC-I candidates")

# Filter 1: Strong MT binding
df_mhc1_strict = df_mhc1[df_mhc1['MT_IC50_nM'] < Config.MHC1_IC50_THRESHOLD].copy()
print(f"[Filter 1] MT IC50 < {Config.MHC1_IC50_THRESHOLD}nM: {len(df_mhc1_strict):,}")

# Filter 2: Weak WT binding (tumor specificity)
df_mhc1_strict = df_mhc1_strict[df_mhc1_strict['WT_IC50_nM'] > Config.MIN_WT_IC50]
print(f"[Filter 2] WT IC50 > {Config.MIN_WT_IC50}nM: {len(df_mhc1_strict):,}")

# Filter 3: High fold change
df_mhc1_strict = df_mhc1_strict[df_mhc1_strict['Fold_Change'] >= Config.MIN_FOLD_CHANGE]
print(f"[Filter 3] Fold change ≥ {Config.MIN_FOLD_CHANGE}x: {len(df_mhc1_strict):,}")

# ============================================================================
# STRICT FILTERING - MHC-II
# ============================================================================

print("\n" + "="*80)
print("APPLYING STRICT FILTERS - MHC-II")
print("="*80)

print(f"\nStarting with: {len(df_mhc2):,} MHC-II candidates")

# Filter 1: Strong MT binding
df_mhc2_strict = df_mhc2[df_mhc2['MT_IC50_nM'] < Config.MHC2_IC50_THRESHOLD].copy()
print(f"[Filter 1] MT IC50 < {Config.MHC2_IC50_THRESHOLD}nM: {len(df_mhc2_strict):,}")

# Filter 2: Weak WT binding
df_mhc2_strict = df_mhc2_strict[df_mhc2_strict['WT_IC50_nM'] > Config.MIN_WT_IC50]
print(f"[Filter 2] WT IC50 > {Config.MIN_WT_IC50}nM: {len(df_mhc2_strict):,}")

# Filter 3: High fold change
df_mhc2_strict = df_mhc2_strict[df_mhc2_strict['Fold_Change'] >= Config.MIN_FOLD_CHANGE]
print(f"[Filter 3] Fold change ≥ {Config.MIN_FOLD_CHANGE}x: {len(df_mhc2_strict):,}")

# ============================================================================
# IMMUNOGENICITY SCORING
# ============================================================================

print("\n" + "="*80)
print("CALCULATING IMMUNOGENICITY SCORES")
print("="*80)

def calculate_immunogenicity_mhc1(row):
    """Calculate comprehensive immunogenicity score"""
    score = 0
    
    # Binding affinity (35%)
    if row['MT_IC50_nM'] < 50:
        binding_score = 35
    elif row['MT_IC50_nM'] < 100:
        binding_score = 30
    elif row['MT_IC50_nM'] < 500:
        binding_score = 20
    else:
        binding_score = 0
    
    # Tumor specificity (35%)
    if row['Fold_Change'] >= 20:
        specificity_score = 35
    elif row['Fold_Change'] >= 10:
        specificity_score = 30
    elif row['Fold_Change'] >= 5:
        specificity_score = 20
    else:
        specificity_score = 0
    
    # Quality score component (30%)
    quality_score = min(row['Quality_Score'] / 100 * 30, 30)
    
    return binding_score + specificity_score + quality_score

def calculate_immunogenicity_mhc2(row):
    """Calculate MHC-II immunogenicity score"""
    score = 0
    
    # Binding affinity (35%)
    if row['MT_IC50_nM'] < 100:
        binding_score = 35
    elif row['MT_IC50_nM'] < 500:
        binding_score = 25
    elif row['MT_IC50_nM'] < 1000:
        binding_score = 15
    else:
        binding_score = 0
    
    # Tumor specificity (35%)
    if row['Fold_Change'] >= 20:
        specificity_score = 35
    elif row['Fold_Change'] >= 10:
        specificity_score = 30
    elif row['Fold_Change'] >= 5:
        specificity_score = 20
    else:
        specificity_score = 0
    
    # Quality score component (30%)
    quality_score = min(row['Quality_Score'] / 100 * 30, 30)
    
    return binding_score + specificity_score + quality_score

df_mhc1_strict['Immunogenicity_Score'] = df_mhc1_strict.apply(calculate_immunogenicity_mhc1, axis=1)
df_mhc2_strict['Immunogenicity_Score'] = df_mhc2_strict.apply(calculate_immunogenicity_mhc2, axis=1)

# ============================================================================
# MUTATION-LEVEL AGGREGATION
# ============================================================================

print("\nAggregating by mutation...")

# MHC-I: Group by mutation and take best peptide
mhc1_by_mutation = df_mhc1_strict.sort_values('Immunogenicity_Score', ascending=False).groupby('mutation_id').first().reset_index()
print(f"✓ MHC-I: {len(mhc1_by_mutation)} unique mutations")

# MHC-II: Group by mutation and take best peptide
mhc2_by_mutation = df_mhc2_strict.sort_values('Immunogenicity_Score', ascending=False).groupby('mutation_id').first().reset_index()
print(f"✓ MHC-II: {len(mhc2_by_mutation)} unique mutations")

# Find mutations with BOTH MHC-I and MHC-II response
common_mutations = set(mhc1_by_mutation['mutation_id']) & set(mhc2_by_mutation['mutation_id'])
print(f"✓ Mutations with BOTH MHC-I and MHC-II: {len(common_mutations)}")

# ============================================================================
# FINAL NEOEPITOPE RANKING
# ============================================================================

print("\n" + "="*80)
print("FINAL NEOEPITOPE RANKING")
print("="*80)

# Add coverage bonus for mutations with both MHC-I and MHC-II
mhc1_by_mutation['Has_MHC2'] = mhc1_by_mutation['mutation_id'].isin(common_mutations)
mhc2_by_mutation['Has_MHC1'] = mhc2_by_mutation['mutation_id'].isin(common_mutations)

# Boost score for dual coverage
mhc1_by_mutation['Final_Score'] = mhc1_by_mutation['Immunogenicity_Score'] + \
                                   mhc1_by_mutation['Has_MHC2'].astype(int) * 15

mhc2_by_mutation['Final_Score'] = mhc2_by_mutation['Immunogenicity_Score'] + \
                                   mhc2_by_mutation['Has_MHC1'].astype(int) * 15

# Sort by final score
mhc1_final = mhc1_by_mutation.sort_values('Final_Score', ascending=False)
mhc2_final = mhc2_by_mutation.sort_values('Final_Score', ascending=False)

# ============================================================================
# SAVE FINAL NEOEPITOPES
# ============================================================================

print("\nSaving final neoepitopes...")

# Top 50 MHC-I neoepitopes
mhc1_final.head(50).to_csv(Config.PHASE5_DIR / "top_50_mhc1_neoepitopes.tsv", sep='\t', index=False)
print(f"✓ Top 50 MHC-I: {Config.PHASE5_DIR / 'top_50_mhc1_neoepitopes.tsv'}")

# Top 50 MHC-II neoepitopes
mhc2_final.head(50).to_csv(Config.PHASE5_DIR / "top_50_mhc2_neoepitopes.tsv", sep='\t', index=False)
print(f"✓ Top 50 MHC-II: {Config.PHASE5_DIR / 'top_50_mhc2_neoepitopes.tsv'}")

# Combined top 20 (best from each)
combined_top = pd.concat([
    mhc1_final.head(20).assign(Source='MHC-I'),
    mhc2_final.head(20).assign(Source='MHC-II')
]).sort_values('Final_Score', ascending=False)

combined_top.to_csv(Config.PHASE5_DIR / "top_20_combined_neoepitopes.tsv", sep='\t', index=False)
print(f"✓ Top 20 Combined: {Config.PHASE5_DIR / 'top_20_combined_neoepitopes.tsv'}")

# Mutations with dual coverage (HIGH PRIORITY)
dual_coverage_mhc1 = mhc1_final[mhc1_final['Has_MHC2']].head(20)
dual_coverage_mhc2 = mhc2_final[mhc2_final['Has_MHC1']].head(20)

dual_coverage = pd.merge(
    dual_coverage_mhc1[['mutation_id', 'gene', 'mt_peptide', 'HLA', 'MT_IC50_nM', 'Fold_Change', 'Final_Score']],
    dual_coverage_mhc2[['mutation_id', 'gene', 'mt_peptide', 'HLA', 'MT_IC50_nM', 'Fold_Change', 'Final_Score']],
    on=['mutation_id', 'gene'],
    suffixes=('_MHC1', '_MHC2')
)

dual_coverage.to_csv(Config.PHASE5_DIR / "dual_coverage_neoepitopes.tsv", sep='\t', index=False)
print(f"✓ Dual Coverage (CD8+CD4): {Config.PHASE5_DIR / 'dual_coverage_neoepitopes.tsv'}")

# ============================================================================
# SUMMARY REPORT
# ============================================================================

print("\n" + "="*80)
print("FINAL NEOEPITOPE SUMMARY")
print("="*80)

print(f"\nMHC-I Neoepitopes:")
print(f"  Total passing strict filters: {len(mhc1_final)}")
print(f"  Strong binders (IC50<50nM): {(mhc1_final['MT_IC50_nM']<50).sum()}")
print(f"  With MHC-II support: {mhc1_final['Has_MHC2'].sum()}")

print(f"\nMHC-II Neoepitopes:")
print(f"  Total passing strict filters: {len(mhc2_final)}")
print(f"  Strong binders (IC50<100nM): {(mhc2_final['MT_IC50_nM']<100).sum()}")
print(f"  With MHC-I support: {mhc2_final['Has_MHC1'].sum()}")

print(f"\nDual Coverage Neoepitopes (HIGHEST PRIORITY):")
print(f"  Mutations with both CD8+ and CD4+ response: {len(common_mutations)}")

print(f"\nTop 10 Genes (MHC-I):")
for gene, count in mhc1_final['gene'].value_counts().head(10).items():
    print(f"  {gene}: {count}")

print(f"\nTop 10 Genes (MHC-II):")
for gene, count in mhc2_final['gene'].value_counts().head(10).items():
    print(f"  {gene}: {count}")

# Generate report
report = Config.PHASE5_DIR / "Phase5_Final_Report.txt"
with open(report, 'w') as f:
    f.write("="*80 + "\n")
    f.write("PHASE 5: FINAL NEOEPITOPE SELECTION REPORT\n")
    f.write("="*80 + "\n\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    f.write("FILTERING SUMMARY:\n")
    f.write(f"MHC-I: {len(df_mhc1):,} → {len(mhc1_final)} final neoepitopes\n")
    f.write(f"MHC-II: {len(df_mhc2):,} → {len(mhc2_final)} final neoepitopes\n")
    f.write(f"Dual Coverage: {len(common_mutations)} mutations\n\n")
    
    f.write("TOP 20 NEOEPITOPES (MHC-I):\n")
    f.write(f"{'Gene':<10}{'Peptide':<15}{'HLA':<15}{'MT IC50':<12}{'Fold':<10}{'Score':<10}\n")
    f.write("-"*80 + "\n")
    for idx, row in mhc1_final.head(20).iterrows():
        f.write(f"{row['gene']:<10}{row['mt_peptide']:<15}{row['HLA']:<15}"
               f"{row['MT_IC50_nM']:<12.2f}{row['Fold_Change']:<10.1f}{row['Final_Score']:<10.1f}\n")

print(f"\n✓ Report: {report}")

print("\n" + "="*80)
print("✓ PHASE 5 COMPLETE - FINAL NEOEPITOPES IDENTIFIED")
print("="*80)
print(f"\nOutputs: {Config.PHASE5_DIR}")
print(f"\nREADY FOR EXPERIMENTAL VALIDATION!")

