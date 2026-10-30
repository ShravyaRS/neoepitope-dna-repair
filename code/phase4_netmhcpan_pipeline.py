"""
===============================================================================
PHASE 4: MHC-I BINDING PREDICTION & IMMUNOGENICITY ASSESSMENT
Industry Standard Pipeline using NetMHCpan 4.1
===============================================================================

This pipeline implements the gold-standard neoepitope prediction workflow:
1. NetMHCpan 4.1 binding predictions (WT and MT peptides)
2. Differential Agretopicity Index (DAI) calculation
3. Strong binder filtering (IC50 < 500 nM)
4. Immunogenicity scoring
5. Candidate prioritization for experimental validation

References:
- Reynisson et al. NetMHCpan-4.1, Nucleic Acids Res. 2020
- Łuksza et al. A neoantigen fitness model, Nature 2017
- Richman et al. Neoantigen dissimilarity, Cell 2019

Author: Neoepitope Analysis Team
Date: November 2025
===============================================================================
"""

import pandas as pd
import numpy as np
import subprocess
import os
from pathlib import Path
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set publication quality
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 11

# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    """Phase 4 Configuration"""
    
    # Paths
    BASE_DIR = Path("/mnt/c/Users/DELL/neo_project")
    PHASE3_DIR = BASE_DIR / "data/output/Phase3"
    PHASE4_DIR = BASE_DIR / "data/output/Phase4"
    
    # Input files from Phase 3
    MUTANT_FASTA = PHASE3_DIR / "peptides/mutant_peptides.fasta"
    WILDTYPE_FASTA = PHASE3_DIR / "peptides/wildtype_peptides.fasta"
    PEPTIDE_TABLE = PHASE3_DIR / "peptides/peptide_table.tsv"
    MUTATION_SUMMARY = PHASE3_DIR / "analysis_data/mutation_summary.tsv"
    
    # Output directories
    PREDICTIONS_DIR = PHASE4_DIR / "predictions"
    FIGURES_DIR = PHASE4_DIR / "figures"
    CANDIDATES_DIR = PHASE4_DIR / "candidate_neoepitopes"
    ANALYSIS_DIR = PHASE4_DIR / "analysis"
    LOGS_DIR = PHASE4_DIR / "logs"
    
    # HLA Alleles (Common alleles for population coverage)
    # Modify based on your cohort or use patient-specific HLA types
    HLA_ALLELES = [
        "HLA-A*02:01",  # 47% Caucasian, 39% Asian - Most common
        "HLA-A*01:01",  # 25% Caucasian
        "HLA-A*03:01",  # 20% Caucasian
        "HLA-A*24:02",  # 40% Asian, 20% Caucasian
        "HLA-B*07:02",  # 25% Caucasian
        "HLA-B*44:02",  # 20% Caucasian
        "HLA-B*08:01",  # 20% Caucasian
        "HLA-C*07:01",  # 25% population
    ]
    
    # Thresholds (Industry Standard)
    STRONG_BINDER_THRESHOLD = 50    # nM (SB: IC50 < 50 nM)
    WEAK_BINDER_THRESHOLD = 500     # nM (WB: IC50 < 500 nM)
    MIN_FOLD_CHANGE = 2.0           # MT should bind 2x better than WT
    MIN_DAI_SCORE = 0.5             # Differential Agretopicity Index
    
    # NetMHCpan settings
    NETMHCPAN_PATH = "netMHCpan"  # Assumes netMHCpan is in PATH
    # If not in PATH, specify full path:
    # NETMHCPAN_PATH = "/path/to/netMHCpan-4.1/netMHCpan"
    
    @classmethod
    def create_directories(cls):
        """Create all output directories"""
        for dir_path in [cls.PHASE4_DIR, cls.PREDICTIONS_DIR, cls.FIGURES_DIR,
                         cls.CANDIDATES_DIR, cls.ANALYSIS_DIR, cls.LOGS_DIR]:
            dir_path.mkdir(parents=True, exist_ok=True)

# ============================================================================
# LOGGER
# ============================================================================

class Logger:
    """Simple logger for pipeline tracking"""
    
    def __init__(self, log_file):
        self.log_file = log_file
        self.start_time = datetime.now()
        
    def log(self, message, level="INFO"):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_entry = f"[{timestamp}] [{level}] {message}"
        print(log_entry)
        with open(self.log_file, 'a', encoding='utf-8') as f:
            f.write(log_entry + "\n")
    
    def section(self, title):
        separator = "=" * 80
        self.log(f"\n{separator}")
        self.log(title)
        self.log(separator)
    
    def subsection(self, title):
        self.log(f"\n--- {title} ---")
    
    def success(self, message):
        self.log(message, level="SUCCESS")
    
    def error(self, message):
        self.log(message, level="ERROR")
    
    def warning(self, message):
        self.log(message, level="WARNING")
    
    def get_runtime(self):
        elapsed = (datetime.now() - self.start_time).total_seconds()
        return f"{elapsed:.2f} seconds"

# ============================================================================
# STEP 1: CHECK NETMHCPAN INSTALLATION
# ============================================================================

def check_netmhcpan_installation(logger):
    """
    Check if NetMHCpan is installed and accessible
    """
    logger.section("STEP 1: CHECKING NETMHCPAN INSTALLATION")
    
    try:
        # Try to run NetMHCpan help
        result = subprocess.run(
            [Config.NETMHCPAN_PATH, "-h"],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if result.returncode == 0:
            logger.success("✓ NetMHCpan is installed and accessible")
            
            # Extract version if possible
            if "NetMHCpan" in result.stdout:
                version_line = [line for line in result.stdout.split('\n') 
                               if 'NetMHCpan' in line][0]
                logger.log(f"Version: {version_line.strip()}")
            
            return True
        else:
            logger.error("✗ NetMHCpan found but returned error")
            return False
            
    except FileNotFoundError:
        logger.error("✗ NetMHCpan not found in PATH")
        logger.error("Please install NetMHCpan 4.1 from:")
        logger.error("https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/")
        logger.error("\nOr specify the full path in Config.NETMHCPAN_PATH")
        return False
    except Exception as e:
        logger.error(f"✗ Error checking NetMHCpan: {str(e)}")
        return False

# ============================================================================
# STEP 2: RUN NETMHCPAN PREDICTIONS
# ============================================================================

def run_netmhcpan_prediction(fasta_file, output_file, peptide_type, logger):
    """
    Run NetMHCpan prediction for given FASTA file
    
    Args:
        fasta_file: Path to peptide FASTA file
        output_file: Path to output file
        peptide_type: "Mutant" or "Wildtype"
    """
    logger.subsection(f"Running NetMHCpan for {peptide_type} Peptides")
    
    # Format HLA alleles for NetMHCpan
    hla_string = ",".join(Config.HLA_ALLELES)
    
    # NetMHCpan command
    cmd = [
        Config.NETMHCPAN_PATH,
        "-p", str(fasta_file),
        "-a", hla_string,
        "-BA",  # Use Binding Affinity mode
        "-xls",  # XLS output format
        "-xlsfile", str(output_file)
    ]
    
    logger.log(f"Command: {' '.join(cmd)}")
    logger.log("This may take 5-15 minutes depending on peptide count...")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800  # 30 minute timeout
        )
        
        if result.returncode == 0:
            logger.success(f"✓ NetMHCpan completed for {peptide_type} peptides")
            
            # Check output file
            if Path(output_file).exists():
                file_size = Path(output_file).stat().st_size / (1024 * 1024)  # MB
                logger.log(f"  Output file size: {file_size:.2f} MB")
                return True
            else:
                logger.error(f"✗ Output file not created: {output_file}")
                return False
        else:
            logger.error(f"✗ NetMHCpan failed with return code: {result.returncode}")
            logger.error(f"Error output: {result.stderr[:500]}")
            return False
            
    except subprocess.TimeoutExpired:
        logger.error("✗ NetMHCpan timed out (>30 minutes)")
        return False
    except Exception as e:
        logger.error(f"✗ Error running NetMHCpan: {str(e)}")
        return False

def parse_netmhcpan_output(xls_file, logger):
    """
    Parse NetMHCpan XLS output file
    
    Returns:
        DataFrame with predictions
    """
    logger.log(f"Parsing NetMHCpan output: {xls_file}")
    
    try:
        # NetMHCpan XLS format has header starting with specific columns
        df = pd.read_csv(xls_file, sep='\t', skiprows=1)
        
        # Standardize column names
        df.columns = df.columns.str.strip()
        
        # Key columns: Pos, HLA, Peptide, Core, Of, Gp, Gl, Ip, Il, Icore, 
        #              Identity, Score_EL, %Rank_EL, Score_BA, %Rank_BA, Aff(nM), BindLevel
        
        logger.success(f"✓ Parsed {len(df):,} predictions")
        logger.log(f"  Columns: {', '.join(df.columns.tolist()[:10])}...")
        
        return df
        
    except Exception as e:
        logger.error(f"✗ Error parsing NetMHCpan output: {str(e)}")
        return None

# ============================================================================
# STEP 3: INTEGRATE PREDICTIONS WITH PHASE 3 DATA
# ============================================================================

def integrate_predictions(mt_predictions, wt_predictions, peptide_table, logger):
    """
    Integrate MT and WT predictions with Phase 3 peptide data
    """
    logger.section("STEP 3: INTEGRATING PREDICTIONS")
    
    # Load peptide table
    logger.log("Loading Phase 3 peptide table...")
    df_peptides = pd.read_csv(peptide_table, sep='\t')
    logger.log(f"  Loaded {len(df_peptides):,} peptides")
    
    # Process MT predictions
    logger.subsection("Processing Mutant Predictions")
    mt_predictions['peptide_type'] = 'MT'
    mt_predictions = mt_predictions.rename(columns={'Peptide': 'mt_peptide'})
    
    # Process WT predictions
    logger.subsection("Processing Wildtype Predictions")
    wt_predictions['peptide_type'] = 'WT'
    wt_predictions = wt_predictions.rename(columns={'Peptide': 'wt_peptide'})
    
    # Merge with peptide table
    logger.subsection("Merging with Phase 3 Data")
    
    # Merge MT predictions
    df_merged = df_peptides.merge(
        mt_predictions[['mt_peptide', 'HLA', 'Aff(nM)', '%Rank_BA', 'BindLevel']],
        on='mt_peptide',
        how='left',
        suffixes=('', '_MT')
    )
    
    # Rename MT columns
    df_merged = df_merged.rename(columns={
        'Aff(nM)': 'MT_Affinity_nM',
        '%Rank_BA': 'MT_Rank',
        'BindLevel': 'MT_BindLevel'
    })
    
    # Merge WT predictions
    df_merged = df_merged.merge(
        wt_predictions[['wt_peptide', 'HLA', 'Aff(nM)', '%Rank_BA', 'BindLevel']],
        left_on=['wt_peptide', 'HLA'],
        right_on=['wt_peptide', 'HLA'],
        how='left',
        suffixes=('', '_WT')
    )
    
    # Rename WT columns
    df_merged = df_merged.rename(columns={
        'Aff(nM)': 'WT_Affinity_nM',
        '%Rank_BA': 'WT_Rank',
        'BindLevel': 'WT_BindLevel'
    })
    
    logger.success(f"✓ Integrated predictions for {len(df_merged):,} peptide-HLA pairs")
    
    return df_merged

# ============================================================================
# STEP 4: CALCULATE IMMUNOGENICITY SCORES
# ============================================================================

def calculate_immunogenicity_scores(df, logger):
    """
    Calculate immunogenicity metrics:
    1. Differential Agretopicity Index (DAI)
    2. Binding fold change
    3. Neoantigen quality score
    """
    logger.section("STEP 4: CALCULATING IMMUNOGENICITY SCORES")
    
    # 1. Fold change (WT_Affinity / MT_Affinity)
    # Higher fold change = better (MT binds stronger than WT)
    logger.subsection("Calculating Binding Fold Change")
    df['Fold_Change'] = df['WT_Affinity_nM'] / df['MT_Affinity_nM']
    df['Fold_Change'] = df['Fold_Change'].fillna(0)
    
    # 2. Differential Agretopicity Index (DAI)
    # DAI = 1 - (MT_Rank / 100) if MT binds and WT doesn't bind strongly
    logger.subsection("Calculating DAI Score")
    
    def calculate_dai(row):
        mt_rank = row['MT_Rank']
        wt_rank = row['WT_Rank']
        
        if pd.isna(mt_rank):
            return 0
        
        # Strong MT binding
        if mt_rank < 2.0:  # Top 2% rank = strong binder
            if pd.isna(wt_rank) or wt_rank > 2.0:
                # WT doesn't bind strongly
                return 1 - (mt_rank / 100)
            else:
                # WT also binds - calculate difference
                return max(0, (wt_rank - mt_rank) / 100)
        else:
            return 0
    
    df['DAI_Score'] = df.apply(calculate_dai, axis=1)
    
    # 3. Neoantigen Quality Score (composite)
    # Combines: binding affinity, fold change, DAI
    logger.subsection("Calculating Neoantigen Quality Score")
    
    def calculate_quality_score(row):
        score = 0
        
        # Component 1: MT binding strength (40% weight)
        if row['MT_Affinity_nM'] < 50:  # Strong binder
            score += 40
        elif row['MT_Affinity_nM'] < 500:  # Weak binder
            score += 20
        
        # Component 2: Fold change (30% weight)
        if row['Fold_Change'] >= 10:
            score += 30
        elif row['Fold_Change'] >= 5:
            score += 20
        elif row['Fold_Change'] >= 2:
            score += 10
        
        # Component 3: DAI (30% weight)
        score += row['DAI_Score'] * 30
        
        return score
    
    df['Quality_Score'] = df.apply(calculate_quality_score, axis=1)
    
    # Statistics
    logger.log(f"\nScore Statistics:")
    logger.log(f"  Fold Change - Mean: {df['Fold_Change'].mean():.2f}, "
              f"Max: {df['Fold_Change'].max():.2f}")
    logger.log(f"  DAI Score - Mean: {df['DAI_Score'].mean():.3f}, "
              f"Max: {df['DAI_Score'].max():.3f}")
    logger.log(f"  Quality Score - Mean: {df['Quality_Score'].mean():.2f}, "
              f"Max: {df['Quality_Score'].max():.2f}")
    
    return df

# ============================================================================
# STEP 5: FILTER AND PRIORITIZE CANDIDATES
# ============================================================================

def filter_candidate_neoepitopes(df, logger):
    """
    Filter for high-quality neoepitope candidates using industry-standard criteria
    """
    logger.section("STEP 5: FILTERING CANDIDATE NEOEPITOPES")
    
    initial_count = len(df)
    logger.log(f"Starting with {initial_count:,} peptide-HLA pairs")
    
    # Filter 1: MT must be a binder (IC50 < 500 nM)
    df_filtered = df[df['MT_Affinity_nM'] < Config.WEAK_BINDER_THRESHOLD].copy()
    logger.log(f"\n[Filter 1] MT binding < 500 nM: {len(df_filtered):,} "
              f"({len(df_filtered)/initial_count*100:.1f}%)")
    
    # Filter 2: WT should not bind (or bind weaker)
    df_filtered = df_filtered[
        (df_filtered['WT_Affinity_nM'] > Config.WEAK_BINDER_THRESHOLD) |
        (df_filtered['Fold_Change'] >= Config.MIN_FOLD_CHANGE)
    ]
    logger.log(f"[Filter 2] WT non-binder or 2x fold change: {len(df_filtered):,} "
              f"({len(df_filtered)/initial_count*100:.1f}%)")
    
    # Filter 3: DAI score threshold
    df_filtered = df_filtered[df_filtered['DAI_Score'] >= Config.MIN_DAI_SCORE]
    logger.log(f"[Filter 3] DAI Score ≥ {Config.MIN_DAI_SCORE}: {len(df_filtered):,} "
              f"({len(df_filtered)/initial_count*100:.1f}%)")
    
    # Categorize by binding strength
    df_filtered['Binder_Category'] = pd.cut(
        df_filtered['MT_Affinity_nM'],
        bins=[0, 50, 500, float('inf')],
        labels=['Strong Binder (SB)', 'Weak Binder (WB)', 'Non-Binder']
    )
    
    logger.subsection("Candidate Summary")
    binder_counts = df_filtered['Binder_Category'].value_counts()
    for category, count in binder_counts.items():
        logger.log(f"  {category}: {count:,}")
    
    # Sort by quality score
    df_filtered = df_filtered.sort_values('Quality_Score', ascending=False)
    
    logger.success(f"\n✓ Final candidates: {len(df_filtered):,} high-quality neoepitopes")
    
    return df_filtered

# ============================================================================
# STEP 6: GENERATE OUTPUTS
# ============================================================================

def generate_candidate_outputs(df_candidates, df_all, logger):
    """
    Generate output files for candidate neoepitopes
    """
    logger.section("STEP 6: GENERATING OUTPUT FILES")
    
    # Output 1: Top candidates (all filters passed)
    logger.subsection("Top Candidate Neoepitopes")
    top_candidates = df_candidates.head(500)
    
    output_file = Config.CANDIDATES_DIR / "top_500_neoepitope_candidates.tsv"
    top_candidates.to_csv(output_file, sep='\t', index=False)
    logger.success(f"✓ Saved: {output_file}")
    
    # Output 2: Strong binders only
    strong_binders = df_candidates[df_candidates['MT_Affinity_nM'] < Config.STRONG_BINDER_THRESHOLD]
    output_file = Config.CANDIDATES_DIR / "strong_binders_IC50_below_50nM.tsv"
    strong_binders.to_csv(output_file, sep='\t', index=False)
    logger.log(f"✓ Strong binders (IC50<50nM): {len(strong_binders):,}")
    
    # Output 3: Gene-specific candidates
    logger.subsection("Gene-Specific Candidate Files")
    top_genes = df_candidates['gene'].value_counts().head(10).index
    
    for gene in top_genes:
        gene_candidates = df_candidates[df_candidates['gene'] == gene]
        output_file = Config.CANDIDATES_DIR / f"candidates_{gene}.tsv"
        gene_candidates.to_csv(output_file, sep='\t', index=False)
        logger.log(f"  {gene}: {len(gene_candidates):,} candidates")
    
    # Output 4: Summary statistics
    logger.subsection("Summary Statistics")
    
    summary_stats = {
        'Total_Peptide_HLA_Pairs': len(df_all),
        'Candidate_Neoepitopes': len(df_candidates),
        'Strong_Binders_SB': len(df_candidates[df_candidates['MT_Affinity_nM'] < 50]),
        'Weak_Binders_WB': len(df_candidates[(df_candidates['MT_Affinity_nM'] >= 50) & 
                                              (df_candidates['MT_Affinity_nM'] < 500)]),
        'Unique_Mutations_with_Candidates': df_candidates['mutation_id'].nunique(),
        'Unique_Genes_with_Candidates': df_candidates['gene'].nunique(),
        'Mean_MT_Affinity_nM': df_candidates['MT_Affinity_nM'].mean(),
        'Mean_Fold_Change': df_candidates['Fold_Change'].mean(),
        'Mean_DAI_Score': df_candidates['DAI_Score'].mean(),
        'Mean_Quality_Score': df_candidates['Quality_Score'].mean()
    }
    
    summary_df = pd.DataFrame([summary_stats])
    output_file = Config.ANALYSIS_DIR / "summary_statistics.tsv"
    summary_df.to_csv(output_file, sep='\t', index=False)
    logger.success(f"✓ Summary statistics: {output_file}")
    
    # Display key metrics
    for key, value in summary_stats.items():
        if isinstance(value, float):
            logger.log(f"  {key}: {value:.2f}")
        else:
            logger.log(f"  {key}: {value:,}")

# ============================================================================
# STEP 7: VISUALIZATIONS
# ============================================================================

def create_phase4_visualizations(df_candidates, df_all, logger):
    """
    Create publication-quality figures for Phase 4
    """
    logger.section("STEP 7: GENERATING VISUALIZATIONS")
    
    # Figure 1: Binding Affinity Distribution
    logger.log("\n[Figure 1] Binding affinity distributions...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Panel A: MT vs WT affinity scatter
    ax1 = axes[0, 0]
    
    # Sample for visualization (too many points)
    df_plot = df_all.sample(min(10000, len(df_all)))
    
    scatter = ax1.scatter(df_plot['WT_Affinity_nM'], df_plot['MT_Affinity_nM'],
                         c=df_plot['Quality_Score'], cmap='viridis',
                         alpha=0.6, s=20, edgecolors='none')
    
    ax1.axhline(y=500, color='red', linestyle='--', linewidth=2, label='IC50 = 500 nM')
    ax1.axvline(x=500, color='red', linestyle='--', linewidth=2)
    ax1.plot([0, 50000], [0, 50000], 'k--', linewidth=1, label='WT = MT')
    
    ax1.set_xlabel('WT Affinity (nM)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('MT Affinity (nM)', fontsize=12, fontweight='bold')
    ax1.set_title('A. MT vs WT Binding Affinity', fontsize=13, fontweight='bold')
    ax1.set_xlim(0, 5000)
    ax1.set_ylim(0, 5000)
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    cbar = plt.colorbar(scatter, ax=ax1)
    cbar.set_label('Quality Score', fontsize=10)
    
    # Panel B: Affinity histogram
    ax2 = axes[0, 1]
    
    ax2.hist(df_candidates['MT_Affinity_nM'], bins=50, color='#2E86AB',
             edgecolor='black', linewidth=1, alpha=0.7, label='Candidates')
    ax2.axvline(x=50, color='red', linestyle='--', linewidth=2, label='SB threshold (50 nM)')
    ax2.axvline(x=500, color='orange', linestyle='--', linewidth=2, label='WB threshold (500 nM)')
    
    ax2.set_xlabel('MT Binding Affinity (nM)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax2.set_title('B. Candidate Neoepitope Binding Distribution', fontsize=13, fontweight='bold')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # Panel C: Fold change distribution
    ax3 = axes[1, 0]
    
    fold_change_filtered = df_candidates[df_candidates['Fold_Change'] < 100]  # Remove outliers
    
    ax3.hist(fold_change_filtered['Fold_Change'], bins=50, color='#F18F01',
             edgecolor='black', linewidth=1, alpha=0.7)
    ax3.axvline(x=Config.MIN_FOLD_CHANGE, color='red', linestyle='--',
                linewidth=2, label=f'Threshold ({Config.MIN_FOLD_CHANGE}x)')
    
    ax3.set_xlabel('Binding Fold Change (WT/MT)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax3.set_title('C. Differential Binding (Fold Change)', fontsize=13, fontweight='bold')
    ax3.legend()
    ax3.grid(axis='y', alpha=0.3)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    # Panel D: Quality score distribution
    ax4 = axes[1, 1]
    
    ax4.hist(df_candidates['Quality_Score'], bins=50, color='#06A77D',
             edgecolor='black', linewidth=1, alpha=0.7)
    ax4.axvline(x=df_candidates['Quality_Score'].median(), color='red',
                linestyle='--', linewidth=2,
                label=f'Median: {df_candidates["Quality_Score"].median():.1f}')
    
    ax4.set_xlabel('Neoantigen Quality Score', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax4.set_title('D. Composite Quality Score Distribution', fontsize=13, fontweight='bold')
    ax4.legend()
    ax4.grid(axis='y', alpha=0.3)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    
    plt.tight_layout()
    output_file = Config.FIGURES_DIR / 'Figure1_Binding_Affinity_Analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.success(f"✓ Saved: {output_file}")
    
    # Figure 2: Gene and HLA Analysis
    logger.log("\n[Figure 2] Gene and HLA analysis...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Panel A: Top genes by candidate count
    ax1 = axes[0, 0]
    
    gene_counts = df_candidates['gene'].value_counts().head(15)
    y_pos = np.arange(len(gene_counts))
    colors_grad = plt.cm.YlOrRd(np.linspace(0.3, 0.9, len(gene_counts)))
    
    bars = ax1.barh(y_pos, gene_counts.values, color=colors_grad, edgecolor='black', linewidth=1.2)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(gene_counts.index, fontsize=11)
    ax1.set_xlabel('Number of Candidate Neoepitopes', fontsize=12, fontweight='bold')
    ax1.set_title('A. Top 15 Genes by Candidate Count', fontsize=13, fontweight='bold')
    ax1.invert_yaxis()
    ax1.grid(axis='x', alpha=0.3, linestyle='--')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    for bar in bars:
        width = bar.get_width()
        ax1.text(width, bar.get_y() + bar.get_height()/2.,
                f' {int(width):,}', va='center', fontsize=9, fontweight='bold')
    
    # Panel B: HLA allele distribution
    ax2 = axes[0, 1]
    
    hla_counts = df_candidates['HLA'].value_counts()
    colors_hla = plt.cm.tab10(np.arange(len(hla_counts)))
    
    bars = ax2.bar(range(len(hla_counts)), hla_counts.values, color=colors_hla,
                   edgecolor='black', linewidth=1.2)
    ax2.set_xticks(range(len(hla_counts)))
    ax2.set_xticklabels(hla_counts.index, rotation=45, ha='right', fontsize=10)
    ax2.set_ylabel('Number of Candidates', fontsize=12, fontweight='bold')
    ax2.set_title('B. Candidate Distribution by HLA Allele', fontsize=13, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    for bar in bars:
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height):,}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Panel C: Binding strength by gene
    ax3 = axes[1, 0]
    
    top_5_genes = df_candidates['gene'].value_counts().head(5).index
    binding_by_gene = []
    
    for gene in top_5_genes:
        gene_data = df_candidates[df_candidates['gene'] == gene]
        sb_count = (gene_data['MT_Affinity_nM'] < 50).sum()
        wb_count = ((gene_data['MT_Affinity_nM'] >= 50) & (gene_data['MT_Affinity_nM'] < 500)).sum()
        binding_by_gene.append({'Gene': gene, 'Strong': sb_count, 'Weak': wb_count})
    
    df_binding = pd.DataFrame(binding_by_gene)
    
    x = np.arange(len(df_binding))
    width = 0.35
    
    ax3.bar(x - width/2, df_binding['Strong'], width, label='Strong Binders (IC50<50nM)',
            color='#06A77D', edgecolor='black', linewidth=1)
    ax3.bar(x + width/2, df_binding['Weak'], width, label='Weak Binders (50-500nM)',
            color='#F18F01', edgecolor='black', linewidth=1)
    
    ax3.set_xticks(x)
    ax3.set_xticklabels(df_binding['Gene'], fontsize=11)
    ax3.set_ylabel('Number of Candidates', fontsize=12, fontweight='bold')
    ax3.set_title('C. Binding Strength by Top 5 Genes', fontsize=13, fontweight='bold')
    ax3.legend()
    ax3.grid(axis='y', alpha=0.3, linestyle='--')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    # Panel D: Quality score by peptide length
    ax4 = axes[1, 1]
    
    length_quality = df_candidates.groupby('peptide_length')['Quality_Score'].apply(list)
    
    bp = ax4.boxplot([length_quality[length] for length in sorted(length_quality.index)],
                     labels=sorted(length_quality.index),
                     patch_artist=True, showmeans=True)
    
    for patch, color in zip(bp['boxes'], ['#2E86AB', '#A23B72', '#F18F01', '#06A77D']):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax4.set_xlabel('Peptide Length', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Quality Score', fontsize=12, fontweight='bold')
    ax4.set_title('D. Quality Score Distribution by Peptide Length', fontsize=13, fontweight='bold')
    ax4.grid(axis='y', alpha=0.3, linestyle='--')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    
    plt.tight_layout()
    output_file = Config.FIGURES_DIR / 'Figure2_Gene_HLA_Analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    logger.success(f"✓ Saved: {output_file}")
    
    # Figure 3: Top Candidates Dashboard
    logger.log("\n[Figure 3] Top candidates summary...")
    
    fig = plt.figure(figsize=(18, 12))
    fig.suptitle('Phase 4 Analysis: Top Neoepitope Candidates Dashboard', 
                 fontsize=18, fontweight='bold', y=0.98)
    
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)
    
    # Key metrics boxes
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.text(0.5, 0.5, f"{len(df_candidates):,}", 
             ha='center', va='center', fontsize=48, fontweight='bold', color='#2E86AB')
    ax1.text(0.5, 0.15, "Total Candidates", 
             ha='center', va='center', fontsize=14, fontweight='bold')
    ax1.axis('off')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    rect = plt.Rectangle((0.05, 0.05), 0.9, 0.9, fill=False, edgecolor='#2E86AB', linewidth=3)
    ax1.add_patch(rect)
    
    ax2 = fig.add_subplot(gs[0, 1])
    sb_count = (df_candidates['MT_Affinity_nM'] < 50).sum()
    ax2.text(0.5, 0.5, f"{sb_count:,}", 
             ha='center', va='center', fontsize=48, fontweight='bold', color='#06A77D')
    ax2.text(0.5, 0.15, "Strong Binders", 
             ha='center', va='center', fontsize=14, fontweight='bold')
    ax2.axis('off')
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    rect = plt.Rectangle((0.05, 0.05), 0.9, 0.9, fill=False, edgecolor='#06A77D', linewidth=3)
    ax2.add_patch(rect)
    
    ax3 = fig.add_subplot(gs[0, 2])
    unique_genes = df_candidates['gene'].nunique()
    ax3.text(0.5, 0.5, f"{unique_genes}", 
             ha='center', va='center', fontsize=48, fontweight='bold', color='#A23B72')
    ax3.text(0.5, 0.15, "Target Genes", 
             ha='center', va='center', fontsize=14, fontweight='bold')
    ax3.axis('off')
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    rect = plt.Rectangle((0.05, 0.05), 0.9, 0.9, fill=False, edgecolor='#A23B72', linewidth=3)
    ax3.add_patch(rect)
    
    # Top 10 candidates table
    ax4 = fig.add_subplot(gs[1:, :])
    ax4.axis('off')
    
    top_10 = df_candidates.nsmallest(10, 'MT_Affinity_nM')[
        ['gene', 'mutation_id', 'mt_peptide', 'HLA', 'MT_Affinity_nM', 
         'WT_Affinity_nM', 'Fold_Change', 'Quality_Score']
    ].copy()
    
    top_10['MT_Affinity_nM'] = top_10['MT_Affinity_nM'].round(2)
    top_10['WT_Affinity_nM'] = top_10['WT_Affinity_nM'].round(2)
    top_10['Fold_Change'] = top_10['Fold_Change'].round(2)
    top_10['Quality_Score'] = top_10['Quality_Score'].round(1)
    
    top_10.columns = ['Gene', 'Mutation', 'Peptide', 'HLA', 'MT IC50', 'WT IC50', 'Fold', 'Quality']
    
    # Truncate mutation_id for display
    top_10['Mutation'] = top_10['Mutation'].str[:20] + '...'
    
    table_data = [top_10.columns.tolist()] + top_10.values.tolist()
    
    table = ax4.table(cellText=table_data, cellLoc='left',
                      loc='center', bbox=[0.05, 0.1, 0.9, 0.8])
    
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    # Style header
    for i in range(len(top_10.columns)):
        table[(0, i)].set_facecolor('#2E86AB')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Alternate row colors
    for i in range(1, len(table_data)):
        for j in range(len(top_10.columns)):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#F0F0F0')
    
    ax4.text(0.5, 0.95, 'Top 10 Neoepitope Candidates by Binding Affinity',
             ha='center', va='top', fontsize=14, fontweight='bold', transform=ax4.transAxes)
    
    plt.savefig(Config.FIGURES_DIR / 'Figure3_Top_Candidates_Dashboard.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    logger.success(f"✓ Saved Figure 3")
    
    logger.success("\n✓ All visualizations completed")

# ============================================================================
# STEP 8: GENERATE FINAL REPORT
# ============================================================================

def generate_final_report(df_candidates, df_all, logger):
    """
    Generate comprehensive Phase 4 analysis report
    """
    logger.section("STEP 8: GENERATING FINAL REPORT")
    
    report_file = Config.PHASE4_DIR / "Phase4_Analysis_Report.txt"
    
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PHASE 4: MHC-I BINDING PREDICTION & IMMUNOGENICITY ASSESSMENT\n")
        f.write("FINAL ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Pipeline Runtime: {logger.get_runtime()}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("1. PREDICTION SUMMARY\n")
        f.write("=" * 80 + "\n")
        f.write(f"Total Peptide-HLA Pairs Analyzed: {len(df_all):,}\n")
        f.write(f"Candidate Neoepitopes Identified: {len(df_candidates):,}\n")
        f.write(f"Filtering Success Rate: {len(df_candidates)/len(df_all)*100:.2f}%\n\n")
        
        f.write("Binding Categories:\n")
        sb = (df_candidates['MT_Affinity_nM'] < 50).sum()
        wb = ((df_candidates['MT_Affinity_nM'] >= 50) & (df_candidates['MT_Affinity_nM'] < 500)).sum()
        f.write(f"  Strong Binders (IC50 < 50 nM): {sb:,} ({sb/len(df_candidates)*100:.1f}%)\n")
        f.write(f"  Weak Binders (IC50 50-500 nM): {wb:,} ({wb/len(df_candidates)*100:.1f}%)\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("2. HLA ALLELE ANALYSIS\n")
        f.write("=" * 80 + "\n")
        f.write("HLA alleles tested:\n")
        for allele in Config.HLA_ALLELES:
            count = (df_candidates['HLA'] == allele).sum()
            f.write(f"  {allele}: {count:,} candidates\n")
        f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("3. TOP 20 GENES BY CANDIDATE COUNT\n")
        f.write("=" * 80 + "\n")
        gene_counts = df_candidates['gene'].value_counts().head(20)
        for i, (gene, count) in enumerate(gene_counts.items(), 1):
            mutations = df_candidates[df_candidates['gene'] == gene]['mutation_id'].nunique()
            sb_gene = (df_candidates[df_candidates['gene'] == gene]['MT_Affinity_nM'] < 50).sum()
            f.write(f"{i:2d}. {gene:10s}: {count:5,} candidates "
                   f"({mutations:4} mutations, {sb_gene:4} SB)\n")
        f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("4. IMMUNOGENICITY METRICS\n")
        f.write("=" * 80 + "\n")
        f.write(f"Mean MT Binding Affinity: {df_candidates['MT_Affinity_nM'].mean():.2f} nM\n")
        f.write(f"Median MT Binding Affinity: {df_candidates['MT_Affinity_nM'].median():.2f} nM\n")
        f.write(f"Mean WT Binding Affinity: {df_candidates['WT_Affinity_nM'].mean():.2f} nM\n")
        f.write(f"Mean Fold Change (WT/MT): {df_candidates['Fold_Change'].mean():.2f}x\n")
        f.write(f"Mean DAI Score: {df_candidates['DAI_Score'].mean():.3f}\n")
        f.write(f"Mean Quality Score: {df_candidates['Quality_Score'].mean():.2f}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("5. TOP 20 CANDIDATE NEOEPITOPES\n")
        f.write("=" * 80 + "\n")
        f.write("Ranked by Quality Score:\n\n")
        
        top_20 = df_candidates.nsmallest(20, 'MT_Affinity_nM')
        f.write(f"{'Rank':<6}{'Gene':<10}{'Peptide':<15}{'HLA':<15}{'MT IC50':<12}{'Fold':<10}{'Quality':<10}\n")
        f.write("-" * 80 + "\n")
        
        for i, (idx, row) in enumerate(top_20.iterrows(), 1):
            f.write(f"{i:<6}{row['gene']:<10}{row['mt_peptide']:<15}{row['HLA']:<15}"
                   f"{row['MT_Affinity_nM']:<12.2f}{row['Fold_Change']:<10.2f}"
                   f"{row['Quality_Score']:<10.1f}\n")
        f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("6. OUTPUT FILES GENERATED\n")
        f.write("=" * 80 + "\n")
        f.write("Candidate Lists:\n")
        f.write(f"  - {Config.CANDIDATES_DIR / 'top_500_neoepitope_candidates.tsv'}\n")
        f.write(f"  - {Config.CANDIDATES_DIR / 'strong_binders_IC50_below_50nM.tsv'}\n")
        f.write(f"  - Gene-specific candidate files (top 10 genes)\n\n")
        
        f.write("Figures:\n")
        f.write(f"  - {Config.FIGURES_DIR / 'Figure1_Binding_Affinity_Analysis.png'}\n")
        f.write(f"  - {Config.FIGURES_DIR / 'Figure2_Gene_HLA_Analysis.png'}\n")
        f.write(f"  - {Config.FIGURES_DIR / 'Figure3_Top_Candidates_Dashboard.png'}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("7. RECOMMENDATIONS FOR EXPERIMENTAL VALIDATION\n")
        f.write("=" * 80 + "\n")
        f.write("Priority 1: Strong Binders (IC50 < 50 nM)\n")
        f.write(f"  - {sb:,} candidates identified\n")
        f.write("  - High likelihood of MHC presentation\n")
        f.write("  - Recommended for in vitro binding assays\n\n")
        
        f.write("Priority 2: High Fold Change (>10x)\n")
        high_fc = (df_candidates['Fold_Change'] > 10).sum()
        f.write(f"  - {high_fc:,} candidates with >10x differential binding\n")
        f.write("  - High tumor specificity\n")
        f.write("  - Recommended for T-cell immunogenicity assays\n\n")
        
        f.write("Priority 3: High Quality Score (>70)\n")
        high_qs = (df_candidates['Quality_Score'] > 70).sum()
        f.write(f"  - {high_qs:,} candidates with Quality Score > 70\n")
        f.write("  - Balanced binding affinity and selectivity\n")
        f.write("  - Recommended for therapeutic vaccine development\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("8. NEXT STEPS\n")
        f.write("=" * 80 + "\n")
        f.write("1. Expression Analysis:\n")
        f.write("   - Correlate candidates with RNA-seq expression data\n")
        f.write("   - Prioritize highly expressed genes\n")
        f.write("   - Filter out low-expression candidates\n\n")
        
        f.write("2. Experimental Validation:\n")
        f.write("   - MHC binding assays (top 50-100 candidates)\n")
        f.write("   - T-cell activation assays (ELISpot, flow cytometry)\n")
        f.write("   - Patient-specific validation with TILs\n\n")
        
        f.write("3. Clinical Translation:\n")
        f.write("   - Vaccine design (peptide or mRNA)\n")
        f.write("   - Manufacturing feasibility assessment\n")
        f.write("   - Regulatory pathway planning\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("PHASE 4 COMPLETE - Ready for Experimental Validation\n")
        f.write("=" * 80 + "\n")
    
    logger.success(f"✓ Final report saved: {report_file}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main pipeline execution"""
    
    print("\n" + "=" * 80)
    print("PHASE 4: MHC-I BINDING PREDICTION & IMMUNOGENICITY ASSESSMENT")
    print("Industry Standard Pipeline (NetMHCpan 4.1)")
    print("=" * 80 + "\n")
    
    # Create directories
    Config.create_directories()
    
    # Initialize logger
    log_file = Config.LOGS_DIR / f"phase4_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    logger = Logger(log_file)
    
    logger.log("Phase 4 Pipeline Started")
    logger.log(f"Output directory: {Config.PHASE4_DIR}")
    
    try:
        # Step 1: Check NetMHCpan installation
        if not check_netmhcpan_installation(logger):
            logger.error("\n" + "=" * 80)
            logger.error("PIPELINE ABORTED: NetMHCpan not installed")
            logger.error("=" * 80)
            logger.error("\nPlease install NetMHCpan 4.1 from:")
            logger.error("https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/")
            logger.error("\nAfter installation, update Config.NETMHCPAN_PATH if needed")
            return
        
        # Step 2: Run NetMHCpan predictions
        logger.section("STEP 2: RUNNING NETMHCPAN PREDICTIONS")
        
        mt_output = Config.PREDICTIONS_DIR / "netmhcpan_mutant_predictions.xls"
        wt_output = Config.PREDICTIONS_DIR / "netmhcpan_wildtype_predictions.xls"
        
        # Run MT predictions
        if not run_netmhcpan_prediction(Config.MUTANT_FASTA, mt_output, "Mutant", logger):
            logger.error("Failed to run mutant predictions. Aborting.")
            return
        
        # Run WT predictions
        if not run_netmhcpan_prediction(Config.WILDTYPE_FASTA, wt_output, "Wildtype", logger):
            logger.error("Failed to run wildtype predictions. Aborting.")
            return
        
        # Parse predictions
        logger.subsection("Parsing Prediction Results")
        mt_predictions = parse_netmhcpan_output(mt_output, logger)
        wt_predictions = parse_netmhcpan_output(wt_output, logger)
        
        if mt_predictions is None or wt_predictions is None:
            logger.error("Failed to parse predictions. Aborting.")
            return
        
        # Step 3: Integrate with Phase 3 data
        df_integrated = integrate_predictions(mt_predictions, wt_predictions, 
                                             Config.PEPTIDE_TABLE, logger)
        
        # Step 4: Calculate immunogenicity scores
        df_scored = calculate_immunogenicity_scores(df_integrated, logger)
        
        # Step 5: Filter candidates
        df_candidates = filter_candidate_neoepitopes(df_scored, logger)
        
        # Step 6: Generate outputs
        generate_candidate_outputs(df_candidates, df_scored, logger)
        
        # Step 7: Create visualizations
        create_phase4_visualizations(df_candidates, df_scored, logger)
        
        # Step 8: Generate report
        generate_final_report(df_candidates, df_scored, logger)
        
        # Final summary
        logger.section("PIPELINE COMPLETED SUCCESSFULLY")
        logger.success(f"Total Runtime: {logger.get_runtime()}")
        logger.success(f"Candidate Neoepitopes: {len(df_candidates):,}")
        logger.success(f"Strong Binders: {(df_candidates['MT_Affinity_nM'] < 50).sum():,}")
        logger.success(f"All outputs saved to: {Config.PHASE4_DIR}")
        
        print("\n" + "=" * 80)
        print("✓ PHASE 4 COMPLETE - Ready for Experimental Validation")
        print("=" * 80 + "\n")
        
    except Exception as e:
        logger.error(f"Pipeline failed with error: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        raise

if __name__ == "__main__":
    main()
