"""
===============================================================================
PHASE 4 MHC-II: CD4+ T-CELL EPITOPE PREDICTION
Complete Pipeline Using YOUR Mutation Data
===============================================================================
"""

import pandas as pd
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime
import requests
import time
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    BASE_DIR = Path("/mnt/c/Users/DELL/neo_project")
    PHASE3_DIR = BASE_DIR / "data/output/Phase3"
    PHASE4_DIR = BASE_DIR / "data/output/Phase4"
    
    MUTATION_SUMMARY = PHASE3_DIR / "analysis_data/mutation_summary.tsv"
    IEDB_MHC2 = Path("/root/software/mhc_ii/mhc_II_binding.py")
    
    MHC2_DIR = PHASE4_DIR / "mhc2"
    MHC2_PEPTIDES = MHC2_DIR / "mhc2_peptides_15mer.tsv"
    MHC2_PREDICTIONS = MHC2_DIR / "mhc2_all_predictions.tsv"
    MHC2_CANDIDATES = MHC2_DIR / "mhc2_candidates.tsv"
    
    # MHC-II alleles (common HLA-DR alleles for population coverage)
    HLA_CLASS2 = [
        "HLA-DRB1*01:01",  # ~20% Caucasian
        "HLA-DRB1*03:01",  # ~20% Caucasian
        "HLA-DRB1*04:01",  # ~25% Caucasian
        "HLA-DRB1*07:01",  # ~25% Caucasian
        "HLA-DRB1*15:01"   # ~20% Caucasian
    ]
    
    PEPTIDE_LENGTH = 15  # Standard for MHC-II
    THRESHOLD_STRONG = 100   # nM
    THRESHOLD_WEAK = 1000    # nM
    
    ENSEMBL_API = "https://rest.ensembl.org"
    
    @classmethod
    def create_directories(cls):
        cls.MHC2_DIR.mkdir(parents=True, exist_ok=True)

Config.create_directories()

print("="*80)
print("PHASE 4 MHC-II: CD4+ T-CELL EPITOPE PREDICTION")
print("Using YOUR Mutation Data from Phase 2")
print("="*80)

# ============================================================================
# STEP 1: LOAD YOUR MUTATION DATA
# ============================================================================

print("\nStep 1: Loading YOUR mutation data...")
df_mutations = pd.read_csv(Config.MUTATION_SUMMARY, sep='\t')
df_mutations = df_mutations[df_mutations['status'] == 'SUCCESS']
print(f"✓ Loaded {len(df_mutations)} mutations from YOUR Phase 2 file")
print(f"  Unique genes: {df_mutations['gene'].nunique()}")

# ============================================================================
# STEP 2: FETCH PROTEIN SEQUENCES & GENERATE 15-MERS
# ============================================================================

print("\nStep 2: Generating 15-mer peptides from YOUR mutations...")

class ProteinFetcher:
    def __init__(self):
        self.cache = {}
    
    def get_sequence(self, gene):
        if gene in self.cache:
            return self.cache[gene]
        
        try:
            # Get gene ID
            url = f"{Config.ENSEMBL_API}/lookup/symbol/homo_sapiens/{gene}"
            response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=10)
            if response.status_code != 200:
                return None
            
            gene_data = response.json()
            gene_id = gene_data['id']
            
            # Get canonical transcript
            url = f"{Config.ENSEMBL_API}/lookup/id/{gene_id}?expand=1"
            response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=10)
            if response.status_code != 200:
                return None
            
            gene_info = response.json()
            
            canonical_transcript = None
            for transcript in gene_info.get('Transcript', []):
                if transcript.get('is_canonical') == 1:
                    canonical_transcript = transcript['id']
                    break
            
            if not canonical_transcript and gene_info.get('Transcript'):
                canonical_transcript = gene_info['Transcript'][0]['id']
            
            if not canonical_transcript:
                return None
            
            # Get protein sequence
            url = f"{Config.ENSEMBL_API}/sequence/id/{canonical_transcript}?type=protein"
            response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=10)
            
            if response.status_code == 200:
                seq = response.json()['seq']
                self.cache[gene] = seq
                return seq
            
            return None
            
        except Exception as e:
            return None

def generate_15mer_peptides(protein_seq, position, ref_aa, alt_aa, mutation_id, gene):
    """Generate 15-mer peptides centered around mutation"""
    peptides = []
    
    if not protein_seq or position is None or position < 1:
        return peptides
    
    pos_idx = position - 1  # Convert to 0-based
    
    if pos_idx < 0 or pos_idx >= len(protein_seq):
        return peptides
    
    # Create mutant sequence
    mutant_seq = protein_seq[:pos_idx] + alt_aa + protein_seq[pos_idx + 1:]
    
    # Generate 15-mers with mutation at different positions (1-15)
    length = Config.PEPTIDE_LENGTH
    
    for offset in range(length):
        start = pos_idx - offset
        end = start + length
        
        if start < 0 or end > len(protein_seq):
            continue
        
        wt_peptide = protein_seq[start:end]
        mt_peptide = mutant_seq[start:end]
        
        if wt_peptide != mt_peptide and len(wt_peptide) == length and len(mt_peptide) == length:
            peptides.append({
                'peptide_id': f"{mutation_id}_15mer_off{offset}",
                'mutation_id': mutation_id,
                'gene': gene,
                'protein_position': position,
                'ref_aa': ref_aa,
                'alt_aa': alt_aa,
                'peptide_length': length,
                'mutation_position_in_peptide': offset + 1,
                'wt_peptide': wt_peptide,
                'mt_peptide': mt_peptide
            })
    
    return peptides

# Generate all 15-mer peptides
fetcher = ProteinFetcher()
all_mhc2_peptides = []

print("Fetching protein sequences and generating peptides...")
unique_genes = df_mutations['gene'].unique()

for gene in unique_genes:
    print(f"  {gene}...", end=" ", flush=True)
    seq = fetcher.get_sequence(gene)
    if seq:
        print(f"✓ ({len(seq)} aa)")
        time.sleep(0.1)  # Rate limiting
    else:
        print("✗ Failed")

for idx, row in df_mutations.iterrows():
    if (idx + 1) % 100 == 0:
        print(f"  Processed {idx + 1}/{len(df_mutations)} mutations...")
    
    gene = row['gene']
    protein_seq = fetcher.get_sequence(gene)
    
    if protein_seq:
        peptides = generate_15mer_peptides(
            protein_seq,
            row['protein_position'],
            row['ref_aa'],
            row['alt_aa'],
            row['mutation_id'],
            gene
        )
        all_mhc2_peptides.extend(peptides)

df_mhc2_peptides = pd.DataFrame(all_mhc2_peptides)

if len(df_mhc2_peptides) == 0:
    print("\n✗ ERROR: No 15-mer peptides generated!")
    exit(1)

df_mhc2_peptides.to_csv(Config.MHC2_PEPTIDES, sep='\t', index=False)
print(f"\n✓ Generated {len(df_mhc2_peptides):,} 15-mer peptides")
print(f"  From {df_mhc2_peptides['mutation_id'].nunique()} mutations")
print(f"  Across {df_mhc2_peptides['gene'].nunique()} genes")

# ============================================================================
# STEP 3: RUN MHC-II BINDING PREDICTIONS
# ============================================================================

print("\n" + "="*80)
print("STEP 3: RUNNING MHC-II BINDING PREDICTIONS")
print("="*80)

print(f"\nProcessing {len(df_mhc2_peptides):,} peptides × {len(Config.HLA_CLASS2)} HLA-II alleles")
print(f"Total predictions: {len(df_mhc2_peptides) * len(Config.HLA_CLASS2):,}")
print("Estimated time: 1-2 hours\n")

mhc2_predictions = []

for hla_idx, allele in enumerate(Config.HLA_CLASS2, 1):
    print(f"[{hla_idx}/{len(Config.HLA_CLASS2)}] {allele}")
    
    # Create temporary FASTA files
    mt_fasta = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta')
    wt_fasta = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta')
    
    try:
        # Write MT peptides in FASTA format
        for idx, row in df_mhc2_peptides.iterrows():
            mt_fasta.write(f">{row['peptide_id']}\n{row['mt_peptide']}\n")
        mt_fasta.close()
        
        # Write WT peptides in FASTA format
        for idx, row in df_mhc2_peptides.iterrows():
            wt_fasta.write(f">{row['peptide_id']}\n{row['wt_peptide']}\n")
        wt_fasta.close()
        
        print(f"  Running MT predictions...", end=" ", flush=True)
        
        # Run MT predictions
        result_mt = subprocess.run(
            ["python3", str(Config.IEDB_MHC2), "netmhciipan", allele, "15", mt_fasta.name],
            capture_output=True, 
            text=True, 
            timeout=3600,
            cwd=str(Config.IEDB_MHC2.parent)
        )
        
        print("✓", end=" ", flush=True)
        print(f"Running WT predictions...", end=" ", flush=True)
        
        # Run WT predictions
        result_wt = subprocess.run(
            ["python3", str(Config.IEDB_MHC2), "netmhciipan", allele, "15", wt_fasta.name],
            capture_output=True, 
            text=True, 
            timeout=3600,
            cwd=str(Config.IEDB_MHC2.parent)
        )
        
        print("✓")
        
        # Parse MT results (column 8 is IC50 for MHC-II)
        mt_results = {}
        for line in result_mt.stdout.strip().split('\n')[1:]:
            parts = line.split('\t')
            if len(parts) >= 9:
                try:
                    seq_num = int(parts[1])
                    ic50 = float(parts[8])  # IC50 column
                    pep_id = df_mhc2_peptides.iloc[seq_num - 1]['peptide_id']
                    mt_results[pep_id] = ic50
                except (ValueError, IndexError):
                    pass
        
        # Parse WT results
        wt_results = {}
        for line in result_wt.stdout.strip().split('\n')[1:]:
            parts = line.split('\t')
            if len(parts) >= 9:
                try:
                    seq_num = int(parts[1])
                    ic50 = float(parts[8])
                    pep_id = df_mhc2_peptides.iloc[seq_num - 1]['peptide_id']
                    wt_results[pep_id] = ic50
                except (ValueError, IndexError):
                    pass
        
        # Combine results
        for idx, row in df_mhc2_peptides.iterrows():
            pep_id = row['peptide_id']
            if pep_id in mt_results:
                mhc2_predictions.append({
                    'peptide_id': pep_id,
                    'gene': row['gene'],
                    'mutation_id': row['mutation_id'],
                    'mt_peptide': row['mt_peptide'],
                    'wt_peptide': row['wt_peptide'],
                    'HLA': allele,
                    'MHC_Class': 'II',
                    'peptide_length': 15,
                    'MT_IC50_nM': mt_results[pep_id],
                    'WT_IC50_nM': wt_results.get(pep_id, 50000),
                    'Fold_Change': wt_results.get(pep_id, 50000) / mt_results[pep_id] if mt_results[pep_id] > 0 else 0
                })
        
        print(f"  ✓ {len([p for p in mhc2_predictions if p['HLA']==allele]):,} predictions")
        
    except subprocess.TimeoutExpired:
        print(f"  ✗ TIMEOUT")
    except Exception as e:
        print(f"  ✗ ERROR: {str(e)}")
    finally:
        Path(mt_fasta.name).unlink(missing_ok=True)
        Path(wt_fasta.name).unlink(missing_ok=True)

print(f"\n✓ MHC-II predictions completed: {len(mhc2_predictions):,}")

if len(mhc2_predictions) == 0:
    print("\n✗ ERROR: No MHC-II predictions generated!")
    exit(1)

# ============================================================================
# STEP 4: PROCESS AND FILTER RESULTS
# ============================================================================

print("\n" + "="*80)
print("STEP 4: PROCESSING MHC-II RESULTS")
print("="*80)

df_mhc2 = pd.DataFrame(mhc2_predictions)

# Calculate quality scores
def calculate_quality_mhc2(row):
    score = 0
    
    # Binding strength (40%)
    if row['MT_IC50_nM'] < Config.THRESHOLD_STRONG:
        score += 40
    elif row['MT_IC50_nM'] < Config.THRESHOLD_WEAK:
        score += 20
    
    # Fold change (30%)
    if row['Fold_Change'] >= 10:
        score += 30
    elif row['Fold_Change'] >= 5:
        score += 20
    elif row['Fold_Change'] >= 2:
        score += 10
    
    # WT non-binding (30%)
    if row['WT_IC50_nM'] > Config.THRESHOLD_WEAK:
        score += 30
    
    return score

df_mhc2['Quality_Score'] = df_mhc2.apply(calculate_quality_mhc2, axis=1)

# Save all predictions
df_mhc2.to_csv(Config.MHC2_PREDICTIONS, sep='\t', index=False)
print(f"✓ All predictions: {Config.MHC2_PREDICTIONS}")

# Filter candidates
df_mhc2_candidates = df_mhc2[
    (df_mhc2['MT_IC50_nM'] < Config.THRESHOLD_WEAK) &
    (df_mhc2['Fold_Change'] >= 2)
].copy()

df_mhc2_candidates = df_mhc2_candidates.sort_values('Quality_Score', ascending=False)

# Save candidates
df_mhc2_candidates.to_csv(Config.MHC2_CANDIDATES, sep='\t', index=False)
print(f"✓ Candidates: {Config.MHC2_CANDIDATES}")

# Save strong binders
df_strong = df_mhc2_candidates[df_mhc2_candidates['MT_IC50_nM'] < Config.THRESHOLD_STRONG]
df_strong.to_csv(Config.MHC2_DIR / "mhc2_strong_binders.tsv", sep='\t', index=False)
print(f"✓ Strong binders: {Config.MHC2_DIR / 'mhc2_strong_binders.tsv'}")

# Gene-specific files
print(f"\nGene-specific files:")
for gene in df_mhc2_candidates['gene'].value_counts().head(10).index:
    gene_df = df_mhc2_candidates[df_mhc2_candidates['gene'] == gene]
    gene_df.to_csv(Config.MHC2_DIR / f"candidates_{gene}_mhc2.tsv", sep='\t', index=False)
    print(f"  {gene}: {len(gene_df):,} candidates")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "="*80)
print("MHC-II RESULTS SUMMARY")
print("="*80)

print(f"\nTotal Predictions: {len(df_mhc2):,}")
print(f"Candidates (IC50<1000nM): {len(df_mhc2_candidates):,} ({len(df_mhc2_candidates)/len(df_mhc2)*100:.1f}%)")
print(f"  Strong Binders (IC50<100nM): {len(df_strong):,}")
print(f"  Weak Binders (100-1000nM): {len(df_mhc2_candidates) - len(df_strong):,}")

print(f"\nUnique Genes with Candidates: {df_mhc2_candidates['gene'].nunique()}")
print(f"Unique Mutations with Candidates: {df_mhc2_candidates['mutation_id'].nunique()}")

print(f"\nTop 10 Genes (MHC-II):")
for gene, count in df_mhc2_candidates['gene'].value_counts().head(10).items():
    print(f"  {gene}: {count:,} candidates")

print(f"\nHLA Allele Distribution:")
for allele, count in df_mhc2_candidates['HLA'].value_counts().items():
    print(f"  {allele}: {count:,} candidates")

# Generate summary report
report_file = Config.MHC2_DIR / "mhc2_analysis_report.txt"
with open(report_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("PHASE 4 MHC-II: CD4+ T-CELL EPITOPE PREDICTION REPORT\n")
    f.write("="*80 + "\n\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    f.write("INPUT DATA:\n")
    f.write(f"  Mutations analyzed: {df_mutations['mutation_id'].nunique()}\n")
    f.write(f"  15-mer peptides generated: {len(df_mhc2_peptides):,}\n")
    f.write(f"  HLA alleles tested: {len(Config.HLA_CLASS2)}\n\n")
    
    f.write("RESULTS:\n")
    f.write(f"  Total predictions: {len(df_mhc2):,}\n")
    f.write(f"  Candidates: {len(df_mhc2_candidates):,}\n")
    f.write(f"  Strong binders: {len(df_strong):,}\n\n")
    
    f.write("TOP 20 MHC-II CANDIDATES:\n")
    f.write(f"{'Gene':<10}{'Peptide':<20}{'HLA':<20}{'MT IC50':<12}{'Fold':<10}{'Quality':<10}\n")
    f.write("-"*80 + "\n")
    for idx, row in df_mhc2_candidates.head(20).iterrows():
        f.write(f"{row['gene']:<10}{row['mt_peptide']:<20}{row['HLA']:<20}"
               f"{row['MT_IC50_nM']:<12.2f}{row['Fold_Change']:<10.2f}{row['Quality_Score']:<10.1f}\n")

print(f"\n✓ Report: {report_file}")

print("\n" + "="*80)
print("✓ PHASE 4 MHC-II COMPLETE")
print("="*80)
print(f"\nAll outputs: {Config.MHC2_DIR}")
print(f"MHC-II Candidates: {len(df_mhc2_candidates):,}")
print(f"Strong Binders: {len(df_strong):,}")
print("\nYour mutation data successfully processed for CD4+ T-cell epitope prediction!")

