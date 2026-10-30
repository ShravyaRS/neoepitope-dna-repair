# ProteinForgeR — Autoregressive Transformer for Protein Sequence Generation

A **comprehensive** protein sequence generator with **structural, functional, and biological validation**. Built in R using `torch`, this pipeline doesn't just generate sequences - it validates them for realistic protein folding, functional motifs, and host interactions.

## Key Features

**STRUCTURAL VALIDATION**
- ESMFold & ColabFold structure prediction
- Folding stability analysis (FoldX, Rosetta, empirical)
- Quality assessment (Ramachandran, radius of gyration, ProSA-like scoring)

**FUNCTIONAL VALIDATION**  
- Active site & catalytic motif detection
- Critical residue conservation analysis
- Signal peptides, localization signals, binding domains
- Influenza-specific validation (HA receptor binding, cleavage sites)

**HOST INTERACTION ANALYSIS**
- Protein-protein interaction prediction (STRING database)
- Domain-based interaction mapping
- Biological pathway analysis (KEGG, Reactome-style)
- Host dependency & pathogenicity risk assessment

**HIGH-THROUGHPUT GENERATION**
- Batched autoregressive transformer 
- Temperature, top-k, nucleus sampling
- Configurable model sizes (2M - 50M+ parameters)

## Architecture

- **Model**: Decoder-only Transformer with causal (autoregressive) masking
- **Tokenization**: Character-level (20 standard amino acids + special tokens)
- **Training Objective**: Next-token prediction (cross-entropy)
- **Validation**: Multi-layered structural, functional, and interaction analysis

## Quick Start

```r
# 1. Install dependencies
source("R/01_install_dependencies.R")

# 2. Place your FASTA file(s) in data/
#    e.g., data/uniprot_sprot.fasta

# 3. Configure validation (edit R/08_run_pipeline.R)
config$protein_family <- "influenza"     # or "general", "enzyme" 
config$organism_type <- "virus"          # or "bacteria", "general"
config$n_structural_analysis <- 50       # Number of structures to predict

# 4. Run the full pipeline
source("R/08_run_pipeline.R")
```

## Validation Pipeline

```
Generated Sequences
        ↓
Sequence Statistics (AA composition, length, k-mers)
        ↓  
Structure Prediction (ESMFold/ColabFold)
        ↓
Stability Assessment (FoldX energy, thermodynamics)
        ↓
Functional Analysis (motifs, active sites, conservation)
        ↓
Host Interactions (PPI, pathways, binding partners)
        ↓
Comprehensive Report (structural quality, functional score, risk assessment)
```

## Project Structure

```
protein-transformer-r/
├── R/
│   ├── 01_install_dependencies.R   # Package installation
│   ├── 02_tokenizer.R              # Amino acid tokenizer
│   ├── 03_data_loader.R            # FASTA parser & batched data loader
│   ├── 04_model.R                  # Transformer architecture (torch)
│   ├── 05_train.R                  # Training loop with logging
│   ├── 06_generate.R               # Sequence generation & sampling
│   ├── 07_evaluate.R               # Quality metrics & analysis
│   ├── 08_run_pipeline.R           # Full end-to-end pipeline
│   ├── 09_structural_prediction.R  # Structure prediction (ESMFold/ColabFold)
│   ├── 10_stability_prediction.R   # Thermodynamic stability analysis
│   ├── 11_functional_validation.R  # Functional motifs & critical residues
│   └── 12_host_interaction_analysis.R # Host-pathogen interactions
├── data/                           # Place your .fasta files here
├── models/                         # Saved model checkpoints
├── outputs/                        # Generated sequences & validation reports
└── README.md
```

## Configuration

Edit the config block in `08_run_pipeline.R`:

```r
config <- list(
  # Model architecture
  d_model = 256, n_heads = 8, n_layers = 6,
  
  # Training  
  n_epochs = 30, batch_size = 32, learning_rate = 3e-4,
  
  # Generation
  n_generate = 500, temperature = 0.8, top_k = 40, top_p = 0.95,
  
  # Validation (NEW!)
  n_structural_analysis = 50,        # Structure prediction count
  n_functional_analysis = 100,       # Functional analysis count
  protein_family = "influenza",      # Family-specific validation
  organism_type = "virus",           # Host interaction analysis
  structure_method = "esmfold"       # Structure predictor
)
```

## Output Files

**Generation:**
- `generated_proteins.fasta` — Generated sequences
- `aa_composition.png` — Amino acid frequency comparison
- `training_curves.png` — Loss and perplexity plots

**Structural Validation:**
- `structures/*.pdb` — Predicted protein structures  
- Quality scores (Ramachandran, compactness, B-factors)

**Functional Analysis:**
- Detected motifs (kinase sites, signal peptides, zinc fingers)
- Active site predictions (catalytic triads, binding pockets)
- Conservation analysis vs. reference sequences

**Host Interactions:**
- Protein-protein interaction predictions
- Biological pathway mappings
- Pathogenicity risk assessment

## External Tool Integration

**Structure Prediction:**
- ESMFold (Meta) — Fast, no installation required
- ColabFold (AlphaFold2) — High accuracy, requires API
- Mock structures for development/testing

**Stability Analysis:**
- FoldX — ΔΔG calculations (optional, requires installation)
- Rosetta — Energy scoring (optional, requires installation)  
- Empirical methods — Built-in, always available

**Databases:**
- STRING — Protein interaction networks
- KEGG/Reactome — Pathway analysis
- Built-in motif databases — No external dependencies

## Why This Matters

Most protein generation models only check sequence statistics (amino acid frequency, length distribution). **ProteinForgeR validates biological realism**:

This catches problems like:
- Sequences that can't fold into stable structures
- Missing critical functional motifs (active sites, binding domains)
- Unrealistic host interaction profiles
- Poor conservation of essential residues

Perfect for computational drug discovery, synthetic biology, and protein engineering research where you need **biologically plausible** sequences, not just statistically similar ones.

## Requirements

- R >= 4.1.0
- torch (R package — installs libtorch automatically)
- Biostrings (Bioconductor, for FASTA parsing)
- ggplot2, dplyr, stringr, bio3d (analysis & visualization)

**Optional External Tools:**
- FoldX (stability prediction)
- Rosetta (energy scoring)
- ESMFold local installation
- HuggingFace API token (for ESMFold API)
