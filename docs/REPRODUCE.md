# Reproducing the analysis

Install the environment first (Python 3.13, versions pinned in `requirements.txt`):

    pip install -r requirements.txt

External predictors (not redistributable here): NetMHCpan 4.1 and MHCflurry 2.0,
orchestrated through pVACseq for MHC class I; the IEDB MHC class II tools
(NetMHCIIpan) for class II. The prediction outputs they produced are deposited
in this record, so every step below can be rerun without them.

## Order of execution

    python phase4_netmhcpan_pipeline.py        # class I predictions  -> 04_mhc1_predictions.tsv
    python phase4_mhc2_complete.py             # class II predictions -> 05_mhc2_predictions.tsv
    python phase5_detailed_neoepitope_ranking.py
    python phase5_final_neoepitope_selection.py
    python Phase12_Advanced.py                 # bipartite network    -> 06_network_degrees.csv
    python PHASE12_WOW_FACTOR.py               # composite score      -> 09_top_ranked_peptides.csv
    python PHASE13_FINAL_WOW.py
    python PHASE10_COMPREHENSIVE_COMBINATION_METAANALYSIS.py

## Composite score

Implemented in `PHASE12_WOW_FACTOR.py`:

    score = 40*s(IC50) + 20*s(FC) + 20*s(share) + 20*s(gene)

    s(IC50) = 1 - (IC50 - 2) / (49328 - 2)          linear min-max
    s(FC)   = min(fold change, 100) / 100            hard cap at 100
    s(share)= min(prediction rows, 100) / 100        hard cap at 100
    s(gene) = TP53 1.0; BRCA1/BRCA2 0.95; ATM 0.9;
              MSH2/MLH1 0.85; POLE/ATR 0.8; all others 0.7

The weights are fixed constants in the script. No calibration step was performed.

## Supplementary tables

`Supplementary_Table_S1.csv` ... `Supplementary_Table_S9.csv` are the tables
cited in the manuscript. S5 carries the per-peptide composite score with its four
components broken out. S9 lists the columns of every file in this record.

## Definitions

- Degree / Connectivity: number of distinct parental sample-level calls that
  generate a peptide. Hubs are peptides of degree > 5 (n = 836 of 33,405).
- `Network_Strength` (file 06) and `epitope_sharing` (file 09) are prediction-row
  counts, i.e. degree x 8 alleles. The manuscript reports degree.
- Mutation counts are recurrence counts: 1,419 sample-level calls = 885 distinct
  genomic variants across 877 samples. 1,417 of the 1,419 produced peptide windows.
