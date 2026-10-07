# SARS-CoV-2 Nucleocapsid Mutation Analysis

Comparative analysis of ~5,600 SARS-CoV-2 Nucleocapsid (N) protein sequences across three lineages — ancestral B.1 (20A), Delta (21J), and Omicron BA.2 — to identify convergent mutations and assess their potential impact on antigen test sensitivity.

*VCU Bioinformatics Capstone*



## Key Findings
- Identified **5 convergent N-protein mutations** arising independently across lineages
- Mapped these mutations to **B-cell epitope regions** targeted by antigen tests

## Pipeline
| Step | Script | Description |
|------|--------|-------------|
| 1 | `fetch_genomes_20A_21J_22C.py` | Queries NCBI Nucleotide for up to 2,000 records per lineage, extracts N-protein translations (≥300 aa, no ambiguous residues), and writes a FASTA plus metadata CSV |
| 2 | MAFFT | Protein multiple sequence alignment → `aligned3.fasta` |
| 3 | `consensus.py` | Builds a majority-rule 20A consensus (excluding sequences >50% gaps) as the ancestral baseline → `20A_consensus.fasta` |
| 4 | `mutation_comparison.py` | Computes per-lineage amino-acid frequencies, flags substitutions vs. the 20A consensus at ≥1% frequency, and plots a frequency heatmap |

Lineages are assigned by NCBI search terms (Pango designations in record metadata), not by re-classification with Pangolin/Nextclade.

## Tech Stack
Python · Biopython · MAFFT · Pandas · NumPy · Matplotlib · Seaborn

## Usage
```bash
git clone https://github.com/ajayman0har/Capstone.git
cd Capstone/scripts
pip install biopython pandas numpy matplotlib seaborn
# MAFFT must be installed separately: https://mafft.cbrc.jp/alignment/software/
# Set your own email in fetch_genomes_20A_21J_22C.py (required by NCBI Entrez)

python fetch_genomes_20A_21J_22C.py
mafft --auto sarscov2_N_proteins_20A_21J_22C_subset.fasta > aligned3.fasta
python consensus.py
python mutation_comparison.py
```

To skip fetching and alignment, run the last two commands on the included `aligned3.fasta`.

## Output
- `20A_consensus.fasta`: ancestral N-protein consensus
- `N_lineage_vs_20A_consensus_heatmap.png`: substitutions per lineage, colored by frequency

## License
MIT
