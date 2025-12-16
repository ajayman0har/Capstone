from Bio import SeqIO
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# ---------- User input ----------
aligned_fasta = "aligned3.fasta"        # your MSA with all sequences
consensus_fasta = "20A_consensus.fasta" # precomputed 20A consensus
threshold = 0.01                        # frequency threshold for mutations
# --------------------------------

# Load 20A consensus as reference
consensus_record = next(SeqIO.parse(consensus_fasta, "fasta"))
ref_seq = str(consensus_record.seq)
print(f"Loaded consensus reference: {consensus_record.id}, length={len(ref_seq)}")

# Read aligned sequences
records = list(SeqIO.parse(aligned_fasta, "fasta"))

# Extract sequences per lineage
seqs_20A = [str(rec.seq) for rec in records if "20A" in rec.description]
seqs_22C = [str(rec.seq) for rec in records if "22C" in rec.description]
seqs_21J = [str(rec.seq) for rec in records if "21J" in rec.description]

print(f"Loaded {len(seqs_20A)} sequences of 20A")
print(f"Loaded {len(seqs_22C)} sequences of 22C")
print(f"Loaded {len(seqs_21J)} sequences of 21J")

# ---------- Functions ----------
def aa_frequencies(seqs):
    """Compute amino acid frequencies per position for a set of sequences."""
    freq_table = []
    seq_len = len(seqs[0])
    for i in range(seq_len):
        col = [s[i] for s in seqs if s[i] != "-"]
        counts = Counter(col)
        total = sum(counts.values())
        freqs = {aa: round(count / total, 3) for aa, count in counts.items()} if total > 0 else {}
        freq_table.append(freqs)
    return freq_table

def difference_matrix_all_freq(ref_seq, freqs_20A, freqs_22C, freqs_21J, threshold=0.01):
    positions = len(ref_seq)
    lineages = ["20A_consensus", "20A", "22C", "21J"]

    aa_matrix = np.empty((4, positions), dtype="<U30")
    color_matrix = np.zeros((4, positions))

    # Row 0 = consensus reference
    aa_matrix[0, :] = list(ref_seq)
    color_matrix[0, :] = 0

    def process_lineage(freqs_list, row_idx):
        for i in range(positions):
            freqs = freqs_list[i] if i < len(freqs_list) else {}
            muts = {aa: f for aa, f in freqs.items() if aa != ref_seq[i] and f >= threshold}
            if muts:
                aa_matrix[row_idx, i] = "\n".join([f"{aa}({round(f,3)})" for aa, f in muts.items()])
                max_f = max(muts.values())
                color_matrix[row_idx, i] = 1 if max_f <= 0.25 else 2 if max_f <= 0.50 else 3 if max_f <= 0.75 else 4
            else:
                aa_matrix[row_idx, i] = ref_seq[i]
                color_matrix[row_idx, i] = 0

    process_lineage(freqs_20A, 1)
    process_lineage(freqs_22C, 2)
    process_lineage(freqs_21J, 3)

    diff_positions = np.where(color_matrix.sum(axis=0) > 0)[0]
    if len(diff_positions) == 0:
        print("No mutations found – skipping heatmap.")
        return None, None, None, None

    return color_matrix[:, diff_positions], aa_matrix[:, diff_positions], lineages, diff_positions + 1

# Compute amino acid frequencies
freqs_20A = aa_frequencies(seqs_20A)
freqs_22C = aa_frequencies(seqs_22C)
freqs_21J = aa_frequencies(seqs_21J)

# Generate heatmap matrices
color_matrix, aa_matrix, lineages, trimmed_positions = difference_matrix_all_freq(
    ref_seq, freqs_20A, freqs_22C, freqs_21J, threshold=threshold
)

# ---------- Plot heatmap ----------
if color_matrix is not None:
    cmap = ListedColormap([
        "#f0f0f0",  # 0 - no mutation
        "#fff7bc",  # 1 - 0–25%
        "#fec44f",  # 2 - 25–50%
        "#fe9929",  # 3 - 50–75%
        "#d95f0e"   # 4 - 75–100%
    ])

    plt.figure(figsize=(max(10, len(trimmed_positions)/12), 4))
    df_plot = pd.DataFrame(color_matrix, index=lineages, columns=trimmed_positions)

    sns.heatmap(
        df_plot,
        cmap=cmap,
        linewidths=0.5,
        linecolor="black",
        annot=aa_matrix,
        fmt="",
        annot_kws={"fontsize": 3, "fontweight": "bold", "ha": "center", "va": "center"},
        yticklabels=lineages,
        xticklabels=trimmed_positions,
        cbar=False
    )
    plt.xticks(fontsize=4)

    legend_elements = [
        Patch(facecolor="#f0f0f0", edgecolor='black', label="0 - No mutation"),
        Patch(facecolor="#fff7bc", edgecolor='black', label="0–25%"),
        Patch(facecolor="#fec44f", edgecolor='black', label="25–50%"),
        Patch(facecolor="#fe9929", edgecolor='black', label="50–75%"),
        Patch(facecolor="#d95f0e", edgecolor='black', label="75–100%")
    ]
    plt.legend(handles=legend_elements, title="Mutation Frequency", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.title("Comparing AA Mutations across 20A, 21J, and 22C")
    plt.xlabel("Amino Acid Position")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig("N_lineage_vs_20A_consensus_heatmap.png", dpi=400)
    plt.show()