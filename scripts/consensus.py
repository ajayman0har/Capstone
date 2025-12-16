# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 16:52:55 2025

@author: ajaym
"""

from Bio import SeqIO
from collections import Counter


aligned_fasta = "aligned3.fasta"   
lineage_tag = "20A"                
max_gap_frac = 0.5                 
min_consensus_freq = 0.5           
output_fasta = "20A_consensus.fasta"
# Read alignment
records = list(SeqIO.parse(aligned_fasta, "fasta"))
print(f"Total sequences in alignment: {len(records)}")

# 20A sequences
seqs_20A_all = [str(rec.seq) for rec in records if lineage_tag in rec.description]
print(f"Found {len(seqs_20A_all)} sequences with tag '{lineage_tag}'")

# too many gaps
seqs_20A = [s for s in seqs_20A_all if s.count("-") / len(s) <= max_gap_frac]
print(f"{len(seqs_20A)} sequences remaining after gap filtering (max {max_gap_frac*100:.0f}%)")

if len(seqs_20A) == 0:
    raise SystemExit("No sequences left after filtering. Check headers or max_gap_frac.")

# Compute consensus
seq_len = len(seqs_20A[0])
assert all(len(s) == seq_len for s in seqs_20A), "Sequences are not the same length!"

consensus = []
for i in range(seq_len):
    col_counts = Counter([s[i] for s in seqs_20A if s[i] != "-"])
    total = sum(col_counts.values())
    if total == 0:
        consensus.append("-")  # fully gapped
        continue
    most_common, count = col_counts.most_common(1)[0]
    freq = count / total
    if freq >= min_consensus_freq:
        consensus.append(most_common)
    else:
        consensus.append(most_common)

consensus_str = "".join(consensus)
print(f"Consensus length: {len(consensus_str)}")
print(f"Number of gaps in consensus: {consensus_str.count('-')}")
print(f"Consensus sequence:\n{consensus_str}")


with open(output_fasta, "w") as f:
    f.write(f">20A_consensus\n{consensus_str}\n")

print(f"Consensus saved to {output_fasta}")