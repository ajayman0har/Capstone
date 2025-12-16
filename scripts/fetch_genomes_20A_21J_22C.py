# -*- coding: utf-8 -*-
"""
Created on Sat Dec  6 16:35:06 2025

@author: ajaym
"""

from Bio import Entrez, SeqIO
import pandas as pd
import time

Entrez.email = "manohara@vcu.edu"

# Lineages to fetch
LINEAGE_TERMS = {"20A": "(20A[All Fields] OR B.1[All Fields])","21J_Delta": "(B.1.617.2[All Fields] OR AY*[All Fields])","22C_Omicron": "(BA.2[All Fields] OR BA.2.*[All Fields])"}

MAX_RESULTS = 2000  # per lineage
BATCH_SIZE = 50    # fetch 50 at a time

def fetch_sequences(lineage_name, lineage_term):
    print(f"Searching {lineage_name}")
    
    base_term = (
        "Severe acute respiratory syndrome coronavirus 2[Organism] "
        "AND human[host] "
        "AND 29000:31000[Sequence Length] "
        f"AND {lineage_term}"
    )

    # Search NCBI for nucleotide records
    handle = Entrez.esearch(db="nucleotide", term=base_term, retmax=MAX_RESULTS)
    result = Entrez.read(handle)
    handle.close()

    ids = result["IdList"]
    print(f"Found {len(ids)} records for {lineage_name}")

    n_proteins, metadata = [], []

    # Fetch GenBank records in batches
    for i in range(0, len(ids), BATCH_SIZE):
        batch_ids = ids[i:i + BATCH_SIZE]
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=batch_ids,
                rettype="gb",
                retmode="text"
            )

            for record in SeqIO.parse(handle, "genbank"):
                for f in record.features:
                    if (
                        f.type == "CDS" and
                        "nucleocapsid" in f.qualifiers.get("product", [""])[0].lower()
                    ):
                        seq = f.qualifiers.get("translation", [""])[0]

                        if not seq or len(seq) < 300:
                            continue
                        if "X" in seq or "?" in seq:
                            continue

                        # Add sequence
                        n_proteins.append(
                            f">{record.id}|{lineage_name} {record.description}\n{seq}\n"
                        )

                        # Add metadata
                        metadata.append({
                            "id": record.id,
                            "organism": record.annotations.get("organism", ""),
                            "collection_date": record.annotations.get("date", ""),
                            "protein_length": len(seq),
                            "lineage": lineage_name,
                            "sequence": seq
                        })

                        break

            handle.close()
            time.sleep(1)

        except Exception as e:
            print(f"Skipped batch {i}-{i + len(batch_ids)}: {e}")

    return n_proteins, metadata


# Run for all lineages
all_proteins, all_metadata = [], []

for lineage_name, lineage_term in LINEAGE_TERMS.items():
    proteins, meta = fetch_sequences(lineage_name, lineage_term)
    all_proteins.extend(proteins)
    all_metadata.extend(meta)

# Save
with open("sarscov2_N_proteins_20A_21J_22C_subset.fasta", "w") as f:
    f.writelines(all_proteins)

# Save metadata
pd.DataFrame(all_metadata).to_csv("sarscov2_N_metadata_20A_21J_22C_subset.csv",index=False)

print("\nSaved sarscov2_N_proteins_20A_21J_22C_subset.fasta")
print("Saved sarscov2_N_metadata_20A_21J_22C_subset.csv")