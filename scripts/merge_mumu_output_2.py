import pandas as pd

# Open input files:
# - MUMU-filtered OTU table
# - Full OTU BLAST annotation table
mumu_file = open(snakemake.input[0], "r")
otu_file = open(snakemake.input[1], "r")

# Open output files:
# - Final merged output (BLAST + MUMU + taxonomy)
# - Filtered MUMU OTU table
# - Filtered BLAST metadata table
# - Duplicates report (all raw duplicates)
output_all = open(snakemake.output[0], "w")
output_table = open(snakemake.output[1], "w")
output_meta = open(snakemake.output[2], "w")
output_dups = open(snakemake.output[3], "w")

# Load input tables as DataFrames
data = pd.read_csv(otu_file, index_col=0)
mumu = pd.read_csv(mumu_file, sep='\t', index_col=0)

# Keep only OTUs that are present in both MUMU and BLAST tables
common_ids = data.index.intersection(mumu.index)
if len(common_ids) == 0:
    raise ValueError("No common OTUs found between MUMU and BLAST input tables.")

# Subset both tables to include only intersecting OTUs
data_inter = data.loc[common_ids]
mumu = mumu.loc[common_ids]

# Identify OTUs that occur more than once (duplicates)
duplicate_ids = data_inter.index[data_inter.index.duplicated(keep=False)].unique()
if len(duplicate_ids) > 0:
    print(f"[natrix2] {len(duplicate_ids)} duplicates saved to 'multi_hits_blast_mumu.csv'")
    # Save all duplicate OTU entries for review
    data.loc[duplicate_ids].to_csv(output_dups)

    # Filter duplicates to keep only those with pident ≥ 90
    valid_dups = data.loc[duplicate_ids]
    valid_dups_filtered = valid_dups[valid_dups["pident"] >= 90.0]

    # Keep only OTUs that have at least one valid (high-pident) hit
    valid_dup_ids = valid_dups_filtered.index.unique()

    # For each duplicate OTU, retain only the best hit (highest pident)
    best_hits = (
        valid_dups_filtered
        .sort_values(by="pident", ascending=False)
        .groupby(level=0)
        .first()
    )

    # Remove all duplicates from the working dataset
    data_inter = data_inter[~data_inter.index.duplicated(keep=False)]

    # Add the selected best hits back to the filtered dataset
    data_inter = pd.concat([data_inter, best_hits])
else:
    print("[natrix2] No duplicate OTUs found.")
    output_dups.write("No duplicate OTUs found.\n")

# Match MUMU abundance table to final list of selected OTUs
mumu = mumu.loc[data_inter.index]

# Merge selected BLAST fields, filtered MUMU abundances, and taxonomy
out = pd.concat([
    data_inter[["sequences", "qlen", "length", "pident", "mismatch", "qstart", "qend", 
                "sstart", "send", "gaps", "evalue"]],
    mumu,
    data_inter["taxonomy"]
], axis=1)

# Export final result table, filtered abundance table, and metadata
out.to_csv(output_all)  # Full merged table
mumu.to_csv(output_table)  # MUMU-filtered OTU counts
data_inter[["sequences", "qlen", "length", "pident", "mismatch", "qstart", "qend", 
            "sstart", "send", "gaps", "evalue", "taxonomy"]].to_csv(output_meta)  # BLAST metadata only
