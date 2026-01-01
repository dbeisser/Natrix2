import pandas as pd
import re
import os

# Load input: OTU_table_mumu.csv
mumu_table_file = open(snakemake.input[0], "r")
# Load input: full_table.csv
blast_full_table_file = open(snakemake.input[1], "r")

# Output files
output_all = open(snakemake.output[0], "w")
output_table = open(snakemake.output[1], "w")
output_meta = open(snakemake.output[2], "w")

# Load tables
blast_full = pd.read_csv(blast_full_table_file, index_col=0)
mumu = pd.read_csv(mumu_table_file, sep="\t", index_col=0)

# Handle empty rows in mumu
empty_rows = mumu[mumu.isna().all(axis=1)]

ids_empty_and_not_in_full = [idx for idx in empty_rows.index if idx not in blast_full.index]
ids_empty_and_in_full = [idx for idx in empty_rows.index if idx in blast_full.index]

mumu_clean = mumu.drop(index=ids_empty_and_not_in_full).copy()

# Debug prints
#print("\nEmpty IDs removed (not present in full_table.csv):")
#print(ids_empty_and_not_in_full)
#print("\nEmpty IDs kept (present in full_table.csv):")
#print(ids_empty_and_in_full)

# Create merge keys (prefix up to underscore)
def make_merge_key(id_value):
    match = re.match(r"(.+?_)", id_value)
    return match.group(1) if match else id_value

mumu_clean["merge_key"] = mumu_clean.index.map(make_merge_key)
blast_full["merge_key"] = blast_full.index.map(make_merge_key)

# Check duplicate merge keys
dup_mumu = mumu_clean["merge_key"][mumu_clean["merge_key"].duplicated(keep=False)].unique().tolist()
dup_blast = blast_full["merge_key"][blast_full["merge_key"].duplicated(keep=False)].unique().tolist()

# Debug prints
#if dup_mumu:
    #print("\nDuplicate merge keys in mumu:", dup_mumu)
#if dup_blast:
    #print("\nDuplicate merge keys in full_table:", dup_blast)

# Right join logic (keep all BLAST keys)
keys_blast = set(blast_full["merge_key"])
common_keys = keys_blast

removed_due_to_key_mumu = mumu_clean[~mumu_clean["merge_key"].isin(common_keys)].index.tolist()
removed_due_to_key_blast = []  # nothing removed here

mumu_filtered  = mumu_clean[mumu_clean["merge_key"].isin(common_keys)].copy()
blast_filtered = blast_full[blast_full["merge_key"].isin(common_keys)].copy()

# Group by merge_key and align
all_keys = sorted(common_keys)
mumu_unique   = mumu_filtered.groupby("merge_key").first().reindex(all_keys)
blast_unique  = blast_filtered.groupby("merge_key").first().reindex(all_keys)
blast_matching = blast_unique

# Identify abundance columns (all non-blast info)
blast_info_cols = [
    "sequences","qlen","length","pident","mismatch",
    "qstart","qend","sstart","send","gaps","evalue"]

abundance_columns = [
    col for col in blast_full.columns
    if col not in blast_info_cols + ["taxonomy","merge_key"]]

# Merge abundances (mumu > blast)
combined_abundances = mumu_unique[abundance_columns].combine_first(
    blast_matching[abundance_columns])

blast_info = blast_matching[blast_info_cols]

# Rebuild OTU IDs using final abundance sums
size_values = combined_abundances.sum(axis=1).astype(int)
new_ids = [f"{key}{size_values.loc[key]}" for key in size_values.index]

# Apply new IDs
combined_abundances.index = new_ids
blast_info.index = new_ids

meta_out = blast_matching[
    blast_info_cols + ["taxonomy"]]
meta_out.index = new_ids

# Combine to final output full-table
taxonomy_col = blast_matching["taxonomy"].rename(index=dict(zip(all_keys, new_ids)))

out = pd.concat([blast_info, combined_abundances, taxonomy_col], axis=1)

# Write merged outputs
out.to_csv(output_all, index_label="seqid")
combined_abundances.to_csv(output_table, index_label="seqid")
meta_out.to_csv(output_meta, index_label="seqid")

# Create list of full_table-only OTUs
unmerged_ids = [idx for idx in blast_full.index if idx not in mumu.index]
unmerged_table = blast_full.loc[unmerged_ids].copy()

unmerged_path = os.path.join(
    os.path.dirname(snakemake.output[0]),
    "unmerged_seqids.csv")
unmerged_table.to_csv(unmerged_path, index_label="seqid")

# Debug prints
#print("Unmerged IDs written to:", unmerged_path)
#print("Number of unmerged IDs:", len(unmerged_ids))
#print("\nIDs removed from mumu (merge_key missing):")
#print(removed_due_to_key_mumu)
#print("\nIDs removed from full_table (merge_key missing):")
#print(removed_due_to_key_blast)
#print("\nIDs only present in full_table (unmerged):")
#print(unmerged_ids)
#print("\n")
