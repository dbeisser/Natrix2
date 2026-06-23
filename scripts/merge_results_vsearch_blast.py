import pandas as pd
import logging

# Function: standardizes the sequence ID format for merging
# Supports:
#  - already-normalized IDs like "N279_63" (kept as-is)
#  - vsearch/swarm style IDs like ">123;456;size=..." (converted to "N123_456")
def new_index(table: pd.DataFrame) -> pd.DataFrame:
    table = table.copy()

    if "seqid" not in table.columns:
        raise ValueError("Expected a 'seqid' column but did not find it.")

    ids = table["seqid"].astype(str).str.strip().str.lstrip(">")

    fixed = []
    for x in ids.tolist():
        # If there is no ';' we assume it's already a stable ID (e.g. N279_63, OTU_1, etc.)
        if ";" not in x:
            fixed.append(x)
            continue

        # If there is ';', convert first two fields into N<field0>_<field1>
        parts = x.replace("size=", "").split(";")
        if len(parts) >= 2:
            fixed.append(f"N{parts[0]}_{parts[1]}")
        else:
            # Fallback: keep whatever we have
            fixed.append(parts[0])

    table["seqid"] = fixed
    table = table.set_index("seqid")
    return table


# Load BLAST results and merged abundance table
blast = pd.read_csv(snakemake.input["blast_result"], sep="\t")
merged = pd.read_csv(snakemake.input["merged"], sep=",")

# Normalize the index in both tables
blast = new_index(blast)
merged = new_index(merged)

# Keep only the best BLAST hit per seqid
blast = blast.sort_values(
    by=["evalue", "pident", "length"],
    ascending=[True, False, False]
)

blast = blast[~blast.index.duplicated(keep="first")]

# Merge both tables by seqid
result = merged.join(blast, how="outer")

# Setup logging
logging.basicConfig(filename=str(snakemake.log), level=logging.INFO)

# Identify rows without taxonomic annotation (taxonomy column may not exist if BLAST file is malformed)
if "taxonomy" in result.columns:
    nas = pd.isna(result["taxonomy"])
else:
    # If taxonomy missing, treat all as NA and still write outputs
    nas = pd.Series([True] * len(result), index=result.index)
    logging.warning("Column 'taxonomy' not found in BLAST results; treating all rows as unassigned.")

nbh_counter = int(nas.sum())
bh_counter = int((~nas).sum())
no_blast_hit = list(result.loc[nas, :].index)

logging.info(
    "%s sequences could be assigned taxonomic information using BLAST, %s could not be assigned taxonomic information using BLAST",
    bh_counter,
    nbh_counter,
)
logging.info("No fitting BLAST hit found for seq_id (first 100 shown):")
logging.info(no_blast_hit[:100])

# Define column groups
cols = set(result.columns.tolist())
cols_include = [
    "sequences", "qlen", "length", "pident", "mismatch",
    "qstart", "qend", "sstart", "send", "gaps",
    "evalue", "taxonomy"
]

# Keep only those include-columns that actually exist in the table
cols_include_existing = [c for c in cols_include if c in cols]

# samples = everything else that is not BLAST/meta columns
samples = list(cols - set(cols_include_existing))

# For "complete" table: include include-cols first, then all remaining
cols_all = cols_include_existing + [c for c in result.columns.tolist() if c not in cols_include_existing]

# Save complete result table
result_complete = result[cols_all]
result_complete.to_csv(snakemake.output["complete"], index_label="seqid")

# Save only rows with valid taxonomy (if taxonomy exists)
if "taxonomy" in result.columns:
    result_filtered = result.loc[~nas, :]
else:
    result_filtered = result.iloc[0:0, :]  # empty
result_filtered.to_csv(snakemake.output["filtered"], index_label="seqid")

# Save OTU/ASV abundance table (only sample columns)
result_otus = result[samples]
result_otus.to_csv(snakemake.output["otus"], index_label="seqid")

# Save metadata table (only BLAST-related information)
result_meta = result[cols_include_existing]
result_meta.to_csv(snakemake.output["metadata"], index_label="seqid")
