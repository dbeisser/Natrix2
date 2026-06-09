import pandas as pd
import re


# Input and output files from Snakemake
#
# snakemake.input[0] = mothur taxonomy assignment file (*.taxonomy)
# snakemake.input[1] = abundance table / ASV table / OTU table
# snakemake.output[0] = full table with counts + taxonomy
# snakemake.output[1] = count table only
# snakemake.output[2] = metadata table with sequence + taxonomy information

TAXONOMY_FILE = snakemake.input[0]
OTU_TABLE_FILE = snakemake.input[1]

FULL_TABLE_OUT = snakemake.output[0]
OTU_TABLE_OUT = snakemake.output[1]
METADATA_TABLE_OUT = snakemake.output[2]


# Helper function for robust ID matching
#
# Mothur may write sequence IDs without the leading ">", while FASTA-derived
# ASV/OTU tables may still contain it. Example:
#   filtered_table.csv:   >1;size=584;
#   mothur taxonomy:      1;size=584;
# These two strings are biologically the same ID, but pandas would not match
# them directly. Therefore, we create a temporary normalized merge key.

def normalize_seqid(seqid):
    """Normalize sequence IDs for merging taxonomy and abundance tables."""
    if pd.isna(seqid):
        return seqid
    return str(seqid).strip().lstrip(">").strip()


# Read input files
#
# Mothur taxonomy output has no header and is tab-separated:
#   seqid <tab> taxonomy
# Keep seqid as string to avoid accidental type conversion.
taxonomy = pd.read_csv(
    TAXONOMY_FILE,
    sep="\t",
    names=["seqid", "taxonomy"],
    dtype={"seqid": str, "taxonomy": str},
)

# Read abundance / ASV / OTU table. Only seqid is forced to string; count columns
# should keep their numeric type.
otu_count = pd.read_csv(
    OTU_TABLE_FILE,
    converters={"seqid": str},
)


# Validate required columns
if "seqid" not in otu_count.columns:
    raise ValueError("Input abundance table must contain a 'seqid' column.")

if "seqid" not in taxonomy.columns or "taxonomy" not in taxonomy.columns:
    raise ValueError("Mothur taxonomy table must contain 'seqid' and 'taxonomy'.")


# Merge taxonomy into the abundance table
# Create temporary normalized IDs only for matching. The original seqid from the
# abundance table is preserved in the final output.
otu_count["merge_id"] = otu_count["seqid"].apply(normalize_seqid)
taxonomy["merge_id"] = taxonomy["seqid"].apply(normalize_seqid)

# Keep only the columns needed for the merge and avoid duplicate taxonomy rows
# causing duplicated abundance rows.
taxonomy_for_merge = taxonomy[["merge_id", "taxonomy"]].drop_duplicates(
    subset="merge_id",
    keep="first",
)

final_file = pd.merge(
    otu_count,
    taxonomy_for_merge,
    on="merge_id",
    how="left",
)

# Remove the temporary merge key from the final table.
final_file = final_file.drop(columns=["merge_id"])

# Write full table: abundance table plus taxonomy column.
final_file.to_csv(FULL_TABLE_OUT, index=False)


# Create metadata table
#
# The metadata table contains sequence ID, sequence string if available, and
# taxonomy. This keeps the previous output structure but is more robust if a
# future table lacks the 'sequences' column.
metadata_columns = [col for col in ["seqid", "sequences", "taxonomy"] if col in final_file.columns]
out2 = final_file[metadata_columns].copy()

# Mothur can include bootstrap/confidence values in parentheses, e.g.:
#   Bacteria(100);Firmicutes(95);
# Add a cleaned taxonomy2 column without these values when taxonomy exists.
if "taxonomy" in out2.columns:
    out2["taxonomy2"] = out2["taxonomy"].apply(
        lambda x: x if pd.isna(x) else re.sub(r"\([^()]*\)", "", str(x))
    )

out2.to_csv(METADATA_TABLE_OUT, index=False)


# Create count table only
#
# Remove non-count metadata columns if present. The original script removed only
# 'sequences'; dropping with errors='ignore' keeps this compatible with ASV and
# OTU tables that may differ slightly in structure.
out3 = otu_count.drop(columns=["sequences", "merge_id"], errors="ignore")
out3.to_csv(OTU_TABLE_OUT, index=False)
